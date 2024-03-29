  subroutine manymd0(emany,lgrad1)
!
!  Subroutine for calculating the many-body energy from the
!  Sutton-Chen potential(s). Requires real space routine to
!  be called first to evaluate the repulsive contribution
!  and the rho parameters.
!
!  Finite clusters only version
!
!  On entry the array scrho must contain the density at each atomic site.
!
!  lfreeze = if .true. only do derivatives for atoms with a degree of
!            freedom
!
!   3/99 Created from many0d.f specifically for first derivative
!        only calculation -> MD
!   3/99 Parallel modifications added
!  10/99 Cubic density function added
!  11/02 Parallel changes made
!  11/03 ndennat/ndentyp replaced
!   7/05 Constant scaling added to sqrt and power functionals
!   7/05 Use of EAM species pointers introduced
!   9/05 rhoderv called to get density derivatives
!   9/05 Call to eamfnderv used to replace EAM function derivative code
!  10/05 Call to eamfnderv added for energy
!   4/06 Modified to handle species specific densities
!   3/07 Printing of EAM densities and energies added as an option
!   3/07 Calculation of emany parallelised
!   5/07 QM/MM schemes added
!   5/07 Call to rhoderv modified by adding rpot
!  12/07 Unused variables removed
!  11/08 Call to rhoderv updated to include the density
!  11/08 Call to rhoderv modified to include x,y,z Cartesian components
!  11/08 rho arrays changed to 2-D to benefit MEAM
!  11/08 call to rhoderv replaced by calls to meamrho/eamrho according to MEAM vs EAM
!  12/08 rho switched back to 1-D array with condensed components
!  12/08 scrho, and defect versions, converted to 2-D array for MEAM
!  12/08 MEAM calculation of density added as option
!   1/09 Total density now printed for MEAM case
!   1/09 Derivatives modified to accommodate MEAM : deriv -> deriv(3)
!   2/09 Derivatives for MEAM added
!   3/09 small replaced by global value smallself from general module
!   3/09 Sign of return arguments from EAM/MEAM function routines corrected for
!   4/09 MEAM screening function derivatives added
!   5/09 MEAM third derivative arguments removed
!   6/09 Site energy and virial added
!  11/09 Region derivatives added
!   4/12 Explicit virial calculation removed as no longer needed
!   4/12 xvir, yvir and zvir removed
!   8/14 MEAM screening made species specific
!   8/14 Taper range passed to eamrho/meamrho
!   8/14 Pair potential derivatives added here so that screening can be included
!   8/14 All MEAM species (i,j,k) now passed to meanscreen
!   2/15 Cycling of loops for zero density removed for lMEAMden case
!   2/15 lvalidij set to be true by baskes potential
!   2/18 Trace added
!
!  Conditions of use:
!
!  GULP is available free of charge to academic institutions
!  and non-commerical establishments only. Copies should be
!  obtained from the author only and should not be distributed
!  in any form by the user to a third party without the express
!  permission of the author. This notice applies to all parts
!  of the program, except any library routines which are
!  distributed with the code for completeness. All rights for
!  such routines remain with the original distributor.
!
!  No claim is made that this program is free from errors and
!  no liability will be accepted for any loss or damage that
!  may result. The user is responsible for checking the validity
!  of their results.
!
!  Copyright Curtin University 2018
!
!  Julian Gale, CIC, Curtin University, February 2018
!
  use configurations, only : nregionno, nregiontype, QMMMmode
  use control
  use current
  use derivatives
  use eam
  use energies,       only : siteenergy
  use general,        only : smallself
  use iochannels,     only : ioout
  use mdlogic
  use optimisation
  use parallel
  use sutton
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  use two
  implicit none
!
!  Passed variables
!
  logical,     intent(in)                      :: lgrad1
  real(dp),    intent(inout)                   :: emany
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: indij
  integer(i4)                                  :: j
  integer(i4)                                  :: k
  integer(i4)                                  :: m
  integer(i4)                                  :: mj
  integer(i4)                                  :: n  
  integer(i4)                                  :: nat1
  integer(i4)                                  :: nat2
  integer(i4)                                  :: nati
  integer(i4)                                  :: natj
  integer(i4)                                  :: neamfnspec2
  integer(i4)                                  :: neamspeci
  integer(i4)                                  :: neamspecj
  integer(i4)                                  :: neamspeck
  integer(i4)                                  :: noff
  integer(i4)                                  :: noffm1
  integer(i4)                                  :: noffset
  integer(i4)                                  :: np
  integer(i4)                                  :: npartial
  integer(i4)                                  :: npot
  integer(i4)                                  :: npots
  integer(i4)                                  :: nregioni
  integer(i4)                                  :: nregionj
  integer(i4)                                  :: nregionk
  integer(i4)                                  :: nregiontypi
  integer(i4)                                  :: nregiontypj
  integer(i4)                                  :: ntyp1  
  integer(i4)                                  :: ntyp2
  integer(i4)                                  :: ntypi
  integer(i4)                                  :: ntypj    
  integer(i4), dimension(:), allocatable       :: npotl
  integer(i4)                                  :: status
  logical                                      :: lnonzeroSij
  logical                                      :: lpartial
  logical                                      :: lopi
  logical                                      :: lopj
  logical                                      :: lopk
  logical                                      :: lQMMMok
  logical                                      :: lvalidij
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: cut2
  real(dp)                                     :: cut2p
  real(dp)                                     :: cut2r
  real(dp)                                     :: deriv(3)
  real(dp)                                     :: derivs(6)
  real(dp)                                     :: drhoij(3,maxmeamcomponent)
  real(dp)                                     :: drhoijs(6,maxmeamcomponent)
  real(dp)                                     :: drhoij2(6,maxmeamcomponent)
  real(dp)                                     :: drhoij2s(21,maxmeamcomponent)
  real(dp)                                     :: drhoij2m(6,3,maxmeamcomponent)
  real(dp)                                     :: drhoji(3,maxmeamcomponent)
  real(dp)                                     :: drhojis(6,maxmeamcomponent)
  real(dp)                                     :: drhoji2(6,maxmeamcomponent)
  real(dp)                                     :: drhoji2s(21,maxmeamcomponent)
  real(dp)                                     :: drhoji2m(6,3,maxmeamcomponent)
  real(dp)                                     :: drhototij(3)
  real(dp)                                     :: drhototijs(6)
  real(dp)                                     :: drhototij2(6)
  real(dp)                                     :: drhototij2s(21)
  real(dp)                                     :: drhototij2m(6,3)
  real(dp)                                     :: drhototij3(10)
  real(dp)                                     :: drhototji(3)
  real(dp)                                     :: drhototjis(6)
  real(dp)                                     :: drhototji2(6)
  real(dp)                                     :: drhototji2s(21)
  real(dp)                                     :: drhototji2m(6,3)
  real(dp)                                     :: drhototji3(10)
  real(dp)                                     :: ebas
  real(dp)                                     :: d1bas
  real(dp)                                     :: d2bas
  real(dp)                                     :: eeam
  real(dp)                                     :: emanytrm
  real(dp)                                     :: oci      
  real(dp)                                     :: ocj  
  real(dp)                                     :: ofct
  real(dp)                                     :: r   
  real(dp)                                     :: r2
  real(dp)                                     :: r2ijmid
  real(dp)                                     :: r2ik
  real(dp)                                     :: r2jk
  real(dp)                                     :: rcut2
  real(dp)                                     :: rcutfactor
  real(dp)                                     :: rhoi
  real(dp)                                     :: rhoj
  real(dp)                                     :: rhoij(maxmeamcomponent)
  real(dp)                                     :: rhoji(maxmeamcomponent)
  real(dp)                                     :: rk
  real(dp)                                     :: rp
  real(dp)                                     :: rscrhoi
  real(dp)                                     :: rscrhoi3
  real(dp)                                     :: rscrhoi5
  real(dp)                                     :: rscrhoj
  real(dp)                                     :: rscrhoj3
  real(dp)                                     :: rscrhoj5
  real(dp)                                     :: scmax
  real(dp)                                     :: Sij
  real(dp)                                     :: Sikj
  real(dp)                                     :: dSikjdr(3)
  real(dp)                                     :: time1
  real(dp)                                     :: time2  
  real(dp)                                     :: xal 
  real(dp)                                     :: yal    
  real(dp)                                     :: zal
  real(dp)                                     :: xcd 
  real(dp)                                     :: ycd    
  real(dp)                                     :: zcd
  real(dp)                                     :: xij0
  real(dp)                                     :: yij0
  real(dp)                                     :: zij0
  real(dp)                                     :: xik0
  real(dp)                                     :: yik0
  real(dp)                                     :: zik0
  real(dp)                                     :: xjk0
  real(dp)                                     :: yjk0
  real(dp)                                     :: zjk0
  type(screening_atoms)                        :: partial
#ifdef TRACE
  call trace_in('manymd0')
#endif
!
  time1 = g_cpu_time()
!
!  For screened MEAM, set scale factor that determines searching cutoff based on which axis of the ellipse is largest.
!  Note: the cutoff is applied to the mid point of the i-j vector and so this is guaranteed to find all distances.
!
  rcutfactor = 1.0_dp
  if (lanyMEAMscreen) then
    neamfnspec2 = neamfnspec*(neamfnspec+1)/2
    do i = 1,neamfnspec
      do j = 1,neamfnspec2
        if (lMEAMscreen(j,i)) then
          rcutfactor = max(1.0_dp,0.25_dp*(1.0_dp + meam_Cmax(j,i)))
        endif
      enddo
    enddo
  endif
!
!  Scale density
!
  call eamscalescrho(1_i4)
  if (lPrintEAM) then
    call mpbarrier
    if (ioproc) then
!
!  Openning banner for energy decomposition
!
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  EAM : Atom No.                Density                 Atom energy (eV) '')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    endif
  endif
!
!  Energy calculation
!
  emany = 0.0_dp
  do i = procid+1,numat,nprocs
    neamspeci = neamfnspecptr(i)
    rhoi = scrho(1,i)
    nregioni = nregionno(nsft+i)
    nregiontypi = nregiontype(nregioni,ncf)
!
!  QM/MM handling : i is a QM atom => exclude
!
    lQMMMok = .true.
    if (QMMMmode(ncf).gt.0) then
      if (nregiontypi.eq.1) lQMMMok = .false.
    endif
    if (neamspeci.gt.0.and.rhoi.gt.1.0d-12.and.lQMMMok) then
      if (lMEAMfn) then
        call meamfnderv(neamfn,neamspeci,scrho(1,i),rhoi,eeam,rscrhoi,rscrhoi3,rscrhoi5,.false.,.false.,.false.)
      else
        call eamfnderv(neamfn,neamspeci,scrho(1,i),eeam,rscrhoi,rscrhoi3,rscrhoi5,.false.,.false.,.false.)
      endif
      emanytrm = occuf(i)*eeam
      emany = emany + emanytrm
      siteenergy(i) = siteenergy(i) + emanytrm
      if (lPrintEAM) then
        write(ioout,'(7x,i8,4x,f24.10,4x,f24.12)') i,rhoi,emanytrm
      endif
    endif
  enddo
  if (lPrintEAM) then
    call mpbarrier
    if (ioproc) then
!
!  Closing banner for energy decomposition
!
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    endif
  endif
!
!  If no forces are needed then we don't need to do loops over atoms, so just return
!
  if (.not.lgrad1) goto 1000
!
!  From here on we can assume that lgrad1  =  .true.
!
!  Allocate local memory
!
  allocate(npotl(npote),stat=status)
  if (status/=0) call outofmemory('manymd0','npotl')
!
!  Find maximum cut-off radius
!
  scmax = 0.0_dp
  do i = 1,npote
    if (nptype(i).eq.19) then
      if (rpot(i).gt.scmax) scmax = rpot(i)
    endif
  enddo
!
!  Set up cutoffs
!
  cut2p = cutp*cutp
  if (lnoreal) goto 999
!
!  Outer loop over sites
!
!  Use Brode-Ahlrichs Algorithm
  noff = numat/2
  noffset  = noff
  if (mod(numat,2_i4).eq.0) then
    noffm1 = noff - 1
  else
    noffm1 = noff
  endif
  iloop: do i = procid+1,numat,nprocs
    if (i.gt.noff) noffset = noffm1
!orig      do 10 i = 2,numat
    lopi = (.not.lfreeze.or.lopf(i))
!     
!  Find EAM species for i
!  
    neamspeci = neamfnspecptr(i)
!
!  Evaluate functional derivatives
!
    if (lMEAMfn) then
      call meamfnderv(neamfn,neamspeci,scrho(1,i),rhoi,eeam,rscrhoi,rscrhoi3,rscrhoi5,.true.,.false.,.false.)
    else
      call eamfnderv(neamfn,neamspeci,scrho(1,i),eeam,rscrhoi,rscrhoi3,rscrhoi5,.true.,.false.,.false.)
      rhoi = scrho(1,i)
    endif
!
    xal = xclat(i)
    yal = yclat(i)
    zal = zclat(i)
    nati = nat(i)
    ntypi = nftype(i)
    nregioni = nregionno(nsft+i)
    nregiontypi = nregiontype(nregioni,ncf)
    oci = occuf(i)
!
!  Start of second atom loop
!
    jloop: do mj = 1,noffset
      j = mod(i+mj-1_i4,numat) + 1
!orig        do 20 j = 1,i-1
      lopj = (.not.lfreeze.or.lopf(j))
      if (.not.lopi.and..not.lopj) cycle jloop
!     
!  Find EAM species for j
!  
      neamspecj = neamfnspecptr(j)
!
!  Find index for i-j
!
      if (neamspeci.gt.neamspecj) then
        indij = neamspeci*(neamspeci - 1)/2 + neamspecj
      else
        indij = neamspecj*(neamspecj - 1)/2 + neamspeci
      endif
!
!  Evaluate functional derivatives
!
      if (lMEAMfn) then
        call meamfnderv(neamfn,neamspecj,scrho(1,j),rhoj,eeam,rscrhoj,rscrhoj3,rscrhoj5,.true.,.false.,.false.)
      else
        call eamfnderv(neamfn,neamspecj,scrho(1,j),eeam,rscrhoj,rscrhoj3,rscrhoj5,.true.,.false.,.false.)
        rhoj = scrho(1,j)
      endif
!
!  Skip if this is not a MEAM density and the densities are zero
!
      if (.not.lMEAMden.and.rhoi.eq.0.0_dp.and.rhoj.eq.0.0_dp) cycle jloop
!
      natj = nat(j)
      ntypj = nftype(j)
      nregionj = nregionno(nsft+j)
      nregiontypj = nregiontype(nregionj,ncf)
!
!  QM/MM handling : i & j are both QM atoms => no forces to compute
!
      if (QMMMmode(ncf).gt.0) then
        if (nregiontypi.eq.1.and.nregiontypj.eq.1) cycle jloop
      endif
!
!  Arrange nat1 to be lower than nat2. If they are equal
!  then arrange that ntyp1 is lower than ntyp2
!
      if (nati.eq.natj) then
        nat1 = nati
        nat2 = natj
        if (ntypi.lt.ntypj) then
          ntyp1 = ntypi
          ntyp2 = ntypj
        else
          ntyp1 = ntypj
          ntyp2 = ntypi
        endif
      elseif (nati.lt.natj) then
        nat1 = nati
        nat2 = nat(j)
        ntyp1 = ntypi
        ntyp2 = nftype(j)
      else
        nat1 = nat(j)
        nat2 = nati
        ntyp1 = nftype(j)
        ntyp2 = ntypi
      endif
      xcd = xclat(j) - xal
      ycd = yclat(j) - yal
      zcd = zclat(j) - zal
      ocj = occuf(j)
      ofct = oci*ocj
!
!  Locate potential number
!  Check whether potential requires specific types
!  Calculate sum of all dispersion terms for pair
!
      npots = 0
      rp = 0.0_dp
      do n = 1,npote
        if (nptype(n).eq.19.or.nptype(n).eq.45.or.nptype(n).eq.55) then
          if (nat1.eq.nspec1(n).and.nat2.eq.nspec2(n)) then
            if ((ntyp1.eq.nptyp1(n).or.nptyp1(n).eq.0).and.(ntyp2.eq.nptyp2(n).or.nptyp2(n).eq.0)) then
              npots = npots + 1
              npotl(npots) = n
              if (rpot(n).gt.rp) rp = rpot(n)
            endif
          endif
        endif
      enddo
!
!  If no valid potentials then there is no need
!  to continue with this pair, unless this is
!  a second derivative calculation, in which
!  case there may be a contribution from triangles
!  of interactions.
!
      if (npots.eq.0) cycle jloop
      cut2r = rp*rp
      if (cut2r.gt.cut2p) cut2r = cut2p
      cut2 = cut2r
      rp = sqrt(cut2)
      r2 = xcd*xcd + ycd*ycd + zcd*zcd
      if (r2.gt.smallself.and.r2.le.cut2) then
!***************************************************
!  Calculate many-body contribution in real space  *
!***************************************************
        deriv(1:3) = 0.0_dp
        if (lmd) derivs(1:6) = 0.0_dp
        r = sqrt(r2)
!***************************************
!  Valid many-body potentials for i-j  *
!***************************************
        if (lMEAM) then
          rhoij(1:maxmeamcomponent) = 0.0_dp
          rhoji(1:maxmeamcomponent) = 0.0_dp
          drhoij(1:3,1:maxmeamcomponent) = 0.0_dp
          drhoji(1:3,1:maxmeamcomponent) = 0.0_dp
          if (lmd) then
            drhoijs(1:6,1:maxmeamcomponent) = 0.0_dp
            drhojis(1:6,1:maxmeamcomponent) = 0.0_dp
          endif
        else
          rhoij(1) = 0.0_dp
          rhoji(1) = 0.0_dp
        endif
        drhototij(1:3) = 0.0_dp
        drhototji(1:3) = 0.0_dp
        if (lmd) then
          drhototijs(1:6) = 0.0_dp
          drhototjis(1:6) = 0.0_dp
        endif
        lvalidij = .false.
        if (npots.gt.0) then
          if (lMEAMden) then
            do m = 1,npots
              npot = npotl(m)
              if (nptype(npot).eq.19) then
                if (r.gt.rpot2(npot).and.r.le.rpot(npot).and.r.le.rpmax) then
                  lvalidij = .true.
!**********************************
!  Calculate density derivatives  *
!**********************************
                  call meamrho(nati,ntypi,natj,ntypj,r,rpot(npot),xcd,ycd,zcd,rhoij,rhoji,drhoij,drhoji, &
                               drhoijs,drhojis,drhoij2,drhoji2,drhoij2s,drhoji2s,drhoij2m,drhoji2m, &
                               1.0_dp,1.0_dp,.true.,lmd,.true.,.false.,twopot(1,npot))
                endif
              endif
            enddo
!*******************************
!  Compute screening function  *
!*******************************
!
!  Initialise screening function to 1
!
            lnonzeroSij = .true.
            Sij = 1.0_dp
!
            if (lanyMEAMscreen) then
!
!  Set cutoffs for possible screening atoms
!
              rcut2 = rcutfactor*r2
!
!  Loop over atoms to search for those that may contribute to the screening
!
              k = 0
              npartial = 0
              do while (k.lt.numat.and.lnonzeroSij)
                k = k + 1
                neamspeck = neamfnspecptr(k)
                if (lMEAMscreen(indij,neamspeck)) then
!
!  Set basic vectors between atoms
!
                  xik0 = xclat(k) - xal
                  yik0 = yclat(k) - yal
                  zik0 = zclat(k) - zal
                  xjk0 = xik0 - xcd
                  yjk0 = yik0 - ycd
                  zjk0 = zik0 - zcd
                  xij0 = 0.5_dp*(xik0 + xjk0)
                  yij0 = 0.5_dp*(yik0 + yjk0)
                  zij0 = 0.5_dp*(zik0 + zjk0)
!
!  Compute square of distance to i-j mid point
!
                  r2ijmid = xij0*xij0 + yij0*yij0 + zij0*zij0
                  if (r2ijmid.lt.rcut2) then
!
!  Complete distances
!
                    r2ik = xik0*xik0 + yik0*yik0 + zik0*zik0
                    r2jk = xjk0*xjk0 + yjk0*yjk0 + zjk0*zjk0
!
!  Compute screening function
!
                    call meamscreen(neamspeci,neamspecj,neamspeck,r2,r2ik,r2jk,Sikj,dSikjdr,lpartial,.true.)
!
!  If screening function contribution is 0, then no need to continue for this pair
!
                    if (Sikj.eq.0.0_dp) then
                      lnonzeroSij = .false.
                      Sij = 0.0_dp
                    else
!
!  Multiply total screening product
!
                      Sij = Sij*Sikj
                      if (lpartial) then
!
!  If this atom has a screening factor between 0 and 1, we need to keep track of it since it will generate non-zero derivatives
!
                        npartial = npartial + 1
                        if (npartial.gt.partial%sa_maxdim) then
                          call changemaxsa(partial,npartial)
                        endif
                        partial%sa_atom(npartial) = k
                        partial%sa_rij(npartial) = sqrt(r2)
                        partial%sa_rik(npartial) = sqrt(r2ik)
                        partial%sa_rjk(npartial) = sqrt(r2jk)
                        partial%sa_xik(npartial) = xik0
                        partial%sa_yik(npartial) = yik0
                        partial%sa_zik(npartial) = zik0
                        partial%sa_xjk(npartial) = xjk0
                        partial%sa_yjk(npartial) = yjk0
                        partial%sa_zjk(npartial) = zjk0
                        partial%sa_Sikj(npartial) = Sikj
                        partial%sa_dSikjdr(1:3,npartial) = dSikjdr(1:3)
                      endif
                    endif
                  endif
!
!  End loop over possible screening atoms
!
                endif
              enddo
!
!  End of screening function
!
            endif
!
!  Only do remainder of work if the screening factor is non-zero
!
            if (lnonzeroSij) then
              do m = 1,npots
                npot = npotl(m)
!
!  Pair potential contribution
!
                if (nptype(npot).eq.45.or.nptype(npot).eq.55) then
                  lvalidij = .true.
                  rk = 1.0_dp/r
                  call baskes(nati,ntypi,neamspeci,natj,ntypj,neamspecj,npot,r,rk,1.0_dp,ebas,d1bas,d2bas,.true.,.false.)
                  ebas  = ebas*ofct
                  d1bas = rk*d1bas*ofct
                  deriv(1) = deriv(1) + d1bas*xcd*Sij
                  deriv(2) = deriv(2) + d1bas*ycd*Sij
                  deriv(3) = deriv(3) + d1bas*zcd*Sij
!
                  call psiscreenderv(i,j,npartial,partial,xcd,ycd,zcd,ebas,Sij,.false.)
                endif
!
!  Density contribution
!
                if (nptype(npot).eq.19) then
                  if (lanyMEAMscreen) then
                    if (npartial.gt.0) then
!
!  Compute derivative contributions of the screening function
!
                      call meamtotalrhoscreenderv(neamspeci,npartial,partial,xcd,ycd,zcd,scrho(1,i), &
                                                  rhoi,rscrhoi,rhoij,Sij,lmd)
!
                      do np = 1,npartial
                        k = partial%sa_atom(np)
                        lopk = (.not.lfreeze.or.lopf(nrelat(k)))
                        nregionk = nregionno(nsft+nrelat(k))
!
!  i-j contribution
!
                        if (lopi) then
                          xdrv(i) = xdrv(i) - partial%sa_drhototij(1,np)*ofct
                          ydrv(i) = ydrv(i) - partial%sa_drhototij(2,np)*ofct
                          zdrv(i) = zdrv(i) - partial%sa_drhototij(3,np)*ofct
                        endif
                        if (lopj) then
                          xdrv(j) = xdrv(j) + partial%sa_drhototij(1,np)*ofct
                          ydrv(j) = ydrv(j) + partial%sa_drhototij(2,np)*ofct
                          zdrv(j) = zdrv(j) + partial%sa_drhototij(3,np)*ofct
                        endif
                        if (nregioni.ne.nregionj) then
                          xregdrv(nregioni) = xregdrv(nregioni) - partial%sa_drhototij(1,np)*ofct
                          yregdrv(nregioni) = yregdrv(nregioni) - partial%sa_drhototij(2,np)*ofct
                          zregdrv(nregioni) = zregdrv(nregioni) - partial%sa_drhototij(3,np)*ofct
                          xregdrv(nregionj) = xregdrv(nregionj) + partial%sa_drhototij(1,np)*ofct
                          yregdrv(nregionj) = yregdrv(nregionj) + partial%sa_drhototij(2,np)*ofct
                          zregdrv(nregionj) = zregdrv(nregionj) + partial%sa_drhototij(3,np)*ofct
                        endif
!
!  i-k contribution
!
                        if (lopi) then
                          xdrv(i) = xdrv(i) - partial%sa_drhototik(1,np)*ofct
                          ydrv(i) = ydrv(i) - partial%sa_drhototik(2,np)*ofct
                          zdrv(i) = zdrv(i) - partial%sa_drhototik(3,np)*ofct
                        endif
                        if (lopk) then
                          xdrv(k) = xdrv(k) + partial%sa_drhototik(1,np)*ofct
                          ydrv(k) = ydrv(k) + partial%sa_drhototik(2,np)*ofct
                          zdrv(k) = zdrv(k) + partial%sa_drhototik(3,np)*ofct
                        endif
                        if (nregioni.ne.nregionk) then
                          xregdrv(nregioni) = xregdrv(nregioni) - partial%sa_drhototik(1,np)*ofct
                          yregdrv(nregioni) = yregdrv(nregioni) - partial%sa_drhototik(2,np)*ofct
                          zregdrv(nregioni) = zregdrv(nregioni) - partial%sa_drhototik(3,np)*ofct
                          xregdrv(nregionk) = xregdrv(nregionk) + partial%sa_drhototik(1,np)*ofct
                          yregdrv(nregionk) = yregdrv(nregionk) + partial%sa_drhototik(2,np)*ofct
                          zregdrv(nregionk) = zregdrv(nregionk) + partial%sa_drhototik(3,np)*ofct
                        endif
!
!  j-k contribution
!
                        if (lopj) then
                          xdrv(j) = xdrv(j) - partial%sa_drhototjk(1,np)*ofct
                          ydrv(j) = ydrv(j) - partial%sa_drhototjk(2,np)*ofct
                          zdrv(j) = zdrv(j) - partial%sa_drhototjk(3,np)*ofct
                        endif
                        if (lopk) then
                          xdrv(k) = xdrv(k) + partial%sa_drhototjk(1,np)*ofct
                          ydrv(k) = ydrv(k) + partial%sa_drhototjk(2,np)*ofct
                          zdrv(k) = zdrv(k) + partial%sa_drhototjk(3,np)*ofct
                        endif
                        if (nregionj.ne.nregionk) then
                          xregdrv(nregionj) = xregdrv(nregionj) - partial%sa_drhototjk(1,np)*ofct
                          yregdrv(nregionj) = yregdrv(nregionj) - partial%sa_drhototjk(2,np)*ofct
                          zregdrv(nregionj) = zregdrv(nregionj) - partial%sa_drhototjk(3,np)*ofct
                          xregdrv(nregionk) = xregdrv(nregionk) + partial%sa_drhototjk(1,np)*ofct
                          yregdrv(nregionk) = yregdrv(nregionk) + partial%sa_drhototjk(2,np)*ofct
                          zregdrv(nregionk) = zregdrv(nregionk) + partial%sa_drhototjk(3,np)*ofct
                        endif
                      enddo

                      call meamtotalrhoscreenderv(neamspecj,npartial,partial,xcd,ycd,zcd,scrho(1,j), &
                                                  rhoj,rscrhoj,rhoji,Sij,lmd)
!
                      do np = 1,npartial
                        k = partial%sa_atom(np)
                        lopk = (.not.lfreeze.or.lopf(nrelat(k)))
                        nregionk = nregionno(nsft+nrelat(k))
!
!  i-j contribution
!
                        if (lopi) then
                          xdrv(i) = xdrv(i) - partial%sa_drhototij(1,np)*ofct
                          ydrv(i) = ydrv(i) - partial%sa_drhototij(2,np)*ofct
                          zdrv(i) = zdrv(i) - partial%sa_drhototij(3,np)*ofct
                        endif
                        if (lopj) then
                          xdrv(j) = xdrv(j) + partial%sa_drhototij(1,np)*ofct
                          ydrv(j) = ydrv(j) + partial%sa_drhototij(2,np)*ofct
                          zdrv(j) = zdrv(j) + partial%sa_drhototij(3,np)*ofct
                        endif
                        if (nregioni.ne.nregionj) then
                          xregdrv(nregioni) = xregdrv(nregioni) - partial%sa_drhototij(1,np)*ofct
                          yregdrv(nregioni) = yregdrv(nregioni) - partial%sa_drhototij(2,np)*ofct
                          zregdrv(nregioni) = zregdrv(nregioni) - partial%sa_drhototij(3,np)*ofct
                          xregdrv(nregionj) = xregdrv(nregionj) + partial%sa_drhototij(1,np)*ofct
                          yregdrv(nregionj) = yregdrv(nregionj) + partial%sa_drhototij(2,np)*ofct
                          zregdrv(nregionj) = zregdrv(nregionj) + partial%sa_drhototij(3,np)*ofct
                        endif
!
!  i-k contribution
!
                        if (lopi) then
                          xdrv(i) = xdrv(i) - partial%sa_drhototik(1,np)*ofct
                          ydrv(i) = ydrv(i) - partial%sa_drhototik(2,np)*ofct
                          zdrv(i) = zdrv(i) - partial%sa_drhototik(3,np)*ofct
                        endif
                        if (lopk) then
                          xdrv(k) = xdrv(k) + partial%sa_drhototik(1,np)*ofct
                          ydrv(k) = ydrv(k) + partial%sa_drhototik(2,np)*ofct
                          zdrv(k) = zdrv(k) + partial%sa_drhototik(3,np)*ofct
                        endif
                        if (nregioni.ne.nregionk) then
                          xregdrv(nregioni) = xregdrv(nregioni) - partial%sa_drhototik(1,np)*ofct
                          yregdrv(nregioni) = yregdrv(nregioni) - partial%sa_drhototik(2,np)*ofct
                          zregdrv(nregioni) = zregdrv(nregioni) - partial%sa_drhototik(3,np)*ofct
                          xregdrv(nregionk) = xregdrv(nregionk) + partial%sa_drhototik(1,np)*ofct
                          yregdrv(nregionk) = yregdrv(nregionk) + partial%sa_drhototik(2,np)*ofct
                          zregdrv(nregionk) = zregdrv(nregionk) + partial%sa_drhototik(3,np)*ofct
                        endif
!
!  j-k contribution
!
                        if (lopj) then
                          xdrv(j) = xdrv(j) - partial%sa_drhototjk(1,np)*ofct
                          ydrv(j) = ydrv(j) - partial%sa_drhototjk(2,np)*ofct
                          zdrv(j) = zdrv(j) - partial%sa_drhototjk(3,np)*ofct
                        endif
                        if (lopk) then
                          xdrv(k) = xdrv(k) + partial%sa_drhototjk(1,np)*ofct
                          ydrv(k) = ydrv(k) + partial%sa_drhototjk(2,np)*ofct
                          zdrv(k) = zdrv(k) + partial%sa_drhototjk(3,np)*ofct
                        endif
                        if (nregionj.ne.nregionk) then
                          xregdrv(nregionj) = xregdrv(nregionj) - partial%sa_drhototjk(1,np)*ofct
                          yregdrv(nregionj) = yregdrv(nregionj) - partial%sa_drhototjk(2,np)*ofct
                          zregdrv(nregionj) = zregdrv(nregionj) - partial%sa_drhototjk(3,np)*ofct
                          xregdrv(nregionk) = xregdrv(nregionk) + partial%sa_drhototjk(1,np)*ofct
                          yregdrv(nregionk) = yregdrv(nregionk) + partial%sa_drhototjk(2,np)*ofct
                          zregdrv(nregionk) = zregdrv(nregionk) + partial%sa_drhototjk(3,np)*ofct
                        endif
                      enddo
                    endif
                  endif
!
!  Scale density and derivatives by screening factor
!
                  drhoij(1:3,1:maxmeamcomponent) = Sij*drhoij(1:3,1:maxmeamcomponent)
                  drhoji(1:3,1:maxmeamcomponent) = Sij*drhoji(1:3,1:maxmeamcomponent)
                  if (lmd) then
                    drhoijs(1:6,1:maxmeamcomponent) = Sij*drhoijs(1:6,1:maxmeamcomponent)
                    drhojis(1:6,1:maxmeamcomponent) = Sij*drhojis(1:6,1:maxmeamcomponent)
                  endif
!
                  call meamtotalrhoderv(neamspeci,scrho(1,i),rhoi,drhoij,drhototij,drhoijs,drhototijs, &
                                        drhoij2,drhototij2,drhoij2s,drhototij2s,drhoij2m,drhototij2m, &
                                        lmd,.false.)
                  call meamtotalrhoderv(neamspecj,scrho(1,j),rhoj,drhoji,drhototji,drhojis,drhototjis, &
                                        drhoji2,drhototji2,drhoji2s,drhototji2s,drhoji2m,drhototji2m, &
                                        lmd,.false.)
                endif
              enddo
            endif
          else
            do m = 1,npots
              npot = npotl(m)
              if (nptype(npot).eq.19) then
                if (r.gt.rpot2(npot).and.r.le.rpot(npot)) then
                  lvalidij = .true.
                  call eamrho(nati,ntypi,natj,ntypj,r,rpot(npot),xcd,ycd,zcd,rhoij,rhoji,drhototij,drhototji, &
                              drhototijs,drhototjis,drhototij2,drhototji2,drhototij2s,drhototji2s, &
                              drhototij2m,drhototji2m,drhototij3,drhototji3,1.0_dp,1.0_dp,lmd,.true.,.false.,.false., &
                              twopot(1,npot))
                endif
              endif
            enddo
          endif
!
!  Combine derivative terms
!
          if (QMMMmode(ncf).gt.0) then
            if (nregiontypi.ne.1) then
              deriv(1:3) = deriv(1:3) + rscrhoi*drhototij(1:3)*ofct
              if (lmd) then
                derivs(1:6) = derivs(1:6) + rscrhoi*drhototijs(1:6)*ofct
              endif
            endif
            if (nregiontypj.ne.1) then
              deriv(1:3) = deriv(1:3) + rscrhoj*drhototji(1:3)*ofct
              if (lmd) then
                derivs(1:6) = derivs(1:6) + rscrhoj*drhototjis(1:6)*ofct
              endif
            endif
          else
            deriv(1:3) = deriv(1:3) + (rscrhoi*drhototij(1:3) + rscrhoj*drhototji(1:3))*ofct
            if (lmd) then
              derivs(1:6) = derivs(1:6) + (rscrhoi*drhototijs(1:6) + rscrhoj*drhototjis(1:6))*ofct
            endif
          endif
        endif
        if (lvalidij) then
!******************************
!  Internal first derivatives *
!******************************
          if (lMEAM) then
            if (lopi) then
              xdrvnr(i) = xdrvnr(i) - deriv(1)
              ydrvnr(i) = ydrvnr(i) - deriv(2)
              zdrvnr(i) = zdrvnr(i) - deriv(3)
            endif
            if (lopj) then
              xdrvnr(j) = xdrvnr(j) + deriv(1)
              ydrvnr(j) = ydrvnr(j) + deriv(2)
              zdrvnr(j) = zdrvnr(j) + deriv(3)
            endif
          else
            if (lopi) then
              xdrv(i) = xdrv(i) - deriv(1)
              ydrv(i) = ydrv(i) - deriv(2)
              zdrv(i) = zdrv(i) - deriv(3)
            endif
            if (lopj) then
              xdrv(j) = xdrv(j) + deriv(1)
              ydrv(j) = ydrv(j) + deriv(2)
              zdrv(j) = zdrv(j) + deriv(3)
            endif
          endif
          if (nregioni.ne.nregionj) then
            xregdrv(nregioni) = xregdrv(nregioni) - deriv(1)
            yregdrv(nregioni) = yregdrv(nregioni) - deriv(2)
            zregdrv(nregioni) = zregdrv(nregioni) - deriv(3)
            xregdrv(nregionj) = xregdrv(nregionj) + deriv(1)
            yregdrv(nregionj) = yregdrv(nregionj) + deriv(2)
            zregdrv(nregionj) = zregdrv(nregionj) + deriv(3)
          endif
        endif
!**************************************
!  End of valid distance i-j section  *
!**************************************
      endif
    enddo jloop
  enddo iloop
!
!  End of real space part - perform general tasks
!
999 continue
!
!  Free local memory
!
  deallocate(npotl,stat=status)
  if (status/=0) call deallocate_error('manymd0','npotl')
!
!  Exit point
!
1000 continue
!
!  Unscale density
!
  call eamscalescrho(-1_i4)
!
!  Timing
!
  time2 = g_cpu_time()
  tmany = tmany + time2 - time1
#ifdef TRACE
  call trace_out('manymd0')
#endif
!
  return
  end
