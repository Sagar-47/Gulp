  subroutine many0dd(emany,lgrad1,lgrad2)
!
!  Subroutine for calculating the many-body energy from the
!  Sutton-Chen potential(s). Requires real space routine to
!  be called first to evaluate the repulsive contribution
!  and the rho parameters.
!
!  Finite clusters only version
!  Distributed memory parallel version.
!
!  On entry the array scrho must contain the density at each atomic site.
!
!  lfreeze = if .true. only do derivatives for atoms with a degree of freedom
!
!   1/17 Created from many0d
!   7/17 I/O handling modified for parallel case
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
!  Julian Gale, CIC, Curtin University, Feb 2018
!
  use configurations, only : nregionno, nregiontype, QMMMmode
  use control
  use current
  use derivatives
  use eam
  use general,        only : smallself
  use iochannels,     only : ioout
  use mdlogic
  use optimisation
  use parallel
  use progress,       only : lduring_opt
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
  logical                                      :: lgrad1
  logical                                      :: lgrad2
  real(dp)                                     :: emany
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: iloc
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: ixc
  integer(i4)                                  :: iyc
  integer(i4)                                  :: izc
  integer(i4)                                  :: ixf
  integer(i4)                                  :: iyf
  integer(i4)                                  :: izf
  integer(i4)                                  :: j
  integer(i4)                                  :: jx
  integer(i4)                                  :: jy
  integer(i4)                                  :: jz
  integer(i4)                                  :: jxc
  integer(i4)                                  :: jyc
  integer(i4)                                  :: jzc
  integer(i4)                                  :: k
  integer(i4)                                  :: m
  integer(i4)                                  :: n
  integer(i4)                                  :: nat1
  integer(i4)                                  :: nat2
  integer(i4)                                  :: nati
  integer(i4)                                  :: natj
  integer(i4)                                  :: natk
  integer(i4)                                  :: neamspeci
  integer(i4)                                  :: neamspecj
  integer(i4)                                  :: neamspeck
  integer(i4)                                  :: nff
  integer(i4)                                  :: npot
  integer(i4)                                  :: npotik
  integer(i4)                                  :: npotjk
  integer(i4), dimension(:), allocatable       :: npotikptr
  integer(i4), dimension(:), allocatable       :: npotjkptr
  integer(i4), dimension(:), allocatable       :: npotl
  integer(i4)                                  :: npots
  integer(i4)                                  :: nregioni
  integer(i4)                                  :: nregionj
  integer(i4)                                  :: nregiontypi
  integer(i4)                                  :: nregiontypj
  integer(i4)                                  :: ntyp1
  integer(i4)                                  :: ntyp2
  integer(i4)                                  :: ntypi
  integer(i4)                                  :: ntypj
  integer(i4)                                  :: ntypk
  integer(i4)                                  :: status
  logical                                      :: lanyvalidik
  logical                                      :: lanyvalidjk
  logical                                      :: lopj
  logical                                      :: lQMMMok
  logical                                      :: lvalidij
  logical                                      :: lvalidik
  logical                                      :: lvalidjk
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: cut2
  real(dp)                                     :: cut2Xscale
  real(dp)                                     :: cut2k
  real(dp)                                     :: cut2p
  real(dp)                                     :: cut2r
  real(dp)                                     :: cut2rk
  real(dp)                                     :: deriv(3)
  real(dp)                                     :: derivs(6)
  real(dp)                                     :: deriv2(6)
  real(dp)                                     :: drhoij(3,maxmeamcomponent)
  real(dp)                                     :: drhoijs(6,maxmeamcomponent)
  real(dp)                                     :: drhoij2(6,maxmeamcomponent)
  real(dp)                                     :: drhoij2s(21,maxmeamcomponent)
  real(dp)                                     :: drhoij2m(6,3,maxmeamcomponent)
  real(dp)                                     :: drhoik(3,maxmeamcomponent)
  real(dp)                                     :: drhoiks(6,maxmeamcomponent)
  real(dp)                                     :: drhoik2(6,maxmeamcomponent)
  real(dp)                                     :: drhoik2s(21,maxmeamcomponent)
  real(dp)                                     :: drhoik2m(6,3,maxmeamcomponent)
  real(dp)                                     :: drhoji(3,maxmeamcomponent)
  real(dp)                                     :: drhojis(6,maxmeamcomponent)
  real(dp)                                     :: drhoji2(6,maxmeamcomponent)
  real(dp)                                     :: drhoji2s(21,maxmeamcomponent)
  real(dp)                                     :: drhoji2m(6,3,maxmeamcomponent)
  real(dp)                                     :: drhojk(3,maxmeamcomponent)
  real(dp)                                     :: drhojks(6,maxmeamcomponent)
  real(dp)                                     :: drhojk2(6,maxmeamcomponent)
  real(dp)                                     :: drhojk2s(21,maxmeamcomponent)
  real(dp)                                     :: drhojk2m(6,3,maxmeamcomponent)
  real(dp)                                     :: drhoki(3,maxmeamcomponent)
  real(dp)                                     :: drhokis(6,maxmeamcomponent)
  real(dp)                                     :: drhoki2(6,maxmeamcomponent)
  real(dp)                                     :: drhoki2s(21,maxmeamcomponent)
  real(dp)                                     :: drhoki2m(6,3,maxmeamcomponent)
  real(dp)                                     :: drhokj(3,maxmeamcomponent)
  real(dp)                                     :: drhokjs(6,maxmeamcomponent)
  real(dp)                                     :: drhokj2(6,maxmeamcomponent)
  real(dp)                                     :: drhokj2s(21,maxmeamcomponent)
  real(dp)                                     :: drhokj2m(6,3,maxmeamcomponent)
  real(dp)                                     :: drhototij(3)
  real(dp)                                     :: drhototijs(6)
  real(dp)                                     :: drhototij2(6)
  real(dp)                                     :: drhototij2s(21)
  real(dp)                                     :: drhototij2m(6,3)
  real(dp)                                     :: drhototij3(10)
  real(dp)                                     :: drhototik(3)
  real(dp)                                     :: drhototiks(6)
  real(dp)                                     :: drhototik2(6)
  real(dp)                                     :: drhototik2s(21)
  real(dp)                                     :: drhototik2m(6,3)
  real(dp)                                     :: drhototik3(10)
  real(dp)                                     :: drhototji(3)
  real(dp)                                     :: drhototjis(6)
  real(dp)                                     :: drhototji2(6)
  real(dp)                                     :: drhototji2s(21)
  real(dp)                                     :: drhototji2m(6,3)
  real(dp)                                     :: drhototji3(10)
  real(dp)                                     :: drhototjk(3)
  real(dp)                                     :: drhototjks(6)
  real(dp)                                     :: drhototjk2(6)
  real(dp)                                     :: drhototjk2s(21)
  real(dp)                                     :: drhototjk2m(6,3)
  real(dp)                                     :: drhototjk3(10)
  real(dp)                                     :: drhototki(3)
  real(dp)                                     :: drhototkis(6)
  real(dp)                                     :: drhototki2(6)
  real(dp)                                     :: drhototki2s(21)
  real(dp)                                     :: drhototki2m(6,3)
  real(dp)                                     :: drhototki3(10)
  real(dp)                                     :: drhototkj(3)
  real(dp)                                     :: drhototkjs(6)
  real(dp)                                     :: drhototkj2(6)
  real(dp)                                     :: drhototkj2s(21)
  real(dp)                                     :: drhototkj2m(6,3)
  real(dp)                                     :: drhototkj3(10)
  real(dp)                                     :: drhototijk2(3,3)
  real(dp)                                     :: drhototjik2(3,3)
  real(dp)                                     :: drhototkij2(3,3)
  real(dp)                                     :: drhototijk2s(6,6)
  real(dp)                                     :: drhototjik2s(6,6)
  real(dp)                                     :: drhototkij2s(6,6)
  real(dp)                                     :: drhototijk2m(6,3)
  real(dp)                                     :: drhototjik2m(6,3)
  real(dp)                                     :: drhototkij2m(6,3)
  real(dp)                                     :: dt1
  real(dp)                                     :: dt2
  real(dp)                                     :: eeam
  real(dp)                                     :: emanytrm
  real(dp)                                     :: oci
  real(dp)                                     :: ocj
  real(dp)                                     :: ock
  real(dp)                                     :: ofct
  real(dp)                                     :: ofctijk
  real(dp)                                     :: r
  real(dp)                                     :: r2
  real(dp)                                     :: rhoi
  real(dp)                                     :: rhoj
  real(dp)                                     :: rhok
  real(dp)                                     :: rhoij(maxmeamcomponent)
  real(dp)                                     :: rhoji(maxmeamcomponent)
  real(dp)                                     :: rhoik(maxmeamcomponent)
  real(dp)                                     :: rhoki(maxmeamcomponent)
  real(dp)                                     :: rhojk(maxmeamcomponent)
  real(dp)                                     :: rhokj(maxmeamcomponent)
  real(dp)                                     :: rik
  real(dp)                                     :: rik2
  real(dp)                                     :: rjk
  real(dp)                                     :: rjk2
  real(dp)                                     :: rp
  real(dp)                                     :: rpik
  real(dp)                                     :: rpjk
  real(dp)                                     :: rpijk
  real(dp)                                     :: rscrhoi
  real(dp)                                     :: rscrhoi3
  real(dp)                                     :: rscrhoi5
  real(dp)                                     :: rscrhoj
  real(dp)                                     :: rscrhoj3
  real(dp)                                     :: rscrhoj5
  real(dp)                                     :: rscrhok
  real(dp)                                     :: rscrhok3
  real(dp)                                     :: rscrhok5
  real(dp)                                     :: scmax
  real(dp)                                     :: time1
  real(dp)                                     :: time2
  real(dp)                                     :: xal
  real(dp)                                     :: yal
  real(dp)                                     :: zal
  real(dp)                                     :: xcd
  real(dp)                                     :: ycd
  real(dp)                                     :: zcd
  real(dp)                                     :: xcd1
  real(dp)                                     :: ycd1
  real(dp)                                     :: zcd1
  real(dp)                                     :: xcd2
  real(dp)                                     :: ycd2
  real(dp)                                     :: zcd2
#ifdef TRACE
  call trace_in('many0dd')
#endif
!
  time1 = g_cpu_time()
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
  do iloc = 1,natomsonnode
    i = node2atom(iloc)
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
      if (lPrintEAM.and.ioproc) then
        write(ioout,'(7x,i8,4x,f24.10,4x,f24.12)') i,rhoi,emanytrm
      endif
    endif
  enddo
  if (lPrintEAM.and.ioproc) then
!
!  Closing banner for energy decomposition
!
    write(ioout,'(''--------------------------------------------------------------------------------'')')
  endif
!
!  If no forces are needed then we don't need to do loops over atoms, so just return
!
  if (.not.lgrad1) goto 1000
!
  if (lduring_opt) then
    cut2Xscale = eamXcutfactor**2
  else
    cut2Xscale = 1.0_dp
  endif
!
!  From here on we can assume that lgrad1  =  .true.
!
!  Find maximum cut - off radius
!
  scmax = 0.0_dp
  do i = 1,npote
    if (nptype(i).eq.19) then
      if (rpot(i).gt.scmax) scmax = rpot(i)
    endif
  enddo
!
!  Find number of unfrozen atoms
!
  if (lfreeze) then
    nff = 0
    do i = 1,numat
      if (lopf(i)) nff = nff + 1
    enddo
  else
    nff = numat
  endif
!
!  Allocate local memory
!
  allocate(npotl(npote),stat=status)
  if (status/=0) call outofmemory('many0dd','npotl')
  if (lgrad2) then
    allocate(npotikptr(npote),stat=status)
    if (status/=0) call outofmemory('many0dd','npotikptr')
    allocate(npotjkptr(npote),stat=status)
    if (status/=0) call outofmemory('many0dd','npotjkptr')
  endif
!
!  Set up cutoffs
!
  cut2p = cutp*cutp
  if (lnoreal) goto 999
!
!  Outer loop over sites
!
  ixc = - 2
  iyc = - 1
  izc =   0
  iloop: do iloc = 1,noptatloc
    i = noptatlocptr(iloc)
!
!  Handle ixc - izc before rho check otherwise values go wrong
!
    ixc = ixc + 3
    iyc = iyc + 3
    izc = izc + 3
!
    ixf = 3*(noptatrptr(i) - 1) + 1
    iyf = ixf + 1
    izf = ixf + 2
!     
!  Find EAM species for i
!  
    neamspeci = neamfnspecptr(i)
!
!  Evaluate functional derivatives
!
    if (lMEAMfn) then
      call meamfnderv(neamfn,neamspeci,scrho(1,i),rhoi,eeam,rscrhoi,rscrhoi3,rscrhoi5,.true.,lgrad2,.false.)
    else
      call eamfnderv(neamfn,neamspeci,scrho(1,i),eeam,rscrhoi,rscrhoi3,rscrhoi5,.true.,lgrad2,.false.)
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
    jxc = - 2
    jyc = - 1
    jzc =   0
    jloop: do j = 1,numat
!
!  Freezing flag  -  need to do this before checking rho values otherwise jxc - jzc go wrong
!
      lopj = (.not.lfreeze.or.lopf(j))
      if (lopj) then
        jxc = jxc + 3
        jyc = jyc + 3
        jzc = jzc + 3
      endif
!
!  Exclude case of i = j
!
      if (i.eq.j) cycle jloop
!     
!  Find EAM species for j
!  
      neamspecj = neamfnspecptr(j)
!
!  Evaluate functional derivatives
!
      if (lMEAMfn) then
        call meamfnderv(neamfn,neamspecj,scrho(1,j),rhoj,eeam,rscrhoj,rscrhoj3,rscrhoj5,.true.,lgrad2,.false.)
      else
        call eamfnderv(neamfn,neamspecj,scrho(1,j),eeam,rscrhoj,rscrhoj3,rscrhoj5,.true.,lgrad2,.false.)
        rhoj = scrho(1,j)
      endif
!
!  If there is no density at either centre and this is not a second derivative run then cycle
!
      if (rhoi.eq.0.0_dp.and.rhoj.eq.0.0_dp.and..not.lgrad2) cycle jloop
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
!
!  Set flags such that if one atom is not being optimised place 3 x 3 
!  second derivative matrix in the on - diagonal block
!
      ix = ixc
      iy = iyc
      iz = izc
      if (lopj) then
        jx = jxc
        jy = jyc
        jz = jzc
      else
        jx = ixf
        jy = iyf
        jz = izf
      endif
      ofct = oci*ocj
!
!  Locate potential number
!  Check whether potential requires specific types
!  Calculate sum of all dispersion terms for pair
!
      npots = 0
      rp = 0.0_dp
      do n = 1,npote
        if (nptype(n).eq.19) then
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
      if (npots.eq.0.and..not.lgrad2) cycle jloop
      if (lgrad2) then
!
!  Need to make cut - off equal to double the maximum to ensure all triangles are included
!
        rp = 2.0_dp*scmax
        cut2r = rp*rp
        if (cut2r.gt.4.0_dp*cut2p) cut2r = cut2p
      else
        cut2r = rp*rp
        if (cut2r.gt.cut2p) cut2r = cut2p
      endif
      cut2 = cut2r
      rp = sqrt(cut2)
      r2 = xcd*xcd + ycd*ycd + zcd*zcd
      if (r2.gt.smallself.and.r2.le.cut2) then
!***************************************************
!  Calculate many-body contribution in real space  *
!***************************************************
        deriv(1:3) = 0.0_dp
        deriv2(1:6) = 0.0_dp
        if (lmd) then
          derivs(1:6) = 0.0_dp
        endif
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
          if (lgrad2) then
            drhoij2(1:6,1:maxmeamcomponent) = 0.0_dp
            drhoji2(1:6,1:maxmeamcomponent) = 0.0_dp
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
        if (lgrad2) then
          drhototij2(1:6) = 0.0_dp
          drhototji2(1:6) = 0.0_dp
        endif
        lvalidij = .false.
        if (npots.gt.0) then
          do m = 1,npots
            npot = npotl(m)
            if (r.gt.rpot2(npot).and.r.le.rpot(npot).and.r.le.rpmax) then
              lvalidij = .true.
!
!  Calculate density derivatives
!
              if (lMEAMden) then
                call meamrho(nati,ntypi,natj,ntypj,r,rpot(npot),xcd,ycd,zcd,rhoij,rhoji,drhoij,drhoji, &
                             drhoijs,drhojis,drhoij2,drhoji2,drhoij2s,drhoji2s,drhoij2m,drhoji2m, &
                             1.0_dp,1.0_dp,.true.,lmd,.true.,lgrad2,twopot(1,npot))
                call meamtotalrhoderv(neamspeci,scrho(1,i),rhoi,drhoij,drhototij,drhoijs,drhototijs, &
                                      drhoij2,drhototij2,drhoij2s,drhototij2s,drhoij2m,drhototij2m, &
                                      lmd,lgrad2)
                call meamtotalrhoderv(neamspecj,scrho(1,j),rhoj,drhoji,drhototji,drhojis,drhototjis, &
                                      drhoji2,drhototji2,drhoji2s,drhototji2s,drhoji2m,drhototji2m, &
                                      lmd,lgrad2)
              else
                call eamrho(nati,ntypi,natj,ntypj,r,rpot(npot),xcd,ycd,zcd,rhoij,rhoji,drhototij,drhototji, &
                            drhototijs,drhototjis,drhototij2,drhototji2,drhototij2s,drhototji2s, &
                            drhototij2m,drhototji2m,drhototij3,drhototji3,1.0_dp,1.0_dp,lmd,.true.,lgrad2,.false., &
                            twopot(1,npot))
              endif
            endif
          enddo
!
!  Combine derivative terms
!
          if (QMMMmode(ncf).gt.0) then
            if (nregiontypi.ne.1) then
              deriv(1:3) = deriv(1:3) + rscrhoi*drhototij(1:3)*ofct
            endif
            if (nregiontypj.ne.1) then
              deriv(1:3) = deriv(1:3) + rscrhoj*drhototji(1:3)*ofct
            endif
          else
            deriv(1:3) = deriv(1:3) + (rscrhoi*drhototij(1:3) + rscrhoj*drhototji(1:3))*ofct
          endif
          if (lmd) then
            derivs(1:6) = derivs(1:6) + (rscrhoi*drhototijs(1:6) + rscrhoj*drhototjis(1:6))*ofct
          endif
          if (lgrad2) then
            deriv2(1:6) = deriv2(1:6) + (rscrhoi*drhototij2(1:6) + rscrhoj*drhototji2(1:6))*ofct
            deriv2(1) = deriv2(1) + ocj*rscrhoi3*drhototij(1)*drhototij(1)*ofct + oci*rscrhoj3*drhototji(1)*drhototji(1)*ofct
            deriv2(2) = deriv2(2) + ocj*rscrhoi3*drhototij(1)*drhototij(2)*ofct + oci*rscrhoj3*drhototji(1)*drhototji(2)*ofct
            deriv2(3) = deriv2(3) + ocj*rscrhoi3*drhototij(1)*drhototij(3)*ofct + oci*rscrhoj3*drhototji(1)*drhototji(3)*ofct
            deriv2(4) = deriv2(4) + ocj*rscrhoi3*drhototij(2)*drhototij(2)*ofct + oci*rscrhoj3*drhototji(2)*drhototji(2)*ofct
            deriv2(5) = deriv2(5) + ocj*rscrhoi3*drhototij(2)*drhototij(3)*ofct + oci*rscrhoj3*drhototji(2)*drhototji(3)*ofct
            deriv2(6) = deriv2(6) + ocj*rscrhoi3*drhototij(3)*drhototij(3)*ofct + oci*rscrhoj3*drhototji(3)*drhototji(3)*ofct
          endif
!
        endif
        if (lvalidij) then
!******************************
!  Internal first derivatives *
!******************************
          if (lMEAM) then
            xdrvnr(i) = xdrvnr(i) - deriv(1)
            ydrvnr(i) = ydrvnr(i) - deriv(2)
            zdrvnr(i) = zdrvnr(i) - deriv(3)
          else
            xdrv(i) = xdrv(i) - deriv(1)
            ydrv(i) = ydrv(i) - deriv(2)
            zdrv(i) = zdrv(i) - deriv(3)
          endif
          if (nregioni.ne.nregionj) then
            xregdrv(nregioni) = xregdrv(nregioni) - deriv(1)
            yregdrv(nregioni) = yregdrv(nregioni) - deriv(2)
            zregdrv(nregioni) = zregdrv(nregioni) - deriv(3)
          endif
!********************************
!  Internal second derivatives  *
!********************************
          if (lgrad2) then
            derv2(jx,ix) = derv2(jx,ix) - deriv2(1)
            derv2(jy,ix) = derv2(jy,ix) - deriv2(2)
            derv2(jz,ix) = derv2(jz,ix) - deriv2(3)
            derv2(jx,iy) = derv2(jx,iy) - deriv2(2)
            derv2(jy,iy) = derv2(jy,iy) - deriv2(4)
            derv2(jz,iy) = derv2(jz,iy) - deriv2(5)
            derv2(jx,iz) = derv2(jx,iz) - deriv2(3)
            derv2(jy,iz) = derv2(jy,iz) - deriv2(5)
            derv2(jz,iz) = derv2(jz,iz) - deriv2(6)
          endif
        endif
        if (lgrad2) then
!******************************************************************
!  Start of third atom loop  -  only needed for second derivatives  *
!******************************************************************
          kloop: do k = 1,numat
!
!  If k  =  i or j then skip
!
            if (k.eq.i.or.k.eq.j) cycle kloop
!
            natk = nat(k)
            ntypk = nftype(k)
            xcd1 = xclat(k) - xal
            ycd1 = yclat(k) - yal
            zcd1 = zclat(k) - zal
            xcd2 = xcd1 - xcd
            ycd2 = ycd1 - ycd
            zcd2 = zcd1 - zcd
            ock = occuf(k)
            ofctijk = ofct*ock
!
!  Check whether there are any potentials between i - k or j - k
!
            npotik = 0
            npotjk = 0
            rpik = 0.0_dp
            rpjk = 0.0_dp
            do n = 1,npote
              if (nptype(n).eq.19) then
                lvalidik = .false.
                if (nati.eq.nspec1(n).and.natk.eq.nspec2(n)) then
                  if (ntypi.eq.nptyp1(n).or.nptyp1(n).eq.0) then
                    if (ntypk.eq.nptyp2(n).or.nptyp2(n).eq.0) lvalidik = .true.
                  endif
                elseif (nati.eq.nspec2(n).and.natk.eq.nspec1(n)) then
                  if (ntypi.eq.nptyp2(n).or.nptyp2(n).eq.0) then
                    if (ntypk.eq.nptyp1(n).or.nptyp1(n).eq.0) lvalidik = .true.
                  endif
                endif
                if (lvalidik) then
                  npotik = npotik + 1
                  npotikptr(npotik) = n
                  if (rpot(n).gt.rpik) rpik = rpot(n)
                endif
!
                lvalidjk = .false.
                if (natj.eq.nspec1(n).and.natk.eq.nspec2(n)) then
                  if (ntypj.eq.nptyp1(n).or.nptyp1(n).eq.0) then
                    if (ntypk.eq.nptyp2(n).or.nptyp2(n).eq.0) lvalidjk = .true.
                  endif
                elseif (natj.eq.nspec2(n).and.natk.eq.nspec1(n)) then
                  if (ntypj.eq.nptyp2(n).or.nptyp2(n).eq.0) then
                    if (ntypk.eq.nptyp1(n).or.nptyp1(n).eq.0) lvalidjk = .true.
                  endif
                endif
                if (lvalidjk) then
                  npotjk = npotjk + 1
                  npotjkptr(npotjk) = n
                  if (rpot(n).gt.rpjk) rpjk = rpot(n)
                endif
              endif
            enddo
            rpijk = max(rpik,rpjk)
!
!  If no valid potentials for i - k or j - k then skip
!
            if ((npotik+npotjk).eq.0) cycle kloop
            cut2rk = rpijk*rpijk
            if (cut2rk.gt.cut2p) cut2rk = cut2p
            cut2k = cut2Xscale*cut2rk
!     
!  Find EAM species for k
!  
            neamspeck = neamfnspecptr(k)
!
!  Evaluate functional derivatives
!
            if (lMEAMfn) then
              call meamfnderv(neamfn,neamspeck,scrho(1,k),rhok,eeam,rscrhok,rscrhok3,rscrhok5,.true.,lgrad2,.false.)
            else
              call eamfnderv(neamfn,neamspeck,scrho(1,k),eeam,rscrhok,rscrhok3,rscrhok5,.true.,lgrad2,.false.)
              rhok = scrho(1,k)
            endif
!
!  If no rho then skip
!
            if (rhoi.eq.0.0_dp.and.rhoj.eq.0.0_dp.and.rhok.eq.0.0_dp) cycle kloop
!
!  Set up constants for k
!
            rik2 = xcd1*xcd1 + ycd1*ycd1 + zcd1*zcd1
            rjk2 = xcd2*xcd2 + ycd2*ycd2 + zcd2*zcd2
!************************************************************
!  Calculate triangular contribution to second derivatives  *
!************************************************************
            lanyvalidik = .false.
            lanyvalidjk = .false.
            if (rik2.gt.smallself.and.rik2.le.cut2k.and.npotik.gt.0) then
!***********************
!  i - k contribution  *
!***********************
!
!  Zero terms
!
              if (lMEAM) then
                rhoik(1:maxmeamcomponent) = 0.0_dp
                rhoki(1:maxmeamcomponent) = 0.0_dp
                drhoik(1:3,1:maxmeamcomponent) = 0.0_dp
                drhoki(1:3,1:maxmeamcomponent) = 0.0_dp
              else
                rhoik(1) = 0.0_dp
                rhoki(1) = 0.0_dp
              endif
              drhototik(1:3) = 0.0_dp
              drhototki(1:3) = 0.0_dp
!
              rik = sqrt(rik2)
!
!  Loop over potentials to find many - body ones
!
              do n = 1,npotik
                m = npotikptr(n)
                lvalidik = .false.
                if (rik.gt.rpot2(m).and.rik.le.rpot(m)) then
                  lvalidik = .true.
                  if (lvalidik) then
!
!  Calculate density derivatives
!
                    if (lMEAMden) then
                      call meamrho(nati,ntypi,natk,ntypk,rik,rpot(m),xcd1,ycd1,zcd1,rhoik,rhoki,drhoik,drhoki, &
                                   drhoiks,drhokis,drhoik2,drhoki2,drhoik2s,drhoki2s,drhoik2m,drhoki2m, &
                                   1.0_dp,1.0_dp,.true.,.false.,.true.,.false.,twopot(1,m))
                      call meamtotalrhoderv(neamspeci,scrho(1,i),rhoi,drhoik,drhototik,drhoiks,drhototiks, &
                                            drhoik2,drhototik2,drhoik2s,drhototik2s,drhoik2m,drhototik2m, &
                                            .false.,.false.)
                      call meamtotalrhoderv(neamspeck,scrho(1,k),rhok,drhoki,drhototki,drhokis,drhototkis, &
                                            drhoki2,drhototki2,drhoki2s,drhototki2s,drhoki2m,drhototki2m, &
                                            .false.,.false.)
                      call meamtotalrhocrossderv(neamspeci,scrho(1,i),rhoi,drhoij,drhototij,drhoik,drhototik, &
                                                 drhoijs,drhoiks,drhototijk2,drhototijk2s,drhototijk2m,.false.)
                    else
                      call eamrho(nati,ntypi,natk,ntypk,rik,rpot(m),xcd1,ycd1,zcd1,rhoik,rhoki,drhototik,drhototki, &
                                  drhototiks,drhototkis,drhototik2,drhototki2,drhototik2s,drhototki2s, &
                                  drhototik2m,drhototki2m,drhototik3,drhototki3,1.0_dp,1.0_dp,.false.,.true.,.false.,.false., &
                                  twopot(1,m))
                    endif
                    lanyvalidik = .true.
                  endif
                endif
              enddo
            endif
            if (rjk2.gt.smallself.and.rjk2.le.cut2.and.npotjk.gt.0) then
!***********************
!  j - k contribution  *
!***********************
!
!  Zero terms
!             
              if (lMEAM) then
                rhojk(1:maxmeamcomponent) = 0.0_dp
                rhokj(1:maxmeamcomponent) = 0.0_dp
                drhojk(1:3,1:maxmeamcomponent) = 0.0_dp
                drhokj(1:3,1:maxmeamcomponent) = 0.0_dp
              else
                rhojk(1) = 0.0_dp
                rhokj(1) = 0.0_dp
              endif
              drhototjk(1:3) = 0.0_dp
              drhototkj(1:3) = 0.0_dp
!
              rjk = sqrt(rjk2)
!
!  Loop over potentials to find many - body ones
!
              do n = 1,npotjk
                m = npotjkptr(n)
                lvalidjk = .false.
                if (rjk.gt.rpot2(m).and.rjk.le.rpot(m)) then
                  lvalidjk = .true.
                  if (lvalidjk) then
!
!  Calculate density derivatives
!
                    if (lMEAMden) then
                      call meamrho(natj,ntypj,natk,ntypk,rjk,rpot(m),xcd2,ycd2,zcd2,rhojk,rhokj,drhojk,drhokj, &
                                   drhojks,drhokjs,drhojk2,drhokj2,drhojk2s,drhokj2s,drhojk2m,drhokj2m, &
                                   1.0_dp,1.0_dp,.true.,.false.,.true.,.false.,twopot(1,m))
                      call meamtotalrhoderv(neamspecj,scrho(1,j),rhoj,drhojk,drhototjk,drhojks,drhototjks, &
                                            drhojk2,drhototjk2,drhojk2s,drhototjk2s,drhojk2m,drhototjk2m, &
                                            .false.,.false.)
                      call meamtotalrhoderv(neamspeck,scrho(1,k),rhok,drhokj,drhototkj,drhokjs,drhototkjs, &
                                            drhokj2,drhototkj2,drhokj2s,drhototkj2s,drhokj2m,drhototkj2m, &
                                            .false.,.false.)
                      call meamtotalrhocrossderv(neamspecj,scrho(1,j),rhoj,drhoji,drhototji,drhojk,drhototjk, &
                                                 drhojis,drhojks,drhototjik2,drhototjik2s,drhototjik2m,.false.)
                    else
                      call eamrho(natj,ntypj,natk,ntypk,rjk,rpot(m),xcd2,ycd2,zcd2,rhojk,rhokj,drhototjk,drhototkj, &
                                  drhototjks,drhototkjs,drhototjk2,drhototkj2,drhototjk2s,drhototkj2s, &
                                  drhototjk2m,drhototkj2m,drhototjk3,drhototkj3,1.0_dp,1.0_dp,.false.,.true.,.false.,.false., &
                                  twopot(1,m))
                    endif
                    lanyvalidjk = .true.
                  endif
                endif
              enddo
            endif
!
!  Cross term derivative for k-i / k-j
!
            if (lanyvalidik.and.lanyvalidjk.and.lMEAMden) then
              call meamtotalrhocrossderv(neamspeck,scrho(1,k),rhok,drhoki,drhototki,drhokj,drhototkj,drhokis,drhokjs, &
                                         drhototkij2,drhototkij2s,drhototkij2m,.false.)
            endif
!*******************************************************
!  Calculate second derivatives for i - k/j - k terms  *
!*******************************************************
!
!  i - k
!
            if (lanyvalidik) then
              if (lMEAM) then
                dt1 = rscrhoi3*ofctijk
                dt2 = rscrhoi*ofctijk
                derv2(jx,ix) = derv2(jx,ix) - dt1*drhototik(1)*drhototij(1) - dt2*drhototijk2(1,1)
                derv2(jy,ix) = derv2(jy,ix) - dt1*drhototik(1)*drhototij(2) - dt2*drhototijk2(2,1)
                derv2(jz,ix) = derv2(jz,ix) - dt1*drhototik(1)*drhototij(3) - dt2*drhototijk2(3,1)
                derv2(jx,iy) = derv2(jx,iy) - dt1*drhototik(2)*drhototij(1) - dt2*drhototijk2(1,2)
                derv2(jy,iy) = derv2(jy,iy) - dt1*drhototik(2)*drhototij(2) - dt2*drhototijk2(2,2)
                derv2(jz,iy) = derv2(jz,iy) - dt1*drhototik(2)*drhototij(3) - dt2*drhototijk2(3,2)
                derv2(jx,iz) = derv2(jx,iz) - dt1*drhototik(3)*drhototij(1) - dt2*drhototijk2(1,3)
                derv2(jy,iz) = derv2(jy,iz) - dt1*drhototik(3)*drhototij(2) - dt2*drhototijk2(2,3)
                derv2(jz,iz) = derv2(jz,iz) - dt1*drhototik(3)*drhototij(3) - dt2*drhototijk2(3,3)
              else
                dt1 = rscrhoi3*ofctijk
                derv2(jx,ix) = derv2(jx,ix) - dt1*drhototik(1)*drhototij(1)
                derv2(jy,ix) = derv2(jy,ix) - dt1*drhototik(1)*drhototij(2)
                derv2(jz,ix) = derv2(jz,ix) - dt1*drhototik(1)*drhototij(3)
                derv2(jx,iy) = derv2(jx,iy) - dt1*drhototik(2)*drhototij(1)
                derv2(jy,iy) = derv2(jy,iy) - dt1*drhototik(2)*drhototij(2)
                derv2(jz,iy) = derv2(jz,iy) - dt1*drhototik(2)*drhototij(3)
                derv2(jx,iz) = derv2(jx,iz) - dt1*drhototik(3)*drhototij(1)
                derv2(jy,iz) = derv2(jy,iz) - dt1*drhototik(3)*drhototij(2)
                derv2(jz,iz) = derv2(jz,iz) - dt1*drhototik(3)*drhototij(3)
              endif
            endif
!
!  j - k
!
            if (lanyvalidjk) then
              if (lMEAM) then
                dt1 = rscrhoj3*ofctijk
                dt2 = rscrhoj*ofctijk
                derv2(jx,ix) = derv2(jx,ix) + dt1*drhototjk(1)*drhototji(1) + dt2*drhototjik2(1,1)
                derv2(jy,ix) = derv2(jy,ix) + dt1*drhototjk(2)*drhototji(1) + dt2*drhototjik2(1,2)
                derv2(jz,ix) = derv2(jz,ix) + dt1*drhototjk(3)*drhototji(1) + dt2*drhototjik2(1,3)
                derv2(jx,iy) = derv2(jx,iy) + dt1*drhototjk(1)*drhototji(2) + dt2*drhototjik2(2,1)
                derv2(jy,iy) = derv2(jy,iy) + dt1*drhototjk(2)*drhototji(2) + dt2*drhototjik2(2,2)
                derv2(jz,iy) = derv2(jz,iy) + dt1*drhototjk(3)*drhototji(2) + dt2*drhototjik2(2,3)
                derv2(jx,iz) = derv2(jx,iz) + dt1*drhototjk(1)*drhototji(3) + dt2*drhototjik2(3,1)
                derv2(jy,iz) = derv2(jy,iz) + dt1*drhototjk(2)*drhototji(3) + dt2*drhototjik2(3,2)
                derv2(jz,iz) = derv2(jz,iz) + dt1*drhototjk(3)*drhototji(3) + dt2*drhototjik2(3,3)
              else
                dt1 = rscrhoj3*ofctijk
                derv2(jx,ix) = derv2(jx,ix) + dt1*drhototjk(1)*drhototji(1)
                derv2(jy,ix) = derv2(jy,ix) + dt1*drhototjk(2)*drhototji(1)
                derv2(jz,ix) = derv2(jz,ix) + dt1*drhototjk(3)*drhototji(1)
                derv2(jx,iy) = derv2(jx,iy) + dt1*drhototjk(1)*drhototji(2)
                derv2(jy,iy) = derv2(jy,iy) + dt1*drhototjk(2)*drhototji(2)
                derv2(jz,iy) = derv2(jz,iy) + dt1*drhototjk(3)*drhototji(2)
                derv2(jx,iz) = derv2(jx,iz) + dt1*drhototjk(1)*drhototji(3)
                derv2(jy,iz) = derv2(jy,iz) + dt1*drhototjk(2)*drhototji(3)
                derv2(jz,iz) = derv2(jz,iz) + dt1*drhototjk(3)*drhototji(3)
              endif
            endif
!
!  i - k/j - k
!
            if (lanyvalidik.and.lanyvalidjk) then
              if (lMEAM) then
                dt1 = rscrhok3*ofctijk
                dt2 = rscrhok*ofctijk
                derv2(jx,ix) = derv2(jx,ix) + dt1*drhototki(1)*drhototkj(1) + dt2*drhototkij2(1,1)
                derv2(jy,ix) = derv2(jy,ix) + dt1*drhototki(1)*drhototkj(2) + dt2*drhototkij2(1,2)
                derv2(jz,ix) = derv2(jz,ix) + dt1*drhototki(1)*drhototkj(3) + dt2*drhototkij2(1,3)
                derv2(jx,iy) = derv2(jx,iy) + dt1*drhototki(2)*drhototkj(1) + dt2*drhototkij2(2,1)
                derv2(jy,iy) = derv2(jy,iy) + dt1*drhototki(2)*drhototkj(2) + dt2*drhototkij2(2,2)
                derv2(jz,iy) = derv2(jz,iy) + dt1*drhototki(2)*drhototkj(3) + dt2*drhototkij2(2,3)
                derv2(jx,iz) = derv2(jx,iz) + dt1*drhototki(3)*drhototkj(1) + dt2*drhototkij2(3,1)
                derv2(jy,iz) = derv2(jy,iz) + dt1*drhototki(3)*drhototkj(2) + dt2*drhototkij2(3,2)
                derv2(jz,iz) = derv2(jz,iz) + dt1*drhototki(3)*drhototkj(3) + dt2*drhototkij2(3,3)
              else
                dt1 = rscrhok3*ofctijk
                derv2(jx,ix) = derv2(jx,ix) + dt1*drhototki(1)*drhototkj(1)
                derv2(jy,ix) = derv2(jy,ix) + dt1*drhototki(1)*drhototkj(2)
                derv2(jz,ix) = derv2(jz,ix) + dt1*drhototki(1)*drhototkj(3)
                derv2(jx,iy) = derv2(jx,iy) + dt1*drhototki(2)*drhototkj(1)
                derv2(jy,iy) = derv2(jy,iy) + dt1*drhototki(2)*drhototkj(2)
                derv2(jz,iy) = derv2(jz,iy) + dt1*drhototki(2)*drhototkj(3)
                derv2(jx,iz) = derv2(jx,iz) + dt1*drhototki(3)*drhototkj(1)
                derv2(jy,iz) = derv2(jy,iz) + dt1*drhototki(3)*drhototkj(2)
                derv2(jz,iz) = derv2(jz,iz) + dt1*drhototki(3)*drhototkj(3)
              endif
            endif
!******************************
!  End of second derivatives  *
!******************************
!***************************
!  End of third atom loop  *
!***************************
          enddo kloop
        endif
!****************************************
!  End of valid distance i - j section  *
!****************************************
      endif
    enddo jloop
  enddo iloop
!
!  End of real space part  -  perform general tasks
!
999 continue
!
!  Free local memory
!
  if (lgrad2) then
    deallocate(npotjkptr,stat=status)
    if (status/=0) call deallocate_error('many0dd','npotjkptr')
    deallocate(npotikptr,stat=status)
    if (status/=0) call deallocate_error('many0dd','npotikptr')
  endif
  deallocate(npotl,stat=status)
  if (status/=0) call deallocate_error('many0dd','npotl')
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
  call trace_out('many0dd')
#endif
!
  return
  end
