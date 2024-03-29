  subroutine manyfcd
!
!  Subroutine for calculating the many-body phonons from the
!  Sutton-Chen potential(s). Requires real space routine to
!  be called first to evaluate the repulsive contribution
!  and the rho parameters. 
!  Unphased version for fc_supercell option. 
!  Distributed memory parallel version.
!
!  On entry the array scrho must contain the density at each atomic site.
!
!   4/17 Created from manyfc
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
  use control
  use current
  use derivatives
  use eam
  use general
  use parallel
  use realvectors,    only : dist, dist2, xtmp, ytmp, ztmp, xtmp2, ytmp2, ztmp2
  use realvectors,    only : cellindex
  use sutton
  use symmetry
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  use two
  use vectors,        only : vector_pair
  implicit none
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ii
  integer(i4)                                  :: iloc
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: j
  integer(i4)                                  :: jj
  integer(i4)                                  :: jx
  integer(i4)                                  :: jy
  integer(i4)                                  :: jz
  integer(i4)                                  :: k
  integer(i4)                                  :: m
  integer(i4)                                  :: n  
  integer(i4)                                  :: nat1
  integer(i4)                                  :: nat2
  integer(i4)                                  :: nati
  integer(i4)                                  :: natj
  integer(i4)                                  :: natk
  integer(i4)                                  :: ncindm
  integer(i4)                                  :: ncindp
  integer(i4)                                  :: neamspeci
  integer(i4)                                  :: neamspecj
  integer(i4)                                  :: neamspeck
  integer(i4)                                  :: nmolonly
  integer(i4)                                  :: nor
  integer(i4)                                  :: nor2
  integer(i4)                                  :: npot
  integer(i4)                                  :: npotijk
  integer(i4)                                  :: npots
  integer(i4)                                  :: ntyp1  
  integer(i4)                                  :: ntyp2
  integer(i4)                                  :: ntypi
  integer(i4)                                  :: ntypj    
  integer(i4)                                  :: ntypk
  integer(i4), dimension(:), allocatable       :: npotl
  integer(i4)                                  :: status
  logical                                      :: lanyvalidik
  logical                                      :: lanyvalidjk
  logical                                      :: lself 
  logical                                      :: lvalidij
  logical                                      :: lvalidik
  logical                                      :: lvalidjk
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: cut2
  real(dp)                                     :: cut2k
  real(dp)                                     :: cut2p
  real(dp)                                     :: cut2r
  real(dp)                                     :: cut2rk
  real(dp)                                     :: d2r(6)
  real(dp)                                     :: deriv(3)
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
  real(dp)                                     :: dt1r
  real(dp)                                     :: dt2r
  real(dp)                                     :: eeam
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
  real(dp)                                     :: rjk
  real(dp)                                     :: rik2
  real(dp)                                     :: rjk2
  real(dp)                                     :: rp
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
  real(dp)                                     :: xcrd 
  real(dp)                                     :: ycrd    
  real(dp)                                     :: zcrd
  real(dp)                                     :: xcrdik
  real(dp)                                     :: ycrdik
  real(dp)                                     :: zcrdik
  real(dp)                                     :: xcrdjk
  real(dp)                                     :: ycrdjk
  real(dp)                                     :: zcrdjk
#ifdef TRACE
  call trace_in('manyfcd')
#endif
!
  time1 = g_cpu_time()
!***************************
!  Set up local variables  *
!***************************
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
  if (lnoreal) goto 1000
!
!  Scale density
!
  call eamscalescrho(1_i4)
!
!  Allocate local memory
!
  allocate(npotl(npote),stat=status)
  if (status/=0) call outofmemory('manyfcd','npotl')
!
!  Outer loop over sites
!
  ix = -2
  iy = -1
  iz =  0
  iloop: do iloc = 1,natomsonnode
    i = node2atom(iloc)
    ix = ix + 3
    iy = iy + 3
    iz = iz + 3
!     
!  Find EAM species for i
!  
    neamspeci = neamfnspecptr(i)
!
!  Evaluate functional derivatives
!
    if (lMEAMfn) then
      call meamfnderv(neamfn,neamspeci,scrho(1,i),rhoi,eeam,rscrhoi,rscrhoi3,rscrhoi5,.true.,.true.,.false.)
    else
      call eamfnderv(neamfn,neamspeci,scrho(1,i),eeam,rscrhoi,rscrhoi3,rscrhoi5,.true.,.true.,.false.)
      rhoi = scrho(1,i)
    endif
!
    xal = xclat(i)
    yal = yclat(i)
    zal = zclat(i)
    nati = nat(i)
    ntypi = nftype(i)
    oci = occuf(i)
!
!  Start of second atom loop
!
    jx = - 2
    jy = - 1
    jz =   0
    jloop: do j = 1,numat
      jx = jx + 3
      jy = jy + 3
      jz = jz + 3
!     
!  Find EAM species for j
!  
      neamspecj = neamfnspecptr(j)
!
!  Evaluate functional derivatives
!
      if (lMEAMfn) then
        call meamfnderv(neamfn,neamspecj,scrho(1,j),rhoj,eeam,rscrhoj,rscrhoj3,rscrhoj5,.true.,.true.,.false.)
      else
        call eamfnderv(neamfn,neamspecj,scrho(1,j),eeam,rscrhoj,rscrhoj3,rscrhoj5,.true.,.true.,.false.)
        rhoj = scrho(1,j)
      endif
!
      natj = nat(j)
      ntypj = nftype(j)
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
      xcrd = xclat(j) - xal
      ycrd = yclat(j) - yal
      zcrd = zclat(j) - zal
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
!  Need to make cut-off equal to double the maximum
!  to ensure all triangles are included
!
      rp = 2.0_dp*scmax
      cut2r = rp*rp
      if (cut2r.gt.4.0_dp*cut2p) cut2r = cut2p
      cut2 = cut2r
!*********************************
!  Find valid vectors for i - j  *
!*********************************
      call rfind(xcrd,ycrd,zcrd,cut2,0.0_dp,0.0_dp,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nmolonly,lself,0_i4,nor)
!
!  Loop over valid vectors
!
      do ii = 1,nor
        r2 = dist(ii)   
        r = sqrt(r2)
        xcd = xtmp(ii)
        ycd = ytmp(ii)
        zcd = ztmp(ii)
!
!  Find cell index
!
        if (abs(cellindex(1,ii)).gt.nd2cell(1).or. &
            abs(cellindex(2,ii)).gt.nd2cell(2).or. &
            abs(cellindex(3,ii)).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
          ncindm = nd2central
          ncindp = nd2central
        else
!
!  Compute index
!
          ncindp = nd2cellptr(nd2cell(1)+1+cellindex(1,ii),nd2cell(2)+1+cellindex(2,ii),nd2cell(3)+1+cellindex(3,ii))
          ncindm = nd2cellptr(nd2cell(1)+1-cellindex(1,ii),nd2cell(2)+1-cellindex(2,ii),nd2cell(3)+1-cellindex(3,ii))
        endif
!***************************************************
!  Calculate many-body contribution in real space  *
!***************************************************
        deriv(1:3) = 0.0_dp
        deriv2(1:6) = 0.0_dp
!***************************************
!  Valid many-body potentials for i-j  *
!***************************************
        if (lMEAM) then
          rhoij(1:maxmeamcomponent) = 0.0_dp
          rhoji(1:maxmeamcomponent) = 0.0_dp
          drhoij(1:3,1:maxmeamcomponent) = 0.0_dp
          drhoji(1:3,1:maxmeamcomponent) = 0.0_dp
          drhoij2(1:6,1:maxmeamcomponent) = 0.0_dp
          drhoji2(1:6,1:maxmeamcomponent) = 0.0_dp
        else
          rhoij(1) = 0.0_dp
          rhoji(1) = 0.0_dp
        endif
        drhototij(1:3) = 0.0_dp
        drhototji(1:3) = 0.0_dp
        drhototij2(1:6) = 0.0_dp
        drhototji2(1:6) = 0.0_dp
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
                             1.0_dp,1.0_dp,.true.,.false.,.true.,.true.,twopot(1,npot))
                call meamtotalrhoderv(neamspeci,scrho(1,i),rhoi,drhoij,drhototij,drhoijs,drhototijs, &
                                      drhoij2,drhototij2,drhoij2s,drhototij2s,drhoij2m,drhototij2m, &
                                      .false.,.true.)
                call meamtotalrhoderv(neamspecj,scrho(1,j),rhoj,drhoji,drhototji,drhojis,drhototjis, &
                                      drhoji2,drhototji2,drhoji2s,drhototji2s,drhoji2m,drhototji2m, &
                                      .false.,.true.)
              else
                call eamrho(nati,ntypi,natj,ntypj,r,rpot(npot),xcd,ycd,zcd,rhoij,rhoji,drhototij,drhototji, &
                            drhototijs,drhototjis,drhototij2,drhototji2,drhototij2s,drhototji2s, &
                            drhototij2m,drhototji2m,drhototij3,drhototji3,1.0_dp,1.0_dp,.false.,.true.,.true.,.false., &
                            twopot(1,npot))
              endif
            endif
          enddo
!
!  Combine derivative terms
!
          deriv(1:3) = deriv(1:3) + (rscrhoi*drhototij(1:3) + rscrhoj*drhototji(1:3))*ofct
          deriv2(1:6) = deriv2(1:6) + (rscrhoi*drhototij2(1:6) + rscrhoj*drhototji2(1:6))*ofct
          deriv2(1) = deriv2(1) + ocj*rscrhoi3*drhototij(1)*drhototij(1)*ofct + oci*rscrhoj3*drhototji(1)*drhototji(1)*ofct
          deriv2(2) = deriv2(2) + ocj*rscrhoi3*drhototij(1)*drhototij(2)*ofct + oci*rscrhoj3*drhototji(1)*drhototji(2)*ofct
          deriv2(3) = deriv2(3) + ocj*rscrhoi3*drhototij(1)*drhototij(3)*ofct + oci*rscrhoj3*drhototji(1)*drhototji(3)*ofct
          deriv2(4) = deriv2(4) + ocj*rscrhoi3*drhototij(2)*drhototij(2)*ofct + oci*rscrhoj3*drhototji(2)*drhototji(2)*ofct
          deriv2(5) = deriv2(5) + ocj*rscrhoi3*drhototij(2)*drhototij(3)*ofct + oci*rscrhoj3*drhototji(2)*drhototji(3)*ofct
          deriv2(6) = deriv2(6) + ocj*rscrhoi3*drhototij(3)*drhototij(3)*ofct + oci*rscrhoj3*drhototji(3)*drhototji(3)*ofct
        endif

        if (lvalidij) then
!
!  Generate products for derivatives between i - j
!
          d2r(1:6) = deriv2(1:6)
!******************************
!  Phased second derivatives  *
!******************************
          d2cell(jx,ix,ncindp) = d2cell(jx,ix,ncindp) - deriv2(1)
          d2cell(jy,ix,ncindp) = d2cell(jy,ix,ncindp) - deriv2(2)
          d2cell(jz,ix,ncindp) = d2cell(jz,ix,ncindp) - deriv2(3)
          d2cell(jx,iy,ncindp) = d2cell(jx,iy,ncindp) - deriv2(2)
          d2cell(jy,iy,ncindp) = d2cell(jy,iy,ncindp) - deriv2(4)
          d2cell(jz,iy,ncindp) = d2cell(jz,iy,ncindp) - deriv2(5)
          d2cell(jx,iz,ncindp) = d2cell(jx,iz,ncindp) - deriv2(3)
          d2cell(jy,iz,ncindp) = d2cell(jy,iz,ncindp) - deriv2(5)
          d2cell(jz,iz,ncindp) = d2cell(jz,iz,ncindp) - deriv2(6)
        endif
!******************************************************************
!  Start of third atom loop - only needed for second derivatives  *
!******************************************************************
        kloop: do k = 1,numat
!
          natk = nat(k)
          ntypk = nftype(k)
          xcrdik = xclat(k) - xal
          ycrdik = yclat(k) - yal
          zcrdik = zclat(k) - zal
          xcrdjk = xcrdik - xcd
          ycrdjk = ycrdik - ycd
          zcrdjk = zcrdik - zcd
          ock = occuf(k)
          ofctijk = oci*ocj*ock
!
!  Check whether there are any potentials between i-k or j-k
!
          npotijk = 0
          rpijk = 0.0_dp
          do n = 1,npote
            if (nptype(n).eq.19) then
              if (natk.eq.nspec1(n).and.(ntypk.eq.nptyp1(n).or.nptyp1(n).eq.0)) then
                npotijk = npotijk + 1
                if (rpot(n).gt.rpijk) rpijk = rpot(n)
              elseif (natk.eq.nspec2(n).and.(ntypk.eq.nptyp2(n).or.nptyp2(n).eq.0)) then
                npotijk = npotijk + 1
                if (rpot(n).gt.rpijk) rpijk = rpot(n)
              endif
            endif
          enddo
!
!  If no valid potentials for i-k or j-k then skip
!
          if (npotijk.eq.0) cycle kloop
          cut2rk = rpijk*rpijk
          if (cut2rk.gt.cut2p) cut2rk = cut2p
          cut2k = cut2rk
!     
!  Find EAM species for k
!  
          neamspeck = neamfnspecptr(k)
!
!  Evaluate functional derivatives
!
          if (lMEAMfn) then
            call meamfnderv(neamfn,neamspeck,scrho(1,k),rhok,eeam,rscrhok,rscrhok3,rscrhok5,.true.,.true.,.false.)
          else
            call eamfnderv(neamfn,neamspeck,scrho(1,k),eeam,rscrhok,rscrhok3,rscrhok5,.true.,.true.,.false.)
            rhok = scrho(1,k)
          endif
!
!  If no rho then skip
!
          if (rhoi.eq.0.0_dp.and.rhoj.eq.0.0_dp.and.rhok.eq.0.0_dp) cycle kloop
!***********************
!  General cell loops  *
!***********************
          call rfindeither(xcrdik,ycrdik,zcrdik,xcrdjk,ycrdjk,zcrdjk,cut2k,cut2k,.false.,nor,nor2)
          do jj = 1,nor2
            rik2 = dist(nor+jj)
            rjk2 = dist2(nor+jj)
            xcd1 = xtmp(nor+jj)
            ycd1 = ytmp(nor+jj)    
            zcd1 = ztmp(nor+jj)  
            xcd2 = xtmp2(nor+jj)
            ycd2 = ytmp2(nor+jj)
            zcd2 = ztmp2(nor+jj)
!************************************************************
!  Calculate triangular contribution to second derivatives  *
!************************************************************
            if (lMEAM) then
              rhoik(1:maxmeamcomponent) = 0.0_dp
              rhoki(1:maxmeamcomponent) = 0.0_dp
              rhojk(1:maxmeamcomponent) = 0.0_dp
              rhokj(1:maxmeamcomponent) = 0.0_dp
              drhoik(1:3,1:maxmeamcomponent) = 0.0_dp
              drhoki(1:3,1:maxmeamcomponent) = 0.0_dp
              drhojk(1:3,1:maxmeamcomponent) = 0.0_dp
              drhokj(1:3,1:maxmeamcomponent) = 0.0_dp
            else
              rhoik(1) = 0.0_dp
              rhoki(1) = 0.0_dp
              rhojk(1) = 0.0_dp
              rhokj(1) = 0.0_dp
            endif
            drhototik(1:3) = 0.0_dp
            drhototki(1:3) = 0.0_dp
            drhototjk(1:3) = 0.0_dp
            drhototkj(1:3) = 0.0_dp
            drhototijk2(1:3,1:3) = 0.0_dp
            drhototjik2(1:3,1:3) = 0.0_dp
            drhototkij2(1:3,1:3) = 0.0_dp
!
            lanyvalidik = .false.
            lanyvalidjk = .false.
            if (rik2.le.cut2k) then
!*********************
!  i-k contribution  *
!*********************
              rik = sqrt(rik2)
!
!  Loop over potentials to find many-body ones
!
              do m = 1,npote
                if (nptype(m).eq.19) then
                  lvalidik = .false.
                  if (rik.gt.rpot2(m).and.rik.le.rpot(m)) then
                    if (nati.eq.nspec1(m).and.natk.eq.nspec2(m)) then
                      if (ntypi.eq.nptyp1(m).or.nptyp1(m).eq.0) then
                        if (ntypk.eq.nptyp2(m).or.nptyp2(m).eq.0) lvalidik = .true.
                      endif
                    elseif (nati.eq.nspec2(m).and.natk.eq.nspec1(m)) then
                      if (ntypi.eq.nptyp2(m).or.nptyp2(m).eq.0) then
                        if (ntypk.eq.nptyp1(m).or.nptyp1(m).eq.0) lvalidik = .true.
                      endif
                    endif
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
                        call meamtotalrhocrossderv(neamspeci,scrho(1,i),rhoi,drhoij,drhototij,drhoik,drhototik,drhoijs,drhoiks, &
                                                   drhototijk2,drhototijk2s,drhototijk2m,lstr)
                      else
                        call eamrho(nati,ntypi,natk,ntypk,rik,rpot(m),xcd1,ycd1,zcd1,rhoik,rhoki,drhototik,drhototki, &
                                    drhototiks,drhototkis,drhototik2,drhototki2,drhototik2s,drhototki2s, &
                                    drhototik2m,drhototki2m,drhototik3,drhototki3,1.0_dp,1.0_dp,.false.,.true.,.false.,.false., &
                                    twopot(1,m))
                      endif
                      lanyvalidik = .true.
                    endif
                  endif
                endif
              enddo
            endif
            if (rjk2.le.cut2) then
!*********************
!  j-k contribution  *
!*********************
              rjk = sqrt(rjk2)
!
!  Loop over potentials to find many-body ones
!
              do m = 1,npote
                if (nptype(m).eq.19) then
                  lvalidjk = .false.
                  if (rjk.gt.rpot2(m).and.rjk.le.rpot(m)) then
                    if (natj.eq.nspec1(m).and.natk.eq.nspec2(m)) then
                      if (ntypj.eq.nptyp1(m).or.nptyp1(m).eq.0) then
                        if (ntypk.eq.nptyp2(m).or.nptyp2(m).eq.0) lvalidjk = .true.
                      endif
                    elseif (natj.eq.nspec2(m).and.natk.eq.nspec1(m)) then
                      if (ntypj.eq.nptyp2(m).or.nptyp2(m).eq.0) then
                        if (ntypk.eq.nptyp1(m).or.nptyp1(m).eq.0) lvalidjk = .true.
                      endif
                    endif
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
                        call meamtotalrhocrossderv(neamspecj,scrho(1,j),rhoj,drhoji,drhototji,drhojk,drhototjk,drhojis,drhojks, &
                                                   drhototjik2,drhototjik2s,drhototjik2m,lstr)
                      else
                        call eamrho(natj,ntypj,natk,ntypk,rjk,rpot(m),xcd2,ycd2,zcd2,rhojk,rhokj,drhototjk,drhototkj, &
                                    drhototjks,drhototkjs,drhototjk2,drhototkj2,drhototjk2s,drhototkj2s, &
                                    drhototjk2m,drhototkj2m,drhototjk3,drhototkj3,1.0_dp,1.0_dp,.false.,.true.,.false.,.false., &
                                    twopot(1,m))
                      endif
                      lanyvalidjk = .true.
                    endif
                  endif
                endif
              enddo
            endif
!
!  Cross term derivative for k-i / k-j
!
            if (lanyvalidik.and.lanyvalidjk.and.lMEAMden) then
              call meamtotalrhocrossderv(neamspeck,scrho(1,k),rhok,drhoki,drhototki,drhokj,drhototkj,drhokis,drhokjs, &
                                         drhototkij2,drhototkij2s,drhototkij2m,lstr)
            endif
!***************************************************
!  Calculate second derivatives for i-k/j-k terms  *
!***************************************************
!
!  i-k
!
            if (lMEAM) then
              dt1r = rscrhoi3*ofctijk
              dt2r = rscrhoi*ofctijk
              d2cell(jx,ix,ncindp) = d2cell(jx,ix,ncindp) - dt1r*drhototik(1)*drhototij(1) - dt2r*drhototijk2(1,1)
              d2cell(jy,ix,ncindp) = d2cell(jy,ix,ncindp) - dt1r*drhototik(1)*drhototij(2) - dt2r*drhototijk2(2,1)
              d2cell(jz,ix,ncindp) = d2cell(jz,ix,ncindp) - dt1r*drhototik(1)*drhototij(3) - dt2r*drhototijk2(3,1)
              d2cell(jx,iy,ncindp) = d2cell(jx,iy,ncindp) - dt1r*drhototik(2)*drhototij(1) - dt2r*drhototijk2(1,2)
              d2cell(jy,iy,ncindp) = d2cell(jy,iy,ncindp) - dt1r*drhototik(2)*drhototij(2) - dt2r*drhototijk2(2,2)
              d2cell(jz,iy,ncindp) = d2cell(jz,iy,ncindp) - dt1r*drhototik(2)*drhototij(3) - dt2r*drhototijk2(3,2)
              d2cell(jx,iz,ncindp) = d2cell(jx,iz,ncindp) - dt1r*drhototik(3)*drhototij(1) - dt2r*drhototijk2(1,3)
              d2cell(jy,iz,ncindp) = d2cell(jy,iz,ncindp) - dt1r*drhototik(3)*drhototij(2) - dt2r*drhototijk2(2,3)
              d2cell(jz,iz,ncindp) = d2cell(jz,iz,ncindp) - dt1r*drhototik(3)*drhototij(3) - dt2r*drhototijk2(3,3)
            else
              dt1r = rscrhoi3*ofctijk
              d2cell(jx,ix,ncindp) = d2cell(jx,ix,ncindp) - dt1r*drhototik(1)*drhototij(1)
              d2cell(jy,ix,ncindp) = d2cell(jy,ix,ncindp) - dt1r*drhototik(1)*drhototij(2)
              d2cell(jz,ix,ncindp) = d2cell(jz,ix,ncindp) - dt1r*drhototik(1)*drhototij(3)
              d2cell(jx,iy,ncindp) = d2cell(jx,iy,ncindp) - dt1r*drhototik(2)*drhototij(1)
              d2cell(jy,iy,ncindp) = d2cell(jy,iy,ncindp) - dt1r*drhototik(2)*drhototij(2)
              d2cell(jz,iy,ncindp) = d2cell(jz,iy,ncindp) - dt1r*drhototik(2)*drhototij(3)
              d2cell(jx,iz,ncindp) = d2cell(jx,iz,ncindp) - dt1r*drhototik(3)*drhototij(1)
              d2cell(jy,iz,ncindp) = d2cell(jy,iz,ncindp) - dt1r*drhototik(3)*drhototij(2)
              d2cell(jz,iz,ncindp) = d2cell(jz,iz,ncindp) - dt1r*drhototik(3)*drhototij(3)
            endif
!
!  j-k
!
            if (lMEAM) then
              dt1r = rscrhoj3*ofctijk
              dt2r = rscrhoj*ofctijk
              d2cell(jx,ix,ncindp) = d2cell(jx,ix,ncindp) + dt1r*drhototjk(1)*drhototji(1) + dt2r*drhototjik2(1,1)
              d2cell(jy,ix,ncindp) = d2cell(jy,ix,ncindp) + dt1r*drhototjk(2)*drhototji(1) + dt2r*drhototjik2(1,2)
              d2cell(jz,ix,ncindp) = d2cell(jz,ix,ncindp) + dt1r*drhototjk(3)*drhototji(1) + dt2r*drhototjik2(1,3)
              d2cell(jx,iy,ncindp) = d2cell(jx,iy,ncindp) + dt1r*drhototjk(1)*drhototji(2) + dt2r*drhototjik2(2,1)
              d2cell(jy,iy,ncindp) = d2cell(jy,iy,ncindp) + dt1r*drhototjk(2)*drhototji(2) + dt2r*drhototjik2(2,2)
              d2cell(jz,iy,ncindp) = d2cell(jz,iy,ncindp) + dt1r*drhototjk(3)*drhototji(2) + dt2r*drhototjik2(2,3)
              d2cell(jx,iz,ncindp) = d2cell(jx,iz,ncindp) + dt1r*drhototjk(1)*drhototji(3) + dt2r*drhototjik2(3,1)
              d2cell(jy,iz,ncindp) = d2cell(jy,iz,ncindp) + dt1r*drhototjk(2)*drhototji(3) + dt2r*drhototjik2(3,2)
              d2cell(jz,iz,ncindp) = d2cell(jz,iz,ncindp) + dt1r*drhototjk(3)*drhototji(3) + dt2r*drhototjik2(3,3)
            else
              dt1r = rscrhoj3*ofctijk
              d2cell(jx,ix,ncindp) = d2cell(jx,ix,ncindp) + dt1r*drhototjk(1)*drhototji(1)
              d2cell(jy,ix,ncindp) = d2cell(jy,ix,ncindp) + dt1r*drhototjk(2)*drhototji(1)
              d2cell(jz,ix,ncindp) = d2cell(jz,ix,ncindp) + dt1r*drhototjk(3)*drhototji(1)
              d2cell(jx,iy,ncindp) = d2cell(jx,iy,ncindp) + dt1r*drhototjk(1)*drhototji(2)
              d2cell(jy,iy,ncindp) = d2cell(jy,iy,ncindp) + dt1r*drhototjk(2)*drhototji(2)
              d2cell(jz,iy,ncindp) = d2cell(jz,iy,ncindp) + dt1r*drhototjk(3)*drhototji(2)
              d2cell(jx,iz,ncindp) = d2cell(jx,iz,ncindp) + dt1r*drhototjk(1)*drhototji(3)
              d2cell(jy,iz,ncindp) = d2cell(jy,iz,ncindp) + dt1r*drhototjk(2)*drhototji(3)
              d2cell(jz,iz,ncindp) = d2cell(jz,iz,ncindp) + dt1r*drhototjk(3)*drhototji(3)
            endif
!
!  i-k/j-k
!
            if (lMEAM) then
              dt1r = rscrhok3*ofctijk
              dt2r = rscrhok*ofctijk
              d2cell(jx,ix,ncindp) = d2cell(jx,ix,ncindp) + dt1r*drhototki(1)*drhototkj(1) + dt2r*drhototkij2(1,1)
              d2cell(jy,ix,ncindp) = d2cell(jy,ix,ncindp) + dt1r*drhototki(1)*drhototkj(2) + dt2r*drhototkij2(1,2)
              d2cell(jz,ix,ncindp) = d2cell(jz,ix,ncindp) + dt1r*drhototki(1)*drhototkj(3) + dt2r*drhototkij2(1,3)
              d2cell(jx,iy,ncindp) = d2cell(jx,iy,ncindp) + dt1r*drhototki(2)*drhototkj(1) + dt2r*drhototkij2(2,1)
              d2cell(jy,iy,ncindp) = d2cell(jy,iy,ncindp) + dt1r*drhototki(2)*drhototkj(2) + dt2r*drhototkij2(2,2)
              d2cell(jz,iy,ncindp) = d2cell(jz,iy,ncindp) + dt1r*drhototki(2)*drhototkj(3) + dt2r*drhototkij2(2,3)
              d2cell(jx,iz,ncindp) = d2cell(jx,iz,ncindp) + dt1r*drhototki(3)*drhototkj(1) + dt2r*drhototkij2(3,1)
              d2cell(jy,iz,ncindp) = d2cell(jy,iz,ncindp) + dt1r*drhototki(3)*drhototkj(2) + dt2r*drhototkij2(3,2)
              d2cell(jz,iz,ncindp) = d2cell(jz,iz,ncindp) + dt1r*drhototki(3)*drhototkj(3) + dt2r*drhototkij2(3,3)
            else
              dt1r = rscrhok3*ofctijk
              d2cell(jx,ix,ncindp) = d2cell(jx,ix,ncindp) + dt1r*drhototki(1)*drhototkj(1)
              d2cell(jy,ix,ncindp) = d2cell(jy,ix,ncindp) + dt1r*drhototki(1)*drhototkj(2)
              d2cell(jz,ix,ncindp) = d2cell(jz,ix,ncindp) + dt1r*drhototki(1)*drhototkj(3)
              d2cell(jx,iy,ncindp) = d2cell(jx,iy,ncindp) + dt1r*drhototki(2)*drhototkj(1)
              d2cell(jy,iy,ncindp) = d2cell(jy,iy,ncindp) + dt1r*drhototki(2)*drhototkj(2)
              d2cell(jz,iy,ncindp) = d2cell(jz,iy,ncindp) + dt1r*drhototki(2)*drhototkj(3)
              d2cell(jx,iz,ncindp) = d2cell(jx,iz,ncindp) + dt1r*drhototki(3)*drhototkj(1)
              d2cell(jy,iz,ncindp) = d2cell(jy,iz,ncindp) + dt1r*drhototki(3)*drhototkj(2)
              d2cell(jz,iz,ncindp) = d2cell(jz,iz,ncindp) + dt1r*drhototki(3)*drhototkj(3)
            endif
          enddo
!***************************
!  End of third atom loop  *
!***************************
        enddo kloop
!**************************************
!  End of valid distance i-j section  *
!**************************************
      enddo
    enddo jloop
  enddo iloop
!
!  Free local memory
!
  deallocate(npotl,stat=status)
  if (status/=0) call deallocate_error('manyfcd','npotl')
!
!  Unscale density
!
  call eamscalescrho(-1_i4)
!
!  Exit point
!
1000 continue
!
!  Timing
!
  time2 = g_cpu_time()
  tmany = tmany + time2 - time1
#ifdef TRACE
  call trace_out('manyfcd')
#endif
!
  return
  end
