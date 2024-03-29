  subroutine real12a2d(eatom,ereal,erecip,ec6,lgrad1,lgrad2,imode)
!
!  Subroutine for calculating defect real energy and up to second derivatives
!
!  Distributed memory parallel version.
!
!  xdefe = region 1 coordinates
!  xclat = region 2 coordinates per unit cell
!
!  imode = 1 => defective region 1 - region 2 - only mode where derivatives are needed
!  imode = 2 => perfect region 1 - region 2
!
!   5/17 Created from real12a2
!   5/17 Commented code removed 
!   5/17 Global summation of dscrho added
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
  use configurations, only : lbsmat,radcfg
  use g_constants
  use control
  use current
  use defects
  use derivatives
  use eam,            only : lMEAMden
  use element,        only : maxele
  use general,        only : cutw, etaw
  use kspace
  use molecule
  use numbers,        only : third
  use parallel
  use region2a
  use shells
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
  integer(i4), intent(in)                      :: imode
  logical,     intent(in)                      :: lgrad1
  logical,     intent(in)                      :: lgrad2
  real(dp),    intent(inout)                   :: eatom
  real(dp),    intent(inout)                   :: ec6
  real(dp),    intent(inout)                   :: ereal
  real(dp),    intent(inout)                   :: erecip
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: iar
  integer(i4)                                  :: iloc
  integer(i4)                                  :: ind
  integer(i4)                                  :: indi
  integer(i4)                                  :: indj
  integer(i4)                                  :: indmi
  integer(i4)                                  :: indmj
  integer(i4)                                  :: indmmi
  integer(i4)                                  :: indmmj
  integer(i4)                                  :: indri
  integer(i4)                                  :: indrif
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: ixx
  integer(i4)                                  :: iyy
  integer(i4)                                  :: izz
  integer(i4)                                  :: ixx2
  integer(i4)                                  :: iyy2
  integer(i4)                                  :: izz2
  integer(i4)                                  :: j
  integer(i4)                                  :: jx
  integer(i4)                                  :: jy
  integer(i4)                                  :: jz
  integer(i4)                                  :: jxx
  integer(i4)                                  :: jyy
  integer(i4)                                  :: jzz
  integer(i4)                                  :: jxx2
  integer(i4)                                  :: jyy2
  integer(i4)                                  :: jzz2
  integer(i4)                                  :: n
  integer(i4)                                  :: n1
  integer(i4)                                  :: nat1
  integer(i4)                                  :: nat2
  integer(i4)                                  :: nati
  integer(i4)                                  :: natj
  integer(i4)                                  :: nbtypeij
  integer(i4)                                  :: nbtypeij2
  integer(i4)                                  :: nloop
  integer(i4)                                  :: nmi
  integer(i4)                                  :: nmj
  integer(i4)                                  :: nmolonly
  integer(i4)                                  :: npsi
  integer(i4)                                  :: npsj
  integer(i4), dimension(:), allocatable       :: npotl
  integer(i4)                                  :: npots
  integer(i4)                                  :: nstep
  integer(i4)                                  :: ntyp1
  integer(i4)                                  :: ntyp2
  integer(i4)                                  :: ntypi
  integer(i4)                                  :: ntypj
  integer(i4)                                  :: status
  logical                                      :: l2bonds
  logical                                      :: l3bonds
  logical                                      :: lbonded
  logical                                      :: lcspair
  logical                                      :: lmatch
  logical                                      :: lmolok
  logical                                      :: lneedmol
  logical                                      :: lorder12loc
  logical                                      :: lptrmol
  logical                                      :: lregion1j
  real(dp)                                     :: c6prod
  real(dp)                                     :: c6self1
  real(dp)                                     :: c6self1d2
  real(dp)                                     :: c6tot
  real(dp)                                     :: cmax
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: cut2
  real(dp)                                     :: cut2e
  real(dp)                                     :: cut2p
  real(dp)                                     :: cut2q
  real(dp)                                     :: cut2r
  real(dp)                                     :: cut2s
  real(dp)                                     :: cut2w
  real(dp)                                     :: d0i
  real(dp)                                     :: d0j
  real(dp)                                     :: d1i
  real(dp)                                     :: d1j
  real(dp)                                     :: d2i2
  real(dp)                                     :: d2ij
  real(dp)                                     :: d2j2
  real(dp)                                     :: deriv
  real(dp)                                     :: deriv2
  real(dp)                                     :: deriv3
  real(dp)                                     :: derive0
  real(dp)                                     :: derive
  real(dp)                                     :: derive2
  real(dp)                                     :: derive3
  real(dp)                                     :: dist
  real(dp),        dimension(:),   allocatable :: dscrholoc
  real(dp),        dimension(:),   allocatable :: dscrhor2dloc
  real(dp),        dimension(:),   allocatable :: dscrhor2ploc
  real(dp),        dimension(:),   allocatable :: dtmp
  real(dp)                                     :: eatom1
  real(dp)                                     :: ec61
  real(dp)                                     :: ereal1
  real(dp)                                     :: erecip1
  real(dp)                                     :: emax
  real(dp)                                     :: esum1(4)
  real(dp)                                     :: esum2(4)
  real(dp)                                     :: eta3
  real(dp)                                     :: factor
  real(dp)                                     :: oci
  real(dp)                                     :: ocj
  real(dp)                                     :: ofct
  real(dp)                                     :: ospfct
  real(dp)                                     :: qli
  real(dp)                                     :: qlj
  real(dp)                                     :: r
  real(dp)                                     :: radi
  real(dp)                                     :: radj
  real(dp)                                     :: radsum
  real(dp)                                     :: rcheck
  real(dp)                                     :: rcheck2
  real(dp)                                     :: rderiv
  real(dp)                                     :: rdiffc
  real(dp)                                     :: rmiddle2
  real(dp)                                     :: rp
  real(dp)                                     :: rpd1
  real(dp)                                     :: rpd2
  real(dp)                                     :: rpd3
  real(dp)                                     :: rpd4
  real(dp)                                     :: rpd5
  real(dp)                                     :: rpd6
  real(dp)                                     :: rtrm1
  real(dp)                                     :: rtrm2
  real(dp)                                     :: rtrm3
  real(dp)                                     :: rtrm32
  real(dp)                                     :: sctrm1
  real(dp)                                     :: sctrm2
  real(dp)                                     :: setrm
  real(dp)                                     :: small
  real(dp)                                     :: time1
  real(dp)                                     :: time2
  real(dp)                                     :: twoeta3
  real(dp)                                     :: xal
  real(dp)                                     :: yal
  real(dp)                                     :: zal
  real(dp)                                     :: xcrd
  real(dp)                                     :: ycrd
  real(dp)                                     :: zcrd
  real(dp)                                     :: xdiffc
  real(dp)                                     :: ydiffc
  real(dp)                                     :: zdiffc
#ifdef TRACE
  call trace_in('real12a2d')
#endif
!
  time1 = g_cpu_time()
!
!  Local variables
!
  small = 1.0d-12
  if (lwolf) then
    twoeta3 = 2.0_dp*etaw*etaw*third
  else
    twoeta3 = 2.0_dp*eta*third
  endif
  if (lc6) then
    eta3 = eta*eta*eta
    c6self1 = 0.5_dp*eta3*third
    if (lgrad2) then
      c6self1d2 = 0.25_dp*eta3*eta
    endif
  endif
!
!  Allocate local memory
!
  allocate(npotl(npote),stat=status)
  if (status/=0) call outofmemory('real12a2d','npotl')
  if (imode.eq.1) then
    allocate(dscrholoc(nreg1),stat=status)
    if (status/=0) call outofmemory('real12a2','dscrholoc')
    allocate(dscrhor2dloc(npreg2),stat=status)
    if (status/=0) call outofmemory('real12a2','dscrhor2dloc')
    dscrholoc(1:nreg1) = 0.0_dp
    dscrhor2dloc(1:npreg2) = 0.0_dp
  else
    allocate(dscrhor2ploc(npreg2),stat=status)
    if (status/=0) call outofmemory('real12a2','dscrhor2ploc')
    dscrhor2ploc(1:npreg2) = 0.0_dp
  endif
!
  eatom1 = 0.0_dp
  ereal1 = 0.0_dp
  erecip1 = 0.0_dp
  ec61 = 0.0_dp
!
!  Set up cutoffs
!
  cut2p = cutp*cutp
  cut2e = rmx2
  cut2s = cuts*cuts
  cut2w = cutw*cutw
  if (lwolf) then        
    cut2q = cut2w        
  else
    cut2q = cut2e
  endif
!
!  Decide whether all interactions lie within unit cell and first
!  neighbours - saves time for large systems
!
  emax = sqrt(cut2q)
  cmax = max(rpmax,emax)
!
  if (imode.eq.1) then
    nloop = nreg1onnode
    n1 = 1
    nstep = 1
  elseif (imode.eq.2) then
    nloop = nreg1old
    n1 = procid + 1
    nstep = nprocs
  endif
  if (lnoreal) goto 999
!
!  Set up pointers to non-region 1 second derivatives
!
  indj = 3*nreg1
  if (ldbsm) indj = indj + nreg1
  jx = indj + 1
  jy = indj + 2
  jz = indj + 3
!**********************************
!  Region 1 - region 1/2a energy  *
!**********************************
!
!  Outer loop over sites
!
  do iloc = n1,nloop,nstep
!
!  Set i attributes according to whether it is defective or perfect region 1
!
    if (imode.eq.1) then
      i = node2reg1(iloc)
      xal = xdefe(i)
      yal = ydefe(i)
      zal = zdefe(i)
      nati = natdefe(i)
      ntypi = ntypdefe(i)
      qli = qdefe(i)
      oci = occdefe(i)
      if (nreldef(i).gt.0) then
        npsi = nreldef(i)
      else
        npsi = 0
      endif
      if (ldefbsmat(i)) then
        radi = radefe(i)
        indri = 3*nloop + iloc
        indrif = 3*nreg1 + i
      else
        radi = 0.0_dp
      endif
      indi = 3*(iloc-1)
    elseif (imode.eq.2) then
      i = iloc
      xal = xperf(i)
      yal = yperf(i)
      zal = zperf(i)
      nati = natp(i)
      ntypi = ntypep(i)
      qli = qp(i)
      oci = occp(i)
      npsi = npsite(i)
      iar = nsft + nrelat(npsi)
      if (lbsmat(iar)) then
        radi = radcfg(iar)
        indri = 3*nloop + i
        indrif = 3*nloop + i
      else
        radi = 0.0_dp
      endif
      indi = 3*(i-1)
    endif
!
    ix = indi + 1
    iy = indi + 2
    iz = indi + 3
!
!  Molecule handling
!
    if (lmol) then
      if (imode.eq.1) then
        nmi = ndefmol(i)
        indmi = ndefind(i)
      elseif (imode.eq.2) then
        nmi = ndefmolp(i)
        indmi = ndefindp(i)
      endif
      call mindtoijk(indmi,ixx,iyy,izz)
    endif
!
!  Find distance from the centre
!
    xdiffc = xal-xdc
    ydiffc = yal-ydc
    zdiffc = zal-zdc
    rdiffc = xdiffc*xdiffc + ydiffc*ydiffc + zdiffc*zdiffc
    rdiffc = sqrt(rdiffc)
    rcheck = rdiffc + cmax + 0.5_dp
    rcheck2 = rcheck*rcheck
!***************************
!  Loop over old region 1  *
!***************************
    lregion1j = .true.
    do j = 1,nreg1old
      natj = natp(j)
      ntypj = ntypep(j)
!
!  Arrange nat1 to be lower than nat2. If they are equal
!  then arrange that ntyp1 is lower than ntyp2
!
      if (nati.eq.natj) then
        nat1 = nati
        nat2 = natj
        if (ntypi.lt.ntypj) then
          lorder12loc = .true.
          ntyp1 = ntypi
          ntyp2 = ntypj
        else
          lorder12loc = .false.
          ntyp1 = ntypj
          ntyp2 = ntypi
        endif
      elseif (nati.lt.natj) then
        lorder12loc = .true.
        nat1 = nati
        nat2 = natj
        ntyp1 = ntypi
        ntyp2 = ntypj
      else
        lorder12loc = .false.
        nat1 = natj
        nat2 = nati
        ntyp1 = ntypj
        ntyp2 = ntypi
      endif
      xcrd = xperf(j) - xal
      ycrd = yperf(j) - yal
      zcrd = zperf(j) - zal
      qlj = qp(j)
      ocj = occp(j)
      ofct = oci*ocj
      factor = qli*qlj*ofct*angstoev
      npsj = npsite(j)
      iar = nsft + nrelat(npsj)
      if (lbsmat(iar)) then
        radj = radcfg(iar)
      else
        radj = 0.0_dp
      endif
      radsum = radi + radj
      lbonded = .false.
      l2bonds = .false.
      l3bonds = .false.
!
!  Possible core-shell flag
!
      lcspair = (abs(nat1-nat2).eq.maxele.or.(oci+ocj).lt.1.0001_dp)
      if (lcspair) then
        ospfct = oci
      else
        ospfct = ofct
      endif
      if (lmol) then
!
!  Molecule handling
!
        nmj = ndefmolp(j)
        indmj = ndefindp(j)
        lmolok = (nmi.eq.nmj.and.nmi.gt.0)
      else
        lmolok = .false.
      endif
!
!  Locate potential number
!  Check whether potential requires specific types
!
      rp = 0.0_dp
      npots = 0
      c6tot = 0.0_dp
      lneedmol = (lmol.and..not.lmolq)
      do n = 1,npote
        if (lmatch(nat1,ntyp1,nspec1(n),nptyp1(n),.true.)) then
          if (lmatch(nat2,ntyp2,nspec2(n),nptyp2(n),.true.)) then
            npots = npots + 1
            npotl(npots)= n
            if (mmexc(n).gt.0.or.(lintra(n).and..not.linter(n).or..not.lintra(n).and.linter(n))) lneedmol = .true.
            if (nptype(n).eq.8.or.nptype(n).eq.33) then
              if (cuts.gt.rp) rp = cuts
            elseif (lc6) then
              if (nptype(n).eq.1.or.nptype(n).eq.7) then
                c6tot = c6tot + twopot(3,n)
                if (repcut(n).gt.rp) rp = repcut(n)
              elseif (nptype(n).eq.2.or.nptype(n).eq.21) then
                c6tot = c6tot + twopot(2,n)
                if (repcut(n).gt.rp) rp = repcut(n)
              elseif (nptype(n).eq.57) then
                c6tot = c6tot + twopot(3,n)*twopot(4,n)*twopot(5,n)**6
                if (repcut(n).gt.rp) rp = repcut(n)
              else
                if (rpot(n).gt.rp) rp = rpot(n)
              endif
            else
              if (rpot(n).gt.rp) rp = rpot(n)
            endif
          endif
        endif
      enddo
!
!  Generate looping indices
!
      cut2r = rp*rp
      if (cut2r.gt.cut2p) cut2r = cut2p
      cut2 = cut2r
      if (cut2e.gt.cut2.and.lewald) cut2 = cut2e
      if (cut2w.gt.cut2.and.lwolf) cut2 = cut2w
      rp = sqrt(cut2)
!
!  If the Coulomb term and the potentials don't need molecule info for
!  for this pair of atoms, then turn off molecule checking to save time
!
      if (.not.lneedmol) lmolok = .false.
!
!  Generate distance squared
!
      r = xcrd*xcrd + ycrd*ycrd + zcrd*zcrd
!
!  Molecule  -  check index
!
      nmolonly = 0
      if (lmolok.and.(r.gt.cut2s.or..not.lcspair)) then
        ind = indmj - indmi
        lptrmol = (ind.eq.0)
        if (.not.lptrmol) then
          call mindtoijk(indmj,jxx,jyy,jzz)
          jxx = jxx - ixx
          jyy = jyy - iyy
          jzz = jzz - izz
          call samemol(lptrmol,nmi,jxx,jyy,jzz,0_i4,0_i4,0_i4)
        endif
        if (lptrmol) then
          if (r.gt.cut2q) nmolonly = 1
          if (npsi.gt.0) then
            call mindtoijk(indmj,jxx,jyy,jzz)
            jxx = jxx - ixx
            jyy = jyy - iyy
            jzz = jzz - izz
            indmmj = nmolind(npsj)
            indmmi = nmolind(npsi)
            call mindtoijk(indmmj,jxx2,jyy2,jzz2)
            call mindtoijk(indmmi,ixx2,iyy2,izz2)
            jxx = jxx + jxx2 - ixx2
            jyy = jyy + jyy2 - iyy2
            jzz = jzz + jzz2 - izz2
            call bonded3(lbonded,l2bonds,l3bonds,nbtypeij,nbtypeij2,npsi,npsj,jxx,jyy,jzz)
          else
            lbonded   = .false.
            l2bonds   = .false.
            l3bonds   = .false.
            nbtypeij  = 0
            nbtypeij2 = 0
          endif
        else
          lbonded   = .false.
          l2bonds   = .false.
          l3bonds   = .false.
          nbtypeij  = 0
          nbtypeij2 = 0
        endif
      else
        lptrmol   = .false.
        lbonded   = .false.
        l2bonds   = .false.
        l3bonds   = .false.
        nbtypeij  = 0
        nbtypeij2 = 0
      endif
      if (r.lt.small) then
!
!  Remove electrostatic self energy
!
        setrm = factor*tweatpi
        erecip1 = erecip1 - setrm
        if (lc6) then
          c6prod = ofct*c6tot
          ec61 = ec61 + c6self1*c6prod
        endif
        if (lgrad2) then
          setrm = setrm*twoeta3
          if (lc6) setrm = setrm - c6prod*c6self1d2
          derv2(jx,ix) = derv2(jx,ix) - setrm
          derv2(jy,iy) = derv2(jy,iy) - setrm
          derv2(jz,iz) = derv2(jz,iz) - setrm
        endif
      elseif (r.le.cut2.or.lptrmol) then
!
!  Store vector
!
        dist = sqrt(r)
        call twobody1(eatom1,ereal1,ec61,lgrad1,lgrad2,.false.,1_i4,1_i4,dist,xcrd,ycrd,zcrd,d0i,d0j, &
                      deriv,deriv2,deriv3,derive0,derive,derive2,derive3,rderiv,npots,npotl,cut2r,cut2q, &
                      cut2s,lptrmol,nmolonly,factor,ofct,ospfct,radsum,rtrm1,rtrm2,rtrm3,rtrm32,sctrm1,sctrm2, &
                      qli,qlj,lcspair,.true.,lregion1j,lbonded,l2bonds,l3bonds,nbtypeij,nbtypeij2, &
                      lregion1j,.false.,d1i,d1j,d2i2,d2ij,d2j2)
!
!  Generate products for derivatives
!
        if (lgrad2) then
          rpd1 = xcrd*xcrd
          rpd2 = ycrd*ycrd
          rpd3 = zcrd*zcrd
          rpd4 = ycrd*zcrd
          rpd5 = xcrd*zcrd
          rpd6 = xcrd*ycrd
        endif
!
!  For mode 2 the terms are to be subtracted  = > change sign
!
        if (imode.eq.2) then
          if (.not.lMEAMden) then
            sctrm1 = - sctrm1
            sctrm2 = - sctrm2
          endif
          if (lgrad1) then
            deriv =  - deriv
            if (radsum.gt.0.0_dp) then
              rtrm1 =  - rtrm1
            endif
          endif
          if (lgrad2) then
            deriv2 =  - deriv2
            if (radsum.gt.0.0_dp) then
              rtrm2 =  - rtrm2
              rderiv =  - rderiv
            endif
          endif
        endif
        if (lsuttonc.and.imode.eq.1) then
!
!  Don't need to add terms for perfect region as rho  =  bulk rho
!
          if (.not.lMEAMden) then
            if (lorder12loc) then
              dscrholoc(i) = dscrholoc(i) + sctrm1*ocj
            else
              dscrholoc(i) = dscrholoc(i) + sctrm2*ocj
            endif
          endif
        endif
!***********************
!  Radial derivatives  *
!***********************
        if (lgrad1) then
          if (radi.gt.0.0_dp) then
            raderv(i) = raderv(i) + rtrm1
            if (lgrad2) then
              derv2(indrif,indri) = derv2(indrif,indri) + rtrm2
            endif
          endif
        endif
!**************************
!  Coordinate Derivatives *
!**************************
!
!  First derivatives
!
        if (lgrad1) then
          xdrv(i) = xdrv(i) - deriv*xcrd
          ydrv(i) = ydrv(i) - deriv*ycrd
          zdrv(i) = zdrv(i) - deriv*zcrd
          if (lgrad2.and.radi.gt.0.0_dp) then
            derv2(indrif,ix) = derv2(indrif,ix) - rderiv*xcrd
            derv2(indrif,iy) = derv2(indrif,iy) - rderiv*ycrd
            derv2(indrif,iz) = derv2(indrif,iz) - rderiv*zcrd
          endif
        endif
!
!  Second derivatives
!
        if (lgrad2) then
          derv2(jx,ix) = derv2(jx,ix) - deriv2*rpd1
          derv2(jx,iy) = derv2(jx,iy) - deriv2*rpd6
          derv2(jx,iz) = derv2(jx,iz) - deriv2*rpd5
          derv2(jy,iy) = derv2(jy,iy) - deriv2*rpd2
          derv2(jy,iz) = derv2(jy,iz) - deriv2*rpd4
          derv2(jz,iz) = derv2(jz,iz) - deriv2*rpd3
          derv2(jx,ix) = derv2(jx,ix) - deriv
          derv2(jy,iy) = derv2(jy,iy) - deriv
          derv2(jz,iz) = derv2(jz,iz) - deriv
        endif
!
!  Symmetrise second derivative matrix
!
        if (lgrad2) then
          derv2(jy,ix) = derv2(jx,iy)
          derv2(jz,ix) = derv2(jx,iz)
          derv2(jz,iy) = derv2(jy,iz)
        endif
      endif
    enddo
!************************
!  Loop over region 2a  *
!************************
    lregion1j = .false.
    rmiddle2 = 0.0_dp
    j = 0
    do while (j.lt.npreg2.and.rmiddle2.lt.rcheck2)
      j = j + 1
      xcrd = xr2a(j) - xal
      ycrd = yr2a(j) - yal
      zcrd = zr2a(j) - zal
      if (.not.lmol) then
!
!  If not a molecular calc and component exceeds maximum
!  cut - off then there is nothing to evaluate
!
        if (abs(xcrd).gt.cmax) goto 10
        if (abs(ycrd).gt.cmax) goto 10
        if (abs(zcrd).gt.cmax) goto 10
      endif
      natj = nr2a(j)
      ntypj = ntr2a(j)
!
!  Arrange nat1 to be lower than nat2. If they are equal
!  then arrange that ntyp1 is lower than ntyp2
!
      if (nati.eq.natj) then
        lorder12loc = .true.
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
        lorder12loc = .true.
        nat1 = nati
        nat2 = natj
        ntyp1 = ntypi
        ntyp2 = ntypj
      else
        lorder12loc = .false.
        nat1 = natj
        nat2 = nati
        ntyp1 = ntypj
        ntyp2 = ntypi
      endif
      qlj = qr2a(j)
      ocj = or2a(j)
      ofct = oci*ocj
      factor = qli*qlj*ofct*angstoev
      npsj = nps(j)
      iar = nsft + nrelat(npsj)
      if (lbsmat(iar)) then
        radj = radcfg(iar)
      else
        radj = 0.0_dp
      endif
      radsum = radi + radj
      lbonded = .false.
      l2bonds = .false.
      l3bonds = .false.
!
!  Distance check relative to centre of region 1
!
      xdiffc = xr2a(j) - xdc
      ydiffc = yr2a(j) - ydc
      zdiffc = zr2a(j) - zdc
      rmiddle2 = xdiffc*xdiffc + ydiffc*ydiffc + zdiffc*zdiffc
!
!  Possible core - shell flag
!
      lcspair = (abs(nat1-nat2).eq.maxele.or.(oci+ocj).lt.1.0001_dp)
      if (lcspair) then
        ospfct = oci
      else
        ospfct = ofct
      endif
      if (lmol) then
!
!  Molecule handling
!
        nmj = nmr2a(j)
        indmj = nmir2a(j)
        lmolok = (nmi.eq.nmj.and.nmi.gt.0)
      else
        lmolok = .false.
      endif
!
!  Locate potential number
!  Check whether potential requires specific types
!
      rp = 0.0_dp
      npots = 0
      c6tot = 0.0_dp
      do n = 1,npote
        if (lmatch(nat1,ntyp1,nspec1(n),nptyp1(n),.true.)) then
          if (lmatch(nat2,ntyp2,nspec2(n),nptyp2(n),.true.)) then
            npots = npots + 1
            npotl(npots) = n
            if (nptype(n).eq.8.or.nptype(n).eq.33) then
              if (cuts.gt.rp) rp = cuts
            elseif (lc6) then
              if (nptype(n).eq.1.or.nptype(n).eq.7) then
                c6tot = c6tot + twopot(3,n)
                if (repcut(n).gt.rp) rp = repcut(n)
              elseif (nptype(n).eq.2.or.nptype(n).eq.21) then
                c6tot = c6tot + twopot(2,n)
                if (repcut(n).gt.rp) rp = repcut(n)
              elseif (nptype(n).eq.57) then
                c6tot = c6tot + twopot(3,n)*twopot(4,n)*twopot(5,n)**6
                if (repcut(n).gt.rp) rp = repcut(n)
              else
                if (rpot(n).gt.rp) rp = rpot(n)
              endif
            else
              if (rpot(n).gt.rp) rp = rpot(n)
            endif
          endif
        endif
      enddo
!
!  Generate looping indices
!
      cut2r = rp*rp
      if (cut2r.gt.cut2p) cut2r = cut2p
      cut2 = cut2r
      if (cut2e.gt.cut2.and.lewald) cut2 = cut2e
      if (cut2w.gt.cut2.and.lwolf) cut2 = cut2w
      rp = sqrt(cut2)
!
!  Generate distance squared
!
      r = xcrd*xcrd + ycrd*ycrd + zcrd*zcrd
!
!  Molecule - check index
!
      if (lmolok.and.(r.gt.cut2s.or..not.lcspair)) then
        ind = indmj - indmi
        lptrmol = (ind.eq.0)
        if (.not.lptrmol) then
          call mindtoijk(indmj,jxx,jyy,jzz)
          call mindtoijk(indmi,ixx,iyy,izz)
          jxx = jxx - ixx
          jyy = jyy - iyy
          jzz = jzz - izz
          call samemol(lptrmol,nmi,jxx,jyy,jzz,0_i4,0_i4,0_i4)
        endif
        if (lptrmol) then
          if (npsi.gt.0) then
            call mindtoijk(indmj,jxx,jyy,jzz)
            call mindtoijk(indmi,ixx,iyy,izz)
            jxx = jxx - ixx
            jyy = jyy - iyy
            jzz = jzz - izz
            indmmj = nmolind(npsj)
            indmmi = nmolind(npsi)
            call mindtoijk(indmmj,jxx2,jyy2,jzz2)
            call mindtoijk(indmmi,ixx2,iyy2,izz2)
            jxx = jxx + jxx2 - ixx2
            jyy = jyy + jyy2 - iyy2
            jzz = jzz + jzz2 - izz2
            call bonded3(lbonded,l2bonds,l3bonds,nbtypeij,nbtypeij2,npsi,npsj,jxx,jyy,jzz)
          else
            lbonded   = .false.
            l2bonds   = .false.
            l3bonds   = .false.
            nbtypeij  = 0
            nbtypeij2 = 0
          endif
        else
          lbonded   = .false.
          l2bonds   = .false.
          l3bonds   = .false.
          nbtypeij  = 0
          nbtypeij2 = 0
        endif
      else
        lptrmol   = .false.
        lbonded   = .false.
        l2bonds   = .false.
        l3bonds   = .false.
        nbtypeij  = 0
        nbtypeij2 = 0
      endif
      if (r.le.cut2.or.lptrmol) then
!
!  Store vector
!
        dist = sqrt(r)
        call twobody1(eatom1,ereal1,ec61,lgrad1,lgrad2,.false.,1_i4,1_i4,dist,xcrd,ycrd,zcrd,d0i,d0j, &
                      deriv,deriv2,deriv3,derive0,derive,derive2,derive3,rderiv,npots,npotl,cut2r,cut2q,cut2s, &
                      lptrmol,nmolonly,factor,ofct,ospfct,radsum,rtrm1,rtrm2,rtrm3,rtrm32,sctrm1,sctrm2, &
                      qli,qlj,lcspair,.true.,lregion1j,lbonded,l2bonds,l3bonds,nbtypeij,nbtypeij2, &
                      lregion1j,.false.,d1i,d1j,d2i2,d2ij,d2j2)
        if (lgrad1) then
          call rpaird(i,j,indri,indrif,ix,iy,iz,jx,jy,jz,lgrad2,xcrd,ycrd,zcrd,deriv,deriv2,rderiv,rtrm1,rtrm2,radi,radj)
        endif
        if (lsuttonc) then
!
!  Don't need to add terms for perfect region as rho = bulk rho
!
          if (.not.lMEAMden) then
            if (lorder12loc) then
              if (imode.eq.1) then
                dscrholoc(i) = dscrholoc(i) + sctrm1*ocj
                dscrhor2dloc(j) = dscrhor2dloc(j) + sctrm2*oci
              else
                dscrhor2ploc(j) = dscrhor2ploc(j) - sctrm2*oci
              endif
            else
              if (imode.eq.1) then
                dscrholoc(i) = dscrholoc(i) + sctrm2*ocj
                dscrhor2dloc(j) = dscrhor2dloc(j) + sctrm1*oci
              else
                dscrhor2ploc(j) = dscrhor2ploc(j) - sctrm1*oci
              endif
            endif
          endif
        endif
        if (lgrad2) then
!
!  Symmetrise second derivative matrix
!
          derv2(jy,ix) = derv2(jx,iy)
          derv2(jz,ix) = derv2(jx,iz)
          derv2(jz,iy) = derv2(jy,iz)
        endif
      endif
10    continue
    enddo
  enddo
!
999 continue
!
!  Globalise energies if in parallel
!
  if (nprocs.gt.1) then
    esum1(1) = eatom1
    esum1(2) = ereal1
    esum1(3) = erecip1
    esum1(4) = ec61
    call sumall(esum1,esum2,4_i4,"real12a2d","esum")
    eatom  = eatom  + esum2(1)
    ereal  = ereal  + esum2(2)
    erecip = erecip + esum2(3)
    ec6    = ec6    + esum2(4)
!
    allocate(dtmp(max(nreg1,npreg2)),stat=status)
    if (status/=0) call outofmemory('real12a2','dtmp')
    if (imode.eq.1) then
      call sumall(dscrholoc,dtmp,nreg1,"real12a2","dscrho")
      dscrho(1,1:nreg1) = dscrho(1,1:nreg1) + dtmp(1:nreg1)
      call sumall(dscrhor2dloc,dtmp,npreg2,"real12a2","dscrhor2d")
      dscrhor2d(1,1:npreg2) = dscrhor2d(1,1:npreg2) + dtmp(1:npreg2)
    else
      call sumall(dscrhor2ploc,dtmp,npreg2,"real12a2","dscrhor2p")
      dscrhor2p(1,1:npreg2) = dscrhor2p(1,1:npreg2) + dtmp(1:npreg2)
    endif
    deallocate(dtmp,stat=status)
    if (status/=0) call deallocate_error('real12a2','dtmp')
  else
    eatom  = eatom  + eatom1
    ereal  = ereal  + ereal1
    erecip = erecip + erecip1
    ec6    = ec6    + ec61
!
    if (imode.eq.1) then
      dscrho(1,1:nreg1) = dscrho(1,1:nreg1) + dscrholoc(1:nreg1)
      dscrhor2d(1,1:npreg2) = dscrhor2d(1,1:npreg2) + dscrhor2dloc(1:npreg2)
    else
      dscrhor2p(1,1:npreg2) = dscrhor2p(1,1:npreg2) + dscrhor2ploc(1:npreg2)
    endif
  endif
!
!  Free local memory
!
  if (imode.eq.1) then
    deallocate(dscrhor2dloc,stat=status)
    if (status/=0) call deallocate_error('real12a2','dscrhor2dloc')
    deallocate(dscrholoc,stat=status)
    if (status/=0) call deallocate_error('real12a2','dscrholoc')
  else
    deallocate(dscrhor2ploc,stat=status)
    if (status/=0) call deallocate_error('real12a2','dscrhor2ploc')
  endif
  deallocate(npotl,stat=status)
  if (status/=0) call deallocate_error('real12a2d','npotl')
!
!  Timing
!
  time2 = g_cpu_time()
  treg1 = treg1 + time2 - time1
#ifdef TRACE
  call trace_out('real12a2d')
#endif
!
  return
  end
