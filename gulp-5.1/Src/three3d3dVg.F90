  subroutine three3d3dVg(i,j,nati,ntypi,natj,ntypj,d3s, &
    xc1,yc1,zc1,x21,y21,z21,nthbk)
!
!  Subroutine for third derivatives of three-body energy for i-j
!  pair with respect to the volume. Called from realrecip3d3dV. 
!  On entry we already have a possible valid i-j pair. 
!  Algorithm modified - less checking now done while searching for 
!  triads and more when a possible triad has been found. 
!  This makes the checking for valid combinations simpler and more 
!  robust (in the case of multiple pivots per triad) at the cost 
!  of more cpu time. Gamma point version.
!
!   1/18 Created from three3d3dV
!   2/18 Trace added
!
!  All d33 type arrays are stored as :
!    1 = xxx 10 = xxy 19 = xxz |
!    2 = yxx 11 = yxy 20 = yxz |
!    3 = zxx 12 = zxy 21 = zxz |for 1-3
!    4 = xyx 13 = xyy 22 = xyz |repeat in 28-54 for 2-3
!    5 = yyx 14 = yyy 23 = yyz |
!    6 = zyx 15 = zyy 24 = zyz |
!    7 = xzx 16 = xzy 25 = xzz |
!    8 = yzx 17 = yzy 26 = yzz |
!    9 = zzx 18 = zzy 27 = zzz |
!
!  ipivot = pointer to i,j,k according to which is the central atom
!           for an asymmetric potential
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
  use g_constants
  use current
  use feworkspace
  use ksample
  use m_three
  use symmetry
  use times
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4)                               :: i
  integer(i4)                               :: j
  integer(i4)                               :: nati
  integer(i4)                               :: natj
  integer(i4)                               :: nthbk
  integer(i4)                               :: ntypi
  integer(i4)                               :: ntypj
  real(dp)                                  :: d3s(3,3,3)
  real(dp)                                  :: x21
  real(dp)                                  :: y21
  real(dp)                                  :: z21
  real(dp)                                  :: xc1
  real(dp)                                  :: yc1
  real(dp)                                  :: zc1
!
!  Local variables
!
  integer(i4)                               :: ii
  integer(i4)                               :: ii1
  integer(i4)                               :: ii2
  integer(i4)                               :: ii3
  integer(i4)                               :: iii
  integer(i4)                               :: imax
  integer(i4)                               :: ind
  integer(i4)                               :: ind2
  integer(i4)                               :: ipivot
  integer(i4)                               :: ipivottriad(3)
  integer(i4)                               :: jj
  integer(i4)                               :: jjmin
  integer(i4)                               :: jkorder(3)
  integer(i4)                               :: jmax
  integer(i4)                               :: k
  integer(i4)                               :: kd3
  integer(i4)                               :: kk
  integer(i4)                               :: kmax
  integer(i4),                         save :: maxvector = 100
  integer(i4)                               :: n
  integer(i4)                               :: n3
  integer(i4)                               :: n3ty
  integer(i4)                               :: natk
  integer(i4)                               :: nmid
  integer(i4)                               :: nsame
  integer(i4)                               :: nt1
  integer(i4)                               :: nt2
  integer(i4)                               :: nt3
  integer(i4)                               :: ntmp
  integer(i4)                               :: ntriad
  integer(i4)                               :: ntyp1
  integer(i4)                               :: ntyp2
  integer(i4)                               :: ntyp3
  integer(i4)                               :: ntypk
  integer(i4)                               :: nvector
  integer(i4)                               :: status
  logical                                   :: bonded2donor
  logical                                   :: bonded2donorJK
  logical                                   :: lkfound
  logical                                   :: lmatch
  logical                                   :: lsymijk
  logical                                   :: lvalid
  logical                                   :: lvalidij
  real(dp)                                  :: ang
  real(dp)                                  :: g_cpu_time
  real(dp)                                  :: cut
  real(dp)                                  :: d0i
  real(dp)                                  :: d0j
  real(dp)                                  :: d0k
  real(dp)                                  :: d1q(3,3)
  real(dp)                                  :: d2q(6)
  real(dp)                                  :: d3l1(3,3,3)
  real(dp)                                  :: d3l2(27)
  real(dp)                                  :: d3l3(27)
  real(dp)                                  :: d3sl(3,3,3)
  real(dp)                                  :: dot
  real(dp)                                  :: e1d1
  real(dp)                                  :: e1d2
  real(dp)                                  :: e1d3
  real(dp)                                  :: e2d(6)
  real(dp)                                  :: e3d(10)
  real(dp)                                  :: ethb
  real(dp)                                  :: ocij
  real(dp)                                  :: ock
  real(dp)                                  :: ofct
  real(dp)                                  :: qli
  real(dp)                                  :: qlj
  real(dp)                                  :: qlk
  real(dp)                                  :: qjk(3)
  real(dp)                                  :: r12
  real(dp)                                  :: r122
  real(dp)                                  :: r13
  real(dp)                                  :: r132
  real(dp)                                  :: r23
  real(dp)                                  :: r232
  real(dp)                                  :: rho1
  real(dp)                                  :: rho2
  real(dp)                                  :: rho3
  real(dp)                                  :: rho4
  real(dp)                                  :: rho5
  real(dp)                                  :: rk32
  real(dp)                                  :: rk33
  real(dp)                                  :: rk34
  real(dp)                                  :: rkthb
  real(dp)                                  :: rkthb3
  real(dp)                                  :: rkthb4
  real(dp)                                  :: ro1
  real(dp)                                  :: ro2
  real(dp)                                  :: ro3
  real(dp)                                  :: ro4
  real(dp)                                  :: ro5
  real(dp)                                  :: symfct
  real(dp)                                  :: the0
  real(dp)                                  :: time1
  real(dp)                                  :: time2
  real(dp)                                  :: tr1
  real(dp)                                  :: tr2
  real(dp)                                  :: tr3
  real(dp)                                  :: tr11
  real(dp)                                  :: tr21
  real(dp)                                  :: tr31
  real(dp)                                  :: trm1
  real(dp)                                  :: trm2
  real(dp)                                  :: trm3
  real(dp)                                  :: trm4
  real(dp)                                  :: trm5
  real(dp)                                  :: trm6
  real(dp)                                  :: trm7
  real(dp)                                  :: trm8
  real(dp)                                  :: trm9
  real(dp)                                  :: trm10
  real(dp)                                  :: trm11
  real(dp)                                  :: trm12
  real(dp)                                  :: trm13
  real(dp)                                  :: trm14
  real(dp)                                  :: trm15
  real(dp)                                  :: trm16
  real(dp)                                  :: trm17
  real(dp)                                  :: trm18
  real(dp)                                  :: trm19
  real(dp)                                  :: trm20
  real(dp)                                  :: trm21
  real(dp)                                  :: trmax
  real(dp)                                  :: ttr11
  real(dp)                                  :: ttr21
  real(dp)                                  :: ttr31
  real(dp)                                  :: v12(3)
  real(dp)                                  :: v13(3)
  real(dp)                                  :: v23(3)
  real(dp)                                  :: x23
  real(dp)                                  :: y23
  real(dp)                                  :: z23
  real(dp)                                  :: x31
  real(dp)                                  :: y31
  real(dp)                                  :: z31
  real(dp)                                  :: xt21
  real(dp)                                  :: yt21
  real(dp)                                  :: zt21
  real(dp)                                  :: xt31
  real(dp)                                  :: yt31
  real(dp)                                  :: zt31
  real(dp), dimension(:), allocatable, save :: xvec
  real(dp), dimension(:), allocatable, save :: yvec
  real(dp), dimension(:), allocatable, save :: zvec
#ifdef TRACE
  call trace_in('three3d3dVg')
#endif
!
  time1 = g_cpu_time()
!
!  Charges are just dummies here since third derivatives are not available for variable charges
!
  qli = 0.0_dp
  qlj = 0.0_dp
  qlk = 0.0_dp
!
!  Allocate local memory
!
  allocate(xvec(maxvector),stat=status)
  if (status/=0) call outofmemory('three3d3dV','xvec')
  allocate(yvec(maxvector),stat=status)
  if (status/=0) call outofmemory('three3d3dV','yvec')
  allocate(zvec(maxvector),stat=status)
  if (status/=0) call outofmemory('three3d3dV','zvec')
!
!  Calculate cartesian K vector components
!
  ocij = occuf(i)*occuf(j)
!
!  Loop over potentials
!
  potloop: do n = 1,nthb
    n3ty = nthrty(n)
    nt1 = ntspec1(n)
    nt2 = ntspec2(n)
    nt3 = ntspec3(n)
    ntyp1 = ntptyp1(n)
    ntyp2 = ntptyp2(n)
    ntyp3 = ntptyp3(n)
    tr11 = thr1(n)
    tr21 = thr2(n)
    tr31 = thr3(n)
    tr1 = tr11*tr11
    tr2 = tr21*tr21
    tr3 = tr31*tr31
    rkthb = thbk(n)
    rkthb3 = 0.0_dp
    rkthb4 = 0.0_dp
    ro1 = 0.0_dp
    ro2 = 0.0_dp
    ro3 = 0.0_dp
    ro4 = 0.0_dp
    ro5 = 0.0_dp
    if (n3ty.eq.2) then
      the0 = theta(n)*degtorad
      ro1 = thrho1(n)
      ro2 = thrho2(n)
      if (ro1.ne.0.0_dp) ro1 = 1.0_dp/ro1
      if (ro2.ne.0.0_dp) ro2 = 1.0_dp/ro2
      lsymijk = .false.
    elseif (n3ty.eq.1) then
      the0 = theta(n)*degtorad
      rkthb3 = thrho2(n)/6.0_dp
      rkthb4 = thrho1(n)/24.0_dp
      lsymijk = .false.
    elseif (n3ty.eq.3) then
      the0 = 0.0_dp
      lsymijk = .true.
    elseif (n3ty.eq.4) then
      the0 = 0.0_dp
      ro1 = theta(n)
      ro2 = thrho1(n)
      ro3 = thrho2(n)
      lsymijk = .true.
    elseif (n3ty.eq.5.or.n3ty.eq.25) then
      the0 = cos(theta(n)*degtorad)
      ro1 = thrho1(n)
      ro2 = thrho2(n)
      ro3 = thrho3(n)
      lsymijk = .false.
    elseif (n3ty.eq.6) then
      the0 = 0.0_dp
      ro1 = thrho1(n)
      ro2 = thrho2(n)
      lsymijk = .false.
    elseif (n3ty.eq.7) then
      the0 = theta(n)
      lsymijk = .false.
    elseif (n3ty.eq.8) then
      the0 = theta(n)*degtorad
      ro1 = thrho1(n)
      ro2 = thrho2(n)
      if (ro1.ne.0.0_dp) ro1 = 1.0_dp/ro1
      if (ro2.ne.0.0_dp) ro2 = 1.0_dp/ro2
      the0 = the0-pi
      the0 = the0*the0
      rkthb = 0.25_dp*rkthb/the0
      lsymijk = .false.
    elseif (n3ty.eq.9) then
      the0 = cos(theta(n)*degtorad)
      rkthb3 = thrho2(n)/6.0_dp
      rkthb4 = thrho1(n)/24.0_dp
      lsymijk = .false.
    elseif (n3ty.eq.10) then
      the0 = theta(n)
      ro1 = thrho1(n)
      ro2 = thrho2(n)
      ro3 = thrho3(n)
      lsymijk = .true.
    elseif (n3ty.eq.11) then
      rkthb3 = thrho3(n)
      ro1 = thrho1(n)
      ro2 = thrho2(n)
      the0 = theta(n)*degtorad
      lsymijk = .false.
    elseif (n3ty.eq.12) then
      the0 = theta(n)
      rkthb3 = nint(thrho1(n))
      lsymijk = .false.
    elseif (n3ty.eq.13) then
      the0 = theta(n)
      rkthb3 = thrho3(n)
      ro1 = thrho1(n)
      ro2 = thrho2(n)  
      lsymijk = .false.
    elseif (n3ty.eq.14) then
      the0 = cos(theta(n)*degtorad)
      ro1 = thrho1(n)
      ro2 = thrho2(n)
      ro3 = thrho3(n)
      lsymijk = .false.
    elseif (n3ty.eq.15) then
      rkthb3 = theta(n)
      ro1 = thrho1(n)
      ro2 = thrho2(n)
      ro3 = thrho3(n)
      lsymijk = .false.
    elseif (n3ty.eq.16) then
      the0 = theta(n)
      ro1 = thrho1(n)
      ro2 = thrho2(n)
      lsymijk = .false.
    elseif (n3ty.eq.17) then
      the0 = theta(n)*degtorad
      rkthb4 = 1.0_dp/(2.0_dp*sin(the0))**2
      rkthb3 = - 4.0_dp*rkthb4*cos(the0)
      the0 = rkthb4*(2.0_dp*cos(the0)**2 + 1.0_dp)
      lsymijk = .false.
    elseif (n3ty.eq.18) then
      rkthb3 = thrho3(n)
      ro1 = thrho1(n)
      ro2 = thrho2(n)
      the0 = cos(theta(n)*degtorad)
      lsymijk = .false.
    elseif (n3ty.eq.19) then
      lsymijk = .false.
    elseif (n3ty.eq.20) then
      the0 = 0.0_dp
      ro1 = theta(n)
      ro2 = thrho2(n)
      ro4 = thrho1(n)
      ro5 = thrho3(n)
      lsymijk = .false.
    elseif (n3ty.eq.21) then
      the0 = theta(n)
      lsymijk = .false.
    elseif (n3ty.eq.22) then
      ro1 = thrho1(n)
      ro2 = thrho2(n)
      the0 = theta(n)*degtorad
      lsymijk = .false.
    elseif (n3ty.eq.23) then
      ro1 = thrho1(n)
      ro2 = thrho2(n)
      the0 = theta(n)
      lsymijk = .false.
    elseif (n3ty.eq.24) then
      the0 = theta(n)*degtorad
      lsymijk = .false.
    endif
!
!  Work out symmetry factor for symmetric potentials
!
    if (lsymijk) then
      nsame = 0
      if (nt1.eq.nt2.and.ntyp1.eq.ntyp2) nsame = nsame + 1
      if (nt1.eq.nt3.and.ntyp1.eq.ntyp3) nsame = nsame + 1
      if (nsame.eq.0) then
        symfct = 1.0_dp
      elseif (nsame.eq.1) then
        symfct = 0.5_dp
      elseif (nsame.eq.2) then
        symfct = 1.0_dp/3.0_dp
      endif
    endif
!
!  Find maximum cut-off
!
    trmax = max(tr1,tr2,tr3)
!
!  Create lattice vectors
!
    cut = trmax
    cut = sqrt(cut)
    call rtlist(nvector,cut,xvec,yvec,zvec,imax,jmax,kmax,nmid,maxvector)
    if (nvector.gt.maxvector) then
!
!  Too many vectors
!
      deallocate(zvec,stat=status)
      if (status/=0) call deallocate_error('three3d3dV','zvec')
      deallocate(yvec,stat=status)
      if (status/=0) call deallocate_error('three3d3dV','yvec')
      deallocate(xvec,stat=status)
      if (status/=0) call deallocate_error('three3d3dV','xvec')
      maxvector = nint(1.1*nvector)
      allocate(xvec(maxvector),stat=status)
      if (status/=0) call outofmemory('three3d3dV','xvec')
      allocate(yvec(maxvector),stat=status)
      if (status/=0) call outofmemory('three3d3dV','yvec')
      allocate(zvec(maxvector),stat=status)
      if (status/=0) call outofmemory('three3d3dV','zvec')
      call rtlist(nvector,cut,xvec,yvec,zvec,imax,jmax,kmax,nmid,maxvector)
    endif
!*****************************
!  Are i and j valid types?  *
!*****************************
    lvalidij = .false.
    if (lmatch(nati,ntypi,nt1,ntyp1,.false.)) then
      if (lmatch(natj,ntypj,nt2,ntyp2,.false.)) then
        lvalidij = .true.
      elseif (lmatch(natj,ntypj,nt3,ntyp3,.false.)) then
        lvalidij = .true.
        ntmp = nt3
        nt3 = nt2
        nt2 = ntmp
        ntmp = ntyp3
        ntyp3 = ntyp2
        ntyp2 = ntmp
      endif
    elseif (lmatch(nati,ntypi,nt2,ntyp2,.false.)) then
      if (lmatch(natj,ntypj,nt1,ntyp1,.false.)) then
        lvalidij = .true.
      elseif (lmatch(natj,ntypj,nt3,ntyp3,.false.)) then
        lvalidij = .true.
        ntmp = nt3
        nt3 = nt1
        nt1 = ntmp
        ntmp = ntyp3
        ntyp3 = ntyp1
        ntyp1 = ntmp
      endif
    elseif (lmatch(nati,ntypi,nt3,ntyp3,.true.)) then
      if (lmatch(natj,ntypj,nt2,ntyp2,.true.)) then
        lvalidij = .true.
        ntmp = nt3
        nt3 = nt1
        nt1 = ntmp
        ntmp = ntyp3
        ntyp3 = ntyp1
        ntyp1 = ntmp
      elseif (lmatch(natj,ntypj,nt1,ntyp1,.true.)) then
        lvalidij = .true.
        ntmp = nt2
        nt2 = nt3
        nt3 = ntmp
        ntmp = ntyp2
        ntyp2 = ntyp3
        ntyp3 = ntmp
      endif
    endif
!
!  Not a valid three body potential based on types  = > skip to next potential
!
    if (.not.lvalidij) cycle potloop
!
!  Dreiding option handling
!
    if (ltdreiding(n)) then
      if (.not.bonded2donor(i)) cycle potloop
    endif
!
!  Check r12 is OK
!  Loop over cell vectors
!
    iiloop: do ii = 1,nvector
      r122 = (xvec(ii) + x21)**2 + (yvec(ii) + y21)**2 + (zvec(ii) + z21)**2
      if (r122.lt.1.0d-10.or.r122.gt.trmax) cycle iiloop
      r12 = sqrt(r122)
!*******************************
!  Inner loop over third site  *
!*******************************
      kloop: do k = 1,numat
        natk = nat(k)
        ntypk = nftype(k)
!
!  Check k is allowed for n, and not equivalent to j
!
        if (.not.lmatch(natk,ntypk,nt3,ntyp3,.true.)) cycle kloop
!
!  Dreiding option handling
!
        if (ltdreiding(n)) then
          if (.not.bonded2donorJK(i,j,k)) cycle kloop
        endif
!
        ock = occuf(k)
        x31 = xclat(k) - xc1
        y31 = yclat(k) - yc1
        z31 = zclat(k) - zc1
        if (j.eq.k) then
          jjmin = ii + 1
        else
          jjmin = 1
        endif
!
!  Check r13 is OK
!  Loop over cell vectors
!
        jjloop: do jj = jjmin,nvector
          r132 = (xvec(jj) + x31)**2 + (yvec(jj) + y31)**2 + (zvec(jj) + z31)**2
          if (r132.lt.1.0d-10.or.r132.gt.trmax) cycle jjloop
!
!  Check r23 is OK
!
          xt31 = x31 + xvec(jj)
          yt31 = y31 + yvec(jj)
          zt31 = z31 + zvec(jj)
          xt21 = x21 + xvec(ii)
          yt21 = y21 + yvec(ii)
          zt21 = z21 + zvec(ii)
          x23 = xt31 - xt21
          y23 = yt31 - yt21
          z23 = zt31 - zt21
          r232 = x23**2 + y23**2 + z23**2
          if (r232.lt.1d-10.or.r232.gt.trmax) cycle jjloop
!
          r13 = sqrt(r132)
          r23 = sqrt(r232)
!*********************************************
!  Possible valid triad - now do full checks *
!*********************************************
          ntriad  =  0
!
!  Pivot  =  i
!
          call validtriad(i,j,k,nati,ntypi,natj,ntypj,natk,ntypk,r12,r13,r23,n, &
            ndim,ii,jj,imax,jmax,kmax,lvalid)
          if (lvalid) then
            ntriad = ntriad + 1
            ipivottriad(ntriad) = 1
            jkorder(ntriad) = 1
            qjk(ntriad) = qf(j)*qf(k)
          else
            call validtriad(i,k,j,nati,ntypi,natk,ntypk,natj,ntypj,r13,r12,r23,n, &
              ndim,jj,ii,imax,jmax,kmax,lvalid)
            if (lvalid) then
              ntriad = ntriad + 1
              ipivottriad(ntriad) = 1
              jkorder(ntriad) = 2
              qjk(ntriad) = qf(j)*qf(k)
            endif
          endif
!
!  Pivot  =  j
!
          call validtriad(j,i,k,natj,ntypj,nati,ntypi,natk,ntypk,r12,r23,r13,n, &
            ndim,nvector-ii+1_i4,nmid-ii+jj,imax,jmax,kmax,lvalid)
          if (lvalid) then
            ntriad = ntriad + 1
            ipivottriad(ntriad) = 2
            jkorder(ntriad) = 1
            qjk(ntriad) = qf(i)*qf(k)
          else
            call validtriad(j,k,i,natj,ntypj,natk,ntypk,nati,ntypi,r23,r12,r13,n, &
              ndim,nmid-ii+jj,nvector-ii+1_i4,imax,jmax,kmax,lvalid)
            if (lvalid) then
              ntriad = ntriad + 1
              ipivottriad(ntriad) = 2
              jkorder(ntriad) = 2
              qjk(ntriad) = qf(i)*qf(k)
            endif
          endif
!
!  Pivot  =  k
!
          call validtriad(k,i,j,natk,ntypk,nati,ntypi,natj,ntypj,r13,r23,r12,n, &
            ndim,nvector-jj+1_i4,nmid-jj+ii,imax,jmax,kmax,lvalid)
          if (lvalid) then
            ntriad = ntriad + 1
            ipivottriad(ntriad) = 3
            jkorder(ntriad) = 1
            qjk(ntriad) = qf(i)*qf(j)
          else
            call validtriad(k,j,i,natk,ntypk,natj,ntypj,nati,ntypi,r23,r13,r12,n, &
              ndim,nmid-jj+ii,nvector-jj+1_i4,imax,jmax,kmax,lvalid)
            if (lvalid) then
              ntriad = ntriad + 1
              ipivottriad(ntriad) = 3
              jkorder(ntriad) = 2
              qjk(ntriad) = qf(i)*qf(j)
            endif
          endif
!
          if (ntriad.gt.0) then
!**************************
!  Valid three-body term  *
!**************************
!
!  Is k already on the list?
!
            lkfound = .false.
            iii = 0
            do while (iii.lt.nthbk.and..not.lkfound)
              iii = iii + 1
              lkfound = (nptrmanyk(iii).eq.k)
            enddo
!
!  If k is not already on the list, initialise d33 for strains
!
!  If the number of terms is too large, then carry on in order to
!  find out number of terms for the error message.
!
            if (.not.lkfound) then
              nthbk = nthbk + 1
              if (nthbk.gt.maxmany) then
                maxmany = nthbk + 40
                call changemaxmany
              endif
              nptrmanyk(nthbk) = k
              kd3 = nthbk
              do kk = 1,108
                d33s(kk,kd3)  = 0.0_dp
              enddo
            else
              kd3 = iii
            endif
!***************************
!  Loop over valid triads  *
!***************************
            do n3 = 1,ntriad
              ipivot = ipivottriad(n3)
!
!  Set potential parameters that depend on atom order
!
              if (ipivot.eq.1) then
                if (jkorder(n3).eq.1) then
                  rho1 = ro1
                  rho2 = ro2
                  rho3 = ro3
                  rho4 = ro4
                  rho5 = ro5
                  ttr11 = tr11
                  ttr21 = tr21
                  ttr31 = tr31
                else
                  rho1 = ro2
                  rho2 = ro1
                  rho3 = ro3
                  rho4 = ro5
                  rho5 = ro4
                  ttr11 = tr21
                  ttr21 = tr11
                  ttr31 = tr31
                endif
              elseif (ipivot.eq.2) then
                if (jkorder(n3).eq.1) then
                  rho1 = ro1
                  rho2 = ro3
                  rho3 = ro2
                  rho4 = ro4
                  rho5 = ro5
                  ttr11 = tr11
                  ttr21 = tr31
                  ttr31 = tr21
                else
                  rho1 = ro2
                  rho2 = ro3
                  rho3 = ro1
                  rho4 = ro5
                  rho5 = ro4
                  ttr11 = tr21
                  ttr21 = tr31
                  ttr31 = tr11
                endif
              else
                if (jkorder(n3).eq.1) then
                  rho1 = ro3
                  rho2 = ro1
                  rho3 = ro2
                  rho4 = ro4
                  rho5 = ro5
                  ttr11 = tr31
                  ttr21 = tr11
                  ttr31 = tr21
                else
                  rho1 = ro3
                  rho2 = ro2
                  rho3 = ro1
                  rho4 = ro5
                  rho5 = ro4
                  ttr11 = tr31
                  ttr21 = tr21
                  ttr31 = tr11
                endif
              endif
!
!  Calculate theta / cos(theta)
!
              if (n3ty.ne.3.and.n3ty.ne.4.and.n3ty.ne.6.and.n3ty.ne.7.and.n3ty.ne.19) then
                if (ipivot.eq.1) then
                  dot = xt21*xt31 + yt21*yt31 + zt21*zt31
                  dot = dot/(r12*r13)
                elseif (ipivot.eq.2) then
                  dot = xt21*x23 + yt21*y23 + zt21*z23
                  dot = - dot/(r12*r23)
                elseif (ipivot.eq.3) then
                  dot = x23*xt31 + y23*yt31 + z23*zt31
                  dot = dot/(r23*r13)
                endif
                if (abs(dot).gt.0.999999999999_dp) dot = sign(1.0_dp,dot)
                ang = acos(dot)
              else
                dot = 0.0_dp
                ang = 0.0_dp
              endif
              ofct = ocij*ock
              if (lsymijk) ofct = ofct*symfct
              if (n3ty.eq.19) then
                rk32 = rkthb*ofct*qjk(n3)
              else
                rk32 = rkthb*ofct
              endif
              if (n3ty.eq.12.or.n3ty.eq.17) then
                rk33 = rkthb3
              else
                rk33 = rkthb3*ofct
              endif
              if (n3ty.eq.17) then
                rk34 = rkthb4
              else
                rk34 = rkthb4*ofct
              endif
              if (n3ty.eq.15) then
                rho1 = thrho1(n)
                rho2 = thrho2(n)
                rho3 = thrho3(n)
              elseif (n3ty.eq.16) then
                rho1 = thrho1(n)
                rho2 = thrho2(n)
              endif
!*****************************************************
!  Calculate derivatives with respect to potentials  *
!*****************************************************
              call threebody(ipivot,n3ty,r12,r13,r23,e1d1,e1d2,e1d3,ethb,e2d,e3d,ttr11,ttr21,ttr31, &
                             rho1,rho2,rho3,rho4,rho5,rk32,rk33,rk34,the0,ang,dot,.true.,.true.,.true., &
                             n,qli,qlj,qlk,d0i,d0j,d0k,d1q,d2q,lthetataper(n),thetatapermin(n), &
                             thetatapermax(n))
!***************************************
!  Calculate phased third derivatives  *
!***************************************
!
!  Set up vector arrays
!
              v12(1) = xt21
              v12(2) = yt21
              v12(3) = zt21
              v13(1) = xt31
              v13(2) = yt31
              v13(3) = zt31
              v23(1) = x23
              v23(2) = y23
              v23(3) = z23
!
!  Zero local arrays
!
              ind = 0
              do ii1 = 1,3
                do ii2 = 1,3
                  do ii3 = 1,3
                    ind = ind + 1
                    d3l1(ii3,ii2,ii1) = 0.0_dp
                    d3l2(ind) = 0.0_dp
                    d3l3(ind) = 0.0_dp
                  enddo
                enddo
              enddo
!
!  Derivatives of i-j with respect to r12
!
!  First, non delta terms
!
              do ii1 = 1,3
                do ii2 = 1,3
                  do ii3 = 1,3
                    trm1 = e3d(1)*v12(ii3)*v12(ii2)*v12(ii1)
                    trm2 = e3d(2)*v13(ii3)*v12(ii2)*v12(ii1)
                    trm3 =  - e3d(3)*v12(ii3)*v23(ii2)*v12(ii1)
                    trm4 =  - e3d(5)*v13(ii3)*v23(ii2)*v12(ii1)
                    d3l1(ii3,ii2,ii1) = d3l1(ii3,ii2,ii1) + trm1
                    d3l1(ii3,ii2,ii1) = d3l1(ii3,ii2,ii1) + trm2
                    d3l1(ii3,ii2,ii1) = d3l1(ii3,ii2,ii1) + trm3
                    d3l1(ii3,ii2,ii1) = d3l1(ii3,ii2,ii1) + trm4
                  enddo
                enddo
              enddo
!
!  Now the delta terms
!
              trm1 = 3.0_dp*e2d(1)*v12(1) + e2d(2)*v13(1) - e2d(3)*v23(1)
              trm2 = e2d(1)*v12(2) + e2d(2)*v13(2)
              trm3 = e2d(1)*v12(3) + e2d(2)*v13(3)
              trm4 = e2d(1)*v12(2) - e2d(3)*v23(2)
              trm5 = e2d(1)*v12(1)
              trm6 = e2d(1)*v12(3) - e2d(3)*v23(3)
              trm7 = e2d(1)*v12(1)
              trm8 = e2d(1)*v12(2)
              trm9 = e2d(1)*v12(1) - e2d(3)*v23(1)
              trm10 = e2d(1)*v12(1) + e2d(2)*v13(1)
              trm11 = 3.0_dp*e2d(1)*v12(2) + e2d(2)*v13(2) - e2d(3)*v23(2)
              trm12 = e2d(1)*v12(3) + e2d(2)*v13(3)
              trm13 = e2d(1)*v12(3) - e2d(3)*v23(3)
              trm14 = e2d(1)*v12(2)
              trm15 = e2d(1)*v12(3)
              trm16 = e2d(1)*v12(1) - e2d(3)*v23(1)
              trm17 = e2d(1)*v12(3)
              trm18 = e2d(1)*v12(2) - e2d(3)*v23(2)
              trm19 = e2d(1)*v12(1) + e2d(2)*v13(1)
              trm20 = e2d(1)*v12(2) + e2d(2)*v13(2)
              trm21 = 3.0_dp*e2d(1)*v12(3) + e2d(2)*v13(3) - e2d(3)*v23(3)
!
              d3l1(1,1,1) = d3l1(1,1,1) + trm1
              d3l1(2,1,1) = d3l1(2,1,1) + trm2
              d3l1(3,1,1) = d3l1(3,1,1) + trm3
              d3l1(1,2,1) = d3l1(1,2,1) + trm4
              d3l1(2,2,1) = d3l1(2,2,1) + trm5
              d3l1(1,3,1) = d3l1(1,3,1) + trm6
              d3l1(3,3,1) = d3l1(3,3,1) + trm7
              d3l1(1,1,2) = d3l1(1,1,2) + trm8
              d3l1(2,1,2) = d3l1(2,1,2) + trm9
              d3l1(1,2,2) = d3l1(1,2,2) + trm10
              d3l1(2,2,2) = d3l1(2,2,2) + trm11
              d3l1(3,2,2) = d3l1(3,2,2) + trm12
              d3l1(2,3,2) = d3l1(2,3,2) + trm13
              d3l1(3,3,2) = d3l1(3,3,2) + trm14
              d3l1(1,1,3) = d3l1(1,1,3) + trm15
              d3l1(3,1,3) = d3l1(3,1,3) + trm16
              d3l1(2,2,3) = d3l1(2,2,3) + trm17
              d3l1(3,2,3) = d3l1(3,2,3) + trm18
              d3l1(1,3,3) = d3l1(1,3,3) + trm19
              d3l1(2,3,3) = d3l1(2,3,3) + trm20
              d3l1(3,3,3) = d3l1(3,3,3) + trm21
!
!  Derivatives of i - j with respect to r13
!
!  First, non delta terms
!
              ind = 0
              do ii1 = 1,3
                do ii2 = 1,3
                  do ii3 = 1,3
                    ind = ind + 1
                    trm1 = e3d(2)*v12(ii3)*v12(ii2)*v13(ii1)
                    trm2 = e3d(4)*v13(ii3)*v12(ii2)*v13(ii1)
                    trm3 = - e3d(5)*v12(ii3)*v23(ii2)*v13(ii1)
                    trm4 = - e3d(8)*v13(ii3)*v23(ii2)*v13(ii1)
                    d3l2(ind) = d3l2(ind) + trm1
                    d3l2(ind) = d3l2(ind) + trm2
                    d3l2(ind) = d3l2(ind) + trm3
                    d3l2(ind) = d3l2(ind) + trm4
                  enddo
                enddo
              enddo
!
!  Now the delta terms
!
              trm1 = e2d(2)*v13(1) + e2d(2)*v12(1) - e2d(5)*v23(1)
              trm2 = e2d(2)*v12(2) - e2d(5)*v23(2)
              trm3 = e2d(2)*v13(1)
              trm4 = e2d(2)*v12(3) - e2d(5)*v23(3)
              trm5 = e2d(2)*v13(1)
              trm6 = e2d(2)*v13(2)
              trm7 = e2d(2)*v12(1) - e2d(5)*v23(1)
              trm8 = e2d(2)*v13(2) + e2d(2)*v12(2) - e2d(5)*v23(2)
              trm9 = e2d(2)*v12(3) - e2d(5)*v23(3)
              trm10 = e2d(2)*v13(2)
              trm11 = e2d(2)*v13(3)
              trm12 = e2d(2)*v12(1) - e2d(5)*v23(1)
              trm13 = e2d(2)*v13(3)
              trm14 = e2d(2)*v12(2) - e2d(5)*v23(2)
              trm15 = e2d(2)*v13(3) + e2d(2)*v12(3) - e2d(5)*v23(3)
!
              d3l2(1) = d3l2(1) + trm1
              d3l2(4) = d3l2(4) + trm2
              d3l2(5) = d3l2(5) + trm3
              d3l2(7) = d3l2(7) + trm4
              d3l2(9) = d3l2(9) + trm5
              d3l2(10) = d3l2(10) + trm6
              d3l2(11) = d3l2(11) + trm7
              d3l2(14) = d3l2(14) + trm8
              d3l2(17) = d3l2(17) + trm9
              d3l2(18) = d3l2(18) + trm10
              d3l2(19) = d3l2(19) + trm11
              d3l2(21) = d3l2(21) + trm12
              d3l2(23) = d3l2(23) + trm13
              d3l2(24) = d3l2(24) + trm14
              d3l2(27) = d3l2(27) + trm15
!
!  Derivatives of i - j with respect to r23
!
!  First, non delta terms
!
              ind = 27
              ind2 = 0
              do ii1 = 1,3
                do ii2 = 1,3
                  do ii3 = 1,3
                    ind = ind + 1
                    ind2 = ind2 + 1
                    trm1 = e3d(3)*v12(ii3)*v12(ii2)*v23(ii1)
                    trm2 = e3d(5)*v13(ii3)*v12(ii2)*v23(ii1)
                    trm3 =  - e3d(6)*v12(ii3)*v23(ii2)*v23(ii1)
                    trm4 =  - e3d(9)*v13(ii3)*v23(ii2)*v23(ii1)
                    d3l3(ind2) = d3l3(ind2) + trm1
                    d3l3(ind2) = d3l3(ind2) + trm2
                    d3l3(ind2) = d3l3(ind2) + trm3
                    d3l3(ind2) = d3l3(ind2) + trm4
                  enddo
                enddo
              enddo
!
!  Now the delta terms
!
              trm1 = e2d(3)*v23(1) - e2d(3)*v12(1) - e2d(5)*v13(1)
              trm2 =  - e2d(3)*v12(2) - e2d(5)*v13(2)
              trm3 =  - e2d(3)*v12(3) - e2d(5)*v13(3)
              trm4 = e2d(3)*v23(1)
              trm5 = e2d(3)*v23(1)
              trm6 = e2d(3)*v23(2)
              trm7 =  - e2d(3)*v12(1) - e2d(5)*v13(1)
              trm8 = e2d(3)*v23(2) - e2d(3)*v12(2) - e2d(5)*v13(2)
              trm9 =  - e2d(3)*v12(3) - e2d(5)*v13(3)
              trm10 = e2d(3)*v23(2)
              trm11 = e2d(3)*v23(3)
              trm12 = e2d(3)*v23(3)
              trm13 =  - e2d(3)*v12(1) - e2d(5)*v13(1)
              trm14 =  - e2d(3)*v12(2) - e2d(5)*v13(2)
              trm15 = e2d(3)*v23(3) - e2d(3)*v12(3) - e2d(5)*v13(3)
!
              d3l3(1) = d3l3(1) + trm1
              d3l3(2) = d3l3(2) + trm2
              d3l3(3) = d3l3(3) + trm3
              d3l3(5) = d3l3(5) + trm4
              d3l3(9) = d3l3(9) + trm5
              d3l3(10) = d3l3(10) + trm6
              d3l3(13) = d3l3(13) + trm7
              d3l3(14) = d3l3(14) + trm8
              d3l3(15) = d3l3(15) + trm9
              d3l3(18) = d3l3(18) + trm10
              d3l3(19) = d3l3(19) + trm11
              d3l3(23) = d3l3(23) + trm12
              d3l3(25) = d3l3(25) + trm13
              d3l3(26) = d3l3(26) + trm14
              d3l3(27) = d3l3(27) + trm15
!*************************
!  Strain contributions  *
!*************************
!
!  Atoms 1 - 2
!
              d3sl(1,1,1) = d3l1(1,1,1)*xt21
              d3sl(2,1,1) = d3l1(2,1,1)*xt21
              d3sl(3,1,1) = d3l1(3,1,1)*xt21
              d3sl(1,2,1) = d3l1(1,2,1)*xt21
              d3sl(2,2,1) = d3l1(2,2,1)*xt21
              d3sl(3,2,1) = d3l1(3,2,1)*xt21
              d3sl(1,3,1) = d3l1(1,3,1)*xt21
              d3sl(2,3,1) = d3l1(2,3,1)*xt21
              d3sl(3,3,1) = d3l1(3,3,1)*xt21
              d3sl(1,1,2) = d3l1(1,1,2)*yt21
              d3sl(2,1,2) = d3l1(2,1,2)*yt21
              d3sl(3,1,2) = d3l1(3,1,2)*yt21
              d3sl(1,2,2) = d3l1(1,2,2)*yt21
              d3sl(2,2,2) = d3l1(2,2,2)*yt21
              d3sl(3,2,2) = d3l1(3,2,2)*yt21
              d3sl(1,3,2) = d3l1(1,3,2)*yt21
              d3sl(2,3,2) = d3l1(2,3,2)*yt21
              d3sl(3,3,2) = d3l1(3,3,2)*yt21
              d3sl(1,1,3) = d3l1(1,1,3)*zt21
              d3sl(2,1,3) = d3l1(2,1,3)*zt21
              d3sl(3,1,3) = d3l1(3,1,3)*zt21
              d3sl(1,2,3) = d3l1(1,2,3)*zt21
              d3sl(2,2,3) = d3l1(2,2,3)*zt21
              d3sl(3,2,3) = d3l1(3,2,3)*zt21
              d3sl(1,3,3) = d3l1(1,3,3)*zt21
              d3sl(2,3,3) = d3l1(2,3,3)*zt21
              d3sl(3,3,3) = d3l1(3,3,3)*zt21
!
              do ii1 = 1,3
                do ii2 = 1,3
                  do ii3 = 1,3
                    d3s(ii3,ii2,ii1) = d3s(ii3,ii2,ii1) + d3sl(ii3,ii2,ii1)
                  enddo
                enddo
              enddo
!
!  Atoms 1 - 3
!
              d3sl(1,1,1) = d3l2(1)*xt31
              d3sl(2,1,1) = d3l2(2)*xt31
              d3sl(3,1,1) = d3l2(3)*xt31
              d3sl(1,2,1) = d3l2(4)*xt31
              d3sl(2,2,1) = d3l2(5)*xt31
              d3sl(3,2,1) = d3l2(6)*xt31
              d3sl(1,3,1) = d3l2(7)*xt31
              d3sl(2,3,1) = d3l2(8)*xt31
              d3sl(3,3,1) = d3l2(9)*xt31
              d3sl(1,1,2) = d3l2(10)*yt31
              d3sl(2,1,2) = d3l2(11)*yt31
              d3sl(3,1,2) = d3l2(12)*yt31
              d3sl(1,2,2) = d3l2(13)*yt31
              d3sl(2,2,2) = d3l2(14)*yt31
              d3sl(3,2,2) = d3l2(15)*yt31
              d3sl(1,3,2) = d3l2(16)*yt31
              d3sl(2,3,2) = d3l2(17)*yt31
              d3sl(3,3,2) = d3l2(18)*yt31
              d3sl(1,1,3) = d3l2(19)*zt31
              d3sl(2,1,3) = d3l2(20)*zt31
              d3sl(3,1,3) = d3l2(21)*zt31
              d3sl(1,2,3) = d3l2(22)*zt31
              d3sl(2,2,3) = d3l2(23)*zt31
              d3sl(3,2,3) = d3l2(24)*zt31
              d3sl(1,3,3) = d3l2(25)*zt31
              d3sl(2,3,3) = d3l2(26)*zt31
              d3sl(3,3,3) = d3l2(27)*zt31
!
              ind = 0
              do ii1 = 1,3
                do ii2 = 1,3
                  do ii3 = 1,3
                    ind = ind + 1
                    d33s(ind,kd3) = d33s(ind,kd3) + d3sl(ii3,ii2,ii1)
                  enddo
                enddo
              enddo
!
!  Atoms 2 - 3
!
              d3sl(1,1,1) = d3l3(1)*x23
              d3sl(2,1,1) = d3l3(2)*x23
              d3sl(3,1,1) = d3l3(3)*x23
              d3sl(1,2,1) = d3l3(4)*x23
              d3sl(2,2,1) = d3l3(5)*x23
              d3sl(3,2,1) = d3l3(6)*x23
              d3sl(1,3,1) = d3l3(7)*x23
              d3sl(2,3,1) = d3l3(8)*x23
              d3sl(3,3,1) = d3l3(9)*x23
              d3sl(1,1,2) = d3l3(10)*y23
              d3sl(2,1,2) = d3l3(11)*y23
              d3sl(3,1,2) = d3l3(12)*y23
              d3sl(1,2,2) = d3l3(13)*y23
              d3sl(2,2,2) = d3l3(14)*y23
              d3sl(3,2,2) = d3l3(15)*y23
              d3sl(1,3,2) = d3l3(16)*y23
              d3sl(2,3,2) = d3l3(17)*y23
              d3sl(3,3,2) = d3l3(18)*y23
              d3sl(1,1,3) = d3l3(19)*z23
              d3sl(2,1,3) = d3l3(20)*z23
              d3sl(3,1,3) = d3l3(21)*z23
              d3sl(1,2,3) = d3l3(22)*z23
              d3sl(2,2,3) = d3l3(23)*z23
              d3sl(3,2,3) = d3l3(24)*z23
              d3sl(1,3,3) = d3l3(25)*z23
              d3sl(2,3,3) = d3l3(26)*z23
              d3sl(3,3,3) = d3l3(27)*z23
!
              ind = 54
              do ii1 = 1,3
                do ii2 = 1,3
                  do ii3 = 1,3
                    ind = ind + 1
                    d33s(ind,kd3) = d33s(ind,kd3) + d3sl(ii3,ii2,ii1)
                  enddo
                enddo
              enddo
!
!  End of loop over valid triads
!
            enddo
          endif
!***********************
!  End of derivatives  *
!***********************
!
!  End of inner loops
!
        enddo jjloop
      enddo kloop
    enddo iiloop
!
!  End of outer loops
!
  enddo potloop
!
!  Free local memory
!
  deallocate(zvec,stat=status)
  if (status/=0) call deallocate_error('three3d3dV','zvec')
  deallocate(yvec,stat=status)
  if (status/=0) call deallocate_error('three3d3dV','yvec')
  deallocate(xvec,stat=status)
  if (status/=0) call deallocate_error('three3d3dV','xvec')
!
!  Timing
!
  time2 = g_cpu_time()
  tthree = tthree + time2 - time1
#ifdef TRACE
  call trace_out('three3d3dVg')
#endif
!
  return
  end
