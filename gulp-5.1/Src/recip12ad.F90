  subroutine recip12ad(erecip,ec6,lgrad1,lgrad2,imode)
!
!  Calculate reciprocal space contribution to region 1 - region 2a
!  interaction.
!
!  Distributed memory parallel version.
!
!  imode = 1 => defective region 1 - region 2a
!
!  imode = 2 => perfect region 1 - region 2a - NB shouldn't be called
!  for derivatives 
!
!   5/17 Created from recip12a
!   5/17 elat renamed to erecip
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
  use g_constants
  use control
  use current
  use defects
  use derivatives
  use kspace
  use parallel
  use species
  use symmetry
  use times
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  use two
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                   :: imode
  logical,     intent(in)                   :: lgrad1
  logical,     intent(in)                   :: lgrad2
  real(dp),    intent(out)                  :: ec6
  real(dp),    intent(out)                  :: erecip
!
!  Local variables
!
  integer(i4)                               :: i
  integer(i4)                               :: idk
  integer(i4)                               :: ii
  integer(i4)                               :: iloc
  integer(i4)                               :: ind
  integer(i4)                               :: iv
  integer(i4)                               :: ix
  integer(i4)                               :: iy
  integer(i4)                               :: iz
  integer(i4)                               :: j
  integer(i4)                               :: jj
  integer(i4)                               :: jx
  integer(i4)                               :: jy
  integer(i4)                               :: jz
  integer(i4)                               :: kk
  integer(i4)                               :: n
  integer(i4)                               :: nati
  integer(i4)                               :: natj
  integer(i4)                               :: ni
  integer(i4)                               :: npsi
  integer(i4)                               :: ntypj
  integer(i4)                               :: ntypi
  integer(i4)                               :: status
  logical                                   :: ldoc6
  logical                                   :: lfound
  logical                                   :: lveck
  real(dp)                                  :: arg
  real(dp)                                  :: argci
  real(dp)                                  :: argck
  real(dp)                                  :: arge
  real(dp)                                  :: c6i
  real(dp)                                  :: c6j
  real(dp)                                  :: c6self2
  real(dp)                                  :: c6t1
  real(dp)                                  :: c6t2
  real(dp)                                  :: c6t3
  real(dp)                                  :: c6t4
  real(dp)                                  :: c6tot
  real(dp)                                  :: cos6
  real(dp)                                  :: cosa
  real(dp)                                  :: cosi
  real(dp)                                  :: cosq
  real(dp)                                  :: csin6
  real(dp)                                  :: csink
  real(dp)                                  :: csink6
  real(dp)                                  :: csinq
  real(dp)                                  :: csprod
  real(dp)                                  :: g_cpu_time
  real(dp)                                  :: d1trm
  real(dp)                                  :: d2trm
  real(dp)                                  :: d21q
  real(dp)                                  :: d22q
  real(dp)                                  :: d23q
  real(dp)                                  :: d24q
  real(dp)                                  :: d25q
  real(dp)                                  :: d26q
  real(dp)                                  :: esum1(2)
  real(dp)                                  :: esum2(2)
  real(dp)                                  :: g_derfc
  real(dp)                                  :: factor
  real(dp), dimension(:), allocatable       :: ktrm6
  real(dp)                                  :: kvv(3)
  real(dp)                                  :: oci
  real(dp)                                  :: ocj
  real(dp)                                  :: phsqk
  real(dp)                                  :: phsqk6
  real(dp)                                  :: phsqksum
  real(dp)                                  :: qfct
  real(dp)                                  :: qli
  real(dp)                                  :: qlj
  real(dp)                                  :: rangstoev
  real(dp)                                  :: rk
  real(dp)                                  :: rk2
  real(dp)                                  :: rketa2
  real(dp)                                  :: rrk2
  real(dp)                                  :: sina
  real(dp)                                  :: sini
  real(dp)                                  :: sinek
  real(dp)                                  :: sinek6
  real(dp)                                  :: sineq
  real(dp)                                  :: sinqx
  real(dp)                                  :: sinqy
  real(dp)                                  :: sinqz
  real(dp)                                  :: time0
  real(dp)                                  :: time1
  real(dp)                                  :: tmps(6)
  real(dp)                                  :: trmk
  real(dp)                                  :: trmk6
  real(dp)                                  :: xal
  real(dp)                                  :: yal
  real(dp)                                  :: zal
  real(dp)                                  :: xci
  real(dp)                                  :: yci
  real(dp)                                  :: zci
  real(dp)                                  :: xd
  real(dp)                                  :: yd
  real(dp)                                  :: zd
  real(dp)                                  :: xrkk
  real(dp)                                  :: yrkk
  real(dp)                                  :: zrkk
  real(dp)                                  :: xpon
#ifdef TRACE
  call trace_in('recip12ad')
#endif
!
  time0 = g_cpu_time()
!
  erecip = 0.0_dp
  ec6 = 0.0_dp
!
  eta4 = 0.25_dp/eta
  rangstoev = 1.0_dp/angstoev
!
!  Check modes
!
  if (imode.eq.2.and.lgrad1) then
    call outerror('recip12ad called with imode 2 - should use recip12a',0_i4)
    call stopnow('recip12ad')
  endif
!
!  Allocate local memory
!
  if (lc6) then
    allocate(ktrm6(nkvec),stat=status)
    if (status/=0) call outofmemory('recip12ad','ktrm6')
  endif
!**********************************
!  Calculate and store k-vectors  *
!**********************************
  if (lc6) then
    c6t1 = vol4pi*rangstoev*sqrtpi/48.0_dp
!
!  Reciprocal space self term
!
    c6self2 = 4.0_dp*c6t1*eta*seta
    if (lra) then
      kvv(1) = kv(1,1)
      kvv(2) = kv(2,2)
      kvv(3) = kv(3,3)
      do i = 1,nkvec
        idk = indk(i)
        ii = (idk/maxindk3) - maxindk
        if (ii.eq.0) then
          factor = 1.0_dp
        else
          factor = 2.0_dp
        endif
        idk = idk - (ii + maxindk)*maxindk3
        jj = (idk/maxindk2) - maxindk
        kk = idk - (jj + maxindk)*maxindk2 - maxindk
        xrk(i) = dble(ii)*kvv(1)
        yrk(i) = dble(jj)*kvv(2)
        zrk(i) = dble(kk)*kvv(3)
        rk2 = xrk(i)*xrk(i) + yrk(i)*yrk(i) + zrk(i)*zrk(i)
        arge = -rk2*eta4
        xpon = exp(arge)
        rrk2 = 1.0_dp/rk2
        ktrm(i) = xpon*vol4pi*factor*rrk2
        rk = sqrt(rk2)
        rketa2 = 0.5_dp*rk/seta
        c6t2 = sqrtpi*g_derfc(rketa2)
        rketa2 = 1.0_dp/rketa2
        c6t3 = 0.5_dp*rketa2**3 - rketa2
        c6t3 = c6t3*xpon
        c6t4 = c6t1*rk2*rk*factor
        ktrm6(i) = c6t4*(c6t2+c6t3)
      enddo
    else
      do i = 1,nkvec
        idk = indk(i)
        ii = (idk/maxindk3) - maxindk
        idk = idk - (ii + maxindk)*maxindk3
        jj = (idk/maxindk2) - maxindk
        kk = idk - (jj + maxindk)*maxindk2 - maxindk
        factor = 2.0_dp
        if (ii.eq.0.and.nkangle.eq.1) then
          factor = 1.0_dp
        elseif (jj.eq.0.and.nkangle.eq.2) then
          factor = 1.0_dp
        elseif (kk.eq.0.and.nkangle.eq.3) then
          factor = 1.0_dp
        elseif (nkangle.eq.0) then
          factor = 1.0_dp
        endif
        xrk(i) = ii*kv(1,1) + jj*kv(1,2) + kk*kv(1,3)
        yrk(i) = ii*kv(2,1) + jj*kv(2,2) + kk*kv(2,3)
        zrk(i) = ii*kv(3,1) + jj*kv(3,2) + kk*kv(3,3)
        rk2 = xrk(i)*xrk(i) + yrk(i)*yrk(i) + zrk(i)*zrk(i)
        arge = -rk2*eta4
        xpon = exp(arge)
        rrk2 = 1.0_dp/rk2
        ktrm(i) = xpon*vol4pi*factor*rrk2
        rk = sqrt(rk2)
        rketa2 = 0.5_dp*rk/seta
        c6t2 = sqrtpi*g_derfc(rketa2)
        rketa2 = 1.0_dp/rketa2
        c6t3 = 0.5_dp*rketa2**3-rketa2
        c6t3 = c6t3*xpon
        c6t4 = c6t1*rk2*rk*factor
        ktrm6(i) = c6t4*(c6t2+c6t3)
      enddo
    endif
  else
    if (lra) then
      kvv(1) = kv(1,1)
      kvv(2) = kv(2,2)
      kvv(3) = kv(3,3)
      do i = 1,nkvec
        idk = indk(i)
        ii = (idk/maxindk3) - maxindk
        if (ii.eq.0) then
          factor = 1.0_dp
        else
          factor = 2.0_dp
        endif
        idk = idk - (ii + maxindk)*maxindk3
        jj = (idk/maxindk2) - maxindk
        kk = idk - (jj + maxindk)*maxindk2 - maxindk
        xrk(i) = ii*kvv(1)
        yrk(i) = jj*kvv(2)
        zrk(i) = kk*kvv(3)
        rk2 = xrk(i)*xrk(i) + yrk(i)*yrk(i) + zrk(i)*zrk(i)
        arge = - rk2*eta4
        xpon = exp(arge)
        rrk2 = 1.0_dp/rk2
        ktrm(i) = xpon*vol4pi*factor*rrk2
      enddo
    else
      do i = 1,nkvec
        idk = indk(i)
        ii = (idk/maxindk3) - maxindk
        idk = idk - (ii + maxindk)*maxindk3
        jj = (idk/maxindk2) - maxindk
        kk = idk - (jj + maxindk)*maxindk2 - maxindk
        factor = 2.0_dp
        if (ii.eq.0.and.nkangle.eq.1) then
          factor = 1.0_dp
        elseif (jj.eq.0.and.nkangle.eq.2) then
          factor = 1.0_dp
        elseif (kk.eq.0.and.nkangle.eq.3) then
          factor = 1.0_dp
        elseif (nkangle.eq.0) then
          factor = 1.0_dp
        endif
        xrk(i) = ii*kv(1,1) + jj*kv(1,2) + kk*kv(1,3)
        yrk(i) = ii*kv(2,1) + jj*kv(2,2) + kk*kv(2,3)
        zrk(i) = ii*kv(3,1) + jj*kv(3,2) + kk*kv(3,3)
        rk2 = xrk(i)*xrk(i) + yrk(i)*yrk(i) + zrk(i)*zrk(i)
        arge = - rk2*eta4
        xpon = exp(arge)
        rrk2 = 1.0_dp/rk2
        ktrm(i) = xpon*vol4pi*factor*rrk2
      enddo
    endif
  endif
!
!  End of set-up section
!
  if (lnorecip) goto 999
  ind = 3*nreg1
  if (ldbsm) ind = ind + nreg1
  jx = ind + 1
  jy = ind + 2
  jz = ind + 3
  lveck = (nkvec.ge.numat)
!
!  Modification added to increase speed - outer loop
!  over k vectors is faster for derivatives as
!  recalculation of cos and sin is avoided.
!
  if (lgrad1.or.lgrad2) lveck = .false.
!
!  If Ewald sum for dispersion then don't use lveck
!  as this requires more vector storage
!
  if (lc6) lveck = .false.
!*****************************
!  Vectorise over k vectors  *
!*****************************
  if (lveck) then
!
!  Start loop over cluster atom - unit cell atom pairs
!
    do iv = 1,nkvec
      csin(iv) = 0.0_dp
      sine(iv) = 0.0_dp
    enddo
    do i = 1,numat
      qli = qf(i)*occuf(i)
      xci = xclat(i)
      yci = yclat(i)
      zci = zclat(i)
      do iv = 1,nkvec
        argc(iv) = xrk(iv)*xci
        argc(iv) = yrk(iv)*yci + argc(iv)
        argc(iv) = zrk(iv)*zci + argc(iv)
      enddo
      do iv = 1,nkvec
        csin(iv) = csin(iv) + qli*cos(argc(iv))
        sine(iv) = sine(iv) + qli*sin(argc(iv))
      enddo
    enddo
    do iloc = 1,nreg1onnode
      i = node2reg1(iloc)
      oci = occdefe(i)
      qli = qdefe(i)*oci
      xal = xdefe(i)
      yal = ydefe(i)
      zal = zdefe(i)
!**********************
!  Lattice energy     *
!**********************
      do iv = 1,nkvec
        argci = xrk(iv)*xal + yrk(iv)*yal + zrk(iv)*zal
        cosi = cos(argci)
        sini = sin(argci)
        erecip = erecip + qli*ktrm(iv)*(cosi*csin(iv)+sini*sine(iv))
      enddo
    enddo
  elseif (lc6.and..not.lc6one) then
!************************************************
!  Algorithm for cases where dispersion cannot  *
!  be factorised into one centre terms          *
!************************************************
    do iloc = 1,nreg1onnode
      i = node2reg1(iloc)
      oci = occdefe(i)
      qli = qdefe(i)*oci
      nati = natdefe(i)
      ntypi = ntypdefe(i)
      xal = xdefe(i)
      yal = ydefe(i)
      zal = zdefe(i)
      ni = 3*(iloc-1)
      ix = ni + 1
      iy = ix + 1
      iz = iy + 1
      do j = 1,numat
        ocj = occuf(j)
        qlj = qf(j)*ocj
        natj = nat(j)
        ntypj = nftype(j)
!
!  Find C6 term for pair
!
        c6tot = 0.0_dp
        do n = 1,npote
          if (nati.eq.nspec1(n).and.natj.eq.nspec2(n)) then
            if ((ntypi.eq.nptyp1(n).or.nptyp1(n).eq.0).and. &
                (ntypj.eq.nptyp2(n).or.nptyp2(n).eq.0)) then
              if (nptype(n).eq.1.or.nptype(n).eq.7) then
                c6tot = c6tot + twopot(3,n)
              elseif (nptype(n).eq.2.or.nptype(n).eq.21) then
                c6tot = c6tot + twopot(2,n)
              elseif (nptype(n).eq.57) then
                c6tot = c6tot + twopot(3,n)*twopot(4,n)*twopot(5,n)**6
              endif
            endif
          elseif (natj.eq.nspec1(n).and.nati.eq.nspec2(n)) then
            if ((ntypj.eq.nptyp1(n).or.nptyp1(n).eq.0).and. &
                (ntypi.eq.nptyp2(n).or.nptyp2(n).eq.0)) then
              if (nptype(n).eq.1.or.nptype(n).eq.7) then
                c6tot = c6tot + twopot(3,n)
              elseif (nptype(n).eq.2.or.nptype(n).eq.21) then
                c6tot = c6tot + twopot(2,n)
              elseif (nptype(n).eq.57) then
                c6tot = c6tot + twopot(3,n)*twopot(4,n)*twopot(5,n)**6
              endif
            endif
          endif
        enddo
        ldoc6 = (abs(c6tot).gt.1.0d-4)
!
!  Find relative vector between atoms
!
        xd = xclat(j) - xal
        yd = yclat(j) - yal
        zd = zclat(j) - zal
        qfct = qli*qlj
        csinq = 0.0_dp
        if (lgrad1) then
          sinqx = 0.0_dp
          sinqy = 0.0_dp
          sinqz = 0.0_dp
          if (lgrad2) then
            d21q = 0.0_dp
            d22q = 0.0_dp
            d23q = 0.0_dp
            d24q = 0.0_dp
            d25q = 0.0_dp
            d26q = 0.0_dp
          endif
        endif
        if (ldoc6) then
          c6tot = c6tot*oci*ocj
          if (lgrad1) then
            csin6 = 0.0_dp
            do iv = 1,nkvec
              xrkk = xrk(iv)
              yrkk = yrk(iv)
              zrkk = zrk(iv)
              arg = xrkk*xd + yrkk*yd + zrkk*zd
              cosa = cos(arg)
              sina = sin(arg)
              cosq = cosa*ktrm(iv)*qfct
              cos6 = cosa*ktrm6(iv)*c6tot
              csinq = csinq + cosq
              csin6 = csin6 + cos6
              d1trm = (ktrm(iv)*qfct - ktrm6(iv)*c6tot)*sina
              sinqx = sinqx + d1trm*xrkk
              sinqy = sinqy + d1trm*yrkk
              sinqz = sinqz + d1trm*zrkk
              if (lgrad2) then
                d2trm = cosq - cos6
                d21q = d21q + d2trm*xrk(iv)*xrk(iv)
                d22q = d22q + d2trm*yrk(iv)*yrk(iv)
                d23q = d23q + d2trm*zrk(iv)*zrk(iv)
                d24q = d24q + d2trm*yrk(iv)*zrk(iv)
                d25q = d25q + d2trm*xrk(iv)*zrk(iv)
                d26q = d26q + d2trm*xrk(iv)*yrk(iv)
              endif
            enddo
          else
            csin6 = 0.0_dp
            do iv = 1,nkvec
              arg = xrk(iv)*xd + yrk(iv)*yd + zrk(iv)*zd
              cosa = cos(arg)
              csinq = csinq + cosa*ktrm(iv)
              csin6 = csin6 + cosa*ktrm6(iv)
            enddo
            csinq = csinq*qfct
            csin6 = csin6*c6tot
          endif
!
!  Lattice energy
!
          erecip = erecip + csinq
          ec6 = ec6 - (csin6+c6self2*c6tot)
        else
          if (lgrad1) then
            do iv = 1,nkvec
              xrkk = xrk(iv)
              yrkk = yrk(iv)
              zrkk = zrk(iv)
              arg = xrkk*xd + yrkk*yd + zrkk*zd
              cosa = cos(arg)*qfct
              sina = sin(arg)*qfct
              cosq = cosa*ktrm(iv)
              csinq = csinq + cosq
              sineq = sina*ktrm(iv)
              sinqx = sinqx + sineq*xrkk
              sinqy = sinqy + sineq*yrkk
              sinqz = sinqz + sineq*zrkk
              if (lgrad2) then
                d21q = d21q + cosq*xrk(iv)*xrk(iv)
                d22q = d22q + cosq*yrk(iv)*yrk(iv)
                d23q = d23q + cosq*zrk(iv)*zrk(iv)
                d24q = d24q + cosq*yrk(iv)*zrk(iv)
                d25q = d25q + cosq*xrk(iv)*zrk(iv)
                d26q = d26q + cosq*xrk(iv)*yrk(iv)
              endif
            enddo
          else
            csinq = 0.0_dp
            do iv = 1,nkvec
              arg = xrk(iv)*xd + yrk(iv)*yd + zrk(iv)*zd
              cosa = cos(arg)
              csinq = csinq + cosa*ktrm(iv)
            enddo
            csinq = csinq*qfct
          endif
!
!  Lattice energy
!
          erecip = erecip + csinq
        endif
!
!  Internal derivatives
!
        if (lgrad1) then
          xdrv(i) = xdrv(i) + sinqx
          ydrv(i) = ydrv(i) + sinqy
          zdrv(i) = zdrv(i) + sinqz
          if (lgrad2) then
            derv2(jx,ix) = derv2(jx,ix) + d21q
            derv2(jx,iy) = derv2(jx,iy) + d26q
            derv2(jx,iz) = derv2(jx,iz) + d25q
            derv2(jy,ix) = derv2(jy,ix) + d26q
            derv2(jy,iy) = derv2(jy,iy) + d22q
            derv2(jy,iz) = derv2(jy,iz) + d24q
            derv2(jz,ix) = derv2(jz,ix) + d25q
            derv2(jz,iy) = derv2(jz,iy) + d24q
            derv2(jz,iz) = derv2(jz,iz) + d23q
          endif
        endif
      enddo
    enddo
  else
!***********************************
!  Vectorise over number of atoms  *
!***********************************
!
!  Start loop over k vectors
!
    do iv = 1,nkvec
      csink = 0.0_dp
      sinek = 0.0_dp
      xrkk = xrk(iv)
      yrkk = yrk(iv)
      zrkk = zrk(iv)
      trmk = ktrm(iv)
      if (lc6one) then
        trmk6 = ktrm6(iv)
        csink6 = 0.0_dp
        sinek6 = 0.0_dp
        do i = 1,numat
          oci = occuf(i)
          qli = qf(i)*oci
          c6i = c6f(i)*oci
          xci = xclat(i)
          yci = yclat(i)
          zci = zclat(i)
          argc(i) = xrkk*xci + yrkk*yci + zrkk*zci
          csin(i) = cos(argc(i))
          sine(i) = sin(argc(i))
          csink6 = csink6 + csin(i)*c6i
          sinek6 = sinek6 + sine(i)*c6i
          csink = csink + csin(i)*qli
          sinek = sinek + sine(i)*qli
        enddo
      else
        do i = 1,numat
          qli = qf(i)*occuf(i)
          xci = xclat(i)
          yci = yclat(i)
          zci = zclat(i)
          argc(i) = xrkk*xci + yrkk*yci + zrkk*zci
          csin(i) = cos(argc(i))*qli
          sine(i) = sin(argc(i))*qli
          csink = csink + csin(i)
          sinek = sinek + sine(i)
        enddo
      endif
!
      if (lgrad2) then
        tmps(1) = xrkk*xrkk
        tmps(2) = yrkk*yrkk
        tmps(3) = zrkk*zrkk
        tmps(4) = yrkk*zrkk
        tmps(5) = xrkk*zrkk
        tmps(6) = xrkk*yrkk
      endif
!************************
!  Internal derivatives *
!************************
!
!  First and second derivatives
!
      do iloc = 1,nreg1onnode
        i = node2reg1(iloc)
        oci = occdefe(i)
        qli = qdefe(i)*oci
        nati = natdefe(i)
        ntypi = ntypdefe(i)
        xal = xdefe(i)
        yal = ydefe(i)
        zal = zdefe(i)
        npsi = 0
        if (lc6one) then
          lfound = .false.
          j = 0
          do while (j.lt.nspec.and..not.lfound)
            j = j + 1
            if (natspec(j).eq.nati.and.(ntypi.eq.ntypspec(j).or.ntypspec(j).eq.0)) lfound = .true.
          enddo
          if (lfound) then
            c6i = c6spec(j)*oci
          else
            c6i = 0.0_dp
          endif
        endif
        argci = xal*xrkk + yal*yrkk + zal*zrkk
        cosi = cos(argci)
        sini = sin(argci)
!**********************
!  Lattice energy     *
!**********************
        erecip = erecip + trmk*qli*(cosi*csink+sini*sinek)
        if (lc6one) ec6 = ec6 - trmk6*c6i*(cosi*csink6+sini*sinek6)
        if (lgrad1) then
!
!  Excursion into second derivatives
!
          if (lgrad2) then
            ni = 3*(iloc-1)
            ix = ni + 1
            iy = ix + 1
            iz = iy + 1
            do j = 1,numat
              ocj = occuf(j)
              qlj = qf(j)*ocj
              csprod = (cosi*csin(j)+sini*sine(j))
              if (lc6one) then
                c6j = c6f(j)*ocj
                argck = csprod*(trmk*qli*qlj-trmk6*c6i*c6j)
              else
                argck = trmk*csprod*qli
              endif
              derv2(jx,ix) = derv2(jx,ix) + argck*tmps(1)
              derv2(jy,iy) = derv2(jy,iy) + argck*tmps(2)
              derv2(jz,iz) = derv2(jz,iz) + argck*tmps(3)
              derv2(jx,iy) = derv2(jx,iy) + argck*tmps(6)
              derv2(jx,iz) = derv2(jx,iz) + argck*tmps(5)
              derv2(jy,iz) = derv2(jy,iz) + argck*tmps(4)
            enddo
          endif
!
!  Return to first derivatives
!
          if (lc6one) then
            phsqk = qli*(cosi*sinek-sini*csink)
            phsqk6 = c6i*(cosi*sinek6-sini*csink6)
            phsqksum = phsqk*trmk - phsqk6*trmk6
          else
            phsqk = qli*(cosi*sinek-sini*csink)
            phsqksum = trmk*phsqk
          endif
          xdrv(i) = xdrv(i) + phsqksum*xrkk
          ydrv(i) = ydrv(i) + phsqksum*yrkk
          zdrv(i) = zdrv(i) + phsqksum*zrkk
        endif
      enddo
!*******************************
!  End of loop over k vectors  *
!*******************************
    enddo
    if (lc6one) then
!
!  Self term
!
      do iloc = 1,nreg1onnode
        i = node2reg1(iloc)
        oci = occdefe(i)
        qli = qdefe(i)*oci
        nati = natdefe(i)
        ntypi = ntypdefe(i)
        xal = xdefe(i)
        yal = ydefe(i)
        zal = zdefe(i)
        npsi = 0
        lfound = .false.
        j = 0
        do while (j.lt.nspec.and..not.lfound)
          j = j + 1
          if (natspec(j).eq.nati.and.(ntypi.eq.ntypspec(j).or.ntypspec(j).eq.0)) lfound = .true.
        enddo
        if (lfound) then
          c6i = c6spec(j)*oci
        else
          c6i = 0.0_dp
        endif
        do j = 1,numat
          ocj = occuf(j)
          c6j = c6f(j)*ocj
          ec6 = ec6 - c6self2*c6i*c6j
        enddo
      enddo
    endif
    if (lgrad2) then
      do iloc = 1,nreg1onnode
        ni = 3*(iloc-1)
        ix = ni + 1
        iy = ix + 1
        iz = iy + 1
        derv2(jy,ix) = derv2(jx,iy)
        derv2(jz,ix) = derv2(jx,iz)
        derv2(jz,iy) = derv2(jy,iz)
      enddo
    endif
  endif
!
!  Convert units to eV
!
  erecip = erecip*angstoev
  ec6 = ec6*angstoev
  if (lgrad1) then
    do i = 1,nreg1
      xdrv(i) = angstoev*xdrv(i)
      ydrv(i) = angstoev*ydrv(i)
      zdrv(i) = angstoev*zdrv(i)
    enddo
    if (lgrad2) then
      do i = 1,3*nreg1onnode
        derv2(jx,i) = angstoev*derv2(jx,i)
        derv2(jy,i) = angstoev*derv2(jy,i)
        derv2(jz,i) = angstoev*derv2(jz,i)
      enddo
    endif
  endif
!
!  Sum terms if this is not lveck algorithm
!
  if (lveck) then
    esum1(1) = erecip
    call sumall(esum1,esum2,1_i4,"recip12ad","esum")
    erecip = esum2(1)
  else
    esum1(1) = erecip
    esum1(2) = ec6
    call sumall(esum1,esum2,2_i4,"recip12ad","esum")
    erecip = esum2(1)
    ec6 = esum2(2)
  endif
!***************
!  Exit point  *
!***************
999 continue
!
!  Free local memory
!
  if (lc6) then
    deallocate(ktrm6,stat=status)
    if (status/=0) call deallocate_error('recip12ad','ktrm6')
  endif
!
!  Timing
!
  time1 = g_cpu_time()
  treg1 = treg1 + time1 - time0
#ifdef TRACE
  call trace_out('recip12ad')
#endif
  return
  end
