  subroutine real1Dfd(matom,erealin,lgrad1)
!
!  This subroutine calculates the electrostatic energy of a 1-D system 
!  in real space based on the algorithm implemented in CRYSTAL. A 
!  neutralising uniform background charge density is applied and then 
!  subtracted again. Because a sum over neutral unit cells is required, 
!  this code is kept separate from the other real space routines. The 
!  conventional real space routines are used to handle the Coulomb
!  subtraction issues in order to keep this routine simple.
!  Finite difference version that focuses on derivatives of matom.
!
!  12/17 Created from real1Dmd
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
!  No claim is made that this program is free from errors and
!  no liability will be accepted for any loss or damage that
!  may result. The user is responsible for checking the validity
!  of their results.
!
!  Copyright Curtin University 2018
!
!  Julian Gale, CIC, Curtin University, Feb 2018
!
  use g_constants,    only : angstoev
  use configurations, only : nregionno, nregiontype, QMMMmode
  use control,        only : lnoreal, lDoQDeriv1
  use current
  use derivatives
  use element,        only : maxele
  use general,        only : nmaxcells, nemorder, smallself
  use kspace,         only : accf1D
  use qmedata,        only : maxloop
  use shells,         only : cuts
  use symmetry,       only : lstr
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed arguments
!
  integer(i4), intent(in)    :: matom
  real(dp),    intent(inout) :: erealin
  logical,     intent(in)    :: lgrad1
!
!  Local variables
!
  integer(i4)                :: i
  integer(i4)                :: j 
  integer(i4)                :: m
  integer(i4)                :: nati
  integer(i4)                :: natj
  integer(i4)                :: nregioni
  integer(i4)                :: nregionj
  integer(i4)                :: nregiontypi
  integer(i4)                :: nregiontypj
  logical                    :: lconverged
  logical                    :: lcspair
  real(dp)                   :: accf
  real(dp)                   :: acell
  real(dp)                   :: g_cpu_time
  real(dp)                   :: cut2s
  real(dp)                   :: d0
  real(dp)                   :: d0i
  real(dp)                   :: d0j
  real(dp)                   :: d0term
  real(dp)                   :: d1
  real(dp)                   :: d1s
  real(dp)                   :: dh1(3)
  real(dp)                   :: dh2(3)
  real(dp)                   :: dh1s
  real(dp)                   :: dh2s
  real(dp)                   :: d2h1(6)
  real(dp)                   :: d2h2(6)
  real(dp)                   :: d2h1m(3)
  real(dp)                   :: d2h2m(3)
  real(dp)                   :: d2h1s
  real(dp)                   :: d2h2s
  real(dp)                   :: d3h1(10)
  real(dp)                   :: d3h1m(6)
  real(dp)                   :: ediff
  real(dp)                   :: elast
  real(dp)                   :: ereal
  real(dp)                   :: esum
  real(dp)                   :: esumem
  real(dp)                   :: esumh
  real(dp)                   :: esum12
  real(dp)                   :: esumem12
  real(dp)                   :: esumh12
  real(dp)                   :: esum2
  real(dp)                   :: esumem2
  real(dp)                   :: esumh2
  real(dp)                   :: e1
  real(dp)                   :: e2
  real(dp)                   :: h1
  real(dp)                   :: h2
  real(dp)                   :: lna
  real(dp)                   :: oci
  real(dp)                   :: ocj
  real(dp)                   :: qi
  real(dp)                   :: qj
  real(dp)                   :: qii
  real(dp)                   :: qij
  real(dp)                   :: r
  real(dp)                   :: rcut
  real(dp)                   :: rr
  real(dp)                   :: t1, t2
  real(dp)                   :: u
  real(dp)                   :: x
  real(dp)                   :: y
  real(dp)                   :: z
!
!  If noreal specified, return
!
  if (lnoreal) then
    erealin = 0.0_dp
    return
  endif
#ifdef TRACE
  call trace_in('real1Dfd')
#endif
!
  t1 = g_cpu_time()
!********************************************************
!  Calculate Coulomb sum converged to desired accuracy  *
!********************************************************
!
!  Loop over number of cells in sum
!
  accf = 10.0**(-accf1D)
  lna = log(a)
  cut2s = cuts*cuts
  m = - 1
  lconverged = .false.
  elast = 0.0_dp
  esum = 0.0_dp
  esum12 = 0.0_dp
  esum2 = 0.0_dp
!
  do while (m.lt.nmaxcells.and..not.lconverged) 
    m = m + 1
!********************************************
!  Direct sum component over neutral cells  *
!********************************************
    acell = dble(m)*a
!
!  Only compute i = matom case
!
    i = matom
    nati = nat(i)
    oci = occuf(i)
    qi = qf(i)*oci
    nregioni = nregionno(nsft+i)
    nregiontypi = nregiontype(nregioni,ncf)
    jloop: do j = 1,numat
!
!  Exclude self term
!
      if (i.eq.j) cycle jloop
!
      natj = nat(j)
      ocj = occuf(j)
      qj = qf(j)*ocj
      nregionj = nregionno(nsft+j)
      nregiontypj = nregiontype(nregionj,ncf)
      lcspair = (abs(nati-natj).eq.maxele.or.(oci+ocj).lt.1.0001_dp)
      if (lcspair) then
        rcut = cut2s
      else
        rcut = smallself
      endif
!  
!  QM/MM handling : i & j are both QM atoms => exclude
!
      if (QMMMmode(ncf).gt.0) then
        if (nregiontypi.eq.1.and.nregiontypj.eq.1) cycle jloop
!
!  QM/MM : Electrostatic embedding : If either i or j are QM atoms => exclude 
!
        if (QMMMmode(ncf).eq.2.and.(nregiontypi.eq.1.or.nregiontypj.eq.1)) cycle jloop
      endif
!
      x = acell + xclat(j) - xclat(i)
      y = yclat(j) - yclat(i)
      z = zclat(j) - zclat(i)
      r = x*x + y*y + z*z
      if (r.gt.rcut) then
        r = sqrt(r)
        rr = 1.0_dp/r
        d0 = qi*qj*rr
        esum = esum + d0
        if (lgrad1) then
          d1 = d0*rr*rr*angstoev
          xdrv(i) = xdrv(i) + d1*x
          ydrv(i) = ydrv(i) + d1*y
          zdrv(i) = zdrv(i) + d1*z
          xdrv(j) = xdrv(j) - d1*x
          ydrv(j) = ydrv(j) - d1*y
          zdrv(j) = zdrv(j) - d1*z
          if (lstr) then
            rstrd(1) = rstrd(1) - d1*x*x
          endif
          if (lDoQDeriv1) then
            d0i = qj*rr*angstoev
            d0j = qi*rr*angstoev
            call d1charge(i,j,.true.,.true.,1_i4,d0i,d0j)
          endif
        endif
      endif
      if (m.gt.0) then
        x = - acell + xclat(j) - xclat(i)
        r = x*x + y*y + z*z
        if (r.gt.rcut) then
          r = sqrt(r)
          rr = 1.0_dp/r
          d0 = qi*qj*rr
          esum = esum + d0
          if (lgrad1) then
            d1 = d0*rr*rr*angstoev
            xdrv(i) = xdrv(i) + d1*x
            ydrv(i) = ydrv(i) + d1*y
            zdrv(i) = zdrv(i) + d1*z
            xdrv(j) = xdrv(j) - d1*x
            ydrv(j) = ydrv(j) - d1*y
            zdrv(j) = zdrv(j) - d1*z
            if (lstr) then
              rstrd(1) = rstrd(1) - d1*x*x
            endif
            if (lDoQDeriv1) then
              d0i = qj*rr*angstoev
              d0j = qi*rr*angstoev
              call d1charge(i,j,.true.,.true.,1_i4,d0i,d0j)
            endif
          endif
        endif
      endif
    enddo jloop
!**********************
!  Self interactions  *
!**********************
    oci = occuf(i)
    qi = qf(i)*oci
    nregioni = nregionno(nsft+i)
    nregiontypi = nregiontype(nregioni,ncf)
!  
!  QM/MM handling : i is a QM atom => exclude
!
    if (QMMMmode(ncf).eq.0.or.nregiontypi.ne.1) then
      r = abs(acell)
      if (r.gt.smallself) then
        rr = 1.0_dp/r
        d0 = qi*qi*rr
!
!  Set region 2 pair flags
!        
        esum = esum + d0
        if (lgrad1.and.lstr) then
          d1 = d0*angstoev
          rstrd(1) = rstrd(1) - d1
        endif
        if (lgrad1.and.lDoQDeriv1) then
          d0i = qi*rr*angstoev
          call d1charge(i,i,.true.,.true.,1_i4,d0i,d0i)
        endif
      endif
    endif
!
!  Neutralising terms
!
!  Background
!
!  and
!
!  Euler-MacLaurin component
!
    esumh = 0.0_dp
    esumh12 = 0.0_dp
    esumh2 = 0.0_dp
    esumem = 0.0_dp
    esumem12 = 0.0_dp
    esumem2 = 0.0_dp
!
    if (m.gt.0) then
      u = (dble(m)+0.5_dp)*a
      jloop2: do j = 1,numat
        ocj = occuf(j)
        qj = qf(j)*ocj
        nregionj = nregionno(nsft+j)
        nregiontypj = nregiontype(nregionj,ncf)
!  
!  QM/MM handling : i & j are both QM atoms => exclude
!
        if (QMMMmode(ncf).gt.0) then
          if (nregiontypi.eq.1.and.nregiontypj.eq.1) cycle jloop2
!
!  QM/MM : Electrostatic embedding : If either i or j are QM atoms => exclude 
!
          if (QMMMmode(ncf).eq.2.and.(nregiontypi.eq.1.or.nregiontypj.eq.1)) cycle jloop2
        endif
!
        qij = qi*qj
        x = xclat(j) - xclat(i)
        y = yclat(j) - yclat(i)
        z = zclat(j) - zclat(i)
        call hfunc(u,+x,1.0_dp,y,z,h1,dh1,d2h1,d3h1,.false.,.false.,.false.)
        call hfunc(u,-x,-1.0_dp,y,z,h2,dh2,d2h2,d3h1,.false.,.false.,.false.)
        esumh = esumh - qij*(h1 + h2 - 2.0_dp*lna)/a
!
        call emfunc(nemorder,u,+x,1.0_dp,y,z,a,e1,dh1,d2h1,dh1s,d2h1s,d2h1m,d3h1,d3h1m,.false.,.false.,.false.)
        call emfunc(nemorder,u,-x,-1.0_dp,y,z,a,e2,dh2,d2h2,dh2s,d2h2s,d2h2m,d3h1,d3h1m,.false.,.false.,.false.)
        esumem = esumem + qij*(e1 + e2)
      enddo jloop2
!
!  QM/MM handling : i is a QM atom => exclude
!
      if (QMMMmode(ncf).eq.0.or.nregiontypi.ne.1) then
        qii = qi*qi
        call hfunc(u,0.0_dp,1.0_dp,0.0_dp,0.0_dp,h1,dh1,d2h1,d3h1,.false.,.false.,.false.)
        esumh = esumh - qii*(h1 - lna)/a
!
        call emfunc(nemorder,u,0.0_dp,1.0_dp,0.0_dp,0.0_dp,a,e1,dh1,d2h1,dh1s,d2h1s,d2h1m, &
          d3h1,d3h1m,.false.,.false.,.false.)
        esumem = esumem + qii*e1
      endif
    endif
!
!  Sum up terms
!
    ereal = esum + esumh + esumem
!
!  Compare to energy with previous number of cells
!  and check for convergence to the required
!  accuracy.
!
    if (abs(ereal).lt.accf) lconverged = .true.
    if (.not.lconverged) then
      ediff = abs((ereal - elast)/ereal)
      lconverged = (ediff.lt.accf)
    endif
    elast = ereal
  enddo
!
!  Divide ereal by number of processors to compensate for the fact that it has already been summed
!
  ereal = ereal*angstoev
  erealin = erealin + ereal
!
!  Save number of cells needed
!
  maxloop(1) = m
!
!  Derivatives of Euler-MacLaurin terms since these are not cumulative
!
  if (lgrad1.and.m.gt.0) then
    u = (dble(m)+0.5_dp)*a
    jloop3: do j = 1,numat
!
!  Exclude self term 
!
      if (i.eq.j) cycle jloop3
!
      ocj = occuf(j)
      qj = qf(j)*ocj
      nregionj = nregionno(nsft+j)
      nregiontypj = nregiontype(nregionj,ncf)
!  
!  QM/MM handling : i & j are both QM atoms => exclude
!
      if (QMMMmode(ncf).gt.0) then
        if (nregiontypi.eq.1.and.nregiontypj.eq.1) cycle jloop3
!
!  QM/MM : Electrostatic embedding : If either i or j are QM atoms => exclude 
!
        if (QMMMmode(ncf).eq.2.and.(nregiontypi.eq.1.or.nregiontypj.eq.1)) cycle jloop3
      endif
      qij = qi*qj
      x = xclat(j) - xclat(i)
      y = yclat(j) - yclat(i)
      z = zclat(j) - zclat(i)
      r = sqrt(x*x + y*y + z*z)
      call hfunc(u,+x,1.0_dp,y,z,h1,dh1,d2h1,d3h1,lgrad1,.false.,.false.)
      call hfunc(u,-x,-1.0_dp,y,z,h2,dh2,d2h2,d3h1,lgrad1,.false.,.false.)
!
      d1 = qij*angstoev/a
      xdrv(i) = xdrv(i) + d1*(dh1(1) + dh2(1))
      ydrv(i) = ydrv(i) + d1*(dh1(2) + dh2(2))
      zdrv(i) = zdrv(i) + d1*(dh1(3) + dh2(3))
      xdrv(j) = xdrv(j) - d1*(dh1(1) + dh2(1))
      ydrv(j) = ydrv(j) - d1*(dh1(2) + dh2(2))
      zdrv(j) = zdrv(j) - d1*(dh1(3) + dh2(3))
      if (lstr) then
        d1s = - (dh1(1)*(u+x)-dh2(1)*(u-x))
        d1s = d1s + (h1+h2)
        d1s = d1s + 2.0_dp*(1.0_dp-lna)
        d1s = d1s*angstoev/a
        d1s = d1s*qij
        rstrd(1) = rstrd(1) + d1s
      endif
!
!  Compute variable charge term = E without charges
!
      d0term = - angstoev*(h1 + h2 - 2.0_dp*lna)/a
!
      call emfunc(nemorder,u,+x,1.0_dp,y,z,a,e1,dh1,d2h1,dh1s, &
                  d2h1s,d2h1m,d3h1,d3h1m,lgrad1,.false.,.false.)
      call emfunc(nemorder,u,-x,-1.0_dp,y,z,a,e2,dh2,d2h2,dh2s, &
                  d2h2s,d2h2m,d3h1,d3h1m,lgrad1,.false.,.false.)
!
      d1 = qij*angstoev
      xdrv(i) = xdrv(i) - d1*(dh1(1) + dh2(1))
      ydrv(i) = ydrv(i) - d1*(dh1(2) + dh2(2))
      zdrv(i) = zdrv(i) - d1*(dh1(3) + dh2(3))
      xdrv(j) = xdrv(j) + d1*(dh1(1) + dh2(1))
      ydrv(j) = ydrv(j) + d1*(dh1(2) + dh2(2))
      zdrv(j) = zdrv(j) + d1*(dh1(3) + dh2(3))
      if (lstr) then
        rstrd(1) = rstrd(1) + d1*(dh1s + dh2s)
      endif
      d0term = d0term + (e1 + e2)*angstoev 
      if (lDoQDeriv1) then
        d0i = qj*d0term
        d0j = qi*d0term
        call d1charge(i,j,.true.,.true.,1_i4,d0i,d0j)
      endif
    enddo jloop3
    if (lstr.or.lDoQDeriv1) then
!***************************
!  Self-interaction terms  *
!***************************
!
!  QM/MM handling : i is a QM atom => exclude
!
      if (QMMMmode(ncf).eq.0.or.nregiontypi.ne.1) then
        qii = qi*qi
        call hfunc(u,0.0_dp,1.0_dp,0.0_dp,0.0_dp,h1,dh1,d2h1,d3h1,lgrad1,.false.,.false.)
!
        if (lstr) then
          d1 = angstoev/a
          d1s = - d1*dh1(1)*u
          d1s = d1s + d1*h1
          d1s = d1s + d1*(1.0_dp-lna)
          d1 = d1*qii
          d1s = d1s*qii
          rstrd(1) = rstrd(1) + d1s
        endif
!
!  Compute variable charge term = E without charges
!
        d0term = - angstoev*(h1 - lna)/a
        call emfunc(nemorder,u,0.0_dp,1.0_dp,0.0_dp,0.0_dp,a,e1,dh1,d2h1, &
                    dh1s,d2h1s,d2h1m,d3h1,d3h1m,lgrad1,.false.,.false.)
!
        if (lstr) then
          rstrd(1) = rstrd(1) + qii*dh1s*angstoev
        endif
        d0term = d0term + e1*angstoev
        if (lDoQDeriv1) then
          d0i = qi*d0term
          call d1charge(i,i,.true.,.true.,1_i4,d0i,d0i)
        endif
      endif
    endif
  endif
!
  t2 = g_cpu_time()
  tatom = tatom + t2 - t1
#ifdef TRACE
  call trace_out('real1Dfd')
#endif
!
  return
  end
