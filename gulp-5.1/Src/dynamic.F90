  subroutine dynamic(xk,yk,zk)
!
!  Master routine for generating phased second derivatives
!  On diagonal block elements (i,i) for element i can no longer be
!  generated by summing off diagonals. Therefore these elements
!  are retained from the normal second derivative matrix and
!  added in afterwards in phonon.
!
!   8/95 Molecule fixing option added
!   3/97 Modified so that only a single K point can be
!        calculated on each call for simplicity
!   5/97 Sutton-Chen potential modifications added. Note that the
!        bulk rho values should be set before calling dynamic
!        from the previous function call.
!   8/97 Geometry set up removed as should never be needed as
!        call to phonon will be preceeded by call to energy.
!   4/01 Calculation of phase factor moved into dynamic from
!        lower level routines.
!   4/01 Modified to allow the option for K points to be specified
!        with reference to the full centred unit cell
!   5/02 Electrostatic contribution for polymers completed
!   8/02 Brenner potential added
!  11/02 Einstein model added
!  11/03 Bond order potential added
!   9/04 Call to bond order charge derivatives added
!   7/06 Sixbody contribution added
!   4/09 Separate call to generate MEAM densities added
!   6/09 Module name changed from three to m_three
!   9/10 Neutron scattering modifications added
!   7/12 Modifications for parallel distributed second derivatives added
!   8/14 Modifications for group velocities added
!   8/14 eatom added to call to density routines
!  10/14 Changed so that fractional k point is passed in
!   2/15 Criterion for calling getBOcharge changed by adding nboQ0
!   2/15 Now calls BOselfp rather than BOself
!   9/16 Unused variable removed
!  12/16 Initialisation of derv2dk modified
!   1/17 Bondorder distributed second derivative calls added
!   2/18 Trace added
!   5/18 Modifications for EEM split bond charge added
!   7/18 emat allocated and passed to dcharge routines if needed
!   7/18 Right-hand dimension of emat corrected for case where dcharge is called
!  11/18 Initialisation of density arrays for manybody potentials added
!   4/19 Density now recomputed for EAM as well as MEAM
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
!  Copyright Curtin University 2019
!
!  Julian Gale, CIC, Curtin University, April 2019
!
  use bondorderdata,     only : nbopot, nboQ, nboQ0
  use control
  use current
  use derivatives
  use eam,               only : lMEAM, lMEAMden, maxmeamcomponent
  use eembonds,          only : neembond
  use four
  use ksample
  use ksample_scatter
  use parallel,          only : nprocs, natomsonnode
  use six
  use sutton
  use m_three
#ifdef TRACE
  use trace,             only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  real(dp),    intent(in)                      :: xk
  real(dp),    intent(in)                      :: yk
  real(dp),    intent(in)                      :: zk
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: j
  integer(i4)                                  :: maxlim
  integer(i4)                                  :: maxlimu
  integer(i4)                                  :: mint
  integer(i4)                                  :: nmax
  integer(i4)                                  :: nmaxu
  integer(i4)                                  :: status
  real(dp)                                     :: eatom
  real(dp)                                     :: ebondorder
  real(dp)                                     :: ebrenner
  real(dp)                                     :: eeinstein
  real(dp),    dimension(:,:), allocatable     :: emat
  real(dp),    dimension(:,:), allocatable     :: ematb
  real(dp)                                     :: kvf(3,3)
  real(dp)                                     :: xkv
  real(dp)                                     :: ykv
  real(dp)                                     :: zkv
#ifdef TRACE
  call trace_in('dynamic')
#endif
!************************************
!  Initialise density for MEAM/EAM  *
!************************************
  if (lsuttonc) then
    if (lMEAM) then
      do i = 1,numat
        scrho(1:maxmeamcomponent,i) = 0.0_dp
        scrho12(1:maxmeamcomponent,i) = 0.0_dp
      enddo
    else
      do i = 1,numat
        scrho(1,i) = 0.0_dp
        scrho12(1,i) = 0.0_dp
      enddo
    endif
  endif
!**********************************************
!  EEM/QEq calculation of charge derivatives  *
!**********************************************
  if (leem) then
!
!  Allocate arrays for matrix equations
!
    nmaxu = numat + 1
    nmax = numat + 1
!
    allocate(emat(nmax,nmaxu),stat=status)
    if (status/=0) call outofmemory('dynamic','emat')
!
    if (leembond) then
      allocate(ematb(nmax,nmaxu),stat=status)
      if (status/=0) call outofmemory('dynamic','ematb')
!
      call dchargesplit(.false.,emat,nmax,ematb,neembond,.true.,.false.)
!
      deallocate(ematb,stat=status)
      if (status/=0) call deallocate_error('dynamic','ematb')
    else
      call dcharge(.false.,emat,nmax,.true.,.false.)
    endif
!
    deallocate(emat,stat=status)
    if (status/=0) call deallocate_error('dynamic','emat')
  endif
!*************************************************
!  Bond Order calculation of charge derivatives  *
!*************************************************
  if ((nboQ+nboQ0).gt.0) then
    call getBOcharge(.true.,.true.)
  endif
!****************************
!  Zero second derivatives  *
!****************************
  mint = 3*numat
  maxlim = mint
  if (nbsmat.gt.0) maxlim = maxlim + numat
  if (nprocs.gt.1) then
    maxlimu = 3_i4*natomsonnode
  else
    maxlimu = maxlim
  endif
  if (maxlimu.gt.maxd2u) then
    maxd2u = maxlimu
    call changemaxd2
  endif
  if (maxlim.gt.maxd2) then
    maxd2 = maxlim
    call changemaxd2
  endif
  do i = 1,maxlimu
    do j = 1,maxlim
      derv2(j,i) = 0.0_dp
      dervi(j,i) = 0.0_dp
    enddo
  enddo
!
!  If group velocities are to be computed then zero k derivatives of dynamical matrix
!
  if (lgroupvelocity) then
    derv2dk(1:3,1:maxlim,1:maxlimu) = 0.0_dpc
  endif
!
!  Select appropriate K vectors
!
  if (lkfull.and.ndim.eq.3) then
    call kvector3Df(kvf)
  else
    kvf(1:3,1:3) = kv(1:3,1:3)
  endif
!***************************
!  Calculate phase factor  *
!***************************
  if (ndim.eq.3) then
    xkv = xk*kvf(1,1) + yk*kvf(1,2) + zk*kvf(1,3)
    ykv = xk*kvf(2,1) + yk*kvf(2,2) + zk*kvf(2,3)
    zkv = xk*kvf(3,1) + yk*kvf(3,2) + zk*kvf(3,3)
  elseif (ndim.eq.2) then
    xkv = xk*kvf(1,1) + yk*kvf(1,2)
    ykv = xk*kvf(2,1) + yk*kvf(2,2)
    zkv = 0.0_dp
  elseif (ndim.eq.1) then
    xkv = xk*kvf(1,1)
    ykv = 0.0_dp
    zkv = 0.0_dp
  endif
!*******************************
!  Reciprocal space component  *
!*******************************
  if (lewald.and.ndim.gt.1) then
    call kindex
    if (ndim.eq.3) then
      if (nprocs.gt.1) then
        call recip3Dpd(xkv,ykv,zkv)
      else
        call recip3Dp(xkv,ykv,zkv)
      endif
    elseif (ndim.eq.2) then
      if (nprocs.gt.1) then
        call recip2Dpd(xkv,ykv)
      else
        call recip2Dp(xkv,ykv)
      endif
    endif
  endif
!*************************
!  Real space component  *
!*************************
  if (nprocs.gt.1) then
    call realpd(xkv,ykv,zkv)
    if (ndim.eq.1) call real1Dpd(xkv)
  else
    call realp(xkv,ykv,zkv)
    if (ndim.eq.1) call real1Dp(xkv)
  endif
!**********************************
!  Bond order charge self-energy  *
!**********************************
  if (nboQ0.gt.0) then
    call BOselfp(xkv,ykv,zkv)
  endif
!*************************
!  Three-body component  *
!*************************
  if (nthb.gt.0) then
    if (nprocs.gt.1) then
      call threepd(xkv,ykv,zkv)
    else
      call threep(xkv,ykv,zkv)
    endif
  endif
!************************
!  Four-body component  *
!************************
  if (nfor.gt.0) then
    if (nprocs.gt.1) then
      call fourpd(xkv,ykv,zkv)
    else
      call fourp(xkv,ykv,zkv)
    endif
  endif
!***********************
!  Six-body component  *
!***********************
  if (nsix.gt.0) then
    if (nprocs.gt.1) then
      call sixpd(xkv,ykv,zkv)
    else
      call sixp(xkv,ykv,zkv)
    endif
  endif
!************************
!  Many-body component  *
!************************
  if (lsuttonc) then
    eatom = 0.0_dp
    call density3(eatom)
!
    if (nprocs.gt.1) then
      call manypd(xkv,ykv,zkv)
    else
      call manyp(xkv,ykv,zkv)
    endif
  endif
!**********************
!  Brenner potential  *
!**********************
  if (lbrenner) then
    ebrenner = 0.0_dp
    if (nprocs.gt.1) then
      call brennerd(ebrenner,xkv,ykv,zkv,.true.,.true.,.true.)
    else
      call brenner(ebrenner,xkv,ykv,zkv,.true.,.true.,.true.)
    endif
  endif
!************************
!  Bondorder potential  *
!************************
  if (nbopot.gt.0) then
    ebondorder = 0.0_dp
    if (nprocs.gt.1) then
      call bondorderd(ebondorder,xkv,ykv,zkv,.true.,.true.,.true.)
    else
      call bondorder(ebondorder,xkv,ykv,zkv,.true.,.true.,.true.)
    endif
  endif
!*****************************
!  Einstein model component  *
!*****************************
  if (leinstein) then
    eeinstein = 0.0_dp
    if (nprocs.gt.1) then
      call einsteinpd(eeinstein,.true.,.true.)
    else
      call einstein(eeinstein,.true.,.true.)
    endif
  endif
#ifdef TRACE
  call trace_out('dynamic')
#endif
!
  return
  end
