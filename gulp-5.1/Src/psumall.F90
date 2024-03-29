  subroutine psumall(eatom,ereal,erecip,ec6,eqeq,eattach,esregion12,esregion2, &
    ethb,efor,eoop,emany,ecmm,ebrenner,epolar,eeinstein,ewolfself,ebondorder, &
    eforce,esix,efield,eradial,ereaxFF,eplane,ecosmo,eone,eedip,eimp,ekim, &
    eplumed,eboQself,lgrad1,lgrad2,lsym)
!
!  Sum all energies and first derivatives across all Nodes
!
!  11/02 Created
!   6/04 Energies now globalised in a single call
!   6/04 ewolfself and ebondorder added
!   8/04 eforce added
!  11/04 esix added
!   3/07 efield added
!   3/07 eradial added
!   7/07 emeta added
!   7/07 ereaxFF added
!   7/07 eplane added
!  10/08 ecosmo added
!   6/09 Virial and site energies added
!  11/09 Region derivatives added
!   1/10 One-body energy added
!   9/10 EDIP energy added
!  11/10 Anisotropic pressure added
!   9/11 Metadynamics internal code replaced with Plumed
!   4/12 Summation of virial removed as no longer needed
!   4/12 xvir, yvir and zvir removed
!   5/12 Atomic stresses added
!  10/12 Support for OpenKIM models added
!   7/13 Improper torsion energy added
!   3/15 Parallel reduction of non-radial forces added
!   3/16 eplumed added
!   9/16 lgrad2 now passed to psumall
!   7/17 eboQself added
!   7/17 eregion2region added 
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
  use configurations, only : maxregion, nregions
  use control,        only : latomicstress, leregion
  use current
  use derivatives
  use eam,            only : lMEAM
  use energies,       only : siteenergy, eregion2region
  use parallel
  use symmetry,       only : lstr
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  logical,     intent(in)                      :: lgrad1
  logical,     intent(in)                      :: lgrad2
  logical,     intent(in)                      :: lsym
  real(dp),    intent(inout)                   :: eatom
  real(dp),    intent(inout)                   :: ereal
  real(dp),    intent(inout)                   :: erecip
  real(dp),    intent(inout)                   :: ec6
  real(dp),    intent(inout)                   :: eqeq
  real(dp),    intent(inout)                   :: eattach
  real(dp),    intent(inout)                   :: esregion12
  real(dp),    intent(inout)                   :: esregion2
  real(dp),    intent(inout)                   :: ethb
  real(dp),    intent(inout)                   :: efor
  real(dp),    intent(inout)                   :: eimp
  real(dp),    intent(inout)                   :: eoop
  real(dp),    intent(inout)                   :: emany
  real(dp),    intent(inout)                   :: ecmm
  real(dp),    intent(inout)                   :: ebrenner
  real(dp),    intent(inout)                   :: epolar
  real(dp),    intent(inout)                   :: eeinstein
  real(dp),    intent(inout)                   :: ewolfself
  real(dp),    intent(inout)                   :: ebondorder
  real(dp),    intent(inout)                   :: eforce
  real(dp),    intent(inout)                   :: esix
  real(dp),    intent(inout)                   :: efield
  real(dp),    intent(inout)                   :: eradial
  real(dp),    intent(inout)                   :: ereaxFF
  real(dp),    intent(inout)                   :: eplane
  real(dp),    intent(inout)                   :: ecosmo
  real(dp),    intent(inout)                   :: eone
  real(dp),    intent(inout)                   :: eedip
  real(dp),    intent(inout)                   :: ekim
  real(dp),    intent(inout)                   :: eplumed
  real(dp),    intent(inout)                   :: eboQself
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ind
  integer(i4)                                  :: j
  integer(i4)                                  :: n
  integer(i4)                                  :: nregs
  integer(i4)                                  :: nrecv
  integer(i4)                                  :: nsend
  integer(i4)                                  :: status
  real(dp)                                     :: g_cpu_time
  real(dp),    dimension(:), allocatable       :: sum
  real(dp),    dimension(:), allocatable       :: sum0
  real(dp)                                     :: tsum0
  real(dp)                                     :: tsuml
!
!  If nprocs = 1 return
!
  if (nprocs.eq.1) return
#ifdef TRACE
  call trace_in('psumall')
#endif
!
!  Allocate workspace array for parallel sums
!
  if (lsym) then
    n = nasym
  else
    n = numat
  endif
  if (leregion) then
    nregs = nregions(ncf)**2
  else
    nregs = 0
  endif
  if (lMEAM) then
    allocate(sum(max(31+nregs,11*n+1+nstrains*(numat+2)+36)),stat=status)
    if (status/=0) call outofmemory('psumall','sum')
    allocate(sum0(max(31+nregs,11*n+1+nstrains*(numat+2)+36)),stat=status)
    if (status/=0) call outofmemory('psumall','sum0')
  else
    allocate(sum(max(31+nregs,8*n+1+nstrains*(numat+1)+36)),stat=status)
    if (status/=0) call outofmemory('psumall','sum')
    allocate(sum0(max(31+nregs,8*n+1+nstrains*(numat+1)+36)),stat=status)
    if (status/=0) call outofmemory('psumall','sum0')
  endif
!****************
!  Global sums  *
!****************
  tsum0 = g_cpu_time()
  sum0(1) = eatom
  sum0(2) = ereal
  sum0(3) = erecip
  sum0(4) = ec6
  sum0(5) = eqeq
  sum0(6) = eattach
  sum0(7) = esregion12
  sum0(8) = esregion2
  sum0(9) = ethb
  sum0(10) = efor
  sum0(11) = eoop
  sum0(12) = emany
  sum0(13) = ecmm
  sum0(14) = ebrenner
  sum0(15) = epolar
  sum0(16) = eeinstein
  sum0(17) = ewolfself
  sum0(18) = ebondorder
  sum0(19) = eforce
  sum0(20) = esix
  sum0(21) = efield
  sum0(22) = eradial
  sum0(23) = ereaxFF
  sum0(24) = eplane
  sum0(25) = ecosmo
  sum0(26) = eone
  sum0(27) = eedip
  sum0(28) = eimp
  sum0(29) = ekim
  sum0(30) = eplumed
  sum0(31) = eboQself
!
  if (leregion) then
    ind = 31
    do i = 1,nregions(ncf)
      do j = 1,nregions(ncf)
        ind = ind + 1
        sum0(ind) = eregion2region(j,i)
      enddo
    enddo
  endif
!
  call sumall(sum0,sum,31_i4+nregs,"psumall","energies")
!
  eatom = sum(1)
  ereal = sum(2)
  erecip = sum(3)
  ec6 = sum(4)
  eqeq = sum(5)
  eattach = sum(6)
  esregion12 = sum(7)
  esregion2 = sum(8)
  ethb = sum(9)
  efor = sum(10)
  eoop = sum(11)
  emany = sum(12)
  ecmm = sum(13)
  ebrenner = sum(14)
  epolar = sum(15)
  eeinstein = sum(16)
  ewolfself = sum(17)
  ebondorder = sum(18)
  eforce = sum(19)
  esix = sum(20)
  efield = sum(21)
  eradial = sum(22)
  ereaxFF = sum(23)
  eplane = sum(24)
  ecosmo = sum(25)
  eone = sum(26)
  eedip = sum(27)
  eimp = sum(28)
  ekim = sum(29)
  eplumed = sum(30)
  eboQself = sum(31)
  if (leregion) then
    ind = 31
    do i = 1,nregions(ncf)
      do j = 1,nregions(ncf)
        ind = ind + 1
        eregion2region(j,i) = sum(ind)
      enddo
    enddo
  endif
!
  if (lgrad1) then
    do i = 1,n
      sum0(i) = xdrv(i)
    enddo
    do i = 1,n
      sum0(n+i) = ydrv(i)
    enddo
    do i = 1,n
      sum0(2*n+i) = zdrv(i)
    enddo
    do i = 1,n
      sum0(3*n+i) = siteenergy(i)
    enddo
    if (lMEAM) then
      do i = 1,n
        sum0(4*n+i) = xdrvnr(i)
      enddo
      do i = 1,n
        sum0(5*n+i) = ydrvnr(i)
      enddo
      do i = 1,n
        sum0(6*n+i) = zdrvnr(i)
      enddo
      nsend = 7*n 
    else
      nsend = 4*n 
    endif
    do i = 1,maxregion
      sum0(nsend+i) = xregdrv(i)
    enddo
    nsend = nsend + maxregion
    do i = 1,maxregion
      sum0(nsend+i) = yregdrv(i)
    enddo
    nsend = nsend + maxregion
    do i = 1,maxregion
      sum0(nsend+i) = zregdrv(i)
    enddo
    nsend = nsend + maxregion
    if (nbsmat.gt.0) then
      do i = 1,n
        sum0(nsend+i) = raderv(i)
      enddo
      nsend = nsend + n
    endif
    if (lstr) then
      do i = 1,nstrains
        sum0(nsend+i) = rstrd(i)
      enddo
      nsend = nsend + nstrains
      if (lMEAM) then
        do i = 1,nstrains
          sum0(nsend+i) = rstrdnr(i)
        enddo
        nsend = nsend + nstrains
      endif
      if (latomicstress) then
        do i = 1,numat
          do j = 1,nstrains
            nsend = nsend + 1
            sum0(nsend) = atomicstress(j,i)
          enddo
        enddo
      endif
      if (lgrad2) then
        do i = 1,nstrains
          do j = 1,nstrains
            nsend = nsend + 1
            sum0(nsend) = sderv2(j,i)
          enddo
        enddo
      endif
    endif
    call sumall(sum0,sum,nsend,"psumall","dervs")
    if (lstr) then
      if (lgrad2) then
        nsend = nsend - nstrains**2
        do i = 1,nstrains
          do j = 1,nstrains
            nsend = nsend + 1
            sderv2(j,i) = sum(nsend)
          enddo
        enddo
        nsend = nsend - nstrains**2
      endif
      if (latomicstress) then
        nsend = nsend - numat*nstrains
        do i = 1,numat
          do j = 1,nstrains
            nsend = nsend + 1
            atomicstress(j,i) = sum(nsend)
          enddo
        enddo
        nsend = nsend - numat*nstrains
      endif
      if (lMEAM) then
        nsend = nsend - nstrains
        do i = 1,nstrains
          rstrdnr(i) = sum(nsend+i)
        enddo
      endif
      nsend = nsend - nstrains
      do i = 1,nstrains
        rstrd(i) = sum(nsend+i)
      enddo
    endif
    do i = 1,n
      xdrv(i) = sum(i)
    enddo
    do i = 1,n
      ydrv(i) = sum(n+i)
    enddo
    do i = 1,n
      zdrv(i) = sum(2*n+i)
    enddo
    do i = 1,n
      siteenergy(i) = sum(3*n+i)
    enddo
    if (lMEAM) then
      do i = 1,n
        xdrvnr(i) = sum(4*n+i)
      enddo
      do i = 1,n
        ydrvnr(i) = sum(5*n+i)
      enddo
      do i = 1,n
        zdrvnr(i) = sum(6*n+i)
      enddo
      nrecv = 7*n
    else
      nrecv = 4*n
    endif
    do i = 1,maxregion
      xregdrv(i) = sum(nrecv+i)
    enddo
    nrecv = nrecv + maxregion
    do i = 1,maxregion
      yregdrv(i) = sum(nrecv+i)
    enddo
    nrecv = nrecv + maxregion
    do i = 1,maxregion
      zregdrv(i) = sum(nrecv+i)
    enddo
    nrecv = nrecv + maxregion
    if (nbsmat.gt.0) then
      do i = 1,n
        raderv(i) = sum(nrecv+i)
      enddo
    endif
  endif
  tsuml = g_cpu_time() - tsum0
  tsum = tsum + tsuml
!
!  Free local memory
!
  deallocate(sum0,stat=status)
  if (status/=0) call deallocate_error('psumall','sum0')
  deallocate(sum,stat=status)
  if (status/=0) call deallocate_error('psumall','sum')
#ifdef TRACE
  call trace_out('psumall')
#endif
!
  return
  end
