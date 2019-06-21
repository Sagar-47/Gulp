  subroutine shellopt
!
!  Subroutine for optimisation of shell positions only.
!
!   5/14 Created from optim
!   6/14 lfreezeok set to false
!   1/17 Call to transmat modified for distributed memory
!   2/17 nmin removed from arguments to minimise
!   3/17 hessian matrix changed to 2-D matrix to accommodate parallel case
!   3/17 Modifications made to allow for new variable order in iopt
!   7/17 Call to setoptptr added
!   2/18 Trace added
!   4/18 Allocation of hessian corrected to use nvaronnode for right-hand side
!   6/18 Parallel handling of nvar corrected
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
!  Julian Gale, CIC, Curtin University, July 2018
!
  use configurations
  use control
  use current
  use element,       only : maxele
  use energies,      only : fcsave
  use gulp_files
  use fitting
  use general
  use iochannels
  use optimisation
  use parallel
  use reallocate
  use symmetry
  use times
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  use xcgc

  implicit none
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: icount
  integer(i4)                                  :: icx
  integer(i4)                                  :: icy
  integer(i4)                                  :: icz
  integer(i4)                                  :: ierror
  integer(i4)                                  :: ifail
  integer(i4)                                  :: ii
  integer(i4)                                  :: ind
  integer(i4)                                  :: indc
  integer(i4), dimension(:), allocatable       :: ioptsave
  integer(i4)                                  :: j
  integer(i4)                                  :: k
  integer(i4)                                  :: mvar
  integer(i4)                                  :: ncellsave
  integer(i4)                                  :: ncellmaxsave
  integer(i4)                                  :: ncellminsave
  integer(i4)                                  :: node
  integer(i4), dimension(:), allocatable       :: node2varsave
  integer(i4)                                  :: nvarsave
  integer(i4)                                  :: nvaronnodesave
  integer(i4), dimension(:), allocatable       :: nvar2localsave
  integer(i4), dimension(:), allocatable       :: nvar2nodesave
  integer(i4),                            save :: ncflast = 0
  integer(i4), dimension(:), allocatable       :: ncount
  integer(i4),                            save :: nhwords = 1
  integer(i4),                            save :: nhuwords = 1
  integer(i4)                                  :: status
  logical                                      :: ldothis
  logical,                                save :: lfirsttime = .true.
  logical,                                save :: lhess2D = .false.
  logical                                      :: lfound
  logical                                      :: lfreezeok
  logical                                      :: lgradloc
  logical                                      :: lgrad2loc
  logical                                      :: lharmloc
  logical                                      :: loptiloc
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: diff
  real(dp)                                     :: fc
  real(dp),    dimension(:,:), pointer,   save :: hess => null()
  real(dp)                                     :: pv
  real(dp)                                     :: time1
  real(dp)                                     :: xf
  real(dp)                                     :: yf
  real(dp)                                     :: zf
#ifdef TRACE
  call trace_in('shellopt')
#endif
!
!  Turn off printing during optimisation
!
  lopprt = .false. 
!
!  Save ncell since we don't want to optimise the cell here
!
  ncellsave = ncell
  ncellmaxsave = ncellmax
  ncellminsave = ncellmin
  ncell = 0
!
!  Save nvar and iopt in case more than shell optimisation was specified in the input
!
  nvarsave = nvar
  nvaronnodesave = nvaronnode
  allocate(ioptsave(nvar),stat=status)
  if (status/=0) call outofmemory('shellopt','ioptsave')
  allocate(node2varsave(nvar),stat=status)
  if (status/=0) call outofmemory('shellopt','node2varsave')
  allocate(nvar2localsave(nvar),stat=status)
  if (status/=0) call outofmemory('shellopt','nvar2localsave')
  allocate(nvar2nodesave(nvar),stat=status)
  if (status/=0) call outofmemory('shellopt','nvar2nodesave')
!
  ioptsave(1:nvar) = iopt(1:nvar)
  node2varsave(1:nvaronnode) = node2var(1:nvaronnode)
  nvar2localsave(1:nvar) = nvar2local(1:nvar)
  nvar2nodesave(1:nvar) = nvar2node(1:nvar)
!
!  Reset nvar and iopt based on shells only
!
  nvar = 0
  do i = 1,nvarsave
    ind = ioptsave(i)
!
!  Exclude strains
!
    if (ind.gt.nstrains) then
      ind = ind - (nstrains+1)
      ind = (ind/3) + 1
      if (iatn(ind).gt.maxele) then
        nvar = nvar + 1
        iopt(nvar) = ioptsave(i)
      endif
    endif
  enddo
!
  if (nprocs.gt.1) then
    nvaronnode = 0
    node = 0
    icount = 0
!
!  If block size hasn't been input then choose a value based on the number of atoms versus processors
!
    if (nblocksizevar.eq.0) then
      nblocksizevar = 3
    endif
!
    do i = 1,nvar
      icount = icount + 1
      nvar2node(i) = node
      if (node.eq.procid) then
        nvaronnode = nvaronnode + 1
        node2var(nvaronnode) = i
        nvar2local(i) = nvaronnode
      else
        nvar2local(i) = 0
      endif
      if (icount.eq.nblocksizevar) then
        icount = 0
        node = node + 1
        if (node.eq.nprocs) node = 0
      endif
    enddo
  else
    nvaronnode = nvar
    do i = 1,nvar
      nvar2node(i) = 0
      node2var(i) = i
      nvar2local(i) = i
    enddo
  endif
!
!  Nullify Hessian pointer and initialise with basic size
!
  if (lfirsttime) then
    lfirsttime = .false.
    if (nprocs.gt.1) lhess2D = .true.
    if (nprocs.gt.1) then
      call realloc(hess,nhwords,nhuwords,ierror)
      if (ierror.ne.0) call outofmemory('shellopt','hess')
    else
      call realloc(hess,nhwords,nhuwords,ierror)
      if (ierror.ne.0) call outofmemory('shellopt','hess')
    endif
  endif
!
  lfreezeok = .false.
  loptiloc = lopt
  lgradloc = lgrad
  lharmloc = lharmrelax
  loptsuccess = .false.
!
!  Set flag to indicate whether second derivatives and therefore tmat will ever be needed
!
  lgrad2loc = (.not.lconj.and..not.llbfgs.and..not.lunit)
  if (lminch.and.mintype.le.2) lgrad2loc = .true.
!
  fc = 0.0_dp
  pv = 0.0_dp
!
!  Transfer all coordinates to x0, including cores
!
  do i = 1,nasym
    call cart2frac(ndimen(ncf),xclat(i),yclat(i),zclat(i),rv,xf,yf,zf,icx,icy,icz)
    x0(3*i+nstrains-2) = xf - icx
    x0(3*i+nstrains-1) = yf - icy
    x0(3*i+nstrains)   = zf - icz
  enddo
!
!  Transfer breathing shell radii
!
  mvar = 3*nasym + nstrains
  if (nbsmat.gt.0) then
    do i = 1,nasym
      x0(mvar+i) = radcfg(nsft+i)
    enddo
  endif
!***********************
!  Set freezing flags  *
!***********************
!
!  Initialise to fixed
!
  do i = 1,nasym
    lopf(i) = .false.
  enddo
  if (lbulknoopt) then
    loptiloc = .false.
  elseif (nvar.gt.0) then
    do i = 1,nvar
      if (iopt(i).gt.nstrains) then
        xc(i) = x0(iopt(i))
      endif
    enddo
    do i = 1,nvar
      ind = iopt(i)
      ind = ind - (nstrains+1)
      ind = (ind/3) + 1
      lopf(ind) = .true.
!
!  Check for constrained atoms
!
      if (ncon.gt.0) then
        do j = 1,ncon
          if (iopt(i).eq.ncvar(j)) then
            indc = ncfix(j)
            if (indc.gt.mvar) then
              indc = indc - mvar
              lopf(indc) = .true.
            elseif (indc.gt.nstrains) then
              indc = indc - (nstrains+1)
              indc = (indc/3) + 1
              lopf(indc) = .true.
            endif
          endif
        enddo
      endif
    enddo
  elseif (nvar.eq.0.and.(lopt.or.lgrad.or.lharmrelax)) then
    nwarn = nwarn + 1
    loptiloc = .false.
    lgradloc = .false.
    lharmloc = .false.
  endif
  if ((loptiloc.or.lharmloc).and.(.not.lconj.or.lminch)) then
!*********************
!  Allocate hessian  *
!*********************
    if (lhess2D) then
!
!  2-D hessian
!
      if (llbfgs) then
        nhwords = nvar*(2*lmbfgsorder + 1) + 2*lmbfgsorder
!
!  Check if minimiser might change and need more memory
!
        if (lminch.and.mintype.le.4) then
          nhuwords = nvaronnode
        else
          nhuwords = 1_i4
        endif
      else
        nhwords = nvar
        nhuwords = nvaronnode
      endif
    else
!
!  1-D triangular half hessian
!
      if (llbfgs) then
        nhwords = nvar*(2*lmbfgsorder + 1) + 2*lmbfgsorder
!
!  Check if minimiser might change and need more memory
!
        if (lminch.and.mintype.le.4) then
          nhwords = max(nhwords,nvar*(nvar+1)/2)
        endif
      else
        nhwords = nvar*(nvar+1)/2
      endif
    endif
    call realloc(hess,nhwords,nhuwords,ierror)
    if (ierror.ne.0) call outofmemory('shellopt','hess')
  endif
  fcsave = fc
!********************************
!       Optimisation            *
!********************************
  if (loptiloc.or.lharmloc) then
    ifail = 1
    lfreeze = lfreezeok
!
!  Set pointers to atoms for optimisation after setting lfreeze
!
    call setoptptr(lfreeze)
!
    if (ioproc) call gflush(ioout)
!
!  Setup transformation matrix if needed
!
    if (ncf.ne.ncflast.and.lgrad2loc) then
      if (nprocs.gt.1) then
        call transmatd
      else
        call transmat
      endif
    endif
!
!  Start minimisation
!
    if (lharmloc) then
!
!  Implicit harmonic relaxation
!
      call harmonicrelax(xc,fc,gc,hess,nhwords,lhess2D,1_i4,.false.)
      ifail = 4
      loptsuccess = .true.
      if (ioproc) call outener
    else
!
!  Minimise static/free energy
!
      call minimise(xc,fc,gc,hess,nhwords,lhess2D,ifail,1_i4,.false.)
      loptsuccess = (ifail.eq.0)
    endif
!
!  Substitute parameters into place
!
    do i = 1,nstrains
      x0(i) = 1.0_dp
    enddo
    do i = 1,nvar
      x0(iopt(i)) = xc(i)
    enddo
!**********************
!  Apply constraints  *
!**********************
    if (ncon.gt.0) then
      do i = 1,ncon
        x0(ncfix(i)) = 0.0_dp
      enddo
      do i = 1,ncon
        x0(ncfix(i)) = x0(ncvar(i))*conco(i) + conadd(i) + x0(ncfix(i))
      enddo
!
!  Handle additive constraints for fractional coordinates - take nearest pair of images
!
      if (ndim.gt.0) then
        allocate(ncount(mvar),stat=status)
        if (status/=0) call outofmemory('shellopt','ncount')
        do i = 1,mvar
          ncount(i) = 0
        enddo
        do i = 1,ncon
          ii = ncfix(i)
          ncount(ii) = ncount(ii) + 1
        enddo
        do i = nstrains+1,mvar
!
!  Select only those coordinates which are fractional
!
          if (ndim.eq.3) then
            ldothis = .true.
          elseif (ndim.eq.2) then
            ldothis = (mod((i-nstrains),3_i4).ne.0)
          elseif (ndim.eq.1) then
            ldothis = (mod((i-nstrains),3_i4).eq.1)
          endif
          if (ncount(i).ge.2.and.ldothis) then
            lfound = .false.
            j = 0
            do while (.not.lfound.and.j.lt.ncon-1)
              j = j + 1
              if (ncfix(j).eq.i) then
                k = j
                do while (.not.lfound.and.k.lt.ncon) 
                  k = k + 1
                  lfound = (ncfix(k).eq.i)
                enddo
              endif
            enddo
            if (lfound) then
              diff = abs(x0(ncvar(j)) - x0(ncvar(k)))
              if (diff.gt.0.5_dp) then
                x0(i) = x0(i) + 0.5_dp
                x0(i) = mod(x0(i),1.0_dp)
              endif
            endif
          endif
        enddo
        deallocate(ncount,stat=status)
        if (status/=0) call deallocate_error('shellopt','ncount')
      endif
    endif
!****************************************
!  Return data to configuration arrays  *
!****************************************
!
!  Atomic positions
!
    if (ndim.ge.1.and.lmodco) then
      do i = 1,nasym
        xcfg(i+nsft) = mod(x0(3*i+(nstrains-2))+1.0_dp,1.0_dp)
      enddo
    else
      do i = 1,nasym
        xcfg(i+nsft) = x0(3*i+(nstrains-2))
      enddo
    endif
    if (ndim.ge.2.and.lmodco) then
      do i = 1,nasym
        ycfg(i+nsft) = mod(x0(3*i+(nstrains-1))+1.0_dp,1.0_dp)
      enddo
    else
      do i = 1,nasym
        ycfg(i+nsft) = x0(3*i+(nstrains-1))
      enddo
    endif
    if (ndim.eq.3.and.lmodco) then
      do i = 1,nasym
        zcfg(i+nsft) = mod(x0(3*i+nstrains)+1.0_dp,1.0_dp)
      enddo
    else
      do i = 1,nasym
        zcfg(i+nsft) = x0(3*i+nstrains)
      enddo
    endif
!
!  Radii
!
    if (nbsmat.gt.0) then
      do i = 1,nasym
        radcfg(i+nsft) = x0(mvar+i)
        rada(i) = x0(mvar+i)
      enddo
    endif
!
!  Copy configuration coordinates back to current arrays
!
    do i = 1,nasym
      xafrac(i) = xcfg(nsft+i)
      yafrac(i) = ycfg(nsft+i)
      zafrac(i) = zcfg(nsft+i)
    enddo
    if (lsymopt) then
      call equpos(.true.,.false.)
    else
      do i = 1,numat
        xfrac(i) = xafrac(i)
        yfrac(i) = yafrac(i)
        zfrac(i) = zafrac(i)
      enddo
    endif
  endif
!************************
!  End of optimisation  *
!************************
  ncflast = ncf
!
!  Timing
!
  time1 = g_cpu_time()
  time1 = time1 - time0
!
!  Restore nvar and iopt in case more than shell optimisation was specified in the input
!
  nvar = nvarsave
  nvaronnode = nvaronnodesave
!
  iopt(1:nvar) = ioptsave(1:nvar)
  node2var(1:nvaronnode) = node2varsave(1:nvaronnode)
  nvar2local(1:nvar) = nvar2localsave(1:nvar)
  nvar2node(1:nvar) = nvar2nodesave(1:nvar)
!
  deallocate(nvar2nodesave,stat=status)
  if (status/=0) call deallocate_error('shellopt','nvar2nodesave')
  deallocate(nvar2localsave,stat=status)
  if (status/=0) call deallocate_error('shellopt','nvar2localsave')
  deallocate(node2varsave,stat=status)
  if (status/=0) call deallocate_error('shellopt','node2varsave')
  deallocate(ioptsave,stat=status)
  if (status/=0) call deallocate_error('shellopt','ioptsave')
!
!  Reset ncell 
!
  ncell = ncellsave
  ncellmax = ncellmaxsave
  ncellmin = ncellminsave
#ifdef TRACE
  call trace_out('shellopt')
#endif
!
  return
  end
