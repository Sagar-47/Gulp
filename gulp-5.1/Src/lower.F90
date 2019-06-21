  subroutine lower(mcv,mcvrptr,freq,nphonatc,nphonatptr,ncfoc,iocptr,maxd2,eigr)
!
!  Lowers the symmetry of a system according to the imaginary eigenvalue eigenvectors.
!
!   3/02 Created from peigen/peigeng
!   7/02 Modified to allow for region 1 phonons only
!   7/02 Probable bug in assignment of shift atoms in xcfg corrected -
!        displacements were being applied to i not j
!   4/04 Lowering of shell position added
!   3/07 Gauss renamed to GULP_gauss
!   8/15 If lower has been applied then set flag to indicate this
!  12/16 mcvrptr added to arguments for the parallel case
!   1/18 Modified to allow for the present of ghost cells
!   2/18 Trace added
!   3/18 Tolerance on imaginary frequencies now can be set from input files
!  11/18 Try to ensure consistent choice of direction
!  11/18 Option added to switch direction of lowering
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
!  Julian Gale, CIC, Curtin University, November 2018
!
  use configurations, only : xcfg, ycfg, zcfg, llowered, nsuperghost
  use control
  use current
  use general,        only : lowerscale, frqtol, lowersign
  use iochannels
  use parallel
  use shells,         only : ncsptr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),      intent(in)   :: iocptr(*)
  integer(i4),      intent(in)   :: maxd2
  integer(i4),      intent(in)   :: mcv
  integer(i4),      intent(in)   :: mcvrptr(mcv)
  integer(i4),      intent(in)   :: ncfoc
  integer(i4),      intent(in)   :: nphonatc
  integer(i4),      intent(in)   :: nphonatptr(*)
  real(dp),         intent(in)   :: freq(*)
  real(dp),         intent(in)   :: eigr(maxd2,*)
!
!  Local variables
!
  integer(i4)                    :: i
  integer(i4)                    :: ixlarge
  integer(i4)                    :: iylarge
  integer(i4)                    :: izlarge
  integer(i4)                    :: iocj
  integer(i4)                    :: j
  integer(i4)                    :: jj
  integer(i4)                    :: jjs
  integer(i4)                    :: ind
  integer(i4)                    :: nimag
  integer(i4)                    :: nmode
  integer(i4)                    :: nghostcell
  integer(i4)                    :: status
  real(dp)                       :: rmat(3,4)
  real(dp)                       :: scale
  real(dp),    allocatable, save :: sum1(:)
  real(dp),    allocatable, save :: sum2(:)
  real(dp),    allocatable, save :: xshift(:)
  real(dp),    allocatable, save :: yshift(:)
  real(dp),    allocatable, save :: zshift(:)
#ifdef TRACE
  call trace_in('lower')
#endif
!*****************************************************
!  Lowering of symmetry according to imaginary modes  *
!******************************************************
  if (index(keyword,'lowe').ne.0) then
!
!  Find number of ghost cells
!
    nghostcell = nsuperghost(1,ncf)*nsuperghost(2,ncf)*nsuperghost(3,ncf)
!
!  Allocate arrays to hold coordinate shifts
!
    allocate(xshift(numat),stat=status)
    if (status/=0) call outofmemory('lower','xshift')
    allocate(yshift(numat),stat=status)
    if (status/=0) call outofmemory('lower','yshift')
    allocate(zshift(numat),stat=status)
    if (status/=0) call outofmemory('lower','zshift')
!
!  Zero shift arrays
!
    xshift(1:numat) = 0.0_dp
    yshift(1:numat) = 0.0_dp
    zshift(1:numat) = 0.0_dp
!
    nimag = 1
    do while (freq(nimag).lt.-frqtol.and.nimag.le.mcv)
      nmode = mcvrptr(nimag)
      if (nmode.ne.0) then
        ind = 0
        scale = lowerscale*abs(freq(nimag))
        do i = 1,ncfoc
          rmat(1,4) = scale*eigr(ind+1,nmode)
          rmat(2,4) = scale*eigr(ind+2,nmode)
          rmat(3,4) = scale*eigr(ind+3,nmode)
          if (ndim.eq.3) then
            do j = 1,3
              rmat(1,j) = rv(1,j)
              rmat(2,j) = rv(2,j)
              rmat(3,j) = rv(3,j)
            enddo
            call GULP_gauss(3_i4,3_i4,1_i4,rmat)
          elseif (ndim.eq.2) then
            do j = 1,2
              rmat(1,j) = rv(1,j)
              rmat(2,j) = rv(2,j)
            enddo
            call GULP_gauss(2_i4,3_i4,1_i4,rmat)
          elseif (ndim.eq.1) then
            rmat(1,4) = rv(1,1)
          endif
          do j = 1,nphonatc
            if (lghost) then
!
!  For ghost cells, equate the atom with it's original image
!
              iocj = (iocptr(j) - 1)/nghostcell + 1
            else
              iocj = iocptr(j)
            endif
            if (iocptr(j).eq.i) then
              jj = nphonatptr(j)
              xshift(jj) = xshift(jj) + rmat(1,4)
              yshift(jj) = yshift(jj) + rmat(2,4)
              zshift(jj) = zshift(jj) + rmat(3,4)
!
!  Move associated shells too
!
              jjs = ncsptr(jj)
              if (jjs.gt.0) then
                xshift(jjs) = xshift(jjs) + rmat(1,4)
                yshift(jjs) = yshift(jjs) + rmat(2,4)
                zshift(jjs) = zshift(jjs) + rmat(3,4)
              endif
            endif
          enddo
          ind = ind + 3
        enddo
!
!  End of condition on mode being local to this processor
!
      endif
      nimag = nimag + 1
    enddo
    nimag = nimag - 1
    if (nimag.gt.0) then
!
!  Only apply shifts if there were imaginary modes
!
      if (nprocs.gt.1) then
!
!  Allocate workspace arrays for communication
!
        allocate(sum1(3*numat),stat=status)
        if (status/=0) call outofmemory('lower','sum1')
        allocate(sum2(3*numat),stat=status)
        if (status/=0) call outofmemory('lower','sum2')
!
!  Transfer shifts to single array
!
        do i = 1,numat
          sum1(i) = xshift(i)
        enddo
        ind = numat
        do i = 1,numat
          sum1(ind+i) = yshift(i)
        enddo
        ind = 2*numat
        do i = 1,numat
          sum1(ind+i) = zshift(i)
        enddo
!
!  Globalise data
!
        call sumall(sum1,sum2,3_i4*numat,"lower","sum1")
!
!  Transfer global data back to shift arrays
!
        do i = 1,numat
          xshift(i) = sum2(i)
        enddo
        ind = numat
        do i = 1,numat
          yshift(i) = sum2(ind+i)
        enddo
        ind = 2*numat
        do i = 1,numat
          zshift(i) = sum2(ind+i)
        enddo
!
!  Free workspace arrays for communication
!
        deallocate(sum2,stat=status)
        if (status/=0) call deallocate_error('lower','sum2')
        deallocate(sum1,stat=status)
        if (status/=0) call deallocate_error('lower','sum1')
      endif
!
!  Try to ensure that a consistent choice of direction is made since this can be
!  subject to numerical noise and the arbitrary sign of the phonon eigenvectors
!
      ixlarge = 0
      i = 0
      do while (ixlarge.eq.0.and.i.lt.numat)
        i = i + 1
        if (abs(xshift(i)).gt.1.0d-3) then
          ixlarge = i
        endif
      enddo
      if (ixlarge.ne.0) then
        scale = sign(1.0_dp,xshift(ixlarge))
      else
        iylarge = 0
        i = 0
        do while (iylarge.eq.0.and.i.lt.numat)
          i = i + 1
          if (abs(yshift(i)).gt.1.0d-3) then
            iylarge = i
          endif
        enddo
        if (iylarge.ne.0) then
          scale = sign(1.0_dp,yshift(iylarge))
        else
          izlarge = 0
          i = 0
          do while (izlarge.eq.0.and.i.lt.numat)
            i = i + 1
            if (abs(zshift(i)).gt.1.0d-3) then
              izlarge = i
            endif
          enddo
          if (izlarge.ne.0) then
            scale = sign(1.0_dp,zshift(izlarge))
          else
            scale = 1.0_dp
          endif
        endif
      endif
!
!  Multiply scale by lowersign to change direction if requested
!
      scale = scale*lowersign
!
!  Add shifts to configuration arrays
!
      do i = 1,numat
        xcfg(nsft+i) = xcfg(nsft+i) + xshift(i)*scale
        ycfg(nsft+i) = ycfg(nsft+i) + yshift(i)*scale
        zcfg(nsft+i) = zcfg(nsft+i) + zshift(i)*scale
      enddo
    endif
!
!  Set flag to indicate whether lowering has occured
!
    if (nimag.ge.1) then
      llowered(ncf) = .true.
    endif
!
!  Output
!
    if (ioproc) then
      if (nimag.gt.1) then
        write(ioout,'(''  Symmetry lowered according to '',i3,'' imaginary mode eigenvectors'',/)') nimag
      elseif (nimag.eq.1) then
        write(ioout,'(''  Symmetry lowered according to one imaginary mode eigenvector'',/)')
      else
        write(ioout,'(''  No imaginary modes present - current symmetry is correct'',/)')
      endif
    endif
!
!  Free workspace arrays
!
    deallocate(zshift,stat=status)
    if (status/=0) call deallocate_error('lower','zshift')
    deallocate(yshift,stat=status)
    if (status/=0) call deallocate_error('lower','yshift')
    deallocate(xshift,stat=status)
    if (status/=0) call deallocate_error('lower','xshift')
  endif
#ifdef TRACE
  call trace_out('lower')
#endif
!
  return
  end
