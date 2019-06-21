  subroutine groupvelocitiesd(nk,mcv,mcvloc,ncfoc,eigc,maxd2,rmass,fscale,lprint)
!
!  Calculates the group velocities for the modes at this k point.
!  Version for distributed memory parallelisation using complex eigenvectors.
!
!  12/16 Created from groupvelocities
!   8/17 Routine wrapped in ifdef since it won't work or be called without MPI
!   2/18 Trace added
!
!  On entry :
!
!  nk          = k point for storing group velocities
!  mcv         = no. of modes
!  mcvloc      = no. of modes on local node
!  ncfoc       = number of condensed core sites
!  eigc        = eigenvectors of dynamical matrix - complex
!  maxd2       = left-hand dimension of eigc
!  rmass       = inverse square root masses for occupied sites
!  fscale      = unit conversion from (eV/Ang**2/amu)^(1/2) to cm^-1
!  lprint      = if true then output group velocities
!
!  On exit :
!
!  groupvelocity = group velocities in cm^-1.Ang
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
#ifdef MPI
  use derivatives,    only : derv2dk
#endif
  use element
  use frequencies
  use iochannels
  use parallel
#ifdef MPI
  use phononatoms,    only : nphonatonnodec, nphonatonnodecptr
#endif
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),      intent(in)    :: maxd2
  integer(i4),      intent(in)    :: mcv
  integer(i4),      intent(in)    :: mcvloc
  integer(i4),      intent(in)    :: nk
  integer(i4),      intent(in)    :: ncfoc
  logical,          intent(in)    :: lprint
  complex(dpc),     intent(in)    :: eigc(maxd2,mcv)
  real(dp),         intent(in)    :: rmass(ncfoc)
  real(dp),         intent(in)    :: fscale
#ifdef MPI
!
!  Local variables
!
  integer(i4)                     :: i
  integer(i4)                     :: ii
  integer(i4)                     :: iloc
  integer(i4)                     :: indi
  integer(i4)                     :: indj
  integer(i4)                     :: j
  integer(i4)                     :: jj
  integer(i4)                     :: m
  integer(i4)                     :: status
  complex(dpc)                    :: cdot
  complex(dpc)                    :: cmassij
  complex(dpc), allocatable, save :: c2tmp(:,:)
  complex(dpc), allocatable, save :: cres(:)
  complex(dpc), external          :: zdotc
  real(dp)                        :: rtmp
!
  integer                         :: idesc(9)
  integer                         :: idese(9)
  integer                         :: idesv(9)
  integer                         :: ifails
  integer                         :: ld
  integer                         :: nb
  integer                         :: nc
#ifdef TRACE
  call trace_in('groupvelocitiesd')
#endif
!
!  Mass weight dynamical matrix derivatives w.r.t. k
!
  do iloc = 1,nphonatonnodec
    i = nphonatonnodecptr(iloc)
    do j = 1,ncfoc
      cmassij = dcmplx(rmass(i)*rmass(j),0.0_dp)
      do ii = 1,3
        indi = ii - 3
        do jj = 1,3
          indj = jj - 3
          derv2dk(1:3,indj+3*j,indi+3*iloc) = cmassij*derv2dk(1:3,indj+3*j,indi+3*iloc)
        enddo
      enddo
    enddo
  enddo
!
!  Zero group velocities
!
  groupvelocity(1:3,1:mcv,nk) = 0.0_dp
!
!  Create workspace arrays
!
  allocate(c2tmp(mcv,mcvloc),stat=status)
  if (status/=0) call outofmemory('groupvelocitiesd','c2tmp')
  allocate(cres(mcv),stat=status)
  if (status/=0) call outofmemory('groupvelocitiesd','cres')
!
!  Set up Blacs descriptors for matrices
!
  nb = nblocksize
  ifails = 0
  nc = mcv
  ld = mcv
  call descinit( idesc, nc, nc, 3*nb, 3*nb, 0, 0, iBlacsContext, ld, ifails )
!
  if (ifails.ne.0) then
    call outerror('initialisation in descinit failed',0_i4)
    call stopnow('groupvelocitiesd')
  endif
!
  call descinit( idesv, nc, 1, 3*nb, 3*nb, 0, 0, iBlacsContext, ld, ifails )
!
  if (ifails.ne.0) then
    call outerror('initialisation in descinit failed',0_i4)
    call stopnow('groupvelocitiesd')
  endif
!
  ld = maxd2
  call descinit( idese, nc, nc, 3*nb, 3*nb, 0, 0, iBlacsContext, ld, ifails )
!
  if (ifails.ne.0) then
    call outerror('initialisation in descinit failed',0_i4)
    call stopnow('groupvelocitiesd')
  endif
!
!  Loop over modes and project according to eigenvectors
!
  do ii = 1,3
    do i = 1,mcvloc
      do j = 1,mcv
        c2tmp(j,i) = derv2dk(ii,j,i)
      enddo
    enddo
    cres = 0.0_dpc
    do m = 1,mcv
      call pzgemv('N',mcv,mcv,(1.0d0,0.0d0),c2tmp,1,1,idesc,eigc,1,m,idese,1, &
                  (0.0d0,0.0d0),cres,1,1,idesv,1)
      call pzdotc(mcv,cdot,cres,1,1,idesv,1,eigc,1,m,idese,1)
      groupvelocity(ii,m,nk) = dble(cdot)
    enddo
  enddo
!
!  Free workspace arrays
!
  deallocate(cres,stat=status)
  if (status/=0) call deallocate_error('groupvelocitiesd','cres')
  deallocate(c2tmp,stat=status)
  if (status/=0) call deallocate_error('groupvelocitiesd','c2tmp')
!
!  Divide by conversion factors and 2 x frequency to convert derivatives
!  of frequency squared to those of frequency w.r.t. k
!
  do m = 1,mcv
    if (abs(freq(m,nk)).gt.1.0d-3) then
      rtmp = 0.5_dp*fscale*fscale/freq(m,nk)
      do ii = 1,3
        groupvelocity(ii,m,nk) = rtmp*groupvelocity(ii,m,nk)
      enddo
    else
      do ii = 1,3
        groupvelocity(ii,m,nk) = 0.0_dp
      enddo
    endif
  enddo
  if (lprint.and.ioproc) then
!
!  Output group velocities
!
    write(ioout,'('' Group velocities (cm-1.Ang) : '',/)')
    write(ioout,'('' Mode '',8x,''X'',13x,''Y'',13x,''Z''/)')
    do m = 1,mcv
      if (freq(m,nk).gt.0.0_dp) then
        write(ioout,'(i6,1x,3f14.6)') m,(groupvelocity(ii,m,nk),ii=1,3)
      else
        write(ioout,'(i6,1x,3f14.6,'' i'')') m,(groupvelocity(ii,m,nk),ii=1,3)
      endif
    enddo
    write(ioout,'(/)')
  endif
#else
  call outerror('groupvelocitiesd called when not compiled with MPI',0_i4)
  call stopnow('groupvelocitiesd')
#endif
#ifdef TRACE
  call trace_out('groupvelocitiesd')
#endif
!
  return
  end
