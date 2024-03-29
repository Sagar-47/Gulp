  subroutine sendall(send,count,sendingnode,caller,desc)
!
!  Uses MPI_bcast to perform a global broadcast for "count" real*8 items
!

!
!  Modules
!
  use parallel
  implicit none
!
!  Passed variables
!
  integer(i4)      :: count
  integer(i4)      :: sendingnode
  real(dp)         :: send(count)
  character(len=*) :: caller, desc
!
!  Other variables
!
#ifdef MPI
  include 'mpif.h'
  integer*4           :: ierror,icode,lenstring,ic,i
  real*8, allocatable :: send_mpi(:)
  character error_string(MPI_max_error_string)
!
  if (nprocs.gt.1) then
!
!  Copy values to MPI compatible variables
!
    allocate(send_mpi(count))
    do i = 1,count
      send_mpi(i) = send(i)
    enddo
!
    call MPI_bcast(send_mpi,count,MPI_double_precision,sendingnode,MPI_comm_GULP,ierror)
!
!  Copy values back from MPI compatible variables
!
    if (procid.ne.sendingnode) then
      do i = 1,count
        send(i) = send_mpi(i)
      enddo
    endif
    deallocate(send_mpi)

!  Check for errors
    if (ierror.ne.MPI_success) then
      write(*,"(5a)")  &
       "Call to MPI_bcast failed in subroutine sendall,", &
       " called from subroutine ", caller, " for ", desc
      icode = ierror
      call MPI_error_string(icode,error_string,lenstring,ierror)
      write(*,"(a)") error_string 
      call MPI_abort(MPI_comm_GULP,1,ierror)
    endif
  endif
#else
#endif
  return
  end

  subroutine isendall(isend,count,sendingnode,caller,desc)
!
!  Uses MPI_bcast to perform a global broadcast for "count" integer*4 items
!

!
!  Modules
!
  use parallel
  implicit none
!
!  Passed variables
!
  integer(i4)      :: count
  integer(i4)      :: sendingnode
  integer(i4)      :: isend(count)
  character(len=*) :: caller, desc
!
!  Other variables
!
#ifdef MPI
  include 'mpif.h'
  integer*4              :: ierror,icode,lenstring,ic,i
  integer*4, allocatable :: isend_mpi(:)
  character error_string(MPI_max_error_string)
!
  if (nprocs.gt.1) then
!
!  Copy values to MPI compatible variables
!
    allocate(isend_mpi(count))
    do i = 1,count
      isend_mpi(i) = isend(i)
    enddo
!
    call MPI_bcast(isend_mpi,count,MPI_integer,sendingnode,MPI_comm_GULP,ierror)
!
!  Copy values back from MPI compatible variables
!
    if (procid.ne.sendingnode) then
      do i = 1,count
        isend(i) = isend_mpi(i)
      enddo
    endif
    deallocate(isend_mpi)

!  Check for errors
    if (ierror.ne.MPI_success) then
      write(*,"(5a)")  &
       "Call to MPI_bcast failed in subroutine isendall,", &
       " called from subroutine ", caller, " for ", desc
      icode = ierror
      call MPI_error_string(icode,error_string,lenstring,ierror)
      write(*,"(a)") error_string 
      call MPI_abort(MPI_comm_GULP,1,ierror)
    endif
  endif
#else
#endif
  return
  end
