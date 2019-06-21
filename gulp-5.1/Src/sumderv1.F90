  subroutine sumderv1(n)
!
!  Completes first derivatives
!  Called by energy.
!
!   3/09 Created from sumderv2
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
  use current,        only : nstrains
  use derivatives
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in) :: n          ! Number of atoms in the asymmetric unit
!
!  Local variables
!
  integer(i4)             :: i          ! Looping index over atoms
#ifdef TRACE
  call trace_in('sumderv1')
#endif
!***********************************************
!  Sum over radial and non-radial derivatives  *
!***********************************************
  do i = 1,n
    xdrv(i) = xdrv(i) + xdrvnr(i)
    ydrv(i) = ydrv(i) + ydrvnr(i)
    zdrv(i) = zdrv(i) + zdrvnr(i)
  enddo
  if (lstr) then
    do i = 1,nstrains
      rstrd(i) = rstrd(i) + rstrdnr(i)
    enddo
  endif
#ifdef TRACE
  call trace_out('sumderv1')
#endif
!
  return
  end
