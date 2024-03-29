  subroutine initmaxnboAdefaults(i)
!
!  Initialises the arrays associated with maxnboA
!
!   9/10 Created from changemax routine
!   1/14 lBOzrlA added
!   9/15 BOccoeffA made into a 2-D array
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
!  Copyright Curtin University 2015
!
!  Julian Gale, CIC, Curtin University, September 2015
!
  use bondorderdata
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in) :: i
!
!  Initialise defaults for new part of array
!
  if (i.ge.1.and.i.le.maxnboA) then
    BOccoeffA(1:5,i) = 0.0_dp
    BOecoeffA(i) = 1.0_dp
    BOhcoeffA(i) = 0.0_dp
    BOlcoeffA(i) = 0.0_dp
    BOmcoeffA(i) = 3.0_dp
    BOncoeffA(i) = 0.0_dp
    nBOtypeA(i) = 1
    lBOzrlA(i) = .false.
  endif
!
  return
  end
