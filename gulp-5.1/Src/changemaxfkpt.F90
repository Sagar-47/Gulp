  subroutine changemaxfkpt
!
!  Alters the size of the arrays associated with maxfkpt
!
!   9/10 Created from changemaxkpt
!   8/14 eigv and groupvelocity added
!   8/16 IR intensity array added
!   1/18 Grueneisen parameters added
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
!  Julian Gale, CIC, Curtin University, January 2018
!
  use control,     only : lgroupvelocity, lgrueneisen
  use frequencies, only : freq, eigv, maxfkpt, maxfqat, lStoreEig, IRintensity
  use frequencies, only : groupvelocity, grueneisen
  use reallocate
  implicit none
!
  integer(i4) :: ierror
!
!  Frequencies
!
  call realloc(freq,3_i4*maxfqat,maxfkpt,ierror)
  if (ierror.ne.0) call outofmemory('changemaxfkpt','freq')
  call realloc(IRintensity,3_i4*maxfqat,maxfkpt,ierror)
  if (ierror.ne.0) call outofmemory('changemaxfkpt','IRintensity')
  if (lStoreEig) then
    call realloc(eigv,3_i4*maxfqat,3_i4*maxfqat,maxfkpt,ierror)
    if (ierror.ne.0) call outofmemory('changemaxfkpt','eigv')
  endif
  if (lgroupvelocity) then
    call realloc(groupvelocity,3_i4,3_i4*maxfqat,maxfkpt,ierror)
    if (ierror.ne.0) call outofmemory('changemaxfkpt','groupvelocity')
  endif
  if (lgrueneisen) then
    call realloc(grueneisen,3_i4*maxfqat,maxfkpt,ierror)
    if (ierror.ne.0) call outofmemory('changemaxfkpt','grueneisen')
  endif
!
  return
  end
