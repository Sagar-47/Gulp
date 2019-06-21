  subroutine oscillatorstrengthc(mcv,mcvloc,mcvptr,nphonatc,nphonatptr,ncfoc,iocptr,eigc,maxd2,oscstrength)
!
!  Calculates the oscillator strengths for the modes.
!  Complex version of oscillatorstrength
!
!  12/16 Created from oscillatorstrength
!   2/18 Trace added
!
!  On entry :
!
!  mcv         = no. of modes
!  mcvloc      = no. of modes on local node
!  mcvptr      = pointer from local to global mode
!  nphonatc    = total number of cores in relevant regions
!  nphonatptr  = pointer from reduce to full atom set
!  ncfoc       = number of condensed core sites
!  iocptr      = pointer that connects sites to condensed sites
!  eigc        = eigenvectors of dynamical matrix - complex
!  maxd2       = left-hand dimension of eigc
!
!  On exit : 
!
!  oscstrength = oscillator strengths for each mode as a 3 x 3 tensor per mode
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
  use current
  use element
  use iochannels
  use parallel,       only : nprocs
  use species,        only : massspec
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),                   intent(in)    :: iocptr(*)
  integer(i4),                   intent(in)    :: maxd2
  integer(i4),                   intent(in)    :: mcv
  integer(i4),                   intent(in)    :: mcvloc
  integer(i4),                   intent(in)    :: mcvptr(mcvloc)
  integer(i4),                   intent(in)    :: nphonatc
  integer(i4),                   intent(in)    :: nphonatptr(*)
  integer(i4),                   intent(in)    :: ncfoc
  complex(dpc),                  intent(in)    :: eigc(maxd2,*)
  real(dp),                      intent(out)   :: oscstrength(3,3,mcv)
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: j
  integer(i4)                                  :: m
  integer(i4)                                  :: mm
  integer(i4)                                  :: status
  complex(dpc)                                 :: osx
  complex(dpc)                                 :: osy
  complex(dpc)                                 :: osz
  real(dp)                                     :: trmj
  real(dp),     dimension(:,:,:), allocatable  :: sum3
#ifdef TRACE
  call trace_in('oscillatorstrengthc')
#endif
!**********************************
!  Oscillator strengths per mode  *
!**********************************
!
!  Zero arrays
!
  oscstrength(1:3,1:3,1:mcv) = 0.0_dp
!
!  Loop over modes
!
  do m = 1,mcvloc
    mm = mcvptr(m)
    osx = 0.0_dp
    osy = 0.0_dp
    osz = 0.0_dp
!
!  Loop over full sites
!
    do i = 1,ncfoc
      ix = 3*(i-1) + 1
      iy = ix + 1
      iz = ix + 2
!
!  Find all cores associated with full site
!
      do j = 1,nphonatc
        if (iocptr(j).eq.i) then
          trmj = occuf(j)/sqrt(massspec(nspecptr(nrelat(nphonatptr(j)))))
!
!  Multiple inverse mass weighted eigenvectors by Born charges
!
          osx = osx + (eigc(ix,m)*bornq(1,1,j) + &
                      eigc(iy,m)*bornq(2,1,j) + &
                      eigc(iz,m)*bornq(3,1,j))*trmj
          osy = osy + (eigc(ix,m)*bornq(1,2,j) + &
                      eigc(iy,m)*bornq(2,2,j) + &
                      eigc(iz,m)*bornq(3,2,j))*trmj
          osz = osz + (eigc(ix,m)*bornq(1,3,j) + &
                      eigc(iy,m)*bornq(2,3,j) + &
                      eigc(iz,m)*bornq(3,3,j))*trmj
        endif
      enddo
    enddo
    oscstrength(1,1,mm) = dble(conjg(osx)*osx)
    oscstrength(2,1,mm) = dble(conjg(osy)*osx)
    oscstrength(3,1,mm) = dble(conjg(osz)*osx)
    oscstrength(1,2,mm) = dble(conjg(osx)*osy)
    oscstrength(2,2,mm) = dble(conjg(osy)*osy)
    oscstrength(3,2,mm) = dble(conjg(osz)*osy)
    oscstrength(1,3,mm) = dble(conjg(osx)*osz)
    oscstrength(2,3,mm) = dble(conjg(osy)*osz)
    oscstrength(3,3,mm) = dble(conjg(osz)*osz)
  enddo
  if (nprocs.gt.1) then
!
!  Globalise oscillator strengths
!
    allocate(sum3(3,3,mcv),stat=status)
    if (status/=0) call outofmemory('oscillatorstrengthc','sum3')
!
    call sumall(oscstrength,sum3,9_i4*mcv,"oscillatorstrengthc","oscstrength")
    oscstrength(1:3,1:3,1:mcv) = sum3(1:3,1:3,1:mcv)
!
    deallocate(sum3,stat=status)
    if (status/=0) call deallocate_error('oscillatorstrengthc','sum3')
  endif
#ifdef TRACE
  call trace_out('oscillatorstrengthc')
#endif
!
  return
  end
