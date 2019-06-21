  subroutine background(ebgd,emad,lgrad1,lgrad2)
!
!  Calculates neutralising background contribution to energy
!  for non-charge neutral solids.
!
!   7/97 Calculation of virial added
!   2/01 Modifications for general dimensionality added
!   6/01 Shifting of the energy for 2-D case added so that
!        defect energies can be calculated
!   5/02 1-D case added
!  11/04 Intent added
!   7/07 Calculation of virial activated for metadynamics
!   3/09 Case of Wolf sum energy trapped
!   9/11 Metadynamics internal code replaced with Plumed
!   9/11 Madelung correction added
!   4/12 Explicit virial calculation removed as no longer needed
!   5/12 Atomic stresses added
!   8/17 Modified for parallel second derivatives to avoid
!        double counting in sderv2
!   1/18 Trace added
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
  use g_constants
  use control,      only : lDoElectrostatics, lwolf, lmadelung, latomicstress
  use current
  use derivatives
  use kspace
  use mdlogic
  use parallel
  use qmedata,      only : maxloop
#ifdef TRACE
  use trace,        only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  real(dp), intent(out)    :: ebgd
  real(dp), intent(out)    :: emad
  logical,  intent(in)     :: lgrad1
  logical,  intent(in)     :: lgrad2
!
!  Local variables
!
  integer(i4)              :: i
  logical                  :: lgrad2loc
  real(dp)                 :: asstress
  real(dp)                 :: demad
  real(dp)                 :: ebgd2
  real(dp)                 :: erealq
  real(dp)                 :: erecipq
  real(dp)                 :: reta
  real(dp)                 :: rvol
  real(dp)                 :: vol
  real(dp)                 :: volume
#ifdef TRACE
  call trace_in('background')
#endif
!
!  Only compute second derivatives on first node otherwise values
!  will be duplicated
!
  lgrad2loc = (lgrad2.and.ioproc)
!
!  Check that electrostatics are being used
!
  emad = 0.0_dp
  if (lDoElectrostatics.and..not.lwolf) then
    if (ndim.eq.3) then
      vol = volume(rv)
      rvol = 1.0_dp/vol
      reta = 1.0_dp/eta
      ebgd = - 0.5_dp*pi*angstoev*(totalcharge**2)*rvol*reta
      if (lmadelung) then
!
!  Check the system is cubic
!
        if (a.ne.b.or.b.ne.c.or.alpha.ne.90.0_dp.or.beta.ne.90.0_dp.or.gamma.ne.90_dp) then
          call outerror('Madelung correction can only be applied to cubic systems',0_i4)
          call stopnow('background')
        endif
        emad = 0.5_dp*2.837297_dp*angstoev*(totalcharge**2)/a
        demad = - emad/a
      else
        demad = 0.0_dp
      endif
      if (lgrad1) then
        strderv(1) = strderv(1) - ebgd - demad
        strderv(2) = strderv(2) - ebgd - demad
        strderv(3) = strderv(3) - ebgd - demad
        if (latomicstress) then
          asstress = (ebgd + demad)/dble(numat)
          do i = 1,numat
            atomicstress(1,i) = atomicstress(1,i) - asstress
            atomicstress(2,i) = atomicstress(2,i) - asstress
            atomicstress(3,i) = atomicstress(3,i) - asstress
          enddo
        endif
        if (lgrad2loc) then
          sderv2(1,1) = sderv2(1,1) + ebgd + emad
          sderv2(2,2) = sderv2(2,2) + ebgd + emad
          sderv2(3,3) = sderv2(3,3) + ebgd + emad
          sderv2(2,1) = sderv2(2,1) + ebgd + emad
          sderv2(3,1) = sderv2(3,1) + ebgd + emad
          sderv2(3,2) = sderv2(3,2) + ebgd + emad
          ebgd2 = 0.5_dp*(ebgd + emad)
          sderv2(4,4) = sderv2(4,4) + ebgd2
          sderv2(5,5) = sderv2(5,5) + ebgd2
          sderv2(6,6) = sderv2(6,6) + ebgd2
        endif
      endif
    elseif (ndim.eq.2) then
      erecipq = 0.0_dp
      erealq = 0.0_dp
      call recip2Dq(erecipq,totalcharge,lgrad1,lgrad2loc)
      call real2Dq(erealq,erecipq,totalcharge,lgrad1,lgrad2loc)
      ebgd = erecipq + erealq
    elseif (ndim.eq.1) then
      call setmaxcell1D(maxloop(1))
      call real1Dq(ebgd,totalcharge,lgrad1,lgrad2loc)
    endif
  else
    ebgd = 0.0_dp
  endif
#ifdef TRACE
  call trace_out('background')
#endif
!
  return
  end
