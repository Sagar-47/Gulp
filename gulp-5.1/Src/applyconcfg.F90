  subroutine applyconcfg
!
!  Apply constraints to configuration arrays to re-symmetrise.
!  Currently needed for simul fitting.
!
!  12/17 Created from setup
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
  use control
  use configurations
  use current
  use datatypes
  use reallocate
  use shells
  use species
  use symmetry
#ifdef TRACE
  use trace,      only : trace_in, trace_out
#endif
  implicit none
!
!  Local variables
!
  integer(i4), dimension(:), allocatable       :: ncount
  integer(i4)                                  :: i
  integer(i4)                                  :: ii
  integer(i4)                                  :: j
  integer(i4)                                  :: k
  integer(i4)                                  :: mvar
  integer(i4)                                  :: nf
  integer(i4)                                  :: nf2
  integer(i4)                                  :: nv
  integer(i4)                                  :: nv2
  integer(i4)                                  :: nvj
  integer(i4)                                  :: nvj2
  integer(i4)                                  :: nvk
  integer(i4)                                  :: nvk2
  integer(i4)                                  :: status
  logical                                      :: lfound
  real(dp)                                     :: diff
  real(dp)                                     :: vcrd
  real(dp)                                     :: vcrdj
  real(dp)                                     :: vcrdk
#ifdef TRACE
  call trace_in('applyconcfg')
#endif
!
!  Find first atom shift
!
  nsft = 0
  if (ncf.gt.1) then
    do i = 1,ncf-1
      nsft = nsft + nascfg(i)
    enddo
  endif
!
!  Set dimensionality
!
  ndim = ndimen(ncf)
!
!  Set number of strains according to the dimensionality and pointer
!
  if (ndim.eq.3) then
    nstrains = 6
  elseif (ndim.eq.2) then
    nstrains = 3
  elseif (ndim.eq.1) then
    nstrains = 1
  else
    nstrains = 0
  endif
!
!  Transfer stored configuration data into working arrays
!
  nasym = nascfg(ncf)
!
!  Set up constraint pointers
!
  if (ncf.eq.ncfg) then
    if (ncontot.gt.0) ncon = ncontot + 1 - n1con(ncf)
  else
    ncon = n1con(ncf+1) - n1con(ncf)
  endif
  if (ncon.gt.maxcon) then
    maxcon = ncon
    call changemaxcon
  endif
  ncfst = n1con(ncf) - 1
!
!  Apply constraints
!
  mvar = 3*nasym + nstrains
  if (ncon.gt.0) then
    do j = 1,ncon
      ncfix(j) = ncfixcfg(ncfst+j)
      ncvar(j) = ncvarcfg(ncfst+j)
      conco(j) = concocfg(ncfst+j)
      conadd(j) = conaddcfg(ncfst+j)
    enddo
    do j = 1,ncon
      nf = ncfix(j)
      if (nf.gt.nstrains.and.nf.le.mvar) then
        nf = nf - nstrains
        nf2 = (nf - 1)/3 + 1
        nf = nf - 3*(nf2 - 1)
        if (nf.eq.1) then
          xcfg(nsft+nf2) = 0.0_dp
        elseif (nf.eq.2) then
          ycfg(nsft+nf2) = 0.0_dp
        else
          zcfg(nsft+nf2) = 0.0_dp
        endif
      endif
    enddo
    do j = 1,ncon
      nv = ncvar(j)
      if (nv.gt.nstrains.and.nf.le.mvar) then
        nv = nv - nstrains
        nv2 = (nv - 1)/3 + 1
        nv = nv - 3*(nv2 - 1)
        if (nv.eq.1) then
          vcrd = xcfg(nsft+nv2)
        elseif (nv.eq.2) then
          vcrd = ycfg(nsft+nv2)
        else
          vcrd = zcfg(nsft+nv2)
        endif
        nf = ncfix(j) - nstrains
        nf2 = (nf - 1)/3 + 1
        nf = nf - 3*(nf2 - 1)
        if (nf.eq.1) then
          xcfg(nsft+nf2) = vcrd*conco(j) + conadd(j) + xcfg(nsft+nf2)
        elseif (nf.eq.2) then
          ycfg(nsft+nf2) = vcrd*conco(j) + conadd(j) + ycfg(nsft+nf2)
        else
          zcfg(nsft+nf2) = vcrd*conco(j) + conadd(j) + zcfg(nsft+nf2)
        endif
      endif
    enddo
!
!  Handle additive constraints for fractional coordinates
!  - take nearest pair of images
!
    if (ndim.eq.3) then
      allocate(ncount(mvar),stat=status)
      if (status/=0) call outofmemory('applyconcfg','ncount')
      do i = 1,mvar
        ncount(i) = 0
      enddo
      do i = 1,ncon
        ii = ncfix(i)
        ncount(ii) = ncount(ii) + 1
      enddo
      do i = nstrains+1,mvar
        if (ncount(i).ge.2) then
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
            nvj = ncvar(j)
            nvj2 = (nvj - (nstrains - 2))/3
            nvj = nvj - nvj2*3 - (nstrains - 3)
            if (nvj.eq.1) then
              vcrdj = xcfg(nsft+nvj2)
            elseif (nvj.eq.2) then
              vcrdj = ycfg(nsft+nvj2)
            else
              vcrdj = zcfg(nsft+nvj2)
            endif
            nvk = ncvar(k)
            nvk2 = (nvk - (nstrains - 2))/3
            nvk = nvk - nvk2*3 - (nstrains - 3)
            if (nvk.eq.1) then
              vcrdk = xcfg(nsft+nvk2)
            elseif (nvk.eq.2) then
              vcrdk = ycfg(nsft+nvk2)
            else
              vcrdk = zcfg(nsft+nvk2)
            endif
            diff = abs(vcrdk - vcrdj)
            if (diff.gt.0.5_dp) then
              nf = ncfix(j)
              nf2 = (nf - (nstrains - 2))/3
              nf = nf - nf2*3 - (nstrains - 3)
              if (nf.eq.1) then
                xcfg(nsft+nf2) = xcfg(nsft+nf2) + 0.5_dp
                xcfg(nsft+nf2) = mod(xcfg(nsft+nf2),1.0_dp)
              elseif (nf.eq.2) then
                ycfg(nsft+nf2) = ycfg(nsft+nf2) + 0.5_dp
                ycfg(nsft+nf2) = mod(ycfg(nsft+nf2),1.0_dp)
              else
                zcfg(nsft+nf2) = zcfg(nsft+nf2) + 0.5_dp
                zcfg(nsft+nf2) = mod(zcfg(nsft+nf2),1.0_dp)
              endif
            endif
          endif
        endif
      enddo
      deallocate(ncount,stat=status)
      if (status/=0) call deallocate_error('applyconcfg','ncount')
    endif
  endif
#ifdef TRACE
  call trace_out('applyconcfg')
#endif
!
  return
  end
