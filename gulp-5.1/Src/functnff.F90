  subroutine functnff(iflag,n,xc,fc,gc,lnonfitcall)
!
!  Supplies the function and derivatives where the second
!  derivatives are determined by finite differences. This
!  routine differs from functn in that it determines the 
!  full second derivative matrices rather than just the 
!  Hessian matrix elements for the optimisation variables.
!
!  Faster version that saves work for bond order potentials.
!  Parallelism here is over atoms in finite difference for
!  Cartesian internal-internal second derivatives
!
!  12/17 Created from functnf
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
  use bondorderdata,  only : nbopot
  use configurations, only : lbsmat
  use control
  use current,        only : numat, nstrains, xclat, yclat, zclat, ndim, radf, nbsmat
  use derivatives,    only : xdrv, ydrv, zdrv, strderv, derv2, derv3, sderv2, raderv
  use gulp_files
  use general
  use m_brenner,      only : setbrennerneigh, unsetbrennerneigh
  use m_bondorder,    only : setbondorderneigh, unsetbondorderneigh
  use m_edip,         only : setedipneigh, unsetedipneigh
  use m_reaxff,       only : setreaxFFneigh, unsetreaxFFneigh
  use parallel
  use symmetry,       only : lstr, lsymderv, lsymderv2
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),            intent(inout) :: iflag            ! Flag to indicate level of derivatives
  integer(i4),            intent(in)    :: n                ! Number of variables
  logical,                intent(in)    :: lnonfitcall      ! Flag to indicate whether this is a call from fitfun or not
  real(dp),               intent(out)   :: fc               ! Energy (on return)
  real(dp),               intent(inout) :: xc(*)            ! Structure in linear array
  real(dp),               intent(out)   :: gc(*)            ! Gradients in linear array (on return)
!
!  Local variables
!
  integer(i4)                           :: i                ! Looping variable over atoms
  integer(i4)                           :: ii               ! Looping variable over atoms on node
  integer(i4)                           :: indi             ! Index for Cartesian coordinate of i in derv2
  integer(i4)                           :: indif            ! Index for Cartesian coordinate of i in derv2
  integer(i4)                           :: ind3             ! 3*numat - used in breathing shell indexing in derv2
  integer(i4)                           :: indj             ! Index for Cartesian coordinate of i in derv2
  integer(i4)                           :: ix               ! Looping variable over Cartesian components for i
  integer(i4)                           :: j                ! Looping variable over atoms
  integer(i4)                           :: jx               ! Looping variable over Cartesian components for j
  integer(i4)                           :: status           ! Error flag for memory allocation
  logical                               :: lbreathing       ! Logical to indicate whether there are any breathing shells
  logical                               :: lgeometryOK
  logical                               :: lgrad1           ! Flag to indicate whether first derivatives are to be calculated
  logical                               :: lgrad2           ! Flag to indicate whether second derivatives are to be calculated
  logical                               :: lsymdervsave     ! Saves value of lsymderv
  logical                               :: lsymderv2save    ! Saves value of lsymderv2
  real(dp)                              :: g_cpu_time
  real(dp)                              :: dfindiffcorr     ! Term used in diagonal strain correction (1 - delta**2)
  real(dp)                              :: rfindiffc        ! 1/(2 x finite difference value for Cartesians)
  real(dp)                              :: rfindiffs        ! 1/(2 x finite difference value for strains)
  real(dp),                        save :: tdmax = 0.0_dp
  real(dp)                              :: t1
  real(dp)                              :: t2
  real(dp)                              :: xci              ! Stores the original value of the variable being finite differenced
  real(dp), dimension(:),   allocatable :: gbb              ! Workspace for backward radial derivatives
  real(dp), dimension(:),   allocatable :: gbf              ! Workspace for forward radial derivatives
  real(dp), dimension(:,:), allocatable :: gcb              ! Workspace for backward Cartesian derivatives
  real(dp), dimension(:,:), allocatable :: gcf              ! Workspace for forward Cartesian derivatives
  real(dp)                              :: gsb(6)           ! Workspace for backward strain derivatives
  real(dp)                              :: gsf(6)           ! Workspace for forward strain derivatives
  real(dp)                              :: x0str(6)         ! Strain workspace
#ifdef TRACE
  call trace_in('functnff')
#endif
!
  t1 = g_cpu_time()
  lgrad1 = (iflag.ge.1)
  lgrad2 = (iflag.ge.2)
  lbreathing = (nbsmat.gt.0) 
!
!  Store values of lsymderv(2) since these must be turned off during finite differencing as it breaks symmetry
!
  lsymdervsave = lsymderv
  lsymderv2save = lsymderv2
  lsymderv = .false.
  lsymderv2 = .false.
!************************************************************
!  Convert optimisation variables to linear structure array *
!************************************************************
  if (lnonfitcall) call xctox0(n,xc,lgeometryOK)
!************************************************************
!  Convert linear structure array to main structure arrays  *
!************************************************************
  if (lx0centroid) then
    call x0tostrcentroid
  else
    call x0tostr
  endif
!
  lfirst = .true.
  if (lgrad2) then
!
!  Allocate local arrays
!
    allocate(gcb(3,numat),stat=status)
    if (status/=0) call outofmemory('functnff','gcb')
    allocate(gcf(3,numat),stat=status)
    if (status/=0) call outofmemory('functnff','gcf')
    if (lbreathing) then
      allocate(gbb(numat),stat=status)
      if (status/=0) call outofmemory('functnff','gbb')
      allocate(gbf(numat),stat=status)
      if (status/=0) call outofmemory('functnff','gbf')
    endif
!
!  Initialise energy setup
!
    if (nbopot.gt.0) then
      call setbondorderneigh
    endif
    if (lEDIP) then
      call setedipneigh
    endif
    if (lreaxFF) then
      call setreaxFFneigh
    endif
    if (lbrenner) then
      call setbrennerneigh
    endif
!
!  Initialise second derivatives
!
    call initdervs(lgrad1,lgrad2)
!*********************************************************************************************
!  Compute second derivatives by finite differences with respect to all Cartesian components *
!*********************************************************************************************
    rfindiffc = 0.5_dp/findiffc
    rfindiffs = 0.5_dp/findiffs
    dfindiffcorr = 0.25_dp*findiffs**2/(1.0_dp - 0.25_dp*findiffs**2)
    ind3 = 3*numat
    do ii = 1,natomsonnode
      i = node2atom(ii)
      indi = 3*(ii-1)
      do ix = 1,3
!
!  Save coordinate
!
        if (ix.eq.1) then
          xci = xclat(i)
        elseif (ix.eq.2) then
          xci = yclat(i)
        else
          xci = zclat(i)
        endif
!
!  Forward step
!
        if (ix.eq.1) then
          xclat(i) = xci + findiffc
        elseif (ix.eq.2) then
          yclat(i) = xci + findiffc
        else
          zclat(i) = xci + findiffc
        endif
        call energyfd(i,ix,findiffc,fc,lgrad1)
        do j = 1,numat
          gcf(1,j) = xdrv(j)
          gcf(2,j) = ydrv(j)
          gcf(3,j) = zdrv(j)
        enddo
        if (lstr) then
          call strfin(.false.)
          gsf(1:nstrains) = strderv(1:nstrains)
        endif
        if (lbreathing) then
          gbf(1:numat) = raderv(1:numat)
        endif
!
!  Backward step
!
        if (ix.eq.1) then
          xclat(i) = xci - findiffc
        elseif (ix.eq.2) then
          yclat(i) = xci - findiffc
        else
          zclat(i) = xci - findiffc
        endif
        call energyfd(i,ix,-findiffc,fc,lgrad1)
        do j = 1,numat
          gcb(1,j) = xdrv(j)
          gcb(2,j) = ydrv(j)
          gcb(3,j) = zdrv(j)
        enddo
        if (lstr) then
          call strfin(.false.)
          gsb(1:nstrains) = strderv(1:nstrains)
        endif
        if (lbreathing) then
          gbb(1:numat) = raderv(1:numat)
        endif
!
!  Calculate Cartesian-Cartesian second derivatives
!
        do j = 1,numat
          indj = 3*(j-1)
          do jx = 1,3
            derv2(indj+jx,indi+ix) = rfindiffc*(gcf(jx,j) - gcb(jx,j))
          enddo
        enddo
        if (lstr) then
!
!  Calculate mixed Cartesian-strain second derivatives
!
          do jx = 1,nstrains
            derv3(indi+ix,jx) = rfindiffc*(gsf(jx) - gsb(jx))
          enddo
        endif
!
!  Restore coordinate
!
        if (ix.eq.1) then
          xclat(i) = xci
        elseif (ix.eq.2) then
          yclat(i) = xci
        else
          zclat(i) = xci
        endif
      enddo
    enddo
    if (lstr) then
!******************************************************************************************
!  Compute second derivatives by finite differences with respect to all strain components *
!******************************************************************************************
!
!  Set initial strain values
!
      x0str(1:nstrains) = 1.0_dp
!
      do ix = 1,nstrains
!
!  Save coordinate
!
        xci = x0str(ix)
!
!  Forward step
!
        x0str(ix) = xci + findiffs
        call x0strainf(x0str,lgeometryOK)
        if (lx0centroid) then
          call x0tostrcentroid
        else
          call x0tostr
        endif
!
        call energy(fc,lgrad1,.false.)
        call strfin(.false.)
!
        gsf(1:nstrains) = strderv(1:nstrains)
        do j = 1,numat
          gcf(1,j) = xdrv(j)
          gcf(2,j) = ydrv(j)
          gcf(3,j) = zdrv(j)
        enddo
        if (lbreathing) then
          gbf(1:numat) = raderv(1:numat)
        endif
!
!  Reset strain back to mid point
!
        if ((ndim.eq.3.and.ix.le.3).or.(ndim.eq.2.and.ix.le.2).or.ndim.eq.1) then
          x0str(ix) = xci/(xci + findiffs)
        else
          x0str(ix) = xci - findiffs
        endif
        call x0strainf(x0str,lgeometryOK)
!
!  If this is an off-diagonal strain then need to recorrect the diagonal
!
        if (ndim.eq.3.and.ix.gt.3) then
          x0str(ix) = xci
          if (ix.eq.4) then
            x0str(2) = xci + dfindiffcorr
            x0str(3) = xci + dfindiffcorr
          elseif (ix.eq.5) then
            x0str(1) = xci + dfindiffcorr
            x0str(3) = xci + dfindiffcorr
          elseif (ix.eq.6) then
            x0str(1) = xci + dfindiffcorr
            x0str(2) = xci + dfindiffcorr
          endif
          call x0strainf(x0str,lgeometryOK)
          if (ix.eq.4) then
            x0str(2) = xci
            x0str(3) = xci
          elseif (ix.eq.5) then
            x0str(1) = xci
            x0str(3) = xci
          elseif (ix.eq.6) then
            x0str(1) = xci
            x0str(2) = xci
          endif
        elseif (ndim.eq.2.and.ix.gt.2) then
          x0str(ix) = xci
          x0str(1) = xci + dfindiffcorr
          x0str(2) = xci + dfindiffcorr
          call x0strainf(x0str,lgeometryOK)
          if (lx0centroid) then
            call x0tostrcentroid
          else
            call x0tostr
          endif
          x0str(1) = xci
          x0str(2) = xci
        endif
!
!  Backward step
!
        x0str(ix) = xci - findiffs
        call x0strainf(x0str,lgeometryOK)
        if (lx0centroid) then
          call x0tostrcentroid
        else
          call x0tostr
        endif
!
        call energy(fc,lgrad1,.false.)
        call strfin(.false.)
!
        gsb(1:nstrains) = strderv(1:nstrains)
        do j = 1,numat
          gcb(1,j) = xdrv(j)
          gcb(2,j) = ydrv(j)
          gcb(3,j) = zdrv(j)
        enddo
        if (lbreathing) then
          gbb(1:numat) = raderv(1:numat)
        endif
!
!  Calculate strain-strain second derivatives
!
        do jx = 1,nstrains
          sderv2(jx,ix) = rfindiffs*(gsf(jx) - gsb(jx))
        enddo
!
        if (lbreathing) then
!
!  Calculate mixed strain-radial second derivatives
!
          do j = 1,numat
            derv3(ind3+j,ix) = rfindiffs*(gbf(j) - gbb(j))
          enddo
        endif
!
!  Reset strain back to mid point
!
        if ((ndim.eq.3.and.ix.le.3).or.(ndim.eq.2.and.ix.le.2).or.ndim.eq.1) then
          x0str(ix) = xci/(xci - findiffs)
        else
          x0str(ix) = xci + findiffs
        endif
        call x0strainf(x0str,lgeometryOK)
!
!  If this is an off-diagonal strain then need to recorrect the diagonal
!
        if (ndim.eq.3.and.ix.gt.3) then
          x0str(ix) = xci
          if (ix.eq.4) then
            x0str(2) = xci + dfindiffcorr
            x0str(3) = xci + dfindiffcorr
          elseif (ix.eq.5) then
            x0str(1) = xci + dfindiffcorr
            x0str(3) = xci + dfindiffcorr
          elseif (ix.eq.6) then
            x0str(1) = xci + dfindiffcorr
            x0str(2) = xci + dfindiffcorr
          endif
          call x0strainf(x0str,lgeometryOK)
          if (ix.eq.4) then
            x0str(2) = xci
            x0str(3) = xci
          elseif (ix.eq.5) then
            x0str(1) = xci
            x0str(3) = xci
          elseif (ix.eq.6) then
            x0str(1) = xci
            x0str(2) = xci
          endif
        elseif (ndim.eq.2.and.ix.gt.2) then
          x0str(ix) = xci
          x0str(1) = xci + dfindiffcorr
          x0str(2) = xci + dfindiffcorr
          call x0strainf(x0str,lgeometryOK)
          x0str(1) = xci
          x0str(2) = xci
        endif
!
!  Restore coordinate
!
        x0str(ix) = xci
      enddo
    endif
    if (lbreathing) then
!**********************************************************************************************
!  Compute second derivatives by finite differences with respect to all breathing shell radii *
!**********************************************************************************************
      do i = 1,numat
        if (lbsmat(i)) then
!
!  Save radius
!
          xci = radf(i)
!
!  Forward step
!
          radf(i) = xci + findiffc
          call energy(fc,lgrad1,.false.)
          do j = 1,numat
            gcf(1,j) = xdrv(j)
            gcf(2,j) = ydrv(j)
            gcf(3,j) = zdrv(j)
          enddo
          gbf(1:numat) = raderv(1:numat)
!
!  Backward step
!
          radf(i) = xci - findiffc
          call energy(fc,lgrad1,.false.)
          do j = 1,numat
            gcb(1,j) = xdrv(j)
            gcb(2,j) = ydrv(j)
            gcb(3,j) = zdrv(j)
          enddo
          gbb(1:numat) = raderv(1:numat)
!
!  Calculate mixed Cartesian-radial second derivatives
!
          do j = 1,numat
            indj = 3*(j-1)
            derv2(indj+1,ind3+i) = rfindiffc*(gcf(1,j) - gcb(1,j))
            derv2(indj+2,ind3+i) = rfindiffc*(gcf(2,j) - gcb(2,j))
            derv2(indj+3,ind3+i) = rfindiffc*(gcf(3,j) - gcb(3,j))
          enddo
!
!  Calculate radial-radial second derivatives
!
          do j = 1,numat
            derv2(ind3+j,ind3+i) = rfindiffc*(gbf(j) - gbb(j))
          enddo
!
!  Restore coordinate
!
          radf(i) = xci
        endif
      enddo
    endif
!****************************************************************************
!  Average duplicate elements in matrices to ensure that they are symmetric *
!****************************************************************************
    if (nprocs.eq.1) then
      do i = 2,numat
        indi = 3*(i-1)
        do j = 1,i-1
          indj = 3*(j-1)
          do ix = 1,3
            do jx = 1,3
              derv2(indj+jx,indi+ix) = 0.5_dp*(derv2(indj+jx,indi+ix) + derv2(indi+ix,indj+jx))
              derv2(indi+ix,indj+jx) = derv2(indj+jx,indi+ix)
            enddo
          enddo
        enddo
      enddo
    endif
    if (lstr) then
      do ix = 2,nstrains
        do jx = 1,ix-1
          sderv2(jx,ix) = 0.5_dp*(sderv2(jx,ix) + sderv2(ix,jx))
          sderv2(ix,jx) = sderv2(jx,ix)
        enddo
      enddo
    endif
    if (lbreathing) then
      do i = 1,numat
        if (lbsmat(i)) then
          do jx = 1,3*numat
            derv2(jx,ind3+i) = 0.5_dp*(derv2(jx,ind3+i) + derv2(ind3+i,jx))
            derv2(ind3+i,jx) = derv2(jx,ind3+i)
          enddo
        endif
      enddo
      do i = 2,numat
        if (lbsmat(i)) then
          do j = 1,i-1
            if (lbsmat(j)) then
              derv2(ind3+j,ind3+i) = 0.5_dp*(derv2(ind3+j,ind3+i) + derv2(ind3+i,ind3+j))
              derv2(ind3+i,ind3+j) = derv2(ind3+j,ind3+i)
            endif
          enddo
        endif
      enddo
    endif
!***********************************************************************************
!  Ensure that translational invariance is obeyed by setting the diagonal elements *
!***********************************************************************************
    do ii = 1,natomsonnode
      i = node2atom(ii)
      indi = 3*(ii-1)
      indif = 3*(i-1)
      do ix = 1,3
        do jx = 1,3
          derv2(indif+jx,indi+ix) = 0.0_dp
        enddo
      enddo
      do j = 1,i-1
        indj = 3*(j-1)
        do ix = 1,3
          do jx = 1,3
            derv2(indif+jx,indi+ix) = derv2(indif+jx,indi+ix) - derv2(indj+jx,indi+ix)
          enddo
        enddo
      enddo
      do j = i+1,numat
        indj = 3*(j-1)
        do ix = 1,3
          do jx = 1,3
            derv2(indif+jx,indi+ix) = derv2(indif+jx,indi+ix) - derv2(indj+jx,indi+ix)
          enddo
        enddo
      enddo
    enddo
!
!  Free local arrays
!
    if (lbreathing) then
      deallocate(gbf,stat=status)
      if (status/=0) call deallocate_error('functnff','gbf')
      deallocate(gbb,stat=status)
      if (status/=0) call deallocate_error('functnff','gbb')
    endif
    deallocate(gcf,stat=status)
    if (status/=0) call deallocate_error('functnff','gcf')
    deallocate(gcb,stat=status)
    if (status/=0) call deallocate_error('functnff','gcb')
!
!  Finalise energy terms
!
    if (lbrenner) then
      call unsetbrennerneigh
    endif
    if (lreaxFF) then
      call unsetreaxFFneigh
    endif
    if (lEDIP) then
      call unsetedipneigh
    endif
    if (nbopot.gt.0) then
      call unsetbondorderneigh
    endif
  endif
!***********************************************************************************************
!  Compute the first derivatives and associated quantities according to normal funct procedure *
!***********************************************************************************************
!
!  Reset structure
!
  if (lx0centroid) then
    call x0tostrcentroid
  else
    call x0tostr
  endif
!
!  Evaluate function and first derivatives
!
  call energy(fc,lgrad1,.false.)
!
!  For surface, get surface energy
!
  if (lseok) call surfaceenergy(fc)
!
!  Complete strain derivatives
!
  if (lstr) call strfin(.false.)
!
!  Output second derivatives 
!
  if (lgrad2) call outderv
!
!  Output energy and derivatives to a file
!
  if (ldrv.and.ioproc) call outdrv(fc,lgrad1,lgrad2)
!
!  Option to write out a .frc file for QMPOT
!
  if (lfrc.and.ioproc) call outfrc(fc,lgrad1,.false.)
!
!  First derivative handling
!
  if (lgrad1.and.lnonfitcall) then
    call getderv1(n,xc,gc,lgrad2,.false.)
  endif
!
!  Restore values of lsymderv(2) before call to obtain first derivatives
!
  lsymderv = lsymdervsave
  lsymderv2 = lsymderv2save
!
  t2 = g_cpu_time()
  if (t2-t1.gt.tdmax) tdmax = t2 - t1
  if (timmax.gt.0.0_dp) then
    if ((timmax-t2).lt.tdmax.and..not.lrelax) iflag = -1
  endif
  if (.not.lgeometryOK) iflag = -2
#ifdef TRACE
  call trace_out('functnff')
#endif
!
  return
  end
