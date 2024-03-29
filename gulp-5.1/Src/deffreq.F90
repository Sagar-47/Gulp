  subroutine deffreq(lprint,fc)
!
!  Calculates defect frequencies.
!
!  The block diagonal elements of derv2 must be stored, as these are
!  not recalculated in dynamic as the summing of off diagonals is no
!  longer applicable.
!
!  leigloc = local flag to indicate whether eigenvectors are 
!            to be generated for this configuration
!  fhenergy= Helmholtz free-energy
!  rtlnz   = R*T*ln(z), where z=partition function
!
!   3/97 Corrected for partially occupied sites
!   3/97 Temperature range introduced for frequency properties
!   4/01 Change to way number of imaginary modes are calculated
!        to try to prevent inclusion of rotations
!   5/02 K point pointer added to references to freq
!   8/02 Channel 51 now closed
!  11/02 Cartesian components of IR intensity added
!   3/03 Code tidied through use of call to compressd2
!  10/04 Eispack call replaced by lapack
!   7/05 Deallocations cleaned
!   5/06 Mass now uses species values
!  11/07 Unused variables removed
!   4/08 Option to generate second derivatives by finite differences added
!   7/08 Output cleaned up for eigenvector case
!  11/08 Modified to handle case where second derivatives are available but
!        finite differences are requested
!   4/09 Use of lfinitediff replaces testing of finite difference value
!   6/09 Modified to allow for the fact that leigen can be set internally
!   3/10 Temporary file on channel 51 deleted when closed.
!   8/10 Equipartition free energy added for clusters
!   8/10 nozero keyword added
!  12/10 Correction added to equipartition free energy 
!   4/11 Use of array mass/rmass replaced with local array massdef/rmassdef
!        to avoid out of bounds error.
!   4/11 Call to changemaxfreqat added
!   1/12 Logic corrected so that frequencies are not printed during fitting.
!   6/12 nobsmodeptr0 and nobsmode added as input.
!   7/12 fobsmodeover added to store overlap information
!   9/12 Removal of translation and rotations added via eckart keyword
!   9/12 Eckart transformation re-written to be in terms of (1-P) matrix
!  10/12 Option added to allow two atoms co-exist on a site with an
!        occupancy greater than one.
!   9/12 Lower option added for finite systems
!  10/13 Raman intensities added for clusters, including from Raman susceptibilities
!   2/14 Projected density of states for clusters corrected
!  12/14 nbsptr renamed to avoid conflict with module version
!  10/15 trmzp now initialised correctly for nozero case
!   7/16 Extra check added to find linear molecules based on the number of atoms
!   8/16 Use of disk channel 51 removed
!   9/16 nphonatptr no longer allocated or deallocated since it is in a module
!  12/16 mcvrptr added to pass to lower
!  12/16 mcvptr added to pass to oscillatorstrengthg
!  12/16 Made exclusively for defect case
!   5/17 nobsmode removed since this no longer used since clusters stopped using deffreq
!  10/17 fhenergy moved to energies module
!  11/17 Format statements tweaked
!   2/18 Trace added
!   3/18 dynam shortened to dyna in keyword check
!   4/18 leigloc now set here for llower and linten
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
!  Julian Gale, CIC, Curtin University, April 2018
!
  use g_constants
  use control
  use current
  use defects
  use derivatives, vectors => dervi
  use element
  use energies,    only : fhenergy
  use gulp_files
  use frequencies
  use general,     only : lfinitediff
  use genetic
  use iochannels
  use parallel
  use projectdos
  use shells
  use species
  use times
#ifdef TRACE
  use trace,       only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  logical,      intent(in)                     :: lprint        ! If true then output results
  real(dp),     intent(in)                     :: fc            ! Internal energy
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4), dimension(:),     allocatable   :: ibocptr
  integer(i4)                                  :: ifail
  integer(i4)                                  :: igroup
  integer(i4)                                  :: ii
  integer(i4)                                  :: imoddb
  integer(i4)                                  :: inat
  integer(i4)                                  :: ind
  integer(i4)                                  :: indi
  integer(i4)                                  :: indj
  integer(i4)                                  :: inert(3)
  integer(i4), dimension(:),     allocatable   :: iocptr
  integer(i4)                                  :: iresid
  integer(i4), dimension(:),     allocatable   :: itmp
  integer(i4)                                  :: itype
  integer(i4)                                  :: j
  integer(i4)                                  :: jj
  integer(i4)                                  :: job
  integer(i4)                                  :: k
  integer(i4)                                  :: l
  integer(i4)                                  :: m
  integer(i4)                                  :: maxlim
  integer(i4)                                  :: mcv
  integer(i4),  dimension(:),     allocatable  :: mcvptr
  integer(i4),  dimension(:),     allocatable  :: mcvrptr
  integer(i4)                                  :: mcvmax
  integer(i4)                                  :: mcvmin
  integer(i4)                                  :: mint
  integer(i4)                                  :: mis
  integer(i4)                                  :: mjs
  integer(i4)                                  :: msv
  integer(i4)                                  :: nati
  integer(i4)                                  :: natj
  integer(i4), dimension(:),     allocatable   :: ndbsptr
  integer(i4)                                  :: nbfoc
  integer(i4)                                  :: ndbs
  integer(i4)                                  :: ncfoc
  integer(i4)                                  :: ni
  integer(i4)                                  :: nimag
  integer(i4)                                  :: np
  integer(i4)                                  :: npc
  integer(i4)                                  :: npfirst
  integer(i4)                                  :: npi
  integer(i4)                                  :: npifirst
  integer(i4)                                  :: npilast
  integer(i4)                                  :: nplast
  integer(i4)                                  :: nproj
  integer(i4)                                  :: nr1
  integer(i4)                                  :: nri
  integer(i4)                                  :: nrj
  integer(i4)                                  :: nsfoc
  integer(i4)                                  :: nsi
  integer(i4)                                  :: nt
  integer(i4)                                  :: ntj
  integer(i4)                                  :: status
  logical                                      :: lbsmi
  logical                                      :: lbsmj
  logical                                      :: lcorei
  logical                                      :: lcorej
  logical                                      :: leigloc
  logical                                      :: lfound
  logical                                      :: lnozero
  logical                                      :: lpocc
  logical                                      :: lproj
  logical                                      :: lprinloc
  real(dp)                                     :: cmfact
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: cv
  real(dp)                                     :: cv2
  real(dp)                                     :: det(2)
  real(dp)                                     :: ent2
  real(dp)                                     :: entropy
  real(dp)                                     :: factor
  real(dp)                                     :: fe_equipartition
  real(dp)                                     :: freqmin
  real(dp)                                     :: fscale
  real(dp),    dimension(:),     allocatable   :: irx
  real(dp),    dimension(:),     allocatable   :: iry
  real(dp),    dimension(:),     allocatable   :: irz
  real(dp),    dimension(:),     allocatable   :: massdef
  real(dp)                                     :: oci
  real(dp)                                     :: qj
  real(dp)                                     :: rkt
  real(dp),    dimension(:),     allocatable   :: rmassdef
  real(dp)                                     :: rmassi
  real(dp)                                     :: rmode
  real(dp)                                     :: rsum
  real(dp)                                     :: rt
  real(dp),    dimension(:),     allocatable   :: rtmp2
  real(dp)                                     :: rtlnz
  real(dp)                                     :: s_equipartition
  real(dp)                                     :: sumr
  real(dp)                                     :: t1d
  real(dp)                                     :: t1i
  real(dp)                                     :: t1t
  real(dp)                                     :: t2d
  real(dp)                                     :: t2i
  real(dp)                                     :: t2t
  real(dp)                                     :: tem
  real(dp)                                     :: trm
  real(dp)                                     :: trm1
  real(dp)                                     :: trmcv
  real(dp)                                     :: trmen
  real(dp)                                     :: trmfe
  real(dp)                                     :: trmfe_eq
  real(dp)                                     :: trmj
  real(dp)                                     :: trms_eq
  real(dp)                                     :: trmzp
  real(dp)                                     :: w
  real(dp),    dimension(:),     allocatable   :: w1
  real(dp),    dimension(:),     allocatable   :: w2
  real(dp),    dimension(:),     allocatable   :: w3
  real(dp)                                     :: wr
  real(dp)                                     :: xd
  real(dp)                                     :: yd
  real(dp)                                     :: zd
  real(dp)                                     :: xi
  real(dp)                                     :: yi
  real(dp)                                     :: zi
  real(dp)                                     :: xir
  real(dp)                                     :: yir
  real(dp)                                     :: zir
  real(dp)                                     :: zpe
#ifdef TRACE
  call trace_in('deffreq')
#endif
!
  t1t = g_cpu_time()
  leigloc = leigen
  if (linten) leigloc = .true.
  if (llower) leigloc = .true.
  lnozero  = (index(keyword,'noze').ne.0)
  lproj = ((nprojcfg(ncf)-nprojdef(ncf)).gt.0) 
  if (lproj) leigloc = .true.
  if (.not.lprint) leigloc = .false.
  lprinloc = (lprint.and.ioproc)
  fscale = sqrt(1.0d23*evtoj*avogadro)
  fscale = fscale/(2.0_dp*pi*speedl)
!
!  If analytic second derivatives are not available (or finite differences are requested) then try finite differences
!
  if (lnoanald2.or.lfinitediff) then
    call dynamicn_defect
  endif
!
!  Set mode for projection of DOS
!
  imoddb = 2
!***********
!  Defect  *
!***********
!
!  Allocate pointer arrays
!
  allocate(ndbsptr(nreg1),stat=status)
  if (status/=0) call outofmemory('deffreq','ndbsptr')
  allocate(ibocptr(nreg1),stat=status)
  if (status/=0) call outofmemory('deffreq','ibocptr')
  allocate(iocptr(nreg1),stat=status)
  if (status/=0) call outofmemory('deffreq','iocptr')
!
  ndbs = 0
  do i = 1,nreg1
    if (ldefbsmat(i)) then
      ndbs = ndbs + 1
      ndbsptr(ndbs) = i
    endif
  enddo
!
!  Setup partial occupancy pointer
!
  ncfoc = 0
  nsfoc = 0
  nbfoc = 0
  do i = 1,nreg1
    ibocptr(i) = 0
    if (occdefe(i).eq.1.0_dp.and..not.lallowgt1) then
!
!  Fully occupied site
!
      if (natdefe(i).gt.maxele) then
        nsfoc = nsfoc + 1
      else
        ncfoc = ncfoc + 1
      endif
      iocptr(i) = ncfoc + nsfoc
      if (ldefbsmat(i)) then
        nbfoc = nbfoc + 1
        ibocptr(i) = nbfoc
      endif
    else
!
!  Partially occupied site
!  Check to see if there is a previous atom on this site
!
      xi = xdefe(i)
      yi = ydefe(i)
      zi = zdefe(i)
      nati = natdefe(i)
      lcorei = (nati.le.maxele)
      lbsmi = (ldefbsmat(i))
      j = 1
      lfound = .false.
      do while (j.lt.i.and..not.lfound)
        lcorej = (natdefe(j).le.maxele)
        if ((lcorei.and.lcorej).or.(.not.lcorei.and..not.lcorej)) then
          xd = xdefe(j) - xi
          yd = ydefe(j) - yi
          zd = zdefe(j) - zi
          rsum = abs(xd) + abs(yd) + abs(zd)
          if (rsum.lt.1.0d-4) then
            lfound = .true.
            iocptr(i) = iocptr(j)
            lbsmj = (ldefbsmat(j))
            if (lbsmi) then
              if (lbsmj) then
                ibocptr(i) = ibocptr(j)
              else
                nbfoc = nbfoc + 1
                ibocptr(i) = nbfoc
              endif
            endif
          endif
        endif
        j = j + 1
      enddo
      if (.not.lfound) then
!
!  Must be new site
!
        if (lcorei) then
          ncfoc = ncfoc + 1
        else
          nsfoc = nsfoc + 1
        endif
        iocptr(i) = ncfoc + nsfoc
        if (ldefbsmat(i)) then
          nbfoc = nbfoc + 1
          ibocptr(i) = nbfoc
        endif
      endif
    endif
  enddo
  lpocc = (nsfoc+ncfoc.ne.nreg1)
!
!  End of partial occupancy pointer generation
!
  mint = 3*nreg1
  maxlim = mint
  if (ndbs.gt.0) maxlim = maxlim + nreg1
  nr1 = nreg1
  msv = 3*nshreg1 + ndbs
  msv = 3*nsfoc + nbfoc
  mcv = 3*ncfoc
!
!  Ensure that frequency array is large enough
!
  call changemaxfreqat(ncfoc)
!
!  Allocate local memory
!
  allocate(mcvptr(mcv),stat=status)
  if (status/=0) call outofmemory('deffreq','mcvptr')
  allocate(mcvrptr(mcv),stat=status)
  if (status/=0) call outofmemory('deffreq','mcvrptr')
  allocate(irx(mcv),stat=status)
  if (status/=0) call outofmemory('deffreq','irx')
  allocate(iry(mcv),stat=status)
  if (status/=0) call outofmemory('deffreq','iry')
  allocate(irz(mcv),stat=status)
  if (status/=0) call outofmemory('deffreq','irz')
  allocate(rtmp2(mcv),stat=status)
  if (status/=0) call outofmemory('deffreq','rtmp2')
  allocate(w1(3*mcv),stat=status)
  if (status/=0) call outofmemory('deffreq','w1')
  allocate(w2(3*mcv),stat=status)
  if (status/=0) call outofmemory('deffreq','w2')
  allocate(w3(mcv),stat=status)
  if (status/=0) call outofmemory('deffreq','w3')
  allocate(massdef(ncfoc),stat=status)
  if (status/=0) call outofmemory('deffreq','massdef')
  allocate(rmassdef(ncfoc),stat=status)
  if (status/=0) call outofmemory('deffreq','rmassdef')
!
  do m = 1,mcv
    mcvptr(m) = m 
    mcvrptr(m) = m 
  enddo
!
!  Calculate inversion square root of masses
!
!  Now modified to handle partial occupancies
!
  do i = 1,ncfoc
    massdef(i) = 0.0_dp
  enddo
  do i = 1,ncoreg1
    ni = natdefe(i)
    nt = ntypdefe(i)
    oci = occdefe(i)
    ii = iocptr(i)
    lfound = .false.
    nsi = 0
    do while (nsi.le.nspec.and..not.lfound)
      nsi = nsi + 1
      lfound = (ni.eq.natspec(nsi).and.nt.eq.ntypspec(nsi))
    enddo
    if (lfound) then
      rmassi = massspec(nsi)*oci
    else
      rmassi = atmass(ni)*oci
    endif
    if (rmassi.eq.0.0) then
      call outerror('mass of element '//atsym(ni)//' is zero',0_i4)
      call stopnow('deffreq')
    endif
    massdef(ii) = massdef(ii) + rmassi
  enddo
  do i = 1,ncfoc
    if (massdef(i).eq.0.0_dp) then
      call outerror('site has total mass of zero in phonon',0_i4)
      call stopnow('deffreq')
    endif
    rmassdef(i) = 1.0_dp/sqrt(massdef(i))
  enddo
!
!  Output frequency header
!
  if (lprinloc) then
    write(ioout,'(/,''  Vibrational Frequency Calculation : '',/)')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
  endif
!*****************************************************************
!  Compress second derivatives according to partial occupancies  *
!*****************************************************************
  if (lpocc) then
    call compressd2(derv2,maxd2,ncfoc,nsfoc,nbfoc,nr1,iocptr,ibocptr) 
  endif
!**************************************************************
!  Compress second derivative matrix w.r.t. breathing shells  *
!**************************************************************
  if (.not.lpocc.and.ndbs.gt.0) then
!
!  Reduce storage of breathing shell data
!
!  Pass 1 : reduce numat to ncfoc+nsfoc & numat to nbfoc in 1D
!
!  Full occupancy
!
    mis = mint
    mjs = mint
    do i = 1,ndbs
      nri = ndbsptr(i)
      do j = 1,ndbs
        nrj = ndbsptr(j)
        derv2(mjs+j,mis+i) = derv2(mjs+nrj,mis+nri)
      enddo
      do j = 1,nr1
        indj = 3*(j-1)
        derv2(indj+1,mis+i) = derv2(indj+1,mis+nri)
        derv2(indj+2,mis+i) = derv2(indj+2,mis+nri)
        derv2(indj+3,mis+i) = derv2(indj+3,mis+nri)
        derv2(mis+i,indj+1) = derv2(mis+nri,indj+1)
        derv2(mis+i,indj+2) = derv2(mis+nri,indj+2)
        derv2(mis+i,indj+3) = derv2(mis+nri,indj+3)
      enddo
    enddo
  endif
!**********************************
!  Eliminate shell contributions  *
!**********************************
  if (msv.gt.0) then
!
!  Invert dynamical shell-shell matrix
!
    do i = 1,msv
      do j = 1,msv
        vectors(j,i) = derv2(mcv+j,mcv+i)
      enddo
    enddo
    job = 1
    t1i = g_cpu_time()
    allocate(itmp(msv),stat=status)
    if (status/=0) call outofmemory('deffreq','itmp')
    call dsifa(vectors,maxd2,msv,itmp,ifail)
    if (ifail.ne.0) then
      call outerror('inversion of shell 2nd derivatives failed',0_i4)
      goto 999
    endif
    call dsidi(vectors,maxd2,msv,itmp,det,inert,w1,job)
    deallocate(itmp,stat=status)
    if (status/=0) call deallocate_error('deffreq','itmp')
    t2i = g_cpu_time()
    tmati = tmati + t2i - t1i
!
!  Resymmetrise and return to derv2
!
    do i = 1,msv
      do j = 1,i
        derv2(mcv+j,mcv+i) = vectors(j,i)
        derv2(mcv+i,mcv+j) = vectors(j,i)
      enddo
    enddo
!***********************************************
!  Corrected second derivatives  =  R - T*S-1*T  *
!***********************************************
!
!  First pass : S-1*T stored in diagonally opposite copy of T
!
    do i = 1,mcv
      do j = 1,msv
        wr = 0.0_dp
        do l = 1,msv
          wr = wr + derv2(mcv+j,mcv+l)*derv2(i,mcv+l)
        enddo
        derv2(mcv+j,i) = wr
      enddo
    enddo
!
!  Second pass : T*(S-1*T) 
!
    do i = 1,mcv
      do j = 1,mcv
        wr = 0.0_dp
        do l = 1,msv
          wr = wr - derv2(j,mcv+l)*derv2(mcv+l,i)
        enddo
        derv2(j,i) = derv2(j,i) + wr
      enddo
    enddo
  endif
!****************************
!  End of shell correction  *
!****************************
!
!  Multiply by mass-factors
!
  do ii = 1,3
    indi = ii - 3
    do jj = 1,3
      indj = jj - 3
      do i = 1,ncfoc
        rmassi = rmassdef(i)
        do j = 1,ncfoc
          derv2(indj+3*j,indi+3*i) = rmassi*rmassdef(j)*derv2(indj+3*j,indi+3*i)
        enddo
      enddo
    enddo
  enddo
!
!  If debugging print out dynamical matrix
!
  if (index(keyword,'dyna').ne.0.and.ioproc) then
    write(ioout,'(/,''  Real Dynamical matrix :'',/)')
    do i = 1,mcv
      write(ioout,'(12f10.5)')(derv2(j,i),j = 1,mcv)
    enddo
  endif
!*********************************
!  Diagonalise dynamical matrix  *
!*********************************
  ifail = 0
  t1d = g_cpu_time()
  do i = 1,mcv
    do j = 1,mcv
      vectors(j,i) = derv2(j,i)
    enddo
  enddo
  if (leigloc) then
    call dsyev('V','U',mcv,vectors,maxd2,freq,w1,3_i4*mcv,ifail)
  else
    call dsyev('N','U',mcv,vectors,maxd2,freq,w1,3_i4*mcv,ifail)
  endif
  t2d = g_cpu_time()
  tdiag = tdiag + t2d - t1d
!
!  Convert frequency units - imaginary freqs denoted by negative no.
!
  do i = 1,mcv
    rt = freq(i,1)
    if (rt.ge.0.0_dp) then
      freq(i,1) = sqrt(rt)*fscale
    else
      rt = abs(rt)
      freq(i,1) = - sqrt(rt)*fscale
    endif
  enddo
!***********************
!  Output frequencies  *
!***********************
  if (leigloc) then
!
!  Normalise eigenvectors
!
    do i = 1,mcv
      sumr = 0.0_dp
      do j = 1,mcv
        sumr = sumr + vectors(j,i)*vectors(j,i)
      enddo
      if (sumr.gt.0.0_dp) then
        sumr = sqrt(sumr)
        rsum = 1.0_dp/sumr
        do j = 1,mcv
          vectors(j,i) = rsum*vectors(j,i)
        enddo
      endif
    enddo
!**********************************
!  Projected densities of states  *
!**********************************
    if (lproj) then
      allocate(itmp(nr1),stat=status)
      if (status/=0) call outofmemory('deffreq','itmp')
!
!  Loop over projections
!
      nproj = nprojcfg(ncf)
      npfirst = 1
      npifirst = 1
      ii = 0
      do i = 1,ncf-1
        npc = nprojcfg(i)
        npfirst = npfirst + npc
        do j = 1,npc
          npifirst = npifirst + nprojit(ii+j)
        enddo
        ii = ii + npc
      enddo
      nplast = npfirst + nproj - 1
      npilast = npifirst
      do i = 1,nproj
        npilast = npilast + nprojit(ii+i)
      enddo
      npilast = npilast - 1
      do np = npfirst,nplast
        if (nprojdb(np).eq.imoddb) then
          do i = 1,nr1
            itmp(i) = 0
          enddo
!
!  Find atoms of projection
!
          do npi = npifirst,npilast
            if (nprojptr(npi).eq.np) then
              if (nprojtyp(npi).gt.99) then
                itmp(nprojnat(npi)) = 1
              else
                inat = nprojnat(npi)
                itype = nprojtyp(npi)
                do i = 1,nr1
                  if (inat.eq.natdefe(i).and.(itype.eq.ntypdefe(i).or.itype.eq.0)) itmp(i) = 1
                enddo
              endif
            endif
          enddo
!
!  Loop over frequencies
!
          do j = 1,mcv
            w = 0.0_dp
            ind = 0
!
!  Loop over atoms
!
            do k = 1,ncfoc
              lfound = .false.
              l = 1
              do while (l.le.ncoreg1.and..not.lfound)
                if (iocptr(l).eq.k) then
                  if (itmp(l).eq.1) then
                    lfound = .true.
                    w = w + vectors(ind+1,j)**2 + vectors(ind+2,j)**2 + vectors(ind+3,j)**2
                  endif
                endif
                l = l + 1
              enddo
              ind = ind + 3
            enddo
            if (ioproc) write(59) w
          enddo
        endif
      enddo
      deallocate(itmp,stat=status)
      if (status/=0) call deallocate_error('deffreq','itmp')
    endif
!*********************************************
!  Evaluate infra-red and Raman intensities  *
!*********************************************
    if (linten) then
!
!  Old algorithm - no use of oscillator strengths
!
      do i = 1,ncfoc
        w3(i) = 0.0_dp
        do j = 1,ncoreg1
          if (iocptr(j).eq.i) then
            qj   = qa(j)
            natj = natdefe(j)
            ntj  = ntypdefe(j)
            trmj = occua(j)/sqrt(atmass(natj))
!
!  Add on any shell charge
!
            do k = 1,nspec
              if ((natspec(k)-maxele).eq.natj) then
                if (ntj.eq.ntypspec(k).or.ntypspec(k).eq.0) then
                  qj = qj + qlspec(k)
                endif
              endif
            enddo
!
!  Term for IR
!
            w3(i) = w3(i) + qj*trmj
          endif
        enddo
      enddo
!
!  Sum eigenvector components multiplied by charge
!
      rkt = boltz*temperature
      if (abs(rkt).gt.1.0d-12) then
        cmfact = planck*speedl/rkt
      endif
      do i = 1,mcv
        xir = 0.0_dp
        yir = 0.0_dp
        zir = 0.0_dp
        ind = 0
        do j = 1,ncfoc
          xir = xir + w3(j)*vectors(ind+1,i)
          yir = yir + w3(j)*vectors(ind+2,i)
          zir = zir + w3(j)*vectors(ind+3,i)
          ind = ind + 3
        enddo
        IRintensity(i,1) = xir*xir + yir*yir + zir*zir
        irx(i) = xir*xir
        iry(i) = yir*yir
        irz(i) = zir*zir
      enddo
    endif
!************************
!  Output eigenvectors  *
!************************
    if (index(keyword,'eige').ne.0.and.ioproc) then
      write(ioout,'(/,''  Frequencies (cm-1) and Eigenvectors : '',/)')
      if (ncfoc.ne.ncoreg1) then
        write(ioout,'(''  Note: eigenvectors in terms of reduced sites due to partial occupancies!'',/)')
      endif
      igroup = mcv/6
      iresid = mcv - igroup*6
      indi = 0
      if (linten) then
        if (igroup.gt.0) then
          do i = 1,igroup
            write(ioout,'(''  Frequency   '',6f10.4)') (freq(indi+j,1),j = 1,6)
            write(ioout,'(''  IR Intensity'',6f10.4)') (IRintensity(indi+j,1),j = 1,6)
            write(ioout,'(''     in X     '',6f10.4)') (irx(indi+j),j = 1,6)
            write(ioout,'(''     in Y     '',6f10.4)') (iry(indi+j),j = 1,6)
            write(ioout,'(''     in Z     '',6f10.4,/)') (irz(indi+j),j = 1,6)
            indj = 0
            do j = 1,ncfoc
              write(ioout,'(i6,'' x '',4x,6f10.6)') j,(vectors(indj+1,indi+k),k = 1,6)
              write(ioout,'(i6,'' y '',4x,6f10.6)') j,(vectors(indj+2,indi+k),k = 1,6)
              write(ioout,'(i6,'' z '',4x,6f10.6)') j,(vectors(indj+3,indi+k),k = 1,6)
              indj = indj + 3
            enddo
            indi = indi + 6
            write(ioout,'(/)')
          enddo
        endif
        if (iresid.gt.0) then
          write(ioout,'(''  Frequency   '',6f10.4)') (freq(indi+j,1),j = 1,iresid)
          write(ioout,'(''  IR Intensity'',6f10.4)') (IRintensity(indi+j,1),j = 1,iresid)
          write(ioout,'(''     in X     '',6f10.4)') (irx(indi+j),j = 1,iresid)
          write(ioout,'(''     in Y     '',6f10.4)') (iry(indi+j),j = 1,iresid)
          write(ioout,'(''     in Z     '',6f10.4)') (irz(indi+j),j = 1,iresid)
          write(ioout,'(/)',advance='no')
          indj = 0
          do j = 1,ncfoc
            write(ioout,'(i6,'' x '',4x,6f10.6)') j,(vectors(indj+1,indi+k),k = 1,iresid)
            write(ioout,'(i6,'' y '',4x,6f10.6)') j,(vectors(indj+2,indi+k),k = 1,iresid)
            write(ioout,'(i6,'' z '',4x,6f10.6)') j,(vectors(indj+3,indi+k),k = 1,iresid)
            indj = indj + 3
          enddo
          write(ioout,'(/)')
        endif
      else
        if (igroup.gt.0) then
          do i = 1,igroup
            write(ioout,'('' Frequency   '',6f10.4,/)') (freq(indi+j,1),j = 1,6)
            indj = 0
            do j = 1,ncfoc
              write(ioout,'(i6,'' x '',4x,6f10.6)') j,(vectors(indj+1,indi+k),k = 1,6)
              write(ioout,'(i6,'' y '',4x,6f10.6)') j,(vectors(indj+2,indi+k),k = 1,6)
              write(ioout,'(i6,'' z '',4x,6f10.6)') j,(vectors(indj+3,indi+k),k = 1,6)
              indj = indj + 3
            enddo
            indi = indi + 6
            write(ioout,'(/)')
          enddo
        endif
        if (iresid.gt.0) then
          write(ioout,'('' Frequency   '',6f10.4)') (freq(indi+j,1),j = 1,iresid)
          write(ioout,'(/)',advance='no')
          indj = 0
          do j = 1,ncfoc
            write(ioout,'(i6,'' x '',4x,6f10.6)') j,(vectors(indj+1,indi+k),k = 1,iresid)
            write(ioout,'(i6,'' y '',4x,6f10.6)') j,(vectors(indj+2,indi+k),k = 1,iresid)
            write(ioout,'(i6,'' z '',4x,6f10.6)') j,(vectors(indj+3,indi+k),k = 1,iresid)
            indj = indj + 3
          enddo
          write(ioout,'(/)')
        endif
      endif
      write(ioout,'(/)')
    elseif (linten.and.ioproc) then
      write(ioout,'(/,''  Frequencies (cm-1) and IR Intensities : '',/)')
      igroup = mcv/6
      iresid = mcv - igroup*6
      indi = 0
      if (igroup.gt.0) then
        do i = 1,igroup
          write(ioout,'(''  Frequency   '',6f10.4)') (freq(indi+j,1),j = 1,6)
          write(ioout,'(''  IR Intensity'',6f10.4)') (IRintensity(indi+j,1),j = 1,6)
          write(ioout,'(''     in X     '',6f10.4)') (irx(indi+j),j = 1,6)
          write(ioout,'(''     in Y     '',6f10.4)') (iry(indi+j),j = 1,6)
          write(ioout,'(''     in Z     '',6f10.4,/)') (irz(indi+j),j = 1,6)
          indi = indi + 6
        enddo
      endif
      if (iresid.gt.0) then
        write(ioout,'(''  Frequency   '',6f10.4)') (freq(indi+j,1),j = 1,iresid)
        write(ioout,'(''  IR Intensity'',6f10.4)') (IRintensity(indi+j,1),j = 1,iresid)
        write(ioout,'(''     in X     '',6f10.4)') (irx(indi+j),j = 1,iresid)
        write(ioout,'(''     in Y     '',6f10.4)') (iry(indi+j),j = 1,iresid)
        write(ioout,'(''     in Z     '',6f10.4)') (irz(indi+j),j = 1,iresid)
        write(ioout,'(/)',advance='no')
      endif
      write(ioout,'(/)')
    elseif (lfreqout.and.lprinloc) then
      write(ioout,'(/,''  Frequencies (cm-1) : '',/)')
      write(ioout,'(9(f8.2))') (freq(j,1),j = 1,mcv)
      write(ioout,'(/)')
    endif
  else
    if (lfreqout.and.lprinloc) then
      write(ioout,'(/,''  Frequencies (cm-1) :'',/)')
      write(ioout,'(9f8.2)')(freq(i,1),i = 1,mcv)
      write(ioout,'(/)')
    endif
  endif
!
!  Output eigenvectors if requested
!
  if (lfreqout.and.lprinloc) then
    write(ioout,'(''--------------------------------------------------------------------------------'',/)')
  endif
!
!  Lower symmetry to remove imaginary modes if selected
!
  call lower(mcv,mcvrptr,freq(1,1),ncfoc,iocptr,ncfoc,iocptr,maxd2,vectors)
!
  mcvmin = 1
!
!  Check for imaginary modes and exclude them
!
  nimag = 0
  do i = 1,mcv
    if (freq(i,1).lt.-0.5_dp) nimag = nimag + 1
  enddo
  mcvmin = mcvmin + max(0,nimag-3)
  mcvmax = mcv
  if (minmode.ne.1) mcvmin = minmode
  if (maxmode.ne.0) mcvmax = maxmode
!*************************************
!  Output phonon related properties  *
!*************************************
  if (lprinloc) then
    tem = temperature - temperaturestep
    do nt = 0,ntemperaturestep
      tem = tem + temperaturestep
!
!  Zero thermodynamic properties
!
      zpe = 0.0_dp
      entropy = 0.0_dp
      fhenergy = 0.0_dp
      rtlnz = 0.0_dp
      cv = 0.0_dp
      fe_equipartition = 0.0_dp
      s_equipartition = 0.0_dp
!***************************************
!  Evaluate phonon related properties  *
!***************************************
      if (tem.gt.1.0d-3) then
!
!  Scale frequencies to hw/kT
!
        rkt = boltz*tem
        cmfact = planck*speedl/rkt
        do i = 1,mcv
          rtmp2(i) = cmfact*freq(i,1)
        enddo
!
!  Store exp(x) in w1, exp(x)-1 in w2 and 1/(exp(x)-1) in w3
!
        do i = 1,mcv
          if (rtmp2(i).lt.12.0_dp) then
            w1(i) = exp(rtmp2(i))
            w2(i) = w1(i) - 1.0_dp
            if (abs(w2(i)).gt.0.0_dp) w3(i) = 1.0_dp/w2(i)
          else
            w3(i) = exp(-rtmp2(i))
          endif
        enddo
!
!  Zero point energy
!
        factor = 0.5_dp*rkt/evtoj
        trmzp = 0.0_dp
        if (.not.lnozero) then
          do i = mcvmin,mcvmax
            if (rtmp2(i).gt.cmfact) trmzp = trmzp + rtmp2(i)
          enddo
          trmzp = factor*trmzp
          zpe = zpe + trmzp
        endif
!
!  Entropy and free energy
!
        trmfe = 0.0_dp
        trmen = 0.0_dp
        freqmin = cmfact
        do i = mcvmin,mcvmax
          trm1 = rtmp2(i)
          if (trm1.gt.freqmin) then
            trm1 = 1.0_dp - exp(-trm1)
            trmfe = trmfe + log(trm1)
            trmen = trmen + rtmp2(i)*w3(i)
          endif
        enddo
!
!  Equipartition free energy and entropy
!
        trmfe_eq = 0.0_dp
        trms_eq = 0.0_dp
        rmode = 0.0_dp
        do i = 1,mcv
          trm1 = rtmp2(i)
          if (trm1.gt.freqmin) then
            trmfe_eq = trmfe_eq + log(trm1)
            trms_eq  = trms_eq + log(trm1) - 1.0_dp
            rmode = rmode + 1.0_dp
          endif
        enddo
!
        factor = 2.0_dp*factor
        trm = factor*trmfe
        fhenergy = fhenergy + trm + trmzp
        rtlnz = rtlnz - trm
        fe_equipartition = fe_equipartition + factor*trmfe_eq
        s_equipartition = s_equipartition + factor*trms_eq
        factor = factor/tem
        entropy = entropy + factor*trmen
!
!  Heat capacity - constant volume
!
        trmcv = 0.0_dp
        do i = mcvmin,mcvmax
          if (rtmp2(i).gt.freqmin) then
            if (rtmp2(i).lt.12.0_dp) then
              trmcv = trmcv + rtmp2(i)*rtmp2(i)*w1(i)*w3(i)*w3(i)
            else
              trmcv = trmcv + rtmp2(i)*rtmp2(i)*w3(i)
            endif
          endif
        enddo
        cv = cv + factor*trmcv
      else
!
!  Zero point energy
!
        factor = 0.5_dp*planck*speedl/evtoj
        if (.not.lnozero) then
          trmzp = 0.0_dp
          do i = mcvmin,mcvmax
            if (freq(i,1).gt.1.0_dp) trmzp = trmzp + freq(i,1)
          enddo
          trmzp = factor*trmzp
          zpe = zpe + trmzp
        endif
      endif
      write(ioout,'(''  Vibrational properties (for region 1): '',''Temperature  =  '',f10.3,'' K'')')tem
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  Zero point energy             =  '',f15.6,'' eV'')')zpe
      if (tem.gt.1.0d-03) then
        trmen = fhenergy - zpe
        entropy = entropy - trmen/tem
        ent2 = entropy*evtoj*avogadro
        cv2 = cv*evtoj*avogadro
        write(ioout,'(''  Entropy                       =  '',f15.6,'' eV/K'')') entropy
        write(ioout,'(''                                =  '',f15.6,'' J/(mol.K)'')') ent2
        write(ioout,'(''  Helmholtz free-energy         =  '',f15.6,'' eV'')') fhenergy + fc
        write(ioout,'(''                                =  '',f15.6,'' kJmol-1'')') (fhenergy+fc)*evtoj*avogadro*0.001_dp
        write(ioout,'(''  Heat capacity - const volume  =  '',f15.6,'' eV/K'')') cv
        write(ioout,'(''                                =  '',f15.6,'' J/(mol.K)'')') cv2
      endif
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    enddo
  endif
!
!  If output if required then call outfreq
!
  if (lprinloc) call outfreq(leigloc,1_i4,freq,IRintensity,ncfoc)
!***************
!  Exit point  *
!***************
999 continue
!
!  Free local memory
!
  deallocate(rmassdef,stat=status)
  if (status/=0) call deallocate_error('deffreq','rmassdef')
  deallocate(massdef,stat=status)
  if (status/=0) call deallocate_error('deffreq','massdef')
  deallocate(w3,stat=status)
  if (status/=0) call deallocate_error('deffreq','w3')
  deallocate(w2,stat=status)
  if (status/=0) call deallocate_error('deffreq','w2')
  deallocate(w1,stat=status)
  if (status/=0) call deallocate_error('deffreq','w1')
  deallocate(rtmp2,stat=status)
  if (status/=0) call deallocate_error('deffreq','rtmp2')
  deallocate(irz,stat=status)
  if (status/=0) call deallocate_error('deffreq','irz')
  deallocate(iry,stat=status)
  if (status/=0) call deallocate_error('deffreq','iry')
  deallocate(irx,stat=status)
  if (status/=0) call deallocate_error('deffreq','irx')
  deallocate(mcvrptr,stat=status)
  if (status/=0) call deallocate_error('deffreq','mcvrptr')
  deallocate(mcvptr,stat=status)
  if (status/=0) call deallocate_error('deffreq','mcvptr')
  deallocate(iocptr,stat=status)
  if (status/=0) call deallocate_error('deffreq','iocptr')
  deallocate(ibocptr,stat=status)
  if (status/=0) call deallocate_error('deffreq','ibocptr')
  deallocate(ndbsptr,stat=status)
  if (status/=0) call deallocate_error('deffreq','ndbsptr')
!
!  Timing
!
  t2t = g_cpu_time()
  tphon = t2t - t1t + tphon
#ifdef TRACE
  call trace_out('deffreq')
#endif
!
  return
  end
