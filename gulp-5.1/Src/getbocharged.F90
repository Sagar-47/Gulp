  subroutine getBOcharged(lgrad1,lgrad2)
!
!  Calculates the charges according to the bond order formalism.
!
!  Distributed memory parallel version.
!
!  On entry :
!
!  lgrad1       = if .true. calculate first derivatives
!  lgrad2       = if .true. calculate second derivatives
!
!  NB: The second charge derivatives with respect to strain exclude the
!      delta terms since these are handled in strfin anyway.
!
!   2/17 Created from getbocharge
!   7/17 Output format changed
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
  use datatypes
  use bondorderdata
  use control,        only : keyword
  use current
  use derivatives,    only : dqdxyz, d2qdxyz2, d2qdxyzs, dqds, d2qds2, nqatoms, nqatomptr, nqatomcell
  use derivatives,    only : maxqatoms, maxqatoms2, maxd2q, maxd2qu, qatomxyz
  use iochannels
  use neighbours
  use parallel,       only : ioproc, natomsonnode, node2atom
  use spatial
  use symmetry,       only : lstr
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  logical,     intent(in)                          :: lgrad1
  logical,     intent(in)                          :: lgrad2
!
!  Local variables
!
  integer(i4)                                      :: i
  integer(i4)                                      :: ic
  integer(i4)                                      :: ii
  integer(i4)                                      :: ijx
  integer(i4)                                      :: ijy
  integer(i4)                                      :: ijz
  integer(i4)                                      :: ijxx
  integer(i4)                                      :: ijxy
  integer(i4)                                      :: ijxz
  integer(i4)                                      :: ijyy
  integer(i4)                                      :: ijyz
  integer(i4)                                      :: ijzz
  integer(i4)                                      :: iloc
  integer(i4)                                      :: imx
  integer(i4)                                      :: imy
  integer(i4)                                      :: imz
  integer(i4)                                      :: ind
  integer(i4)                                      :: ind2
  integer(i4)                                      :: indn
  integer(i4), dimension(:,:,:), allocatable, save :: ineigh
  integer(i4)                                      :: itmp
  integer(i4)                                      :: ix
  integer(i4)                                      :: ixyz
  integer(i4)                                      :: iy
  integer(i4)                                      :: iz
  integer(i4)                                      :: j
  integer(i4)                                      :: jc
  integer(i4)                                      :: jix
  integer(i4)                                      :: jiy
  integer(i4)                                      :: jiz
  integer(i4)                                      :: jj
  integer(i4)                                      :: jx
  integer(i4)                                      :: jy
  integer(i4)                                      :: jz
  integer(i4)                                      :: k
  integer(i4)                                      :: kk
  integer(i4)                                      :: kl
  integer(i4)                                      :: ks
  integer(i4)                                      :: kt
  integer(i4)                                      :: m
  integer(i4)                                      :: maxqatomsloc
  integer(i4)                                      :: maxqatomsloc2
  integer(i4)                                      :: maxxy
  integer(i4)                                      :: maxx
  integer(i4)                                      :: n1i
  integer(i4)                                      :: n1j
  integer(i4)                                      :: nati
  integer(i4)                                      :: natj
  integer(i4)                                      :: nqatoms2
  integer(i4)                                      :: nboij
  integer(i4)                                      :: ni
  integer(i4)                                      :: nj
  integer(i4)                                      :: nmin
  integer(i4)                                      :: nn
  integer(i4)                                      :: nn2
  integer(i4)                                      :: nptr
  integer(i4), dimension(:,:),   allocatable, save :: neighno
  integer(i4), dimension(:),     allocatable, save :: nneigh
  integer(i4)                                      :: nsplower(3)
  integer(i4)                                      :: nspupper(3)
  integer(i4)                                      :: ntypi
  integer(i4)                                      :: ntypj
  integer(i4)                                      :: status
  logical                                          :: lfound
  logical                                          :: lfound1
  logical                                          :: lmaxneighok
  logical                                          :: lok
  real(dp)                                         :: bR22
  real(dp)                                         :: g_cpu_time
  real(dp)                                         :: d1trm
  real(dp)                                         :: d2trm
  real(dp)                                         :: dfdr
  real(dp)                                         :: d2fdr2
  real(dp)                                         :: d3fdr3
  real(dp)                                         :: f
  real(dp),    dimension(:),     allocatable, save :: rBOcutmax
  real(dp)                                         :: rij
  real(dp)                                         :: r2
  real(dp)                                         :: rpdji(6)
  real(dp)                                         :: rrij
  real(dp)                                         :: rtmp
  real(dp)                                         :: sBOq0
  real(dp)                                         :: t1
  real(dp)                                         :: t2
  real(dp),    dimension(:,:),   allocatable, save :: rneigh
  real(dp),    dimension(:,:),   allocatable, save :: xneigh
  real(dp),    dimension(:,:),   allocatable, save :: yneigh
  real(dp),    dimension(:,:),   allocatable, save :: zneigh
  real(dp)                                         :: xdiff
  real(dp)                                         :: ydiff
  real(dp)                                         :: zdiff
  real(dp)                                         :: xi
  real(dp)                                         :: yi     
  real(dp)                                         :: zi
  real(dp)                                         :: xji
  real(dp)                                         :: yji
  real(dp)                                         :: zji
  real(dp)                                         :: xji0
  real(dp)                                         :: yji0
  real(dp)                                         :: zji0
#ifdef TRACE
  call trace_in('getbocharged')
#endif
!
  t1 = g_cpu_time()
!
!  Allocate memory that does not depend on maxneigh
!
  allocate(nneigh(numat),stat=status)
  if (status/=0) call outofmemory('getBOcharged','nneigh')
  allocate(rBOcutmax(numat),stat=status)
  if (status/=0) call outofmemory('getBOcharged','rBOcutmax')
!
!  Reinitialisation point should maxneigh be increased
!
100 continue
  lmaxneighok = .true.
  if (allocated(rneigh)) then
    deallocate(zneigh,stat=status)
    if (status/=0) call deallocate_error('getBOcharged','zneigh')
    deallocate(yneigh,stat=status)
    if (status/=0) call deallocate_error('getBOcharged','yneigh')
    deallocate(xneigh,stat=status)
    if (status/=0) call deallocate_error('getBOcharged','xneigh')
    deallocate(rneigh,stat=status)
    if (status/=0) call deallocate_error('getBOcharged','rneigh')
    deallocate(ineigh,stat=status)
    if (status/=0) call deallocate_error('getBOcharged','ineigh')
    deallocate(neighno,stat=status)
    if (status/=0) call deallocate_error('getBOcharged','neighno')
  endif
!
!  Initialise charges
!
  do i = 1,numat
    qf(i) = 0.0_dp
  enddo
!
!  Allocate local memory
!
  allocate(neighno(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('getBOcharged','neighno')
  allocate(ineigh(3_i4,maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('getBOcharged','ineigh')
  allocate(rneigh(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('getBOcharged','rneigh')
  allocate(xneigh(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('getBOcharged','xneigh')
  allocate(yneigh(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('getBOcharged','yneigh')
  allocate(zneigh(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('getBOcharged','zneigh')
!*************************************
!  Find cut-off radii for all atoms  *
!*************************************
  do i = 1,numat
    nati = nat(i)
    ntypi = nftype(i)
    rBOcutmax(i) = 0.0_dp
!
!  Check twobody potentials
!
    do j = 1,nboQ
      if (nati.eq.nBOspecQ1(j).and.(ntypi.eq.nBOtypQ1(j).or.nBOtypQ1(j).eq.0)) then
        rBOcutmax(i) = max(rBOcutmax(i),rBOmaxQ(j))
      endif
      if (nati.eq.nBOspecQ2(j).and.(ntypi.eq.nBOtypQ2(j).or.nBOtypQ2(j).eq.0)) then
        rBOcutmax(i) = max(rBOcutmax(i),rBOmaxQ(j))
      endif
    enddo
  enddo
!****************************************
!  Calculate neighbour lists for atoms  *
!****************************************
  if (lspatialok) then
    maxxy = nspcell(1)*nspcell(2)
    maxx  = nspcell(1)
!  
!  Loop over all spatial cells
!     
    do ixyz = 1,ncellpernode
      ind = ncellnodeptr(ixyz)
      ind2 = ind - 1
      iz = ind2/maxxy
      ind2 = ind2 - maxxy*iz
      iy = ind2/maxx
      ix = ind2 - maxx*iy + 1
      iy = iy + 1
      iz = iz + 1
!
!  Set cell search bounds
!
      nspupper(1) = min(ix+ncellsearch(1),nspcell(1))
      nspupper(2) = min(iy+ncellsearch(2),nspcell(2))
      nspupper(3) = min(iz+ncellsearch(3),nspcell(3))
      nsplower(1) = max(ix-ncellsearch(1),1)
      nsplower(2) = max(iy-ncellsearch(2),1)
      nsplower(3) = max(iz-ncellsearch(3),1)
!     
!  Get number of atoms in this cell
!       
      ni = nspcellat(ind)
      n1i = nspcellat1ptr(ind)
!  
!  Loop over atoms in the cell finding neighbours
!     
      do ii = 1,ni
        i = nspcellatptr(n1i+ii)
        ic = nspcellatptrcell(n1i+ii)
        nneigh(i) = 0
        nati = nat(i)
        ntypi = nftype(i)
!       
!  Set coordinates
!
        xi = xinbox(i) + xvec2cell(ic)
        yi = yinbox(i) + yvec2cell(ic)
        zi = zinbox(i) + zvec2cell(ic)
!                             
!  Compute square of cut-off for distance checking
!  
        bR22 = rBOcutmax(i)**2
!
!  Loop over neighbouring cells
!
        do imz = nsplower(3),nspupper(3)
          do imy = nsplower(2),nspupper(2)
            do imx = nsplower(1),nspupper(1)
              indn = (imz-1)*maxxy + (imy-1)*maxx + imx
!                         
!  Loop over atoms within neighbouring cells  
!                         
              nj = nspcellat(indn)
              n1j = nspcellat1ptr(indn)
              do jj = 1,nj
                j = nspcellatptr(n1j+jj)
                jc = nspcellatptrcell(n1j+jj)
!                     
!  Exclude self term    
!                         
                if (i.ne.j.or.ind.ne.indn) then
!                             
!  Set centre cell coordinate differences
!  
                  xji = xvec2cell(jc) + xinbox(j) - xi
                  yji = yvec2cell(jc) + yinbox(j) - yi
                  zji = zvec2cell(jc) + zinbox(j) - zi
!  
                  r2 = xji*xji + yji*yji + zji*zji
                  if (r2 .lt. bR22) then
!
!  Check whether j is within two body cut-off
!
                    natj = nat(j)
                    ntypj = nftype(j)
                    m = 0
                    lok = .false.
                    do while (m.lt.nboQ.and..not.lok)
                      m = m + 1
                      if (nati.eq.nBOspecQ1(m).and.(ntypi.eq.nBOtypQ1(m).or.nBOtypQ1(m).eq.0).and. &
                          natj.eq.nBOspecQ2(m).and.(ntypj.eq.nBOtypQ2(m).or.nBOtypQ2(m).eq.0)) then
                        lok = (r2.lt.rBOmaxQ(m)**2)
                      elseif (nati.eq.nBOspecQ2(m).and.(ntypi.eq.nBOtypQ2(m).or.nBOtypQ2(m).eq.0).and. &
                          natj.eq.nBOspecQ1(m).and.(ntypj.eq.nBOtypQ1(m).or.nBOtypQ1(m).eq.0)) then
                        lok = (r2.lt.rBOmaxQ(m)**2)
                      endif
                    enddo
                    if (lok) then
                      if (nneigh(i).ge.maxneigh.or..not.lmaxneighok) then
                        lmaxneighok = .false.
                        nneigh(i) = nneigh(i) + 1
                      else        
                        rij = sqrt(r2)
                        nneigh(i) = nneigh(i) + 1
                        neighno(nneigh(i),i) = j
                        ineigh(1,nneigh(i),i) = ivec2cell(1,jc) - ivec2cell(1,ic)
                        ineigh(2,nneigh(i),i) = ivec2cell(2,jc) - ivec2cell(2,ic)
                        ineigh(3,nneigh(i),i) = ivec2cell(3,jc) - ivec2cell(3,ic)
                        rneigh(nneigh(i),i) = rij
                        xneigh(nneigh(i),i) = xji
                        yneigh(nneigh(i),i) = yji
                        zneigh(nneigh(i),i) = zji
                      endif
                    endif
                  endif
                endif
!                             
              enddo
            enddo
          enddo
        enddo             
!                                   
      enddo                 
!                                 
    enddo                       
  else
    do i = 1,numat
      nneigh(i) = 0
      nati = nat(i)
      ntypi = nftype(i)
!     
!  Compute square of cut-off for distance checking   
!     
      bR22 = rBOcutmax(i)**2
!
!  Loop over atoms
!
      do j = 1,numat
        natj = nat(j)
        ntypj = nftype(j)
!
!  Set centre cell coordinate differences
!
        xji0 = xclat(j) - xclat(i)
        yji0 = yclat(j) - yclat(i)
        zji0 = zclat(j) - zclat(i)
!
!  Loop over unit cells
!
        do ii = 1,iimax
!
!  Exclude self term
!
          if (i.ne.j.or.ii.ne.iimid) then
            xji = xji0 + xvec1cell(ii)
            yji = yji0 + yvec1cell(ii)
            zji = zji0 + zvec1cell(ii)
            r2 = xji*xji + yji*yji + zji*zji
            if (r2 .lt. bR22) then
              m = 0
              lok = .false.   
              do while (m.lt.nboQ.and..not.lok)
                m = m + 1   
                if (nati.eq.nBOspecQ1(m).and.(ntypi.eq.nBOtypQ1(m).or.nBOtypQ1(m).eq.0).and. &
                    natj.eq.nBOspecQ2(m).and.(ntypj.eq.nBOtypQ2(m).or.nBOtypQ2(m).eq.0)) then
                  lok = (r2.lt.rBOmaxQ(m)**2)
                elseif (nati.eq.nBOspecQ2(m).and.(ntypi.eq.nBOtypQ2(m).or.nBOtypQ2(m).eq.0).and. &
                    natj.eq.nBOspecQ1(m).and.(ntypj.eq.nBOtypQ1(m).or.nBOtypQ1(m).eq.0)) then
                  lok = (r2.lt.rBOmaxQ(m)**2)
                endif
              enddo
              if (lok) then
                if (nneigh(i).ge.maxneigh.or..not.lmaxneighok) then
                  lmaxneighok = .false.
                  nneigh(i) = nneigh(i) + 1
                else
                  rij = sqrt(r2)
                  nneigh(i) = nneigh(i) + 1
                  neighno(nneigh(i),i) = j
                  ineigh(1,nneigh(i),i) = ivec1cell(1,ii)
                  ineigh(2,nneigh(i),i) = ivec1cell(2,ii)
                  ineigh(3,nneigh(i),i) = ivec1cell(3,ii)
                  rneigh(nneigh(i),i) = rij
                  xneigh(nneigh(i),i) = xji
                  yneigh(nneigh(i),i) = yji
                  zneigh(nneigh(i),i) = zji
                endif
              endif
            endif
          endif
        enddo
      enddo
    enddo
  endif
!**********************************************************************
!  If maxneigh has been exceeded, increase value and return to start  *
!**********************************************************************
  if (.not.lmaxneighok) then
    do i = 1,numat
      if (nneigh(i).gt.maxneigh) maxneigh = nneigh(i)
    enddo
    if (ioproc.and.index(keyword,'verb').ne.0) then
      write(ioout,'(/,''  Increasing maxneigh to '',i6)') maxneigh
    endif
    goto 100
  endif
!*******************************
!  Sort neighbours into order  *
!******************************* 
  if (lspatialok) then
    do i = 1,numat
!               
!  Build pointer
!               
      do nn = 1,nneigh(i)
        nmin = numat + 1 
        do nn2 = nn,nneigh(i) 
          if (neighno(nn2,i).lt.nmin) then
            nmin = neighno(nn2,i)
            nptr = nn2  
          endif
        enddo       
!         
!  Sort quantities
!
        if (nptr.ne.nn) then
          itmp = neighno(nptr,i)
          neighno(nptr,i) = neighno(nn,i)
          neighno(nn,i)  = itmp
          itmp = ineigh(1,nptr,i)
          ineigh(1,nptr,i) = ineigh(1,nn,i)
          ineigh(1,nn,i)  = itmp
          itmp = ineigh(2,nptr,i)
          ineigh(2,nptr,i) = ineigh(2,nn,i)
          ineigh(2,nn,i)  = itmp
          itmp = ineigh(3,nptr,i)
          ineigh(3,nptr,i) = ineigh(3,nn,i)
          ineigh(3,nn,i)  = itmp
          rtmp = rneigh(nptr,i)
          rneigh(nptr,i) = rneigh(nn,i)
          rneigh(nn,i)  = rtmp
          rtmp = xneigh(nptr,i)
          xneigh(nptr,i) = xneigh(nn,i)
          xneigh(nn,i)  = rtmp
          rtmp = yneigh(nptr,i)
          yneigh(nptr,i) = yneigh(nn,i)
          yneigh(nn,i)  = rtmp
          rtmp = zneigh(nptr,i)
          zneigh(nptr,i) = zneigh(nn,i)
          zneigh(nn,i)  = rtmp
        endif  
      enddo         
    enddo
  endif
!*********************************************
!  Calculate numbers of neighbours for each  *
!*********************************************
  if (ioproc.and.index(keyword,'debu').ne.0) then
    write(ioout,'(/,''  Number of neighbours for atoms :'',/)')
    write(ioout,'(''  i'',2x,''Nat'',3x,''N'')')
    do i = 1,numat
      write(ioout,'(i4,1x,i3,1x,i6)') i,nat(i),nneigh(i)
    enddo
    write(ioout,'(/,''  Neighbours of atoms :'',/)')
    do i = 1,numat
      write(ioout,'(i4,8(1x,i4))') i,(neighno(nn,i),nn=1,nneigh(i))
    enddo
  endif
!*******************************************************
!  Copy neighbour information into charge atom arrays  *
!*******************************************************
  if (lgrad1) then
    maxqatomsloc = 0
    do i = 1,numat
      maxqatomsloc = max(maxqatomsloc,nneigh(i))
    enddo
    if (maxqatomsloc.gt.maxqatoms) then
      maxqatoms = maxqatomsloc
      call changemaxqatoms
    endif
    do iloc = 1,natomsonnode
      i = node2atom(iloc)
      nqatoms(iloc) = nneigh(i)
      do nn = 1,nneigh(i)
        nqatomptr(nn,iloc) = neighno(nn,i)
        nqatomcell(1:3,nn,iloc) = ineigh(1:3,nn,i)
        qatomxyz(1,nn,iloc) = xneigh(nn,i)
        qatomxyz(2,nn,iloc) = yneigh(nn,i)
        qatomxyz(3,nn,iloc) = zneigh(nn,i)
      enddo
    enddo
  endif
!*****************************
!  Set up derivative arrays  *
!*****************************
  if (lgrad1) then
!
!  Check dimensions
!
    if (natomsonnode.gt.maxd2qu) then                   
      maxd2qu = natomsonnode
      call changemaxd2q                          
    endif
    if (3*numat.gt.maxd2q) then                  
      maxd2q = 3*numat                           
      call changemaxd2q                          
    endif
    if (lgrad2) then
      maxqatomsloc2 = (3*maxqatomsloc + 3)*(3*maxqatomsloc + 6)/2
      if (maxqatomsloc2.gt.maxqatoms2) then
        maxqatoms2 = maxqatomsloc2
        call changemaxqatoms2
      endif
    endif
!   
!  Zero arrays
!  
    do i = 1,natomsonnode
      do j = 1,3*numat
        dqdxyz(j,i) = 0.0_dp
      enddo
    enddo
    if (ndim.gt.0.and.lstr) then
      do i = 1,natomsonnode
        do j = 1,nstrains
          dqds(j,i) = 0.0_dp 
        enddo 
      enddo
    endif
    if (lgrad2) then
      do i = 1,natomsonnode
        nqatoms2 = (3*nqatoms(i) + 3)*(3*nqatoms(i) + 6)/2
        do j = 1,nqatoms2
          d2qdxyz2(j,i) = 0.0_dp
        enddo
      enddo
      if (ndim.gt.0.and.lstr) then
        do i = 1,natomsonnode
          do j = 1,3*nqatoms(i)
            do k = 1,nstrains
              d2qdxyzs(k,j,i) = 0.0_dp
            enddo
          enddo
          do j = 1,nstrains*(nstrains + 1)/2
            d2qds2(j,i) = 0.0_dp 
          enddo 
        enddo
      endif
    endif
  endif
!*****************************
!  Loop over pairs of atoms  *
!*****************************
  do iloc = 1,natomsonnode
    i = node2atom(iloc)
!
!  Set variables relating to i
!
    nati = nat(i)
    ntypi = nftype(i)
    if (lgrad1) then
      ix = 3*(i-1) + 1
      iy = ix + 1
      iz = ix + 2
    endif
!
!  Loop over neighbours of i 
!
    ni = 1
    do while (ni.le.nneigh(i))
!
      j = neighno(ni,i)
!
!  Set variables relating to j
!
      natj = nat(j)
      ntypj = nftype(j)
      if (lgrad1) then
        jx = 3*(j-1) + 1
        jy = jx + 1
        jz = jx + 2
      endif
!
!  Set up i-j quantities
!
      rij = rneigh(ni,i)
      xji = xneigh(ni,i)
      yji = yneigh(ni,i)
      zji = zneigh(ni,i)
!
!  Find i in neighbour list for j
!
      nj = 1
      lfound = .false.
      do while (nj.lt.nneigh(j).and..not.lfound)
        if (neighno(nj,j).eq.i) then
          xdiff = xneigh(nj,j) + xji
          ydiff = yneigh(nj,j) + yji
          zdiff = zneigh(nj,j) + zji
          lfound = ((abs(xdiff)+abs(ydiff)+abs(zdiff)).lt.1.0d-6)
        endif
        if (.not.lfound) nj = nj + 1
      enddo
!
!  Find two-body bond order potential between i and j
!
      lfound = .false.
      lfound1 = .false.
      nboij = 0
      do while (.not.lfound.and.nboij.lt.nboQ) 
        nboij = nboij + 1
        if (nBOspecQ1(nboij).eq.nati.and.nBOspecQ2(nboij).eq.natj) then
          if ((nBOtypQ1(nboij).eq.ntypi.or.nBOtypQ1(nboij).eq.0).and. &
              (nBOtypQ2(nboij).eq.ntypj.or.nBOtypQ2(nboij).eq.0)) then
            lfound = .true.
            lfound1 = .true.
          endif
        elseif (nBOspecQ1(nboij).eq.natj.and.nBOspecQ2(nboij).eq.nati) then
          if ((nBOtypQ1(nboij).eq.ntypj.or.nBOtypQ1(nboij).eq.0).and. &
              (nBOtypQ2(nboij).eq.ntypi.or.nBOtypQ2(nboij).eq.0)) then
            lfound = .true.
          endif
        endif
      enddo
      if (lfound) then
!***********************************************
!  Valid two-body bond order charge potential  *
!***********************************************
!
!  Calculate fij
!
        if (nBOtaperQ(nboij).eq.2) then
          call staper(rij,rBOminQ(nboij),rBOmaxQ(nboij),f,dfdr,d2fdr2,d3fdr3,lgrad1,lgrad2,.false.)
        else
          call ctaper(rij,rBOminQ(nboij),rBOmaxQ(nboij),f,dfdr,d2fdr2,d3fdr3,lgrad1,lgrad2,.false.)
        endif
!
!  Add contribution to charge
!
        if (lfound1) then
          qf(i) = qf(i) + BOq0(nboij)*f
        else
          qf(i) = qf(i) - BOq0(nboij)*f
        endif
        if (lgrad1) then
!
!  First derivatives
!
          rrij = 1.0_dp/rij
          if (lfound1) then
            sBOq0 = BOq0(nboij)
          else
            sBOq0 = - BOq0(nboij)
          endif
          d1trm = sBOq0*rrij*dfdr
!
          dqdxyz(ix,iloc) = dqdxyz(ix,iloc) - d1trm*xji
          dqdxyz(iy,iloc) = dqdxyz(iy,iloc) - d1trm*yji
          dqdxyz(iz,iloc) = dqdxyz(iz,iloc) - d1trm*zji
          dqdxyz(jx,iloc) = dqdxyz(jx,iloc) + d1trm*xji
          dqdxyz(jy,iloc) = dqdxyz(jy,iloc) + d1trm*yji
          dqdxyz(jz,iloc) = dqdxyz(jz,iloc) + d1trm*zji
          if (ndim.gt.0.and.lstr) then
            rpdji(1) = xji*xji
            rpdji(2) = yji*yji
            rpdji(3) = zji*zji
            rpdji(4) = yji*zji
            rpdji(5) = xji*zji
            rpdji(6) = xji*yji
            do kl = 1,nstrains
              ks = nstrptr(kl)
              dqds(kl,iloc) = dqds(kl,iloc) + d1trm*rpdji(ks)
            enddo
          endif
          if (lgrad2) then
!***********************
!  Second derivatives  *
!***********************
            jix = 3*(ni - 1) + 1
            jiy = jix + 1
            jiz = jiy + 1
!
            ijx = 3*(nj - 1) + 1
            ijy = ijx + 1
            ijz = ijy + 1
!
            if (i.ne.j) then
              d2trm = sBOq0*rrij*rrij*(d2fdr2 - rrij*dfdr)
!
              ijxx = jix*(jix + 1)/2
              ijyy = jiy*(jiy + 1)/2
              ijzz = jiz*(jiz + 1)/2
              ijxy = ijyy - 1
              ijxz = ijzz - 2
              ijyz = ijzz - 1
!
              d2qdxyz2(ijxx,iloc) = d2qdxyz2(ijxx,iloc) - d2trm*xji*xji
              d2qdxyz2(ijxy,iloc) = d2qdxyz2(ijxy,iloc) - d2trm*xji*yji
              d2qdxyz2(ijxz,iloc) = d2qdxyz2(ijxz,iloc) - d2trm*xji*zji
              d2qdxyz2(ijyy,iloc) = d2qdxyz2(ijyy,iloc) - d2trm*yji*yji
              d2qdxyz2(ijyz,iloc) = d2qdxyz2(ijyz,iloc) - d2trm*yji*zji
              d2qdxyz2(ijzz,iloc) = d2qdxyz2(ijzz,iloc) - d2trm*zji*zji
              d2qdxyz2(ijxx,iloc) = d2qdxyz2(ijxx,iloc) - d1trm
              d2qdxyz2(ijyy,iloc) = d2qdxyz2(ijyy,iloc) - d1trm
              d2qdxyz2(ijzz,iloc) = d2qdxyz2(ijzz,iloc) - d1trm
            endif
!
!  Strain derivatives
!
            if (ndim.gt.0.and.lstr) then
              ind = 0
              do kk = 1,nstrains
                ks = nstrptr(kk)
                do kl = 1,kk
                  kt = nstrptr(kl)
                  ind = ind + 1
                  d2qds2(ind,iloc) = d2qds2(ind,iloc) + d2trm*rpdji(ks)*rpdji(kt)
                enddo
!
                d2qdxyzs(kk,jix,iloc) = d2qdxyzs(kk,jix,iloc) + d2trm*rpdji(ks)*xji
                d2qdxyzs(kk,jiy,iloc) = d2qdxyzs(kk,jiy,iloc) + d2trm*rpdji(ks)*yji
                d2qdxyzs(kk,jiz,iloc) = d2qdxyzs(kk,jiz,iloc) + d2trm*rpdji(ks)*zji
              enddo
              if (ndim.eq.3) then
                d2qdxyzs(1,jix,iloc) = d2qdxyzs(1,jix,iloc) + 2.0_dp*d1trm*xji
                d2qdxyzs(2,jiy,iloc) = d2qdxyzs(2,jiy,iloc) + 2.0_dp*d1trm*yji
                d2qdxyzs(3,jiz,iloc) = d2qdxyzs(3,jiz,iloc) + 2.0_dp*d1trm*zji
                d2qdxyzs(5,jix,iloc) = d2qdxyzs(5,jix,iloc) + d1trm*zji
                d2qdxyzs(6,jix,iloc) = d2qdxyzs(6,jix,iloc) + d1trm*yji
                d2qdxyzs(4,jiy,iloc) = d2qdxyzs(4,jiy,iloc) + d1trm*zji
                d2qdxyzs(6,jiy,iloc) = d2qdxyzs(6,jiy,iloc) + d1trm*xji
                d2qdxyzs(4,jiz,iloc) = d2qdxyzs(4,jiz,iloc) + d1trm*yji
                d2qdxyzs(5,jiz,iloc) = d2qdxyzs(5,jiz,iloc) + d1trm*xji
              elseif (ndim.eq.2) then
                d2qdxyzs(1,jix,iloc) = d2qdxyzs(1,jix,iloc) + 2.0_dp*d1trm*xji
                d2qdxyzs(2,jiy,iloc) = d2qdxyzs(2,jiy,iloc) + 2.0_dp*d1trm*yji
                d2qdxyzs(3,jix,iloc) = d2qdxyzs(3,jix,iloc) + d1trm*yji
                d2qdxyzs(3,jiy,iloc) = d2qdxyzs(3,jiy,iloc) + d1trm*xji
              elseif (ndim.eq.1) then
                d2qdxyzs(1,jix,iloc) = d2qdxyzs(1,jix,iloc) + 2.0_dp*d1trm*xji
              endif
            endif
          endif
        endif
      endif
!
!  End of loop over neighbours of i
!
      ni = ni + 1
    enddo
  enddo
!
!  Global sum of charges
!
  call sumall(qf,qa,numat,"getbocharged","qf")
!
!  Set up charges for asymmetric unit
!
  qf(1:numat) = qa(1:numat)
  do i = 1,nasym
    qa(i) = qf(nrel2(i))
  enddo
!********************************
!  Debugging output of charges  *
!********************************
  if (ioproc.and.index(keyword,'debu').ne.0) then
    write(ioout,'(/,''  BO charges of atoms :'',/)')
    write(ioout,'(''  Atom no.       Charge'')')
    do i = 1,nasym
      write(ioout,'(i8,5x,f12.8)') i,qa(i)
    enddo
  endif
!
!  Free local memory
!
  deallocate(zneigh,stat=status)
  if (status/=0) call deallocate_error('getBOcharged','zneigh')
  deallocate(yneigh,stat=status)
  if (status/=0) call deallocate_error('getBOcharged','yneigh')
  deallocate(xneigh,stat=status)
  if (status/=0) call deallocate_error('getBOcharged','xneigh')
  deallocate(rneigh,stat=status)
  if (status/=0) call deallocate_error('getBOcharged','rneigh')
  deallocate(ineigh,stat=status)
  if (status/=0) call deallocate_error('getBOcharged','ineigh')
  deallocate(neighno,stat=status)
  if (status/=0) call deallocate_error('getBOcharged','neighno')
  deallocate(rBOcutmax,stat=status)
  if (status/=0) call deallocate_error('getBOcharged','rBOcutmax')
  deallocate(nneigh,stat=status)
  if (status/=0) call deallocate_error('getBOcharged','nneigh')
!
  t2 = g_cpu_time()
  tbondorder = tbondorder + t2 - t1
#ifdef TRACE
  call trace_out('getbocharged')
#endif
!
  return
  end
