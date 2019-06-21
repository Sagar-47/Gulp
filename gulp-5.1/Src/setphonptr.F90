  subroutine setphonptr
!
!  Sets up pointers for atoms involved in phonon calculations
!
!   7/12 Created from code in phonon
!  11/12 Parallel distribution code separated
!  11/16 Distributed second derivative modifications made
!  12/16 Definition of pointers for cores & shells on node corrected
!   4/17 Parallel distribution added for 2D case
!   7/17 Unused module variables removed
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
!  Copyright Curtin University 2017
!
!  Julian Gale, CIC, Curtin University, July 2017
!
  use configurations,  only : nregionno, lbsmat
  use control
  use current
  use element
  use parallel
  use phononatoms
  use shells
  implicit none
!
!  Local variables
!
  integer(i4)                                       :: i
  integer(i4)                                       :: ii
!************************************************
!  Find number of atoms for phonon calculation  *
!************************************************
!
!  nphonat  = number of atoms involved in phonon calculation
!  nphonatb = number of breathing shells involved in phonon calculation
!  nphonatc = number of cores involved in phonon calculation
!  nphonats = number of shells involved in phonon calculation
!  nphonatptr = pointer from nphonat to numat reference
!  nphonatcptr = pointer from nphonatc to ncore reference
!  nphonatsptr = pointer from nphonats to nshel reference
!  nphonatrptr = pointer from numat to nphonat reference
!  nphonatrcptr = pointer from ncore to nphonatc reference
!  nphonatrsptr = pointer from nshel to nphonats reference
!
!  Names as above with "onnode" refer to the local node component
!
  if (ndim.eq.2) then
! DEBUG - need to handle parallel regions case
    nphonat = 0
    nphonatb = 0
    nphonatc = 0
    nphonats = 0
    do i = 1,numat
      if (nregionno(nsft+nrelat(i)).eq.1) then
        nphonat = nphonat + 1
        if (nat(i).gt.maxele) then
          nphonats = nphonats + 1
          nphonatsptr(nphonats) = i
          nphonatrsptr(i) = nphonats
        else
          nphonatc = nphonatc + 1
          nphonatcptr(nphonatc) = i
          nphonatrcptr(i) = nphonatc
        endif
        nphonatptr(nphonat) = i
        nphonatrptr(i) = nphonat
        if (lbsmat(nsft+nrelat(i))) then
          nphonatb = nphonatb + 1
        endif
      else
        nphonatrptr(i) = 0
      endif
    enddo
    if (nprocs.gt.1) then
      nphonatonnode = 0
      nphonatonnodeb = 0
      do ii = 1,natomsonnode
        i = node2atom(ii)
        if (nregionno(nsft+nrelat(i)).eq.1) then
          nphonatonnode = nphonatonnode + 1
          nphonatonnodeptr(nphonatonnode) = i
          nphonatonnoderptr(i) = nphonatonnode
          if (lbsmat(nsft+nrelat(i))) then
            nphonatonnodeb = nphonatonnodeb + 1
          endif
        endif
      enddo
      nphonatonnodec = 0
      do ii = 1,ncoreonnode
        i = node2atom(ncoonnodeptr(ii))
        if (nregionno(nsft+nrelat(i)).eq.1) then
          nphonatonnodec = nphonatonnodec + 1
          nphonatonnodecptr(nphonatonnodec) = i
          nphonatonnodercptr(i) = nphonatonnodec
        endif
      enddo
      nphonatonnodes = 0
      do ii = 1,nshellonnode
        i = node2atom(nshonnodeptr(ii))
        if (nregionno(nsft+nrelat(i)).eq.1) then
          nphonatonnodes = nphonatonnodes + 1
          nphonatonnodesptr(nphonatonnodes) = i
          nphonatonnodersptr(i) = nphonatonnodes
        endif
      enddo
    else
      nphonatonnode  = nphonat
      nphonatonnodeb = nphonatb
      nphonatonnodec = nphonatc
      nphonatonnodes = nphonats
      nphonatonnodeptr(1:numat) = nphonatptr(1:numat)
      nphonatonnoderptr(1:numat) = nphonatrptr(1:numat)
      nphonatonnodecptr(1:ncore) = nphonatcptr(1:ncore)
      nphonatonnodercptr(1:ncore) = nphonatrcptr(1:ncore)
      nphonatonnodesptr(1:nshell) = nphonatsptr(1:nshell)
      nphonatonnodersptr(1:numat) = nphonatrsptr(1:numat)
    endif
  else
    nphonat = numat
    nphonatb = 0
    do i = 1,numat
      nphonatptr(i) = i
      nphonatrptr(i) = i
      if (lbsmat(nsft+nrelat(i))) then
        nphonatb = nphonatb + 1
      endif
    enddo
    nphonatc = ncore
    nphonats = nshell
    do i = 1,ncore
      nphonatcptr(i) = i
      nphonatrcptr(i) = i
    enddo
    do i = 1,nshell
      nphonatsptr(i) = i + ncore
      nphonatrsptr(ncore+i) = i
    enddo
    if (nprocs.gt.1) then
      nphonatonnode = natomsonnode
      nphonatonnodeb = 0
      do i = 1,natomsonnode
        nphonatonnodeptr(i) = node2atom(i)
        nphonatonnoderptr(node2atom(i)) = i
        if (lbsmat(nsft+nrelat(node2atom(i)))) then
          nphonatonnodeb = nphonatonnodeb + 1
        endif
      enddo
      nphonatonnodec = ncoreonnode
      nphonatonnodes = nshellonnode
      do i = 1,ncoreonnode
        nphonatonnodecptr(i) = node2atom(ncoonnodeptr(i))
        nphonatonnodercptr(node2atom(ncoonnodeptr(i))) = i
      enddo
      do i = 1,nshellonnode
        nphonatonnodesptr(i) = node2atom(nshonnodeptr(i))
        nphonatonnodersptr(node2atom(nshonnodeptr(i))) = i
      enddo
    else
      nphonatonnode  = nphonat
      nphonatonnodeb = nphonatb
      nphonatonnodec = nphonatc
      nphonatonnodes = nphonats
      nphonatonnodeptr(1:numat) = nphonatptr(1:numat)
      nphonatonnoderptr(1:numat) = nphonatrptr(1:numat)
      nphonatonnodecptr(1:ncore) = nphonatcptr(1:ncore)
      nphonatonnodercptr(1:ncore) = nphonatrcptr(1:ncore)
      nphonatonnodesptr(1:nshell) = nphonatsptr(1:nshell)
      nphonatonnodersptr(1:numat) = nphonatrsptr(1:numat)
    endif
  endif
!
  return
  end
