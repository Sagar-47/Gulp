  subroutine einsteinpd(eeinstein,lgrad1,lgrad2)
!
!  Subroutine for calculating the energy and derivatives due to the Einstein model
!  Distributed memory version. 
!
!  11/12 Created from einstein
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
  use current
  use configurations
  use derivatives
  use energies,       only : siteenergy, eregion2region
  use optimisation,   only : lopf, lfreeze
  use parallel
  use symmetry
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  logical,  intent(in)  :: lgrad1
  logical,  intent(in)  :: lgrad2
  real(dp), intent(out) :: eeinstein
!
!  Local variables
!
  integer(i4)           :: i
  integer(i4)           :: iloc
  integer(i4)           :: ix
  integer(i4)           :: iy
  integer(i4)           :: iz
  integer(i4)           :: ni
  integer(i4)           :: nr
  integer(i4)           :: nregioni
  real(dp)              :: deltae
  real(dp)              :: ki
  real(dp)              :: r2
  real(dp)              :: xf
  real(dp)              :: yf
  real(dp)              :: zf
  real(dp)              :: xi
  real(dp)              :: yi
  real(dp)              :: zi
!
!  Initialise integral of force x distance
!
  eeinstein = 0.0_dp
!
!  If there are no Einstein model atoms return
!
  if (.not.leinstein) return
#ifdef TRACE
  call trace_in('einsteinpd')
#endif
!
!  Initialise second derivative indices
!
  if (lgrad2) then
    ix = -2
    iy = -1
    iz = 0
  endif
!
!  Loop over asymmetric unit performing integral
!
  do iloc = 1,natomsonnode
    i = node2atom(iloc)
    ni = nrelat(i)
    nregioni = nregionno(nsft+ni)
    if (lopf(ni).or..not.lfreeze) then
      if (lgrad2) then
        ix = ix + 3
        iy = iy + 3
        iz = iz + 3
      endif
      if (leinsteinat(nsft+ni)) then
!
!  Find nearest image in fraction coordinates - must use image of symmetry
!  related atom to avoid problems where lattice points are only defined for
!  the asymmetric unit.
!
        xf = xafrac(ni) - xeinsteinat(nsft + ni)
        yf = yafrac(ni) - yeinsteinat(nsft + ni)
        zf = zafrac(ni) - zeinsteinat(nsft + ni)
        xf = mod(xf+10.0_dp,1.0_dp)
        yf = mod(yf+10.0_dp,1.0_dp)
        zf = mod(zf+10.0_dp,1.0_dp)
        if (xf.gt.0.5_dp) xf = xf - 1.0_dp
        if (yf.gt.0.5_dp) yf = yf - 1.0_dp
        if (zf.gt.0.5_dp) zf = zf - 1.0_dp
!
!  Convert coordinate difference into Cartesian frame
!
        xi = xf*rv(1,1) + yf*rv(1,2) + zf*rv(1,3)
        yi = xf*rv(2,1) + yf*rv(2,2) + zf*rv(2,3)
        zi = xf*rv(3,1) + yf*rv(3,2) + zf*rv(3,3)
!
!  If symmetry is being used, then rotate difference
!
        if (i.ne.ni) then
          xf = xi
          yf = yi
          zf = zi
          nr = nrotop(i)
          xi = xf*rop(1,1,nr) + yf*rop(1,2,nr) + zf*rop(1,3,nr)
          yi = xf*rop(2,1,nr) + yf*rop(2,2,nr) + zf*rop(2,3,nr)
          zi = xf*rop(3,1,nr) + yf*rop(3,2,nr) + zf*rop(3,3,nr)
        endif
!
!  Evaluate energy and derivatives
!
        ki = keinsteinat(nsft + ni)
        r2 = xi*xi + yi*yi + zi*zi
        deltae = 0.5_dp*ki*r2
        eeinstein = eeinstein + deltae
!
        eregion2region(nregioni,nregioni) = eregion2region(nregioni,nregioni) + deltae
        siteenergy(i) = siteenergy(i) + deltae
        if (lgrad1) then
          xdrv(i) = xdrv(i) + ki*xi
          ydrv(i) = ydrv(i) + ki*yi
          zdrv(i) = zdrv(i) + ki*zi
          xregdrv(nregioni) = xregdrv(nregioni) + ki*xi
          yregdrv(nregioni) = yregdrv(nregioni) + ki*yi
          zregdrv(nregioni) = zregdrv(nregioni) + ki*zi
          if (lgrad2) then
            derv2d(ix) = derv2d(ix) + ki
            derv2d(iy) = derv2d(iy) + ki
            derv2d(iz) = derv2d(iz) + ki
          endif
        endif
      endif
    endif
  enddo
#ifdef TRACE
  call trace_out('einsteinpd')
#endif
!
  return
  end
