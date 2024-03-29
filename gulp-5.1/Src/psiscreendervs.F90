  subroutine psiscreendervs(i,j,npartial,partial,xij,yij,zij,energy,Sij,strfct,lstr)
!
!  Calculates the derivatives from the psibaskes contribution allowing for the
!  screening factors. Symmetry adapted version called from manysd
!
!  On entry :
!
!  i         = asymmetric unit atom that is the focus of derivatives
!  j         = full unit atom whose pair with i is being considered
!  energy    = energy contribution from potential
!  Sij       = total screening factor for i-j pair
!  npartial  = number of atoms contributing to partial screening of i-j pair
!  partial   = data type containing the information regarding the i-k-j trio for screening
!  strfct    = scale factor to correct strain derivatives
!  xij       = x component of i-j vector
!  yij       = y component of i-j vector
!  zij       = z component of i-j vector
!
!  On exit :
!
!  Derivatives added to appropriate arrays
!
!   7/18 Created from psiscreenderv
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
!  Julian Gale, CIC, Curtin University, July 2018
!
  use control,        only : latomicstress
  use current
  use derivatives
  use eam
  use numbers,        only : third
  use optimisation,   only : lopf, lfreeze
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),           intent(in)    :: i
  integer(i4),           intent(in)    :: j
  integer(i4),           intent(in)    :: npartial
  logical,               intent(in)    :: lstr
  real(dp),              intent(in)    :: energy
  real(dp),              intent(in)    :: Sij
  real(dp),              intent(in)    :: strfct
  real(dp),              intent(in)    :: xij
  real(dp),              intent(in)    :: yij
  real(dp),              intent(in)    :: zij
  type(screening_atoms), intent(inout) :: partial
!
!  Local variables
!
  integer(i4)                          :: k                   ! Atom index for k
  integer(i4)                          :: ka                  ! Atom index for k in asymmetric unit
  integer(i4)                          :: kl                  ! Looping index over strains
  integer(i4)                          :: klp                 ! Pointer from strain number to strain
  integer(i4)                          :: np                  ! Looping index over partial screening atom
  logical                              :: lopi                ! Optimisation flag for i
  logical                              :: lopk                ! Optimisation flag for k
  real(dp)                             :: derivk(3)           ! Workspace for derivatives of k
  real(dp)                             :: Sij_rest            ! Sij product, but without current Sikj contribution (and Silj for second derivatives)
  real(dp)                             :: strloc(6)           ! Local strain workspace array
  real(dp)                             :: xik                 ! Local copy of partial%sa_xik(np) for brevity
  real(dp)                             :: yik                 ! Local copy of partial%sa_yik(np) for brevity
  real(dp)                             :: zik                 ! Local copy of partial%sa_zik(np) for brevity
  real(dp)                             :: xjk                 ! Local copy of partial%sa_xjk(np) for brevity
  real(dp)                             :: yjk                 ! Local copy of partial%sa_yjk(np) for brevity
  real(dp)                             :: zjk                 ! Local copy of partial%sa_zjk(np) for brevity
#ifdef TRACE
  call trace_in('psiscreendervs')
#endif
!
!  Set optimisation flag for i
!
  lopi = (.not.lfreeze.or.lopf(i))
!
!  Is screening term = 0 or 1? If so then there are no derivatives.
!
  if (Sij.ne.0.0_dp.and.Sij.ne.1.0_dp) then
!
!  Now loop over partially screening atoms to construct full derivative
!
    do np = 1,npartial
      k = partial%sa_atom(np)
      ka = nrelat(k)
      lopk = (.not.lfreeze.or.lopf(ka))
!
!  Set local scalars to keep code compact
!
      xik = partial%sa_xik(np)
      yik = partial%sa_yik(np)
      zik = partial%sa_zik(np)
      xjk = partial%sa_xjk(np)
      yjk = partial%sa_yjk(np)
      zjk = partial%sa_zjk(np)
!
!  Compute Sij without the current Sikj : Since all Sikj must be > 0 then we can just divide
!
      Sij_rest = Sij/partial%sa_Sikj(np)
!
      partial%sa_drhototik(1,np) = energy*Sij_rest*partial%sa_dSikjdr(2,np)*xik
      partial%sa_drhototik(2,np) = energy*Sij_rest*partial%sa_dSikjdr(2,np)*yik
      partial%sa_drhototik(3,np) = energy*Sij_rest*partial%sa_dSikjdr(2,np)*zik
      partial%sa_drhototjk(1,np) = energy*Sij_rest*partial%sa_dSikjdr(3,np)*xjk
      partial%sa_drhototjk(2,np) = energy*Sij_rest*partial%sa_dSikjdr(3,np)*yjk
      partial%sa_drhototjk(3,np) = energy*Sij_rest*partial%sa_dSikjdr(3,np)*zjk
!
!  i-j contribution
!
      if (lopi) then
        xdrv(i) = xdrv(i) - energy*Sij_rest*partial%sa_dSikjdr(1,np)*xij
        ydrv(i) = ydrv(i) - energy*Sij_rest*partial%sa_dSikjdr(1,np)*yij
        zdrv(i) = zdrv(i) - energy*Sij_rest*partial%sa_dSikjdr(1,np)*zij
      endif
!
!  i-k contribution
!
      if (lopi) then
        xdrv(i) = xdrv(i) - energy*Sij_rest*partial%sa_dSikjdr(2,np)*xik
        ydrv(i) = ydrv(i) - energy*Sij_rest*partial%sa_dSikjdr(2,np)*yik
        zdrv(i) = zdrv(i) - energy*Sij_rest*partial%sa_dSikjdr(2,np)*zik
      endif
      if (lopk) then
!
!  Manipulate the contribution to the derivatives of k to those acting on the image of k in the asymmetric unit
!
        derivk(1) = energy*Sij_rest*partial%sa_dSikjdr(2,np)*xik
        derivk(2) = energy*Sij_rest*partial%sa_dSikjdr(2,np)*yik
        derivk(3) = energy*Sij_rest*partial%sa_dSikjdr(2,np)*zik
!
        call symdervrot(k,nrel2(ka),derivk)
!
        xdrv(ka) = xdrv(ka) + derivk(1)
        ydrv(ka) = ydrv(ka) + derivk(2)
        zdrv(ka) = zdrv(ka) + derivk(3)
      endif
!
!  Strains
!
      if (lstr) then
        strloc(1) = energy*Sij_rest*partial%sa_dSikjdr(1,np)*xij*xij + &
                    energy*Sij_rest*partial%sa_dSikjdr(2,np)*xik*xik + &
                    energy*Sij_rest*partial%sa_dSikjdr(3,np)*xjk*xjk
        strloc(2) = energy*Sij_rest*partial%sa_dSikjdr(1,np)*yij*yij + &
                    energy*Sij_rest*partial%sa_dSikjdr(2,np)*yik*yik + &
                    energy*Sij_rest*partial%sa_dSikjdr(3,np)*yjk*yjk
        strloc(3) = energy*Sij_rest*partial%sa_dSikjdr(1,np)*zij*zij + &
                    energy*Sij_rest*partial%sa_dSikjdr(2,np)*zik*zik + &
                    energy*Sij_rest*partial%sa_dSikjdr(3,np)*zjk*zjk
        strloc(4) = energy*Sij_rest*partial%sa_dSikjdr(1,np)*yij*zij + &
                    energy*Sij_rest*partial%sa_dSikjdr(2,np)*yik*zik + &
                    energy*Sij_rest*partial%sa_dSikjdr(3,np)*yjk*zjk
        strloc(5) = energy*Sij_rest*partial%sa_dSikjdr(1,np)*xij*zij + &
                    energy*Sij_rest*partial%sa_dSikjdr(2,np)*xik*zik + &
                    energy*Sij_rest*partial%sa_dSikjdr(3,np)*xjk*zjk
        strloc(6) = energy*Sij_rest*partial%sa_dSikjdr(1,np)*xij*yij + &
                    energy*Sij_rest*partial%sa_dSikjdr(2,np)*xik*yik + &
                    energy*Sij_rest*partial%sa_dSikjdr(3,np)*xjk*yjk
!
        do kl = 1,nstrains
          klp = nstrptr(kl)
          rstrd(kl) = rstrd(kl) + strloc(klp)*strfct
        enddo
! DEBUG - atomic stresses not computed by manysd at present
!        if (latomicstress) then
!          do kl = 1,nstrains
!            klp = nstrptr(kl)
!            atomicstress(kl,i) = atomicstress(kl,i) + third*strloc(klp)
!            atomicstress(kl,j) = atomicstress(kl,j) + third*strloc(klp)
!            atomicstress(kl,k) = atomicstress(kl,k) + third*strloc(klp)
!          enddo
!        endif
      endif
    enddo
  endif
#ifdef TRACE
  call trace_out('psiscreendervs')
#endif
!
  return
  end
