  subroutine psiscreenderv(i,j,npartial,partial,xij,yij,zij,energy,Sij,lstr)
!
!  Calculates the derivatives from the psibaskes contribution allowing for the
!  screening factors
!
!  On entry :
!
!  energy    = energy contribution from potential
!  Sij       = total screening factor for i-j pair
!  npartial  = number of atoms contributing to partial screening of i-j pair
!  partial   = data type containing the information regarding the i-k-j trio for screening
!  xij       = x component of i-j vector
!  yij       = y component of i-j vector
!  zij       = z component of i-j vector
!
!  On exit :
!
!  Derivatives added to appropriate arrays
!
!   8/14 Created from meamtotalscreenderv
!   9/14 Bug fixed in referencing of partial S derivatives
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
  use configurations, only : nregionno
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
  real(dp),              intent(in)    :: xij
  real(dp),              intent(in)    :: yij
  real(dp),              intent(in)    :: zij
  type(screening_atoms), intent(inout) :: partial
!
!  Local variables
!
  integer(i4)                          :: k                   ! Atom index for k
  integer(i4)                          :: kl                  ! Looping index over strains
  integer(i4)                          :: klp                 ! Pointer from strain number to strain
  integer(i4)                          :: np                  ! Looping index over partial screening atom
  integer(i4)                          :: nregioni            ! Region number for i
  integer(i4)                          :: nregionj            ! Region number for j
  integer(i4)                          :: nregionk            ! Region number for k
  logical                              :: lopi                ! Optimisation flag for i
  logical                              :: lopj                ! Optimisation flag for j
  logical                              :: lopk                ! Optimisation flag for k
  real(dp)                             :: Sij_rest            ! Sij product, but without current Sikj contribution (and Silj for second derivatives)
  real(dp)                             :: strloc(6)           ! Local strain workspace array
  real(dp)                             :: xik                 ! Local copy of partial%sa_xik(np) for brevity
  real(dp)                             :: yik                 ! Local copy of partial%sa_yik(np) for brevity
  real(dp)                             :: zik                 ! Local copy of partial%sa_zik(np) for brevity
  real(dp)                             :: xjk                 ! Local copy of partial%sa_xjk(np) for brevity
  real(dp)                             :: yjk                 ! Local copy of partial%sa_yjk(np) for brevity
  real(dp)                             :: zjk                 ! Local copy of partial%sa_zjk(np) for brevity
#ifdef TRACE
  call trace_in('psiscreenderv')
#endif
!
!  Set terms for i & j
!
  lopi = (.not.lfreeze.or.lopf(nrelat(i)))
  lopj = (.not.lfreeze.or.lopf(nrelat(j)))
  nregioni = nregionno(nsft+nrelat(i))
  nregionj = nregionno(nsft+nrelat(j))
!
!  Is screening term = 0 or 1? If so then there are no derivatives.
!
  if (Sij.ne.0.0_dp.and.Sij.ne.1.0_dp) then
!
!  Now loop over partially screening atoms to construct full derivative
!
    do np = 1,npartial
      k = partial%sa_atom(np)
      lopk = (.not.lfreeze.or.lopf(nrelat(k)))
      nregionk = nregionno(nsft+nrelat(k))
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
      if (lopj) then
        xdrv(j) = xdrv(j) + energy*Sij_rest*partial%sa_dSikjdr(1,np)*xij
        ydrv(j) = ydrv(j) + energy*Sij_rest*partial%sa_dSikjdr(1,np)*yij
        zdrv(j) = zdrv(j) + energy*Sij_rest*partial%sa_dSikjdr(1,np)*zij
      endif
      if (nregioni.ne.nregionj) then
        xregdrv(nregioni) = xregdrv(nregioni) - energy*Sij_rest*partial%sa_dSikjdr(1,np)*xij
        yregdrv(nregioni) = yregdrv(nregioni) - energy*Sij_rest*partial%sa_dSikjdr(1,np)*yij
        zregdrv(nregioni) = zregdrv(nregioni) - energy*Sij_rest*partial%sa_dSikjdr(1,np)*zij
        xregdrv(nregionj) = xregdrv(nregionj) + energy*Sij_rest*partial%sa_dSikjdr(1,np)*xij
        yregdrv(nregionj) = yregdrv(nregionj) + energy*Sij_rest*partial%sa_dSikjdr(1,np)*yij
        zregdrv(nregionj) = zregdrv(nregionj) + energy*Sij_rest*partial%sa_dSikjdr(1,np)*zij
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
        xdrv(k) = xdrv(k) + energy*Sij_rest*partial%sa_dSikjdr(2,np)*xik
        ydrv(k) = ydrv(k) + energy*Sij_rest*partial%sa_dSikjdr(2,np)*yik
        zdrv(k) = zdrv(k) + energy*Sij_rest*partial%sa_dSikjdr(2,np)*zik
      endif
      if (nregioni.ne.nregionk) then
        xregdrv(nregioni) = xregdrv(nregioni) - energy*Sij_rest*partial%sa_dSikjdr(2,np)*xik
        yregdrv(nregioni) = yregdrv(nregioni) - energy*Sij_rest*partial%sa_dSikjdr(2,np)*yik
        zregdrv(nregioni) = zregdrv(nregioni) - energy*Sij_rest*partial%sa_dSikjdr(2,np)*zik
        xregdrv(nregionk) = xregdrv(nregionk) + energy*Sij_rest*partial%sa_dSikjdr(2,np)*xik
        yregdrv(nregionk) = yregdrv(nregionk) + energy*Sij_rest*partial%sa_dSikjdr(2,np)*yik
        zregdrv(nregionk) = zregdrv(nregionk) + energy*Sij_rest*partial%sa_dSikjdr(2,np)*zik
      endif
!
!  j-k contribution
!
      if (lopj) then
        xdrv(j) = xdrv(j) - energy*Sij_rest*partial%sa_dSikjdr(3,np)*xjk
        ydrv(j) = ydrv(j) - energy*Sij_rest*partial%sa_dSikjdr(3,np)*yjk
        zdrv(j) = zdrv(j) - energy*Sij_rest*partial%sa_dSikjdr(3,np)*zjk
      endif
      if (lopk) then
        xdrv(k) = xdrv(k) + energy*Sij_rest*partial%sa_dSikjdr(3,np)*xjk
        ydrv(k) = ydrv(k) + energy*Sij_rest*partial%sa_dSikjdr(3,np)*yjk
        zdrv(k) = zdrv(k) + energy*Sij_rest*partial%sa_dSikjdr(3,np)*zjk
      endif
      if (nregionj.ne.nregionk) then
        xregdrv(nregionj) = xregdrv(nregionj) - energy*Sij_rest*partial%sa_dSikjdr(3,np)*xjk
        yregdrv(nregionj) = yregdrv(nregionj) - energy*Sij_rest*partial%sa_dSikjdr(3,np)*yjk
        zregdrv(nregionj) = zregdrv(nregionj) - energy*Sij_rest*partial%sa_dSikjdr(3,np)*zjk
        xregdrv(nregionk) = xregdrv(nregionk) + energy*Sij_rest*partial%sa_dSikjdr(3,np)*xjk
        yregdrv(nregionk) = yregdrv(nregionk) + energy*Sij_rest*partial%sa_dSikjdr(3,np)*yjk
        zregdrv(nregionk) = zregdrv(nregionk) + energy*Sij_rest*partial%sa_dSikjdr(3,np)*zjk
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
          rstrd(kl) = rstrd(kl) + strloc(klp)
        enddo
        if (latomicstress) then
          do kl = 1,nstrains
            klp = nstrptr(kl)
            atomicstress(kl,i) = atomicstress(kl,i) + third*strloc(klp)
            atomicstress(kl,j) = atomicstress(kl,j) + third*strloc(klp)
            atomicstress(kl,k) = atomicstress(kl,k) + third*strloc(klp)
          enddo
        endif
      endif
    enddo
  endif
#ifdef TRACE
  call trace_out('psiscreenderv')
#endif
!
  return
  end
