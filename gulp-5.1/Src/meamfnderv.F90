  subroutine meamfnderv(neamfn,m,scrho,rhototal,eeam,rscrho,rscrho3,rscrho5,lgrad1,lgrad2,lgrad3)
!
!  Calculates the derivatives of the MEAM function.
!
!  On entry :
!
!  neamfn   = pointer to function type
!  m        = EAM function species number
!  scrho    = density components at atom 
!  lgrad1   = if .true. calculate first derivatives 
!  lgrad2   = if .true. calculate second derivatives 
!  lgrad3   = if .true. calculate third derivatives 
!
!  On exit :
!
!  eeam     = local contribution to energy
!  rhototal = total rho from sum of component contributions
!  rscrho   = term needed for energy/first derivatives
!  rscrho3  = term needed for second derivatives
!  rscrho5  = term needed for third derivatives
!
!  NB: In the above terms the negative sign of the embedding energy is ignored
!      since this is handled in the calling routine.
!
!  12/08 Created from eamfnderv
!  12/08 rho0 added to Baskes form of functional
!   1/09 Total rho now returned in argument list
!   3/09 l replaced as argument directly by m
!   3/09 Return arguments sign changed so that they directly represent the energy
!        and its derivatives.
!  12/11 Sign of Baskes functional reversed to be consistent with papers.
!   6/13 Igarashi functional added
!   8/13 Spline functional added
!   8/14 Baskes functional modified to allow for Ec and A being separated
!  12/14 Derivatives corrected for alloy parameter
!   2/15 Bug in second and third derivatives of Johnson function corrected
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
!  Julian Gale, CIC, Curtin University, February 2018
!
  use datatypes
  use eam,           only : eamfnpar, neampower, eamalloy
  use eam,           only : eamfnnumeric, eamfnnumericdrho, neamfnnumeric
  use eam,           only : eamfnnumeric1, eamfnnumeric2, eamfnnumeric3
  use eam,           only : eamfnnumeric4, eamfnnumeric5, eamfnnumeric6
  use eam,           only : eamfnnumeric7, eamfnnumeric8
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: m
  integer(i4), intent(in)  :: neamfn
  logical,     intent(in)  :: lgrad1
  logical,     intent(in)  :: lgrad2
  logical,     intent(in)  :: lgrad3
  real(dp),    intent(out) :: eeam
  real(dp),    intent(in)  :: scrho(*)
  real(dp),    intent(out) :: rhototal
  real(dp),    intent(out) :: rscrho
  real(dp),    intent(out) :: rscrho3
  real(dp),    intent(out) :: rscrho5
!
!  Local variables
!
  integer(i4)              :: k
  integer(i4)              :: kk
  real(dp)                 :: AE0
  real(dp)                 :: alpha
  real(dp)                 :: beta
  real(dp)                 :: gamma
  real(dp)                 :: delta
  real(dp)                 :: atrm
  real(dp)                 :: btrm
  real(dp)                 :: c0
  real(dp)                 :: c1
  real(dp)                 :: c2
  real(dp)                 :: c3
  real(dp)                 :: c4
  real(dp)                 :: drho
  real(dp)                 :: drho2
  real(dp)                 :: drho3
  real(dp)                 :: drho4
  real(dp)                 :: dxdrho
  real(dp)                 :: d1rootrho
  real(dp)                 :: d2rootrho
  real(dp)                 :: d2xdrho2
  real(dp)                 :: d3xdrho3
  real(dp)                 :: dydrho
  real(dp)                 :: d2ydrho2
  real(dp)                 :: d3ydrho3
  real(dp)                 :: exptrm
  real(dp)                 :: f0
  real(dp)                 :: f1
  real(dp)                 :: logrho
  real(dp)                 :: phi0
  real(dp)                 :: pp
  real(dp)                 :: rho13
  real(dp)                 :: rho23
  real(dp)                 :: rho43
  real(dp)                 :: rho53
  real(dp)                 :: rhoe
  real(dp)                 :: rhoi
  real(dp)                 :: rhon
  real(dp)                 :: rhorr
  real(dp)                 :: rhorrn
  real(dp)                 :: rhorrnx
  real(dp)                 :: rhorrny
  real(dp)                 :: rhorrn1
  real(dp)                 :: rhorrn2
  real(dp)                 :: rhorrn3
  real(dp)                 :: rhorrn1x
  real(dp)                 :: rhorrn2x
  real(dp)                 :: rhorrn3x
  real(dp)                 :: rhorrn1y
  real(dp)                 :: rhorrn2y
  real(dp)                 :: rhorrn3y
  real(dp)                 :: rk
  real(dp)                 :: rlogtrm
  real(dp)                 :: rnp
  real(dp)                 :: rnpx
  real(dp)                 :: rnpy
  real(dp)                 :: rootk
  real(dp)                 :: rootrho
  real(dp)                 :: trm1
#ifdef TRACE
  call trace_in('meanfnderv')
#endif
!
  if (m.gt.0) then
    if (scrho(1).gt.1.0d-12) then
      f0 = eamfnpar(1,m)
!
!  Compute total density from MEAM components
!
      call meamtotalrho(m,scrho,rhototal)
!
!  Evaluate functional derivatives
!
      if (neamfn.eq.1) then
        rootrho = sqrt(rhototal)
        eeam = - f0*rootrho
        if (lgrad1) then
          rscrho = - 0.5_dp/rootrho
          if (lgrad2) then
            rscrho3 = - 2.0_dp*f0*rscrho**3
            if (lgrad3) then
              rscrho5 = - 1.5_dp*rscrho3/rhototal
            endif
          endif
          rscrho = f0*rscrho
        endif
      elseif (neamfn.eq.2) then
        rnp = 1.0_dp/dble(neampower)
        rhon = (rhototal)**(rnp - 1.0_dp)
        eeam = - f0*rhon*rhototal
        if (lgrad1) then
          rscrho = - f0*rnp*rhon
          if (lgrad2) then
            rscrho3 = rscrho*(rnp - 1.0_dp)/rhototal
            if (lgrad3) then
              rscrho5 = rscrho3*(rnp - 2.0_dp)/rhototal
            endif
          endif
        endif
      elseif (neamfn.eq.3) then
        rnp = 1.0_dp/dble(neampower)
        rhoi = rhototal
        f1 = eamfnpar(2,m)
        rhoe = 1.0_dp/eamfnpar(3,m)
        rhorr = rhoi*rhoe
        rhorrn3 = rhorr**(rnp - 3.0_dp)
        rhorrn2 = rhorrn3*rhorr
        rhorrn1 = rhorrn2*rhorr
        rhorrn = rhorrn1*rhorr
        rlogtrm = log(rhorrn)
        eeam = - (f0*(1.0_dp - rlogtrm)*rhorrn - f1*rhorr)
        if (lgrad1) then
          rscrho = f0*rnp*rhoe*rlogtrm*rhorrn1 + f1*rhoe
          if (lgrad2) then
            rscrho3 = f0*rnp*rhoe*rhoe*rhorrn2*((rnp - 1.0_dp)*rlogtrm + rnp)
            if (lgrad3) then
              rscrho5 = f0*rnp*rhoe*rhoe*rhoe*rhorrn3*((rnp - 1.0_dp)*(rnp - 2.0_dp)*rlogtrm - rnp*(2.0_dp*rnp - 3.0_dp))
            endif
          endif
        endif
      elseif (neamfn.eq.4) then
        pp = rhototal/eamfnnumericdrho(m) + 1.0_dp 
        kk = nint(pp)
        kk = max(1,min(kk,neamfnnumeric(m)-1))
        pp = pp - dble(kk)
        eeam = (((eamfnnumeric3(kk,m)*pp + eamfnnumeric2(kk,m))*pp + eamfnnumeric1(kk,m))*pp + eamfnnumeric(kk,m))
        if (lgrad1) then
          rscrho = ((eamfnnumeric6(kk,m)*pp + eamfnnumeric5(kk,m))*pp + eamfnnumeric4(kk,m))
          if (lgrad2) then
            rscrho3 = (eamfnnumeric8(kk,m)*pp + eamfnnumeric7(kk,m))
            if (lgrad3) then
              rscrho5 = 0.0_dp
            endif
          endif
        endif
      elseif (neamfn.eq.5) then
!
!  Johnson
!
        alpha = eamfnpar(4,m)
        beta = eamfnpar(5,m)
        gamma = eamfnpar(6,m)
        rhoi = rhototal
        f1 = eamfnpar(2,m)
        rhoe = 1.0_dp/eamfnpar(3,m)
        rhorr = rhoi*rhoe
!
        rnpx = alpha/beta
        rhorrn3x = rhorr**(rnpx - 3.0_dp)
        rhorrn2x = rhorrn3x*rhorr
        rhorrn1x = rhorrn2x*rhorr
        rhorrnx = rhorrn1x*rhorr
        rlogtrm = log(rhorrnx)
!
        rnpy = gamma/beta
        rhorrn3y = rhorr**(rnpy - 3.0_dp)
        rhorrn2y = rhorrn3y*rhorr
        rhorrn1y = rhorrn2y*rhorr
        rhorrny = rhorrn1y*rhorr
!
        eeam = - f0*(1.0_dp - rlogtrm)*rhorrnx - f1*rhorrny
        if (lgrad1) then
          dxdrho = rnpx*rhorrn1x*rhoe
          dydrho = rnpy*rhorrn1y*rhoe
          rscrho = f0*rlogtrm*dxdrho - f1*dydrho
          if (lgrad2) then
            d2xdrho2 = rnpx*(rnpx - 1.0_dp)*rhorrn2x*rhoe*rhoe
            d2ydrho2 = rnpy*(rnpy - 1.0_dp)*rhorrn2y*rhoe*rhoe
            rscrho3 = f0*(rlogtrm*d2xdrho2 + (dxdrho**2)/rhorrnx) - f1*d2ydrho2
            if (lgrad3) then
              d3xdrho3 = rnpx*(rnpx - 1.0_dp)*(rnpx - 2.0_dp)*rhorrn3x*rhoe*rhoe*rhoe
              d3ydrho3 = rnpy*(rnpy - 1.0_dp)*(rnpy - 2.0_dp)*rhorrn3y*rhoe*rhoe*rhoe
              rscrho5 = f0*(rlogtrm*d3xdrho3 + 3.0_dp*d2xdrho2*dxdrho/rhorrnx - (dxdrho**3/rhorrnx**2))  &
                        - f1*d3ydrho3
            endif
          endif
        endif
      elseif (neamfn.eq.6) then
!
!  Glue
!
!  First find correct density range and then calculate appropriate polynomial functional
!
        rhoi = rhototal
        if (rhoi.lt.eamfnpar(1,m)) then
!
!  Low density range
!
          rhoe = eamfnpar(1,m)
          drho = rhoi - rhoe
          drho2 = drho**2
          drho3 = drho2*drho
          drho4 = drho3*drho
          c4 = eamfnpar(3,m)
          c3 = eamfnpar(4,m)
          c2 = eamfnpar(5,m)
          c1 = eamfnpar(6,m)
          c0 = eamfnpar(7,m)
          eeam = - (c4*drho4 + c3*drho3 + c2*drho2 + c1*drho + c0)
          if (lgrad1) then
            rscrho = - (4.0_dp*c4*drho3 + 3.0_dp*c3*drho2 + 2.0_dp*c2*drho + c1)
            if (lgrad2) then
              rscrho3 = - (12.0_dp*c4*drho2 + 6.0_dp*c3*drho + 2.0_dp*c2)
              if (lgrad3) then
                rscrho5 = - (24.0_dp*c4*drho + 6.0_dp*c3)
              endif 
            endif 
          endif 
        elseif (rhoi.ge.eamfnpar(2,m)) then
!
!  High density range
!
          rhoe = eamfnpar(2,m)
          drho = rhoi - rhoe
          drho2 = drho**2
          drho3 = drho2*drho
          c3 = eamfnpar(13,m)
          c2 = eamfnpar(14,m)
          c1 = eamfnpar(15,m)
          c0 = eamfnpar(16,m)
          eeam = - (c3*drho3 + c2*drho2 + c1*drho + c0)
          if (lgrad1) then
            rscrho = - (3.0_dp*c3*drho2 + 2.0_dp*c2*drho + c1)
            if (lgrad2) then
              rscrho3 = - (6.0_dp*c3*drho + 2.0_dp*c2)
              if (lgrad3) then
                rscrho5 = - 6.0_dp*c3
              endif 
            endif 
          endif 
        else
!
!  Mid density range
!
          rhoe = eamfnpar(2,m)
          drho = rhoi - rhoe
          drho2 = drho**2
          drho3 = drho2*drho
          drho4 = drho3*drho
          c4 = eamfnpar(8,m)
          c3 = eamfnpar(9,m)
          c2 = eamfnpar(10,m)
          c1 = eamfnpar(11,m)
          c0 = eamfnpar(12,m)
          eeam = - (c4*drho4 + c3*drho3 + c2*drho2 + c1*drho + c0)
          if (lgrad1) then
            rscrho = - (4.0_dp*c4*drho3 + 3.0_dp*c3*drho2 + 2.0_dp*c2*drho + c1)
            if (lgrad2) then
              rscrho3 = - (12.0_dp*c4*drho2 + 6.0_dp*c3*drho + 2.0_dp*c2)
              if (lgrad3) then
                rscrho5 = - (24.0_dp*c4*drho + 6.0_dp*c3)
              endif 
            endif 
          endif 
        endif
      elseif (neamfn.eq.7) then
!
!  Foiles
!
        c1 = eamfnpar(1,m)
        c2 = eamfnpar(2,m)
        c3 = eamfnpar(3,m)
        c4 = eamfnpar(4,m)
        rho53 = rhototal**(5.0_dp/3.0_dp)
        eeam = - c1*rhototal**2 - c2*rhototal - c3*rho53/(c4 + rhototal)
        if (lgrad1) then
          rho23 = rhototal**(2.0_dp/3.0_dp)
          rscrho = - 2.0_dp*c1*rhototal - c2 - c3*((5.0_dp/3.0_dp)*rho23*(c4 + rhototal) - rho53)/(c4 + rhototal)**2
          if (lgrad2) then
            rho13 = rhototal**(-1.0_dp/3.0_dp)
            rscrho3 = - (2.0_dp*c1 + c3*((10.0_dp/9.0_dp)*rho13*(c4 + rhototal)**2 - &
                      (10.0_dp/3.0_dp)*rho23*(c4 + rhototal) + 2.0_dp*rho53)/(c4 + rhototal)**3)
            if (lgrad3) then
              rho43 = rhototal**(-1.0_dp/3.0_dp)
              rscrho5 = - (c3*(-(10.0_dp/27.0_dp)*rho43*(c4 + rhototal)**3 - (10.0_dp/3.0_dp)*rho13*(c4 + rhototal)**2 + &
                        10.0_dp*rho23*(c4 + rhototal) - 6.0_dp*rho53)/(c4 + rhototal)**4)
            endif
          endif
        endif
      elseif (neamfn.eq.8) then
!
!  Mei-Davenport - note sign change wrapped into f0 to keep consistency with paper
!
        f0 = - eamfnpar(1,m)
        alpha = eamfnpar(2,m)
        beta = eamfnpar(3,m)
        gamma = eamfnpar(4,m)
        delta = eamfnpar(5,m)
        phi0 = eamfnpar(6,m)
!
        rnpx = alpha/beta
        rhorrn3x = rhototal**(rnpx - 3.0_dp)
        rhorrn2x = rhorrn3x*rhototal
        rhorrn1x = rhorrn2x*rhototal
        rhorrnx = rhorrn1x*rhototal
        rlogtrm = log(rhototal)
!
        eeam = f0*(1.0_dp - rnpx*rlogtrm)*rhorrnx
        if (lgrad1) then
          rscrho = - f0*rnpx*rnpx*rlogtrm*rhorrn1x
          if (lgrad2) then
            rscrho3 = - f0*rnpx*rnpx*rhorrn2x*((rnpx - 1.0_dp)*rlogtrm + 1.0_dp)
            if (lgrad3) then
              rscrho5 = - f0*rnpx*rnpx*rhorrn3x*((rnpx - 1.0_dp)*(rnpx - 2.0_dp)*rlogtrm + 2.0_dp*rnpx - 3.0_dp)
            endif
          endif
        endif
        do k = 1,3
          rk = dble(k)
          rootk = sqrt(rk)
          exptrm = exp(-(rootk - 1.0_dp)*gamma)
          f0 = - 0.5_dp*phi0*eamfnpar(6+k,m)*exptrm
          atrm = 1.0_dp + (rootk - 1.0_dp)*delta
          btrm = rootk*delta/beta
!
          rnpx = rootk*gamma/beta
          rhorrn3x = rhototal**(rnpx - 3.0_dp)
          rhorrn2x = rhorrn3x*rhototal
          rhorrn1x = rhorrn2x*rhototal
          rhorrnx = rhorrn1x*rhototal
!
          eeam = eeam - f0*(atrm - btrm*rlogtrm)*rhorrnx
          if (lgrad1) then
            trm1 = (rnpx*atrm - rnpx*btrm*rlogtrm - btrm)
            rscrho = rscrho - f0*rhorrn1x*trm1
            if (lgrad2) then
              rscrho3 = rscrho3 - f0*rhorrn2x*(rnpx*(rnpx - 1.0_dp)*(atrm - btrm*rlogtrm) - btrm*(2.0_dp*rnpx - 1.0_dp))
              if (lgrad3) then
                rscrho5 = rscrho5 - f0*rhorrn3x*(rnpx*(rnpx - 1.0_dp)*(rnpx - 2.0_dp)*(atrm - btrm*rlogtrm) - &
                                    btrm*(3.0_dp*rnpx*rnpx - 6.0_dp*rnpx + 2.0_dp))
              endif
            endif
          endif
        enddo
      elseif (neamfn.eq.9) then
!
!  Baskes
!
        AE0 = f0*eamfnpar(2,m)
        rhoe = 1.0_dp/eamfnpar(3,m)
        rhorr = rhototal*rhoe
        logrho = log(rhorr)
        eeam = AE0*rhorr*logrho
        if (lgrad1) then
          rscrho = AE0*(logrho + 1.0_dp)*rhoe
          if (lgrad2) then
            rscrho3 = AE0*rhoe/rhototal
            if (lgrad3) then
              rscrho5 = rscrho3/rhototal
            endif
          endif
        endif
      elseif (neamfn.eq.12) then
!
!  Igarashi
!
        f1 = eamfnpar(2,m)
        rhoe = abs(rhototal)
        rootrho = sqrt(rhoe*(1.0_dp + f1*rhoe))
        eeam = - f0*rootrho
        if (lgrad1) then
          d1rootrho = 1.0_dp + 2.0_dp*f1*rhoe
          rscrho = - 0.5_dp*f0*d1rootrho/rootrho
          if (lgrad2) then
            d2rootrho = 2.0_dp*f1
            rscrho3 = 0.25_dp*d1rootrho*d1rootrho/rootrho**3 - 0.5_dp*d2rootrho/rootrho
            if (lgrad3) then
              rscrho5 = f0*(-0.375_dp*(d1rootrho**3)/rootrho**5 + 0.75_dp*d2rootrho*d1rootrho/rootrho**3)
            endif
            rscrho3 = f0*rscrho3
          endif
        endif
      elseif (neamfn.eq.13) then
!
!  Spline
!
        if (scrho(1).ge.eamfnpar(5,m).and.scrho(1).lt.eamfnpar(6,m)) then
          rhon = (scrho(1) - eamfnpar(5,m))
          eeam = - (eamfnpar(1,m)*rhon**3 + eamfnpar(2,m)*rhon**2 + eamfnpar(3,m)*rhon + eamfnpar(4,m))
          if (lgrad1) then
            rscrho = - (3.0_dp*eamfnpar(1,m)*rhon**2 + 2.0_dp*eamfnpar(2,m)*rhon + eamfnpar(3,m))
            if (lgrad2) then
              rscrho3 = - (6.0_dp*eamfnpar(1,m)*rhon + 2.0_dp*eamfnpar(2,m))
              if (lgrad3) then
                rscrho5 = - 6.0_dp*eamfnpar(1,m)
              endif
            endif
          endif
        endif
      endif
!
!  Correct for alloy parameter
!
      if (lgrad1) then
        rscrho = rscrho/eamalloy(1,m)
        if (lgrad2) then
          rscrho3 = rscrho3/eamalloy(1,m)
          if (lgrad3) then
            rscrho5 = rscrho5/eamalloy(1,m)
          endif
        endif
      endif
    else
      eeam = 0.0_dp
      if (lgrad1) then
        rscrho = 0.0_dp
        if (lgrad2) then
          rscrho3 = 0.0_dp
          if (lgrad3) then
            rscrho5 = 0.0_dp
          endif
        endif
      endif
    endif
  else
    eeam = 0.0_dp
    if (lgrad1) then
      rscrho = 0.0_dp
      if (lgrad2) then
        rscrho3 = 0.0_dp
        if (lgrad3) then
          rscrho5 = 0.0_dp
        endif
      endif
    endif
  endif
#ifdef TRACE
  call trace_out('meanfnderv')
#endif
!
  return
  end
