  subroutine outener
!
!  Output energy components
!
!   6/95 Dipole energy added
!   3/97 Output in KJmol-1 added
!  12/97 Self energy term added to total for EEM/QEq
!   3/01 Output of surface energy added
!  11/01 Attachment energy added
!   5/02 Scaling of shift added
!   5/02 Brenner energy added
!   8/02 External force modifications added
!   9/02 Polymer energy added
!   9/02 Region 1-2 interaction energy output
!  10/02 Region 2 energy output
!  11/02 Einstein energy output
!   1/03 Wolf self energy output
!   6/03 XML modifications added
!   7/03 Bond order potential energy added
!  11/04 Six body energy added
!   7/05 Labelling of QEq terms changed to SM when needed
!   7/05 Charge coupled potential energy added
!   7/06 Six body energy added
!   3/07 Electric field energy added
!   3/07 Radial force added
!   3/07 Chemshell energy added
!   7/07 ReaxFF energy added
!   7/07 Plane potential energy added
!  12/07 Name of Coulomb correction changed to reaxFF for this 
!        case to reflect option specified
!  10/08 COSMO energy added
!   2/09 Modified to accommodate new version of FoX and gulp_cml
!   6/09 Module name changed from three to m_three
!   1/10 One-body potentials added
!   9/10 EDIP added
!  11/10 Energy from anisotropic pressure added
!   1/11 Energy from anisotropic pressure removed
!   9/11 Madelung correction added
!  11/11 Region breakdown of energies added
!  10/12 OpenKIM energy added
!   7/13 Improper torsion energy added
!   3/14 Harmonic relaxation energy added
!   3/15 kcal keyword added
!   3/16 eplumed added
!   4/17 ChemShell interaction modified
!   6/17 Old ChemShell restored as a compile option
!  11/18 Output of attachment energies disabled if using KIM models since
!        values will not be correct
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
!  Julian Gale, CIC, Curtin University, November 2018
!
  use bondorderdata,  only : nbopot, nboQ
  use cellmultipole
  use chargecoupled,  only : nCCspec
  use configurations
  use control
  use current
  use element,        only : lSandM
  use energies
  use four
  use g_constants,    only : kjmtoev, kcaltoev
  use gulp_cml,       only : gulp_cml_outener, lcml
  use iochannels
  use kim_models,     only : lkim_model
  use m_plumed,       only : lplumed
  use m_three
  use one,            only : none
  use plane,          only : nplanepot
  use shifts
  use six
  use sutton
  use symmetry
  implicit none
!
!  Local variables
!
  character(len=10) :: cmmword(4)
  character(len=6)  :: unitword
  integer(i4)       :: i
  integer(i4)       :: icf
  integer(i4)       :: j
  logical           :: lkcal
  logical           :: lkjmol
  logical           :: lpress
  real(dp)          :: etot
  real(dp)          :: evtokjmol
  real(dp)          :: sft
  real(dp)          :: units
!
  data cmmword/'monopole  ','dipole    ','quadrupole','octopole  '/
!
  if (lcml) then
    call gulp_cml_outener 
  endif
!
!  Set units for output
!
  lkcal  = (index(keyword,'kcal').eq.1.or.index(keyword,' kcal').ne.0)
  lkjmol = (index(keyword,'kjmo').eq.1.or.index(keyword,' kjmo').ne.0)
!
  if (lkcal) then
    units = 1.0_dp/kcaltoev
    unitword = 'kcal  '
  elseif (lkjmol) then
    units = 1.0_dp/kjmtoev
    unitword = 'kJ/mol'
  else
    units = 1.0_dp
    unitword = 'eV    '
  endif
!
  evtokjmol = 23.0604_dp*4.184_dp
  lpress = ((abs(press).gt.0.0_dp.and.ndim.eq.3).or.lanisotropicpress)
  sft = shift(nshcfg(ncf))*shscalecfg(ncf)
  etot = erecip + eatom + sft + ethb + ereal + efor + epv + ecmm + ec6 + eoop + emany + &
         edipole + ebgd + eself + eqeq + evib + epolar + ebrenner + eforce + esix + &
         esregion12 + esregion2 + eeinstein + ewolfself + ebondorder + eboQself +  &
         echargecoupled + efield + eradial + ereaxFF + eplane + ecosmo + eone + &
         eedip + emad + eimp + ekim + erelax + eplumed
!
  if (lfree) then
    write(ioout,'(/,''  Components of free energy : '',/)')
  elseif (lpress) then
    write(ioout,'(/,''  Components of enthalpy : '',/)')
  else
    if (icmm.gt.0) then
      write(ioout,'(/,''  Components of energy : (CMM to '',a10,'')'',/)') cmmword(icmm)
    else
      write(ioout,'(/,''  Components of energy : '',/)')
    endif
  endif
  write(ioout,'(''--------------------------------------------------------------------------------'')')
  write(ioout,'(''  Interatomic potentials     = '',f20.8,1x,a6)') eatom*units,unitword
  if (none.gt.0) then
    write(ioout,'(''  One-body potentials        = '',f20.8,1x,a6)') eone*units,unitword
  endif
  if (nthb.gt.0) then
    write(ioout,'(''  Three-body potentials      = '',f20.8,1x,a6)') ethb*units,unitword
  endif
  if (nfor.gt.0) then
    write(ioout,'(''  Four-body potentials       = '',f20.8,1x,a6)') efor*units,unitword
    write(ioout,'(''  Improper torsions          = '',f20.8,1x,a6)') eimp*units,unitword
    write(ioout,'(''  Out of plane potentials    = '',f20.8,1x,a6)') eoop*units,unitword
  endif
  if (nsix.gt.0) then
    write(ioout,'(''  Six-body potentials        = '',f20.8,1x,a6)') esix*units,unitword
  endif
  if (lsuttonc) then
    write(ioout,'(''  Many-body potentials       = '',f20.8,1x,a6)') emany*units,unitword
  endif
  if (lbrenner) then
    write(ioout,'(''  Brenner potentials         = '',f20.8,1x,a6)') ebrenner*units,unitword
  endif
  if (nbopot.gt.0) then
    write(ioout,'(''  Bond-order potentials      = '',f20.8,1x,a6)') ebondorder*units,unitword
  endif
  if (lreaxFF) then
    write(ioout,'(''  ReaxFF force field         = '',f20.8,1x,a6)') ereaxFF*units,unitword
  endif
  if (lEDIP) then
    write(ioout,'(''  EDIP force field           = '',f20.8,1x,a6)') eedip*units,unitword
  endif
  if (nboQ.gt.0) then
    write(ioout,'(''  Bond-order self energy     = '',f20.8,1x,a6)') eboQself*units,unitword
  endif
  if (nCCspec.gt.0) then
    write(ioout,'(''  Charge-coupled potentials  = '',f20.8,1x,a6)') echargecoupled*units,unitword
  endif
  if (lkim_model) then
    write(ioout,'(''  OpenKIM potentials         = '',f20.8,1x,a6)') ekim*units,unitword
  endif
  if (lplumed) then
    write(ioout,'(''  Plumed                     = '',f20.8,1x,a6)') eplumed*units,unitword
  endif
  write(ioout,'(''  Monopole - monopole (real) = '',f20.8,1x,a6)') ereal*units,unitword
  if (icmm.gt.0) then
    write(ioout,'(''  Monopole - multipole       = '',f20.8,1x,a6)') ecmm*units,unitword
  endif
  if (ndim.gt.0) then
    write(ioout,'(''  Monopole - monopole (recip)= '',f20.8,1x,a6)') erecip*units,unitword
    write(ioout,'(''  Monopole - monopole (total)= '',f20.8,1x,a6)') (erecip + ereal)*units,unitword
  endif
  if (lc6.and.ndim.eq.3) then
    write(ioout,'(''  Dispersion (real+recip)    = '',f20.8,1x,a6)') ec6*units,unitword
  endif
  if (abs(epolar).gt.1d-8) then
    write(ioout,'(''  Polarisation energy        = '',f20.8,1x,a6)') epolar*units,unitword
  endif
  if (abs(erelax).gt.1d-8) then
    write(ioout,'(''  Harmonic relaxation energy = '',f20.8,1x,a6)') erelax*units,unitword
  endif
  if (abs(ecosmo).gt.1d-8) then
    write(ioout,'(''  Solvation energy           = '',f20.8,1x,a6)') ecosmo*units,unitword
  endif
  if (leinstein) then
    write(ioout,'(''  Einstein energy            = '',f20.8,1x,a6)') eeinstein*units,unitword
  endif
  if (nplanepot.gt.0) then
    write(ioout,'(''  Plane potential energy     = '',f20.8,1x,a6)') eplane*units,unitword
  endif
  if (abs(sft).gt.1d-8) then
    write(ioout,'(''  Energy shift               = '',f20.8,1x,a6)') sft*units,unitword
  endif
  if (edipole.gt.1.0d-8) then
    write(ioout,'(''  Dipole correction energy   = '',f20.8,1x,a6)') edipole*units,unitword
  endif
  if (abs(ebgd).gt.1.0d-8) then
    write(ioout,'(''  Neutralising energy        = '',f20.8,1x,a6)') ebgd*units,unitword
  endif
  if (abs(emad).gt.1.0d-8) then
    write(ioout,'(''  Madelung neutrality energy = '',f20.8,1x,a6)') emad*units,unitword
  endif
  if (abs(ewolfself).gt.1.0d-8) then
    write(ioout,'(''  Wolf sum self energy       = '',f20.8,1x,a6)') ewolfself*units,unitword
  endif
  if (abs(eself).gt.1.0d-8) then
    write(ioout,'(''  Self energy (EEM/QEq/SM)   = '',f20.8,1x,a6)') eself*units,unitword
  endif
  if (abs(eqeq).gt.1.0d-8) then
    if (lSandM) then
      write(ioout,'(''  SM Coulomb correction      = '',f20.8,1x,a6)') eqeq*units,unitword
    else
      write(ioout,'(''  QEq Coulomb correction     = '',f20.8,1x,a6)') eqeq*units,unitword
    endif
  endif
  if (abs(esregion12).gt.1.0d-8) then
    write(ioout,'(''  Region 1-2 interaction     = '',f20.8,1x,a6)') esregion12*units,unitword
  endif
  if (abs(esregion2).gt.1.0d-8) then
    write(ioout,'(''  Region 2-2 interaction     = '',f20.8,1x,a6)') esregion2*units,unitword
  endif
  if (lpress) then
    write(ioout,'(''  Pressure*volume            = '',f20.8,1x,a6)') epv*units,unitword
  endif
  if (abs(eradial).gt.1.0d-8) then
    write(ioout,'(''  Radial                     = '',f20.8,1x,a6)') eradial*units,unitword
  endif
  if (abs(eforce).gt.1.0d-8) then
    write(ioout,'(''  External_force*distance    = '',f20.8,1x,a6)') eforce*units,unitword
  endif
  if (abs(efield).gt.1.0d-8) then
    write(ioout,'(''  Electric_field*distance    = '',f20.8,1x,a6)') efield*units,unitword
  endif
  if (abs(echemsh).gt.1.0d-8) then
    write(ioout,'(''  Chemshell contribution     = '',f20.8,1x,a6)') echemsh*units,unitword
  endif
  if (abs(evib).gt.1.0d-8) then
    write(ioout,'(''  Vibrational contribution   = '',f20.8,1x,a6)') evib*units,unitword
  endif
  write(ioout,'(''--------------------------------------------------------------------------------'')')
  icf = icentfct(ncbl)
  if (lsymopt.and.icf.gt.1) then
    if (lfree) then
      write(ioout,'(''  Total free energy : '')')
    elseif (lpress) then
      write(ioout,'(''  Total lattice enthalpy : '')')
    else
      write(ioout,'(''  Total lattice energy : '')')
    endif
    write(ioout,'(''    Primitive unit cell      = '',f20.8,1x,a6)') etot*units,unitword
    write(ioout,'(''    Non-primitive unit cell  = '',f20.8,1x,a6)') etot*dble(icf)*units,unitword
  else
    if (lfree) then
      write(ioout,'(''  Total free energy          = '',f20.8,1x,a6)') etot*units,unitword
    elseif (lpress) then
      write(ioout,'(''  Total lattice enthalpy     = '',f20.8,1x,a6)') etot*units,unitword
    else
      write(ioout,'(''  Total lattice energy       = '',f20.8,1x,a6)') etot*units,unitword
    endif
  endif
  write(ioout,'(''--------------------------------------------------------------------------------'')')
  if (lsymopt.and.icf.gt.1) then
    if (lfree) then
      write(ioout,'(''  Total free energy (in kJmol-1): '')')
    elseif (lpress) then
      write(ioout,'(''  Total lattice enthalpy (in kJmol-1): '')')
    else
      write(ioout,'(''  Total lattice energy (in kJmol-1): '')')
    endif
    write(ioout,'(''    Primitive unit cell      = '',f20.4,'' kJ/(mole unit cells)'')') etot*evtokjmol
    write(ioout,'(''    Non-primitive unit cell  = '',f20.4,'' kJ/(mole unit cells)'')') etot*dble(icf)*evtokjmol
  else
    if (lfree) then
      if (ndim.gt.0) then
        write(ioout,'(''  Total free energy          = '',f20.4,'' kJ/(mole unit cells)'')') etot*evtokjmol
      else
        write(ioout,'(''  Total free energy          = '',f20.4,'' kJ/mol'')') etot*evtokjmol
      endif
    elseif (lpress) then
      write(ioout,'(''  Total lattice enthalpy     = '',f20.4,'' kJ/(mole unit cells)'')') etot*evtokjmol
    else
      if (ndim.gt.0) then
        write(ioout,'(''  Total lattice energy       = '',f20.4,'' kJ/(mole unit cells)'')') etot*evtokjmol
      else
        write(ioout,'(''  Total lattice energy       = '',f20.4,'' kJ/mol'')') etot*evtokjmol
      endif
    endif
  endif
  write(ioout,'(''--------------------------------------------------------------------------------'')')
!
!  Surface energy calculation
!
  if (lseok) then
    if (ndim.eq.2) then
      write(ioout,'(''  Bulk energy                = '',f20.6,1x,a6)') sbulkecfg(ncf)*units,unitword
      if (nregions(ncf).gt.1) then
        write(ioout,'(''  Surface energy (region 1)  = '',f20.6,'' J/m**2'')') esurface
      else
        write(ioout,'(''  Surface energy (full slab) = '',f20.6,'' J/m**2'')') esurface
        write(ioout,'(''  Surface energy (per  side) = '',f20.6,'' J/m**2'')') esurface*0.5_dp
      endif
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    else
      if (nregions(ncf).gt.1) then
        write(ioout,'(''  Polymer energy (region 1)  = '',f20.6,'' eV/Ang'')') esurface
        write(ioout,'(''--------------------------------------------------------------------------------'')')
      endif
    endif
  endif
!
!  Attachment energy
!
  if (ndim.eq.2.and.abs(eattach).gt.1.0d-8) then
    write(ioout,'(''  Attachment energy          = '',f20.8,1x,a6)') eattach*units,unitword
    write(ioout,'(''  Attachment energy/unit     = '',f20.8,1x,a6)') (eattach/dble(nzmol))*units,unitword
    write(ioout,'(''--------------------------------------------------------------------------------'')')
  endif
!
!  Attachment energy - not output if using OpenKIM since eattach will not be correct
!
  if (.not.lkim_model) then
    if (ndim.eq.2.and.abs(eattach).gt.1.0d-8) then
      write(ioout,'(''  Attachment energy          = '',f20.8,1x,a6)') eattach*units,unitword
      write(ioout,'(''  Attachment energy/unit     = '',f20.8,1x,a6)') (eattach/dble(nzmol))*units,unitword
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    endif
  endif
!
!  Region - region energies
!
  if (leregion) then
    write(ioout,'(/,''  Region - region interaction energies (eV):'',/)') 
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    write(ioout,'(8x,4(8x,i4,4x))') (i,i=1,min(nregions(ncf),4_i4))
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    do i = 1,nregions(ncf)
      write(ioout,'(2x,i4,2x,4(f16.4))') i,(eregion2region(j,i),j=1,min(i,4_i4))
    enddo
    write(ioout,'(''--------------------------------------------------------------------------------'')')
  endif
!
  return
  end
