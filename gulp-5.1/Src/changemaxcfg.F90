  subroutine changemaxcfg
!
!  Alters the size of the arrays associated with maxcfg
!
!   9/03 lregionrigid flag added
!   1/05 ndistancereset added
!  11/05 hmssg initialised as (unknown)
!   2/06 omegadirtype added
!   3/06 omega damping factor added
!  11/06 NEB modifications added
!  11/06 lfcborn, lfcphon and lfcprop initialised to false
!   2/07 Electric field arrays added
!   3/07 Radial force added
!   3/07 Initialisation of further variables added
!   5/07 nregiontype added
!   5/07 QMMMmode added
!   7/07 Metadynamics arrays added
!   3/08 New MD integrator arrays added
!   4/08 freaction added
!   8/08 Metadynamics arrays modified
!   8/08 pr_conscfg added
!  10/08 COSMO modifications added
!   6/09 PDF changes added
!   9/10 lfcscatter added
!   9/10 Initialisations now performed in a subroutine
!  11/10 Anisotropic pressure added
!   9/11 Metadynamics internal code replaced with Plumed
!   6/12 nobsmodecfg added
!  10/12 Arrays added for OpenKIM models
!  12/12 Time-dependent field added
!  12/12 Modified to allow for multiple time-dependent fields
!   3/13 nreg1cfg, nvacacfg and nintecfg  added
!   3/13 lindvacptr added
!   7/13 Symmetry number array added
!   8/13 Thermal conductivity arrays added
!   9/13 nsuper changed to be a 2-D array
!   9/13 Raman directions added
!   4/14 lbornkin added
!  10/14 nd2cellcfg added
!   4/15 Ghost supercell array added
!   8/15 New lower options added
!   9/15 translate noise added
!   5/16 SPME arrays added
!   8/16 OpenKIM changes added
!   9/16 hmssgio added
!   3/17 fix_atom option added
!   4/17 qonsas option added
!  10/17 Initial coordinates added to restart info
!   1/18 Reallocation for pkim_model removed for F03 version
!   2/18 ltrantherm added
!   2/18 nxks, nyks, nzks converted to a single array
!   2/18 Arrays for handling primitive cell added
!   3/18 Sign option added to translate
!   4/18 Arrays for scan_cell added
!   4/18 Twist option added
!   5/18 lcscanstrain added
!   6/18 straincfg added
!  11/18 Modified for version 2 of OpenKIM
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
  use configurations
  use cosmic
  use current,        only : maxat
  use defects
  use distances,      only : ndistancereset
  use field
  use freeze
  use genetic,        only : xmaxcfg, xmincfg
  use kim_models
  use ksample
  use moldyn
  use m_pr
  use m_pdfneutron,   only : changemaxpdfcfg
  use g_neb
  use observables,    only : freaction, maxobs, nobsmodecfg
  use potentialgrid
  use potentialpoints
  use projectdos
  use radial
  use reallocate
  use scan
  use shifts
  use spme,           only : nqkgrid
  use symmetry
  use thermalcond
  implicit none
!
!  Local variables
!
  integer(i4)       :: ierror, i
  integer(i4), save :: oldmaxcfg = 0
!
!  Configuration data
!
  call realloc(anisotropicpresscfg,6_i4,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','anisotropicpresscfg')
  call realloc(lbornkin,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','lbornkin')
  call realloc(bornk,3_i4,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','bornk')
  call realloc(cosmoeigen,3_i4,3_i4,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','cosmoeigen')
  call realloc(dhklcfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','dhklcfg')
  call realloc(lcosmoeigin,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','lcosmoeigin')
  call realloc(cosmoepsilon,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','cosmoepsilon')
  call realloc(cosmodrsolv,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','cosmodrsolv')
  call realloc(cosmorsolv,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','cosmorsolv')
  call realloc(energycfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','energycfg')
  call realloc(freaction,maxcfg,maxobs,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','freaction')
  call realloc_ch1(hmssg,16_i4,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','hmssg')
  call realloc_ch1(hmssgio,16_i4,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','hmssgio')
  call realloc(ifhr,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','ifhr')
  call realloc(iflags,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','iflags')
  call realloc(ifso,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','ifso')
  call realloc(iperm,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','iperm')
  call realloc(iufree,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','iufree')
  call realloc(ivso,3_i4,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','ivso')
  call realloc(lanisotropicpresscfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','lanisotropicpresscfg')
  call realloc(lcscanstrain,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','lcscanstrain')
  call realloc(ldeflin,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','ldeflin')
  call realloc(lfcborn,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','lfcborn')
  call realloc(lfcphon,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','lfcphon')
  call realloc(lfcprop,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','lfcprop')
  call realloc(lfcscatter,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','lfcscatter')
  call realloc(lindintptr,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','lindintptr')
  call realloc(lindvacptr,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','lindvacptr')
  call realloc(linitcfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','linitcfg')
  call realloc(llowered,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','llowered')
  call realloc(lmdconstrain,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','lmdconstrain')
  call realloc(lomega,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','lomega')
  call realloc(lopfc,6_i4*maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','lopfc')
  call realloc(lopfreg,3_i4*maxregion,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','lopfreg')
  call realloc(lregionrigid,maxregion,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','lregionrigid')
  call realloc(lreldin,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','lreldin')
  call realloc(lsymset,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','lsymset')
  call realloc(ltrannoise,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','ltrannoise')
  call realloc(ltrantherm,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','ltrantherm')
  call realloc(lufree,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','lufree')
  call realloc(lvecin,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','lvecin')
  call realloc(maxmodecfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','maxmodecfg')
  call realloc(minmodecfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','minmodecfg')
  call realloc(ncellmaxcfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','ncellmaxcfg')
  call realloc(ncellmincfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','ncellmincfg')
  call realloc(ncorepcfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','ncorepcfg')
  call realloc(ngocfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','ngocfg')
  call realloc(nqkgrid,3_i4,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nqkgrid')
  call realloc(lnebvaryspring,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','lnebvaryspring')
  call realloc(nebspring,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nebspring')
  call realloc(nebspringmin,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nebspringmin')
  call realloc(nebfinalcell,6_i4,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nebfinalcell')
  call realloc(nebfinalradius,maxat,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nebfinalradius')
  call realloc(nebfinalxyz,3_i4,maxat,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nebfinalxyz')
  call realloc(nfixatomtypecfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nfixatomtypecfg')
  call realloc(nfixatomcfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nfixatomcfg')
  call realloc(ngocfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','ngocfg')
  call realloc(ninternalmaxcfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','ninternalmaxcfg')
  call realloc(ninternalmincfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','ninternalmincfg')
  call realloc(nnebreplica,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nnebreplica')
  call realloc(nreg1cfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nreg1cfg')
  call realloc(nintecfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nintecfg')
  call realloc(nvacacfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nvacacfg')
  call realloc(nregions,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nregions')
  call realloc(nummodecfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nummodecfg')
  call realloc(n1con,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','n1con')
  call realloc(n1var,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','n1var')
  call realloc_ch80(names,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','names')
  call realloc(nascfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nascfg')
  call realloc(nbornstep,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nbornstep')
  call realloc(ndcentyp,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','ndcentyp')
  call realloc(nd2cellcfg,3_i4,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nd2cellcfg')
  call realloc(ndimen,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','ndimen')
  call realloc(ndistancereset,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','ndistancereset')
  call realloc(neiglow,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','neiglow')
  call realloc(neighigh,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','neighigh')
  call realloc(nensemble,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nensemble')
  call realloc(nmdconstrainatom,2_i4,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nmdconstrainatom')
  call realloc(nmdconstraindist,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nmdconstraindist')
  call realloc(nmdeq,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nmdeq')
  call realloc(nmdprod,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nmdprod')
  call realloc(nmdsamp,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nmdsamp')
  call realloc(nmdvelmode,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nmdvelmode')
  call realloc(nmdvelmodp,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nmdvelmodp')
  call realloc(nmdwrite,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nmdwrite')
  call realloc(nobsmodecfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nobsmodecfg')
  call realloc(nomegastep,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nomegastep')
  call realloc(norigkpt,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','norigkpt')
  call realloc(npotptcfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','npotptcfg')
  call realloc(nprojcfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nprojcfg')
  call realloc(nprojdef,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nprojdef')
  call realloc(nsasexcludemax,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nsasexcludemax')
  call realloc(nsasexcludemin,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nsasexcludemin')
  call realloc(nshcfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nshcfg')
  call realloc(nccscfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nccscfg')
  call realloc(nregiontype,maxregion,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nregiontype')
  call realloc(nspcg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nspcg')
  call realloc(nsregion2,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nsregion2')
  call realloc(nsuper,3_i4,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nsuper')
  call realloc(nsuperghost,3_i4,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nsuperghost')
  call realloc(ntempstp,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','ntempstp')
  call realloc(ntempstpstart,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','ntempstpstart')
  call realloc(ncscan,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','ncscan')
  call realloc(ntran,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','ntran')
  call realloc(ntwistcfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','ntwistcfg')
  call realloc(nvarcfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nvarcfg')
  call realloc(nzmolcfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nzmolcfg')
  call realloc(nks,3_i4,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nks')
  call realloc(nksala,3_i4,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nksala')
  call realloc(nxpg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nxpg')
  call realloc(nypg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nypg')
  call realloc(nzpg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','nzpg')
  call realloc(QMMMmode,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','QMMMmode')
  call realloc(omega,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','omega')
  call realloc(omegadamping,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','omegadamping')
  call realloc(omegadir,6_i4,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','omegadir')
  call realloc(omegadirtype,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','omegadirtype')
  call realloc(omegastep,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','omegastep')
  call realloc(presscfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','presscfg')
  call realloc(qonsastarget,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','qonsastarget')
  call realloc(ramandir,6_i4,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','ramandir')
  call realloc(ramandirtype,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','ramandirtype')
  call realloc(qpres,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','qpres')
  call realloc(qtemp,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','qtemp')
  call realloc(reg1,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','reg1')
  call realloc(reg1last,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','reg1last')
  call realloc(reg2,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','reg2')
  call realloc(reg2a1,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','reg2a1')
  call realloc(ropcfg,3_i4,3_i4,maxsymop,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','ropcfg')
  call realloc(rufree,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','rufree')
  call realloc(rvcfg,3_i4,3_i4,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','rvcfg')
  call realloc(rvpcfg,3_i4,3_i4,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','rvpcfg')
  call realloc(shift,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','shift')
  call realloc(shscalecfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','shscalecfg')
  call realloc(sbulkecfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','sbulkecfg')
  call realloc(straincfg,6_i4,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','straincfg')
  call realloc(stresscfg,6_i4,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','stresscfg')
  call realloc(symnocfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','symnocfg')
  call realloc(taubcfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','taubcfg')
  call realloc(tautcfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','tautcfg')
  call realloc(pr_conscfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','pr_conscfg')
  call realloc(tempcfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','tempcfg')
  call realloc(tempstp,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','tempstp')
  call realloc(tmdeq,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','tmdeq')
  call realloc(tmdforcestart,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','tmdforcestart')
  call realloc(tmdforcestop,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','tmdforcestop')
  call realloc(tmdprod,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','tmdprod')
  call realloc(tmdsamp,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','tmdsamp')
  call realloc(tmdscale,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','tmdscale')
  call realloc(tmdscint,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','tmdscint')
  call realloc(tmdwrite,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','tmdwrite')
  call realloc(trannoise,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','trannoise')
  call realloc(trantherm,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','trantherm')
  call realloc(tstep,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','tstep')
  call realloc(totalchargecfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','totalchargecfg')
  call realloc(vitcfg,3_i4,maxsymop,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','vitcfg')
  call realloc(xdcent,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','xdcent')
  call realloc(ydcent,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','ydcent')
  call realloc(zdcent,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','zdcent')
  call realloc(xmaxpg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','xmaxpg')
  call realloc(ymaxpg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','ymaxpg')
  call realloc(zmaxpg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','zmaxpg')
  call realloc(xminpg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','xminpg')
  call realloc(yminpg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','yminpg')
  call realloc(zminpg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','zminpg')
  call realloc(cscan,6_i4,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','cscan')
  call realloc(xtran,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','xtran')
  call realloc(ytran,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','ytran')
  call realloc(ztran,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','ztran')
  call realloc(xmaxcfg,3_i4,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','xmaxcfg')
  call realloc(xmincfg,3_i4,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','xmincfg')
  call realloc(xufree,3_i4,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','xufree')
!
  call realloc(lfieldcfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','lfieldcfg')
  call realloc(ntdfieldcfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','ltdfieldcfg')
  call realloc(fieldcfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','fieldcfg')
  call realloc(fielddirectioncfg,3_i4,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','fielddirectioncfg')
  call realloc(td_fieldcfg,3_i4,maxtdfield,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','td_fieldcfg')
  call realloc(td_fielddirectioncfg,3_i4,maxtdfield,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','td_fielddirectioncfg')
  call realloc(tmdfieldstart,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','tmdfieldstart')
  call realloc(tmdfieldstop,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','tmdfieldstop')
!
  call realloc(lradialcfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','lradialcfg')
  call realloc(radialKcfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','radialKcfg')
  call realloc(radialXYZcfg,3_i4,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','radialXYZcfg')
!
!  OpenKIM arrays
!
  call realloc(lkim_model_cfg_OK,maxcfg,maxkimmodel,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','lkim_model_cfg_OK')
!
  call realloc(lomega_af_in,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','lomega_af_in')
  call realloc(omega_af,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','omega_af')
  call realloc(B_pr_cfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','B_pr_cfg')
  call realloc(v_s_cfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','v_s_cfg')
  call realloc(v_p_cfg,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcfg','v_p_cfg')
!
!  Reallocate PDF related variables
!
  call changemaxpdfcfg
!
!  Initialise defaults for new part of array
!
  if (maxcfg.gt.oldmaxcfg) then
    do i = oldmaxcfg+1,maxcfg
      call initmaxcfgdefaults(i)
    enddo
  endif
!
!  Save current value of maxcfg for next call
!
  oldmaxcfg = maxcfg
!
  return
  end
