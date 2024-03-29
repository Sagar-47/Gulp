  subroutine nullpointer
!
!  This subroutine nullifies all pointers in the uses so that
!  they may be subsequently passed to realloc without problems
!  even if the pointer hasn't been previously assigned.
!
!   9/06 Arrays added for literal symbols
!  10/06 Neutron arrays nullified
!  11/06 NEB arrays added
!  11/06 x/y/zfracimage arrays added
!   2/07 nbondedtype added
!   2/07 Electric field arrays added
!   3/07 Radial force added
!   3/07 Chemshell changes added
!   4/07 UFF data arrays added
!   5/07 xfsave/yfsave/zfsave added
!   5/07 UFFtor added
!   5/07 nregiontype added
!   5/07 Bond increment charge arrays added
!   5/07 Partial occupancy arrays added
!   5/07 ltdreiding added
!   7/07 symbolUFFspec added
!   7/07 GCMC molecule flag added
!   7/07 Metadynamics variables added
!   7/07 ReaxFF data arrays added
!   7/07 Plane potential added
!  12/07 Extra ReaxFF array added
!   1/08 ltrialatom added
!   3/08 numofspec added
!   4/08 Extra ReaxFF arrays added
!   4/08 nobptr3 added
!   4/08 freaction added
!   5/08 New arrays for out of plane term added
!   8/08 pr_conscfg added
!  10/08 COSMO arrays added
!  10/08 MEAM modifications added
!  11/08 Logical array to flag out of plane potentials added
!  11/08 lmeamspec added
!  12/08 Module input renamed to gulpinput
!  12/08 eammeamcoeff array changed from density to function related
!   1/09 swap move added to Monte Carlo
!   1/09 Arrays apot, bpot, cpot & dpot replaced by new 2-D array twopot
!   1/09 leshift/lgshift logicals introduced to replace use of tpot(3/4,) for energy
!        and gradient shifts
!   1/09 Core-shell vector array added
!   1/09 lreaxFFpboOK added 
!   3/09 lkptdispersion added
!   6/09 Module name changed from three to m_three
!   6/09 EVB arrays added
!   6/09 New molecule indexing arrays added
!   6/09 Neutron related variables from ers29 removed
!   7/09 EAM symbol arrays added
!   1/10 One-body potentials added
!   3/10 Arrays that flag rule generated potentials added
!   6/10 Modified to handle optional second atom type for bond-order
!        attractive and repulsive terms.
!   9/10 Neutron scattering modifications added
!   9/10 EDIP arrays added
!   1/11 Reference cell vectors added
!   8/11 nreaxFFfixQspecptr added
!  11/11 eregion2region added
!   4/12 xvir, yvir and zvir removed
!   5/12 Atomic stress array added
!   6/12 Eigenvector input for mode observables added
!   7/12 nphonatptr and nphonatrptr added
!  10/12 Support for OpenKIM models added
!  12/12 Time-dependent field added
!  12/12 Modified to allow for multiple time-dependent fields
!   3/12 New defect arrays added
!   7/13 Symmetry number array added
!   8/13 Thermal conductivity arrays added
!   8/13 symspec added
!   9/13 Raman susceptibilities added
!  12/13 Strain derivative arrays separated from stresses
!   1/14 lBOzrl arrays added
!   1/14 EDIP2p and lEDIPpi added
!   4/14 lbornkin added
!   7/14 lEDIP3mod added
!   8/14 eigv and groupvelocity added
!   8/14 MEAM screening parameters made species specific
!  10/14 nd2cellcfg added
!  12/14 rtrm1 changed from scalar to array
!   2/15 nqatomcell and qatomxyz added
!   2/15 nlibnobo added
!   2/15 lmm3se added
!   3/15 bond type arrays added
!   4/15 Ghost supercell array added
!   7/15 External potential added
!   8/15 llowered added
!   9/15 translate noise added
!  12/15 nBOtapertype added
!   1/16 Species specific VDW radius added for cosmo
!   3/16 BO coordination potential added
!   5/16 Multiple mcswaps added
!   5/16 SPME arrays added
!   7/16 nullifying of pkim_model removed since this is a c_ptr now
!   8/16 IR intensity array added
!   8/16 OpenKIM changes added
!   9/16 hmssgio added
!   9/16 ncoshptr added
!  11/16 phonon pointers for local node added
!   1/17 optimisation atom pointer arrays added
!   1/17 neemlocptr added
!   2/17 parallel variable arrays added
!   3/17 fix_atom option added
!   4/17 ChemShell interaction modified
!   4/17 qonsas option added
!   5/17 Region 1 parallel distribution pointers added
!   6/17 Old ChemShell restored as a compile option
!  10/17 Initial coordinates added to restart info
!  11/17 msd arrays added
!   1/18 icosxsp, icosysp, icoszsp added
!   1/18 Grueneisen parameters added
!   1/18 nullifying of pkim_model removed since this is a c_ptr now - again
!   1/18 Trace added
!   2/18 nxks, nyks, nzks converted to a single array
!   2/18 Configuration arrays for tracking primitive cell added
!   3/18 Sign option added to translate
!   4/18 scan_cell arrays added
!   4/18 Split bonds added
!   5/18 nbondqb added
!   5/18 lcscanstrain added
!   5/18 qrange added
!   6/18 e0range added
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
  use bondcharge
  use bondorderdata
  use bondvectors
  use cellmultipole
  use chargecoupled
  use configurations
  use cosmic
  use cosmicpwtloc
  use current
  use defects
  use derivatives
  use dispersion
  use distances
  use eam
  use EDIPdata
  use eembonds
  use eemdata
  use energies,          only : siteenergy, eregion2region
  use feworkspace
  use field
  use fitting
  use four
  use freeze
  use frequencies
  use gaconf
  use general
  use genetic
  use gulpchemsh
  use gulpinput
  use kim_cellimages
  use kim_models
  use ksample
  use ksample_scatter
  use kspace
  use library
  use moldyn
  use molecule
  use montecarlo
  use m_pr
  use m_three
  use g_neb
  use observables
  use one
  use optimisation
  use parallel
  use partial
  use phononatoms
  use plane
  use polarise
  use potchange
  use potentialgrid
  use potentialinterpolation
  use potentialpoints
  use potentialsites
  use potentialxyz
  use projectdos
  use properties
  use radial
  use reallocate
  use realvectors
  use reaxFFdata
  use region2a
  use scan
  use scatterdata
  use shells
  use shellextrapolation
  use shifts
  use six
  use spatial
  use spatialbo
  use species
  use splinedata
  use spme,            only : nqkgrid
  use sutton
  use symmetry
  use thermalcond
#ifdef TRACE
  use trace,           only : trace_in, trace_out
#endif
  use transform
  use two
  use velocities
  use xcgc
  use uffdata
!
!  Local variables
!
  integer(i4)       :: ierror
#ifdef TRACE
  call trace_in('nullpointer')
#endif
!
#ifdef OLDCS
  if (ichemsh_qm < 0) then
#else
  if (ichemsh_link.eq.0) then
#endif
!     
!  GULP standalone application (not a ChemShell call)
! 
    nullify(symbolbondQ)
    nullify(nbondQspec1)
    nullify(nbondQspec2)
    nullify(nbondQtyp1)
    nullify(nbondQtyp2)
    nullify(bondQincrement)
!
    nullify(nBOspec1)
    nullify(nBOspec2)
    nullify(nBOspecA1)
    nullify(nBOspecA2)
    nullify(nBOspecR1)
    nullify(nBOspecR2)
    nullify(nBOspecQ0)
    nullify(nBOspecQ1)
    nullify(nBOspecQ2)
    nullify(nBOspecZ)
    nullify(nBOtaperQ)
    nullify(nBOtapertype)
    nullify(nBOtyp1)
    nullify(nBOtyp2)
    nullify(nBOtypA1)
    nullify(nBOtypA2)
    nullify(nBOtypR1)
    nullify(nBOtypR2)
    nullify(nBOtypQ0)
    nullify(nBOtypQ1)
    nullify(nBOtypQ2)
    nullify(nBOtypZ)
    nullify(nBOtypeA)
    nullify(nBOtypeR)
    nullify(nBOtypeQ)
    nullify(nBOtypeQ0)
    nullify(nBOtypeT)
    nullify(lBOzrlA)
    nullify(lBOzrlR)
    nullify(BOcombi)
    nullify(BOacoeff)
    nullify(BObcoeff)
    nullify(BOccoeffA)
    nullify(BOccoeffR)
    nullify(BOccoeffZ)
    nullify(BOchiA)
    nullify(BOchiR)
    nullify(BOecoeffA)
    nullify(BOecoeffR)
    nullify(BOecoeffZ)
    nullify(BOhcoeffA)
    nullify(BOhcoeffR)
    nullify(BOlcoeffA)
    nullify(BOlcoeffR)
    nullify(BOmcoeffA)
    nullify(BOmcoeffR)
    nullify(BOncoeffA)
    nullify(BOncoeffR)
    nullify(BOzcoeffZ)
    nullify(BOq0)
    nullify(BOq0pot)
    nullify(BOq0ref)
    nullify(BOq0rho)
    nullify(BOzacoeff)
    nullify(BOzbcoeff)
    nullify(rBOmax)
    nullify(rBOmin)
    nullify(rBOmaxQ)
    nullify(rBOminQ)
!
    nullify(natCCspec)
    nullify(ntypCCspec)
    nullify(nCCparNb)
    nullify(CCbeta)
    nullify(CCeta)
    nullify(CClambda)
    nullify(CCmu)
    nullify(CCparA)
    nullify(CCparB)
    nullify(CCparC)
    nullify(CCparD)
    nullify(CCparDL)
    nullify(CCparDU)
    nullify(CCparH)
    nullify(CCparM)
    nullify(CCparN)
    nullify(CCparIE)
    nullify(CCparAE)
    nullify(CCparQL)
    nullify(CCparQU)
    nullify(CCvdwC)
    nullify(rCCmaxL)
    nullify(rCCmaxS)
    nullify(rCCminL)
    nullify(rCCminS)
!
    nullify(anisotropicpresscfg)
    nullify(atom2local)
    nullify(atom2locala)
    nullify(atom2localv)
    nullify(atom2node)
    nullify(atom2nodea)
    nullify(atom2nodev)
    nullify(reg12local)
    nullify(reg12node)
    nullify(node2atom)
    nullify(node2atoma)
    nullify(node2atomv)
    nullify(node2pts)
    nullify(node2reg1)
    nullify(node2var)
    nullify(npts2local)
    nullify(npts2node)
    nullify(nvar2local)
    nullify(nvar2node)
    nullify(nboxat)
    nullify(ncoonnodeptr)
    nullify(nshonnodeptr)
    nullify(names)
    nullify(ioptcfg)
    nullify(maxmodecfg)
    nullify(minmodecfg)
    nullify(n1var)
    nullify(nd2cellcfg)
    nullify(ndimen)
    nullify(nascfg)
    nullify(natcfg)
    nullify(nbornstep)
    nullify(n1con)
    nullify(ncellmaxcfg)
    nullify(ncellmincfg)
    nullify(ncfixcfg)
    nullify(nconcfg)
    nullify(ncorepcfg)
    nullify(ncvarcfg)
    nullify(neamfnspecptr)
    nullify(neamspecptr)
    nullify(neembonded)
    nullify(neemptr)
    nullify(neemrptr)
    nullify(neemlocptr)
    nullify(neemlocrptr)
    nullify(neiglow)
    nullify(neighigh)
    nullify(nfixatomtypecfg)
    nullify(nfixatomcfg)
    nullify(ninternalmaxcfg)
    nullify(ninternalmincfg)
    nullify(nomegastep)
    nullify(nregionno)
    nullify(nregions)
    nullify(nregiontype)
    nullify(nrotop)
    nullify(nspecptrcfg)
    nullify(nsregion2)
    nullify(nsuper)
    nullify(nsuperghost)
    nullify(ntempstp)
    nullify(ntempstpstart)
    nullify(ntypcfg)
    nullify(ntwistcfg)
    nullify(nummodecfg)
    nullify(nvarcfg)
    nullify(nzmolcfg)
    nullify(QMMMmode)
    nullify(lanisotropicpresscfg)
    nullify(leinsteinat)
    nullify(lfcborn)
    nullify(lfcphon)
    nullify(lfcprop)
    nullify(lfcscatter)
    nullify(lgcmcmol)
    nullify(linitcfg)
    nullify(llowered)
    nullify(lomega)
    nullify(lopfc)
    nullify(lopfi)
    nullify(lopfreg)
    nullify(lqmatom)
    nullify(lregionrigid)
    nullify(ltdforcecfg)
    nullify(lvecin)
    nullify(cncfg)
    nullify(conaddcfg)
    nullify(concocfg)
    nullify(dhklcfg)
    nullify(energycfg)
    nullify(extpotcfg)
    nullify(forcecfg)
    nullify(keinsteinat)
    nullify(occucfg)
    nullify(omega)
    nullify(omegadamping)
    nullify(omegadir)
    nullify(omegadirtype)
    nullify(omegastep)
    nullify(oxcfg)
    nullify(qlcfg)
    nullify(presscfg)
    nullify(radcfg)
    nullify(ramandir)
    nullify(ramandirtype)
    nullify(rvcfg)
    nullify(rvpcfg)
    nullify(sbulkecfg)
    nullify(straincfg)
    nullify(stresscfg)
    nullify(symnocfg)
    nullify(symboleamfnspec)
    nullify(symboleamspec)
    nullify(tdforcecfg)
    nullify(tempcfg)
    nullify(tempstp)
    nullify(totalchargecfg)
    nullify(xcfg)
    nullify(ycfg)
    nullify(zcfg)
    nullify(xinitcfg)
    nullify(yinitcfg)
    nullify(zinitcfg)
    nullify(xeinsteinat)
    nullify(yeinsteinat)
    nullify(zeinsteinat)
    nullify(nbonds)
    nullify(nbonded)
    nullify(nbondedtype)
    nullify(nbondind)
    nullify(nbondqb)
    nullify(icosx)
    nullify(icosy)
    nullify(icosz)
    nullify(icosxsp)
    nullify(icosysp)
    nullify(icoszsp)
    nullify(iopt)
    nullify(iatn)
    nullify(nat)
    nullify(natype)
    nullify(nftype)
    nullify(ncfix)
    nullify(ncvar)
    nullify(neqv)
    nullify(nrel2)
    nullify(nrelat)
    nullify(nspecptr)
    nullify(lbsmat)
    nullify(lsliceatom)
    nullify(c6a)
    nullify(c6f)
    nullify(cna)
    nullify(cnf)
    nullify(conadd)
    nullify(conco)
    nullify(mass)
    nullify(msdx)
    nullify(msdy)
    nullify(msdz)
    nullify(occua)
    nullify(occuf)
    nullify(oxa)
    nullify(oxf)
    nullify(qa)
    nullify(qbond)
    nullify(qf)
    nullify(rada)
    nullify(radf)
    nullify(rmass)
    nullify(xalat)
    nullify(yalat)
    nullify(zalat)
    nullify(xclat)
    nullify(yclat)
    nullify(zclat)
    nullify(xafrac)
    nullify(yafrac)
    nullify(zafrac)
    nullify(xfrac)
    nullify(yfrac)
    nullify(zfrac)
    nullify(xfracimage)
    nullify(yfracimage)
    nullify(zfracimage)
    nullify(xfsave)
    nullify(yfsave)
    nullify(zfsave)
    nullify(xinitial)
    nullify(yinitial)
    nullify(zinitial)
    nullify(xstore)
    nullify(ystore)
    nullify(zstore)
    nullify(rstore)
    nullify(x0)
    nullify(nbondsdef)
    nullify(nbondeddef)
    nullify(nreldef)
    nullify(nreldefcfg)
    nullify(idopt)
    nullify(idoptcfg)
    nullify(inddfix)
    nullify(inddeffix)
    nullify(natdefe)
    nullify(natp)
    nullify(natdvacptrcfg)
    nullify(ntypdefe)
    nullify(ntypep)
    nullify(ntypdvacptrcfg)
    nullify(ncdfix)
    nullify(ncdvar)
    nullify(ndcentyp)
    nullify(ndefcfg)
    nullify(ndefind)
    nullify(ndefindp)
    nullify(ndefmol)
    nullify(ndefmolp)
    nullify(ndefnat)
    nullify(ndeftp)
    nullify(ndeftyp)
    nullify(ndeqv)
    nullify(ndintptr)
    nullify(ndvacptr)
    nullify(ndintptrcfg)
    nullify(ndvacptrcfg)
    nullify(ndrel)
    nullify(ndrelop)
    nullify(ndsptr)
    nullify(npsite)
    nullify(nptrr1)
    nullify(nreg1cfg)
    nullify(natdefecfg)
    nullify(ntypdefecfg)
    nullify(nintecfg)
    nullify(nvacacfg)
    nullify(xdefecfg)
    nullify(ydefecfg)
    nullify(zdefecfg)
    nullify(qdefecfg)
    nullify(occdefecfg)
    nullify(radefecfg)
    nullify(ndefmolcfg)
    nullify(ndefindcfg)
    nullify(lbrdvacptrcfg)
    nullify(ldefbsmatcfg)
    nullify(ldqmatomcfg)
    nullify(ldefbsmat)
    nullify(ldeffix)
    nullify(ldeflin)
    nullify(ldfix)
    nullify(ldqmatom)
    nullify(lr1created)
    nullify(lreldin)
    nullify(lindintptr)
    nullify(lindvacptr)
    nullify(dconco)
    nullify(dscrho)
    nullify(dscrhor2d)
    nullify(dscrhor2p)
    nullify(occdefe)
    nullify(occp)
    nullify(qdefe)
    nullify(qp)
    nullify(radefe)
    nullify(reg1)
    nullify(reg1last)
    nullify(reg2)
    nullify(reg2a1)
    nullify(xdefe)
    nullify(ydefe)
    nullify(zdefe)
    nullify(xperf)
    nullify(yperf)
    nullify(zperf)
    nullify(xdef)
    nullify(ydef)
    nullify(zdef)
    nullify(xdcent)
    nullify(ydcent)
    nullify(zdcent)
    nullify(xyzdvacptrcfg)
!
    nullify(atomicstress)
    nullify(derv2dk)
    nullify(derv2)
    nullify(derv2d)
    nullify(dervi)
    nullify(derv3)
    nullify(diagblock)
    nullify(dqds)
    nullify(d2qds2)
    nullify(dqdxyz)
    nullify(d2qdxyz2)
    nullify(d2qdxyzs)
    nullify(qatomxyz)
    nullify(raderv)
    nullify(eregion2region)
    nullify(siteenergy)
    nullify(sumatomicstress)
    nullify(xdrv)
    nullify(ydrv)
    nullify(zdrv)
    nullify(xregdrv)
    nullify(yregdrv)
    nullify(zregdrv)
!
    nullify(nqatoms)
    nullify(nqatomcell)
    nullify(nqatomptr)
    nullify(ndde)
    nullify(ndds)
    nullify(ndispcfg)
    nullify(ndstart)
    nullify(ndend)
    nullify(xdisp)
    nullify(ydisp)
    nullify(zdisp)
!
    nullify(ndenfn)
    nullify(ndenfncomp)
    nullify(neamfnnumeric)
    nullify(neamnat)
    nullify(neamnat2)
    nullify(neamfnnat)
    nullify(neamtyp)
    nullify(neamtyp2)
    nullify(neamfnmeamorder)
    nullify(neamfntyp)
    nullify(neammeamorder)
    nullify(denpar)
    nullify(eamalloy)
    nullify(eamdenfile)
    nullify(eamfnfile)
    nullify(eamfnnumeric) 
    nullify(eamfnnumeric1)
    nullify(eamfnnumeric2)
    nullify(eamfnnumeric3)
    nullify(eamfnnumeric4)
    nullify(eamfnnumeric5)
    nullify(eamfnnumeric6)
    nullify(eamfnnumeric7)
    nullify(eamfnnumeric8)
    nullify(eamfnnumericdrho)
    nullify(eamfnpar)
    nullify(eamfnmeamcoeff)
    nullify(eamtaperdrho)
    nullify(eamtaperrho)
    nullify(lmeamspec)
    nullify(lMEAMscreen)
    nullify(meam_Cmin)
    nullify(meam_Cmax)
!
    nullify(nfatyp)
    nullify(nfcfix)
    nullify(nfcotyp)
    nullify(nfcvar)
    nullify(nfcfg)
    nullify(nfpot)
    nullify(nfpot2)
    nullify(nfpot3)
    nullify(nftyp)
    nullify(nfvar)
    nullify(nfvar2)
    nullify(nfvar3)
    nullify(nfitptr)
    nullify(fconadd)
    nullify(fconco)
    nullify(fconpower)
    nullify(fres)
    nullify(scale)
    nullify(icell41)
    nullify(icell42)
    nullify(icell43)
    nullify(ilind)
    nullify(ilnum)
    nullify(jkind)
    nullify(symbol4)
    nullify(mmfexc)
    nullify(n4botype)
    nullify(nforptr)
    nullify(nforty)
    nullify(nfptyp1)
    nullify(nfptyp2)
    nullify(nfptyp3)
    nullify(nfptyp4)
    nullify(nfspec1)
    nullify(nfspec2)
    nullify(nfspec3)
    nullify(nfspec4)
    nullify(npfor)
    nullify(lfdreiding)
    nullify(lfintra)
    nullify(lfinter)
    nullify(lgenerated4)
    nullify(lonly3oop)
    nullify(loutofplane)
    nullify(fork)
    nullify(for1)
    nullify(for2)
    nullify(for3)
    nullify(for4)
    nullify(for1min)
    nullify(for2min)
    nullify(for3min)
    nullify(for4min)
    nullify(forpoly)
    nullify(iltor)
    nullify(ilftor)
    nullify(ilxtor)
    nullify(lopiltor)
    nullify(liltorswitch)
    nullify(ljktorswitch)
    nullify(lsurfiltor)
    nullify(nctor)
    nullify(neqiltor)
    nullify(nfortor)
    nullify(oiltor)
    nullify(riltor)
    nullify(xiltor)
    nullify(yiltor)
    nullify(ziltor)
!
    nullify(lfieldcfg)
    nullify(ntdfieldcfg)
    nullify(fieldcfg)
    nullify(fielddirectioncfg)
    nullify(td_fieldcfg)
    nullify(td_fielddirectioncfg)
    nullify(tmdfieldstart)
    nullify(tmdfieldstop)
!
    nullify(ibocptr)
    nullify(ibocshptr)
    nullify(iocptr)
    nullify(iocshptr)
    nullify(nbsptr)
!
    nullify(nphonatptr)
    nullify(nphonatcptr)
    nullify(nphonatsptr)
    nullify(nphonatrptr)
    nullify(nphonatrcptr)
    nullify(nphonatrsptr)
!
    nullify(nphonatonnodeptr)
    nullify(nphonatonnodecptr)
    nullify(nphonatonnodesptr)
    nullify(nphonatonnoderptr)
    nullify(nphonatonnodercptr)
    nullify(nphonatonnodersptr)
!
    nullify(nqrange)
    nullify(nqrangetype)
    nullify(nqrnow)
    nullify(qrangemax)
    nullify(qrangemin)
    nullify(chirange)
    nullify(murange)
    nullify(e0range)
    nullify(q0range)
    nullify(radrange)
    nullify(zetarange)
    nullify(znucrange)
!
    nullify(icell61)
    nullify(icell62)
    nullify(icell63)
    nullify(icell64)
    nullify(icell65)
    nullify(ijind)
    nullify(klind)
    nullify(mnind)
    nullify(symbol6)
    nullify(mmsexc)
    nullify(n6botype)
    nullify(nsixptr)
    nullify(nsixty)
    nullify(nsptyp1)
    nullify(nsptyp2)
    nullify(nsptyp3)
    nullify(nsptyp4)
    nullify(nsptyp5)
    nullify(nsptyp6)
    nullify(nsspec1)
    nullify(nsspec2)
    nullify(nsspec3)
    nullify(nsspec4)
    nullify(nsspec5)
    nullify(nsspec6)
    nullify(npsix)
    nullify(lsintra)
    nullify(lsinter)
    nullify(sixk)
    nullify(six1)
    nullify(six2)
    nullify(six3)
    nullify(six4)
    nullify(six5)
!
    nullify(eigv)
    nullify(freq)
    nullify(groupvelocity)
    nullify(grueneisen)
    nullify(IRintensity)
    nullify(titleword)
    nullify(ndiscret)
    nullify(xmax)
    nullify(xmaxcfg)
    nullify(xmin)
    nullify(xmincfg)
    nullify(words)
    nullify(norigkpt)
    nullify(nkptcfg)
    nullify(nlorder)
    nullify(floats)
    nullify(nks)
    nullify(nksala)
    nullify(lkptdispersion)
    nullify(wkpt)
    nullify(xkpt)
    nullify(ykpt)
    nullify(zkpt)
    nullify(wskpt)
    nullify(xskpt)
    nullify(yskpt)
    nullify(zskpt)
    nullify(indk)
    nullify(argc)
    nullify(csin)
    nullify(kmod)
    nullify(ktrm)
    nullify(ktrms)
    nullify(sine)
    nullify(xrk)
    nullify(yrk)
    nullify(zrk)
    nullify(xrk0)
    nullify(yrk0)
    nullify(zrk0)
    nullify(libname)
    nullify(libspec)
    nullify(lmdconstrain)
    nullify(nensemble)
    nullify(nmdconstrainatom)
    nullify(nmdconstraindist)
    nullify(nmdeq)
    nullify(nmdprod)
    nullify(nmdsamp)
    nullify(nmdvelmode)
    nullify(nmdvelmodp)
    nullify(nmdwrite)
    nullify(lfix)
    nullify(qpres)
    nullify(qtemp)
    nullify(tmdeq)
    nullify(tmdforcestart)
    nullify(tmdforcestop)
    nullify(tmdprod)
    nullify(tmdsamp)
    nullify(tmdscale)
    nullify(tmdscint)
    nullify(tmdwrite)
    nullify(tstep)
    nullify(xabsco)
    nullify(yabsco)
    nullify(zabsco)
    nullify(moldim)
    nullify(moldimi)
    nullify(molgcmc)
    nullify(ixshift)
    nullify(iyshift)
    nullify(izshift)
    nullify(natmol)
    nullify(n1connect)
    nullify(n2connect)
    nullify(nconnectcfg)
    nullify(nconnectind)
    nullify(nconnecttype)
    nullify(nmolatom)
    nullify(nmolind)
    nullify(nmollist)
    nullify(nmolptr)
    nullify(natbondtype)
    nullify(ntypbondtype)
    nullify(nbondtypeptr)
    nullify(nobond)
    nullify(nobotyp)
    nullify(nfgracfg)
    nullify(nfgrat)
    nullify(nfstraincfg)
    nullify(nfstraint)
    nullify(nobcfg)
    nullify(nobptr)
    nullify(nobptr2)
    nullify(nobptr3)
    nullify(nobtyp)
    nullify(nobsmodeat)
    nullify(nobsmodecfg)
    nullify(fcalc)
    nullify(fcalcoriginal)
    nullify(fgrad)
    nullify(fgradweight)
    nullify(finfo)
    nullify(fstrain)
    nullify(fstrainweight)
    nullify(fobs)
    nullify(fobsmode)
    nullify(fobsmodefreq)
    nullify(fobsmodeover)
    nullify(freaction)
    nullify(weight)
    nullify(lopf)
    nullify(rmode)
    nullify(natpolspec)
    nullify(ntyppolspec)
    nullify(noptatptr)
    nullify(noptatrptr)
    nullify(noptatlocptr)
    nullify(noptatlocrptr)
    nullify(dpolar)
    nullify(dpolspec)
    nullify(qpolar)
    nullify(qpolspec)
    nullify(npchng)
    nullify(nxpg)
    nullify(nypg)
    nullify(nzpg)
    nullify(xmaxpg)
    nullify(xminpg)
    nullify(ymaxpg)
    nullify(yminpg)
    nullify(zmaxpg)
    nullify(zminpg)
    nullify(npotsitecfg)
    nullify(xpotsite)
    nullify(ypotsite)
    nullify(zpotsite)
    nullify(vpotsite)
    nullify(npotptcfg)
    nullify(xpotpt)
    nullify(ypotpt)
    nullify(zpotpt)
    nullify(vpotpt)
    nullify(v2xyz)
    nullify(vx)
    nullify(vy)
    nullify(vz)
    nullify(v2xyz12)
    nullify(vx12)
    nullify(vy12)
    nullify(vz12)
    nullify(nprojit)
    nullify(nprojcfg)
    nullify(nprojnat)
    nullify(nprojtyp)
    nullify(nprojdb)
    nullify(nprojdef)
    nullify(nprojptr)
    nullify(nr2a)
    nullify(ntr2a)
    nullify(nmr2a)
    nullify(nmir2a)
    nullify(nps)
    nullify(ndsptr2a)
    nullify(ndeqv2a)
    nullify(ndrel2a)
    nullify(ndrelop2a)
    nullify(ldbr2a)
    nullify(xdis)
    nullify(ydis)
    nullify(zdis)
    nullify(xr2a)
    nullify(yr2a)
    nullify(zr2a)
    nullify(qr2a)
    nullify(or2a)
    nullify(rr2a)
    nullify(ncscan)
    nullify(lcscanstrain)
    nullify(ntran)
    nullify(ltranat)
    nullify(ltranatminus)
    nullify(ltrannoise)
    nullify(ltrantherm)
    nullify(trannoise)
    nullify(trantherm)
    nullify(cscan)
    nullify(xtran)
    nullify(ytran)
    nullify(ztran)
    nullify(ncsptr)
    nullify(ncoptr)
    nullify(ncoshptr)
    nullify(nshptr)
    nullify(natratiomspec)
    nullify(ntypratiomspec)
    nullify(ratiom)
    nullify(ratiomspec)
    nullify(csvector)
    nullify(nshcfg)
    nullify(shift)
    nullify(shscalecfg)
    nullify(symspec)
    nullify(natspec)
    nullify(numofspec)
    nullify(ntypspec)
    nullify(lbrspec)
    nullify(ldefshspec)
    nullify(lgastinlibspec)
    nullify(lgastinspec)
    nullify(linspec)
    nullify(lmassinspec)
    nullify(lnmrinspec)
    nullify(lvdwinspec)
    nullify(lqinspec)
    nullify(lmask)
    nullify(c6spec)
    nullify(qlspec)
    nullify(gastspec)
    nullify(massspec)
    nullify(nmrspec)
    nullify(radspec)
    nullify(vdwspec)
    nullify(nsplpt)
    nullify(nsplty)
    nullify(d1f)
    nullify(d2f)
    nullify(splf)
    nullify(splr)
    nullify(scrho)
    nullify(scrho12)
    nullify(hmssg)
    nullify(hmssgio)
    nullify(ifhr)
    nullify(iflags)
    nullify(ifso)
    nullify(iperm)
    nullify(ngocfg)
    nullify(nccscfg)
    nullify(nspcg)
    nullify(nspcgp)
    nullify(ivso)
    nullify(lsymset)
    nullify(ropcfg)
    nullify(vitcfg)
    nullify(icell31)
    nullify(icell32)
    nullify(i3ind)
    nullify(j3ind)
    nullify(k3ind)
    nullify(mmtexc)
    nullify(n3botype)
    nullify(n3bondno)
    nullify(n3bondnono)
    nullify(nthbptr)
    nullify(nthrty)
    nullify(ntptyp1)
    nullify(ntptyp2)
    nullify(ntptyp3)
    nullify(ntspec1)
    nullify(ntspec2)
    nullify(ntspec3)
    nullify(lgenerated3)
    nullify(lthetataper)
    nullify(symbol3)
    nullify(ltdreiding)
    nullify(ltintra)
    nullify(ltinter)
    nullify(thbk)
    nullify(theta)
    nullify(thetatapermax)
    nullify(thetatapermin)
    nullify(thr1)
    nullify(thr2)
    nullify(thr3)
    nullify(thr1min)
    nullify(thr2min)
    nullify(thr3min)
    nullify(thrho1)
    nullify(thrho2)
    nullify(thrho3)
    nullify(threepoly)
    nullify(tmat)
    nullify(ipot)
    nullify(mmexc)
    nullify(mmexcse)
    nullify(n2botype)
    nullify(natse)
    nullify(ncombipower)
    nullify(ntypse)
    nullify(nattab)
    nullify(ntypab)
    nullify(symbol2)
    nullify(nptype)
    nullify(nptyp1)
    nullify(nptyp2)
    nullify(nspec1)
    nullify(nspec2)
    nullify(lcombine)
    nullify(lgenerated2)
    nullify(lintra)
    nullify(linter)
    nullify(leshift)
    nullify(lgshift)
    nullify(lmm3se)
    nullify(atoma)
    nullify(atomb)
    nullify(twopot)
    nullify(epsilon)
    nullify(eshift)
    nullify(gshift)
    nullify(repcut)
    nullify(rhopot)
    nullify(rpot)
    nullify(rpot2)
    nullify(scale14)
    nullify(sigma)
    nullify(tapergrad)
    nullify(taperpot)
    nullify(tpot)
    nullify(iufree)
    nullify(lufree)
    nullify(rufree)
    nullify(xufree)
    nullify(velx)
    nullify(vely)
    nullify(velz)
    nullify(x2)
    nullify(y2)
    nullify(z2)
    nullify(x3)
    nullify(y3)
    nullify(z3)
    nullify(x4)
    nullify(y4)
    nullify(z4)
    nullify(x5)
    nullify(y5)
    nullify(z5)
!
    nullify(cellindex)
    nullify(nbotype)
    nullify(nbotype2)
    nullify(lbonded)
    nullify(l2bonds)
    nullify(l3bonds)
    nullify(lptrmol)
    nullify(deriv)
    nullify(deriv2)
    nullify(deriv3)
    nullify(derive0)
    nullify(derive)
    nullify(derive2)
    nullify(derive3)
    nullify(dist)
    nullify(dist2)
    nullify(dist3)
    nullify(d0i)
    nullify(d0j)
    nullify(d1i)
    nullify(d1j)
    nullify(d2i2)
    nullify(d2ij)
    nullify(d2j2)
    nullify(rderiv)
    nullify(rpd)
    nullify(rtrm1)
    nullify(rtrm2)
    nullify(rtrm3)
    nullify(rtrm32)
    nullify(xtmp)
    nullify(ytmp)
    nullify(ztmp)
    nullify(xtmp2)
    nullify(ytmp2)
    nullify(ztmp2)
    nullify(xtmp3)
    nullify(ytmp3)
    nullify(ztmp3)
!
    nullify(gc)
    nullify(xc)
    nullify(xcother)
!
    nullify(ithbest)
    nullify(xbest)
    nullify(xconf)
    nullify(fconf)
!
!  One-body pointers
!
    nullify(symbol1)
    nullify(nptyp11)
    nullify(nspec11)
    nullify(onepot)
!
!  MC pointers
!
    nullify(ngcmcnat)
    nullify(ngcmctype)
    nullify(ngcmcmolat)
    nullify(ngcmcmolnat)
    nullify(ngcmcmoltype)
    nullify(nptrdestroyable)
    nullify(nptrmoveable)
    nullify(nptrrotateable)
    nullify(nptrstrainable)
    nullify(nptrswapable)
    nullify(nptrtrialatom)
    nullify(nmcswapnat)
    nullify(nmcswappair)
    nullify(nmcswapspec)
    nullify(nmcswaptype)
    nullify(nswapable)
    nullify(ltrialatom)
    nullify(lmcswapany)
    nullify(pswap)
    nullify(xgcmcmol)
    nullify(ygcmcmol)
    nullify(zgcmcmol)
!
!  COSMO arrays
!
    nullify(nallnearsegptr)
    nullify(nallnearsegrptr)
    nullify(nar)
    nullify(nnearseg)
    nullify(nnearsegptr)
    nullify(nnearsegptrcell)
    nullify(npwt)
    nullify(npwtptr)
    nullify(npwtloc)
    nullify(npwtptrloc)
    nullify(nsasexcludemax)
    nullify(nsasexcludemin)
    nullify(nsasparticleptr)
    nullify(nsasparticlepartptr)
    nullify(nset)
    nullify(nsetf)
    nullify(atsrad)
    nullify(cosmoatomptr)
    nullify(cosmoA)
    nullify(cosmoBq)
    nullify(cosmoeigen)
    nullify(cosmoepsilon)
    nullify(cosmodrsolv)
    nullify(cosmorsolv)
    nullify(cosmotm)
    nullify(cosmopwt)
    nullify(cosmowt)
    nullify(lcosmoeigin)
    nullify(qonsas)
    nullify(qonsastarget)
    nullify(qsasparticles)
    nullify(sas)
    nullify(segweight)
    nullify(sphere1h)
    nullify(sphere1)
    nullify(sphere2)
    nullify(spxyz)
    nullify(spxyzouter)
!
    nullify(lbornkin)
    nullify(bornk)
    nullify(bornq)
!
    nullify(lbuffercell)
    nullify(natomcell)
    nullify(natomnodeptr)
    nullify(ncellnodeptr)
    nullify(nspcellat)
    nullify(nspcellat1ptr)
    nullify(nspcell2atptr)
    nullify(nspcellatptr)
    nullify(nspcellatptrcell)
    nullify(xinbox)
    nullify(yinbox)
    nullify(zinbox)
!
    nullify(lnebvaryspring)
    nullify(nebspring)
    nullify(nebspringmin)
    nullify(nebfinalcell)
    nullify(nebfinalradius)
    nullify(nebfinalxyz)
    nullify(nebreplicacell)
    nullify(nebreplicaradius)
    nullify(nebreplicacfgptr)
    nullify(nebreplicaxyz)
    nullify(nnebreplica)
    nullify(nnebreplicano)
!
    nullify(lbuffercellbo)
    nullify(natomcellbo)
    nullify(natomnodeptrbo)
    nullify(ncellnodeptrbo)
    nullify(nspcellatbo)
    nullify(nspcellat1ptrbo)
    nullify(nspcellatptrbo)
    nullify(nspcellatptrcellbo)
    nullify(xinboxbo)
    nullify(yinboxbo)
    nullify(zinboxbo)
!
    nullify(nptrfork)
    nullify(nptrforl)
    nullify(nptrmanyk)
    nullify(d33)
    nullify(d33r)
    nullify(d33i)
    nullify(d33s)
    nullify(d33rs)
    nullify(d33is)
    nullify(d34)
    nullify(d34r)
    nullify(d34i)
    nullify(d34s)
    nullify(d34rs)
    nullify(d34is)
!
    nullify(nbondvecind)
    nullify(nbtypevec)
    nullify(nbtype2vec)
    nullify(lbondedvec)
    nullify(l2bondsvec)
    nullify(l3bondsvec)
!
    nullify(icosxs)
    nullify(icosys)
    nullify(icoszs)
    nullify(ndistance)
    nullify(ndistancecell)
    nullify(ndistanceij)
    nullify(ndistanceind)
    nullify(ndistancemolonly)
    nullify(ndistanceptr)
    nullify(ndistancereset)
    nullify(ndistbotype)
    nullify(ndistbotype2)
    nullify(distl1bond)
    nullify(distl2bond)
    nullify(distl3bond)
    nullify(distlptrmol)
    nullify(distlself)
    nullify(distance)
    nullify(distance2)
    nullify(distancexyz)
!
    nullify(dRinterpolate)
    nullify(rdRinterpolate)
    nullify(FofR)
    nullify(dFofR)
!
    nullify(xshellsave)
    nullify(yshellsave)
    nullify(zshellsave)
!
    nullify(lradialcfg)
    nullify(radialKcfg)
    nullify(radialXYZcfg)
!
!  Thermal conductivity
!
    nullify(lomega_af_in)
    nullify(omega_af)
    nullify(B_pr_cfg)
    nullify(v_s_cfg)
    nullify(v_p_cfg)
!
!  ReaxFF arrays
!
    nullify(symbolreaxFFspec)
    nullify(natreaxFFspec)
    nullify(ntypreaxFFspec)
    nullify(nreaxFFfixQspecptr)
    nullify(nreaxFFval3)
    nullify(lreaxFFbocorrect)
    nullify(lreaxFFmorseinput)
    nullify(lreaxFFpboOK)
    nullify(lreaxFFtorsinput)
    nullify(qreaxFF)
    nullify(reaxFFr)
    nullify(reaxFFalpha)
    nullify(reaxFFeps)
    nullify(reaxFFrvdw)
    nullify(reaxFFgammaw)
    nullify(reaxFFpover)
    nullify(reaxFFpunder)
    nullify(reaxFFhincrement)
    nullify(reaxFFlp)
    nullify(reaxFFmorse)
    nullify(reaxFFoc1)
    nullify(reaxFFoc2)
    nullify(reaxFFuc1)
    nullify(reaxFFpboc)
    nullify(reaxFFval)
    nullify(reaxFFval1)
    nullify(reaxFFval3)
    nullify(reaxFFconj3)
    nullify(reaxFFhb3)
    nullify(reaxFFpen2)
    nullify(reaxFFpen3)
    nullify(reaxFFtor4)
    nullify(reaxFFDe)
    nullify(reaxFFpbe)
    nullify(reaxFFpbo)
    nullify(reaxFFrmax)
    nullify(reaxFFrmaxpair)
    nullify(reaxFFDeVDW)
    nullify(reaxFFalphaVDW)
    nullify(reaxFFr0VDW)
    nullify(reaxFFgammaVDW)
    nullify(reaxFFgammaQ)
!
!  UFF arrays
!
    nullify(symbolUFFspec)
    nullify(natUFFspec)
    nullify(ntypUFFspec)
    nullify(nUFFtype)
    nullify(UFFr)
    nullify(UFFtheta)
    nullify(UFFtor)
    nullify(UFFchi)
    nullify(UFFx)
    nullify(UFFd)
    nullify(UFFzeta)
    nullify(UFFZeff)
    nullify(UFFKoop)
    nullify(UFFthetaoop)
!
!  Plane potential arrays
!
    nullify(planepotsymbol)
    nullify(nplanepotpower)
    nullify(nplanepottype)
    nullify(natplanepot)
    nullify(ntypplanepot)
    nullify(planepot)
    nullify(planepotrmin)
    nullify(planepotrmax)
!
!  MD integrator related arrays
!
    nullify(pekin)
    nullify(taubcfg)
    nullify(tautcfg)
    nullify(pr_conscfg)
!
!  Scatter arrays
!
    nullify(n_qpoint)
    nullify(q_ordinate)
    nullify(Hold_Q)
    nullify(Hold_smq)
    nullify(Hold_tau)
    nullify(Qvector)
    nullify(scatlencoh)
    nullify(scatleninc)
    nullify(sofomega)
    nullify(sofomega_fit)
    nullify(tauvector)
!
!  EDIP arrays
!
    nullify(symbolEDIPspec)
    nullify(natEDIPspec)
    nullify(ntypEDIPspec)
    nullify(lEDIPpairOK)
    nullify(lEDIPtriadOK)
    nullify(lEDIPpi)
    nullify(lEDIP3mod)
    nullify(lEDIP3orig)
    nullify(EDIPrmaxpair)
    nullify(EDIPrmax)
    nullify(EDIPfhigh)
    nullify(EDIPflow)
    nullify(EDIPphigh)
    nullify(EDIPplow)
    nullify(EDIPalpha)
    nullify(EDIPZdih)
    nullify(EDIPZrep)
    nullify(EDIPc0)
    nullify(EDIP2epsilon)
    nullify(EDIP2a)
    nullify(EDIP2aprime)
    nullify(EDIP2B)
    nullify(EDIP2p)
    nullify(EDIP2beta)
    nullify(EDIP2sigma)
    nullify(EDIP3lambda0)
    nullify(EDIP3lambdap)
    nullify(EDIP3gamma0)
    nullify(EDIP3gammap)
    nullify(EDIP3q)
    nullify(EDIP3kq2)
    nullify(EDIP3Z0)
!
!  KIM related arrays
!
    nullify(lkim_model_cfg_OK)
    nullify(kim_cellvector)
    nullify(kim_coord)
    nullify(kim_cutoffs)
    nullify(kim_forces)
    nullify(kim_ncontributing)
    nullify(kim_ncontributing_only)
    nullify(kim_nspec)
    nullify(kim_nspecies)
    nullify(kim_particle)
!
!  Property arrays
!
    nullify(ramanasus)
!
!  SPME arrays
!
    nullify(nqkgrid)
  else
!        
!  ChemShell call
!        
    maxcfg = 1
    call realloc(nascfg,maxcfg,ierror)
    if (ierror.ne.0) call outofmemory('nullpointer','nascfg')
    nascfg = 0
#ifdef OLDCS
    nullify(shell_force)
#endif

    nboA = 0
    nboR = 0
    nboQ = 0
    nboQ0 = 0
    nbopot = 0
    nlibnboA = 0
    nlibnboR = 0
    nlibnboQ = 0
    nlibnboQ0 = 0
    nlibnbopot = 0
    nlibnobo = 0
     
  endif
#ifdef TRACE
  call trace_out('nullpointer')
#endif
!
  return
  end
