*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:40  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_NUCGEO_INC
#define CERNLIB_GEANT321_NUCGEO_INC
*
*
* nucgeo.inc
*
*
*=== nucgeo ===========================================================*
*
*----------------------------------------------------------------------*
*                                                                      *
*     NUClear GEOmetry common:                                         *
*                                                                      *
*     Created on  20  july  1991   by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on 24-mar-93     by    Alfredo Ferrari               *
*                                                                      *
*     Included in the following routines:                              *
*                                                                      *
*                        BIMSEL                                        *
*                        NUCNUC                                        *
*                        NWISEL                                        *
*                        PHDSET                                        *
*                        PIOABS                                        *
*                        PIONUC                                        *
*                        PRENUC                                        *
*                        PREPRE                                        *
*                        RSTSEL                                        *
*                        SBCOMP                                        *
*                        SIGFER                                        *
*                        UMOFIN                                        *
*                                                                      *
*     Description of the variables (NUCGEO):                           *
*                                                                      *
*                        Radtot = total radius of the nucleus          *
*                        Radiu0 = radius of the nucleus constant       *
*                                 density core                         *
*                        Radiu1 = radius at the nucleus skin depth     *
*                                 end                                  *
*                        Rad1o2 = half density radius of the nucleus   *
*                        Skindp = Skin depth of the nucleus ( where    *
*                                 density decreases linearly with the  *
*                                 radius from rhocen to rhoskn,        *
*                                 Radiu1 = Radiu0 + Skindp and         *
*                                 Rad1o2 = Radiu0 + 1/(2 Omalhl)       *
*                                        x Skindp )                    *
*                                 Rhoskn = Rhocen * Alphal             *
*                        Halodp = Halo depth of the nucleus ( where    *
*                                 density decreases linearly with the  *
*                                 radius from rhoskn to 0,             *
*                                 Radtot = Radiu0 + Skindp + Halodp    *
*                        Alphal = fraction of the central density the  *
*                                 transition from skin to halo occurs  *
*                                 at                                   *
*                        Omalhl = 1 - Alphal                           *
*                        Radskn = Radius at which the density would be *
*                                 zero if the skin depth behaviour is  *
*                                 continued (Radskn = Radiu0 + Skindp  *
*                                 / Omalhl)                            *
*                        Skneff = "effective" skin depth corresponding *
*                                 to Radskn, Skneff = Skindp / Omalhl  *
*                        Radpro = equivalent radius of the projectile  *
*                        Bimptr = "true" impact parameter of the proj- *
*                                 ectile (referred to the centre of    *
*                                 the projectile)                      *
*                        Rimptr = "true" radius of the interaction     *
*                                 point (referred to the centre of     *
*                                 the projectile)                      *
*                    X,Y,Zimptr = "true" position of the interaction   *
*                                 point (referred to the centre of     *
*                                 the projectile)                      *
*                        Rhocen = central density of the nucleus       *
*                        Rhocor = density at the transition point from *
*                                 core to skin                         *
*                        Rhoskn = density at the transition point from *
*                                 skin to halo, Rhoskn = Alphal Rhocen *
*                        Rhoimp = density of the nucleus at the "ef-   *
*                                 fective interaction point"           *
*                        Rhoimt = density of the nucleus at the "true" *
*                                 position of the interaction point    *
*                    Pfr,Ekfcen = Maximum Fermi momentum/energy in the *
*                                 central core                         *
*                    Pfr,Ekfpro = Maximum Fermi momentum/energy at the *
*                                 true interaction point for the proj- *
*                                 ectile (they are computed as for nu- *
*                                 cleons without any reduction factor) *
*                    Pfr,Ekfimp = Maximum Fermi momentum/energy at the *
*                                 "effective" interaction point for    *
*                                 the target nucleon                   *
*                    Pfr,Ekfbim = Maximum Fermi momentum/energy at     *
*                                 r = " effective" impact parameter    *
*                                 (Bimpct) for the nucleon with the    *
*                                 deepest well                         *
*                        Vprwll = well depth for the present projectile*
*                                 at the "true" interaction point      *
*                        Bimpct = "effective" impact parameter of the  *
*                                 projectile (referred to the centre   *
*                                 of the target nucleon)               *
*                        Rimpct = "effective" radius of the interac-   *
*                                 tion point (referred to the centre   *
*                                 of the target nucleon)               *
*                    X,Y,Zimpct = "effective" position of the interac- *
*                                 tion point (referred to the centre   *
*                                 of the target nucleon)               *
*                        Wllred = reduction factor to be applied to    *
*                                 the Ipwell well to get the proper    *
*                                 well for the projectile              *
*                        Clmbbr = Coulomb barrier for the present pro- *
*                                 jectile                              *
*                        Rdclmb = radius corresponding to the Coulomb  *
*                                 barrier at which Coulomb effects are *
*                                 supposed to be overcome by the nuc-  *
*                                 lear potential:                      *
*                                     Rdclmb = Clmbbr / (zZe^2)        *
*                        Bfclmb = correction factor for the impact pa- *
*                                 rameter, for boo such that the actual*
*                                 b =< Rdclmb:                         *
*                                  Bfclmb = sqrt ( 1 - Clmbbr/Ekproj ) *
*                        Bfceff = actual correction factor for the imp-*
*                                 act parameter, b = boo / Bfceff      *
*                                 for boo =< Rdclmb x Bfclmb:          *
*                                        Bfceff = Bfclmb               *
*                                 for boo > Rdclmb x Bfclmb:           *
*                                   Bfceff = 1 / ( x + sqrt (1+x^2) )  *
*                                  x = Clmbbr x Rdclmb / (2 Ekproj boo)*
*                        Ipwell = index of the target nucleon well to  *
*                                 be used in computing the one for the *
*                                 projectile                           *
*                                 ( 1 = proton, 2 = neutron )          *
*                        Itncmx = index of the target nucleon with     *
*                                 largest Fermi momentum               *
*                                 ( 1 = proton, 2 = neutron )          *
*                        Kprin  = particle index of the projectile     *
*                        Ntargt = number of target nucleons (2 at max) *
*                        Knucim = particle index of the target nucleon *
*                                 ( 1 = proton, 8 = neutron )          *
*                        Knuci2 = particle index of the 2nd target     *
*                                 nucleon for absorption on a couple   *
*                                 of nucleons                          *
*                                                                      *
*     Description of the variables (NUCPWI):                           *
*                                                                      *
*                        Almbar = Reduced De Broglie wavelength        *
*                        Bimmax = maximum impact parameter (at oo from *
*                                 the nucleus)                         *
*                        Siggeo = Geometrical cross section summed     *
*                                 over all partial waves (assuming     *
*                                 opacity=1 for any l)                 *
*                                 Siggeo = pi ( Almbar(lmax+1) )^2     *
*                        Lllmax = highest partial wave, it corresponds *
*                                 to: Almbar Lllmax >= Bimmax, where   *
*                                 the >= means that the smallest       *
*                                 integer >= Bimmax / Almbar is used   *
*                        Lllact = partial wave index of the present    *
*                                 interaction                          *
*                                                                      *
*----------------------------------------------------------------------*
*
      PARAMETER ( PI    = 3.14159265358979322702D+00 )
      PARAMETER ( PISQ  = 9.86960440108935854694D+00 )
      PARAMETER ( SQRT12 = 3.464101615137755D+00  )
*  This is log(11)/(4log(3)), it is ok for alphal = 0.1, bethal = 0.01 )
      PARAMETER ( SKTOHL = 0.5456645846610345D+00 )
*  This is log(99/19)/(4log(3)), it is ok for alphal = 0.05,
*  bethal = 0.01 )
*     PARAMETER ( SKTOHL = 0.3756286198494407D+00 )
*  This is log(99/4)/(4log(3)), it is ok for alphal = 0.2,
*  bethal = 0.01 )
*     PARAMETER ( SKTOHL = 0.7301997078753058D+00 )
*  This is log(99/4)/(4log(3)), it is ok for alphal = 0.02,
*  bethal = 0.001 )
*     PARAMETER ( SKTOHL = 0.7301997078753058D+00 )
      PARAMETER ( RZNUCL = 1.12        D+00 )
      PARAMETER ( RMSPRO = 0.8         D+00 )
      PARAMETER ( R0PROT = RMSPRO / SQRT12  )
      PARAMETER ( ARHPRO = 1.D+00 / 8.D+00 / PI / R0PROT / R0PROT
     &          / R0PROT )
      PARAMETER ( RLLE04 = RZNUCL )
      PARAMETER ( RLLE16 = RZNUCL )
      PARAMETER ( RLGT16 = RZNUCL )
      PARAMETER ( RCLE04 = 0.75D+00 / PI / RLLE04 / RLLE04 / RLLE04 )
      PARAMETER ( RCLE16 = 0.75D+00 / PI / RLLE16 / RLLE16 / RLLE16 )
      PARAMETER ( RCGT16 = 0.75D+00 / PI / RLGT16 / RLGT16 / RLGT16 )
      PARAMETER ( SKLE04 = 1.4D+00 )
      PARAMETER ( SKLE16 = 1.9D+00 )
      PARAMETER ( SKGT16 = 2.4D+00 )
      PARAMETER ( HLLE04 = SKTOHL * SKLE04 )
      PARAMETER ( HLLE16 = SKTOHL * SKLE16 )
      PARAMETER ( HLGT16 = SKTOHL * SKGT16 )
      PARAMETER ( ALPHA0 = 0.1D+00 )
*2    PARAMETER ( ALPHA0 = 0.05D+00 )
*0    PARAMETER ( ALPHA0 = 0.2D+00 )
*3    PARAMETER ( ALPHA0 = 0.02D+00 )
*4    PARAMETER ( ALPHA0 = 0.25D+00 )
      PARAMETER ( OMALH0 = 1.D+00 - ALPHA0 )
      PARAMETER ( GAMSK0 = 0.9D+00 )
*0    PARAMETER ( GAMSK0 = 0.8D+00 )
*2    PARAMETER ( GAMSK0 = 0.9D+00 )
*3    PARAMETER ( GAMSK0 = 0.9D+00 )
*4    PARAMETER ( GAMSK0 = 0.75D+00 )
      PARAMETER ( OMGAS0 = 1.D+00 - GAMSK0 )
      PARAMETER ( POTME0 = 0.6666666666666667D+00 )
      PARAMETER ( POTBA0 = 1.D+00 )
*  This parameter is the Panofsky ratio
      PARAMETER ( PNFRAT = 1.533D+00 )
*  This parameter set the branching ratio for radiative pi- capture
*  at rest in complex nuclei (it is a bit larger than the experimental
*  one to compensate for Pauli blocking etc etc ). This value is
*  for A=oo
      PARAMETER ( RADPIM = 0.035D+00 )
      PARAMETER ( RDPMHL = 14.D+00   )
*  Probability for pi- absorption to have a second proton in the couple
*  of nucleons ( the first one must be a proton )
      PARAMETER ( APMRST = 4.D+00 / 44.D+00 )
*  Probability for pi- absorption to have a second proton in the couple
*  of nucleons ( the first one must be a proton )
      PARAMETER ( APMPRO = 1.D+00 / 6.D+00 )
*  Probability for pi+ absorption to have a proton in the couple
*  of nucleons ( the first one must be a neutron )
      PARAMETER ( APPPRO = 1.D+00 / 6.D+00 )
*  Probability for pi0 absorption to have a p in the couple
*  as first particle
      PARAMETER ( AP0PFS = 0.5D+00 )
*  Probability for pi0 absorption to have a p in the couple
*  as second particle for a p first particle
      PARAMETER ( AP0PFP = 1.D+00 / 3.D+00 )
*  Probability for pi0 absorption to have a p in the couple
*  as second particle for a n first particle
      PARAMETER ( AP0NFP = 2.D+00 / 3.D+00 )
      PARAMETER ( MXSCIN = 50     )
*
      LOGICAL LABRST, LELSTC, LINELS, LCHEXC, LABSRP, LABSTH
*  NUClear Geometry Input data
      COMMON / FKNUGI / RHOTAB (2:260), RHATAB (2:260), ALPTAB (2:260),
     &                  RADTAB (2:260), SKITAB (2:260), HALTAB (2:260),
     &                  SK3TAB (2:260), SK4TAB (2:260), HABTAB (2:260),
     &                  CWSTAB (2:260), EKATAB (2:260), PFATAB (2:260),
     &                  PFRTAB (2:260)
*  NUClear GEOmetry
      COMMON / FKNUGE / RADTOT, RADIU1, RADIU0, RAD1O2, SKINDP, HALODP,
     &                  ALPHAL, OMALHL, RADSKN, SKNEFF, CPARWS, RADPRO,
     &                  RADCOR, RADCO2, RADMAX, BIMPTR, RIMPTR, XIMPTR,
     &                  YIMPTR, ZIMPTR, RHOIMT, EKFPRO, PFRPRO, RHOCEN,
     &                  RHOCOR, RHOSKN, EKFCEN (2), PFRCEN (2), EKFBIM,
     &                  PFRBIM, RHOIMP, EKFIMP, PFRIMP, EKFIM2, PFRIM2,
     &                  VPRWLL, RIMPCT, BIMPCT, XIMPCT, YIMPCT, ZIMPCT,
     &                  XBIMPC, YBIMPC, ZBIMPC, CXIMPC, CYIMPC, CZIMPC,
     &                  SQRIMP, SIGMAP, SIGMAN, SIGMAA, RHORED, R0TRAJ,
     &                  R1TRAJ, SBUSED, SBTOT , SBRES , RHOAVE, EKFAVE,
     &                  PFRAVE, AVEBIN, ACOLL , ZCOLL , RADSIG, OPACTY,
     &                  EKECON, PNUCCO, EKEWLL, PPRWLL, PXPROJ, PYPROJ,
     &                  PZPROJ, EKFERM, PNFRMI, PXFERM, PYFERM, PZFERM,
     &                  EKFER2, PNFRM2, PXFER2, PYFER2, PZFER2, EKFER3
      COMMON / FKNUGE /
     &                  PNFRM3, PXFER3, PYFER3, PZFER3, RHOMEM, EKFMEM,
     &                  BIMMEM, WLLRED, VPRBIM, POTINC, POTOUT, EEXMIN,
     &                  EEXDEL, EEXANY, CLMBBR, RDCLMB, BFCLMB, BFCEFF,
     &                  BNPROJ, BNDNUC, DEBRLM, SK4PAR, UBIMPC, VBIMPC,
     &                  WBIMPC, BNDPOT, SIGMAT, WLLRES, POTBAR, POTMES,
     &                  HHLP   (2), FORTOT (2), BNENRG (3), DEFNUC (2),
     &                  SIGMPR (3), SIGMNU (3), SIGMAB (3), IPWELL,
     &                  ITNCMX, KPRIN , NTARGT, KNUCIM, KNUCI2, KNUCI3,
     &                  IEVPRE, ISFCOL
*  NUClear Partial Waves Informations
      COMMON / FKNUPW / ALMBAR, BIMMAX, SIGGEO, LLLMAX, LLLACT
*  NUClear Geometry Interaction Informations
      COMMON / FKNUII / HOLEXP (MXSCIN), XEXPIN (MXSCIN),
     &                  YEXPIN (MXSCIN), ZEXPIN (MXSCIN),
     &                  RHOEXP, EKFEXP, EHLFIX,
     &                  NHLEXP, NHLFIX, IPRTYP, ISCTYP (MXSCIN),
     &                  NUSCIN, NEXPEM, LABRST, LELSTC, LINELS, LCHEXC,
     &                  LABSRP, LABSTH
      DIMENSION AWSTAB (2:260)
      EQUIVALENCE ( DEFPRO, DEFNUC (1) )
      EQUIVALENCE ( DEFNEU, DEFNUC (2) )
      EQUIVALENCE ( OMALHL, SK3PAR )
      EQUIVALENCE ( ALPHAL, HABPAR )
      EQUIVALENCE ( ALPTAB (2), AWSTAB (2) )
      EQUIVALENCE ( SIGMPE, SIGMPR (1) )
      EQUIVALENCE ( SIGMPC, SIGMPR (2) )
      EQUIVALENCE ( SIGMPI, SIGMPR (3) )
      EQUIVALENCE ( SIGMNE, SIGMNU (1) )
      EQUIVALENCE ( SIGMNC, SIGMNU (2) )
      EQUIVALENCE ( SIGMNI, SIGMNU (3) )
      EQUIVALENCE ( SIGMA2, SIGMAB (1) )
      EQUIVALENCE ( SIGMA3, SIGMAB (2) )
      EQUIVALENCE ( SIGMAS, SIGMAB (3) )
 

#endif
