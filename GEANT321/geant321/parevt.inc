*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:41  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_PAREVT_INC
#define CERNLIB_GEANT321_PAREVT_INC
*
*
* parevt.inc
*
*                                                                      *
*=== parevt ===========================================================*
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     Created on 28-3-90 by A. Ferrari, INFN - Milan                   *
*                                                                      *
*     Last change  on  20 april  1992    by    Alfredo Ferrari         *
*                                                                      *
*           included in:                                               *
*                          fluka                                       *
*                          bdnopt                                      *
*                          corevt                                      *
*                          distnu                                      *
*                          epilog                                      *
*                          eventv                                      *
*                          evevap                                      *
*                          evvini                                      *
*                          evxtes                                      *
*                          ferevv                                      *
*                          kaskad                                      *
*                          incini                                      *
*                          nucevv                                      *
*                          pmprab                                      *
*           and in     :                                               *
*                          prenuc                                      *
*                          umofin                                      *
*                          bimsel                                      *
*                          nucnuc                                      *
*                          nwisel                                      *
*                          pioabs                                      *
*                          prepre                                      *
*                          rstsel                                      *
*                                                                      *
*           frdiff = fraction of the total number of events which are  *
*                    diffractive ones                                  *
*           pthdff = momentum threshold for diffractive events         *
*                    diffractive ones                                  *
*           fsprd0 = asymptotic reduction factor for cascade particle  *
*                    energy slope parameter                            *
*           fshpfn = modifying factor for bamjet momentum sharing jet  *
*                    function                                          *
*           rn1gsc = correlation factor between shower and grey part-  *
*                    ticle multiplicities                              *
*           rn2gsc = correlation factor between shower and grey part-  *
*                    ticle multiplicities                              *
*           dpower = power exponent used in the binomial distributions *
*                    for grey particles                                *
*           lpower = flag to signal if dpower is different from the    *
*                    default option (dpower=2)                         *
*           ldiffr = flag for particles undergoing diffractive events  *
*           linctv = flag to switch on/off cascade nucleons and exci-  *
*                    tation energy (it is important for testing the    *
*                    generator                                         *
*           levprt = flag for nuclear evaporation                      *
*           ldeexg = flag for nuclear deexcitation gammas              *
*           lgdhpr = flag for geometry dependent calculations in the   *
*                    preequilibrium part                               *
*           lpreex = flag for explicit treatment of the first inte-    *
*                    ctions in the preequilibrium part                 *
*           lhlfix = flag for "freezing" the hole depth of the holes   *
*                    produced in the explicitly treated interactions   *
*                    in the preequilibrium part                        *
*           lprfix = flag for "freezing" the already sampled particle  *
*                    exciton energy in the preequilibrium part         *
*           lparwv = flag to apply the Heisenberg principle according  *
*                    to the reduced DeBroglie wavelength in the calcu- *
*                    lation of the impact parameter and in computing   *
*                    the particle position after interactions          *
*           ilvmod = flag to select the level density option           *
*                                                                      *
*----------------------------------------------------------------------*
*
*
*  30% of events, randomly choosen, are diffractive events!!
*  Now 20% (J. Ranft & A. Ferrari, 25-3-90)
*
      PARAMETER ( FRDIFF = 0.2D+00 )
      PARAMETER ( PTHDFF = 8.5D+00 )
      PARAMETER ( ETHSEA = 1.0D+00 )
 
      LOGICAL LDIFFR, LINCTV, LEVPRT, LHEAVY, LDEEXG, LGDHPR, LPREEX,
     &        LHLFIX, LPRFIX, LPARWV, LPOWER
      COMMON / FKPARE / DPOWER, FSPRD0, FSHPFN, RN1GSC, RN2GSC,
     &                  LDIFFR (39),    LPOWER, LINCTV, LEVPRT, LHEAVY,
     &                  LDEEXG, LGDHPR, LPREEX, LHLFIX, LPRFIX, LPARWV,
     &                  ILVMOD
 

#endif
