*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:41  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_PARNUC_INC
#define CERNLIB_GEANT321_PARNUC_INC
*
*
* parnuc.inc
*
*
*=== parnuc ===========================================================*
*
*----------------------------------------------------------------------*
*                                                                      *
*     Created on  2 august 1991    by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on 25-mar-93     by    Alfredo Ferrari               *
*                                                                      *
*     included in :                                                    *
*                                                                      *
*                   Couset                                             *
*                   Nucnuc                                             *
*                   Nwisel                                             *
*                   Phdset                                             *
*                   Phdwll                                             *
*                   Pioabs                                             *
*                   Pionuc                                             *
*                   Pmprab                                             *
*                   Pncset                                             *
*                   Pncset                                             *
*                   Prenuc                                             *
*                   Prepre                                             *
*                   Sbcomp                                             *
*                   Umofin                                             *
*                                                                      *
*----------------------------------------------------------------------*
*
      PARAMETER ( PIGRK  = PIPIPI )
      PARAMETER ( ALEVEL = 8.D-03 )
*  Standard parameter for central density: it is exactly r0nucl of
*  Parevt
      PARAMETER ( RCNUCL = 1.12D+00 )
*  "Naive" r0 for sigmas
      PARAMETER ( R0SIG  = 1.3D+00 )
*  r0 for sigmas according to Kikuchi & Kawai
      PARAMETER ( R0SIGK = 1.5D+00 )
*  r0 for Coulomb barrier according to Kikuchi & Kawai
      PARAMETER ( RCOULB = 1.5D+00 )
*  Parameter for Coulomb barrier according to Hermes (Evap-5)
      PARAMETER ( COULBH = 0.88235D-03 )
*  Central density in nuclei
      PARAMETER ( RHONU0 = 0.75D+00 / PIGRK / RCNUCL / RCNUCL / RCNUCL )
*  For Taufor, Aefrav & Rhonuc the combinations 1 and 3 are ok. If the
*  Hannes expression for Taufor is used we must use Taufor from
*  combination 3 and Aefrav and Rhonuc from combination 1
*  Combination 3:
*3    PARAMETER ( TAUFOR = 3.0D+00 )
*  Combination 1:
*     PARAMETER ( TAUFOR = 5.0D+00 )
      PARAMETER ( TAUFOR = 10.0D+00 )
*  Thresholds for explicit interactions:
*    Primary projectile:
      PARAMETER ( EKEEXP = 0.03D+00 )
*    Secondary nucleons:
      PARAMETER ( EKREXP = 0.05D+00 )
*    Just for Pauli checks:
      PARAMETER ( EKEMNM = 0.01D+00 )
*
      PARAMETER ( NCPMX = 50 )
      COMMON / FKPARN / EKORI , PXORI , PYORI , PZORI , PTORI ,
     &                  ENNUC  (NCPMX), PNUCL  (NCPMX), EKFNUC (NCPMX),
     &                  XSTNUC (NCPMX), YSTNUC (NCPMX), ZSTNUC (NCPMX),
     &                  PXNUCL (NCPMX), PYNUCL (NCPMX), PZNUCL (NCPMX),
     &                  RSTNUC (NCPMX), RHNUCL (NCPMX), FREEPA (NCPMX),
     &                  CRRPAN (NCPMX), CRRPAP (NCPMX), BSTNUC (NCPMX),
     &                  WLRNUC (NCPMX), BNDGAV, DEFMIN, KPNUCL (NCPMX),
     &                  KRFNUC (NCPMX), KPORI , IBORI , IBNUCL, NPNUC ,
     &                  INUCTS (NCPMX), NNUCTS, ISFNUC (NCPMX)
 

#endif
