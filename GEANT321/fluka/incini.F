*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:57  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.43  by  S.Giani
*-- Author :
*=== incini ===========================================================*
*                                                                      *
      SUBROUTINE INCINI
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*                                                                      *
*     Created on  10  june  1990   by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on 14-apr-93     by    Alfredo Ferrari               *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*
      PARAMETER ( PI = PIPIPI )
#include "geant321/fheavy.inc"
#include "geant321/inpdat2.inc"
#include "geant321/inpflg.inc"
#include "geant321/nucdat.inc"
#include "geant321/parevt.inc"
      COMMON / FKNUCO / HELP (2), HHLP (2), FTVTH (2), FINCX (2),
     &                  EKPOLD (2), BBOLD, ZZOLD, SQROLD, ASEASQ,
     &                  FSPRED, FEX0RD
      BBOLD  = - 1.D+10
      ZZOLD  = - 1.D+10
      SQROLD = - 1.D+10
      APFRMX = PLABRC * ( 9.D+00 * PI / 8.D+00 )**0.3333333333333333D+00
     &       / R0NUCL
      AMNUCL (1) = AMPROT
      AMNUCL (2) = AMNEUT
      AMNUSQ (1) = AMPROT * AMPROT
      AMNUSQ (2) = AMNEUT * AMNEUT
      AMNHLP = 0.5D+00 * ( AMNUCL (1) + AMNUCL (2) )
      ASQHLP = AMNHLP**2
*     ASQHLP = 0.5D+00 * ( AMNUSQ (1) + AMNUSQ (2) )
      AEFRMX = SQRT ( ASQHLP + APFRMX**2 ) - AMNHLP
      AEFRMA = 0.3D+00 * APFRMX**2 / AMNHLP * ( 1.D+00 - APFRMX**2 /
     &         ( 5.6D+00 * ASQHLP ) )
      AV0WEL = AEFRMX + EBNDAV
      EBNDNG (1) = EBNDAV
      EBNDNG (2) = EBNDAV
      AEXC12 = 0.001D+00 * FKENER ( ONEONE*12, SIXSIX )
      CEXC12 = 0.001D+00 * ENRG   ( ONEONE*12, SIXSIX )
      AMMC12 = 12.D+00 * AMUAMU + AEXC12
      AMNC12 = AMMC12 - 6.D+00 * AMELEC +
     &         FERTHO * 6.D+00**EXPEBN
      AEXO16 = 0.001D+00 * FKENER ( ONEONE*16, EIGEIG )
      CEXO16 = 0.001D+00 * ENRG   ( ONEONE*16, EIGEIG )
      AMMO16 = 16.D+00 * AMUAMU + AEXO16
      AMNO16 = AMMO16 - 8.D+00 * AMELEC +
     &         FERTHO * 8.D+00**EXPEBN
      AEXS28 = 0.001D+00 * FKENER ( ONEONE*28, ONEONE*14 )
      CEXS28 = 0.001D+00 * ENRG   ( ONEONE*28, ONEONE*14 )
      AMMS28 = 28.D+00 * AMUAMU + AEXS28
      AMNS28 = AMMS28 - 14.D+00 * AMELEC +
     &         FERTHO * 14.D+00**EXPEBN
      AEXC40 = 0.001D+00 * FKENER ( ONEONE*40, ONEONE*20 )
      CEXC40 = 0.001D+00 * ENRG   ( ONEONE*40, ONEONE*20 )
      AMMC40 = 40.D+00 * AMUAMU + AEXC40
      AMNC40 = AMMC40 - 20.D+00 * AMELEC +
     &         FERTHO * 20.D+00**EXPEBN
      AEXF56 = 0.001D+00 * FKENER ( ONEONE*56, ONEONE*26 )
      CEXF56 = 0.001D+00 * ENRG   ( ONEONE*56, ONEONE*26 )
      AMMF56 = 56.D+00 * AMUAMU + AEXF56
      AMNF56 = AMMF56 - 26.D+00 * AMELEC +
     &         FERTHO * 26.D+00**EXPEBN
      AEX107 = 0.001D+00 * FKENER ( ONEONE*107, ONEONE*47 )
      CEX107 = 0.001D+00 * ENRG   ( ONEONE*107, ONEONE*47 )
      AMM107 = 107.D+00 * AMUAMU + AEX107
      AMN107 = AMM107 - 47.D+00 * AMELEC +
     &         FERTHO * 47.D+00**EXPEBN
      AEX132 = 0.001D+00 * FKENER ( ONEONE*132, ONEONE*54 )
      CEX132 = 0.001D+00 * ENRG   ( ONEONE*132, ONEONE*54 )
      AMM132 = 132.D+00 * AMUAMU + AEX132
      AMN132 = AMM132 - 54.D+00 * AMELEC +
     &         FERTHO * 54.D+00**EXPEBN
      AEX181 = 0.001D+00 * FKENER ( ONEONE*181, ONEONE*73 )
      CEX181 = 0.001D+00 * ENRG   ( ONEONE*181, ONEONE*73 )
      AMM181 = 181.D+00 * AMUAMU + AEX181
      AMN181 = AMM181 - 73.D+00 * AMELEC +
     &         FERTHO * 73.D+00**EXPEBN
      AEX208 = 0.001D+00 * FKENER ( ONEONE*208, ONEONE*82 )
      CEX208 = 0.001D+00 * ENRG   ( ONEONE*208, ONEONE*82 )
      AMM208 = 208.D+00 * AMUAMU + AEX208
      AMN208 = AMM208 - 82.D+00 * AMELEC +
     &         FERTHO * 82.D+00**EXPEBN
      AEX238 = 0.001D+00 * FKENER ( ONEONE*238, ONEONE*92 )
      CEX238 = 0.001D+00 * ENRG   ( ONEONE*238, ONEONE*92 )
      AMM238 = 238.D+00 * AMUAMU + AEX238
      AMN238 = AMM238 - 92.D+00 * AMELEC +
     &         FERTHO * 92.D+00**EXPEBN
#if defined(CERNLIB_FDEBUG)
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Maximum Fermi momentum  : ',REAL(APFRMX),
     &                  ' GeV/c ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Maximum Fermi energy    : ',REAL(AEFRMX),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Average Fermi energy    : ',REAL(AEFRMA),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Average binding energy  : ',REAL(EBNDAV),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Nuclear well depth      : ',REAL(AV0WEL),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Excess  mass  for 12-C  : ',REAL(AEXC12),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Cameron E. m. for 12-C  : ',REAL(CEXC12),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Atomic  mass  for 12-C  : ',REAL(AMMC12),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Nuclear mass  for 12-C  : ',REAL(AMNC12),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Excess  mass  for 16-O  : ',REAL(AEXO16),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Cameron E. m. for 16-O  : ',REAL(CEXO16),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Atomic  mass  for 16-O  : ',REAL(AMMO16),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Nuclear mass  for 16-O  : ',REAL(AMNO16),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Excess  mass  for 40-Ca : ',REAL(AEXC40),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Cameron E. m. for 40-Ca : ',REAL(CEXC40),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Atomic  mass  for 40-Ca : ',REAL(AMMC40),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Nuclear mass  for 40-Ca : ',REAL(AMNC40),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Excess  mass  for 56-Fe : ',REAL(AEXF56),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Cameron E. m. for 56-Fe : ',REAL(CEXF56),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Atomic  mass  for 56-Fe : ',REAL(AMMF56),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Nuclear mass  for 56-Fe : ',REAL(AMNF56),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Excess  mass  for 107-Ag: ',REAL(AEX107),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Cameron E. m. for 107-Ag: ',REAL(CEX107),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Atomic  mass  for 107-Ag: ',REAL(AMM107),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Nuclear mass  for 107-Ag: ',REAL(AMN107),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Excess  mass  for 132-Xe: ',REAL(AEX132),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Cameron E. m. for 132-Xe: ',REAL(CEX132),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Atomic  mass  for 132-Xe: ',REAL(AMM132),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Nuclear mass  for 132-Xe: ',REAL(AMN132),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Excess  mass  for 181-Ta: ',REAL(AEX181),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Cameron E. m. for 181-Ta: ',REAL(CEX181),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Atomic  mass  for 181-Ta: ',REAL(AMM181),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Nuclear mass  for 181-Ta: ',REAL(AMN181),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Excess  mass  for 208-Pb: ',REAL(AEX208),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Cameron E. m. for 208-Pb: ',REAL(CEX208),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Atomic  mass  for 208-Pb: ',REAL(AMM208),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Nuclear mass  for 208-Pb: ',REAL(AMN208),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Excess  mass  for 238-U : ',REAL(AEX238),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Cameron E. m. for 238-U : ',REAL(CEX238),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Atomic  mass  for 238-U : ',REAL(AMM238),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
      WRITE ( LUNOUT,* )' **** Nuclear mass  for 238-U : ',REAL(AMN238),
     &                  ' GeV   ****'
      WRITE ( LUNOUT,* )
#endif
      AMHEAV (1) = AMUAMU + 1.D-03 * FKENER ( ONEONE, ZERZER )
      AMHEAV (2) = AMUAMU + 1.D-03 * FKENER ( ONEONE, ONEONE )
      AMHEAV (3) = 2.D+00 * AMUAMU + 1.D-03 * FKENER ( TWOTWO, ONEONE )
      AMHEAV (4) = 3.D+00 * AMUAMU + 1.D-03 * FKENER ( THRTHR, ONEONE )
      AMHEAV (5) = 3.D+00 * AMUAMU + 1.D-03 * FKENER ( THRTHR, TWOTWO )
      AMHEAV (6) = 4.D+00 * AMUAMU + 1.D-03 * FKENER ( FOUFOU, TWOTWO )
      ELBNDE (0) = 0.D+00
      DO 2000 IZ = 1, 100
         ELBNDE ( IZ ) = FERTHO *  IZ **EXPEBN
2000  CONTINUE
      IF ( LEVPRT ) THEN
#if defined(CERNLIB_FDEBUG)
         WRITE ( LUNOUT, * )' **** Evaporation from residual nucleus',
     &                      ' activated **** '
         IF ( LDEEXG ) WRITE ( LUNOUT, * )' **** Deexcitation gamma',
     &                      ' production activated **** '
         IF ( LHEAVY ) WRITE ( LUNOUT, * )' **** Evaporated "heavies"',
     &                      ' transport activated **** '
         IF ( IFISS .GT. 0 )
     &                 WRITE ( LUNOUT, * )' **** High Energy fission ',
     &                      ' requested & activated **** '
#endif
      ELSE
         LDEEXG = .FALSE.
         LHEAVY = .FALSE.
         IFISS  = 0
      END IF
      RETURN
*=== End of subroutine incini =========================================*
      END
