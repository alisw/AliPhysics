*
* $Id$
*
* $Log$
* Revision 1.1.1.7  1997/07/13 13:25:18  cernlib
* Import version 7.09
*
* Revision 1.7  1997/07/02 15:03:08  pdflib
* potential bug removed in GRSGALO
*
* Revision 1.6  1997/04/08 15:53:26  pdflib
* Lambda_QCD for CTEQ4 sets corrected
*
* Revision 1.5  1996/12/10 10:49:31  pdflib
* Version 7.07: consistency change for u,d-valence in CTEQ4
*
* Revision 1.4  1996/11/15 15:07:21  pdflib
* Version 7.06: tabs removed and bug fixed in grsgalo
*
* Revision 1.3  1996/11/07 09:36:12  pdflib
* New version 7.05: cut routine Ctq4RTb in 9 pieces
*
* Revision 1.1.1.1  1996/04/12 15:29:26  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
      SUBROUTINE PDFVERS
C
#include "pdf/w5051p7.inc"
#include "pdf/w50516.inc"
#include "pdf/w50517.inc"
      CHARACTER*80 CARD
      IDATQQ = 0
      ITIMQQ = 0
#if (!defined(CERNLIB_CRAY))&&(!defined(CERNLIB_CONVEX))
      IQCHAW = 4
#endif
#if defined(CERNLIB_CRAY)||defined(CERNLIB_CONVEX)
      IQCHAW = 8
#endif
      CALL CFILL(' ',CARD,1,80)
      CALL UHTOC(
     + 78H PDF      7.09      970702 16.05 in the CERN Computer Program 
     +Library  W5051  ,IQCHAW,CARD,80)
C
C print here selected fields of the CARD string
         IF(N6.LE.0) N6 = L6
         IF(FIRST) THEN
            WRITE(N6,*)  ' '
            WRITE(N6,*)  ' *****  PDFLIB Version:  ',CARD(10:14),
     +      '  Released on  ',CARD(20:26),'  at  ',CARD(27:32),
     +      ' ',CARD(33:76),'  *****'
C     +      '  in the  ',CARD(33:62),'  *****'
C            WRITE(N6,*)  ' *****  Library compiled on ',IDATQQ,
            WRITE(N6,*)  ' *****  Library compiled on ',CARD(20:26),
     +      ' at ',CARD(27:32),'  *****'
C     +      ' at ',ITIMQQ,'  *****'
         ENDIF
      RETURN
      END
