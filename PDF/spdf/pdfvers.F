*
* $Id$
*
* $Log$
* Revision 1.1.1.8  2000/04/19 13:29:31  mclareni
* Import version 8.04
*
* Revision 1.12  2000/04/17 18:54:47  plothow
* version 8.04: small mods, dlog to log
*
* Revision 1.11  2000/04/15 14:10:46  plothow
* Version 8.03: Nuclear PDFs introduced
*
* Revision 1.10  2000/04/09 14:47:30  plothow
* version 8.02 - add CTEQ4 ref2 PDFs
*
* Revision 1.9  2000/04/09 10:06:26  plothow
* new version 8.00 - MRS 98 and 99, CTEQ5 and GRV98 added
*
* Revision 1.8  2000/03/23 11:39:42  pdflib
* beta function changed to pdf_beta function
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
#include "pdf/w505190.inc"
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
C      CALL UHTOC(
C     + 78H PDF      7.09      970702 16.05 in the CERN Computer Program 
C     +Library  W5051  ,IQCHAW,CARD,80)
C
       PDFVER(1) = '8.04      '
       PDFVER(2) = '2000-04-17'
       PDFVER(3) = '12.24     '
C
C print here selected fields of the CARD string
         IF(N6.LE.0) N6 = L6
         IF(FIRST) THEN
C            WRITE(N6,*)  ' '
C            WRITE(N6,*)  ' *****  PDFLIB Version:  ',CARD(10:14),
C     +      '  Released on  ',CARD(20:26),'  at  ',CARD(27:32),
C     +      ' ',CARD(33:76),'  *****'
C     +      '  in the  ',CARD(33:62),'  *****'
C
C            WRITE(N6,*)  ' *****  Library compiled on ',IDATQQ,
C     +      ' at ',ITIMQQ,'  *****'
C
C            WRITE(N6,*)  ' *****  Library compiled on ',CARD(20:26),
C     +      ' at ',CARD(27:32),'  *****'
C
            WRITE(N6,*)  ' '
C
            WRITE(N6,*)  
     +      ' *****  CERN Computer Program Library     -    ',
     +      ' Reference:  W5051  *****'
C
            WRITE(N6,*)  ' *****  PDFLIB Version:  ',PDFVER(1)(1:4),
     +      '  Released on  ',PDFVER(2)(1:10),
     +      '  at  ',PDFVER(3)(1:5),'  *****'
C
         ENDIF
      RETURN
      END
