*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:55  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.30  by  S.Giani
*-- Author :
      SUBROUTINE GSCVOL
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    ROUTINE TO STORE THE RELEVANT CONTENTS OF THE COMMON BLOCK  *
C.    *    /GCVOLU/ IN THE COMMON /GCVOL2/                             *
C.    *    THE COMMON BLOCK CAN BE RELOADED FROM ARRAY BY THE ROUTINE  *
C.    *    GFCVOL.                                                     *
C.    *                                                                *
C.    *    ==>Called by : GDRAW,GDRVOL,GINVOL,GMEDIA,GTMEDI            *
C.    *         Authors : R.Brun, F.Bruyant  *********                 *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcvol2.inc"
C.
C.          ----------------------------------------
C.
      NLEVE2 = NLEVEL
      IF (NLEVE2.LE.0) GO TO 99
*
      NLEVMX=IQ(LQ(JVOLUM)-1)
      IF (JGPAR2.EQ.0) THEN
         CALL MZBOOK (IXCONS, JGPAR2, JGPAR2, 1,'GPA2',
     +                                NLEVMX+1, 0, NLEVMX+1, 2, 0)
      ENDIF
*
      INFRO2 = INFROM
      IF(NLEVEL+1.GT.IQ(JGPAR2-1)) THEN
         CALL MZPUSH(IXCONS,JGPAR2,NLEVEL+1-IQ(JGPAR2-1),
     +                             NLEVEL+1-IQ(JGPAR2-1),' ')
      ENDIF
      DO 40 I = 1,NLEVEL
*
         NAMES2(I) = NAMES(I)
         NUMB2(I)  = NUMBER(I)
         LVOL2(I)  = LVOLUM(I)
         LIND2(I)  = LINDEX(I)
         GONLY2(I) = GONLY(I)
         NLDEV2(I) = NLDEV(I)
         LINMX2(I) = LINMX(I)
*
         DO 10 J = 1,3
            GTRAN2(J,I) = GTRAN(J,I)
   10    CONTINUE
*
         DO 20 J = 1,10
            GRMAT2(J,I) = GRMAT(J,I)
   20    CONTINUE
*
         IQ(JGPAR2+I) = IQ(JGPAR+I)
         LQ(JGPAR2-I) = LQ(JGPAR-I)
*
   40 CONTINUE
*                                                             END GSCVOL
   99 END
