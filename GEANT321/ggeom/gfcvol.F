*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:48  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.28  by  S.Giani
*-- Author :
      SUBROUTINE GFCVOL
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    ROUTINE TO RELOAD THE RELEVANT CONTENTS OF THE COMMON       *
C.    *    BLOCK /GCVOLU/ FROM THE COMMON /GCVOL2/                     *
C.    *    THE ROUTINE ASSUMES THAT THE CONTENTS WERE STORED USING     *
C.    *    THE ROUTINE GSCVOL.                                         *
C.    *                                                                *
C.    *    ==>Called by : GDRAW, GDRVOL, GINVOL, GMEDIA, GTMEDI        *
C.    *         Authors : A.McPherson, F.Bruyant                       *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcvol2.inc"
C.
C.          ----------------------------------------
C.
      NLEVEL = NLEVE2
      IF (NLEVEL.LE.0) GO TO 99
*
      INFROM = INFRO2
      DO 40 I = 1,NLEVEL
*
         NAMES(I)  = NAMES2(I)
         NUMBER(I) = NUMB2(I)
         LVOLUM(I) = LVOL2(I)
         LINDEX(I) = LIND2(I)
         GONLY(I)  = GONLY2(I)
         NLDEV(I)  = NLDEV2(I)
         LINMX(I)  = LINMX2(I)
*
         DO 10 J = 1,3
            GTRAN(J,I) = GTRAN2(J,I)
   10    CONTINUE
*
         DO 20 J = 1,10
            GRMAT(J,I) = GRMAT2(J,I)
   20    CONTINUE
*
         IQ(JGPAR+I) = IQ(JGPAR2+I)
         LQ(JGPAR-I) = LQ(JGPAR2-I)
*
   40 CONTINUE
*                                                             END GFCVOL
   99 END
