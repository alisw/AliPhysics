*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:12  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.20  by  S.Giani
*-- Author :
      SUBROUTINE GSDETU(IUSET,IUDET,NUPAR,UPAR)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *      Storing and Retrieving User Detector parameters           *
C.    *      -----------------------------------------------           *
C.    *                                                                *
C.    *   Stores user parameters for detector IUDET of set IUSET.      *
C.    * IUSET     user set identifier                                  *
C.    * IUDET     user detector identifier                             *
C.    * NUPAR     number of user parameters                            *
C.    * UPAR      array of NUPAR user floating point parameters.       *
C.    *   The  routine  is  used at  initialisation  time  once  the   *
C.    * geometrical volumes have been defined.                         *
C.    *                                                                *
C.    *                                                                *
C.    *       JS = LQ(JSET-ISET)                                       *
C.    *       JD = LQ(JS-IDET)                                         *
C.    *       JDU= LQ(JD-1)                                            *
C.    *       Q(JDU+1) = UPAR(1)  1st user parameter,etc               *
C.    *                                                                *
C.    *    ==>Called by : <USER>, UGEOM                                *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
      DIMENSION UPAR(1)
      CHARACTER*4 IUSET,IUDET
C.
C.    ------------------------------------------------------------------
C.
      IF(JSET.LE.0)GO TO 90
      NSET=IQ(JSET-1)
      IF(NSET.LE.0)GO TO 90
      CALL GLOOK(IUSET,IQ(JSET+1),NSET,ISET)
      IF(ISET.LE.0)GO TO 90
      JS=LQ(JSET-ISET)
      NDET=IQ(JS-1)
      IF(NDET.LE.0)GO TO 90
      CALL GLOOK(IUDET,IQ(JS+1),NDET,IDET)
      IF(IDET.LE.0)GO TO 90
      JD=LQ(JS-IDET)
C
      IF(NUPAR.LE.0)GO TO 90
      CALL MZBOOK(IXCONS,JDU,JD,-3,'SJDU',0,0,NUPAR,3,0)
C
      CALL UCOPY(UPAR,Q(JDU+1),NUPAR)
      GO TO 99
C
C              Error
C
  90  WRITE(CHMAIL,1000)IUSET,IUDET
      CALL GMAIL(0,0)
 1000 FORMAT(' ***** GSDETU ERROR FOR SET ',A4,' OR DETECTOR ',A4)
C
  99  RETURN
      END
