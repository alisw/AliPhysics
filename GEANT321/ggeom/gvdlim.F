*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1999/05/18 15:55:17  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:57  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.31  by  S.Giani
*-- Author :
      SUBROUTINE GVDLIM(JVO,IN,IAXIS,CLOW,CHIGH,IERR)
C.
C.    **********************************************************
C.    *                                                        *
C.    *    ROUTINE TO FIND THE COORDINATE LIMITS ALONG AXIS    *
C.    *    IAXIS FOR THE INTH CONTENT OF THE VOLUME POINTED    *
C.    *    TO BY JVO. THE LIMITS ARE RETURNED AS CLOW AND      *
C.    *    CHIGH.IF THEY CANNOT BE FOUND THEN IERR IS SET TO 1.*
C.    *                                                        *
C.    *    ==>Called by : GGCLOS                               *
C.    *         Author  S.Giani  *********                     *
C.    *                                                        *
C.    **********************************************************
C.
#include "geant321/gcbank.inc"
      DIMENSION PAR(100),ATT(20)
C.
C.         -----------------------------------------
C.
      JIN=LQ(JVO-IN)
      IROT=Q(JIN+4)
      IVOT=Q(JIN+2)
      JVOT=LQ(JVOLUM-IVOT)
      ISH=Q(JVOT+2)
C
      CALL GFIPAR(JVOT,JIN,IN,NPAR,NATT,PAR,ATT)
C
C              X,Y OR Z COORDINATES.
C
      IF(IAXIS.GT.3) GO TO 20
C
      CALL GVDCAR(IAXIS,ISH,IROT,PAR,CL,CH,IERR)
      IF(IERR.NE.0) GO TO 999
C
      CLOW=CL+Q(JIN+IAXIS+4)
      CHIGH=CH+Q(JIN+IAXIS+4)
C
      GO TO 999
C
   20 CONTINUE
      IF(IAXIS.GT.5) GO TO 30
C
C          R COORDINATES.
C
      CALL GVDRAD(IAXIS,ISH,IROT,Q(JIN+5),PAR,CLOW, CHIGH,IERR)
 
C
      GO TO 999
C
   30 CONTINUE
C
C          PHI OR THETA COORDINATES.
C
      IF(IAXIS.EQ.6) CALL GVDPHI(ISH,IROT,Q(JIN+5), PAR,CLOW,CHIGH,
     +IERR)
C
      IF(IAXIS.EQ.7) CALL GVDTHE(ISH,IROT,Q(JIN+5), PAR,CLOW,CHIGH,
     +IERR)
C
  999 CONTINUE
      END
