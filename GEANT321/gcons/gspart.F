*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:17  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.19  by  S.Giani
*-- Author :
      SUBROUTINE GSPART(IPART,NAPART,ITRTYP,AMASS,CHARGE,TLIFE,
     +            UBUF,NWBUF)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Store particle parameters                                *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GPART                                *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcmzfo.inc"
#include "geant321/gcunit.inc"
      DIMENSION UBUF(1)
      CHARACTER*(*) NAPART
      CHARACTER*20 NAME
C.
C.    ------------------------------------------------------------------
C.
      IF(IPART.LE.0)GO TO 99
      IF(JPART.LE.0)THEN
         CALL MZBOOK(IXCONS,JPART,JPART,1,'PART',NPART,NPART,0,3,0)
         IQ(JPART-5)=0
      ENDIF
      IF(IPART.GT.NPART)THEN
         CALL MZPUSH(IXCONS,JPART,IPART-NPART,0,'I')
         NPART=IPART
         JPA1=0
      ELSE
         JPA1=LQ(JPART-IPART)
         IF(JPA1.GT.0) THEN
            WRITE(CHMAIL,10000)
            CALL GMAIL(1,0)
            CALL GPPART(IPART)
            CALL MZDROP(IXCONS,LQ(JPART-IPART),' ')
         ENDIF
      ENDIF
      CALL MZBOOK(IXCONS,JPA,JPART,-IPART,'PART',2,2,NWBUF+9,IOPART,0)
C
      NAME=NAPART
      NCH=LNBLNK(NAME)
      IF(NCH.GT.0)THEN
         IF(NAME(NCH:NCH).EQ.'$')NAME(NCH:NCH)=' '
      ENDIF
      CALL UCTOH(NAME,IQ(JPA+1),4,20)
C
      Q(JPA + 6) = ITRTYP
      Q(JPA + 7) = AMASS
      Q(JPA + 8) = CHARGE
      Q(JPA + 9) = TLIFE
      IF(NWBUF.GT.0)CALL UCOPY(UBUF,Q(JPA+10),NWBUF)
C
      IF(JPA1.GT.0) THEN
         CALL GPPART(-IPART)
      ENDIF
C
  99  RETURN
10000 FORMAT(' *** GSPART ***: Warning, particle redefinition:')
      END
