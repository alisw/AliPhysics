*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:11  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.20  by  S.Giani
*-- Author :
      SUBROUTINE GSDETA (IUSET, IUDET, IUALI, NWHI, NWDI, IALI)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *                 Handling Detector Aliases                      *
C.    *                 --------------------------                     *
C.    *                                                                *
C.    * Detector  'aliases'  can  be  specified  for  any  sensitive   *
C.    * detector for  which the  user needs to  store more  than one   *
C.    * type of hit.                                                   *
C.    *   Defines an  alias IUALI for  detector IUDET of  set IUSET.   *
C.    * Enters it in the JSET structure as an additional detector in   *
C.    * the corresponding set,  at the position IALI.  Copies to the   *
C.    * link  position  IALI  the  GSDET  parameter  bank  from  the   *
C.    * original detector  IUDET,  with empty  links to  the GSDETH,   *
C.    * GSDETD and GSDETU parameter banks.    The user can therefore   *
C.    * call  these   three  routines   again  with   the  arguments   *
C.    * appropriate to the  detector IUALI.  Several aliases  can be   *
C.    * defined  for the  same  detector through  as  many calls  to   *
C.    * GSDETA.                                                        *
C.    *                                                                *
C.    *    ==>Called by : <USER>, UGEOM                                *
C.    *       Author    F.Bruyant  *********                           *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcmzfo.inc"
      CHARACTER*4 IUSET,IUDET,IUALI
C.
C.    ------------------------------------------------------------------
C.
      IF (JSET.LE.0) GO TO 90
      NSET = IQ(JSET-1)
      IF (NSET.LE.0) GO TO 90
      CALL GLOOK(IUSET,IQ(JSET+1),NSET,ISET)
      IF (ISET.LE.0) GO TO 90
      JS = LQ(JSET-ISET)
      NDET = IQ(JS-1)
      IF (NDET.LE.0) GO TO 90
      CALL GLOOK(IUDET,IQ(JS+1),NDET,IDET)
      IF (IDET.LE.0) GO TO 90
      CALL GLOOK(IUALI,IQ(JS+1),NDET,IALI)
      IF (IALI.NE.0) GO TO 95
      JD = LQ(JS-IDET)
      NDATA = IQ(JD-1)
C
      CALL MZPUSH(IXCONS, JS, 1, 1,'I')
      NDET = NDET +1
      IALI = NDET
      CALL UCTOH(IUALI,IQ(JS+IALI),4,4)
C
      CALL MZBOOK(IXCONS,JD2,JS,-IALI, 'SEJD',4,4, NDATA,IOSEJD,0)
      IQ(JD2-5)=10*ISET+IALI
      JD = LQ(JS-IDET)
      CALL UCOPY (IQ(JD+1), IQ(JD2+1), NDATA)
      IQ(JD2+7) = NWHI
      IQ(JD2+8) = NWDI
      IQ(JD2+10) = IDET
      GO TO 99
C
   90 WRITE (CHMAIL, 1000) IUSET, IUDET
      CALL GMAIL(0,0)
 1000 FORMAT (' ***** GSDETA ERROR FOR SET ',A4,' OR DETECTOR ',A4)
      GO TO 99
   95 WRITE (CHMAIL, 2000) IUSET, IUALI
      CALL GMAIL(0,0)
 2000 FORMAT (' ***** GSDETA ERROR FOR SET ', A4, ' ALIAS NAME '
     +, A4, ' ALREADY USED')
C
  99  RETURN
      END
