*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:27  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZUSED(N,IUSED)
*
************************************************************************
*
*        To mark used records
* Input:
*    N      Number of consecutive records used
*    IUSED  First of the N records used
*
* Called by RZCOP1,RZEXPD,RZMDIR,RZOUT,RZVOUT
*
*  Author  : R.Brun DD/US/PD
*  Written : 02.04.86
*  Last mod: 04.05.86
*
************************************************************************
#include "zebra/rzcl.inc"
*
*-----------------------------------------------------------------------
*
*           First. Remove used records from the list of available records
*
      NFREE=IQ(KQSP+LFREE+1)
      IF(NFREE.EQ.0)GO TO 99
*
      I=1
  10  IF(I.LE.NFREE)THEN
         L1=LFREE+2*I
         IR1=IQ(KQSP+L1)
         IRL=IQ(KQSP+L1+1)
         IF(IUSED.EQ.IR1)THEN
            IF(IR1+N.LE.IRL)THEN
               IQ(KQSP+L1)=IR1+N
            ELSE
               NFREE=NFREE-1
               IQ(KQSP+LFREE+1)=NFREE
               IF(NFREE.GT.0)THEN
                  CALL UCOPY2(IQ(KQSP+L1+2),IQ(KQSP+L1),2*(NFREE-I+1))
               ENDIF
            ENDIF
            GO TO 20
         ENDIF
         I=I+1
         GO TO 10
      ENDIF
*
*         Second. Mark used records in list of used records
*
  20  NUSED=IQ(KQSP+LUSED+1)
      IUSEL=IUSED+N-1
      IUS=2
      IF(NUSED.EQ.0)THEN
         IQ(KQSP+LUSED+1)=1
         IQ(KQSP+LUSED+IUS)=IUSED
         IQ(KQSP+LUSED+IUS+1)=IUSEL
         GO TO 99
      ENDIF
  30  IF(IUSEL.LT.IQ(KQSP+LUSED+IUS))THEN
         IF(IUSEL.EQ.IQ(KQSP+LUSED+IUS)-1)THEN
            IQ(KQSP+LUSED+IUS)=IUSED
         ELSE
            IF(IUS.GT.2.AND.IUSEL.LE.IQ(KQSP+LUSED+IUS-1))GO TO 99
            NDATA=IQ(KQSP+LUSED-1)
            IF(NDATA.LT.2*NUSED+3)THEN
               CALL MZPUSH(JQPDVS,LUSED,0,10,'I')
            ENDIF
            NLEFT=2*NUSED-IUS+2
            IF(NLEFT.GT.0)THEN
               CALL UCOPY2(IQ(KQSP+LUSED+IUS),
     +                     IQ(KQSP+LUSED+IUS+2),NLEFT)
            ENDIF
            NUSED=NUSED+1
            IQ(KQSP+LUSED+1)=NUSED
            IQ(KQSP+LUSED+IUS)=IUSED
            IQ(KQSP+LUSED+IUS+1)=IUSEL
         ENDIF
         GO TO 99
      ENDIF
*
      IF(IUSED.EQ.IQ(KQSP+LUSED+IUS+1)+1)THEN
         IQ(KQSP+LUSED+IUS+1)=IUSEL
         IF(IUS+2.LT.2*NUSED)THEN
            IF(IQ(KQSP+LUSED+IUS+1).EQ.IQ(KQSP+LUSED+IUS+2))THEN
               IQ(KQSP+LUSED+IUS+1)=IQ(KQSP+LUSED+IUS+3)
               NLEFT=2*NUSED-IUS-2
               IF(NLEFT.GT.0)THEN
                  CALL UCOPY2(IQ(KQSP+LUSED+IUS+4),
     +                        IQ(KQSP+LUSED+IUS+2),NLEFT)
               ENDIF
               NUSED=NUSED-1
               IQ(KQSP+LUSED+1)=NUSED
            ENDIF
         ENDIF
         GO TO 99
      ENDIF
*
      IUS=IUS+2
      IF(IUS.LE.2*NUSED)GO TO 30
*
      NDATA=IQ(KQSP+LUSED-1)
      IF(NDATA.LT.2*NUSED+3)THEN
         CALL MZPUSH(JQPDVS,LUSED,0,10,'I')
      ENDIF
      IQ(KQSP+LUSED+IUS)=IUSED
      IQ(KQSP+LUSED+IUS+1)=IUSEL
      NUSED=NUSED+1
      IQ(KQSP+LUSED+1)=NUSED
*
  99  RETURN
      END
