*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:58  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:25  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZKEYS(MAXDIM,MAXKEY,KEYS,NKEYS)
*
************************************************************************
*
*        Returns the list of keys currently in the CWD
* Input:
*   MAXDIM   The actual first dimension of output array KEYS.   It should
*            in principle be at least equal to the number of key elements
*            NWKEY as declared to RZMDIR.
*   MAXKEY   The actual second dimension of output array KEYS
* Output:
*   KEYS*    A 2-dimensional array  dimensioned KEYS(MAXDIM,MAXKEY).   It
*            will contain the key vectors associated with the CWD.
*            Its first index runs over the  key elements for a given key,
*            while its second index runs over the different keys.
*   NKEYS*   Actual number of keys present for the CWD.
*
* Called by <USER>
*
*  Author  : R.Brun DD/US/PD
*  Written : 05.04.86
*  Last mod: 09.02.93 JDS. Return NK and NWK in IQUEST in all cases
*
************************************************************************
#include "zebra/rzcl.inc"
#include "zebra/rzk.inc"
      DIMENSION KEYS(MAXDIM,MAXKEY)
*
*-----------------------------------------------------------------------
*

#include "zebra/q_jbyt.inc"

      IQUEST(1)=0
      NKEYS=0
      IF(LQRS.EQ.0)GO TO 99
      IF(LCDIR.EQ.0)GO TO 99
      LK = IQ(KQSP+LCDIR+KLK)
      NK = IQ(KQSP+LCDIR+KNKEYS)
      NWK= IQ(KQSP+LCDIR+KNWKEY)
      IF(NK.LE.MAXKEY)THEN
         NKEYS=NK
         IQUEST(11)=NK
         IQUEST(12)=NWK
      ELSE
         NKEYS=MAXKEY
         IQUEST(1)=1
         IQUEST(11)=NK
         IQUEST(12)=NWK
      ENDIF
      IF(NWK.GT.MAXDIM)THEN
         IQUEST(1)=1
         IQUEST(11)=NKEYS
         IQUEST(12)=NWK
      ENDIF
*
      DO 20 I=1,NKEYS
         K=LK+(NWK+1)*(I-1)
         DO 10 J=1,NWK
            IKDES=(J-1)/10
            IKBIT1=3*J-30*IKDES-2
            IF(JBYT(IQ(KQSP+LCDIR+KKDES+IKDES),IKBIT1,3).LT.3)THEN
               KEYS(J,I)=IQ(KQSP+LCDIR+K+J)
            ELSE
               CALL ZITOH(IQ(KQSP+LCDIR+K+J),KEYS(J,I),1)
            ENDIF
  10     CONTINUE
  20  CONTINUE
*
  99  RETURN
      END
