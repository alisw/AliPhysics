*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:52  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:23  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZFRFZ(LUNFZ,CHOPT)
*
************************************************************************
*
*        To read the sequential file LUNFZ into the CWD
*     NB. A call to FZOPEN must preceede this call
* Input:
*   LUNFZ   Logical unit number of the FZ sequential access file
*   CHOPT   default, read all cycles for path CHPATH
*           'H' read only the highest cycle
*
* Called by <USER>
*
*  Author  : R.Brun DD/US/PD
*  Written : 14.05.86
*  Last mod: 21.10.91
*
************************************************************************
#include "zebra/zbcdch.inc"
#include "zebra/rzcl.inc"
#include "zebra/rzch.inc"
#include "zebra/rzk.inc"
#include "zebra/rzckey.inc"
#include "zebra/rzclun.inc"
      CHARACTER*(*) CHOPT
      CHARACTER*128 CHFORM
      CHARACTER*1  BSLASH
      DIMENSION IHDIR(4),ICDIR(KNMAX)
      EQUIVALENCE (ICDIR(1),KEY(1))
      LOGICAL RZSAME
*
*-----------------------------------------------------------------------
*
      BSLASH=CQALLC(61:61)
      IQUEST(1)=0
      IF(LQRS.EQ.0)GO TO 99
*
*           Check permission
*
      IFLAG=0
      CALL RZMODS('RZFRFZ',IFLAG)
      IF(IFLAG.NE.0)GO TO 99
      CALL UOPTC(CHOPT,'H',IOPTH)
*
*           Save CWD name
*
*
      CALL RZCDIR(CHWOLD,'R')
      ITIME=0
*
*           Read general header and find next RZ construct
*
  10  NH=KNSIZE
      CALL FZIN(LUNFZ,JQPDVS,0,0,'S',NH,IHEAD)
      IF(IQUEST(1).NE.0)GO TO 99
      IF(NH.NE.2)GO TO 10
      IF(IHEAD(1).NE.12345)GO TO 10
      NLPI0=IHEAD(2)
      NLPI=NLPI0
*
*           Read next directory
*
  20  CONTINUE
      ITIME=ITIME+1
      NH=KNSIZE
      CALL FZIN(LUNFZ,JQPDVS,0,0,'S',NH,IHEAD)
      IF(IQUEST(1).NE.0)GO TO 90
      IF(NH.EQ.1.AND.IHEAD(1).EQ.99.AND.ITIME.NE.0)GO TO 90
      IF(NH.NE.KNSIZE)GO TO 20
      IF(IHEAD(1).NE.1)GO TO 20
      IF(IHEAD(2).EQ.NLPI0)GO TO 70
*
*          Go back levels
*
      IF(IHEAD(2).LE.NLPI)THEN
         CALL MZDROP(JQPDVS,LCDIR,' ')
         CHL=BSLASH
         ICHL=1
         DO 30 I=NLPI-1,IHEAD(2),-1
            CHFORM=CHL(1:ICHL)//BSLASH
            CHL=CHFORM
            ICHL=ICHL+1
  30     CONTINUE
         CALL RZCDIR(CHL,' ')
      ENDIF
*
*           New subdirectory. Check if directory does not exist already
*
      LS=IQ(KQSP+LCDIR+KLS)
      NSDIR=IQ(KQSP+LCDIR+KNSD)
      CALL ZITOH(ICDIR,IHDIR,4)
      CALL UHTOC(IHDIR,4,CHL,16)
      NWKEY=ICDIR(KNWKEY)
      KTAGS=KKDES+(NWKEY-1)/10+1
      DO 40 I=1,NSDIR
         IF(RZSAME(ICDIR,IQ(KQSP+LCDIR+LS+7*(I-1)),4))GO TO 60
  40  CONTINUE
*
*           Create subdirectory
*
      CALL RZMDIR(CHL,NWKEY,'?',' ')
      IF(IQUEST(1).NE.0)GO TO 90
*
*           Set CWD to new branch
*
  60  CALL RZCDIR(CHL,' ')
*
*         Is directory big enough ?
*
      IF(IQ(KQSP+LCDIR-1).LT.ICDIR(KLE))THEN
         NM=ICDIR(KLE)-IQ(KQSP+LCDIR-1)
         CALL RZEXPD('RZFRFZ',NM)
         IF(IQUEST(1).NE.0) GO TO 90
      ENDIF
      CALL UCOPY(ICDIR(KKDES),IQ(KQSP+LCDIR+KKDES),2*NWKEY+KTAGS-KKDES)
      CALL UCOPY(ICDIR(KDATEC),IQ(KQSP+LCDIR+KDATEC),2)
*
  70  NLPI=IHEAD(2)
*
*           Copy keys from sequential file to CWD
*
      CALL SBIT1(IQ(KQSP+LTOP),2)
      CALL RZFRF1(LUNFZ,IOPTH)
      IF(IQUEST(1).EQ.0) GO TO 20
*
*           Set CWD to original value
*
  90  ISAVE = 2
      IQ1=IQUEST(1)
      CALL RZSAVE
      CALL RZCDIR(CHWOLD,' ')
      IQUEST(1)=IQ1
*
  99  RETURN
      END
