*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:27:02  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:25  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZMDIR(CHDIR,NWKEY,CHFORM,CHTAG)
*
************************************************************************
*
*         To create a subdirectory of the CWD
* Input:
*   CHDIR   Character variable specifying the name of the directory to be
*           created.  All characters, but / \ * ~ ?  are allowed in a
*           directory name.
*   NWKEY   Number of words associated to a key (maximum 5)
*   CHFORM  Character variable describing each element  of the key vector
*           'B' Bit string but not zero
*           'H' Hollerith (4 characters)
*           'A' same as 'H' but RZLDIR
*           'I' Integer (nonzero)
*           Ex: CHFORM='IIH' for NWKEY=3 and the 2 first keys are integer
*               and the third one is Hollerith
*   CHTAG   Character array defined as CHARACTER*8 CHTAG(NWKEY).
*           Each  element of  the  array allows  the  description of  the
*           corresponding element in the key vector with a tag of up to 8
*           characters.
*
* Called by <USER>
*
*  Author  : R.Brun DD/US/PD
*  Written : 03.04.86
*  Last mod: 24.06.92 JDS - protection against invalid directory names
*          : 04.03.94 S.Banerjee (Change in cycle structure)
*
************************************************************************
#include "zebra/zunit.inc"
#include "zebra/rzcl.inc"
#include "zebra/rzclun.inc"
#include "zebra/rzk.inc"
#include "zebra/rzcycle.inc"
      CHARACTER*(*) CHDIR,CHFORM,CHTAG(*)
      DIMENSION IHDIR(4)
      LOGICAL RZSAME
*
*-----------------------------------------------------------------------
*
#include "zebra/q_jbyt.inc"
      IQUEST(1)=0
      IF(LQRS.EQ.0)GO TO 99
*
      NCD=LENOCC(CHDIR)
      IF(NCD.GT.16)NCD=16
*
*     Check directory name
*
      IF(NCD.EQ.0) THEN
         IQUEST(1) = 2
         GOTO 99
      ENDIF
 
      JX = ICLUNS(CHDIR,1,NCD)
      IF(JX.GT.0) THEN
         IQUEST(1) = 3
         GOTO 99
      ENDIF
 
      CALL VBLANK(IHDIR,4)
      CALL UCTOH(CHDIR,IHDIR,4,NCD)
      CALL ZHTOI(IHDIR,IHDIR,4)
      LOGLV=JBYT(IQ(KQSP+LTOP),15,3)-3
*
*                Check NWKEY.LE.KNMAX
*
      IF(NWKEY.LE.0.OR.NWKEY.GT.KNMAX)THEN
         IF(LOGLV.GE.-2) WRITE(IQLOG,1000)
 1000    FORMAT(' RZMDIR. NWKEY input value is invalid')
         IQUEST(1) =1
         IQUEST(11)=NWKEY
         GO TO 99
      ENDIF
*
*                Check if subdirectory already exists
*
      LS = IQ(KQSP+LCDIR+KLS)
      NSDIR=IQ(KQSP+LCDIR+KNSD)
      IF(NSDIR.GT.0)THEN
         DO 5 I=1,NSDIR
            IF(RZSAME(IHDIR,IQ(KQSP+LCDIR+LS+7*(I-1)),4))THEN
               IF(LOGLV.GE.-2) WRITE(IQLOG,2000) CHDIR
 2000          FORMAT(' RZMDIR. Already existing directory ',A)
               IQUEST(1)=4
               GO TO 99
            ENDIF
   5     CONTINUE
      ENDIF
*
*           Check if WRITE permission
*
      IFLAG=0
      CALL RZMODS('RZMDIR',IFLAG)
      IF(IFLAG.NE.0)GO TO 99
*
*           Check if enough space. Increase directory size by LREC
*           if required.
*
      NWFREE=IQ(KQSP+LCDIR+KNFREE)
      IF(NWFREE.LT.8)THEN
         CALL RZEXPD('RZMDIR',40)
      ENDIF
*
*     Allocate record for new directory
*
      CALL RZALLO('RZMDIR',1,IRNEXT)
      IF(IRNEXT.LE.0)GO TO 99
*
*           Insert new directory name into S,  Move K
*
      LD = IQ(KQSP+LCDIR+KLD)
      LK = IQ(KQSP+LCDIR+KLK)
      LF = IQ(KQSP+LCDIR+KLF)
      CALL UCOPY2(IQ(KQSP+LCDIR+LK),IQ(KQSP+LCDIR+LK+7),LF-LK)
      CALL UCOPY(IHDIR,IQ(KQSP+LCDIR+LK),4)
      IQ(KQSP+LCDIR+LK+4)=NCD
      IF (KVSCYC.EQ.0) THEN
         IQ(KQSP+LCDIR+LK+5)=IRNEXT+2**18
      ELSE
         IQ(KQSP+LCDIR+LK+5)=IRNEXT
      ENDIF
      IDTIME=0
      CALL RZDATE(IDTIME,IDATE,ITIME,2)
      IQ(KQSP+LCDIR+LK+6)=IDTIME
      LK=LK+7
      LF=LF+7
*
      NSDIR=NSDIR+1
      IQ(KQSP+LCDIR+KNFREE)=IQ(KQSP+LCDIR+KNFREE)-7
      IQ(KQSP+LCDIR+KNSD)=NSDIR
      IQ(KQSP+LCDIR+KLK)=LK
      IQ(KQSP+LCDIR+KLF)=LF
*
*            Create the subdirectory bank
*
      CALL MZBOOK(JQPDVS,LSDIR,LCDIR,-1,'RZ  ',2,1,LREC,2,0)
*
      CALL UCOPY(IHDIR,IQ(KQSP+LSDIR+1),4)
      KTAGS=KKDES+(NWKEY-1)/10+1
      LDS=KTAGS+2*NWKEY
      LBS=LDS+2
      IQ(KQSP+LSDIR+KUP)    = IQ(KQSP+LCDIR+LD+1)
      IQ(KQSP+LSDIR+KPW1)   = IHPWD(1)
      IQ(KQSP+LSDIR+KPW1+1) = IHPWD(2)
      IQ(KQSP+LSDIR+KPW1+2) = NCD
      CALL SBYT(NHPWD,IQ(KQSP+LSDIR+KPW1+2),6,5)
      IF(IMODEX.GT.0)THEN
         CALL SBIT1(IQ(KQSP+LSDIR+KPW1+2),12)
      ENDIF
      IQ(KQSP+LSDIR+KDATEC) = IDTIME
      IQ(KQSP+LSDIR+KDATEM) = IDTIME
      IQ(KQSP+LSDIR+KQUOTA) = IQ(KQSP+LCDIR+KQUOTA)-5
      IQ(KQSP+LSDIR+KRUSED) = IQ(KQSP+LSDIR+13)+1
      IQ(KQSP+LSDIR+KWUSED) = IQ(KQSP+LSDIR+14)+LREC
      IQ(KQSP+LSDIR+KRZVER) = KVSCYC
      IQ(KQSP+LSDIR+KIP1)   = 1
      IQ(KQSP+LSDIR+KNFREE) = LREC-LBS
      IQ(KQSP+LSDIR+KLD)    = LDS
      IQ(KQSP+LSDIR+KLB)    = LBS
      IQ(KQSP+LSDIR+KLS)    = LBS
      IQ(KQSP+LSDIR+KLK)    = LBS
      IQ(KQSP+LSDIR+KLF)    = LBS
      IQ(KQSP+LSDIR+KLC)    = LREC+1
      IQ(KQSP+LSDIR+KLE)    = LREC
      IQ(KQSP+LSDIR+KNWKEY) = NWKEY
      IQ(KQSP+LSDIR+LDS)    = 1
      IQ(KQSP+LSDIR+LDS+1)  = IRNEXT
*
      IF(CHFORM(1:1).NE.'?')THEN
         NCHF=LEN(CHFORM)
         NCH =LEN(CHTAG(1))
         IF(NCH.GT.8)NCH=8
         DO 30 I=1,NWKEY
            IF(NCH.LT.8)CALL VBLANK(IHDIR,2)
            CALL UCTOH(CHTAG(I),IHDIR,4,NCH)
            CALL UCOPY(IHDIR,IQ(KQSP+LSDIR+KTAGS+2*(I-1)),2)
            IF(I.LE.NCHF)THEN
               IF(CHFORM(I:I).EQ.'B')IFORM=1
               IF(CHFORM(I:I).EQ.'I')IFORM=2
               IF(CHFORM(I:I).EQ.' ')IFORM=2
               IF(CHFORM(I:I).EQ.'H')IFORM=3
               IF(CHFORM(I:I).EQ.'A')IFORM=4
            ENDIF
            IKDES=(I-1)/10
            IKBIT1=3*I-30*IKDES-2
            CALL SBYT(IFORM,IQ(KQSP+LSDIR+KKDES+IKDES),IKBIT1,3)
  30     CONTINUE
         CALL ZHTOI(IQ(KQSP+LSDIR+KTAGS),IQ(KQSP+LSDIR+KTAGS),2*NWKEY)
      ENDIF
*
*            Write new directory on file
*
      CALL RZIODO(LUN,LREC,IRNEXT,IQ(KQSP+LSDIR+1),2)
      IF(IQUEST(1).NE.0)THEN
        CALL SBIT0(IQ(KQSP+LCDIR),2)
        CALL MZDROP(JQPDVS,LSDIR,' ')
        GO TO 99
      ENDIF
*
      CALL SBIT1(IQ(KQSP+LSDIR),2)
      CALL SBIT1(IQ(KQSP+LSDIR),IQDROP)
*
*            Mark records effectively used
*
      CALL RZUSED(1,IRNEXT)
*
  99  RETURN
      END
