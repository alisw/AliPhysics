*
* $Id$
*
* $Log$
* Revision 1.4  1997/11/24 14:47:21  jamie
* set IOPTRR in /RZCOUT/
*
* Revision 1.3  1996/04/24 17:27:21  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.2  1996/03/08 08:08:06  jamie
* Bug fixes for opt R handling in rz(v)out
*
* Revision 1.1.1.1  1996/03/06 10:47:27  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZVOUT(V,N,KEYU,ICYCLE,CHOPT)
*
************************************************************************
*
*           User FORTRAN array V of length N is output on a RZ file
* Input:
*   V       Fortran array of length N
*   KEYU    Keyword vector of length NWKEY as specified by RZMDIR.
*   ICYCLE  only used as Input parameter if 'A' option is given (see below)
*   CHOPT   Character variable specifying the selected options.
*           mode
*             default
*                   Array V contains only floating points
*             'I'   Array V contains integers
*             'B'   Array V contains bit-patterns
*             'H'   Array V contains Holleriths
*             'D'   Array V contains Double precision words
*             'A'   Key will not be visible by RZLDIR
*             'R'   Replace option.
*             'S'   Used for sequential operation. Application
*                   guarantees that keys are unique and objects
*                   are simply appended without checking all existing keys
* Output:
*   ICYCLE  Cycle number associated to the key entered
*           ICYCLE is 1 if KEY was  not already present in the directory,
*           and one larger than the previous  cycle associated to the key
*           otherwise.
*           ICYCLE is only Input parameter when option 'A' is given
*
* Called by <USER>
*
*  Author  : R.Brun DD/US/PD
*  Written : 03.04.86
*  Last mod: 25.06.93 - JDS Add IEVENT to RZWRT call (see RZWRT)
*          : 04.03.94 - S.Banerjee (Change in cycle structure)
*          : 05.09.94 - Add option S
*          : 23.03.95 J.Shiers - key # in cycles block is KEY(1)
*          : 10.04.95 J.Shiers - move IKYV lines to after definition
*          : 08.03.96 J.Shiers - use ICOLD (not 1) if opt R specified
*          : 24.11.97 M.Brun - set IOPTRR in /RZCOUT/
*
************************************************************************
#include "zebra/zunit.inc"
#include "zebra/rzcl.inc"
#include "zebra/rzclun.inc"
#include "zebra/rzcout.inc"
#include "zebra/rzk.inc"
#include "zebra/rzckey.inc"
#include "zebra/rzcycle.inc"
#include "zebra/mzioc.inc"
      CHARACTER*(*) CHOPT
      DIMENSION KEYU(*),V(*)
      DIMENSION    IOPTV(7)
      EQUIVALENCE (IOPTA,IOPTV(1)), (IOPTD,IOPTV(2))
     +,           (IOPTI,IOPTV(3)), (IOPTB,IOPTV(4))
     +,           (IOPTH,IOPTV(5)), (IOPTR,IOPTV(6))
     +,           (IOPTS,IOPTV(7))
*
*-----------------------------------------------------------------------
*
#include "zebra/q_jbyt.inc"
*
      IQUEST(1)=0
      IEVENT   =0
*
*     Loglevel
*
      LOGLV = JBYT(IQ(KQSP+LTOP),15,3)-3
*
      CALL UOPTC(CHOPT,'ADIBHRS',IOPTV)
*
*     Set also /RZCOUT/ common variable
*
      IOPTRR = IOPTR
*
*     Options R & S are incompatible
*
      IF(IOPTR.NE.0.AND.IOPTS.NE.0) THEN
         IF(LOGLV.GE.3) WRITE(IQPRNT,10000)
10000    FORMAT(' RZVOUT. ERROR - options R and S are incompatible')
         IQUEST(1)=5
         GOTO 999
      ENDIF
*
*         Check if WRITE permission on file and directory
*
      IF(LQRS.EQ.0)GOTO 999
      IF(N.LE.0)THEN
         IQUEST(1)=4
         GOTO 999
      ENDIF
      IFLAG=0
      CALL RZMODS('RZVOUT',IFLAG)
      IF(IFLAG.NE.0)GOTO 999
*
*         Write current buffer if not same directory
*         Get last record written in that directory
*         Create buffer bank
*
      LROUT  = LQ(KQSP+LTOP-6)
      IROUT  = IQ(KQSP+LTOP+KIROUT)
      IRLOUT = IQ(KQSP+LCDIR+KRLOUT)
      IP1    = IQ(KQSP+LCDIR+KIP1)
      IF(LROUT.EQ.0)THEN
         CALL MZBOOK(JQPDVS,LROUT,LTOP,-6,'RZOU',0,0,LREC+1,2,-1)
         IQ(KQSP+LROUT-5)=LUN
         IROUT=0
         IP1=1
      ENDIF
      IF(IROUT.NE.IRLOUT.AND.IRLOUT.NE.0)THEN
         CALL RZIODO(LUN,LREC,IRLOUT,IQ(KQSP+LROUT+1),1)
#if defined(CERNLIB_QMVAX)
         IF(IRELAT.NE.0)UNLOCK(UNIT=LUN)
#endif
         IROUT=IRLOUT
         IQ(KQSP+LTOP+KIROUT)=IROUT
         IP1=IQ(KQSP+LCDIR+KIP1)
         IF(IQ(KQSP+LTOP+KIRIN).EQ.IROUT)IQ(KQSP+LTOP+KIRIN) = 0
      ENDIF
*
*         Is directory big enough to accomodate new cycle ?
*
      IF(IOPTR.EQ.0) THEN
         IF(IQ(KQSP+LCDIR+KNFREE).LT.IQ(KQSP+LCDIR+KNWKEY)+4*KLCYCL+1)
     +                                                              THEN
            CALL RZEXPD('RZVOUT',10*(IQ(KQSP+LCDIR+KNWKEY)+KLCYCL+1))
            IF(IQUEST(1).NE.0)GOTO 999
         ENDIF
      ENDIF
 
      LK = IQ(KQSP+LCDIR+KLK)
      LF = IQ(KQSP+LCDIR+KLF)
      LC = IQ(KQSP+LCDIR+KLC)
      NWFREE=IQ(KQSP+LCDIR+KNFREE)
 
      IF(IOPTR.EQ.0) THEN
*
*         Compute how many records and how many words
*         are necessary to write array V.
*
         IR1=IRLOUT
         IALLOC=0
         NDATA=N
         NLEFT=LREC-IP1+1
         IF(NDATA.LE.NLEFT)THEN
            NR=0
         ELSE
            NR=(NDATA-NLEFT-1)/LREC + 1
         ENDIF
         IF(IRLOUT.EQ.0)NR=NR+1
         IF(NR.GT.0)THEN
            CALL RZALLO('RZVOUT',NR,IALLOC)
            IF(IALLOC.EQ.0)GOTO 999
         ENDIF
         IF(IRLOUT.EQ.0.OR.IP1.GT.LREC)THEN
            IP1=1
            NLEFT=LREC
            IRLOUT=IALLOC
            IR1=IALLOC
            IR2=IALLOC+1
            IR3=IALLOC+NR-1
            IF(NR.EQ.1)THEN
               IR2=0
               IR3=0
            ENDIF
         ELSE
            IR2=IALLOC
            IR3=IALLOC+NR-1
         ENDIF
 
      ENDIF
 
      NKEYS = IQ(KQSP+LCDIR+KNKEYS)
      NWKEY = IQ(KQSP+LCDIR+KNWKEY)
      IQUEST(7)=NKEYS
      IQUEST(8)=NWKEY
*
*         Convert input key vector to internal format
*
         DO 10 I=1,NWKEY
            IKDES=(I-1)/10
            IKBIT1=3*I-30*IKDES-2
            IF(JBYT(IQ(KQSP+LCDIR+KKDES+IKDES),IKBIT1,3).LT.3)THEN
               KEY(I)=KEYU(I)
            ELSE
               CALL ZHTOI(KEYU(I),KEY(I),1)
            ENDIF
   10    CONTINUE
*
*         Search if KEY is already entered
*
      IF(IOPTS.EQ.0) THEN
         IF(NKEYS.GT.0)THEN
            DO 30 I=1,NKEYS
               DO 20 K=1,NWKEY
                  LKC=LK+(NWKEY+1)*(I-1)
                  IF(IQ(KQSP+LCDIR+LKC+K).NE.KEY(K))GOTO 30
   20          CONTINUE
               LCOLD=IQ(KQSP+LCDIR+LKC)
               IF (KVSCYC.NE.0) THEN
*                 IF (IQ(KQSP+LCDIR+LCOLD+KKYCYC).NE.I) THEN
                  IF (IQ(KQSP+LCDIR+LCOLD+KKYCYC).NE.
     +                IQ(KQSP+LCDIR+LKC+1)) THEN
                     IQUEST(1) = 11
                     GO TO 999
                  ENDIF
               ENDIF
               ICOLD = JBYT(IQ(KQSP+LCDIR+LCOLD+KCNCYC),21,12)
*              IKYV = I
               IKYV = IQ(KQSP+LCDIR+LKC+1)
 
               IF(IOPTR.EQ.0) THEN
                  ICYCLE=ICOLD+1
               ELSE
*                 ICYCLE=1
                  ICYCLE=ICOLD
                  IF (KVSCYC.EQ.0) THEN
                     IP1 = JBYT(IQ(KQSP+LCDIR+LCOLD+KORCYC), 1,16)
                     IR1 = JBYT(IQ(KQSP+LCDIR+LCOLD+KFRCYC),17,16)
                     IR2 = JBYT(IQ(KQSP+LCDIR+LCOLD+KSRCYC),17,16)
                     NWORDS = JBYT(IQ(KQSP+LCDIR+LCOLD+KNWCYC), 1,20)
                  ELSE
                     IP1 = JBYT(IQ(KQSP+LCDIR+LCOLD+KORCYC), 1,20)
                     IR1 = IQ(KQSP+LCDIR+LCOLD+KFRCYC)
                     IR2 = IQ(KQSP+LCDIR+LCOLD+KSRCYC)
                     NWORDS = IQ(KQSP+LCDIR+LCOLD+KNWCYC)
                  ENDIF
                  NLEFT=LREC-IP1+1
                  IF(NWORDS.LE.NLEFT)THEN
                     NR=0
                  ELSE
                     NR=(NWORDS-NLEFT-1)/LREC + 1
                  ENDIF
                  IF(LOGLV.GE.3) WRITE(IQPRNT,10100) IP1,IR1,NWORDS,NR
10100          FORMAT(' RZVOUT. object starts at word ',I6,
     +                ' in record ',I6,' nwords = ',I6,' nrecs = ',I6)
                  IF(LOGLV.GE.3.AND.IR2.NE.0) WRITE(IQPRNT,10200) IR2
10200          FORMAT(' RZVOUT. object continues in record ',I6)
                  IRSAVE = IRLOUT
                  IF(IR1.NE.IRLOUT) THEN
                     CALL RZIODO(LUN,LREC,IR1,IQ(KQSP+LROUT+1),1)
                     IF(IQUEST(1).NE.0)GOTO 999
                     IROUT = IR1
                     IRLOUT = IR1
                  ENDIF
               ENDIF
 
               IQUEST(20)=I
               GOTO 50
   30       CONTINUE
         ENDIF
*
*     Object must already exist if R option is given
*
         IF(IOPTR.NE.0) THEN
            IQUEST(1) = 6
            IF(LOGLV.GE.-2) WRITE(IQLOG,10300)
10300    FORMAT(' RZVOUT. Error - object does not exist')
            GOTO 999
         ENDIF
 
      ENDIF
*
*         New KEY, append to the list
*
      IQUEST(20)=NWKEY+1
      NWFREE=NWFREE-NWKEY-1
      IQ(KQSP+LCDIR+KNKEYS)=IQ(KQSP+LCDIR+KNKEYS)+1
      LCOLD = 0
      LKC   = LF
      LF    = LF+NWKEY+1
      ICYCLE= 1
      DO 40 I=1,NWKEY
         IQ(KQSP+LCDIR+LKC+I)=KEY(I)
   40 CONTINUE
*     IKYV  = IQ(KQSP+LCDIR+KNKEYS)
      IKYV  = IQ(KQSP+LCDIR+LKC+1)
*
*         Write user array
*
   50 IQ3=IRLOUT
      IQ4=IP1
 
#if !defined(CERNLIB_FQXISN)
      IF(IMODEX.GT.0)THEN
         MFO(1)= 3
         IF(IOPTI.NE.0)MFO(1)=2
         IF(IOPTB.NE.0)MFO(1)=1
         IF(IOPTH.NE.0)MFO(1)=5
         MFO(2)=-1
         JFOEND= 2
      ENDIF
#endif
      CALL RZWRT(V,NDATA,IOPTB,IEVENT)
      IF(IQUEST(1).NE.0)THEN
         IF(ICYCLE.EQ.1.AND.IOPTR.EQ.0)THEN
            IQ(KQSP+LCDIR+KNKEYS)=IQ(KQSP+LCDIR+KNKEYS)-1
         ENDIF
         IROUT=IQ(KQSP+LTOP+KIROUT)
         IF(IROUT.GT.0)THEN
            IQ1=IQUEST(1)
            CALL RZIODO(LUN,LREC,IROUT,IQ(KQSP+LROUT+1),1)
#if defined(CERNLIB_QMVAX)
            IF(IRELAT.NE.0)UNLOCK(UNIT=LUN)
#endif
            IQUEST(1)=IQ1
         ENDIF
         GOTO 999
      ENDIF
 
      IF(IOPTR.EQ.0) THEN
         IF(IRLOUT.EQ.IR1) GOTO 60
         IF(IRLOUT.GE.IR2.AND.IRLOUT.LE.IR3) GOTO 60
         IF(IROUT.EQ.IRLOUT)IROUT=0
         IRLOUT = 0
         IP1 = 1
   60    CONTINUE
*
*         Create a new cycle
*
         LC    = LC-KLCYCL
         NWFREE= NWFREE-KLCYCL
         IFORM = 3
         IF(IOPTI.NE.0) IFORM = 2
         IF(IOPTB.NE.0) IFORM = 1
         IF(IOPTH.NE.0) IFORM = 5
         IQ(KQSP+LCDIR+LKC)       = LC
         IQ(KQSP+LCDIR+LC+KPPCYC) = LCOLD
         IQ(KQSP+LCDIR+LC+KFLCYC) = IFORM
         CALL RZDATE(IQ(KQSP+LCDIR+LC+KFLCYC),IDATE,ITIME,2)
         IF(IOPTA.NE.0) CALL SBIT1(IQ(KQSP+LCDIR+LC+KFLCYC),4)
         IQ(KQSP+LCDIR+LC+KORCYC) = IQ4
         IQ(KQSP+LCDIR+LC+KNWCYC) = NDATA
         CALL SBYT(ICYCLE,IQ(KQSP+LCDIR+LC+KCNCYC),21,12)
         IF (KVSCYC.EQ.0) THEN
            IF (NLEFT.LT.NDATA)
     +         CALL SBYT(IR2,IQ(KQSP+LCDIR+LC+KSRCYC),17,16)
            CALL SBYT(IR1,IQ(KQSP+LCDIR+LC+KFRCYC),17,16)
         ELSE
            IF (NLEFT.LT.NDATA) THEN
               IQ(KQSP+LCDIR+LC+KSRCYC) = IR2
            ELSE
               IQ(KQSP+LCDIR+LC+KSRCYC) = 0
            ENDIF
            IQ(KQSP+LCDIR+LC+KFRCYC) = IR1
            IQ(KQSP+LCDIR+LC+KKYCYC) = IKYV
         ENDIF
      ENDIF
*
*         Update internal pointers in the directory
*
      IQUEST(3)=IQ3
      IQUEST(4)=IQ4
      IQUEST(5)=0
      IQUEST(6)=ICYCLE
      IQUEST(11)=NDATA
 
      IF(IOPTR.EQ.0) THEN
 
         IQ(KQSP+LTOP+KIROUT)=IROUT
         IQ(KQSP+LCDIR+KRUSED)=IQ(KQSP+LCDIR+KRUSED)+NR
         NWUSED=IQ(KQSP+LCDIR+KWUSED)+NDATA
         IF(NWUSED.GT.1000000)THEN
            IQ(KQSP+LCDIR+KMEGA)=IQ(KQSP+LCDIR+KMEGA)+1
            IQ(KQSP+LCDIR+KWUSED)=NWUSED-1000000
         ELSE
            IQ(KQSP+LCDIR+KWUSED)=NWUSED
         ENDIF
         IQ(KQSP+LCDIR+KRLOUT)=IRLOUT
         IQ(KQSP+LCDIR+KIP1)=IP1
         IQ(KQSP+LCDIR+KNFREE)=NWFREE
         IQ(KQSP+LCDIR+KLF)=LF
         IQ(KQSP+LCDIR+KLC)=LC
 
      ENDIF
 
      IQUEST(9)=IQ(KQSP+LCDIR+KQUOTA)-IQ(KQSP+LCDIR+KRUSED)
*
*         Mark used records
*
      IF(IOPTR.EQ.0.AND.NR.GT.0)THEN
         CALL RZUSED(NR,IALLOC)
      ENDIF
*
*     Restore last record written if required
*
      IF(IOPTR.NE.0.AND.IRSAVE.NE.IRLOUT) THEN
         CALL RZIODO(LUN,LREC,IRLOUT,IQ(KQSP+LROUT+1),2)
         IF(IQUEST(1).NE.0)GOTO 999
         IROUT  = IRSAVE
         IRLOUT = IRSAVE
         CALL RZIODO(LUN,LREC,IRLOUT,IQ(KQSP+LROUT+1),1)
         IF(IQUEST(1).NE.0)GOTO 999
      ENDIF
*
  999 END
