*
* $Id$
*
* $Log$
* Revision 1.3  1996/04/24 17:27:03  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.2  1996/03/08 08:08:04  jamie
* Bug fixes for opt R handling in rz(v)out
*
* Revision 1.1.1.1  1996/03/06 10:47:26  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZOUT(IXDIV,LSUP,KEYU,ICYCLE,CHOPT)
*
************************************************************************
*
*       To write data structure pointed by LSUP
* Input:
*   IXDIV   Index of the division(s)
*           May be zero if the 'D' option is not selected
*           May be a compound  index (see MZIXCO on page 24)   if the 'D'
*           option is selected
*   LSUP    Supporting address of the data structure  (may be zero if the
*           'D' option is selected)
*   KEYU    Keyword vector of length NWKEY as specified by RZMDIR.
*   ICYCLE  only used as Input parameter if 'A' option is given (see below)
*   CHOPT   Character variable specifying the selected options.
*           data structure
*             default
*                   The data structure  supported by the bank  at LSUP is
*                   written out (link 0 is not followed)
*             'D'   Complete division(s)
*                   default:    Dropped banks are squeezed out
*                         (slower but maybe more economic than 'DI')
*             'DI'  Immediate  dump  of  divisions   with  dropped  banks
*                   included
*             'L'   Write  the data  structure  supported  by the  linear
*                   structure at LSUP (link 0 is followed)
*             'S'   Single bank at LSUP
*             'R'   Replace option.
*             'Q'   seQuential mode - no check made to see if key already
*                   exists (option S in RZVOUT etc.)
*           mode
*             default
*                   Keep banks available after output
*             'N'   No links, i.e. linkless handling
*             'W'   Drop data structure or wipe division(s) after output
*             'A'   Key will not be visible by RZLDIR
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
*  Written : 04.04.86
*  Last mod: 14.09.93 - RB. Protect against directories > 64K
*          : 04.03.94 S.Banerjee (Change in cycle structure)
*          : 27.09.94 J.Shiers - add option Q (cf S in RZVOUT)
*          : 17.02.95 J.Shiers - move definition of logl
*          : 21.02.95 J.Shiers - improve warning message for
*                                big directories (objects auto-deleted)
*          : 23.03.95 J.Shiers - key # in cycles block is KEY(1)
*          : 10.04.95 J.Shiers - move IKYV lines to after definition
*          : 08.03.96 J.Shiers - only delete previous cycle if opt R
*                                not specified. If opt R is specified,
*                                use previous cycle number (and not 1)
*
*  In RZ file format version 0, RZ pointers are stored in 16 bits.
*  Thus, pointers to previous cycle etc. must not exceed 65536.
*  An attempt to save a new cycle under such conditions will result
*  in the previous cycle being deleted with a warning message.
*
*  For files created with RZ version 1, this is not necessary as
*  pointers are 32 bit.
*
************************************************************************
#include "zebra/zunit.inc"
#include "zebra/rzcl.inc"
#include "zebra/rzclun.inc"
#include "zebra/eqlqrrf.inc"
#include "zebra/rzcout.inc"
#include "zebra/rzk.inc"
#include "zebra/rzckey.inc"
#include "zebra/rzcycle.inc"
#include "zebra/mzct.inc"
#include "zebra/fzcx.inc"
      CHARACTER*(*) CHOPT
      DIMENSION KEYU(*)
      DIMENSION LSUP(1)
      DIMENSION    IOPTV(9)
*
*-----------------------------------------------------------------------
*
#include "zebra/q_jbyt.inc"
*
      IQUEST(1)=0
      CALL MZSDIV(IXDIV,1)
      IXDIVX=IXDIV
      LQRRF(KQT+1)=LSUP(1)
      IEVFLX=1
      NWTXX =0
      NWSEGX=0
      NWTABX=0
      NWBKX =0
      NWUHCX=0
      CALL UOPTC(CHOPT,'ADILSNWRQ',IOPTV)
      IOPTXA=IOPTV(1)
      IOPTXD=IOPTV(2)
      IOPTXI=IOPTV(3)
      IOPTXL=IOPTV(4)
      IOPTXS=IOPTV(5)
      IOPTXN=IOPTV(6)
      IOPTXW=IOPTV(7)
      IOPTXR=IOPTV(8)
      IOPTRR=IOPTV(8)
      IOPTXQ=IOPTV(9)
      IOPTXN=IOPTXN+IOPTXS
      IOPTXM=0
*
*     Loglevel
*
      LOGLV = JBYT(IQ(KQSP+LTOP),15,3)-3
*
*     Options R & Q are incompatible
*
      IF(IOPTXR.NE.0.AND.IOPTXQ.NE.0) THEN
         IF(LOGLV.GE.3) WRITE(IQPRNT,10000)
10000    FORMAT(' RZOUT. ERROR - options R and Q are incompatible')
         IQUEST(1)=5
         GOTO 999
      ENDIF
*
*         Check if WRITE permission on file and directory
*
      IF(LQRS.EQ.0)GOTO 999
      IFLAG=0
      CALL RZMODS('RZOUT ',IFLAG)
      IF(IFLAG.NE.0) GOTO 999
*
*         Write current buffer if not same directory
*         Get last record written in that directory
*         Create buffer bank
*
   10 LROUT  = LQ(KQSP+LTOP-6)
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
         IF(IQUEST(1).NE.0)GOTO 999
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
      NKEYS = IQ(KQSP+LCDIR+KNKEYS)
      NWKEY = IQ(KQSP+LCDIR+KNWKEY)
 
      IF(IOPTXR.EQ.0) THEN
         IF(IQ(KQSP+LCDIR+KNFREE).LT.NWKEY+4*KLCYCL+1)THEN
            CALL RZEXPD('RZOUT ',10*(NWKEY+KLCYCL+1))
            IF(IQUEST(1).NE.0)GOTO 999
         ENDIF
      ENDIF
 
      LK = IQ(KQSP+LCDIR+KLK)
      LF = IQ(KQSP+LCDIR+KLF)
      LC = IQ(KQSP+LCDIR+KLC)
      LE = IQ(KQSP+LCDIR+KLE)
      NWFREE=IQ(KQSP+LCDIR+KNFREE)
*
*         Convert input key vector to internal format
*
      IQUEST(7)=NKEYS
      IQUEST(8)=NWKEY
*
      DO 20 I=1,NWKEY
         IKDES=(I-1)/10
         IKBIT1=3*I-30*IKDES-2
         IF(JBYT(IQ(KQSP+LCDIR+KKDES+IKDES),IKBIT1,3).LT.3)THEN
            KEY(I)=KEYU(I)
         ELSE
            CALL ZHTOI(KEYU(I),KEY(I),1)
         ENDIF
   20 CONTINUE
*
*         Search if KEY is already entered
*
      IF(IOPTXQ.EQ.0) THEN
         IF(NKEYS.GT.0)THEN
            DO 40 I=1,NKEYS
               DO 30 K=1,NWKEY
                  LKC=LK+(NWKEY+1)*(I-1)
                  IF(IQ(KQSP+LCDIR+LKC+K).NE.KEY(K))GOTO 40
   30          CONTINUE
*
*     Protect against directories > 65536
*
               IF(IOPTXR.EQ.0.AND.KVSCYC.EQ.0.AND.LE.GT.65536)THEN
                  CALL RZDELK(KEYU,0,'C')
                  IF(LOGLV.GE.-2) THEN
                     WRITE(IQPRNT,10100)
10100 FORMAT(' RZOUT: current RZ file cannot support > 64K records ',/,
     +       '        or individual directories > 64K')
                     WRITE(IQPRNT,10200)KEY(1)
10200 FORMAT(' RZOUT: previous cycle(s) for this key (',I8,')',
     +       '  deleted')
                     WRITE(IQPRNT,10300)
10300 FORMAT(' RZOUT: please consult ZEBRA manual for further details')
                  ENDIF
                  GOTO 10
               ENDIF
 
               LCOLD=IQ(KQSP+LCDIR+LKC)
               IF (KVSCYC.NE.0) THEN
*                 IF (IQ(KQSP+LCDIR+LCOLD+KKYCYC).NE.I) THEN
                  LKC=LK+(NWKEY+1)*(I-1)
                  IF (IQ(KQSP+LCDIR+LCOLD+KKYCYC).NE.
     +                IQ(KQSP+LCDIR+LKC+1)) THEN
                     IQUEST(1) = 11
                     GO TO 999
                  ENDIF
               ENDIF
               ICOLD = JBYT(IQ(KQSP+LCDIR+LCOLD+KCNCYC),21,12)
*              IKYV = I
               IKYV = IQ(KQSP+LCDIR+LKC+1)
               IF(IOPTXR.EQ.0) THEN
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
                  IF(LOGLV.GE.3) WRITE(IQPRNT,10400) IP1,IR1,NWORDS,NR
10400          FORMAT(' RZOUT. object starts at word ',I6,
     +                ' in record ',I6,' nwords = ',I6,' nrecs = ',I6)
                  IF(LOGLV.GE.3.AND.IR2.NE.0) WRITE(IQPRNT,10500) IR2
10500          FORMAT(' RZOUT. object continues in record ',I6)
                  IRSAVE = IRLOUT
                  IF(IR1.NE.IRLOUT) THEN
                     CALL RZIODO(LUN,LREC,IR1,IQ(KQSP+LROUT+1),1)
                     IF(IQUEST(1).NE.0)GOTO 999
                     IROUT = IR1
                     IRLOUT = IR1
                  ENDIF
               ENDIF
               IQUEST(20)=I
               GOTO 60
   40       CONTINUE
         ENDIF
*
*     Object must already exist if R option is given
*
         IF(IOPTXR.NE.0) THEN
            IQUEST(1) = 6
            IF(LOGLV.GE.-2) WRITE(IQLOG,10600)
10600    FORMAT(' RZOUT. Error - object does not exist')
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
      DO 50 I=1,NWKEY
         IQ(KQSP+LCDIR+LKC+I)=KEY(I)
   50 CONTINUE
*     IKYV  = IQ(KQSP+LCDIR+KNKEYS)
      IKYV  = IQ(KQSP+LCDIR+LKC+1)
*
*         Construct table of material to be written
*
   60 CALL MZSDIV(IXDIVX,0)
      LENTRX = LQRRF(KQT+1)
      JQSTMV = -1
      MODTBX = 0
      JFLGAX = 0
      CALL FZOTAB
      IF(IQUEST(1).NE.0)THEN
         IF(LOGLV.GE.-2) WRITE(IQLOG,10700)
10700    FORMAT(' RZOUT. Unable to construct table')
         GOTO 80
      ENDIF
      IF(IOPTXN.NE.0)NWTABX=0
      NDATA=NWTABX+NWBKX+3
*
*     Replace option - record allocation not required
*
      IF(IOPTXR.NE.0) THEN
*
*     Check if exactly the same number of words are required
*
         IF(NDATA.NE.NWORDS) THEN
            IQUEST(1) = 1
            IF(LOGLV.GE.-2) WRITE(IQLOG,10800) NWORDS,NDATA
10800    FORMAT(' RZOUT. Error - existing object required ',I10,
     +          ' words. New object requires ',I10)
            GOTO 90
         ENDIF
 
      ELSE
*
*         Compute how many records and how many words
*         are necessary to write data structure.
*
         IR1=IRLOUT
         IALLOC=0
         NLEFT=LREC-IP1+1
         IF(NDATA.LE.NLEFT)THEN
            NR=0
         ELSE
            NR=(NDATA-NLEFT-1)/LREC + 1
         ENDIF
         IF(IRLOUT.EQ.0)NR=NR+1
         IF(NR.GT.0)THEN
            CALL RZALLO('RZOUT ',NR,IALLOC)
            IF(IALLOC.EQ.0)GOTO 80
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
      IQ3 =IRLOUT
      IQ4 =IP1
*
*         Write data structure at LSUP according to table
*
      CALL RZOBKN
      IF(IQUEST(1).NE.0)THEN
         IROUT=IQ(KQSP+LTOP+KIROUT)
         IF(IROUT.GT.0)THEN
            IQ1=IQUEST(1)
            CALL RZIODO(LUN,LREC,IROUT,IQ(KQSP+LROUT+1),1)
#if defined(CERNLIB_QMVAX)
            IF(IRELAT.NE.0)UNLOCK(UNIT=LUN)
#endif
            IQUEST(1)=IQ1
         ENDIF
         GOTO 80
      ENDIF
 
      IF(IOPTXR.EQ.0) THEN
 
         IF(IRLOUT.EQ.IR1) GOTO 70
         IF(IRLOUT.GE.IR2.AND.IRLOUT.LE.IR3) GOTO 70
         IF(IROUT.EQ.IRLOUT)IROUT=0
         IRLOUT = 0
         IP1 = 1
   70    CONTINUE
*
*         Create a new cycle
*
         LC     = LC-KLCYCL
         NWFREE = NWFREE-KLCYCL
         IQ(KQSP+LCDIR+LKC)       = LC
         IQ(KQSP+LCDIR+LC+KPPCYC) = LCOLD
         IQ(KQSP+LCDIR+LC+KFLCYC) = 0
         CALL RZDATE (IQ(KQSP+LCDIR+LC+KFLCYC),IDATE,ITIME,2)
         IF(IOPTXA.NE.0)CALL SBIT1(IQ(KQSP+LCDIR+LC+KFLCYC),4)
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
 
      IF(IOPTXR.EQ.0) THEN
 
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
*         Drop data structure or wipe division (option W)
*
      IF (IOPTXW.NE.0) THEN
         LSUP(1) = LQRRF(KQT+1)
         IF (LSUP(1).NE.0) THEN
            CALL MZDROP(IXDIV,LSUP,' ')
            LSUP(1)=0
         ENDIF
      ENDIF
*
*         Mark used records
*
      IF(IOPTXR.EQ.0.AND.NR.GT.0)THEN
         CALL RZUSED(NR,IALLOC)
      ENDIF
      GOTO 90
   80 IF(ICYCLE.EQ.1.AND.IOPTXR.EQ.0)THEN
         IQ(KQSP+LCDIR+KNKEYS)=IQ(KQSP+LCDIR+KNKEYS)-1
      ENDIF
 
   90 CONTINUE
*
*     Restore last record written if required
*
      IF(IOPTXR.NE.0.AND.IRSAVE.NE.IRLOUT) THEN
         CALL RZIODO(LUN,LREC,IRLOUT,IQ(KQSP+LROUT+1),2)
         IF(IQUEST(1).NE.0)GOTO 999
         IROUT  = IRSAVE
         IRLOUT = IRSAVE
         CALL RZIODO(LUN,LREC,IRLOUT,IQ(KQSP+LROUT+1),1)
         IF(IQUEST(1).NE.0)GOTO 999
      ENDIF
*
  999 END
