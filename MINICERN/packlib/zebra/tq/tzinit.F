*
* $Id$
*
* $Log$
* Revision 1.3  1996/04/24 17:27:35  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.2  1996/04/18 16:14:53  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:27  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE TZINIT (LUNP,IXDIV)

C--   Master routine for title input

#include "zebra/zmach.inc"
#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
#include "zebra/mqsys.inc"
#include "zebra/eqlqt.inc"
#include "zebra/mzcwk.inc"
#include "zebra/tzuc.inc"
#include "zebra/tzc1.inc"
C--------------    END CDE                             -----------------
      COMMON /SLATE/ NDSLAT,NESLAT,NFSLAT,NGSLAT,DUMMY(36)
      CHARACTER    COL(LGL)*1
      EQUIVALENCE (COL(1), LINE(1:1))
      CHARACTER    CHWORK*8
#if defined(CERNLIB_EQUHOLCH)
      INTEGER      IWORK(2)
      EQUIVALENCE (CHWORK,IWORK)
#endif

      DIMENSION    LUNP(1)
      DIMENSION    IFLAGS(5)
      EQUIVALENCE (IFLAGS(1),IFLLOG)

      PARAMETER    (NCC = 7)
      CHARACTER    CCTEXT(NCC)*6

#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HTZIN, 4HIT   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HTZINIT /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'TZINIT  ')
#endif
      DATA CCTEXT  / 'LOG   ', 'PR*INT', 'US*ER ', 'KI*LL '
     +,              'ANY*WA', 'DO    ', 'FIN*IS' /

#include "zebra/q_sbyt.inc"


#include "zebra/qtraceq.inc"

C--       decide division

      CALL MZSDIV (IXDIV,7)
      L = LQT(KQT+1)
      IF (L.NE.0)  JQDIVI = MZFDIV (-7, L)
      IF (JQDIVI.EQ.0)  JQDIVI = 2
      JSTOR  = JQSTOR
      IXTITL = MSBYT (JQSTOR,JQDIVI,27,4)

C--       invert the order of the title banks,
C-        so that the new ones can be created at the beginning

      IF (L.NE.0)  CALL ZTOPSY (IXTITL,L)

      CALL VZERO (NPARA,23)
      CALL VZERO (NREADY,12)
      IF (NQLOGL.GE.-1)  IFLLOG = 1
      BLANK = ' '

C--       get the input unit

      LUNTQ = LUNP(1)
      IF (LUNTQ.LT.0)   LUNTQ = IQTTIN
      IF (LUNTQ.EQ.0)   LUNTQ = IQREAD
#if defined(CERNLIB_QPRINT)
      IF (IFLLOG.NE.0)  THEN
          JDV = JQDIVI
          IF (JDV.EQ.JQDVSY)  JDV = 24
          WRITE (IQLOG,9000) LUNTQ,JQSTOR,JDV
        ENDIF
 9000 FORMAT (1X/' TZINIT.  Read title banks from LUNTQ = ',I3/
     F10X,'for store',I3,' into division',I3/1X)
#endif

C----              Look for next control line

   12 IF (NREADY.NE.0)             GO TO 14

#include "zebra/tzread1.inc"

   14 NREADY = 0
      IF (NCHORG.EQ.0)             GO TO 12
      IF (LINE(1:2).EQ.'*.')       GO TO 12
      IF (LINE(1:4).EQ.'*CMZ')     GO TO 12

      NCHLN  = NCHORG
      J = INDEX (LINE(1:NCHLN),' #.')
      IF (J.NE.0)  NCHLN = LNBLNK (LINE(1:J))
      IF (NCHLN.EQ.0)              GO TO 12

      IF (COL(1).NE.'*')           GO TO 91
      IF (IFLLOG.NE.0)  WRITE (IQLOG,9012) LINE(1:NCHORG)
 9012 FORMAT ('    > ',A)

      IF (LINE(1:2).EQ.'*-')       GO TO 12

      NHEAD = MIN (NCHORG,80)
      LHEAD(1:NHEAD) = LINE(1:NHEAD)
      CALL CLTOU (LINE(1:NCHLN))
      IGNORE = 0

C----              Analyse control line

      JA  = ICNEXT (LINE,1,NCHLN) + 1
      JE  = NESLAT - 1
      NN  = NDSLAT

      IF (NN.LT.3)                 GO TO 92
      JCC = ICNTH (LINE(2:NN),CCTEXT,NCC)

C-        log   print user  kill  anyw  do    finis
C-        1     2     3     4     5     6     7

      IF (JCC.EQ.0)                GO TO 92
      IF (JCC.EQ.NCC)              GO TO 78

      JA = ICNEXT (LINE,JE+1,NCHLN)
      JE = NESLAT - 1
      NN = NDSLAT
      IF (JCC.EQ.NCC-1)            GO TO 21

C--          *LOG, *PRINT, *USER, *KILL, *ANYWAY (OFF)

      IF (NN.EQ.0)  THEN
          IFLAGS(JCC) = 1
          IF (JCC.EQ.2)  IFLLOG = 1
          GO TO 12
        ENDIF

      IF (NN.NE.3)                 GO TO 92
      CHWORK = LINE(JA:JA+2)
      IF (CHWORK(1:3).NE.'OFF')    GO TO 92
      IFLAGS(JCC) = 0
      IF (JCC.EQ.1)  IFLPRI = 0
      GO TO 12

C--------          Analyse  *DO idh idn -opt x       --------

   21 IF (NN.EQ.0)                 GO TO 92
      NN = MIN(NN,4)
      CHWORK = ' '
      CHWORK(1:NN) = LINE(JA:JA+NN-1)
      CALL VZERO (NAME,5)
      CALL VZERO (IFLEXA,16)
#if defined(CERNLIB_EQUHOLCH)
      NAME(1) = IWORK(1)
#endif
#if !defined(CERNLIB_EQUHOLCH)
      CALL UCTOH (CHWORK,NAME,4,4)
#endif
      IDNUM   = -1

   24 JA = ICNEXT (LINE,JE+1,NCHLN)
      IF (JA.GT.NCHLN)             GO TO 61
      JE = NESLAT - 1
      NN = NDSLAT

      IF (COL(JA).NE.'-')  THEN
          IF (IDNUM.GE.0)          GO TO 92
          IDNUM = ICDECI (LINE,JA,JE)
          IF (IDNUM.LT.0)          GO TO 92
          IF (NGSLAT.NE.0)         GO TO 92
          GO TO 24
        ENDIF

      IDNUM = MAX (IDNUM,0)
      IF (NN.EQ.1)                 GO TO 92

      JA = JA + 1
      JOPT = INDEX ('FISENUCA',COL(JA))
      IF (JOPT.EQ.0)               GO TO 92

      JA = JA + 1
      NN = NN - 2
      GO TO (31, 34, 36, 37, 39, 41, 44, 47), JOPT
C-            F   I   S   E   N   U   C   A
C-            1   2   3   4   5   6   7   8

C--       option  -F(format), Fortran FORMAT

   31 IF (COL(JA).NE.'(')          GO TO 92
      IF (COL(JE).NE.')')          GO TO 92
      JFMTC1 = JA
      JFMTC2 = JE
      GO TO 24

C--       option -Ii or -If ... or -I(text), I/O characteristic

   34 JOPT = INDEX ('(BIFD', COL(JA))
C-                   12345
      IF (JOPT.EQ.0)               GO TO 92
      IF (JOPT.NE.1)  THEN
          IF (NN.NE.1)             GO TO 92
          NAME(5) = JOPT - 1
          GO TO 24
        ENDIF

      JE = ICFIND (')', LINE,JA,NCHLN)
      IF (NGSLAT.EQ.0)             GO TO 92

      CALL MZIOBK (NAME,20, LINE(JA+1:JE-1))
      GO TO 24

C--       option -S[n], true size of the bank

   36 IFLSIZ = 1
      IF (NN.EQ.0)                 GO TO 24
      GO TO 39

C--       option -E[n], exact number of data words

   37 IFLEXA = 1
      IF (NN.EQ.0)                 GO TO 24

C--       option -Nn, bank size

   39 N  = ICDECI (LINE,JA,JE)
      IF (NGSLAT.NE.0)             GO TO 92
      IF (N.LE.0)                  GO TO 92
      IF (NAME(4).NE.0)  THEN
          IF (NAME(4).NE.N)        GO TO 94
        ENDIF
      NAME(4) = N
      GO TO 24

C--       option -U[n], call TZUSER

   41 IFLTZU = 1
      IF (NN.EQ.0)                 GO TO 24
      IVALUS = ICDECI (LINE,JA,JE)
      IF (NGSLAT.NE.0)             GO TO 92
      GO TO 24

C--       option -C[a][/e], column usage

   44 JCOLA = 1
      IF (COL(JA).EQ.'/')          GO TO 45
      JCOLA = ICDECI (LINE,JA,JE)
      IF (JCOLA.LE.0)              GO TO 92
      IF (JCOLA.GE.LGL-3)          GO TO 92
      IF (NGSLAT.EQ.0)             GO TO 24
      JA = NESLAT
      IF (COL(JA).NE.'/')          GO TO 92
   45 JCOLE = ICDECI (LINE,JA+1,JE)
      IF (NGSLAT.NE.0)             GO TO 92
      JCOLE = MIN (JCOLE,LGL)
      IF (JCOLE.LT.JCOLA+2)        GO TO 92
      GO TO 24

C--       option -A[n][C][W]

   47 CALL TZACW (LINE(JA-1:JE))
      IF (NCHPW.GE.0)              GO TO 24
      NCHPW = 0
      GO TO 92

C--------          Read data of next title bank      --------

   61 IDNUM = MAX (IDNUM,0)
      IF (NAME(4).NE.0)            GO TO 63

C--                read the data through the scratch area

      IF (JFMTC1.NE.0)             GO TO 93
      LPUTA = LQWKFZ
      LPUTE = LPUTA + NQWKTT
      GO TO 66

C--                pre-lift the bank for filling

   63 IFLPRE = 7
      CALL MZLIFT (IXTITL,LOLD,LQT(KQT+1),1,NAME,-1)
      LPUTA = KQS + LOLD + 9
      LPUTE = LPUTA + NAME(4)
      IF (JFMTC1.EQ.0)             GO TO 66

C--                Read data en block with Fortran FORMAT

      NWOCC = NAME(4)
#include "zebra/tzread2.inc"
      GO TO 68

C--                Read data in free-field format

   66 CALL TZFREE
      IF (NFAULT.NE.0)  THEN
          NFATAL = NFATAL + NFAULT
          IF (IFLPRE.EQ.0)         GO TO 77
          LNEW = -1
          GO TO 69
        ENDIF
      NWOCC = LPUTX - LPUTA

C--       lift the bank and copy the data

      IF (IFLPRE.EQ.0)  THEN
          NAME(4) = NWOCC
          CALL MZLIFT (IXTITL,LOLD,LQT(KQT+1),1,NAME,-1)
          CALL UCOPY (LQ(LPUTA),IQ(KQS+LOLD+1),NAME(4))
        ENDIF

C--       handle 'exact' or 'size'

      IF (IFLEXA.NE.0)  THEN
          IF (NWOCC.NE.NAME(4))    GO TO 97
        ENDIF

      IF (IFLSIZ.EQ.0)             GO TO 68
      N = NAME(4) - NWOCC
      IF (N.NE.0)  CALL VZERO (IQ(KQS+LOLD+1+NWOCC),N)
      NWOCC = NAME(4)

C--       calling TZUSER

   68 IQ(KQS+LOLD-5) = IDNUM
      IF (IFLTZU+IFLUSE.EQ.0)      GO TO 71
      NPARA = IVALUS

      LNEW = 0
      CALL TZUSER (LQT(KQT+1))
      IF (LNEW.LT.-1)              GO TO 98
      IF (LNEW.EQ.0)               GO TO 71
      NWOCC = NAME(4)
   69 CALL MZDROP (IXTITL, LQT(KQT+1), '.')
      IF (LNEW.EQ.-1)              GO TO 77

C--       possibly reduce the bank

   71 N = NAME(4) - NWOCC
      IF (N.GT.0)  THEN
          L = LQT(KQT+1)
          CALL MZPUSH (IXTITL,L,0,-N,'I')
        ENDIF

   77 IF (NREADY.GE.0)             GO TO 12

C-------           Finished

   78 L = LQT(KQT+1)
      CALL ZTOPSY (IXTITL,L)
      IQUEST(1) = NFATAL
      IF (NFATAL.NE.0)  THEN
          IF (IFLANY.EQ.0)  CALL ZFATAM ('TZINIT fails.')
        ENDIF
#include "zebra/qtrace99.inc"
      RETURN

C----              Error handling

C--       unheaded data line
   91 IF (IGNORE.NE.0)             GO TO 12
      WRITE (IQLOG,9091) LINE(1:NCHORG)
 9091 FORMAT (' !!i> ',A)
      IF (IFLKIL.GT.0)  NFATAL = NFATAL + 1
      GO TO 12

C--       faulty control line *xxx
   92 IF (IFLLOG.EQ.0)  WRITE (IQLOG,9012) LHEAD(1:NHEAD)
      WRITE (IQLOG,9092) BLANK(1:JA)
 9092 FORMAT (' !!f ',A,'^->  !!! fault')
      IGNORE = 7
      NFATAL = NFATAL + 1
      GO TO 12

C--       faults with *DO
   93 IF (IFLLOG.EQ.0)  WRITE (IQLOG,9012) LHEAD(1:NHEAD)
      WRITE (IQLOG,9093)
 9093 FORMAT (' !!f   fault :  -F(...) requires -Nn')
      IGNORE = 7
      NFATAL = NFATAL + 1
      GO TO 12

   94 IF (IFLLOG.EQ.0)  WRITE (IQLOG,9012) LHEAD(1:NHEAD)
      WRITE (IQLOG,9094)
 9094 FORMAT (' !!f   fault :  contradictory bank sizes')
      IGNORE = 7
      NFATAL = NFATAL + 1
      GO TO 12

C--       premature EoF
   96 IF (IFLLOG.EQ.0)  WRITE (IQLOG,9012) LHEAD(1:NHEAD)
      WRITE (IQLOG,9096)
 9096 FORMAT (' !!f   fault :  premature EoF')
      NFATAL = NFATAL + 1
      GO TO 78

C--       exact fault
   97 IF (IFLLOG.EQ.0)  WRITE (IQLOG,9012) LHEAD(1:NHEAD)
      WRITE (IQLOG,9097) NAME(4),NWOCC
 9097 FORMAT (' !!f   fault :  # of words expected / read=',2I6)
      NFATAL = NFATAL + 1
      GO TO 77

C--       user kill
   98 IF (IFLLOG.EQ.0)  WRITE (IQLOG,9012) LHEAD(1:NHEAD)
      WRITE (IQLOG,9098)
 9098 FORMAT (' !!f   fault :  kill requested by TZUSER')
      NFATAL = NFATAL + 1
      GO TO 77
      END
*      ==================================================
#include "zebra/qcardl.inc"
