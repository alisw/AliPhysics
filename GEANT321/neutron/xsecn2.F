*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:22:00  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/04 23/02/95  14.46.01  by  S.Giani
*-- Author :
      SUBROUTINE XSECN2(ICOM,IREC,IUNIT,IGAMS,LGAM,ELTOL,INABS,LNAB,
     + ITHRMS,LTHRM,IDICTS,LDICT,NTX,NTS,IGCBS,LGCB,AWR,Q,LR,QLR,
     + BUF,IBUF,LIM,LAST,INEL)
C       THIS ROUTINE READS THE REMAINDER OF INPUT I/O UNIT(s),
C       SELECTS THE ELEMENTS NEEDED FOR THE CALCULATIONS,
C       AND STORES THE CROSS SECTION DATA IN CORE
#include "geant321/minput.inc"
#include "geant321/mconst.inc"
#include "geant321/mmicab.inc"
      CHARACTER*4 MARK
      DIMENSION BUF(*),IBUF(*),ICOM(*),IGAMS(*),LGAM(*),INABS(*),
     +LNAB(*),ITHRMS(*),LTHRM(*),AWR(*),IDICTS(NNR,NNUC),ELTOL(*),
     +LDICT(NNR,NNUC),Q(NQ,NNUC),NTX(*),NTS(*),IGCBS(NGR,NNUC),
     +LGCB(NGR,NNUC),IREC(*),LR(NQ,NNUC),QLR(NQ,NNUC)
      DIMENSION INEL(*),IUNIT(*)
C       ASSIGN THE DEFAULT VALUES
      LEN=0
C       INITIALIZE THE COUNTERS FOR THE LOOP
C       NISR EQUALS THE NUMBER OF ISOTOPES READ
C       IRECNO EQUALS THE NEXT RECORD NUMBER TO BE READ ON INPUT
C  I/O UNIT (NUNIT)
C       LAST EQUALS THE LAST CORE POSITION USED IN THE CALLING
CROUTINE (INPUT1)
C       LST EQUALS THE LAST POSITION USED IN THE BUF ARRAY
C   (I.E. (BUF(LST) = D(LAST)))
      NISR=0
      IRECNO=1
      LST=0
C       PRINT OUT THE CROSS SECTION DIRECTORY IF CALLED FOR
   10 CONTINUE
C       START LOOP TO READ IN THE DATA ON INPUT I/O UNIT
      DO 370 II=1,NI
         IR   = IREC(II)
         IF(NUNIT.NE.IUNIT(II)) IRECNO = 1
         NUNIT= IUNIT(II)
         IF(NUNIT.LE.0) THEN
           WRITE(IOUT,'(/,'' XSECN2 : Wrong unit number '',I10)') NUNIT
           GOTO 370
         ENDIF
         IF(NISR.GE.NNUC)GO TO 370
         IF(IR.EQ.0)GO TO 370
C       LOOP TO LOCATE THE I CONTROL BLOCK RECORD (IR=IREC(II))
CZ x-section endmark = 'ENDE'
CZ file endmark ='ENDF'
         MARK = '   '
   20    IF(MARK.EQ.'ENDE') IRECNO = IRECNO + 1
         IF(MARK.EQ.'ENDF') GOTO 50
         IF(IR.EQ.IRECNO) GOTO 30
         READ(NUNIT,'(A)') MARK
         GO TO 20
C       CHECK TO DETERMINE THE ISOTOPE NUMBER FOR THE RANDOM WALK
   30    DO 40 I=1,NNUC
            IF(ICOM(I).EQ.II)GO TO 60
   40    CONTINUE
   50    WRITE(IOUT,10000)II
10000 FORMAT('0',10X,'ERROR IN ROUTINE XSECN2, II=',I6,/)
         GO TO 390
C       READ I CONTROL BLOCK RECORD OFF INPUT I/O UNIT (NUNIT) FOR
C       THE ELEMENT CORRESPONDING TO IREC(II) AND ICOM(I)
   60    IJK=I
         READ(NUNIT,'(I10,4G13.7,1I10,/,6I10)') IBUF(LST+1),(BUF(LST+
     +   IK),IK=2,5),(IBUF(LST+IJ),IJ=6,12)
         NISR=NISR+1
C       ASSIGN VALUES TO ARRAYS NEEDED FOR THE RANDOM WALK
         ISO=IJK
         NEL=INEL(II)
         AWR(ISO)=BUF(LST+2)
CZ store accuracy of xs
         ELTOL(ISO) = BUF(LST+4)
         IFLAGU=IBUF(LST+6)
         LGAM(ISO)=IBUF(LST+7)
         NTX(ISO)=IBUF(LST+8)
         NTS(ISO)=IBUF(LST+9)
         LTHRM(ISO)=IBUF(LST+11)
         LNAB(ISO)=IBUF(LST+12)
C       READ IN THE ISOTOPE DICTIONARY (IDICT ARRAY)
C       FROM INPUT I/O UNIT (NUNIT)
         READ(NUNIT,'((8I10))')(LDICT(J,ISO),J=1,NNR)
   70    CONTINUE
C       READ IN ENDF/B FILE3 CROSS SECTION DATA
C       READ IN ENDF/B FILE4 ANGULAR DISTRIBUTION DATA
C       READ IN ENDF/B FILE5 SECONDARY ENERGY DISTRIBUTION DATA
         DO 190 I2=1,NNR
            LZ=LDICT(I2,ISO)
            IF(LZ.EQ.0)GO TO 190
            LEN=LIM-LAST
            IF(LEN.LT.LZ)GO TO 380
            IDICTS(I2,ISO)=LAST+1-LMOX2
CZ changed in order to read ASCII input file
C I2 < 67  -> x-section data
C I2 < 123 -> angular distribution
C I2 < 134 -> secondary energy distribution
C I2 = 134 ->
            IF(I2.LT.67) THEN
               READ(NUNIT,'((6G13.7))')(BUF(LST+I),I=1,LZ)
            ELSE IF(I2.LT.123) THEN
C ------------------- I2 = 67 -----------------------------
               READ(NUNIT,'((8I10))') (IBUF(LST+I),I=1,2), (IBUF(LST+
     +         J+2),J=1,2*IBUF(LST+1))
               K = 2*IBUF(LST+1) + 2 + 1
               DO 80 J=1,IBUF(LST+2)
                  READ(NUNIT,'(G13.7,I10,/,(6G13.7))') BUF(LST+K),
     +            IBUF(LST+K+1), (BUF(LST+IK+K+1),IK=1,IBUF(LST+K+1)*2)
                  K = K + 2 + IBUF(LST+K+1)*2
   80          CONTINUE
            ELSE IF(I2.LT.134) THEN
C-------------------- I2 = 123 ----------------------------
               READ(NUNIT,'(2I10,G13.7,2I10,/,(8I10))') (IBUF(LST+I),
     +         I=1,2),BUF(LST+3),(IBUF(LST+J),J=4,5), (IBUF(LST+K+5),K=
     +         1,2*IBUF(LST+4))
               ID = 2*IBUF(LST+4) + 5
               LF = IBUF(LST+2)
               NP2 = 2*IBUF(LST+5)
               READ(NUNIT,'((6G13.7))') (BUF(LST+ID+I),I=1,NP2)
               ID = ID + NP2
               KEND = 1
               IF(LF.EQ.5.OR.LF.EQ.11) KEND = 2
               DO 100 K=1,KEND
                  READ(NUNIT,'((8I10))') (IBUF(LST+ID+I),I=1,2)
                  NR2 = 2*IBUF(LST+ID+1)
                  NE = IBUF(LST+ID+2)
                  ID = ID + 2
                  READ(NUNIT,'((8I10))') (IBUF(LST+ID+I),I=1,NR2)
                  ID = ID + NR2
                  IEND = NE
                  IF(LF.EQ.5.OR.LF.EQ.11) IEND = 1
                  IF(LF.EQ.7.OR.LF.EQ.9) IEND = 1
                  DO 90 I=1,IEND
                     IF(LF.EQ.1) THEN
                        READ(NUNIT,'(G13.7,2I10)') BUF(LST+ID+1),
     +                  (IBUF(LST+ID+J),J=2,3)
                        NR2 = 2*IBUF(LST+ID+2)
                        NP2 = 2*IBUF(LST+ID+3)
                        ID = ID + 3
                        READ(NUNIT,'((8I10))') (IBUF(LST+ID+J),J=1,
     +                  NR2)
                        ID = ID + NR2
                     ELSE
                        NP2 = 2*NE
                     ENDIF
                     READ(NUNIT,'((6G13.7))') (BUF(LST+ID+J),J=1,NP2)
                     ID = ID + NP2
   90             CONTINUE
  100          CONTINUE
            ELSE
C ------------------ I2 = 134 --------------------------------------
               READ(NUNIT,'(I10)') IBUF(LST+1)
               LNU = IBUF(LST+1)
               IF(LNU.NE.2) THEN
                  READ(NUNIT,'(I10,/,(6G13.7))') IBUF(LST+2), (BUF(LST
     +            +I+2),I=1,IBUF(LST+2))
               ELSE
                  READ(NUNIT,'((8I10))') (IBUF(LST+I),I=2,3)
                  NR2 = IBUF(LST+2)*2
                  READ(NUNIT,'((8I10))') (IBUF(LST+3+J),J=1,NR2)
                  NP2 = IBUF(LST+3)*2
                  READ(NUNIT,'((6G13.7))') (BUF(LST+3+NR2+J),J=1,NP2)
               ENDIF
            ENDIF
CZ end of change
            IF(I2.GT.66)GO TO 120
  110       CONTINUE
            GO TO 180
  120       IF(I2.GT.122)GO TO 150
  130       CONTINUE
            CALL ANGCDF(BUF(LST+1),BUF(LST+1),LZ)
  140       CONTINUE
            GO TO 180
  150       IF(I2.GT.133)GO TO 170
  160       CONTINUE
            GO TO 180
  170       CONTINUE
  180       CONTINUE
            LAST=LAST+LZ
            LST=LST+LZ
  190    CONTINUE
C       READ IN THE AVERAGE PHOTON PRODUCTION ARRAY
         LZ=LGAM(ISO)
         IF(LZ.EQ.0)GO TO 210
         LEN=LIM-LAST
         IF(LEN.LT.LZ)GO TO 380
         IGAMS(ISO)=LAST+1-LMOX2
         READ(NUNIT,'((6G13.7))')(BUF(LST+I),I=1,LZ)
  200    CONTINUE
         LAST=LAST+LZ
         LST=LST+LZ
  210    CONTINUE
C       READ IN THE TOTAL NEUTRON DISAPPERANCE ARRAY
         LZ=LNAB(ISO)
         IF(LZ.EQ.0)GO TO 230
         LEN=LIM-LAST
         IF(LEN.LT.LZ)GO TO 380
         INABS(ISO)=LAST+1-LMOX2
         READ(NUNIT,'((6G13.7))')(BUF(LST+I),I=1,LZ)
  220    CONTINUE
         LAST=LAST+LZ
         LST=LST+LZ
  230    CONTINUE
C       READ IN THE Q VALUE ARRAY
         READ(NUNIT,'((6G13.7))')(Q(I,ISO),I=1,NQ)
  240    CONTINUE
C       READ IN THE LR VALUE ARRAY
         READ(NUNIT,'((8I10))')(LR(I,ISO),I=1,NQ)
  250    CONTINUE
C       READ IN THE QLR VALUE ARRAY
         READ(NUNIT,'((6G13.7))')(QLR(I,ISO),I=1,NQ)
  260    CONTINUE
C       READ IN THE PHOTON DATA DICTIONARY (GCB ARRAY)
C       FROM INPUT I/O UNIT (NUNIT)
C       CURRENT STORAGE IS SET TO ACCOMODATE UP TO 30 INTERACTIONS
C       (I.E. (2*NTX(ISO)+2*NTS(ISO)).LE.NGR)
         L=2*NTX(ISO)+2*NTS(ISO)
         IF(L.EQ.0)GO TO 350
         L1=2*NTX(ISO)
         L2=L1+1
         READ(NUNIT,'((8I10))')(LGCB(J,ISO),J=1,L)
  270    CONTINUE
C       READ IN ENDF/B FILE12 PHOTON MULTIPLICATION DATA
C       READ IN ENDF/B FILE13 PHOTON CROSS SECTION DATA
         NNTX=NTX(ISO)
         DO 300 I2=1,NNTX
            LZ=LGCB(2*I2,ISO)
            IF(LZ.EQ.0)GO TO 300
            LEN=LIM-LAST
            IF(LEN.LT.LZ)GO TO 380
            IGCBS(2*I2-1,ISO)=LGCB(2*I2-1,ISO)
            IGCBS(2*I2,ISO)=LAST+1-LMOX2
CZ changed in order to read ASCII xsection file
            READ(NUNIT,'((8I10))') (IBUF(LST+I),I=1,2)
            READ(NUNIT,'((6G13.7))') (BUF(LST+J+2),J=1,IBUF(LST+2))
            ID = IBUF(LST+2) + 2 + LST
            DO 280 K = 1, IBUF(LST+1)
               READ(NUNIT,'(2(G13.7,I10))') BUF(ID+1),IBUF(ID+2),
     +         BUF(ID+3),IBUF(ID+4)
               ID = ID + 4
               READ(NUNIT,'((6G13.7))') (BUF(ID + J),J=1,IBUF(LST+2))
               ID = ID + IBUF(LST+2)
  280       CONTINUE
CZ end of change
  290       CONTINUE
            LAST=LAST+LZ
            LST=LST+LZ
  300    CONTINUE
C       READ IN ENDF/B FILE15 PHOTON SECONDARY ENERGY DISTRIBUTIONS
         NNTS=NTS(ISO)
         IF(NNTS.EQ.0)GO TO 350
         DO 340 I2=1,NNTS
            LZ=LGCB(L1+2*I2,ISO)
            IF(LZ.EQ.0)GO TO 340
            LEN=LIM-LAST
            IF(LEN.LT.LZ)GO TO 380
            IGCBS(L1+2*I2-1,ISO)=LGCB(L1+2*I2-1,ISO)
            IGCBS(L1+2*I2,ISO)=LAST+1-LMOX2
CZ changed in order to read ASCII xsection file
            READ(NUNIT,'(2I10,G13.7,2I10,/,(8I10))') (IBUF(LST+I),I=1,
     +      2),BUF(LST+3), (IBUF(LST+J),J=4,5), (IBUF(LST+K+5),K=1,2*
     +      IBUF(LST+4))
            ID = 2*IBUF(LST+4) + 5
            LF = IBUF(LST+2)
            NP2 = 2*IBUF(LST+5)
            READ(NUNIT,'((6G13.7))') (BUF(LST+ID+I),I=1,NP2)
            ID = ID + NP2
            KEND = 1
            IF(LF.EQ.5.OR.LF.EQ.11) KEND = 2
            DO 320 K=1,KEND
               READ(NUNIT,'((8I10))') (IBUF(LST+ID+I),I=1,2)
               NR2 = 2*IBUF(LST+ID+1)
               NE = IBUF(LST+ID+2)
               ID = ID + 2
               READ(NUNIT,'((8I10))') (IBUF(LST+ID+I),I=1,NR2)
               ID = ID + NR2
               IEND = NE
               IF(LF.EQ.5.OR.LF.EQ.11) IEND = 1
               IF(LF.EQ.7.OR.LF.EQ.9) IEND = 1
               DO 310 I=1,IEND
                  IF(LF.EQ.1) THEN
                     READ(NUNIT,'(G13.7,2I10)') BUF(LST+ID+1), (IBUF(L
     +               ST+ID+J),J=2,3)
                     NR2 = 2*IBUF(LST+ID+2)
                     NP2 = 2*IBUF(LST+ID+3)
                     ID = ID + 3
                     READ(NUNIT,'((8I10))') (IBUF(LST+ID+J),J=1,NR2)
                     ID = ID + NR2
                  ELSE
                     NP2 = 2*NE
                  ENDIF
                  READ(NUNIT,'((6G13.7))') (BUF(LST+ID+J),J=1,NP2)
                  ID = ID + NP2
  310          CONTINUE
  320       CONTINUE
CZ end of change
  330       CONTINUE
            LAST=LAST+LZ
            LST=LST+LZ
  340    CONTINUE
  350    CONTINUE
C       READ IN THE THERMAL CROSS SECTION DATA ARRAY
         LZ=LTHRM(ISO)
         IF(LZ.EQ.0)GO TO 360
         LEN=LIM-LAST
         IF(LEN.LT.LZ)GO TO 380
         ITHRMS(ISO)=LAST+1
         READ(NUNIT,'((6G13.7))')(BUF(LST+I),I=1,LZ)
         LAST=LAST+LZ
         LST=LST+LZ
  360    CONTINUE
  370 CONTINUE
      GO TO 400
  380 WRITE(IOUT,10100)LZ,LEN
10100 FORMAT('0','NOT ENOUGH SPACE TO READ IN RECORD',/,5X,
     +'LENGTH OF RECORD=',I10,/,5X,'SPACE AVAILABLE=',I10)
  390 PRINT '('' CALOR: ERROR in XSECN2 ====> STOP '')'
      STOP
  400 RETURN
      END
