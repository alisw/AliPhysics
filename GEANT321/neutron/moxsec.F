*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:57  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/04 23/02/95  14.46.01  by  S.Giani
*-- Author :
      SUBROUTINE MOXSEC
C************************************************************
C
C  setup cross-section tables for MICAP
C
C  Called by: MORINI
C
C  INPUT: MICAP element IDs in KE  = LD(LFP11)
C         element densities in RHO = D (LFP12)
C
C  Author : C.Zeitnitz
C
C  See USER's GUIDE TO MICAP ORNL/TM-10340
C  for details and pointer description (MPOINT)
C
C************************************************************
#include "geant321/mmicap.inc"
#include "geant321/minput.inc"
#include "geant321/mpoint.inc"
#include "geant321/mconst.inc"
#include "geant321/cmagic.inc"
C
      CHARACTER*80 XSFILE
      CHARACTER*70 CCOMM
      CHARACTER*11 CGEANT
      CHARACTER*4  MARK
      CHARACTER*20 MATNAM
      INTEGER INEL(134)
      LOGICAL NOHED,XSTOP
C
C       CALCULATE THE NUMBER OF ELEMENTS (NNUC)
C       AND GENERATE THE ISOTOPE NUMBER ARRAY (IN(NMIX))
      NWTOT = 0
      DO 10 I=1,NMIX
         LD(LFP17+I-1)=0
         LD(LFP16+I-1)=0
   10 CONTINUE
C       INITIALIZE THE NUMBER OF ELEMENTS (NNUC)
      NNUC=0
      DO 30 I=1,NMIX
         IF(LD(LFP16+I-1).GT.0)GO TO 30
         NNUC=NNUC+1
         LD(LFP16+I-1)=NNUC
         DO 20 J=I+1,NMIX
            IF(LD(LFP11+I-1).NE.LD(LFP11+J-1))GO TO 20
            LD(LFP16+J-1)=NNUC
   20    CONTINUE
   30 CONTINUE
C get number of isoptopes from xsection file(s)
      LT = LTEMP
      NII = 0
   40 CONTINUE
        NUNIT = IQ(LT+NTUNIT)
        READ(NUNIT,'(I10)') NIS
        NII = NII + NIS
        IQ(LT+NTMPNI) = NIS
        LT = LQ(LT)
      IF(LT.GT.0) GOTO 40
C allocate needed memory for x-section
      NW = 2*NII+13*NNUC+2*NNR*NNUC+4*NGR*NNUC+3*NQ*NNUC+26*MEDIA + 2
      NI = NII
      NWTOT = NWTOT + NW
      CALL CHKZEB(NW,IXCONS)
      CALL MZBOOK(IXCONS,LMOX1,0,2,'MOX1',0,0,NW,0,-1)
C       SET UP THE B CONTROL BLOCK LOCATION NUMBER ARRAY ICOM(NNUC)
C LFP170 points to length of x-section data
      LFP170 = LMOX1 + 2
      LFP18=LFP170+NNUC
      LFP18A=LFP18+NII
      LFP19=LFP18A+NII
      LFP20=LFP19+NMIX
C       SET UP THE ARRAY (IREC(NII))
      CALL XSECN1(NII,D(LFP11),D(LFP16),D(LFP17),
     +                D(LFP18),D(LFP18A),D(LFP170),D(LFP19),
     +                D(LFP20),D(LFP20),INEL)
C check if all isotopes have been found in the x-section file(s)
      XSTOP = .FALSE.
      DO 50  I=1,NMIX
        IF(LD(LFP19+I-1).EQ.0) THEN
         WRITE(IOUT,10100)LD(LFP19+I-1)
10000    FORMAT(' MICAP: Could not find x-section of element ',I8)
         XSTOP = .TRUE.
        ENDIF
   50 CONTINUE
      IF(XSTOP) THEN
        PRINT '('' CALOR : Neutron x-section not found ===> STOP '')'
        STOP
      ENDIF
      LFP21=LFP20+NNUC
C store xs accuracy at LFP210 (used for thinning  in XSECN3)
      LFP210 = LFP21 + NNUC
      LFP22=LFP210+NNUC
      LFP23=LFP22+NNUC
      LFP24=LFP23+NNUC
      LFP25=LFP24+NNUC
      LFP26=LFP25+NNUC
      LFP27=LFP26+NNR*NNUC
      LFP28=LFP27+NNR*NNUC
      LFP29=LFP28+NNUC
      LFP30=LFP29+NNUC
      LFP31=LFP30+NGR*NNUC
      LFP32=LFP31+NGR*NNUC
      LFP33=LFP32+MEDIA
      LFP34=LFP33+MEDIA
      LFP35=LFP34+NNUC
      LFP36=LFP35+3*NQ*NNUC
C       CLEAR THE STORAGE LOCATIONS FOR THE DICTIONARIES, ETC.
      CALL CLEAR(D,LFP20,LFP36-1)
C       ESTABLISH THE RANDOM WALK STORAGE LOCATIONS
      LFP41=LFP36
      LFP42=LFP41+2*NNUC
      LFP45=LFP42+24*MEDIA
      LFP46=LFP45+NGR*NNUC
      NW = 0
      DO 60 INUC=1,NNUC
         NW = NW + LD(LFP170+INUC-1)
   60 CONTINUE
      NW = NW + 2
      NWTOT = NWTOT + NW
      CALL CHKZEB(NW,IXCONS)
      CALL MZBOOK(IXCONS,LMOX2,0,2,'MOX2',0,0,NW,0,-1)
      LFP43 = LMOX2 + 2
      LAST = LFP43 - 1
      MAXD = LMOX2 + NW
C       PLACE THE MICROSCOPIC CROSS SECTION DATA INTO THE CORE
      CALL XSECN2(D(LFP17),D(LFP18),D(LFP18A),
     +            D(LFP20),D(LFP21),D(LFP210),D(LFP22),D(LFP23),
     +            D(LFP24),D(LFP25),D(LFP26),D(LFP27),D(LFP28),
     +            D(LFP29),D(LFP30),D(LFP31),D(LFP34),D(LFP35),
     +            D(LFP35+NQ*NNUC),D(LFP35+2*NQ*NNUC),
     +            D(LFP43),D(LFP43),MAXD,LAST,INEL)
C determine length needed for macroscopic xs and mixing
      NW = 0
      DO 90  IM=1,MEDIA
         NM = 0
         LZ = 0
         DO 70 IN=1,NMIX
            IF(LD(LFP10+IN-1).NE.IM) GOTO 70
            NM = NM+1
            II = LD(LFP16+IN-1)
            LZ = MAX0(LD(LFP27+NNR*(II-1)),LZ)
C           LZ = MAX0(LDICT(1,II),LZ)
   70    CONTINUE
         IF(NM.GT.1) LZ = 4*LZ
         DO 80 J=1,NMIX
            IF(LD(LFP10+J-1).NE.IM) GOTO 80
            II = LD(LFP16+J-1)
   80    CONTINUE
         NW = NW + LZ
   90 CONTINUE
      NW = NW + 2
      NWTOT = NWTOT + NW
      CALL CHKZEB(NW,IXCONS)
      CALL MZBOOK(IXCONS,LMOX3,0,2,'MOX3',0,0,NW,0,-1)
      LAST = LMOX3 + 1
      LFP44=LAST+1
      MAXD = LMOX3+NW
C       SET, MIX AND THIN THE TOTAL CROSS SECTIONS
C       ACCORDING TO THE MIXING TABLE
      CALL XSECN3(D(LFP10),D(LFP11),D(LFP12),D(LFP16),D(LFP26),
     +            D(LFP27),D(LFP32),D(LFP33),D(LFP44),D(LFP44),
     +            D,MAXD,LAST)
C       ESTABLISH THE PHOTON TOTAL CROSS SECTION DATA DICTIONARY
C       STORAGE LOCATIONS
C determine number of words needed for photon production xs
      NW = 0
      DO 110 I=1,NNUC
         DO 100 J=1,LD(LFP28+I-1)
            LZ = LD(LFP31 + 2*J - 1 + NGR*(I-1))
C           LZ = LGCB(2*J,I)
            NW = NW + LZ
  100    CONTINUE
  110 CONTINUE
      NW = NW + 2*NGR*NNUC+2
      NWTOT = NWTOT + NW + 1
      CALL CHKZEB(NW,IXCONS)
      CALL MZBOOK(IXCONS,LMOX4,0,2,'MOX4',0,0,NW,0,-1)
      LMAG2 = LMOX4 + 1
      LD(LMAG2) = NMAGIC
      LFP45 = LMAG2 + 1
      LFP46 = LFP45 + NGR*NNUC
      LFP47 = LFP46 + NGR*NNUC
      LAST = LFP47 - 1
      MAXD = LMOX4 + NW
C       CLEAR THE STORAGE LOCATIONS FOR THE PHOTON DICTIONARIES
C       OF THE TOTAL PHOTON PRODUCTION CROSS SECTIONS
      CALL CLEAR(D,LFP45,LFP47-1)
C       SUM THE PHOTON PARTIAL DISTRIBUTIONS OF THE ENDF/B-V
C       FILE 12 AND FILE 13 DATA (BY MT NUMBER) AND PLACE THESE
C       MICROSCOPIC MULTIPLICITIES TIMES CROSS SECTIONS IN CORE
      CALL XSECN5(D(LFP28),D(LFP30),D(LFP31),D(LFP45),D(LFP46),
     +           D(LFP47),D(LFP47),D,D,MAXD,LAST)
C
C print out media to print unit IOUT
C      WRITE(IOUT,10000)
10100 FORMAT(23X,'MICAP Material Parameters',/,
     +       23X,'-------------------------',/)
      WRITE(IOUT,10200)
10200 FORMAT(8X,'GEANT Material Parameters',10X,
     +        6X,'MICAP Material Parameters',/,
     +       8X,25('-'),10X,6X,25('-'))
      WRITE(IOUT,10300)
10300 FORMAT(1X,'Material',16X,'No/Iso',4X,'A',5X,'Z',2X,'|',
     +        4X,'A',5X,'Z',3X,'Density',
     +        3X,'Coll.Len',/,44('-'),'+',33('-'))
      MFLAG = 0
      KMED  = 0
      NISO  = 1
      DO 130 I=0,NMIX-1
C get GEANT name of material
         MARK = '/   '
         IF(LD(LFP11+I)/1000.NE.LD(LFP13+I)) THEN
            MARK = '/  *'
            MFLAG=1
         ENDIF
         K1 = LD(LFP16+I)-1
         LS1 = LD(LFP26+NNR*K1)+LMOX2
         LEN = LD(LFP27+NNR*K1)/2
         EN = 1.E6
         CALL TBSPLT(D(LS1),EN,LEN,XSEC)
         XSEC = 1./XSEC/D(LFP12+I)
         IF(LD(LFP140+I).NE.0.) THEN
            WRITE(CGEANT,'(F6.1,I5)')  D(LFP140+I),LD(LFP13+I)
         ELSE
            WRITE(CGEANT,'(A11)')  '    -     -'
         ENDIF
         IF(KMED.NE.LD(LFP10+I)) THEN
          NISO = 1
          CALL GFMATE(LD(LGE2MO+LD(LFP10+I)),MATNAM,AA,ZZ,DENS,
     +       RADL,ABSL,UB,NW)
          NBLK = LNBLNK(MATNAM)
          DO 120 JC=NBLK+1,20
             WRITE(MATNAM(JC:JC),'(A1)') '.'
  120     CONTINUE
          WRITE(MARK(2:3),'(I2)') NISO
          WRITE(IOUT,10400) MATNAM,LD(LGE2MO+LD(LFP10+I)),MARK,
     +     CGEANT,
     +     D(LFP34+LD(LFP16+I)-1)*1.008665,
     +     LD(LFP11+I)/1000,D(LFP12+I),XSEC
          KMED = LD(LFP10+I)
         ELSE
          WRITE(MARK(2:3),'(I2)') NISO
          WRITE(IOUT,10500) LD(LGE2MO+LD(LFP10+I)),MARK,CGEANT,
     +     D(LFP34+LD(LFP16+I)-1)*1.008665,
     +     LD(LFP11+I)/1000,D(LFP12+I),XSEC
         ENDIF
10400    FORMAT(1X,A20,I6,A4,A11,'  |',F6.1,I5,1X,E11.4,1X,E9.3)
10500    FORMAT(1X,20X,I6,A4,A11,'  |',F6.1,I5,1X,E11.4,1X,E9.3)
         LD(LFP13+I) = LD(LFP11+I)/1000
         NISO = NISO + 1
  130 CONTINUE
      WRITE(IOUT,'(78(''-''),/,48X,''Density in (Atoms/barn/cm)'')')
      WRITE(IOUT,'(36X,
     +    ''Collision Length for 1 MeV neutron in (cm)'',/)')
      IF(MFLAG.EQ.1) WRITE(IOUT,'(/,
     + 15X,''*******************************************'',/,
     + 15X,''*               W A R N I N G             *'',/,
     + 15X,''*   Marked isotopes (*) not found in the  *'',/,
     + 15X,''*        cross-section file(s)            *'',/,
     + 15X,''*    Cross-sections of the isotope with   *'',/,
     + 15X,''*    the closest Z will be used instead   *'',/,
     + 15X,''*******************************************'',/)')
C which x-section files have been used?
      LT = LTEMP
      LCI = LCISO
      LC = 0
      NOHED=.TRUE.
  140 CONTINUE
C first check if x-section file has been used!
        NUNIT = IQ(LT+NTUNIT)
        DO 150 I=0,NMIX-1
           KISO = LD(LFP16+I)
           MISO = LD(LFP17+KISO-1)
           IF(NUNIT.EQ.LD(LFP18A+MISO-1)) GOTO 160
  150   CONTINUE
C unit never used !
        GOTO 190
  160   CONTINUE
C search for comments for selected isotopes
        NCOM = IQ(LCI+1)
        DO 180 J=1,NCOM
           K = (J-1)*81 + 2
           JZ = IQ(LCI+K)
           JA = IQ(LCI+K+1)
           CCOMM = ' '
           CALL UHTOC(IQ(LCI+K+2),4,CCOMM,70)
           DO 170 I=0,NMIX-1
              KISO = LD(LFP16+I)
              IA = NINT(D(LFP34+KISO-1)*1.008665)
              IZ = LD(LFP11+I)/1000
              MISO = LD(LFP17+KISO-1)
C print the comment, if the isotope is correct and has been read from
C the current x-section file
              IF(IA.EQ.JA .AND. IZ.EQ.JZ .AND.
     +           NUNIT.EQ.LD(LFP18A+MISO-1)) THEN
                IF(NOHED) THEN
                 WRITE(IOUT,'(/,23X,''COMMENTS ABOUT ISOTOPE DATA'')')
                 WRITE(IOUT,'(  23X,''---------------------------'',/)')
                 NOHED = .FALSE.
                ENDIF
                LC = LC + 1
                WRITE(IOUT,'(I4,'') '',A70)') LC,CCOMM
                GOTO 180
              ENDIF
  170      CONTINUE
  180   CONTINUE
  190   LT = LQ(LT)
        LCI = LQ(LCI)
      IF(LT.GT.0.AND.LCI.GT.0) GOTO 140
C print the x-section file names and comments
      WRITE(IOUT,'(/,20X,''USED NEUTRON CROSS-SECTION FILES'')')
      WRITE(IOUT,'(  20X,''--------------------------------'',/)')
      LT = LTEMP
  200 CONTINUE
C first check if x-section file has been used!
        NUNIT = IQ(LT+NTUNIT)
        DO 210 I=0,NMIX-1
           KISO = LD(LFP16+I)
           MISO = LD(LFP17+KISO-1)
           IF(NUNIT.EQ.LD(LFP18A+MISO-1)) GOTO 220
  210   CONTINUE
C unit never used !
        GOTO 230
  220   CONTINUE
C get file name of x-section file
        XSFILE = ' '
        COMMEN = ' '
        DATSTR = ' '
        CALL UHTOC(IQ(LT+NTNAME+1),4,XSFILE,IQ(LT+NTNAME))
        CALL UHTOC(IQ(LT+NTCOMM+1),4,COMMEN,IQ(LT+NTCOMM))
        CALL UHTOC(IQ(LT+NTDATS+1),4,DATSTR,IQ(LT+NTDATS))
        WRITE(IOUT,'('' File      : '',A66)') XSFILE
        WRITE(IOUT,'('' Generated : '',A24,/,
     +               '' Comment   : '',A66,/)') DATSTR,COMMEN
  230   LT = LQ(LT)
      IF(LT.GT.0) GOTO 200
      WRITE(IOUT,'(/,'' MICAP :'',I10,
     + '' words used in GCBANK for neutron x-section tables''/)') NWTOT
      RETURN
      END
