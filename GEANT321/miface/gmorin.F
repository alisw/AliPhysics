*
* $Id$
*
* $Log$
* Revision 1.3  1996/09/30 14:26:05  ravndal
* Windows NT related modifications
*
* Revision 1.2  1996/04/26 12:21:12  ravndal
* NAP character*20 declaration included
*
* Revision 1.1.1.1  1995/10/24 10:21:52  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/04 23/02/95  14.53.54  by  S.Giani
*-- Author :    Christian Zeitnitz   21/07/92
      SUBROUTINE GMORIN
C**************************************************************
C                  Initialize MICAP
C                  ================
C Called by : CALINI
C
C Purpose : setup cross-section tables and initialize pointer
C           print flags etc.
C
C Author : C.Zeitnitz
C
C last modification: Changed in order to read new x-section file
C
C
C for details see MICAP manual ORNL/TM-10340
C*************************************************************
C MICAP commons
#include "geant321/mmicap.inc"
#include "geant321/mpoint.inc"
#include "geant321/minput.inc"
#include "geant321/mmass.inc"
#include "geant321/mconst.inc"
#include "geant321/cmagic.inc"
#include "geant321/cerrcm.inc"
#include "geant321/camass.inc"
C GEANT common
#include "geant321/gccuts.inc"
#include "geant321/gcflag.inc"
#include "geant321/gcunit.inc"
C pointer to material/mixture bank (NMATE,JMATE)
#include "geant321/gcnum.inc"
C
      COMMON / QUEST / IQUEST(100)
C
C  Avogadro number multiplied by 1.E-24
      PARAMETER(XNAVO = 0.60221367)
C
      DIMENSION A(100),AGEA(100),Z(100),DEN(100),MID(100,2),IDI(20,2)
      DIMENSION IPID(0:11)
      CHARACTER*100 XSFILE
      CHARACTER*4   CNAME
      CHARACTER*20  NAP
      CHARACTER*70 CCOMM
#if defined(CERNLIB_UNIX)||defined(CERNLIB_CRAY)||defined(CERNLIB_VAX)
      CHARACTER*100 CHROOT
#endif
      LOGICAL OPENED,EXISTS,FMIST,FSINGL,FMIFL
C GEANT Particle IDs used to extract masses from GEANT
      DATA IPID /14 , 13 , 8 , 7 , 9 , 5 , 6 , 45 , 46 , 49 , 47 , 1/
C
C   Initialization flag of GEANT
C   Loop for XMASS extracted from GCALOR routine CALINI  K.L-P
C
      IFINIT(7) = 1
      DO 5 I=0,11
         CALL GFPART(IPID(I),NAP,ITR,AM,CH,TL,UB,NW)
         XMASS(I)=AM
    5 CONTINUE
C
C  neutron energy cut (eV)
      ECUT = CUTNEU * 1.E9
C get time cut off from GEANT
      TCUT = TOFMAX
C temperature for thermal neutron xsection (Kelvin)
C only temporary constant !!!
      TEMP = 300.0
C xsection file unit
      MICROS = 31
      IOUT = LOUT
      INN  = LIN
C open MICAP I/O units
      INQUIRE(UNIT=MICROS,OPENED=OPENED)
      IF(OPENED) THEN
         REWIND MICROS
      ELSE
#if defined(CERNLIB_VAX)
         XSFILE='xsneut95.dat'
         INQUIRE(FILE=XSFILE,EXIST=EXISTS)
         IF(.NOT.EXISTS) THEN
            ISTAT = LIB$SYS_TRNLOG('CERN_ROOT',NALL,CHROOT,,,%VAL(0))
            IF(ISTAT.EQ.1) XSFILE = 'CERN_ROOT:[LIB]xsneut95.dat'
         ENDIF
         INQUIRE(FILE=XSFILE,EXIST=EXISTS)
         IF(.NOT.EXISTS) THEN
           PRINT*,'**********************************'
           PRINT*,'*        G C A L O R             *'
           PRINT*,'*        -----------             *'
           PRINT*,'*  File XSNEUT95.DAT not found   *'
           PRINT*,'*         Program STOP           *'
           PRINT*,'*   Check CERN_ROOT environment  *'
           PRINT*,'*           variable             *'
           PRINT*,'**********************************'
           STOP
         ENDIF
         OPEN(UNIT=MICROS,FILE=XSFILE, STATUS='OLD',READONLY)
#endif
#if defined(CERNLIB_UNIX)||defined(CERNLIB_CRAY)
         XSFILE = 'xsneut95.dat'
         INQUIRE(FILE=XSFILE,EXIST=EXISTS)
         IF(.NOT.EXISTS) THEN
            CHROOT=' '
            CALL GETENV('CERN_ROOT',CHROOT)
            LNROOT = LNBLNK(CHROOT)
            IF(LNROOT.GT.0)
     +      XSFILE = CHROOT(1:LNROOT)//'/lib/xsneut95.dat'
         ENDIF
         INQUIRE(FILE=XSFILE,EXIST=EXISTS)
         IF(.NOT.EXISTS) THEN
           PRINT*,'**********************************'
           PRINT*,'*        G C A L O R             *'
           PRINT*,'*        -----------             *'
           PRINT*,'*  File XSNEUT95.DAT not found   *'
           PRINT*,'*         Program STOP           *'
           PRINT*,'*   Check CERN_ROOT environment  *'
           PRINT*,'*           variable             *'
           PRINT*,'**********************************'
           STOP
         ENDIF
         OPEN(UNIT=MICROS,FILE=XSFILE,STATUS='OLD')
#endif
#if defined(CERNLIB_MSDOS)||defined(CERNLIB_WINNT)
         CHROOT=' '
         CALL GETENVF('CERN_ROOT',CHROOT)
         LNROOT = LNBLNK(CHROOT)
         IF(LNROOT.LE.0) THEN
            XSFILE = 'xsneut95.dat'
         ELSE
            XSFILE = CHROOT(1:LNROOT)//'\\lib\\xsneut95.dat'
            INQUIRE(FILE=XSFILE,EXIST=EXISTS)
            IF(.NOT.EXISTS) XSFILE='xsneut95.dat'
         ENDIF
         OPEN(UNIT=MICROS,FILE=XSFILE)
#endif
#if defined(CERNLIB_IBMVM)
         XSFILE = '\XSNEUT95 DAT *'
         OPEN(UNIT=MICROS,FILE=XSFILE,STATUS='OLD')
#endif
      ENDIF
C setup the link areas needed for x-section banks
      CALL MZLINK(IXCONS,'MICTMP',LTEMP,LTEMP,LTEMP)
      CALL MZLINK(IXCONS,'MMICAP',LMAG2,LMOX4,LMAG2)
      CALL MZLINK(IXCONS,'MPOINT',LMAG1,LFP210,LMAG1)
C
      LSUP = 0
      LCSUP = 0
      NUNIT = MICROS
C pointers into TEMP bank
      NTUNIT = 1
      NTNAME = NTUNIT + 1
      NTMPNI = NTNAME + 1 + 80/4
      NTCOMM = NTMPNI + 1
      NTDATS = NTCOMM + 1 + 80/4
      NTLIST = NTDATS + 1 + 24/4
   10 CONTINUE
C read comment and date of xsection file
        READ(NUNIT,'(A80,/,A24)') COMMEN,DATSTR
C read in material definition array
        READ(NUNIT,'(I10)') NISO
        NWW = NISO * 3 + 12 + NTLIST
C get temporary buffer
        CALL CHKZEB(NWW,IXCONS)
        IF(LSUP.EQ.0) Then
C create a top level bank for the list of isotopes
          CALL MZBOOK(IXCONS,LTEMP,LSUP,1,'TEMP',3,0,NWW,0,-1)
          LT = LTEMP
        ELSE
C create an additional bank in the linear structure TEMP
          CALL MZBOOK(IXCONS,LT,LSUP,0,'TEMP',3,0,NWW,0,-1)
          LSUP = LT
        ENDIF
        NREC = NISO * 3 / 12
        NN = 0
C store the unit number of the file in bank TEMP
        IQ(LT + NTUNIT) = NUNIT
C store the file name in bank TEMP
        CALL UCTOH(XSFILE,IQ(LT+NTNAME+1),4,LNBLNK(XSFILE))
        IQ(LT + NTNAME) = LNBLNK(XSFILE)
C store the comment and date string in bank TEMP
        IQ(LT + NTCOMM) = LNBLNK(COMMEN)
        CALL UCTOH(COMMEN,IQ(LT+NTCOMM+1),4,LNBLNK(COMMEN))
        IQ(LT + NTDATS) = LNBLNK(DATSTR)
        CALL UCTOH(DATSTR,IQ(LT+NTDATS+1),4,LNBLNK(DATSTR))
        DO 20 I=1,NREC
           LL = (I-1)*12 + LT + NTLIST
           READ(NUNIT,'(12I6)') (IQ(L),L=LL,LL+11)
   20   CONTINUE
C
C get number of comment lines for different isotopes
        READ(NUNIT,'(I10)') NCOM
        NWW = NCOM * 80 + 2
C get CISO bank
        CALL CHKZEB(NWW,IXCONS)
        IF(LCSUP.EQ.0) Then
C create a top level bank for the isotope comments
          CALL MZBOOK(IXCONS,LCISO,LCSUP,1,'CISO',3,0,NWW,0,-1)
          LC = LCISO
        ELSE
C create an additional bank in the linear structure CISO
          CALL MZBOOK(IXCONS,LC,LCSUP,0,'CISO',3,0,NWW,0,-1)
          LCSUP = LC
        ENDIF
        IQ(LC+1) = NCOM
        DO 30 I=1,NCOM
           J = (I-1)*81 + 2
           READ(NUNIT,'(I4,I4,A70)') IQ(LC+J),IQ(LC+J+1),
     +                CCOMM
           CALL UCTOH(CCOMM,IQ(LC+J+2),4,70)
   30   CONTINUE
C
C---------------------------------------------------------------------
C check the existence of secondary x-section files stored in bank MIFL
C real messy code !!! But its fortran after all !!! CZ Jan 95
        XSFILE = ' '
        IF(NUNIT.EQ.MICROS) THEN
          FMIFL = .FALSE.
          CALL MZINQD(IXCONS)
          IF(LMIFIL.GE.IQUEST(3) .AND. LMIFIL.LE.IQUEST(4)) THEN
             CALL UHTOC(IQ(LMIFIL-4),4,CNAME,4)
             IF(CNAME.EQ.'MIFL') FMIFL = .TRUE.
          ENDIF
          IXSF=LMIFIL
        ENDIF
        IF(FMIFL) THEN
   40     CONTINUE
C get the file name
          CALL UHTOC(IQ(IXSF+2),4,XSFILE,IQ(IXSF+1))
C
          INQUIRE(FILE=XSFILE,EXIST=EXISTS)
          IF(.NOT.EXISTS) THEN
             PRINT '(70(''*''))'
             PRINT*,' * MICAP : x-section file not found : ',XSFILE
             PRINT '(70(''*''))'
          ELSE
            IXSF = IXSF + 101
C last name in the list ?
            IF(IXSF-LMIFIL .GE. IQ(LMIFIL-1) ) FMIFL = .FALSE.
C find a free unit number (greater 31), and use it
            DO 50 I=NUNIT+1,99
              INQUIRE(UNIT=I,OPENED=OPENED)
              IF(.NOT.OPENED) THEN
                 NUNIT = I
#if defined(CERNLIB_UNIX)||defined(CERNLIB_CRAY)||defined(CERNLIB_IBMVM)
                 OPEN(UNIT=I,FILE=XSFILE,STATUS='OLD')
#endif
#if defined(CERNLIB_VAX)
                 OPEN(UNIT=I,FILE=XSFILE,STATUS='OLD',READONLY)
#endif
                 GOTO 10
              ENDIF
   50       CONTINUE
            PRINT '(70(''*''))'
            PRINT *,'* MICAP : No more free units available !'
            PRINT '(70(''*''))'
          ENDIF
        ENDIF
C---------------------------------------------------------------------
      CALL VZERO(MATIDS,4000)
      LT = LTEMP
   60 CONTINUE
        NUNIT = IQ(LT + NTUNIT)
        KK = LT + NTLIST
        DO 90 I=1,100
           NIS = IQ(KK)
           KK = KK + 1
           IF(NIS.EQ.0) GOTO 100
           IF(MATIDS(I,1,1).EQ.0) THEN
              MATIDS(I,1,1) = NIS
              MATIDS(I,1,2) = NUNIT
C is the Z of the element correct?
           ELSE IF(IQ(KK)/1000.EQ.I) THEN
C overwrite existing element with the one stored in new file
              DO 70 J=2,MATIDS(I,1,1)+1
                MATIDS(I,J,1) = 0
                MATIDS(I,J,2) = 0
   70         CONTINUE
              MATIDS(I,1,1) = NIS
              MATIDS(I,1,2) = NUNIT
           ELSE
C no action
              KK = KK + 2 * NIS
              GOTO 90
           ENDIF
C maximal 20 isotopes per element
           NIS = MIN(NIS,20)
           DO 80 J=2,NIS+1
              MATIDS(I,J,1) = IQ(KK)
              MATIDS(I,J,2) = IQ(KK+1)
              KK = KK + 2
   80      CONTINUE
   90   CONTINUE
  100   CONTINUE
        LT = LQ(LT)
      IF(LT.GT.0) GOTO 60
C
C       DEFINE CROSS SECTION DIMENSIONING VARIABLES
C         NNR EQUALS THE NUMBER OF NEUTRON RECORDS
C         NQ EQUALS THE NUMBER OF Q VALUES
C         NGR EQUALS THE NUMBER OF GAMMA RECORDS
C       SET THE DEFAULT VALUES FOR THE CURRENT CROSS SECTION DATA
      NNR=134
      NQ=66
      NGR=60
C
C       SET THE DEFAULT VALUES FOR THE NEUTRON, PROTON, DEUTERON,
C       TRITON, HELIUM-3, AND ALPHA PARTICLE MASSES (IN EV)
      ZN=XMASS(1)*1.E9
      ZP=XMASS(0)*1.E9
      ZD=XMASS(7)*1.E9
      ZT=XMASS(8)*1.E9
      ZHE3=XMASS(9)*1.E9
      ZA=XMASS(10)*1.E9
C       SET THE DEFAULT VALUES FOR THE NEUTRON, PROTON, DEUTERON,
C       TRITON, HELIUM-3, AND ALPHA PARTICLE MASSES (IN AMU)
      XAMU=0.93149432*1.E9
      AN=ZN/XAMU
      AP=ZP/XAMU
      AD=ZD/XAMU
      AT=ZT/XAMU
      AHE3=ZHE3/XAMU
      AA=ZA/XAMU
C now preprocess all materials xs
      MEDIA = 0
      NMIX = 0
      NMAT = 0
#if defined(CERNLIB_MDEBUG)
      PRINT *,' MICAP-INI : setup materials '
      PRINT '('' NMATE='',I20,'' JMATE='',I20)',NMATE,JMATE
      PRINT '('' NTMED='',I20,'' JTMED='',I20)',NTMED,JTMED
#endif
C Check if material option bank MIST exists
      FMIST = .FALSE.
      CALL MZINQD(IXCONS)
      IF(LMIST.GE.IQUEST(3) .AND. LMIST.LE.IQUEST(4)) THEN
         CALL UHTOC(IQ(LMIST-4),4,CNAME,4)
         IF(CNAME.EQ.'MIST') FMIST = .TRUE.
      ENDIF
C 1. loop over tracking media -> get NMIX,MEDIA
      DO 140 I=1,NTMED
         JTM = LQ(JTMED - I)
         IF(JTM.LE.0) GOTO 140
C valid tracking medium found get material number
C and get corresponding material parameters from JMATE structure
         IMA = INT(Q(JTM+6))
         IF(IMA.LE.0 .OR. IMA.GT.NMATE) GOTO 140
C count number of elements and number of mixing operations
         JMA = LQ(JMATE-IMA)
         IF(JMA.LE.0) GOTO 140
         IF(Q(JMA+6) .LE. 1.0 .OR. Q(JMA+6) .GE. 240.) GOTO 140
C Check if for material IMA single isotopes are selected
         FSINGL = .FALSE.
         IF(FMIST) THEN
           DO 110 KIM=1,IQ(LMIST-1),2
              IF(IMA.EQ.IQ(LMIST+KIM).AND.IQ(LMIST+KIM+1).EQ.0) THEN
                 FSINGL = .TRUE.
                 GOTO 120
              ENDIF
  110      CONTINUE
  120      CONTINUE
         ENDIF
C get number of elements in material max = 100
         KK = MIN(ABS(Q(JMA+11)),100.)
C relation between MICAP and GEANT material number
         MEDIA = MEDIA + 1
C mixture ?
         KK1 = KK
         IF(KK.GT.1) THEN
            JMIXT = LQ(JMA - 5)
C
C check if more than one isotope has to taken into account for all
C elements in the mixture
            DO 130 K=1,KK
               IA = NINT(Q(JMIXT+K))
               IZ = NINT(Q(JMIXT+K+KK))
               CALL MATISO(IZ,IA,NNI,IDI,FSINGL,NUNIT)
               KK1 = KK1 + NNI - 1
  130       CONTINUE
         ELSE
            IA  = NINT(Q(JMA+6))
            IZ  = NINT(Q(JMA+7))
            CALL MATISO(IZ,IA,NNI,IDI,FSINGL,NUNIT)
            KK1 = KK1 + NNI - 1
         ENDIF
         NMIX = NMIX + KK1
  140 CONTINUE
C allocate ZEBRA bank for material information
      NW = 9 * NMIX + MEDIA + 10
C define link area for MICAP banks in GCBANK
      CALL CHKZEB(NW,IXCONS)
      CALL MZBOOK(IXCONS,LMOMA,0,2,'MOME',0,0,NW,0,-1)
      LMAG1 = LMOMA + 1
      IQ(LMAG1) = NMAGIC
      LGE2MO = LMAG1  + 1
      LFP10  = LGE2MO + MEDIA + 1
      LFP11  = LFP10  + NMIX
      LFP12  = LFP11  + NMIX
      LFP13  = LFP12  + NMIX
      LFP14  = LFP13  + NMIX
      LFP140 = LFP14  + NMIX
      LFP16  = LFP140 + NMIX
      LFP17  = LFP16  + NMIX
C 2. loop over tracking media
      MEDIA1 = 0
      NMIX1 = 0
      DO 230 I=1,NTMED
         JTM = LQ(JTMED - I)
         IF(JTM.LE.0) GOTO 230
C valid tracking medium found get material number
C and get corresponding material parameters from JMATE structure
         IMA = INT(Q(JTM+6))
#if defined(CERNLIB_MDEBUG)
         PRINT '('' IMATE ='',I10)',IMA
#endif
         IF(IMA.LE.0 .OR. IMA.GT.NMATE) GOTO 230
C count number of elements and number of mixing operations
         JMA = LQ(JMATE-IMA)
         IF(JMA.LE.0) GOTO 230
         IF(Q(JMA+6) .LE. 1.0 .OR. Q(JMA+6) .GE. 240.) GOTO 230
C Check if for material IMA single isotopes are selected
         FSINGL = .FALSE.
         IF(FMIST) THEN
           DO 150 KIM=1,IQ(LMIST-1),2
              IF(IMA.EQ.IQ(LMIST+KIM).AND.IQ(LMIST+KIM+1).EQ.0) THEN
                 FSINGL = .TRUE.
                 GOTO 160
              ENDIF
  150      CONTINUE
  160      CONTINUE
         ENDIF
         NMAT = NMAT + 1
C get number of elements in material max = 100
         RHO1 = Q(JMA+8)
         KK = MIN1(ABS(Q(JMA+11)),100.)
C relation between MICAP and GEANT material number
C check if medium IMA already stored (multiple tracking media)
         CALL VZERO(AGEA,100)
         DO 180 KMI=1,MEDIA1
            IF(IQ(LGE2MO+KMI).EQ.IMA) THEN
               IF(KK.EQ.1) THEN
                  IA  = NINT(Q(JMA+6))
                  IZ  = NINT(Q(JMA+7))
                  CALL MATISO(IZ,IA,NNI,IDI,FSINGL,NUNIT)
                  NMIX = NMIX - NNI
               ELSE
                  JMIXT = LQ(JMA - 5)
                  DO 170 K=1,KK
                     IA = NINT(Q(JMIXT+K))
                     IZ = NINT(Q(JMIXT+K+KK))
                     CALL MATISO(IZ,IA,NNI,IDI,FSINGL,NUNIT)
                     NMIX  = NMIX - NNI
  170             CONTINUE
               ENDIF
               MEDIA = MEDIA - 1
               GOTO 230
            ENDIF
  180    CONTINUE
         MEDIA1 = MEDIA1 + 1
         IQ(LGE2MO+MEDIA1) = IMA
C mixture ?
         KK2 = KK
         IF(KK.GT.1) THEN
            JMIXT = LQ(JMA - 5)
            KPOS = 1
            DO 200 K=1,KK
               AMOL = Q(LQ(JMIXT-1) + 2)
               XMOLCM = RHO1/AMOL*XNAVO
               IA = NINT(Q(JMIXT+K))
               IZ = NINT(Q(JMIXT+K+KK))
               CALL MATISO(IZ,IA,NNI,IDI,FSINGL,NUNIT)
               KK2 = KK2 + NNI - 1
               DO 190 KJ=1,NNI
                  KKPOS = KPOS + KJ - 1
                  IF(KJ.EQ.1) THEN
                     AGEA(KKPOS) = Q(JMIXT+K)
                  ELSE
                     AGEA(KKPOS) = 0.
                  ENDIF
                  MID(KKPOS,1) = IDI(KJ,1)
                  MID(KKPOS,2) = NUNIT
                  IIZ = IDI(KJ,1)/1000
                  IIA = IDI(KJ,1) - IIZ * 1000
                  IF(IIA.NE.0 .AND. NNI.GT.1.) THEN
                    A(KKPOS) = FLOAT(IIA)
                  ELSE
                    A(KKPOS) = Q(JMIXT+K)
                  ENDIF
                  Z(KKPOS) =  Q(JMIXT+K+KK)
                  WISO = FLOAT(IDI(KJ,2))/100.
                  WI = Q(JMIXT+K+2*KK)*AMOL/A(KKPOS)*WISO
                  DEN(KKPOS) = XMOLCM * WI
#if defined(CERNLIB_MDEBUG)
                  PRINT '('' MIXT: El #'',I3,'' A,Z :'',2F10.2,
     +                    '' Rho='',F10.2,'' Den ='',F10.5)',
     +                       KJ,A(KJ),Z(KJ),RHO1,DEN(KJ)
#endif
  190          CONTINUE
               KPOS = KPOS + NNI
  200       CONTINUE
C element or compound
         ELSE
            IA  = NINT(Q(JMA+6))
            IZ  = NINT(Q(JMA+7))
            CALL MATISO(IZ,IA,NNI,IDI,FSINGL,NUNIT)
            KK2 = KK2 + NNI - 1
            DO 210 KJ=1,NNI
               IF(KJ.EQ.1) THEN
                  AGEA(KJ) = Q(JMA+6)
               ELSE
                  AGEA(KJ) = 0.
               ENDIF
               MID(KJ,1) = IDI(KJ,1)
               MID(KJ,2) = NUNIT
               IIZ = IDI(KJ,1)/1000
               IIA = IDI(KJ,1) - IIZ * 1000
               IF(IIA.NE.0 .AND. NNI.GT.1.) THEN
                  A(KJ) = FLOAT(IIA)
               ELSE
                  A(KJ) = Q(JMA+6)
               ENDIF
               Z(KJ) =  Q(JMA+7)
               WISO = FLOAT(IDI(KJ,2))/100.
               DEN(KJ) = RHO1/A(KJ) * WISO *XNAVO
#if defined(CERNLIB_MDEBUG)
               PRINT '('' ELEM: Iso #'',I3,'' A,Z :'',2F10.2,
     +                 '' Rho='',F10.2,'' Den ='',F10.5)',
     +                  KJ,A(KJ),Z(KJ),RHO1,DEN(KJ)
#endif
  210       CONTINUE
         ENDIF
C
C fill MICAP material arrays
C actual number of isotopes given by KK2
C
         DO 220 J = NMIX1 + 1, NMIX1 + KK2
            IQ(LFP10+J-1) = MEDIA1
            IQ(LFP11+J-1) = MID(J-NMIX1,1)
C check if bound hydrogen has been selected
            IF(NINT(Z(J-NMIX1)).EQ.1.AND.KK.GT.1) IQ(LFP11+J-1) = 1000
            Q(LFP12+J-1) = DEN(J-NMIX1)
            IQ(LFP13+J-1) = NINT(Z(J-NMIX1))
            Q(LFP14+J-1) = A(J-NMIX1)
            Q(LFP140+J-1) = AGEA(J-NMIX1)
  220    CONTINUE
         NMIX1 = NMIX1 + KK2
  230 CONTINUE
      IF(NMIX.LE.0) THEN
         PRINT *,' GCALOR: NO tracking media found ===> STOP '
         STOP
      ENDIF
C read cross-sections and perform mixing and thinning
      CALL MOXSEC
C close MICAP cross-section file(s)
      LT = LTEMP
  240 CONTINUE
        CLOSE(UNIT=IQ(LT+NTUNIT))
        LT = LQ(LT)
      IF(LT.GT.0) GOTO 240
C Drop temporary linear structures
      CALL MZDROP(IXCONS,LTEMP,'L')
      CALL MZDROP(IXCONS,LCISO,'L')
      RETURN
      END
