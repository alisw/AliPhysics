*CMZ :          19/11/98  21.32.03  by  Federico Carminati
*-- Author :
C***********************************************************************
C
      SUBROUTINE ZDC_ADDANG(THETA1,PHI1,THETA2,PHI2,THETA3,PHI3)
C
C  This subroutine takes spherical polar angles THETA1 and PHI1
C  and adds to them the direction PHI2 on a cone of half-opening angle
C  THETA2 (= angle from cone axis to cone surface) to produce
C  the new space angles THETA3 and PHI3. All angles are in
C  degrees, with 0<=THETA<=180 and 0<=PHI<360. There is no
C  other restriction on THETA1, PHI1, THETA2 or PHI2.
C
C             27-AUG-1992          S. Coutu
C
C***********************************************************************
#undef CERNLIB_GEANT321_GCONSP_INC
#include "geant321/gconsp.inc"
      DOUBLE PRECISION THETA1,PHI1,THETA2,PHI2,THETA3,PHI3
      DOUBLE PRECISION TEMP,CX,CY,CZ,CT1,ST1,CT2,ST2,
     +     CP1,SP1,CP2,SP2,RTHETA3
C
      CT1=COS(THETA1*DEGRAD)
      ST1=SIN(THETA1*DEGRAD)
      CP1=COS(PHI1*DEGRAD)
      SP1=SIN(PHI1*DEGRAD)
      CT2=COS(THETA2*DEGRAD)
      ST2=SIN(THETA2*DEGRAD)
      CP2=COS(PHI2*DEGRAD)
      SP2=SIN(PHI2*DEGRAD)
      CX=CT1*CP1*ST2*CP2+ST1*CP1*CT2-SP1*ST2*SP2
      CY=CT1*SP1*ST2*CP2+ST1*SP1*CT2+CP1*ST2*SP2
      CZ=CT1*CT2-ST1*ST2*CP2
C
      RTHETA3=ACOS(CZ)
      THETA3=RTHETA3*RADDEG
      IF(THETA3.EQ.0.OR.THETA3.EQ.180.)THEN
         PHI3=0.
      ELSE
         TEMP=CX/SIN(RTHETA3)
         IF(TEMP.GT.1.)TEMP=1.
         IF(TEMP.LT.-1.)TEMP=-1.
         PHI3=ACOS(TEMP)*RADDEG
         IF(CY.LT.0.)PHI3=360.-PHI3
      ENDIF
*
      END

*CMZ :          21/11/98  17.33.14  by  Federico Carminati
*-- Author :
C==================================================================
      SUBROUTINE ZDC_EFERMI(ID)
*
#include "zdc_common.inc"
#undef CERNLIB_GEANT321_GCONSP_INC
#include "geant321/gconsp.inc"
      DIMENSION RNDM(3)
*
      CALL GRNDM(RNDM,3)
*
      IF(ID.EQ.14)THEN
         DO I=1,200
            IF(RNDM(3).GE.PROBINTP(I-1).AND.RNDM(3).LT.PROBINTP(I))GO TO
     $           10
         END DO
      ELSE IF(ID.EQ.13) THEN
         DO I=1,200
            IF(RNDM(3).GE.PROBINTN(I-1).AND.RNDM(3).LT.PROBINTN(I))GO TO
     $           10
         END DO
      ENDIF
*
   10 PIPPO=PP(I)+0.001
      PHI=TWOPI*RNDM(1)
      COST=1.-2*RNDM(2)
      TET=ACOS(COST)
      DDP(1)=PIPPO*SIN(TET)*COS(PHI)
      DDP(2)=PIPPO*SIN(TET)*SIN(PHI)
      DDP(3)=PIPPO*COST
*
      END

*CMZ :          23/02/99  07.53.48  by  Federico Carminati
*CMZ :  2.03/01 06/08/98  10.13.59  by  Federico Carminati
*-- Author :   E.SCOMPARIN, 13/05/1996.
      SUBROUTINE ZDC_END
C
C *** TERMINATION OF THE ZDC SIMULATION AT THE END OF A RUN ***
C
C CALLED BY : SXEND
C ORIGIN    : E. SCOMPARIN
C
C
*      CALL HROUT(0,ICYCLE,' ')
*      CALL HREND('ALIZ')
*      CLOSE(1)

      PRINT 10000
10000 FORMAT(1H /1H ,36('*'),' ZDC_END  ',36('*')/
     $ 1H ,20X,'ZDC simulation package terminated'/
     $ 1H ,80('*'))
*
      END

*CMZ :          23/02/99  07.53.48  by  Federico Carminati
*CMZ :  2.03/01 07/08/98  09.48.28  by  Federico Carminati
*-- Author :   E.SCOMPARIN, 13/05/1996.
      SUBROUTINE ZDC_EVE
C
C *** TERMINATION OF THE ZDC SIMULATION AFTER EACH EVENT ***
C
C CALLED BY : SXEVE
C ORIGIN    : E. SCOMPARIN
C
C --  RETRIEVE HITS INFORMATION FOR THE CURRENT EVENT
#undef CERNLIB_GEANT321_GCFLAG_INC
#include "geant321/gcflag.inc"
#include "zdc_common.inc"
*
      INTEGER NUMVSP(NVZP),NUMVSN(NVZN),
     +        NUMBP(NVZP,1000),NUMBN(NVZN,1000),
     +        ITRP(1000),ITRN(1000)
      REAL HITP(NHZ,1000),HITN(NHZ,1000)
*
      REAL ENEH(1000),XZPH(1000),YZPH(1000),EZPH(1000,NZPTX,NZPTY),
     +                XZNH(1000),YZNH(1000),EZNH(1000,NZNTX,NZNTY)
      INTEGER NPAH(1000),SFLH(1000),NPPH(1000,NZPTX,NZPTY),
     +                              NPNH(1000,NZNTX,NZNTY)
*
      DATA NUMVSP/0,0/
      DATA NUMVSN/0,0/
*
      DO JZERO=1,1000
         DO KZERO=1,NHZ
            HITP(KZERO,JZERO)=0
            HITN(KZERO,JZERO)=0
         ENDDO
      ENDDO
*
      DO I=1,1000
         NPAH(I)=0
         XZPH(I)=99.
         YZPH(I)=99.
         XZNH(I)=99.
         YZNH(I)=99.
      ENDDO
*
      ITOLD=0
*      CALL GTREVE
*
C --  RETRIEVE HITS INFORMATION FOR THE CURRENT EVENT
      IF(ISECF.EQ.1.OR.(ISECF.EQ.0.AND.MOD(IEVENT,1).EQ.0))THEN
         WRITE(6,*)'EVENT NO. ',IEVENT
         CALL GPHITS('*','*')
         WRITE(6,*)
      ENDIF
*
      NHMAXP=1000
      NHMAXN=1000
*
      CALL GFHITS(ZSET,ZP1,NVZP,NHZ,NHMAXP,0,
     +            NUMVSP,ITRP,NUMBP,HITP,NHITP)
      CALL GFHITS(ZSET,ZN1,NVZN,NHZ,NHMAXN,0,
     +            NUMVSN,ITRN,NUMBN,HITN,NHITN)
*
* Fill n-tuple
*
      ITCTP=0
      DO I=1,NHITP
         IF(ITRP(I).NE.ITCTP)THEN
            ITCTP=ITRP(I)
            NPAH(ITCTP)=HITP(1,I)
            ENEH(ITCTP)=HITP(2,I)
            XZPH(ITCTP)=HITP(3,I)
            YZPH(ITCTP)=HITP(4,I)
            SFLH(ITCTP)=INT(HITP(5,I))
         ENDIF
         J1=NUMBP(1,I)
         J2=NUMBP(2,I)
         EZPH(ITCTP,J1,J2)=HITP(7,I)
         NPPH(ITCTP,J1,J2)=HITP(6,I)
      ENDDO
      ITCTN=0
      DO I=1,NHITN
         IF(ITRN(I).NE.ITCTN)THEN
            ITCTN=ITRN(I)
            NPAH(ITCTN)=HITN(1,I)
            ENEH(ITCTN)=HITN(2,I)
            XZNH(ITCTN)=HITN(3,I)
            YZNH(ITCTN)=HITN(4,I)
            SFLH(ITCTN)=INT(HITN(5,I))
         ENDIF
         J1=NUMBN(1,I)
         J2=NUMBN(2,I)
         EZNH(ITCTN,J1,J2)=HITN(7,I)
         NPNH(ITCTN,J1,J2)=HITN(6,I)
      ENDDO
*
      DO I=1,1000
         IF(NPAH(I).NE.0)THEN
            NPA=NPAH(I)
            ENE=ENEH(I)
            XZP=XZPH(I)
            YZP=YZPH(I)
            SFL=SFLH(I)
            DO J=1,NZPTX
               DO K=1,NZPTY
                  EZP(J,K)=EZPH(I,J,K)
                  NPP(J,K)=NPPH(I,J,K)
               ENDDO
            ENDDO
            XZN=XZNH(I)
            YZN=YZNH(I)
            DO J=1,NZNTX
               DO K=1,NZNTY
                  EZN(J,K)=EZNH(I,J,K)
                  NPN(J,K)=NPNH(I,J,K)
               ENDDO
            ENDDO
*            CALL HFNT(10)
         ENDIF
      ENDDO
*
      END

*CMZ :          21/11/98  17.33.14  by  Federico Carminati
*-- Author :
C==================================================================
      SUBROUTINE ZDC_GFERMI(A,Z)
*
#include "zdc_common.inc"
#undef CERNLIB_GEANT321_GCONSP_INC
#include "geant321/gconsp.inc"
*     CALCOLO DELLA DISTRIBUZIONE DEI MOMENTI SECONDO
*     LA DISTRIBUZIONE DOPPIA GAUSSIANA (ILINOV), UGUALE
*     PER NEUTRONI E PROTONI
*
      PROBINTP(0)=0.0
      PROBINTN(0)=0.0
      SIG1=0.113
      SIG2=0.250
      ALFA=0.18*((A/12.)**0.333)
      XK=(4.*PI)/((1.+ALFA)*TWOPI**1.5)
*
      DO I=1,200
         P=FLOAT(I)*0.005
         PP(I)=P
         E1=(P*P)/(2.*SIG1*SIG1)
         E2=(P*P)/(2.*SIG2*SIG2)
         F1=EXP(-(E1))
         F2=EXP(-(E2))
         PROBP=XK*P*P*(F1/(SIG1**3)+F2/(SIG2**3))*0.005
         PROBINTP(I)=PROBINTP(I-1)+PROBP
         PROBINTN(I)=PROBINTP(I)
      END DO
*
      END

*CMZ :          23/02/99  07.53.48  by  Federico Carminati
*CMZ :  2.03/01 25/08/98  23.39.54  by  Federico Carminati
*-- Author :
C-- Author :    E.SCOMPARIN, 07/05/1996
      SUBROUTINE ZDC_INIT
C
C *** INITIALISATION FOR THE ZDC  SIMULATION ***
C
C CALLED BY : SXINIT
C ORIGIN    : E. SCOMPARIN
C
C
#undef CERNLIB_GEANT321_GCFLAG_INC
#include "geant321/gcflag.inc"
#undef CERNLIB_GEANT321_GCKINE_INC
#include "geant321/gckine.inc"
C
#include "zdc_common.inc"
*
* Initialize COMMON block ZDC_CGEOM
*
      DATA HDZN/4.0,4.0,50.0/
      DATA HDZP/8.0,8.0,75.0/
C Coordinates of the center of the ZDC front face in the MRS
      DATA ZNPOS/-0.5,0.,11613./
      DATA ZPPOS/-19.0,0.,11563./
      DATA FIZN/0.,0.01825,50.0/
      DATA FIZP/0.,0.01825,75.0/
      DATA GRZN/0.025,0.025,50.0/
      DATA GRZP/0.040,0.040,75.0/
      DATA NCEN/11,11,0/
      DATA NCEP/10,10,0/
*
* Initialize COMMON block ZDC_HITS
*
* Data for set ZDC
      DATA ZSET,ZN1,ZP1 /'ZDC ', 'ZN1 ','ZP1 '/
      DATA CHNMSZN/'ZNTX','ZN1 '/
      DATA CHNMSZP/'ZPTX','ZP1 '/
      DATA NBITSZN /2,2/
      DATA NBITSZP /2,2/
      DATA IDTYPZN,IDTYPZP /1, 2/
      DATA NWHI,NWDI /1000, 1000/
*
* Data for Detectors ZN1,ZP1
*
      DATA CHNAMH /'IPAR','EPAR','XPAR','YPAR','SFLG','FELE','ENER'/
      DATA NBITSH /8,32,32,32,1,32,32/
      DATA ORIG /0.,0.,10.,10.,0.,0.,0./
      DATA FACT /1.,1.,10000.,10000.,1.,10000.,10000./
      DATA NHSUM /2/
*
* Read tables of Cerenkov light produced in the fibers
*
      OPEN(67,FILE='ZDC/light22620362207s',
     +    FORM='FORMATTED', STATUS='OLD')
      OPEN(68,FILE='ZDC/light22620362208s',
     +    FORM='FORMATTED', STATUS='OLD')
      OPEN(69,FILE='ZDC/light22620362209s',
     +    FORM='FORMATTED', STATUS='OLD')
      OPEN(70,FILE='ZDC/light22620362210s',
     +    FORM='FORMATTED', STATUS='OLD')
      DO 10 K=1,NALFA
         READ (67,10000)(TABLE(1,K,J),J=1,NBE)
         READ (68,10000)(TABLE(2,K,J),J=1,NBE)
         READ (69,10000)(TABLE(3,K,J),J=1,NBE)
         READ (70,10000)(TABLE(4,K,J),J=1,NBE)
   10 CONTINUE
10000 FORMAT(7(1X,F10.5))
      CLOSE (67)
      CLOSE (68)
      CLOSE (69)
      CLOSE (70)

      IFERMI=1
*
*  Set Fermi flag according to SWITCH 3
      IF(ISWIT(3).EQ.1) IFERMI=1
*  define the LHC energy/nucleon
      EPERN=2760.
*
*  Open data file to read in particles from event generators
      IF(IHIJ.GE.1)THEN
         OPEN(27,FILE=FILEH,STATUS='OLD')
      ELSE IF(IVNU.GE.1)THEN
         OPEN(27,FILE=FILEV,STATUS='OLD')
      ENDIF
C  GET FERMI MOMENTUM DISTRIBUTIONS
      CALL ZDC_GFERMI(207.,82.)
C
*      CALL HROPEN (1,'ALIZ',FILEA,'N',1024,ISTAT)
*
*      CALL HBNT(10,'ALIZDC','D')
*      CALL HBNAME(10,'ALI',NPA,
*     +'NPA[0,112]:U,ENE,SFL[0,1]:U,'//
*     +'XZP:R*4:16:[-99,99],YZP:R*4:16:[-99,99],'//
*     +'EZP(2,2),NPP(2,2),'//
*     +'XZN:R*4:16:[-99,99],YZN:R*4:16:[-99,99],'//
*     +'EZN(2,2),NPN(2,2)')
*
      PRINT 10100
10100 FORMAT(1H /1H ,36('*'),' ZDC_INIT ',36('*')/
     + 1H ,20X,'ZDC simulation package initialised'/
     + 1H ,80('*'))

      END

*CMZ :          21/11/98  17.33.14  by  Federico Carminati
*-- Author :
      SUBROUTINE ZDC_KINE(NT)
*     Generation of event kinematics                              *
* ----------------------------------------------------------------*
#undef CERNLIB_GEANT321_GCKINE_INC
#include "geant321/gckine.inc"
#undef CERNLIB_GEANT321_GCFLAG_INC
#include "geant321/gcflag.inc"
#include "zdc_common.inc"
*
      DIMENSION KATT(10000,2),PATT(10000,4)
*
      INTEGER ITA(49),ITB(49)
      REAL PLAB(3),UBUF(1),UB(1)
      DIMENSION VERTEX(3)
      REAL BALP(4),DDDP(4),PLABS(4)
*
      DATA UBUF/0./
      DATA ITA/22,-11,11,0,-13,13,111,211,-211,130,321,-321,
     +2112,2212,-2212,310,221,3122,3222,3212,3112,3322,3312,
     +3334,-2112,-3122,-3112,-3212,-3222,-3322,-3312,-3334,
     +-15,15,411,-411,421,-421,431,-431,4122,0,0,0,0,0,0,0,0/
      DATA ITB/10,-12,12,0,-14,14,110,120,-120,230,130,-130,
     +1220,1120,-1120,0,220,2130,1130,1230,2230,1330,2330,3331,
     +-1220,-2130,-1130,-1230,-2230,-1330,-2330,-3331,-16,16,
     +-240,240,-140,140,0,0,2140,0,0,0,0,0,0,0,0/
      DATA IVR/0/
*
      DIMENSION RNDM(2)
      CHARACTER*20 PARTNAME
*
* --- BEAMTY=1 Gaussian beam (width sigx sigy, mean values fx fy)
* --- BEAMTY=2 Uniform beam  (width sigx sigy, mean values fx fy)
*     Data cards BMTY, INBM
*
      DO KZ=1,3
         VERTEX(KZ)=0
      ENDDO
      IF(BEAMTY.EQ.1) THEN
         CALL SXGAUS(RNDM,2)
         VERTEX(1)= VERTEX(1)+FX+RNDM(1)*SIGX
         VERTEX(2)= VERTEX(2)+FY+RNDM(2)*SIGY
         VERTEX(3)= 0.05
      ELSEIF(BEAMTY.EQ.2) THEN
         CALL GRNDM(RNDM,2)
         VERTEX(1)= VERTEX(1)+FX+(2.*RNDM(1)-1.)*SIGX
         VERTEX(2)= VERTEX(2)+FY+(2.*RNDM(2)-1.)*SIGY
         VERTEX(3)= 0.05
      ENDIF
*
*  Fill the vertex bank
*
      CALL GSVERT(VERTEX,0,0,0.,0,NVTX)
      IF(NVTX.EQ.0)WRITE(6,*)'Error defining vertex'
*
* Read HIJING event file
*
      IF(IHIJ.GE.1)THEN
   10    READ(27,*,END=70) NMUL,EATT,JATT,NT,NP,N0,N01,N10,N11,B,NATT
         IF(IHIJ.EQ.2)THEN
            WRITE(6,*)'HIJING EVENT HEADER= ',NMUL,EATT,JATT,NT,NP,
     +      N0,N01,N10,N11,B,NATT
         ENDIF
         READ(27,*) (KATT(I,1),KATT(I,2),(PATT(I,J),J=1,4),I=1,NATT)
         IVR=IVR+1
         IF(IVR.LT.IHIJF) GOTO 10
*
         BIMP=B
         PARTI=FLOAT(NT)+FLOAT(NP)
         RMULTI=FLOAT(NATT)
*
         DO I=1,NATT
* --  Remove particles with negative z-direction
            IF(PATT(I,3).LT.0.0) GOTO 30
* --  Remove (if required) the spectators
            IF(IHIJSP.EQ.1)THEN
               IF(KATT(I,2).EQ.0..OR.KATT(I,2).EQ.10.) GOTO 30
            ENDIF
* --  Define participant flag : UBUF(1)
            IF(KATT(I,2).EQ.0..OR.KATT(I,2).EQ.10.) THEN
               UBUF(1)=0.
            ELSE
               UBUF(1)=1.
            ENDIF
* --  Translate the particle ID from HIJING to GEANT
            DO J=1,49
               IF(KATT(I,1).EQ.ITA(J)) THEN
                  IPART=J
                  GOTO 20
               ENDIF
            ENDDO
   20       CONTINUE
*
            DO JJ=1,3
               PLAB(JJ)=PATT(I,JJ)
            ENDDO
*
* Rotations to account for beam divergence and crossing angle
            IF(BMDIV.NE.0.) CALL ZDC_ROTP(PLAB,0)
            IF(BMCRA.NE.0.) CALL ZDC_ROTP(PLAB,1)
*
* Fill the tracks bank
            CALL GSKINE(PLAB,IPART,NVTX,UBUF,1,NT)
            IF(IHIJ.EQ.2)THEN
               WRITE(6,*)'IPART= ',IPART,'PLAB= ',PLAB
            ENDIF
            IF(NT.EQ.0)WRITE(6,*)'Error defining track'
*
   30       CONTINUE
         ENDDO
*
* Read VENUS event file
*
      ELSE IF(IVNU.GE.1)THEN
   40    READ(27,*,END=80)LEVT,NREVT,NMUL,B,KOLEVT,COLEVT, PMXEVT,
     +   EGYEVT,NP,NT,NATT
         IF(LEVT.NE.1) THEN
            WRITE(6,*) ' Error in reading VENUS data'
            STOP
         ENDIF
         IF(IVNU.EQ.2)THEN
            WRITE(6,*)'VENUS EVENT HEADER= ',LEVT,NREVT,NATT,B,KOLEVT,
     +      COLEVT,PMXEVT,EGYEVT,NP,NT
         ENDIF
*
         DO I=1,NATT
*--  read produced particles
            READ(27,*,END=80) LPTL,NREVT,NRPTL,KATT(I,1), (PATT(I,J),J=
     +      1,4),AMASSS
*
*         lptl ......... record label (lptl=3)
*         nrevt ........ event number
*         nrptl ........ particle number
*         katt(1,1)..... particle id
*         patt ......... 4-momentum (px,py,pz,en) in lab
*         amasss ........ mass of particle
*
            IF(LPTL.NE.3) THEN
               WRITE(6,*) ' ERROR IN READING VENUS DATA'
               STOP
            ENDIF
         ENDDO
*
         IVR=IVR+1
         IF(IVR.LT.IVNUF) GOTO 40
*
         BIMP=B
         PARTI=(FLOAT(NT)+FLOAT(NP))/2.
         RMULTI=FLOAT(NATT)
*
         DO I=1,NATT
* --  remove particles in negative direction
            IF(PATT(I,3).LT.0.0) GOTO 60
* --  if required remove the spectators
            IF(IVNSP.EQ.1)THEN
               IF((KATT(I,1).EQ.1120.OR.KATT(I,1).EQ.1220).AND. (PATT(I
     +         ,4).LT.(EPERN*1.005).AND. PATT(I,4).GT.(EPERN*0.995)))
     +         GOTO 60
            ENDIF
* --  flag the spectators (UBUF(1)=0)
            IF((KATT(I,1).EQ.1120.OR.KATT(I,1).EQ.1220).AND. (PATT(I,4)
     +      .GT.(EPERN*0.995).AND. PATT(I,4).LT.(EPERN*1.005))) THEN
               UBUF(1)=0.
            ELSE
               UBUF(1)=1.
            ENDIF
* --  Translate the particle ID from VENUS to GEANT
            DO J=1,49
               IF(KATT(I,1).EQ.ITB(J))THEN
                  IPART=J
                  GOTO 50
               ENDIF
            ENDDO
   50       CONTINUE
            DO JJ=1,3
               PLAB(JJ)=PATT(I,JJ)
            ENDDO
*
* Rotations to account for beam divergence and crossing angle
            IF(BMCRA.NE.0.)CALL ZDC_ROTP(PLAB,1)
            IF(BMDIV.NE.0.)CALL ZDC_ROTP(PLAB,0)
*
* Fill the tracks bank
            CALL GSKINE(PLAB,IPART,NVTX,UBUF,1,NT)
            IF(IVNU.EQ.2)THEN
               WRITE(6,*)'IPART= ',IPART,'PLAB= ',PLAB
            ENDIF
            IF(NT.EQ.0)WRITE(6,*)'Error defining track'
   60       CONTINUE
         ENDDO
*
      ELSE
*
* Generation of test particles
*
* Definition of the incident particle
*
         IPART = NINT(PKINE(1))
         UBUF(1)=1
         CALL GFPART(IPART,PARTNAME,ITRTYP,AMASS,CHARGE,TLIFE,UB,NB)
         AAP=1.0
*
* --- In case of a ion , find out energy and momentum
* --- NAAP=A=Z+N - atomic number of the initial particle
         IF(IPART.GE.45) THEN
            NAAP = INT(0.5 + AMASS/0.93149432)
            AAP = NAAP
            UBUF(1)=0.
         ENDIF
         IF((IPART.GE.45).OR.(IPART.EQ.13).OR.(IPART.EQ.14)) THEN
            UBUF(1)=0.
         ENDIF
*
         IF(IPART.EQ.14)THEN
            AMA=0.938
         ELSE IF(IPART.EQ.13)THEN
            AMA=0.939
         ELSE
            AMA=AMASS
         ENDIF
*
         IF(PKINE(6).EQ.0.)THEN
            PLAB(1)=PKINE(2)*PKINE(3)
            PLAB(2)=PKINE(2)*PKINE(4)
            PLAB(3)=PKINE(2)*PKINE(5)
         ELSE
            SCANG=2*ATAN(EXP(-(PKINE(6))))
            PLAB(1)=-PKINE(2)*SIN(SCANG)
            PLAB(2)=0.
            PLAB(3)=PKINE(2)*COS(SCANG)
         ENDIF
*
* Rotations to account for beam divergence and crossing angle
         IF(BMCRA.NE.0.)CALL ZDC_ROTP(PLAB,1)
         IF(BMDIV.NE.0.)CALL ZDC_ROTP(PLAB,0)
*
* --  in case of a nucleon, if required, apply the Fermi momentum
         IF(UBUF(1).EQ.0.AND.IFERMI.EQ.1) THEN
            IF(IPART.EQ.13.OR.IPART.EQ.14)THEN
               CALL ZDC_EFERMI(IPART)
            ELSE IF(IPART.GE.45)THEN
               CALL ZDC_EFERMI(13)
               GOLDHABER=SQRT(AAP*(207-AAP)/(207-1.))
               DO JJ=1,3
                  DDP(JJ)=DDP(JJ)*GOLDHABER
               ENDDO
            ENDIF
            DO II=1,3
               BALP(II)=-PLAB(II)
            ENDDO
            BALP(4)=SQRT(BALP(1)**2+BALP(2)**2+BALP(3)**2+AMA**2)
            DO II=1,3
               DDDP(II)=DDP(II)
            ENDDO
            DDDP(4)=SQRT(DDDP(1)**2+DDDP(2)**2+DDDP(3)**2+AMA**2)
            CALL LORENF(AMA,BALP,DDDP,PLABS)
            DO II=1,3
               PLAB(II)=PLABS(II)
            ENDDO
c           WRITE(6,*)'BALP ',BALP,'DDDP ',DDDP,'PLABS ',PLABS
         ENDIF
*
         BIMP=999.
*
*  Fill the tracks bank
         CALL GSKINE(PLAB,IPART,NVTX,UBUF,1,NT)
         IF(NT.EQ.0)WRITE(6,*)'Error defining track'

      ENDIF
*
* Some debug printings
      IF(IHIJ.EQ.2.OR.IVNU.EQ.2)THEN
         CALL GPVERT(NVTX)
         DO KK=1,NT
            CALL GPKINE(KK)
         ENDDO
      ENDIF
*
C
      RETURN
   70 WRITE(6,*) 'END OF HIJING FILE'
      IEORUN=1
      GOTO 999
   80 WRITE(6,*) 'END OF VENUS FILE'
      IEORUN=1
*
  999 END

*CMZ :          19/11/98  09.54.52  by  Federico Carminati
*-- Author :
      SUBROUTINE  ZDC_ROTP(PLAB,ICROSS)
#undef CERNLIB_GEANT321_GCKINE_INC
#include "geant321/gckine.inc"
#undef CERNLIB_GEANT321_GCFLAG_INC
#include "geant321/gcflag.inc"
C
#undef CERNLIB_GEANT321_GCONSP_INC
#include "geant321/gconsp.inc"
#include "zdc_common.inc"
C
      DOUBLE PRECISION TETPART,FIPART,TETDIV,FIDIV,TETSUM,FISUM
      DOUBLE PRECISION DPLAB(3)
      DIMENSION PLAB(3)
      DIMENSION RVEC(1)
      DIMENSION RNDM(1)
C
C TAKE NON-ZERO DIVERGENCY
c        write(6,*)'PLAB (NODIV) = ',plab
      PMOD=0.0
      DO IJ=1,3
         DPLAB(IJ)=PLAB(IJ)
         PMOD=PMOD+PLAB(IJ)*PLAB(IJ)
      ENDDO
      PMOD=SQRT(PMOD)
      IF(ICROSS.EQ.0)THEN
C---  TETDIV is given in radiants and it is read by data cards
         CALL RNORML(RVEC,1)
c          WRITE(6,*)RVEC
         TETDIV=PKINE(10)*ABS(RVEC(1))
C---  PHIDIV is randomly chosen between 0 and 360
         CALL GRNDM(RNDM,1)
         FIDIV=RNDM(1)*TWOPI
      ELSE IF(ICROSS.EQ.1)THEN
         IF(ISWIT(7).EQ.0)THEN
            TETDIV=0.
            FIDIV=0.
         ELSE IF(ISWIT(7).EQ.1)THEN
            TETDIV=PKINE(4)
            FIDIV=0.
         ELSE IF(ISWIT(7).EQ.2)THEN
            TETDIV=PKINE(4)
            FIDIV=PIBY2
         ENDIF
      ENDIF
      TETPART=DATAN(DSQRT(DPLAB(1)**2+DPLAB(2)**2)/DPLAB(3))
      IF(DPLAB(2).NE.0..OR.DPLAB(1).NE.0.) THEN
         FIPART=DATAN2(DPLAB(2),DPLAB(1))
      ELSE
         FIPART=0.
      ENDIF
      IF (FIPART.LT.0)FIPART=FIPART+TWOPI
C        write(6,*)' TETPART= ',tetpart,' FIPART= ',fipart
      TETDIV=TETDIV*RADDEG
      FIDIV=FIDIV*RADDEG
      TETPART=TETPART*RADDEG
      FIPART=FIPART*RADDEG
c        CALL ZDC_ADDANG(TETDIV,FIDIV,TETPART,FIPART,TETSUM,FISUM)
      CALL ZDC_ADDANG(TETPART,FIPART,TETDIV,FIDIV,TETSUM,FISUM)
C        write(6,*)' TETSUM= ',tetsum,'FISUM= ',fisum
      TETSUM=TETSUM*DEGRAD
      FISUM=FISUM*DEGRAD
      PLAB(1)=PMOD*DSIN(TETSUM)*DCOS(FISUM)
      PLAB(2)=PMOD*DSIN(TETSUM)*DSIN(FISUM)
      PLAB(3)=PMOD*DCOS(TETSUM)
C        write(6,*)'PLAB (DIV) = ',plab
      RETURN
      END

*CMZ :          19/11/98  09.54.52  by  Federico Carminati
*CMZ :  2.03/01 25/08/98  23.40.04  by  Federico Carminati
*-- Author :   E.SCOMPARIN, 13/05/1996.
      SUBROUTINE ZDC_STEP
C
C *** RECORDING OF ZDC HITS AFTER EACH TRACKING STEP ***
C
C CALLED BY : SXSTEP
C ORIGIN    : E. SCOMPARIN
C
#undef CERNLIB_GEANT321_GCKINE_INC
#undef CERNLIB_GEANT321_GCONSP_INC
#undef CERNLIB_GEANT321_GCVOLU_INC
#undef CERNLIB_GEANT321_GCKMAX_INC
#undef CERNLIB_GEANT321_GCKING_INC
#undef CERNLIB_GEANT321_GCNUM_INC
#undef CERNLIB_GEANT321_GCTRAK_INC
#undef CERNLIB_GEANT321_GCFLAG_INC
#undef CERNLIB_GEANT321_GCTMED_INC
#undef CERNLIB_GEANT321_GCMATE_INC
#undef CERNLIB_GEANT321_GCSETS_INC
#include "geant321/gckine.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcking.inc"
#include "geant321/gcnum.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcflag.inc"
#include "geant321/gctmed.inc"
#include "geant321/gcmate.inc"
#include "geant321/gcsets.inc"
#include "zdc_common.inc"
*
      DIMENSION LNAM(15),LNUM(15)
      DIMENSION UB(1)
*
      DIMENSION VERT0(3),PVERT0(4)
      CHARACTER*20 CHNPART0
*
      DIMENSION XLOCM(3),XANGM(3),XLOCL(3),XANGL(3),FNULA(3),ZENTR(3)
      DATA ALFMAX /110.0/
*
* New track
*
      IF(ITRA.NE.ITOLD)THEN
         IF(INWVOL.EQ.1.AND.IDTYPE.GT.0)THEN
*
* Impact point on the ZDC (either ZPRO or ZNEU)
*
            NLEV0=NLEVEL
            CALL GFPATH(ISET,IDET,NUMBV,NLEV,LNAM,LNUM)
            CALL UCTOH('ZNEU',IZNEU,4,4)
            CALL UCTOH('ZPRO',IZPRO,4,4)
            DO I=NLEV,1,-1
               IF(LNAM(I).EQ.IZNEU.OR.LNAM(I).EQ.IZPRO)THEN
                  NLEVEL=0
                  CALL GLVOLU(I,LNAM,LNUM,IER)
                  IF(IDTYPE.EQ.1)THEN
                     CALL GMTOD(VECT,HITSBUPN(3),1)
                     HITSBUPP(3)=99.
                     HITSBUPP(4)=99.
                  ELSE IF(IDTYPE.EQ.2)THEN
                     CALL GMTOD(VECT,HITSBUPP(3),1)
                     HITSBUPN(3)=99.
                     HITSBUPN(4)=99.
                  ENDIF
                  NLEVEL=NLEV0
                  GOTO 10
               ENDIF
            ENDDO
   10       CONTINUE
*
* Particle type and energy for the primary track
*
            CALL GFKINE(ITRA,VERT0,PVERT0,IPART0,NVERT0,UB,NB)
            CALL GFPART(IPART0,CHNPART0,ITRTYP0,AMASS0,CHARGE0,TLIFE0,
     +      UB,NB)
            PVERTM=VMOD(PVERT0,3)
            GEKIN0=SQRT(PVERTM**2+AMASS0**2)
            IF(IPART.NE.IPART0.OR.GEKIN.LT.GEKIN0*0.95)THEN
               HITSBUPP(1)=FLOAT(IPART0)
               HITSBUPP(2)=GEKIN0
               HITSBUPP(5)=1.
               HITSBUPN(1)=FLOAT(IPART0)
               HITSBUPN(2)=GEKIN0
               HITSBUPN(5)=1.
            ELSE
               HITSBUPP(1)=FLOAT(IPART)
               HITSBUPP(2)=GEKIN
               HITSBUPP(5)=0.
               HITSBUPN(1)=FLOAT(IPART)
               HITSBUPN(2)=GEKIN
               HITSBUPN(5)=0.
            ENDIF
            HITSBUPP(6)=0.
            HITSBUPP(7)=0.
            HITSBUPN(6)=0.
            HITSBUPN(7)=0.
            IF(IDTYPE.EQ.1)THEN
               CALL GSAHIT(ISET,IDET,ITRA,NUMBV,HITSBUPN,IHIT)
            ELSE IF(IDTYPE.EQ.2)THEN
               CALL GSAHIT(ISET,IDET,ITRA,NUMBV,HITSBUPP,IHIT)
            ENDIF
            ITOLD=ITRA
         ENDIF
      ENDIF
*
* Energy deposition in the passive material
*
      IF(INWVOL.EQ.0.AND.IDTYPE.GT.0)THEN
         HITS(6)=0.
         IF(ISTOP.EQ.2)THEN
            HITS(7)=GEKIN
         ELSE
            HITS(7)=DESTEP
         ENDIF
         IF(IDTYPE.EQ.1)THEN
            DO KCOPY=1,5
               HITSBUPN(KCOPY)=HITS(KCOPY)
            ENDDO
         ELSE IF(IDTYPE.EQ.2)THEN
            DO KCOPY=1,5
               HITSBUPP(KCOPY)=HITS(KCOPY)
            ENDDO
         ENDIF
         CALL GSCHIT(ISET,IDET,ITRA,NUMBV,HITS,NHSUM,IHIT)
      ENDIF
*
* Select charged relativistic particles
*
      IF(ABS(CHARGE).EQ.0)RETURN
      IF(GEKIN.GT.0.00001) THEN
         BETA  = VECT(7)/GETOT
      ELSE
         RETURN
      ENDIF
      IF(BETA.LT.0.67) RETURN

*
* Particle inside fiber
*
      IF(NUMED.EQ.803.AND.IDTYPE.GT.0) THEN
*
         XLOCM(1)= VECT(1)
         XLOCM(2)= VECT(2)
         XLOCM(3)= VECT(3)
         XANGM(1)= VECT(4)
         XANGM(2)= VECT(5)
         XANGM(3)= VECT(6)
*
* --- Calculating the position of the centre of the fiber in MARS
         DO KZERO=1,3
            FNULA(KZERO)=0
            ZENTR(KZERO)=0
         ENDDO
         CALL GDTOM(FNULA,ZENTR,1)
*
* --- Angle between the particle trajectory and the fiber axis
*
         CALL GMTOD(XANGM,XANGL,2)
         IF(XANGL(3).GT.1.0) XANGL(3)=SIGN(1.0,XANGL(3))
         ALFAR=ACOS(XANGL(3))
         ALFAD=ALFAR*RADDEG
*
         IF (ALFAD.GT.ALFMAX) RETURN
         IALFA = INT(1.+ALFAD/2.)
*
* --- Distance between the particle trajectory and the fiber axis
*
         CALL GMTOD(XLOCM,XLOCL,1)
         IF(ABS(XANGL(1)).GT.0.0001) THEN
            DCOEFF = XANGL(2)/XANGL(1)
            BE=ABS(XLOCL(2)-DCOEFF*XLOCL(1))/SQRT(DCOEFF*DCOEFF+1.)
         ELSE
            BE = ABS(XLOCL(1))
         ENDIF
         IF(IDTYPE.EQ.1)THEN
            RADIUS=FIZN(2)
         ELSEIF(IDTYPE.EQ.2)THEN
            RADIUS=FIZP(2)
         ENDIF
         IF(BE.GT.RADIUS) THEN
            PRINT*,'DISTANCE GREATER THAN FIBER RADIUS'
            RETURN
         ENDIF
         IBE = 1+INT(1000.*BE)
         IF(IBE.GT.NBE)IBE=NBE
*
* --- Finding out proper index for beta
         IF(BETA.GE.0.67 .AND. BETA.LE.0.75) IBETA=1
         IF(BETA.GT.0.75 .AND. BETA.LE.0.85) IBETA=2
         IF(BETA.GT.0.85 .AND. BETA.LE.0.95) IBETA=3
         IF(BETA.GT.0.95 .AND. BETA.LE.1.00) IBETA=4
*
* --- Looking into the tables
         OUT = CHARGE*CHARGE*TABLE(IBETA,IALFA,IBE)
         CALL RNPSSN(OUT,N,IERROR)
         HITS(6)=FLOAT(N)
         HITS(7)=0.
         IF(IDTYPE.EQ.1)THEN
            CALL UCOPY(HITSBUPN,HITS,4)
         ELSE IF(IDTYPE.EQ.2)THEN
            CALL UCOPY(HITSBUPP,HITS,4)
         ENDIF
         CALL GSCHIT(ISET,IDET,ITRA,NUMBV,HITS,NHSUM,IHIT)
      ENDIF
*
  999 END

*CMZ :          18/11/98  16.55.53  by  Federico Carminati
*CMZ :  2.03/01 06/08/98  10.13.59  by  Federico Carminati
*-- Author :   E. SCOMPARIN
      SUBROUTINE ZDC_TRKI(IFLAG)
C
C *** DEFINITION OF THE GEOMETRY OF THE PHOS ***
C *** NVE 24-SEP-1990 CERN GENEVA ***
C
C CALLED BY : SXTRKI
C ORIGIN    : E. SCOMPARIN
C
#undef CERNLIB_GEANT321_GCKINE_INC
#undef CERNLIB_GEANT321_GCKMAX_INC
#undef CERNLIB_GEANT321_GCKING_INC
#undef CERNLIB_GEANT321_GCTRAK_INC
#undef CERNLIB_GEANT321_GCTMED_INC
#undef CERNLIB_GEANT321_GCFLAG_INC
#undef CERNLIB_GEANT321_GCMATE_INC
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gctrak.inc"
#include "geant321/gctmed.inc"
#include "geant321/gcflag.inc"
#include "geant321/gcmate.inc"
#include "zdc_common.inc"

      REAL UBUF(1)
      DATA ITOLDG/0/

* Fetch new track (VENUS or HIJING)
      IF(ITRA.NE.ITOLDG)THEN
         CALL GFKINE(ITRA,VERT,PVERT,IPART,NVERT,UBUF,NUBUF)
         IFLPART=INT(UBUF(1))
         IF(IHIJ.EQ.2.OR.IVNU.EQ.2)WRITE(6,*)'ITRA,IFLPART ',ITRA,
     +   IFLPART
         ITOLDG=ITRA
      ENDIF
*
      END

*CMZ :          15/02/99  14.40.30  by  Federico Carminati
*-- Author :  Rene Brun
      SUBROUTINE ZDC_DEFS
C
C *** DEFINITION OF THE ZDC parameters
C     CALLED BY : SXDEFS
C
      end

      subroutine zdc_setbeam(beam, fxoff, fyoff, sx, sy,
     +               div, angle, cross)
* beam  : 1 = gaussian beam
*       : 2 = uniform beam
* fxoff : x-coordinate of beam offset
* fyoff : y-coordinate of beam offset
* sx    : sigma-x of the beam (gaussian or uniform)
* sy    : sigma-y of the beam (gaussian or uniform)
* div   : divergency of the beam (32*10**-6 rad for LHC)
* angle : beam crossing angle (100*10**-6 rad for LHC)
* cross : 1 = horizontal beam crossing
*       : 2 = vertical beam crossing
#include "zdc_common.inc"
      integer beam, cross
      real fxoff,fyoff,sx,sy,div,angle
      beamty = type
      fx     = fxoff
      fy     = fyoff
      sigx   = sx
      sigy   = sy
      bmdiv  = div
      bmcra  = angle
      iflcr  = cross
      end

      subroutine zdc_sethijing(hij, hijf, hijsp, file)
* IHIJ  : 1 = read HIJING event file
*       : 2 =  "     "      "    "    + debug
* IHIJF : event number of the first event to be read from file
* IHIJSP: 0 = read all particles
*       : 1 = remove spectator nucleons
#include "zdc_common.inc"
      integer hij,hijf,hijsp
      character *(*) file
      ihij   = hij
      ihijf  = hijf
      ihijsp = hijsp
      fileh  = file
      end

      subroutine zdc_setvenus(hiv, hivf, hivsp, file)
* IHIV  : 1 = read VENUS event file
*       : 2 =  "     "      "    "    + debug
* IHIVF : event number of the first event to be read from file
* IHIVSP: 0 = read all particles
*       : 1 = remove spectator nucleons
#include "zdc_common.inc"
      integer hiv,hivf,hivsp
      character *(*) file
      ihiv   = hiv
      ihivf  = hivf
      ihivsp = hivsp
      filev  = file
      end

      subroutine zdc_setkine(code, pmom, cx, cy, cz, type, fermi)
* code     : GEANT code of the test particle
* pmom     : absolute value of particle momentum
* cx,cy,cz : director cosines of the track (if type)
* type     :  0 = take director cosines from cx,cy,cz)
*          : <>0 = pseudorapidity of the test particle
* fermi    : 0 = no Fermi motion for the spectator nucleons
*          : 1 = Fermi motion for the spectator nucleons
#undef CERNLIB_GEANT321_GCKINE_INC
#include "geant321/gckine.inc"
#include "zdc_common.inc"
      integer code,type,fermi
      real pmom,cx,cy,cz
      pkine(1)  = code
      pkine(2)  = pmom
      pkine(3)  = cx
      pkine(4)  = cy
      pkine(5)  = cz
      pkine(6)  = type
      ifermi    = fermi
      end
*CMZ :          19/11/98  09.31.36  by  Federico Carminati
*-- Author :    Federico Carminati   19/11/98
      SUBROUTINE SXGAUS(RVEC,N)
#undef CERNLIB_GEANT321_GCONSP_INC
#include "geant321/gconsp.inc"
*KEND.
*
* Vector of random numbers distributed like a gaussian
      DIMENSION RNDM(2), RVEC(N)
      DO K=1,(N-1)/2+1
         CALL GRNDM(RNDM,2)
         PHI = TWOPI*RNDM(1)
         RHO = SQRT(-2.*LOG(RNDM(2)))
         RVEC(2*K-1) = RHO*SIN(PHI)
         IF(2*K.LE.N) RVEC(2*K  ) = RHO*COS(PHI)
      ENDDO
*
      END
