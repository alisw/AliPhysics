*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:58  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.38  by  S.Giani
*-- Author :
      SUBROUTINE NUCREC(NOPT,IREC)
C
C *** NUCLEAR REACTION KINEMATICS AT LOW ENERGIES ***
C *** NVE 18-MAY-1988 CERN GENEVA ***
C
C CALLED BY : GHEISH, GNSLWD
C ORIGIN    : H.FESEFELDT (12-FEB-1987)
C
C NOPT=1   N M(A,Z) --> G (G) M(A+1,Z  )    NEUTRON CAPTURE
C NOPT=2   N M(A,Z) --> N (G) M(A  ,Z  )    INELASTIC NEUTRON SCATT.
C NOPT=3   N M(A,Z) --> P (G) M(A  ,Z-1)
C NOPT=4   N M(A,Z) --> D (G) M(A-1,Z-1)
C NOPT=5   N M(A,Z) --> T (G) M(A-2,Z-1)
C NOPT=6   N M(A,Z) --> ALP.  M(A-3,Z-2)
C NOPT=7   N M(A,Z) --> N N   M(A-1,Z  )
C NOPT=8   N M(A,Z) --> N P   M(A-1,Z-1)
C NOPT=9   N M(A,Z) --> P P   M(A-1,Z-2)
C NOPT=11  P M(A,Z) --> G (G) M(A+1,Z+1)    PROTON CAPTURE
C NOPT=12  P M(A,Z) --> N (G) M(A  ,Z  )    INELASTIC PROTON SCATT.
C NOPT=13  P M(A,Z) --> P (G) M(A  ,Z+1)
C NOPT=14  P M(A,Z) --> D (G) M(A-1,Z  )
C NOPT=15  P M(A,Z) --> T (G) M(A-2,Z  )
C NOPT=16  P M(A,Z) --> ALP.  M(A-3,Z-1)
C NOPT=17  P M(A,Z) --> N N   M(A-1,Z+1)
C NOPT=18  P M(A,Z) --> N P   M(A-1,Z  )
C NOPT=19  P M(A,Z) --> P P   M(A-1,Z-1)
C SIMILAR FOR D,T,ALPHA SCATTERING ON NUCLEI
C
C NOTE : DOUBLE PRECISION CALCULATIONS ARE VITAL FOR THESE LOW
C        ENERGY PROCESSES
C        THEREFORE THE VARS OF /GENIO/ ARE DECLARED DOUBLE PRECISION
C        ALSO A DOUBLE PRECISION VERSION OF THE PHASE SPACE PACKAGE
C        "PHPNUC" HAS BEEN INTRODUCED
C *** HMF 29-AUG-1989 RWTH AACHEN ***
C
#include "geant321/s_defcom.inc"
#include "geant321/s_nucio.inc"
C
      DIMENSION QVAL(10),TCH(10)
      DIMENSION RNDM(2)
C
C** PROGRAM RETURNS WITH NOPT=0, IF INELASTIC SCATTERING ENERGETICALLY
C** NOT POSSIBLE, OR IF WRONG PARTICLES ENTER THIS ROUTINE: ONLY FOR
C** PROTONS,NEUTRONS, DEUTERIUM, TRITIUM AND ALPHAS.
C** IF EK > 100 MEV, THIS ROUTINE IS CERTAINLY NOT ADEQUATE.
C
      NOPT=0
      IF (IREC .EQ. 0) GO TO 9999
C
      IF (NPRT(9) .AND. (EK .GT. 0.1)) PRINT 9000,EK,IPART
 9000 FORMAT(' *NUCREC* ENERGY TOO HIGH EK = ',G12.5,' GEV ',
     $ ' KPART = ',I3)
      IF (EK .GT. 0.1) GO TO 9999
C
C%%%      IF(IPART.EQ.16) GOTO 2
C%%%      IF(IPART.EQ.14) GOTO 3
C%%%      IF(IPART.EQ.30) GOTO 4
C%%%      IF(IPART.EQ.31) GOTO 5
C%%%      IF(IPART.EQ.32) GOTO 6
C%%%      GO TO 9999
C%%%    2 AMAS = ATOMAS(1.,0.)
C%%%      GOTO 8
C%%%    3 AMAS = ATOMAS(1.,1.)
C%%%      GOTO 8
C%%%    4 AMAS = ATOMAS(2.,1.)
C%%%      GOTO 8
C%%%    5 AMAS = ATOMAS(3.,1.)
C%%%      GOTO 8
C%%%    6 AMAS = ATOMAS(4.,2.)
C
      IF (IPART .EQ. 16) GO TO 8
      IF (IPART .EQ. 14) GO TO 8
      IF (IPART .EQ. 30) GO TO 8
      IF (IPART .EQ. 31) GO TO 8
      IF( IPART .EQ. 32) GO TO 8
      GO TO 9999
C** SET BEAM PARTICLE, TAKE EK AS FUNDAMENTAL QUANTITY
C** DUE TO THE DIFFICULT KINEMATIC, ALL MASSES HAVE TO BE ASSIGNED
C** THE BEST MEASURED VALUES.
 8    CONTINUE
      CALL VZERO(QVAL,10)
      CALL VZERO(TCH ,10)
C --- GET MASS WHICH MATCHES GEANT ---
      AMAS=RMASS(IPART)
      EN=EK+ABS(AMAS)
      P =SQRT(ABS(EN*EN-AMAS*AMAS))
      PP=SQRT(PX*PX+PY*PY+PZ*PZ)
      IF (PP .GT. 1.0E-6) GO TO 8000
      CALL GRNDM(RNDM,2)
      PHINVE=TWPI*RNDM(1)
      COST=-1.+2.*RNDM(2)
      IF (COST .LE. -1.) COST=-1.
      IF (COST .GE.  1.) COST= 1.
      RTHNVE=ACOS(COST)
      PX=SIN(RTHNVE)*COS(PHINVE)
      PY=SIN(RTHNVE)*SIN(PHINVE)
      PZ=COS(RTHNVE)
      PP=1.
 8000 CONTINUE
      PX=PX/PP
      PY=PY/PP
      PZ=PZ/PP
      CALL VZERO(PV,10*MXGKPV)
      PV(1,1) =PX*P
      PV(2,1) =PY*P
      PV(3,1) =PZ*P
      PV(4,1) =EN
      PV(5,1) =AMAS
      PV(6,1) =NCH
      PV(7,1) =TOF
      PV(8,1) =IPART
      PV(9,1) =0.
      PV(10,1)=USERW
      PV(1,2) =0.
      PV(2,2) =0.
      PV(3,2) =0.
      PV(4,2) =0.
      PV(5,2) =ATOMAS(ATNO2,ZNO2)
      PV(6,2) =ZNO2
      PV(7,2) =TOF
      PV(8,2) =0.
      PV(9,2) =0.
      PV(10,2)=0.
C** CALCULATE Q-VALUE OF REACTIONS
      IF(IPART.EQ.16) GOTO 20
      IF(IPART.EQ.14) GOTO 30
      IF(IPART.EQ.30) GOTO 40
      IF(IPART.EQ.31) GOTO 50
      IF(IPART.EQ.32) GOTO 60
   20 PV(5,11)=ATOMAS(ATNO2+1.,ZNO2   )
      PV(6,11)=ZNO2
      PV(5,21)=0.
      PV(6,21)=0.
      PV(8,21)=1.
      PV(5,31)=0.
      PV(6,31)=0.
      PV(8,31)=1.
C
      PV(5,12)=PV(5,2)
      PV(6,12)=PV(6,2)
      PV(5,22)=RMASS(16)
      PV(6,22)=0.
      PV(8,22)=16.
      PV(5,32)=0.
      PV(6,32)=0.
      PV(8,32)=1.
C
      PV(5,13)=ATOMAS(ATNO2   ,ZNO2-1.)
      PV(6,13)=ZNO2-1.
      PV(5,23)=RMASS(14)
      PV(6,23)=1.
      PV(8,23)=14.
      PV(5,33)=0.
      PV(6,33)=0.
      PV(8,33)=1.
C
      PV(5,14)=ATOMAS(ATNO2-1.,ZNO2-1.)
      PV(6,14)=ZNO2-1.
      PV(5,24)=RMASS(30)
      PV(6,24)=1.
      PV(8,24)=30.
      PV(5,34)=0.
      PV(6,34)=0.
      PV(8,34)=1.
C
      PV(5,15)=ATOMAS(ATNO2-2.,ZNO2-1.)
      PV(6,15)=ZNO2-1.
      PV(5,25)=RMASS(31)
      PV(6,25)=1.
      PV(8,25)=31.
      PV(5,35)=0.
      PV(6,35)=0.
      PV(8,35)=1.
C
      PV(5,16)=ATOMAS(ATNO2-3.,ZNO2-2.)
      PV(6,16)=ZNO2-2.
      PV(5,26)=RMASS(32)
      PV(6,26)=2.
      PV(8,26)=32.
      PV(5,36)=0.
      PV(6,36)=0.
      PV(8,36)=1.
C
      PV(5,17)=ATOMAS(ATNO2-1.,ZNO2   )
      PV(6,17)=ZNO2
      PV(5,27)=PV(5,22)
      PV(6,27)=0.
      PV(8,27)=16.
      PV(5,37)=PV(5,22)
      PV(6,37)=0.
      PV(8,37)=16.
C
      PV(5,18)=PV(5,14)
      PV(6,18)=PV(6,14)
      PV(5,28)=PV(5,22)
      PV(6,28)=0.
      PV(8,28)=16.
      PV(5,38)=PV(5,23)
      PV(6,38)=1.
      PV(8,38)=14.
C
      PV(5,19)=ATOMAS(ATNO2-1.,ZNO2-2.)
      PV(6,19)=ZNO2-2.
      PV(5,29)=PV(5,23)
      PV(6,29)=1.
      PV(8,29)=14.
      PV(5,39)=PV(5,23)
      PV(6,39)=1.
      PV(8,39)=14.
C
      GOTO 70
   30 PV(5,11)=ATOMAS(ATNO2+1.,ZNO2+1.)
      PV(6,11)=ZNO2+1.
      PV(5,21)=0.
      PV(6,21)=0.
      PV(8,21)=1.
      PV(5,31)=0.
      PV(6,31)=0.
      PV(8,31)=1.
C
      PV(5,12)=ATOMAS(ATNO2   ,ZNO2+1.)
      PV(6,12)=ZNO2+1.
      PV(5,22)=RMASS(16)
      PV(6,22)=0.
      PV(8,22)=16.
      PV(5,32)=0.
      PV(6,32)=0.
      PV(8,32)=1.
C
      PV(5,13)=PV(5,2)
      PV(6,13)=PV(6,2)
      PV(5,23)=RMASS(14)
      PV(6,23)=1.
      PV(8,23)=14.
      PV(5,33)=0.
      PV(6,33)=0.
      PV(8,33)=1.
C
      PV(5,14)=ATOMAS(ATNO2-1.,ZNO2   )
      PV(6,14)=ZNO2
      PV(5,24)=RMASS(30)
      PV(6,24)=1.
      PV(8,24)=30.
      PV(5,34)=0.
      PV(6,34)=0.
      PV(8,34)=1.
C
      PV(5,15)=ATOMAS(ATNO2-2.,ZNO2   )
      PV(6,15)=ZNO2
      PV(5,25)=RMASS(31)
      PV(6,25)=1.
      PV(8,25)=31.
      PV(5,35)=0.
      PV(6,35)=0.
      PV(8,35)=1.
C
      PV(5,16)=ATOMAS(ATNO2-3.,ZNO2-1.)
      PV(6,16)=ZNO2-1.
      PV(5,26)=RMASS(32)
      PV(6,26)=2.
      PV(8,26)=32.
      PV(5,36)=0.
      PV(6,36)=0.
      PV(8,36)=1.
C
      PV(5,17)=ATOMAS(ATNO2-1.,ZNO2+1.)
      PV(6,17)=ZNO2+1.
      PV(5,27)=PV(5,22)
      PV(6,27)=0.
      PV(8,27)=16.
      PV(5,37)=PV(5,22)
      PV(6,37)=0.
      PV(8,37)=16.
C
      PV(5,18)=PV(5,14)
      PV(6,18)=PV(6,14)
      PV(5,28)=PV(5,22)
      PV(6,28)=0.
      PV(8,28)=16.
      PV(5,38)=PV(5,23)
      PV(6,38)=1.
      PV(8,38)=14.
C
      PV(5,19)=ATOMAS(ATNO2-1.,ZNO2-1.)
      PV(6,19)=ZNO2-1.
      PV(5,29)=PV(5,23)
      PV(6,29)=1.
      PV(8,29)=14.
      PV(5,39)=PV(5,23)
      PV(6,39)=1.
      PV(8,39)=14.
C
      NOPT=10
      GOTO 70
   40 PV(5,11)=ATOMAS(ATNO2+2.,ZNO2+1.)
      PV(6,11)=ZNO2+1.
      PV(5,21)=0.
      PV(6,21)=0.
      PV(8,21)=1.
      PV(5,31)=0.
      PV(6,31)=0.
      PV(8,31)=1.
C
      PV(5,12)=ATOMAS(ATNO2+1.,ZNO2+1.)
      PV(6,12)=ZNO2+1.
      PV(5,22)=RMASS(16)
      PV(6,22)=0.
      PV(8,22)=16.
      PV(5,32)=0.
      PV(6,32)=0.
      PV(8,32)=1.
C
      PV(5,13)=ATOMAS(ATNO2+1.,ZNO2   )
      PV(6,13)=ZNO2
      PV(5,23)=RMASS(14)
      PV(6,23)=1.
      PV(8,23)=14.
      PV(5,33)=0.
      PV(6,33)=0.
      PV(8,33)=1.
C
      PV(5,14)=PV(5,2)
      PV(6,14)=PV(6,2)
      PV(5,24)=RMASS(30)
      PV(6,24)=1.
      PV(8,24)=30.
      PV(5,34)=0.
      PV(6,34)=0.
      PV(8,34)=1.
C
      PV(5,15)=ATOMAS(ATNO2-1.,ZNO2   )
      PV(6,15)=ZNO2
      PV(5,25)=RMASS(31)
      PV(6,25)=1.
      PV(8,25)=31.
      PV(5,35)=0.
      PV(6,35)=0.
      PV(8,35)=1.
C
      PV(5,16)=ATOMAS(ATNO2-2.,ZNO2-1.)
      PV(6,16)=ZNO2-1.
      PV(5,26)=RMASS(32)
      PV(6,26)=2.
      PV(8,26)=32.
      PV(5,36)=0.
      PV(6,36)=0.
      PV(8,36)=1.
C
      PV(5,17)=ATOMAS(ATNO2   ,ZNO2+1.)
      PV(6,17)=ZNO2+1.
      PV(5,27)=PV(5,22)
      PV(6,27)=0.
      PV(8,27)=16.
      PV(5,37)=PV(5,22)
      PV(6,37)=0.
      PV(8,37)=16.
C
      PV(5,18)=PV(5,14)
      PV(6,18)=PV(6,14)
      PV(5,28)=PV(5,22)
      PV(6,28)=0.
      PV(8,28)=16.
      PV(5,38)=PV(5,23)
      PV(6,38)=1.
      PV(8,38)=14.
C
      PV(5,19)=ATOMAS(ATNO2   ,ZNO2-1.)
      PV(6,19)=ZNO2-1.
      PV(5,29)=PV(5,23)
      PV(6,29)=1.
      PV(8,29)=14.
      PV(5,39)=PV(5,23)
      PV(6,39)=1.
      PV(8,39)=14.
C
      NOPT=20
      GOTO 70
   50 PV(5,11)=ATOMAS(ATNO2+3.,ZNO2+1.)
      PV(6,11)=ZNO2+1.
      PV(5,21)=0.
      PV(6,21)=0.
      PV(8,21)=1.
      PV(5,31)=0.
      PV(6,31)=0.
      PV(8,31)=1.
C
      PV(5,12)=ATOMAS(ATNO2+2.,ZNO2+1.)
      PV(6,12)=ZNO2+1.
      PV(5,22)=RMASS(16)
      PV(6,22)=0.
      PV(8,22)=16.
      PV(5,32)=0.
      PV(6,32)=0.
      PV(8,32)=1.
C
      PV(5,13)=ATOMAS(ATNO2+2.,ZNO2   )
      PV(6,13)=ZNO2
      PV(5,23)=RMASS(14)
      PV(6,23)=1.
      PV(8,23)=14.
      PV(5,33)=0.
      PV(6,33)=0.
      PV(8,33)=1.
C
      PV(5,14)=ATOMAS(ATNO2+1.,ZNO2   )
      PV(6,14)=ZNO2
      PV(5,24)=RMASS(30)
      PV(6,24)=1.
      PV(8,24)=30.
      PV(5,34)=0.
      PV(6,34)=0.
      PV(8,34)=1.
C
      PV(5,15)=PV(5,2)
      PV(6,15)=PV(6,2)
      PV(5,25)=RMASS(31)
      PV(6,25)=1.
      PV(8,25)=31.
      PV(5,35)=0.
      PV(6,35)=0.
      PV(8,35)=1.
C
      PV(5,16)=ATOMAS(ATNO2-1.,ZNO2-1.)
      PV(6,16)=ZNO2-1.
      PV(5,26)=RMASS(32)
      PV(6,26)=2.
      PV(8,26)=32.
      PV(5,36)=0.
      PV(6,36)=0.
      PV(8,36)=1.
C
      PV(5,17)=ATOMAS(ATNO2+1.,ZNO2+1.)
      PV(6,17)=ZNO2+1.
      PV(5,27)=PV(5,22)
      PV(6,27)=0.
      PV(8,27)=16.
      PV(5,37)=PV(5,22)
      PV(6,37)=0.
      PV(8,37)=16.
C
      PV(5,18)=PV(5,14)
      PV(6,18)=PV(6,14)
      PV(5,28)=PV(5,22)
      PV(6,28)=0.
      PV(8,28)=16.
      PV(5,38)=PV(5,23)
      PV(6,38)=1.
      PV(8,38)=14.
C
      PV(5,19)=ATOMAS(ATNO2+1.,ZNO2-1.)
      PV(6,19)=ZNO2-1.
      PV(5,29)=PV(5,23)
      PV(6,29)=1.
      PV(8,29)=14.
      PV(5,39)=PV(5,23)
      PV(6,39)=1.
      PV(8,39)=14.
C
      NOPT=30
      GOTO 70
   60 PV(5,11)=ATOMAS(ATNO2+4.,ZNO2+2.)
      PV(6,11)=ZNO2+2.
      PV(5,21)=0.
      PV(6,21)=0.
      PV(8,21)=1.
      PV(5,31)=0.
      PV(6,31)=0.
      PV(8,31)=1.
C
      PV(5,12)=ATOMAS(ATNO2+3.,ZNO2+2.)
      PV(6,12)=ZNO2+2.
      PV(5,22)=RMASS(16)
      PV(6,22)=0.
      PV(8,22)=16.
      PV(5,32)=0.
      PV(6,32)=0.
      PV(8,32)=1.
C
      PV(5,13)=ATOMAS(ATNO2+3.,ZNO2+1.)
      PV(6,13)=ZNO2+1.
      PV(5,23)=RMASS(14)
      PV(6,23)=1.
      PV(8,23)=14.
      PV(5,33)=0.
      PV(6,33)=0.
      PV(8,33)=1.
C
      PV(5,14)=ATOMAS(ATNO2+2.,ZNO2+1.)
      PV(6,14)=ZNO2+1.
      PV(5,24)=RMASS(30)
      PV(6,24)=1.
      PV(8,24)=30.
      PV(5,34)=0.
      PV(6,34)=0.
      PV(8,34)=1.
C
      PV(5,15)=ATOMAS(ATNO2+1.,ZNO2+1.)
      PV(6,15)=ZNO2+1.
      PV(5,25)=RMASS(31)
      PV(6,25)=1.
      PV(8,25)=31.
      PV(5,35)=0.
      PV(6,35)=0.
      PV(8,35)=1.
C
      PV(5,16)=PV(5,2)
      PV(6,16)=PV(6,2)
      PV(5,26)=RMASS(32)
      PV(6,26)=2.
      PV(8,26)=32.
      PV(5,36)=0.
      PV(6,36)=0.
      PV(8,36)=1.
C
      PV(5,17)=ATOMAS(ATNO2+2.,ZNO2+2.)
      PV(6,17)=ZNO2+2.
      PV(5,27)=PV(5,22)
      PV(6,27)=0.
      PV(8,27)=16.
      PV(5,37)=PV(5,22)
      PV(6,37)=0.
      PV(8,37)=16.
C
      PV(5,18)=PV(5,14)
      PV(6,18)=PV(6,14)
      PV(5,28)=PV(5,22)
      PV(6,28)=0.
      PV(8,28)=16.
      PV(5,38)=PV(5,23)
      PV(6,38)=1.
      PV(8,38)=14.
C
      PV(5,19)=ATOMAS(ATNO2+2.,ZNO2   )
      PV(6,19)=ZNO2
      PV(5,29)=PV(5,23)
      PV(6,29)=1.
      PV(8,29)=14.
      PV(5,39)=PV(5,23)
      PV(6,39)=1.
      PV(8,39)=14.
C
      NOPT=40
   70 QV     =EK+PV(5,2)+PV(5,1)
      TC     =   PV(6,2)+PV(6,1)
      QVAL(1)=QV - PV(5,11)
      TCH (1)=TC - PV(6,11)
      QVAL(2)=QV - PV(5,12) - PV(5,22)
      TCH (2)=TC - PV(6,12) - PV(6,22)
      QVAL(3)=QV - PV(5,13) - PV(5,23)
      TCH (3)=TC - PV(6,13) - PV(6,23)
      QVAL(4)=QV - PV(5,14) - PV(5,24)
      TCH (4)=TC - PV(6,14) - PV(6,24)
      QVAL(5)=QV - PV(5,15) - PV(5,25)
      TCH (5)=TC - PV(6,15) - PV(6,25)
      QVAL(6)=QV - PV(5,16) - PV(5,26)
      TCH (6)=TC - PV(6,16) - PV(6,26)
      QVAL(7)=QV - PV(5,17) - PV(5,27) - PV(5,37)
      TCH (7)=TC - PV(6,17) - PV(6,27) - PV(6,37)
      QVAL(8)=QV - PV(5,18) - PV(5,28) - PV(5,38)
      TCH (8)=TC - PV(6,18) - PV(6,28) - PV(6,38)
      QVAL(9)=QV - PV(5,19) - PV(5,29) - PV(5,39)
      TCH (9)=TC - PV(6,19) - PV(6,29) - PV(6,39)
   74 QV = 0
      IF(IREC.EQ.2) QVAL(1)=0.
      IF(IPART.NE.16) GOTO 75
      CALL GRNDM(RNDM,2)
      IF(RNDM(1).GT.((ATNO2-1.)/230.)**2) QVAL(1)=0.
      EKA=7.9254/ATNO2
      IF(RNDM(2).LT.EK/EKA) GOTO 75
      QVAL(3)=0.
      QVAL(4)=0.
      QVAL(5)=0.
      QVAL(6)=0.
      QVAL(9)=0.
   75 DO 71 I=1,9
      IF(PV(5,10+I).LT.0.5) QVAL(I)=0.
      IF(QVAL(I).LT.0.    ) QVAL(I)=0.
      IF(ABS(TCH(I)-0.1).GT.0.5 ) QVAL(I)=0.
      QV=QV+QVAL(I)
   71 CONTINUE
      CALL GRNDM(RNDM,1)
      RAN=RNDM(1)
      QV1=0.
      DO 72 I=1,9
      IF(QVAL(I).EQ.0.) GOTO 72
      QV1=QV1+QVAL(I)/QV
      IF(RAN.LE.QV1) GOTO 73
   72 CONTINUE
C** REACTION KINEMATICALLY NOT POSSIBLE
      NOPT=0
      GO TO 9999
   73 NOPT=NOPT+I
      PV(5,3)=PV(5,10+I)
      PV(6,3)=PV(6,10+I)
      PV(8,3)=0.
      PV(5,4)=PV(5,20+I)
      PV(6,4)=PV(6,20+I)
      PV(8,4)=PV(8,20+I)
      PV(5,5)=PV(5,30+I)
      PV(6,5)=PV(6,30+I)
      PV(8,5)=PV(8,30+I)
      NT=2
      RAN=EK*10.
      IF(RAN.GT.0.5) RAN=0.5
      CALL GRNDM(RNDM,1)
      IF(RNDM(1).LT.RAN) NT=3
      IF(MOD(NOPT,10).GE.7) NT=3
C** CALCULATE CMS ENERGY
   80 PV(4,2)=PV(5,2)
      CALL ADD(1,2,MXGKPV)
      PV(1,MXGKPV)=-PV(1,MXGKPV)
      PV(2,MXGKPV)=-PV(2,MXGKPV)
      PV(3,MXGKPV)=-PV(3,MXGKPV)
C** SET QUANTITIES FOR PHASE SPACE ROUTINE IN CMS
      TECM=PV(5,MXGKPV)
      NPG=NT
      KGENEV=1
      DO 81 I=1,NPG
   81 AMASS(I)=PV(5,2+I)
C --- Invoke double precision version of the phase space package ---
      CALL PHPNUC
      DO 83 I=1,NPG
      DO 82 J=1,4
   82 PV(J,2+I)=PCM(J,I)
C** TRANSFORM INTO LAB.SYSTEM
      CALL LOR(2+I,MXGKPV,2+I)
      PV(7,2+I)=TOF
   83 CONTINUE
C** SET CHARGES AND PARTICLE INDEX FOR LOW MASS FRAGMENTS
      IF (ABS(PV(5,3)-RMASS(14)) .LT. 0.0001) GO TO 84
      IF (ABS(PV(5,3)-RMASS(16)) .LT. 0.0001) GO TO 85
      IF (ABS(PV(5,3)-RMASS(30)) .LT. 0.0001) GO TO 86
      IF (ABS(PV(5,3)-RMASS(31)) .LT. 0.0001) GO TO 87
      IF (ABS(PV(5,3)-RMASS(32)) .LT. 0.0001) GO TO 88
      GOTO 89
   84 PV(6,3)=1.
      PV(8,3)=14.
      GOTO 89
   85 PV(6,3)=0.
      PV(8,3)=16.
      GOTO 89
   86 PV(6,3)=1.
      PV(8,3)=30.
      GOTO 89
   87 PV(6,3)=1.
      PV(8,3)=31.
      GOTO 89
   88 PV(6,3)=2.
      PV(8,3)=32.
   89 NTT=2+NT
      DO 90 I=1,NTT
      IPP=IFIX(PV(8,I)+0.01)
      IF(IPP.EQ.0) GOTO 90
      EK=PV(4,I)-PV(5,I)
      IF(I.LT.3) GOTO 92
      IF(IPP.LT.30) GOTO 92
      CALL GRNDM(RNDM,1)
      EK=EK*0.5*RNDM(1)
   92 IF(EK.LT.1.E-6) EK=1.E-6
      PV(5,I)=RMASS(IPP)
      PV(4,I)=EK+PV(5,I)
      P=SQRT(ABS(PV(4,I)**2-PV(5,I)**2))
      PP=SQRT(PV(1,I)**2+PV(2,I)**2+PV(3,I)**2)
      IF(PP.GT.1.E-6) GOTO 91
      CALL GRNDM(RNDM,2)
      PHINVE=TWPI*RNDM(1)
      COST=-1.+2.*RNDM(2)
      IF (COST .LE. -1.) COST=-1.
      IF (COST .GE.  1.) COST= 1.
      RTHNVE=ACOS(COST)
      PV(1,I)=SIN(RTHNVE)*COS(PHINVE)
      PV(2,I)=SIN(RTHNVE)*SIN(PHINVE)
      PV(3,I)=COS(RTHNVE)
      PP=1.
   91 PV(1,I)=PV(1,I)*P/PP
      PV(2,I)=PV(2,I)*P/PP
      PV(3,I)=PV(3,I)*P/PP
   90 CONTINUE
      IF(.NOT.NPRT(4)) GOTO 100
      WRITE(NEWBCD,1000) XEND,YEND,ZEND,IND,NOPT
 1000 FORMAT(1H ,'Nuclear reaction at (X,Y,Z) ',3(G12.5,1X)
     +,' Material ',I5,' NOPT ',I5)
      DO 95 I=1,NTT
         WRITE(NEWBCD,1001) I,(PV(J,I),J=1,10)
   95 CONTINUE
 1001 FORMAT(1H ,I3,1X,10(G10.3,1X))
  100 INTCT=INTCT+1.
C** SET INTERACTION MODE ACCORDING TO GHEISHA-CONVENTION
C** N-CAPTURE
      IF(PV(8,3).GT.0.) GOTO 110
      CALL SETCUR(4)
      NTK=NTK+1
      IF(NT.EQ.3) CALL SETTRK(5)
      GO TO 9999
 110  CONTINUE
      CALL SETCUR(4)
      NTK=NTK+1
      CALL SETTRK(3)
      IF(NT.EQ.3) CALL SETTRK(5)
      CALL SETTRK(3)
C
 9999 CONTINUE
      END
