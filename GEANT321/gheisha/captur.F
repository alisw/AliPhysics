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
      SUBROUTINE CAPTUR(NOPT)
C
C *** ROUTINE FOR CAPTURE OF NEUTRAL BARYONS ***
C *** NVE 04-MAR-1988 CERN GENEVA ***
C
C ORIGIN : H.FESEFELDT (02-DEC-1986)
C
#include "geant321/s_defcom.inc"
      DIMENSION RNDM(3)
C
      NOPT=1
      IER(81)=IER(81)+1
      PV(1,1)=PX*P
      PV(2,1)=PY*P
      PV(3,1)=PZ*P
      PV(4,1)=EN
      PV(5,1)=ABS(AMAS)
      PV(6,1)=NCH
      PV(7,1)=TOF
      PV(8,1)=IPART
      PV(9,1)=0.
      PV(10,1)=USERW
      ND=IND+1
      PV(1,2)=0.
      PV(2,2)=0.
      PV(3,2)=0.
      PV(4,2)=ATOMAS(ATNO(ND),ZNO(ND))
      PV(5,2)=PV(4,2)
      PV(6,2)=ZNO(ND)
      PV(7,2)=TOF
      PV(8,2)=0.
      PV(9,2)=0.
      PV(10,2)=0.
      CALL ADD(1,2,MXGKPV)
      PV(1,MXGKPV)=-PV(1,MXGKPV)
      PV(2,MXGKPV)=-PV(2,MXGKPV)
      PV(3,MXGKPV)=-PV(3,MXGKPV)
      CALL NORMAL(RAN)
      P=0.0065+RAN*0.0010
      CALL GRNDM(RNDM,3)
      COST=-1.+RNDM(1)*2.
      SINT=SQRT(ABS(1.-COST*COST))
      PHI=TWPI*RNDM(2)
      PV(1,3)=P*SINT*SIN(PHI)
      PV(2,3)=P*SINT*COS(PHI)
      PV(3,3)=P*COST
      PV(4,3)=P
      PV(5,3)=0.
      PV(6,3)=0.
      PV(8,3)=1.
      PV(9,3)=0.
      PV(10,3)=0.
      RAN=RNDM(3)
      TOF=TOF-480.*LOG(RAN)
      PV(7,3)=TOF
      CALL LOR(3,MXGKPV,3)
      NT=3
      XP=0.008-P
      IF(XP.LT.0.) GOTO 9
      NT=4
      CALL GRNDM(RNDM,2)
      COST=-1.+RNDM(1)*2.
      SINT=SQRT(ABS(1.-COST*COST))
      PHI=TWPI*RNDM(2)
      PV(1,4)=XP*SINT*SIN(PHI)
      PV(2,4)=XP*SINT*COS(PHI)
      PV(3,4)=XP*COST
      PV(4,4)=XP
      PV(5,4)=0.
      PV(6,4)=0.
      PV(7,4)=TOF
      PV(8,4)=1.
      PV(9,4)=0.
      PV(10,4)=0.
      CALL LOR(4,MXGKPV,4)
    9 INTCT=INTCT+1.
      CALL SETCUR(3)
      NTK=NTK+1
      IF(NT.EQ.4) CALL SETTRK(4)
      IF(NPRT(4))
     *WRITE(NEWBCD,1002) XEND,YEND,ZEND,EN,XP
1002  FORMAT(' *CAPTUR* NEUTRON CAPTURE   POSITION',3(2X,F8.2),2X,
     * 'PHOTON ENERGIES',2X,2F8.4)
C
      END
