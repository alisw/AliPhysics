*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:08  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.20  by  S.Giani
*-- Author :
      SUBROUTINE GCDERR (ICD,ERP,ERS)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       ROUTINE TO CALCULATE THE ERROR ON THE CURRENT DIVISION   *
C.    *       INFORMATION AS OBTAINED BY "GCDRIF"                      *
C.    *                                                                *
C.    *       INPUT:  ICD = DIGITIZED CURRENT DIVISION INFORMATION     *
C.    *                     ( 0 ... 1000 )                             *
C.    *               ERP = VARIANCE OF GAUSSIAN DISTRIBUTED PEDESTAL  *
C.    *                     ERRORS ON THE MEASURED PULSE HEIGHTS       *
C.    *                     RELATIVE TO THE SUM OF THE PULSE HEIGHTS   *
C.    *               ERS = VARIANCE OF GAUSSIAN DISTRIBUTED SLOPE     *
C.    *                     ERRORS ON THE MEASURED PULSE HEIGHTS       *
C.    *                     RELATIVE TO EACH PULSE HEIGHT              *
C.    *       OUTPUT: ICD = VALUE UPDATED ACCORDING TO RANDOM ERRORS   *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GUDIGI                               *
C.    *       Author    D.Mitaroff  *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
      COMMON  /GCYDR/ ZL
      DIMENSION RNDM(4)
C
      Z = ICD
C
C---- CALCULATE PEDESTAL ERRORS FOR VARIANCE BEING 1.
      CALL GRNDM(RNDM,4)
*          call rannor(eps1,eps2)
      RY=RNDM(1)
      RZ=RNDM(2)
      RX=6.283185*RZ
      A1=SQRT(-2.*LOG(RY))
      EPS1=A1*SIN(RX)
      EPS2=A1*COS(RX)
      DZP = - EPS1 * Z + EPS2 * (ZL - Z)
C
C---- CALCULATE SLOPE ERRORS FOR VARIANCE BEING 1.
*          call rannor(eps1,eps2)
      RY=RNDM(3)
      RZ=RNDM(4)
      RX=6.283185*RZ
      A1=SQRT(-2.*LOG(RY))
      EPS1=A1*SIN(RX)
      EPS2=A1*COS(RX)
      DZS = (EPS2 - EPS1) * Z * (ZL - Z) / ZL
C
C---- CALCULATE NEW VALUE OF ICD FOR VARIANCES ERP, ERS.
      Z = Z + ERP * DZP + ERS * DZS
      IF (Z .LT. 0.)   Z = 0.
      IF (Z .GT. ZL)   Z = ZL
C
      ICD = Z
C
      END
