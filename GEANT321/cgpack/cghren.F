*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:44  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.31  by  S.Giani
*-- Author :
      SUBROUTINE CGHREN(NT,NOLD,NNEW,NFACE,XYZ,IP,IFACE,NXYZ)
************************************************************************
*                                                                      *
*     Name: CGHREN                                                     *
*     Author: E. Chernyaev                       Date:    04.08.88     *
*                                                Revised:              *
*                                                                      *
*     Function: Transform coordinates to screen system and             *
*               make renumeration                                      *
*                                                                      *
*     References: none                                                 *
*                                                                      *
*     Input:    NT - number of transformation                          *
*             NOLD - number of old nodes                               *
*             NNEW - number of new nodes                               *
*            NFACE - number of new faces                               *
*         XYZ(3,*) - node coordinates                                  *
*          IP(2,*) - work array for renumbering                        *
*         IFACE(*) - faces                                             *
*                                                                      *
*     Output: NXYZ - total number of nodes after renumeration          *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/cgdelt.inc"
#include "geant321/cgctra.inc"
      REAL      XYZ(3,*),SXYZ(3)
*SG
      INTEGER  IP(2,*),IFACE(*)
*SG
*
**          T R A N S F O R M   T O   S C R E E N   COORDINATES
**          D I S C R E T I S A T I O N
*
      DO 200 I=1,NNEW
        DO 100 K=1,3
          SXYZ(K) = TSCRN(1,K,NT)*XYZ(1,I) + TSCRN(2,K,NT)*XYZ(2,I) +
     +              TSCRN(3,K,NT)*XYZ(3,I) + TSCRN(4,K,NT)
*          IF (SXYZ(K) .GE. 0.) KK = (SXYZ(K) + DESCR) * DELSCR
*          IF (SXYZ(K) .LT. 0.) KK = (SXYZ(K) - DESCR) * DELSCR
*          SXYZ(K) = KK * EESCR
  100     CONTINUE
        IP(1,I) = I
        XYZ(1,I) = SXYZ(1)
        XYZ(2,I) = SXYZ(2)
        XYZ(3,I) = SXYZ(3)
  200   CONTINUE
*
**         S H E L L   S O R T   O F   C O O R D N A T E S
*
*      ISTEP  = 1
*  290 ISTEP  = ISTEP*3 + 1
*      IF (ISTEP*2 .LT. NNEW)    GOTO 290
*  300 ISTEP  = ISTEP/3
*      DO 500 M=1,NNEW-ISTEP
*        IF(XYZ(1,M)-XYZ(1,M+ISTEP))     500,310,350
*  310   IF(XYZ(2,M)-XYZ(2,M+ISTEP))     500,320,350
*  320   IF(XYZ(3,M)-XYZ(3,M+ISTEP))     500,500,350
**
*  350   SXYZ(1)    = XYZ(1,M+ISTEP)
*        SXYZ(2)    = XYZ(2,M+ISTEP)
*        SXYZ(3)    = XYZ(3,M+ISTEP)
*        IPCUR      =  IP(1,M+ISTEP)
*        I          = M
*  400   XYZ(1,I+ISTEP) = XYZ(1,I)
*        XYZ(2,I+ISTEP) = XYZ(2,I)
*        XYZ(3,I+ISTEP) = XYZ(3,I)
*        IP (1,I+ISTEP) =  IP(1,I)
*        I          = I - ISTEP
*        IF (I .LE. 0)           GOTO 450
*        IF (XYZ(1,I)-SXYZ(1))   450,410,400
*  410   IF (XYZ(2,I)-SXYZ(2))   450,420,400
*  420   IF (XYZ(3,I)-SXYZ(3))   450,450,400
*  450   XYZ(1,I+ISTEP) = SXYZ(1)
*        XYZ(2,I+ISTEP) = SXYZ(2)
*        XYZ(3,I+ISTEP) = SXYZ(3)
*        IP (1,I+ISTEP) = IPCUR
**
*  500   CONTINUE
*      IF (ISTEP .NE. 1)          GOTO 300
*
**          N O D E   R E N U M E R A T I O N
*
      NN           = 1
      NIP11=IP(1,1)
      IP(2,NIP11)= NN + NOLD
      DO 650 I=2,NNEW
*        IF (XYZ(1,I) .NE. XYZ(1,I-1)) GOTO 610
*        IF (XYZ(2,I) .NE. XYZ(2,I-1)) GOTO 610
*        IF (XYZ(3,I) .EQ. XYZ(3,I-1)) GOTO 620
  610   NN = NN + 1
        XYZ(1,NN) = XYZ(1,I)
        XYZ(2,NN) = XYZ(2,I)
        XYZ(3,NN) = XYZ(3,I)
  620   NIP=IP(1,I)
        IP(2,NIP)= NN + NOLD
  650   CONTINUE
*
**          S E T   N E W   N O D E   N U M B E R S   I N   F A C E S
*
      JF     = 1
      DO 800 NF=1,NFACE
        NEDGE  = IFACE(JF)
        DO 700 NE=1,NEDGE*2
          NIF=IFACE(JF+NE)
          IFACE(JF+NE) = IP(2,NIF)
  700     CONTINUE
        JF = JF + 1 + NEDGE*2
  800   CONTINUE
*
      NXYZ   = NN + NOLD
      RETURN
      END
