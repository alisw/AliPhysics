*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:45  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.32  by  S.Giani
*-- Author :
      SUBROUTINE CGWIRE(CG,NTRAN,IFHIDE,NMAX,WIRE,ISHAPE,SHADE)
************************************************************************
*                                                                      *
*     Name: CGWIRE                                                     *
*     Author: E. Chernyaev                       Date:    15.04.88     *
*     Revised by: S.Giani                        Date:    24.04.91     *
*                                                                      *
*     Function: Transfer CG-object to WIRE-object                      *
*                                                                      *
*     References: CGTSTR, CGFVIS, CGWSOR                               *
*                                                                      *
*     Input:  CG - CG-object                                           *
*             NT - number of transformation to screen coordinates      *
*         IFHIDE - flag ( 0 - all edges,-1 - only not hidden)          *
*          NMAX  - max-length of WIRE array                            *
*                                                                      *
*     Output: WIRE - WIRE-object                                       *
*                    WIRE(1) - length of WIRE-object                   *
*                        = 0 if error in CG-object structure           *
*                        < 0 if no space in WIRE array or error        *
*                            in parameters                             *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/cggpar.inc"
#include "geant321/cgctra.inc"
*SG
#include "geant321/gcdraw.inc"
#include "geant321/gcgobj.inc"
#include "geant321/gcspee.inc"
*SG
      REAL      CG(*),WIRE(*)
      INTEGER   SHADE(*)
*-
      LLEP=ABS(LEP)
      WIRE(1)= 0.
      IF(LLEP.NE.1)SHADE(1)=0.
      MMM=1
*           T E S T   P A R A M E T E R S
      IF (NTRAN .LE. 0)                 GOTO 999
      IF (NTRAN .GT. NTMAX)             GOTO 999
      IF (NMAX  .LE. 0)                 GOTO 999
      CALL CGTSTR(CG,IREP)
      IF (IREP .LT. 0)                  GOTO 999
      NFACE  = CG(KCGNF)
*SG
      IF (NFACE .EQ. 0)THEN
         WIRE(KCGNF)=0
         GOTO 999
      ENDIF
*SG
*           I N I T I A L I S A T I O N
      NWIRE  = 0
      JCG    = LCGHEA
      JWR    = LCGHEA
*
**          L O O P   A L O N G   F A C E S
*SG
      NTIM=0
*SG
      DO 300 NF=1,NFACE
        JCGFAC = JCG
        NEDGE  = CG(JCG+KCGNE)
*           D E F I N E   V A L U E   O F  "IFVIS"
        IFVIS  = -1
        IF (IFHIDE .NE. 0)              GOTO 200
        JCG    = JCG + LCGFAC
        DO 100 NE=1,NEDGE
          IF (CG(JCG+KCGAE) .LT. 0.)    GOTO 200
          JCG   = JCG + LCGEDG
  100     CONTINUE
        IFVIS = 0
        GOTO 210
*
**           T E S T   F A C E   V I S I B I L I T Y
**           M O V E   E D G E S   T O   W I R E - O B J E C T
*
  200   CALL CGFARE(NTRAN,CG(JCGFAC+1),IFVIS,ISHAPE)
  210   IF (IFHIDE.NE.0 .AND. IFVIS.LT.0)       GOTO 290
        JCG    = JCGFAC
        AFACE  = CG(JCG+KCGAF)
        JCG    = JCG + LCGFAC
        DO 250 NE=1,NEDGE
          IF (CG(JCG+KCGAE) .LT. 0.)            AEDGE =-AFACE - 1
          IF (CG(JCG+KCGAE) .GE. 0.)            AEDGE = AFACE
          IF (AEDGE.LT.0. .AND. IFVIS.LT.0)     GOTO 240
          IF (NCLAS3+JWR+LCGEDG .GT. NMAX)             GOTO 220
          WIRE(JWR+KCGAE) = AEDGE
          WIRE(JWR+KCGX1) = CG(JCG+KCGX1)
          WIRE(JWR+KCGY1) = CG(JCG+KCGY1)
          WIRE(JWR+KCGZ1) = CG(JCG+KCGZ1)
          WIRE(JWR+KCGX2) = CG(JCG+KCGX2)
          WIRE(JWR+KCGY2) = CG(JCG+KCGY2)
          WIRE(JWR+KCGZ2) = CG(JCG+KCGZ2)
  220     NWIRE  = NWIRE + 1
          JWR    = JWR + LCGEDG
  240     JCG    = JCG + LCGEDG
  250     CONTINUE
        MMM=MMM+1
        IF(LLEP.NE.1)SHADE(MMM)=NWIRE
  290   JCG    = JCGFAC + LCGFAC + NEDGE*LCGEDG
  300   CONTINUE
*
**          D E L E T E  V E R Y   S M A L L   W I R E S
**          S E T   N E E D E D   D I R E C T I O N
**          S O R T   W I R E S
*
*SG
*        Finding the total number of words needed to build
*        the whole Wire Structure.
*
      NCLAS3=NCLAS3+JWR
      IF(NCLAS3.GE.NMAX)THEN
         KCGST=-8
         GOTO 999
      ENDIF
*SG
      IF (JWR .GT. NMAX)        GOTO 998
*SG
      LINFIL=IBITS(LINATT,13,3)
      IF(LINFIL.EQ.0)CALL CGWSOR(NWIRE,WIRE(KCGNF+1))
*SG
      WIRE(KCGSIZ) = LCGHEA + NWIRE*LCGEDG
      IF(LLEP.NE.1)SHADE(1)=MMM
      WIRE(KCGATT) = CG(KCGATT)
      WIRE(KCGNF)  = NWIRE
      GOTO 999
*
  998 WIRE(KCGSIZ) = NMAX - JWR
  999 RETURN
      END
