*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:50  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.29  by  S.Giani
*-- Author :
      SUBROUTINE GGVCHK (JVOM, IN, NVAR, LVAR)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    SUBR. GGVCHK (JVOM, IN, NVAR*, LVAR*)                       *
C.    *                                                                *
C.    *   Checks volume parameters of IN'th content of mother volume   *
C.    *    at address JVOM or, when IN = 0, of mother volume itself.   *
C.    *   Returns NVAR = 0, when no variable parameters                *
C.    *       or  NVAR =-1, when no variable parameters AND not a      *
C.    *        possible source of local development                    *
C.    *       otherwise,                                               *
C.    *        the list LVAR of the NVAR variable parameters positions *
C.    *                                                                *
C.    *   Called by : GGDVLP                                           *
C.    *   Author    : F.Bruyant                                        *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcshno.inc"
C.
      DIMENSION LVAR(*)
      INTEGER   LDP(7,12), NDP(12)
      SAVE LDP,NDP
C.
      DATA NDP/ 3, 4, 5, 7, 3, 3, 5, 5, 4, 3, 2, 2 /
      DATA LDP/ 1, 2, 3, 4*0,  1, 2, 3, 4, 3*0,  1, 2, 3, 4, 5, 2*0,
     +        1, 4, 5, 6, 8, 9, 10,  1, 2, 3, 4*0,  1, 2, 3, 4*0,
     +        1, 2, 3, 4, 5, 2*0,  1, 2, 3, 4, 5, 2*0,  1, 2, 3, 4, 3*0,
     +        1, 2, 3, 4*0,  2, 3, 5*0,  2, 3, 5*0 /
C.
C.    ------------------------------------------------------------------
*
      NVAR = 0
      NSP  = 0
      IF (IN.EQ.0) THEN
         NPAR = Q(JVOM+5)
         IF (NPAR.LE.0) THEN
            NIN = Q(JVOM+3)
            IF (NIN.LT.0) NVAR = -1
            GO TO 999
         ENDIF
         JPAR = JVOM +6
         ISH = Q(JVOM+2)
      ELSE
         JIN  = LQ(JVOM-IN)
         IVO  = Q(JIN+2)
         JVO  = LQ(JVOLUM-IVO)
         NPAR = Q(JVO+5)
         IF (NPAR.GT.0) THEN
            JPAR = JVO +6
         ELSE
            NIN  = Q(JVOM+3)
            IF (NIN.LT.0) THEN
               JPAR = 0
            ELSE
               JPAR = JIN +9
               NIN  = Q(JVO+3)
               IF (NIN.LT.0) THEN
                  NSP = -1
               ELSE
*
                  DO 9 IIN = 1, NIN
                     JJIN = LQ(JVO-IIN)
                     IIVO = Q(JJIN+2)
                     JJVO = LQ(JVOLUM-IIVO)
                     NNPAR= Q(JJVO+5)
                     IF (NNPAR.GT.0) THEN
                        JJPAR= JJVO + 6
                     ELSE
                        NSP=-1
                        GO TO 10
                     ENDIF
*
                     IISH = Q(JJVO+2)
                     NND  = NDP(IISH)
                     IF (IISH.LE.10) THEN
*
*   *                  Volumes other than PGON or PCON
*
                        DO 6 ID = 1,NND
                           IIDP = LDP(ID,IISH)
                           IF (Q(JJPAR+IIDP).LT.0.) THEN
                              NSP = -1
                              GO TO 10
                           ENDIF
    6                   CONTINUE
*
                     ELSE IF (IISH.LE.12) THEN
*
*   *                  PGON and PCON volumes
*
                        IPZ  = 15 -IISH
                        JJPAR = JJPAR +IPZ
                        NNZ   = 2*Q(JJPAR)
                        DO 8 IZ = 1, NNZ
                           DO 7 ID = 1, NND
                              IIDP = LDP(ID,IISH)
                              IF (Q(JJPAR+IIDP).LT.0.) THEN
                                 NSP = -1
                                 GO TO 10
                              ENDIF
    7                      CONTINUE
    8                   CONTINUE
                     ENDIF
    9             CONTINUE
*
               ENDIF
*
            ENDIF
         ENDIF
  10     ISH = Q(JVO+2)
      ENDIF
*        Shape 28 will not be supported in the future
      IF (ISH.EQ.28.OR.ISH.EQ.13.OR.ISH.EQ.NSCTUB.OR.ISH.EQ.14)THEN
         NVAR = 0
         GO TO 999
      ENDIF
*
      ND = NDP(ISH)
*
      IF (JPAR.LE.0) THEN
         IF (ISH.LE.10) THEN
*
*  *       Volumes other PGON or PCON
*
            NVAR = ND
            DO 11 ID = 1, ND
   11       LVAR(ID) = LDP(ID,ISH)
         ELSE IF (ISH.LE.12) THEN
*
*  *       PGON and PCON volumes
*
            IPZ  = 15 - ISH
            NVAR = 6
            DO 12 ID = 1, NVAR
   12       LVAR(ID) = IPZ + ID
         ELSE
*
            GO TO 900
         ENDIF
         GO TO 999
      ENDIF
*
      IF (ISH.LE.10) THEN
*
*   *   Volumes other than PGON or PCON
*
         DO 19 ID = 1,ND
            IDP = LDP(ID,ISH)
            IF (Q(JPAR+IDP).GE.0.) GO TO 19
            NVAR = NVAR +1
            LVAR(NVAR) = IDP
   19    CONTINUE
*
      ELSE IF (ISH.LE.12) THEN
*
*   *   PGON and PCON volumes
*
         IPZ  = 15 -ISH
         JPAR = JPAR +IPZ
         NZ   = 2*Q(JPAR)
         INC  = IPZ
         DO 29 IZ = 1,NZ
            DO 28 ID = 1,ND
               IDP = LDP(ID,ISH)
               IF (Q(JPAR+IDP).GE.0.) GO TO 28
               NVAR = NVAR +1
               LVAR(NVAR) = INC +IDP
   28       CONTINUE
            INC = INC +3
   29    CONTINUE
      ELSE
*
         GO TO 900
      ENDIF
      IF (NVAR.EQ.0) NVAR = NSP
      GO TO 999
*
 900  WRITE (CHMAIL, 1001) ISH
      CALL GMAIL (0, 0)
 1001 FORMAT (' GGVCHK : No code for shape ISH=',I5)
*                                                             END GGVCHK
 999  END
