*
* $Id$
*
* $Log$
* Revision 1.4  1997/10/17 10:11:47  mclareni
* Negative square root protection
*
* Revision 1.3  1996/05/03 07:31:39  cernlib
* Replace in Max( .., 0.0) the 0.0 by ZERO to have the proper variable type.
*
* Revision 1.2  1996/04/26 12:19:59  ravndal
* neg. SQRT protection
*
* Revision 1.1.1.1  1995/10/24 10:19:58  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*$ CREATE STALIN.FOR
*COPY STALIN
*                                                                      *
*=== stalin ===========================================================*
*                                                                      *
      SUBROUTINE STALIN
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
       PARAMETER (ZERO=0)
*----------------------------------------------------------------------*
*                                                                      *
*     STAbility LINe calculation:                                      *
*                                                                      *
*     Created on 04 december 1992  by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on 04-dec-92     by    Alfredo Ferrari               *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*
#include "geant321/isotop.inc"
*
      DIMENSION ZNORM (260)
*  +-------------------------------------------------------------------*
*  |
      DO 1000 IZ=1,100
         DO 500 J=1,2
            ASTLIN (J,IZ) = 0.D+00
  500    CONTINUE
 1000 CONTINUE
*  |
*  +-------------------------------------------------------------------*
*  +-------------------------------------------------------------------*
*  |
      DO 2000 IA=1,260
         ZNORM (IA) = 0.D+00
         DO 1500 J=1,2
            ZSTLIN (J,IA) = 0.D+00
 1500    CONTINUE
 2000 CONTINUE
*  |
*  +-------------------------------------------------------------------*
*  +-------------------------------------------------------------------*
*  |  Loop on the Atomic Number
      DO 3000 IZ=1,100
*  |  +----------------------------------------------------------------*
*  |  |    Loop on the stable isotopes
         DO 2500 IS = ISONDX (1,IZ), ISONDX (2,IZ)
            IA = ISOMNM (IS)
            ASTLIN (1,IZ) = ASTLIN (1,IZ) + ABUISO (IS) * IA
            ASTLIN (2,IZ) = ASTLIN (2,IZ) + ABUISO (IS) * IA**2
            ZNORM    (IA) = ZNORM (IA) + ABUISO (IS)
            ZSTLIN (1,IA) = ZSTLIN (1,IA) + ABUISO (IS) * IZ
            ZSTLIN (2,IA) = ZSTLIN (2,IA) + ABUISO (IS) * IZ**2
 2500    CONTINUE
*  |  |
*  |  +----------------------------------------------------------------*
*  |  Normalize and print A_stab versus Z data:
         ASTLIN (2,IZ) = MAX ( SQRT (
     &                   MAX ( ASTLIN(2,IZ)-ASTLIN(1,IZ)**2,ZERO)),
     &                   HLFHLF )
 3000 CONTINUE
*  |
*  +-------------------------------------------------------------------*
*  +-------------------------------------------------------------------*
*  |  Normalize and print Z_stab versus A data:
      TMP10 = 1.D-10
      TMP03 = 0.3D+00
      DO 4000 IA=1,260
         ZSTLIN (1,IA) = ZSTLIN (1,IA) / MAX ( ZNORM (IA), TMP10 )
         ZSTLIN (2,IA) = ZSTLIN (2,IA) / MAX ( ZNORM (IA), TMP10 )
         ZSTLIN (2,IA) = MAX ( ZSTLIN (2,IA), ZSTLIN (1,IA)**2 )
         IF ( ZNORM (IA) .GT. ANGLGB )
     &     ZSTLIN (2,IA) = MAX ( SQRT (
     &                     MAX ( ZSTLIN(2,IA)-ZSTLIN(1,IA)**2,ZERO)),
     &                     TMP03 )
 4000 CONTINUE
*  |
*  +-------------------------------------------------------------------*
*  +-------------------------------------------------------------------*
*  |  Normalize and print Z_stab versus A data:
      DO 5000 IA=1,260
         IF ( ZNORM (IA) .LE. ANGLGB ) THEN
            DO 4200 JA = IA-1,1,-1
               IF ( ZNORM (JA) .GT. ANGLGB ) THEN
                  IA1 = JA
                  GO TO 4300
               END IF
 4200       CONTINUE
 4300       CONTINUE
            DO 4400 JA = IA+1,260
               IF ( ZNORM (JA) .GT. ANGLGB ) THEN
                  IA2 = JA
                  GO TO 4500
               END IF
 4400       CONTINUE
            IA2 = IA1
            IA1 = IA1 - 1
 4500       CONTINUE
            DIAIA1 = IA-IA1
            DIA2IA = IA2-IA1
            ZSTLIN (1,IA) = DIAIA1 / DIA2IA
     &                    * ( ZSTLIN (1,IA2) - ZSTLIN (1,IA1) )
     &                    + ZSTLIN (1,IA1)
            ZSTLIN (2,IA) = DIAIA1 / DIA2IA
     &                    * ( ZSTLIN (2,IA2) - ZSTLIN (2,IA1) )
     &                    + ZSTLIN (2,IA1)
         END IF
         IZ = MIN ( 100, NINT (ZSTLIN(1,IA)) )
         ATOZ = IZ / ASTLIN (1,IZ)
         ZSTLIN (2,IA) = MAX ( ZSTLIN (2,IA), ATOZ * ASTLIN (2,IZ) )
 5000 CONTINUE
*  |
*  +-------------------------------------------------------------------*
      RETURN
      END
