*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:58  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 25/07/94  15.04.20  by  S.Giani
*-- Author :
*$ CREATE NUDISV.FOR
*COPY NUDISV
*
*=== nudisv ===========================================================*
*
      INTEGER FUNCTION NUDISV ( ANU, KB, EXPLAM, ASEASQ, APOWER, PRZERO)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*                                                                      *
*     Last change  on  05-may-92    by  Alfredo Ferrari, INFN-Milan    *
*                                                                      *
*     Nudist91:                                                        *
*      multiplicity distribution of chains produced in hadron-nucleus  *
*      collisions                                                      *
*                                                                      *
*      Input variables:                                                *
*                       anu    = average number of collisions          *
*                       kb     = baryon number of the projectile       *
*                       explam = exp ( - (anu-1) )                     *
*      Output variables:                                               *
*                       nudisv = number of high energy collisions      *
*----------------------------------------------------------------------*
*
      PARAMETER ( INMAX  = 30 )
      PARAMETER ( IAMAX  = 13 )
      PARAMETER ( ANUMAX = 6.0D+00 )
      PARAMETER ( ANUMED = 0.5D+00 * ( 3.75D+00 + ANUMAX ) )
*
C      AB,(KB.NE.0) BARYON NUCLEUS
C
      DIMENSION AB (INMAX,IAMAX), ANUB (IAMAX), ANUM (9),
     &          ANSESQ (IAMAX), ANPOWR (IAMAX), NULAST (IAMAX)
      DIMENSION AAB (IAMAX), DIST (0:INMAX), IEX (3,5)
      REAL RNDM(1)
      LOGICAL LFIRST
      SAVE LFIRST, ANUB, ANUM, AB, ANSESQ, NULAST, ANPOWR
      DATA ANUB/1.D0,1.1D0,1.25D0,1.45D0,1.65D0,1.74D0,2.09D0,2.66D0,
     &          3.0D0,3.1D0,3.75D0,ANUMED,ANUMAX/
      DATA ANUM/1.D0,1.45D0,1.52D0,1.77D0,2.19D0,2.42D0,2.47D0,
     *2.90D0,5.D0/
      DATA AB  /   1.D+00, 29*0.D+00, 90*0.D+00,
     & 0.5771D+00, 0.2565D+00, 0.1159D+00, 0.0397D+00,
     & 0.0095D+00, 0.0013D+00, 0.0001D+00, 23*0.D+00,
     & 0.5457D+00, 0.2548D+00, 0.1303D+00, 0.0515D+00,
     & 0.0146D+00, 0.0028D+00, 0.0003D+00, 23*0.D+00,
     & 0.4462D+00, 0.2445D+00, 0.1599D+00, 0.0919D+00, 0.0406D+00,
     & 0.0132D+00, 0.0032D+00, 0.0005D+00, 0.0001D+00, 21*0.D+00,
     & 0.3387D+00, 0.2105D+00, 0.1670D+00, 0.1291D+00, 0.0849D+00,
     & 0.0441D+00, 0.0180D+00, 0.0059D+00, 0.0014D+00, 0.0003D+00,
     & 0.0001D+00, 19*0.D+00,
     & 0.2959D+00, 0.1885D+00, 0.1673D+00, 0.1344D+00, 0.1044D+00,
     & 0.0653D+00, 0.0336D+00, 0.0144D+00, 0.0046D+00, 0.0012D+00,
     & 0.0003D+00, 19*0.D+00,
     & 0.2835D+00, 0.1818D+00, 0.1572D+00, 0.1357D+00, 0.1082D+00,
     & 0.0715D+00, 0.0380D+00, 0.0161D+00, 0.0057D+00, 0.0018D+00,
     & 0.0004D+00, 0.0001D+00, 18*0.D+00,
     & 0.2247D+00, 0.1481D+00, 0.1342D+00, 0.1328D+00, 0.1211D+00,
     & 0.0984D+00, 0.0700D+00, 0.0394D+00, 0.0198D+00, 0.0077D+00,
     & 0.0028D+00, 0.0008D+00, 0.0002D+00, 17*0.D+00,  60*0.D+00/
      DATA IEX / 1, 3, 6, 1, 2, 3, 2, 4, 5, 10, 12, 11, 11, 13, 12 /
      DATA LFIRST / .TRUE. /
*
*  +-------------------------------------------------------------------*
*  |      First call: perform the initialization
      IF ( LFIRST ) THEN
         LFIRST = .FALSE.
         IPOWER = NINT (-APOWER)
*  |  +----------------------------------------------------------------*
*  |  |
         DO 100 IA = 1, IAMAX
            AAB (IA) = 0.D+00
*  |  |  +-------------------------------------------------------------*
*  |  |  |          This loop computes the normalization factor
            DO 80 IN = 1, INMAX
               AAB (IA) = AAB (IA) + AB (IN,IA)
   80       CONTINUE
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF ( AAB (IA) .GT. 0.D+00 ) THEN
               ANUAVE = 0.D+00
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               DO 90 IN = 1, INMAX
                  AB (IN,IA) = AB (IN,IA) / AAB (IA)
                  ANUAVE = ANUAVE + IN * AB (IN,IA)
   90          CONTINUE
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
               ANUB (IA) = ANUAVE
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
  100    CONTINUE
*  |  |
*  |  +----------------------------------------------------------------*
*  |  +----------------------------------------------------------------*
*  |  |
         DO 200 IE = 1,5
            IA  = IEX (2,IE)
            IR1 = IEX (1,IE)
            IR2 = IEX (3,IE)
            ANUEN1 = 0.D+00
            ANUEN2 = 0.D+00
            AAB (IA) = 0.D+00
            AEXTR1 = 0.5D+00 * MIN ( 1, IR1 - 1 )
            AEXTR2 = 0.5D+00 * MIN ( 1, IR2 - 1 )
*  |  |  +-------------------------------------------------------------*
*  |  |  |          This loop fills the ab(i,ia) array,
*  |  |  |          extrapolating from the data for <nu>s,
*  |  |  |          anub (ir1,r2)
            DO 130 IN = 1, INMAX
               ANUBE1 = ANUEN1
               ANUBE2 = ANUEN2
               ANUEN1 = ( AEXTR1 + IN  ) * ANUB (IA)
     &                / ANUB (IR1)
               ANUEN2 = ( AEXTR2 + IN  ) * ANUB (IA)
     &                / ANUB (IR2)
               NABE1  = INT (ANUBE1)
               NABE2  = INT (ANUBE2)
               NAEN1  = NINT (ANUEN1)
               NAEN2  = NINT (ANUEN2)
               IF ( AB (IN,IR1) .GT. 0.D+00 ) THEN
                  DO 110 INN = MAX (NABE1,1), NAEN1
                     IF (INN .GT. 1) THEN
                        ANUBE0 = 0.5D+00 + (INN-1)
                        IF (INN .GT. INMAX ) GO TO 110
                     ELSE
                        ANUBE0 = 0.D+00
                     END IF
                     ANUEN0 = 0.5D+00 + INN
                     ANUBEG = MAX ( ANUBE1, ANUBE0 )
                     ANUEND = MIN ( ANUEN1, ANUEN0 )
                     WEIGHT = AB (IN,IR1) * MAX ( ZERZER, ( ANUEND
     &                      - ANUBEG ) / ( ANUEN1 - ANUBE1 ) )
                     AB (INN,IA) = AB (INN,IA) + WEIGHT
                     AAB (IA) = AAB (IA) + WEIGHT
  110             CONTINUE
               END IF
               IF ( AB (IN,IR2) .GT. 0.D+00 ) THEN
                  DO 120 INN = MAX (NABE2,1), NAEN2
                     IF (INN .GT. 1) THEN
                        ANUBE0 = 0.5D+00 + (INN-1)
                        IF (INN .GT. INMAX ) GO TO 120
                     ELSE
                        ANUBE0 = 0.D+00
                     END IF
                     ANUEN0 = 0.5D+00 + INN
                     ANUBEG = MAX ( ANUBE2, ANUBE0 )
                     ANUEND = MIN ( ANUEN2, ANUEN0 )
                     WEIGHT = AB (IN,IR2) * MAX ( ZERZER, ( ANUEND
     &                      - ANUBEG ) / ( ANUEN2 - ANUBE2 ) )
                     AB (INN,IA) = AB (INN,IA) + WEIGHT
                     AAB (IA) = AAB (IA) + WEIGHT
  120             CONTINUE
               END IF
  130       CONTINUE
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            DO 140 IN =1, INMAX
               AB (IN,IA) = AB (IN,IA) / AAB (IA)
  140       CONTINUE
*  |  |  |
*  |  |  +-------------------------------------------------------------*
  200    CONTINUE
*  |  |
*  |  +----------------------------------------------------------------*
*  |  +----------------------------------------------------------------*
*  |  |             This loop simply creates a cumulative distribution
         DO 12 IA = 1, IAMAX
            AAB (IA) = 0.D+00
            NULAST (IA) = INMAX
*  |  |  +-------------------------------------------------------------*
*  |  |  |          This loop computes the normalization factor
            DO 18 IN = INMAX, 1, -1
               AAB (IA) = AAB (IA) + AB (IN,IA)
               IF ( AB (IN,IA) .LE. 0.D+00 ) NULAST (IA) = IN - 1
   18       CONTINUE
*  |  |  |
*  |  |  +-------------------------------------------------------------*
            ANUAVE = AB (1,IA)
            ANSESQ (IA) = 0.D+00
            AB   (1,IA) = AB (1,IA) / AAB (IA)
*  |  |  +-------------------------------------------------------------*
*  |  |  |          Create the cumulative distribution
            DO 13 IN = 2, INMAX
               AB (IN,IA) = AB (IN,IA) / AAB (IA) + AB (IN-1,IA)
               ANUAVE = ANUAVE + IN * ( AB (IN,IA)
     &                - AB (IN-1,IA) )
               ANSESQ (IA) = ANSESQ (IA) + (IN-1)*(IN-1)
     &                     * ( AB (IN,IA) - AB (IN-1,IA) )
   13       CONTINUE
*  |  |  |
*  |  |  +-------------------------------------------------------------*
            ANUB   (IA) = ANUAVE
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF ( IPOWER .LT. 7 ) THEN
               ANPOWR (IA) = AB (1,IA)
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            ELSE IF ( IPOWER .LT. 11 ) THEN
               ANPOWR (IA) = AB (1,IA) * FPOWER ( IPOWER, 1, ANUAVE )
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            ELSE
               ANPOWR (IA) = AB (1,IA) * FPOWER ( IPOWER, 1, ANUAVE )
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  +-------------------------------------------------------------*
*  |  |  |          Compute < nu**y(nu) >:
            DO 11 IN = 2, INMAX
               IF ( IPOWER .LT. 7 ) THEN
                  IF ( APOWER .GT. 0.D+00 ) THEN
                     DPOWER = APOWER
                  ELSE
                     DPOWER = FPOWER ( IPOWER, IN, ANUAVE )
                  END IF
                  ANPOWR (IA) = ANPOWR (IA) + IN**DPOWER
     &                        * ( AB (IN,IA) - AB (IN-1,IA) )
               ELSE IF ( IPOWER .LT. 11 ) THEN
                  ANPOWR (IA) = ANPOWR (IA) + IN
     &                        * FPOWER ( IPOWER, IN, ANUAVE )
     &                        * ( AB (IN,IA) - AB (IN-1,IA) )
               ELSE
                  ANPOWR (IA) = ANPOWR (IA) + IN
     &                        **FPOWER ( IPOWER, IN, ANUAVE )
     &                        * ( AB (IN,IA) - AB (IN-1,IA) )
               END IF
   11       CONTINUE
*  |  |  |
*  |  |  +-------------------------------------------------------------*
   12    CONTINUE
*  |  |
*  |  +----------------------------------------------------------------*
*  |  +----------------------------------------------------------------*
*  |  |
         DO 17 IA = 1, IAMAX
            ANUAC = AB (1,IA)
            ANUSQ = AB (1,IA)
            ANUTH = AB (1,IA)
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            DO 15 I = 2, NULAST (IA)
               ANUAC = ANUAC + I     * ( AB (I,IA) - AB(I-1,IA) )
               ANUSQ = ANUSQ + I*I   * ( AB (I,IA) - AB(I-1,IA) )
               ANUTH = ANUTH + I*I*I * ( AB (I,IA) - AB(I-1,IA) )
   15       CONTINUE
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF ( IA .GT. 1 ) THEN
               ANUSE2 = ( ANUSQ - 2.D+00 * ANUB (IA) + 1.D+00 ) /
     &                ( ANUB (IA) - 1.D+00 )**2
               ANUSE3 = ( ANUTH - 3.D+00 * ANUSQ + 3.D+00 * ANUB (IA)
     &                - 1.D+00 ) / ( ANUB (IA) - 1.D+00 )**3
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            ELSE
               ANUSE2 = 1.D+00
               ANUSE3 = 1.D+00
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
   17    CONTINUE
*  |  |
*  |  +----------------------------------------------------------------*
         NUDISV = -1
         RETURN
      END IF
*  |
*  +-------------------------------------------------------------------*
*  +-------------------------------------------------------------------*
*  |                Select the proper distributions for interpolations
      DO 40 I = 1, IAMAX
         IF (ANU .LT. ANUB(I)) GO TO 41
   40 CONTINUE
*  |
*  +-------------------------------------------------------------------*
      I=IAMAX + 1
   41 CONTINUE
      NANU  = I - 1
      IF (NANU.GE.IAMAX) GO TO 50
      IF (NANU.LE.0) NANU = 1
      NANUN = NANU + 1
      WEIGH1 = ( ANU - ANUB (NANU)  ) / ( ANUB (NANUN) - ANUB (NANU) )
      WEIGH2 = ( ANUB (NANUN) - ANU ) / ( ANUB (NANUN) - ANUB (NANU) )
      DIST (0) = 0.D+00
      ASEASQ = WEIGH1 * ANSESQ (NANUN) + WEIGH2 * ANSESQ (NANU)
      APOWER = WEIGH1 * ANPOWR (NANUN) + WEIGH2 * ANPOWR (NANU)
*  +-------------------------------------------------------------------*
*  |   Create the proper cumulative distribution by interpolation
*  |   ( note that since weigh1 + weigh2 = 1 it will be automatically
*  |   normalized )
      DO 20 IN = 1, NULAST (NANUN)
          DIST (IN) = WEIGH1 * AB (IN,NANUN) + WEIGH2 * AB (IN,NANU)
          WEIGHT = (IN-1) * ( DIST (IN) - DIST (IN-1) )
   20 CONTINUE
*  |
*  +-------------------------------------------------------------------*
      PRZERO = DIST (1)
      CALL GRNDM(RNDM,1)
      X=RNDM(1)
*  +-------------------------------------------------------------------*
*  |   Compute the proper nu from the cumulative distribution
      DO 30 IN = 1, NULAST (NANUN)
         IF ( X .LE. DIST (IN) ) GO TO 31
   30 CONTINUE
*  |
*  +-------------------------------------------------------------------*
      IN = NULAST (NANUN)
   31 CONTINUE
      NUDISV = IN
      RETURN
*  Come here for <nu> larger than 8
   50 CONTINUE
      APOWER = -1.D+00
      ASEASQ = -1.D+00
      NUD    = KPOIS (EXPLAM)
      NUDISV = NUD + 1
      PRZERO = EXPLAM
      WRITE (LUNERR,*)' *** Nudisv called with <nu> >= 8 !!',REAL(ANU),
     &                ' ***'
      RETURN
*=== end of function Nudisv ===========================================*
      END
