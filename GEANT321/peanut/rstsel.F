*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:22:03  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.46  by  S.Giani
*-- Author :
*$ CREATE RSTSEL.FOR
*COPY RSTSEL
*
*=== rstsel ===========================================================*
*
      SUBROUTINE RSTSEL ( JPROJ )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*
#include "geant321/balanc.inc"
#include "geant321/nucdat.inc"
#include "geant321/nucgeo.inc"
#include "geant321/part.inc"
#include "geant321/parevt.inc"
#include "geant321/resnuc.inc"
*
      REAL RNDM(1)
      SAVE ABTAR , ZZTAR , PACORE, PASKIN,
     &     PAHALO, PACHLP, XHLP  , AUSFL , ZUSFL , ONUSFL, KPROJ ,
     &     ITFRMI, ITFRM2
*
      KPROJ  = JPROJ
      ABTAR  = IBTAR
      ZZTAR  = ICHTAR
      AUSFL  = ABTAR
      ZUSFL  = ZZTAR
      ONUSFL = 0.D+00
      CALL RACO ( CXIMPC, CYIMPC, CZIMPC )
      UBIMPC = CXIMPC
      VBIMPC = CYIMPC
      WBIMPC = CZIMPC
      PACORE = FRHINC (RADIU0) / ABTAR
      PASKIN = FRHINC (RADIU1) / ABTAR - PACORE
      PAHALO = 1.D+00 - PACORE - PASKIN
      GO TO 500
      ENTRY RSTNXT
  500 CONTINUE
      LABRST = .TRUE.
      IF ( KPROJ .EQ. 2 ) THEN
         NTARGT = 1
         STOP 'pbar_rest'
      ELSE IF ( KPROJ .EQ. 9 ) THEN
         NTARGT = 1
         STOP 'nbar_rest'
      ELSE IF ( KPROJ .EQ. 14 ) THEN
         CALL GRNDM(RNDM,1)
         IF ( RNDM (1) .LT. RADPIM * ( ABTAR + RDPMHL ) / ABTAR )
     &      THEN
            NTARGT = 1
            KNUCIM = 1
            ITFRMI = 1
            IPRTYP = - ( 14 * 10 + KNUCIM )
         ELSE
            NTARGT = 2
            KNUCIM = 1
            ITFRMI = 1
            PRCOL0 = APMRST
            PRCOLP = PRCOL0 * ( ZUSFL - ONUSFL )
            PRCOLP = PRCOLP / ( PRCOLP + ( 1.D+00 - PRCOL0 ) * ( AUSFL
     &             - ZUSFL ) )
            CALL GRNDM(RNDM,1)
            IF ( RNDM (1) .LT. PRCOLP ) THEN
               KNUCI2 = 1
               ITFRM2 = 1
            ELSE
               KNUCI2 = 8
               ITFRM2 = 2
            END IF
            IPRTYP = - ( 14 * 100 + KNUCIM * 10 + KNUCI2 )
         END IF
      ELSE IF ( KPROJ .EQ. 15 ) THEN
         NTARGT = 1
         NTARGT = 2
         STOP 'K-_rest'
      ELSE
         STOP '???_rest'
      END IF
      CALL GRNDM(RNDM,1)
      RNDMAB = RNDM (1)
      IF ( RNDMAB .LT. PACORE ) THEN
         RNDMAB = RNDMAB / PACORE
         RADSAM = RADIU0
         RHOSAM = RHOCEN
         XHLP   = 0.D+00
      ELSE IF ( RNDMAB - PACORE .LT. PASKIN ) THEN
         RNDMAB = ( RNDMAB - PACORE ) / PASKIN
         RADSAM = RADIU1
         RHOSAM = RHOCOR
         XHLP   = RADIU0 / RADIU1
      ELSE
         RNDMAB = ( RNDMAB - PACORE - PASKIN ) / PAHALO
         RHOSAM = RHOSKN
         RADSAM = RADTOT
         XHLP   = RADIU1 / RADTOT
      END IF
      XHLP3  = XHLP * XHLP * XHLP
 1000 CONTINUE
         BIMPCT = RADSAM * ( RNDMAB * ( 1.D+00 - XHLP3 ) + XHLP3 )
     &          **0.3333333333333333D+00
         RHOIMP = FRHONC (BIMPCT)
         CALL GRNDM(RNDM,1)
         RNDMAB = RNDM (1)
         RHORAT = RHOIMP / RHOSAM
         IF ( RNDMAB .GE. RHORAT ) THEN
            RNDMAB = ( RNDMAB - RHORAT ) / ( 1.D+00 - RHORAT )
            GO TO 1000
         END IF
      PFRIMP = FPFRNC ( RHOIMP, ITFRMI )
      EKFIMP = FEKFNC ( PFRIMP, ITFRMI )
      RHOIMT = RHOIMP
      IF ( NTARGT .EQ. 2 ) THEN
         IF ( RHOIMP .GE. RHOCEN ) THEN
            PFRIM2 = PFRCEN (ITFRM2)
            EKFIM2 = EKFCEN (ITFRM2)
            EKFPRO = EKFCEN (IPWELL)
            PFRPRO = PFRCEN (IPWELL)
         ELSE IF ( ITFRM2 .NE. ITFRMI ) THEN
            PFRIM2 = PFRIMP * PFRCEN (ITFRM2) / PFRCEN (ITFRMI)
            EKFIM2 = SQRT ( PFRIM2 * PFRIM2 + AMNUSQ (ITFRM2) )
     &             - AMNUCL (ITFRM2)
            IF ( IPWELL .EQ. ITFRMI ) THEN
               EKFPRO = EKFIMP
               PFRPRO = PFRIMP
            ELSE
               PFRPRO = PFRIM2
               EKFPRO = EKFIM2
            END IF
         ELSE
            PFRIM2 = PFRIMP
            EKFIM2 = EKFIMP
            IF ( IPWELL .EQ. ITFRMI ) THEN
               PFRPRO = PFRIMP
               EKFPRO = EKFIMP
            ELSE
               PFRPRO = PFRIMP * PFRCEN (IPWELL) / PFRCEN (ITFRMI)
               EKFPRO = SQRT ( PFRPRO * PFRPRO + AMNUSQ (IPWELL) )
     &                - AMNUCL (IPWELL)
            END IF
         END IF
      ELSE
         IF ( RHOIMP .GE. RHOCEN ) THEN
            PFRPRO = PFRCEN (IPWELL)
            EKFPRO = EKFCEN (IPWELL)
         ELSE IF ( IPWELL .EQ. ITFRMI ) THEN
            PFRPRO = PFRIMP
            EKFPRO = EKFIMP
         ELSE
            PFRPRO = PFRIMP * PFRCEN (IPWELL) / PFRCEN (ITFRMI)
            EKFPRO = SQRT ( PFRPRO * PFRPRO + AMNUSQ (IPWELL) )
     &             - AMNUCL (IPWELL)
         END IF
      END IF
      VPRWLL = WLLRED * ( EKFPRO + BNDNUC )
      BIMPTR = BIMPCT
      BIMMEM = BIMPTR
      RIMPCT = BIMPCT
      RIMPTR = BIMPTR
      XBIMPC = UBIMPC * BIMPCT
      YBIMPC = VBIMPC * BIMPCT
      ZBIMPC = WBIMPC * BIMPCT
      XIMPCT = XBIMPC
      YIMPCT = YBIMPC
      ZIMPCT = ZBIMPC
      XIMPTR = XIMPCT
      YIMPTR = YIMPCT
      ZIMPTR = ZIMPCT
      EKEWLL = EKECON + VPRWLL
      PPRWLL = SQRT ( EKEWLL * ( EKEWLL + 2.D+00 * AM (KPROJ) ) )
      PXPROJ = PPRWLL * CXIMPC
      PYPROJ = PPRWLL * CYIMPC
      PZPROJ = PPRWLL * CZIMPC
      PNFRMI = PFNCLV ( ITFRMI, .TRUE. )
      IF ( PNFRMI .LT. -100.D+00 ) GO TO 500
      CALL RACO ( PXFERM, PYFERM, PZFERM )
      PXFERM = PXFERM * PNFRMI
      PYFERM = PYFERM * PNFRMI
      PZFERM = PZFERM * PNFRMI
      EKFIMP = MAX ( EKFERM, EKFIMP )
      IF ( NTARGT .EQ. 2 ) THEN
         PNFRM2 = PFNCLV ( ITFRM2, .FALSE. )
         IF ( PNFRM2 .LT. -100.D+00 ) GO TO 500
         CALL RACO ( PXFER2, PYFER2, PZFER2 )
         PXFER2 = PXFER2 * PNFRM2
         PYFER2 = PYFER2 * PNFRM2
         PZFER2 = PZFER2 * PNFRM2
         EKFIM2 = MAX ( EKFER2, EKFIM2 )
         ERES   = EKEWLL + AM (KPROJ) + AM (KNUCIM) + EKFERM
     &          + AM (KNUCI2) + EKFER2
         PXRES  = PXPROJ + PXFERM + PXFER2
         PYRES  = PYPROJ + PYFERM + PYFER2
         PZRES  = PZPROJ + PZFERM + PZFER2
         PTRES2 = PXRES**2 + PYRES**2 + PZRES**2
      ELSE
         ERES   = EKEWLL + AM (KPROJ) + AM (KNUCIM) + EKFERM
         PXRES  = PXPROJ + PXFERM
         PYRES  = PYPROJ + PYFERM
         PZRES  = PZPROJ + PZFERM
         PTRES2 = PXRES**2 + PYRES**2 + PZRES**2
      END IF
      RETURN
*=== End of subroutine rstsel =========================================*
      END
