*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:04  cernlib
* Geant
*
*
#include "geant321/pilot.h"
#if defined(CERNLIB_OLDNAME)
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*=== zeroin ===========================================================*
*
      SUBROUTINE ZEROIN
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*======================================================================*
*                                                                      *
*     Created by Alfredo Ferrari on 11 july 1990   INFN - Milan        *
*                                                                      *
*     Last change  on    23-apr-93  by  Alfredo Ferrari                *
*                                                                      *
*     This routine zeroes a few arrays, and initialize some data       *
*                                                                      *
*======================================================================*
*
#include "geant321/paprop.inc"
#include "geant321/part2.inc"
*
*  +-------------------------------------------------------------------*
*  |               Set the "effective" masses for particles
      DO 8000 I = 1, IDMAX8
*  |  +----------------------------------------------------------------*
*  |  |            Baryons: "effective" mass = actual one - baryon
*  |  |                      number x proton mass (final state of any
*  |  |                      baryon)
         IF ( IIBAR (I) .GT. 0 ) THEN
            AAMDSC (I) = AAM (I) - IIBAR (I) * AAM (1)
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |            Antibaryons: "effective mass" = actual mass +
*  |  |                          |baryon number| x proton mass
         ELSE IF ( IIBAR (I) .LT. 0 ) THEN
            AAMDSC (I) = AAM (I) - IIBAR (I) * AAM (1)
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |            Mesons: "effective" mass = actual mass
         ELSE
            AAMDSC (I) = AAM (I)
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
 8000 CONTINUE
*  |
*  +-------------------------------------------------------------------*
*    Electron: "effective" mass = 0
      AAMDSC (3) = ZERZER
*    Positron: "effective" mass = 2 x mass
      AAMDSC (4) = 2.D+00 * AAM (4)
*    Neutron: "effective" mass = 0 (it should be taken into account the
*                         neutrino energy when decaying)
      AAMDSC (8) = ZERZER
*    Muon+: "effective" mass = muon mass + positron mass
      AAMDSC (10) = AAM (10) + AAM (4)
*    Muon-: "effective" mass = muon mass - electron mass
*           (final stable state)
      AAMDSC (11) = AAM (11) - AAM (3)
*    Tau+: "effective" mass = Tau mass + positron mass
      AAMDSC (131) = AAM (131) + AAM (4)
*    Tau-: "effective" mass = Tau mass - electron mass
*           (final stable state)
      AAMDSC (132) = AAM (132) - AAM (3)
*  +-------------------------------------------------------------------*
*  |               Set the "effective" masses for particles in Paprop
      DO 8500 I = 1, NALLWP
         KP = IPTOKP (I)
*  |  +----------------------------------------------------------------*
*  |  |
         IF ( KP .GT. 0 ) THEN
            AMDISC (I) = AAMDSC (KP)
            AM     (I) = AAM  (KP)
            ICHRGE (I) = IICH (KP)
            THALF  (I) = TAU  (KP)
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE IF ( I .NE. 30 ) THEN
            AMDISC (I) = ZERZER
            AM     (I) = ZERZER
            ICHRGE (I) = 0
            THALF  (I) = ZERZER
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE
            AMDISC (I) = AAMDSC (1)
            AM     (I) = AAM  (1)
            ICHRGE (I) = IICH (1)
            THALF  (I) = TAU  (1)
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
*  |  +----------------------------------------------------------------*
*  |  |   Set the particle spin (in 1/2 units)
*  |  |   If the particle is a hadron :
*  |  |            mesons         --> 0 spin
*  |  |       baryons,antibaryons --> 1/2 spin
*  |  |   otherwise :
*  |  |       e-/+,mu-/+,vu,vubar --> 1/2 spin
*  |  |             photons       --> 1 spin
         IF ( LHADRO (I) .AND. IIBAR (KP) .EQ. 0 ) THEN
            JSPINP (I) = 0
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE
            IF ( I .EQ. 7 ) THEN
               JSPINP (I) = 2
            ELSE
               JSPINP (I) = 1
            END IF
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
 8500 CONTINUE
*  |
*  +-------------------------------------------------------------------*
*=== end of subroutine zeroin =========================================*
      RETURN
      END
#endif
