*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:46  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/03 07/10/94  18.41.50  by  S.Giani
*-- Author :
      SUBROUTINE GUHADR
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       User routine to generate one hadronic interaction        *
C.    *                                                                *
C.    *    ==>Called by : GTHADR,GTNEUT                                *
C.    *                                                                *
C.    ******************************************************************
C.
*======================================================================*
*                                                                      *
*   An interface with a part of the FLUKA shower code is available in  *
*   GEANT 3.21. The following conditions and warnings must be taken    *
*   into account when using the FLUKA routines.                        *
*                                                                      *
*======================================================================*
*======================================================================*
*                                                                      *
*                   FFFF   L      U  U  K  K   AA                      *
*                   F      L      U  U  K K   A  A                     *
*                   FFF    L      U  U  KK    AAAA                     *
*                   F      L      U  U  K K   A  A                     *
*                   F      LLLL    UU   K  K  A  A                     *
*                                                                      *
*                   (C) Copyright of the authors                       *
*                                                                      *
*             A. Fasso'*, A. Ferrari#, J. Ranft$, P.R. Sala#           *
*                                                                      *
*               *: CERN, #: INFN -Milan, $: CERN/Frascati              *
*                                                                      *
*                (e-mail: FERRARIA@CERNVM.CERN.CH)                     *
*                                                                      *
*======================================================================*
*                                                                      *
*   - All the rights concerning FLUKA or parts of it are only of the   *
*     authors and are independent from those of the GEANT code         *
*                                                                      *
*   - FLUKA [1-6] is a standalone code capable of simulating the inter-*
*     action and transport of all components of EM and hadronic cas-   *
*     cades up to several TeV. However, only cross sections and        *
*     models for hadronic elastic and inelastic interactions (end 1992 *
*     status) are included in this GEANT version.                      *
*                                                                      *
*   - The most recent FLUKA model [4,6] for nucleon and pion interac-  *
*     tions in the intermediate energy range is not fully implemented  *
*     in GEANT. Only a simplified version, limited to p and n below    *
*     250 MeV, is available in GEANT 3.21. A coarser model is used for *
*     other projectiles in this energy range. However the implemented  *
*     parts should be adequate for most detector simulations and sim-  *
*     ilar applications for which GEANT is generally used. Their accu- *
*     racy could be insufficient for some nuclear physics studies or   *
*     demanding simulations at low energies, where the more sophistic- *
*     ated models [4,6] could be required.                             *
*                                                                      *
*   - The performances of GEANT-FLUKA are therefore not representative *
*     of those of FLUKA standalone and should be referenced as such    *
*     rather than simply GEANT or FLUKA.                               *
*                                                                      *
*   - The authors reserve the right of publishing about the physical   *
*     models developed for FLUKA. Running the FLUKA routines in isol-  *
*     ation for benchmarks (or equivalent use) is not permitted,       *
*     except after consultations or in collaboration with the authors. *
*                                                                      *
*   - The FLUKA routines are supposed to be included and used in       *
*     GEANT only. Any other use must be authorized by the authors.     *
*                                                                      *
*   - References: at least reference [5] should be always quoted when  *
*                 reporting results obtained with GEANT-FLUKA          *
*                                                                      *
* [1] A. Fasso', A. Ferrari, J. Ranft, P. R. Sala, G. R. Stevenson and *
*     J. M. Zazula, "FLUKA92", presented at the workshop on "Simulat-  *
*     ing Accelerator Radiation Environment", SARE, Santa Fe, 11-15    *
*     january (1993), Proceedings in press.                            *
*                                                                      *
* [2] P.A. Aarnio, A. Fasso', A. Ferrari, J.-H. Moehring, J. Ranft,    *
*     P.R. Sala, G.R. Stevenson and J.M. Zazula, Proc. MC93 Int. Conf. *
*     on Monte-Carlo Simulation in High-Energy and Nuclear Physics,    *
*     Tallahassee, Florida, 22-26 february (1993), World Scientific,   *
*     p. 88 (1994).                                                    *
*                                                                      *
* [3] P.A. Aarnio, A. Fasso', A. Ferrari, J.-H. Moehring, J. Ranft,    *
*     P.R. Sala, G.R. Stevenson and J.M. Zazula, ibidem, p. 100 (1994) *
*                                                                      *
* [4] A. Ferrari and P.R. Sala, ibidem, p. 277 (1994).                 *
*                                                                      *
* [5] A. Fasso', A. Ferrari, J. Ranft and P.R. Sala, "FLUKA: present   *
*     status and future developments", presented at the IV Int. Conf.  *
*     on Calorimetry in High Energy Physics, La Biodola (Elba),        *
*     September 19-25 1993, Proceedings in press.                      *
*                                                                      *
* [6] A. Fasso', A. Ferrari, J. Ranft, and P.R. Sala, "FLUKA: Perf-    *
*     ormances and Applications in the Intermediate Energy Range",     *
*     presented at the "Specialists' Meeting on Shielding Aspects of   *
*     Accelerators, Targets & Irradiation Facilities", Arlington,      *
*     April 28-29 1994, Proceedings in press.                          *
*                                                                      *
*======================================================================*
C.
C.    ------------------------------------------------------------------
C.
C          GHEISHA only if IHADR<3 (default)
C
      CALL GHEISH
C
      END
