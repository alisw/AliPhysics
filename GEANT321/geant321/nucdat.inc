*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:40  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_NUCDAT_INC
#define CERNLIB_GEANT321_NUCDAT_INC
*
*
* nucdat.inc
*
*                                                                     *
*=== nucdat ==========================================================*
*                                                                     *
*---------------------------------------------------------------------*
*                                                                     *
*     Include file Nucdat                                             *
*                                                                     *
*     Created on 20 april 1990  by            Alfredo Ferrari         *
*                                               INFN Milan            *
*                                                                     *
*     Last change on 09-jun-92                                        *
*                                                                     *
*     Included in the following routines:                             *
*                                                                     *
*            COREVT                                                   *
*            CORRIN                                                   *
*            DISTNU                                                   *
*            DRES                                                     *
*            EVDEEX                                                   *
*            EVENTV                                                   *
*            EVEVAP                                                   *
*            FEREVV                                                   *
*            FERHAV                                                   *
*            INCINI                                                   *
*            NUCEVV                                                   *
*            NUCRIV                                                   *
*            RAKEKV                                                   *
*            RBKEKV                                                   *
*            EVENTQ                                                   *
*            NUCRIN                                                   *
*                                                                     *
*---------------------------------------------------------------------*
*
      PARAMETER ( AMUAMU = 0.93149432  D+00 )
*     PARAMETER ( AMPROT = 0.9382796   D+00 )
      PARAMETER ( AMPROT = 0.93827231  D+00 )
*     PARAMETER ( AMNEUT = 0.9395731   D+00 )
      PARAMETER ( AMNEUT = 0.93956563  D+00 )
      PARAMETER ( AMELEC = AMELCT )
      PARAMETER ( R0NUCL = 1.12        D+00 )
      PARAMETER ( RCCOUL = 1.7         D+00 )
      PARAMETER ( COULPR = 0.001439965 D+00 )
*     PARAMETER ( FERTHO = 15.73       D-09 )
*     PARAMETER ( EXPEBN = 7.D+00 / 3.D+00  )
*     PARAMETER ( BEXC12 = FERTHO * 65.41634134195703D+00 )
      PARAMETER ( FERTHO = 14.33       D-09 )
      PARAMETER ( EXPEBN = 2.39        D+00 )
      PARAMETER ( BEXC12 = FERTHO * 72.40715579499394D+00 )
      PARAMETER ( AMUC12 = AMUAMU - AMELEC / 2.D+00 - BEXC12 / 12.D+00 )
      PARAMETER ( AMHYDR = AMPROT + AMELEC  )
      PARAMETER ( AMHTON = AMHYDR - AMNEUT  )
      PARAMETER ( AMNTOU = AMNEUT - AMUC12  )
      PARAMETER ( AMUCSQ = AMUC12 * AMUC12 )
      PARAMETER ( EBNDAV = 0.5D+00 * (AMPROT + AMNEUT) - AMUC12 )
*   Gammin : threshold for deexcitation gammas production, set to 1 keV
*   (this means that up to 1 keV of energy unbalancing can occur
*    during an event)
      PARAMETER ( GAMMIN = 1.0D-06 )
      PARAMETER ( GAMNSQ = 2.0D+00 * GAMMIN * GAMMIN )
*   Tvepsi : "epsilon" for excitation energy, set to gammin / 100
      PARAMETER ( TVEPSI = GAMMIN / 100.D+00 )
*
      COMMON /FKNUDA/ AV0WEL,     APFRMX,     AEFRMX,     AEFRMA,
     &                RDSNUC,     V0WELL (2), PFRMMX (2), EFRMMX (2),
     &                EFRMAV (2), AMNUCL (2), AMNUSQ (2), EBNDNG (2),
     &                VEFFNU (2), ESLOPE (2), PKMNNU (2), EKMNNU (2),
     &                PKMXNU (2), EKMXNU (2), EKMNAV (2), EKINAV (2),
     &                EXMNAV (2), EKUPNU (2), EXMNNU (2), EXUPNU (2),
     &                ERCLAV (2), ESWELL (2), FINCUP (2), AMRCAV    ,
     &                AMRCSQ    , ATO1O3    , ZTO1O3    , ELBNDE (0:100)
 

#endif
