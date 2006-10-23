#include "cfortran.h"
#include "Rtypes.h"

#include "Fdimpar.h"

extern "C" {
/*$ CREATE DPDXCM.ADD
*COPY DPDXCM
*
*=== dpdxcm ===========================================================*
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 1989-2006         by        Alfredo Ferrari        *
*     All Rights Reserved.                                             *
*                                                                      *
*                                                                      *
*     Include file: dpdxcm  (DP/DX CoMmon)                             *
*                                                                      *
*     Created  on  10 february 1991   by        Alfredo Ferrari        *
*                                                INFN - Milan          *
*                                                                      *
*     Last change on  21-may-06       by        Alfredo Ferrari        *
*                                                                      *
*     Included in the following routines:                              *
*                                                                      *
*              blockmvax/bdtrns.f                                      *
*              dedxmvax/dedsdv.f                                       *
*              dedxmvax/dedx.f                                         *
*              dedxmvax/dedxfl.f                                       *
*              dedxmvax/deltad.f                                       *
*              dedxmvax/deltar.f                                       *
*              dedxmvax/deltas.f                                       *
*              dedxmvax/delthr.f                                       *
*              dedxmvax/delths.f                                       *
*              dedxmvax/dpdx.f                                         *
*              dedxmvax/dpdxio.f                                       *
*              dedxmvax/enion.f                                        *
*              dedxmvax/enionf.f                                       *
*              dedxmvax/gdedxc.f                                       *
*              dedxmvax/heabre.f                                       *
*              dedxmvax/hosufl.f                                       *
*              dedxmvax/hvbrem.f                                       *
*              dedxmvax/hvpair.f                                       *
*              dedxmvax/t0zffc.f                                       *
*              elsmvax/sigtab.f                                        *
*              emfmvax/ededxf.f                                        *
*              emfmvax/emfin.f                                         *
*              emfmvax/emfret.f                                        *
*              emfmvax/emfsco.f                                        *
*              emfmvax/emfstp.f                                        *
*              emfmvax/pdedxf.f                                        *
*              kaskadmvax/hmsnsc.f                                     *
*              kaskadmvax/kashea.f                                     *
*              kaskadmvax/kaskad.f                                     *
*              kaskadmvax/mulhad.f                                     *
*              kaskadmvax/omegah.f                                     *
*              kaskadmvax/sgttot.f                                     *
*              kaskadmvax/stepop.f                                     *
*              mainmvax/deflts.f                                       *
*              mainmvax/dltcrd.f                                       *
*              mainmvax/matcrd.f                                       *
*              mainmvax/zeroin.f                                       *
*                                                                      *
*           Avionp (m) = average ionization potential (eV) of medium m *
*           Ccster (m) = Sternheimer cbar   parameter for medium m     *
*           X0ster (m) = Sternheimer x0     parameter for medium m     *
*           Xester (m) = Sternheimer x1     parameter for medium m     *
*           Amster (m) = Sternheimer m      parameter for medium m     *
*           Aaster (m) = Sternheimer a      parameter for medium m     *
*           D0ster (m) = Sternheimer delta0 parameter for medium m     *
*           Aviont (m) = auxiliary ionization potential of medium m    *
*           T0dpdx (m) = delta ray production threshold of medium m    *
*                        (all particle but e+/e-)                      *
*           Tedpdx (m) = delta ray production threshold of medium m    *
*                        (electrons and positrons)                     *
*           Gaspfl (m) = pressure (atm) if a gas                       *
*           Pthrmx     = maximum momentum of the tabulations           *
*           Anpicm (m) = average number of primary ionization per cm   *
*                        for a mip for medium m (at NTP for a gas)     *
*           Frstip (m) = first ionization potential for medium m (GeV) *
*           Faltmt (m) = density modifying factor for a possible alt-  *
*                        ernate material for medium m                  *
*           Maltmt (m) = alternate material for medium m               *
*           Msdpdx (m) = possible "special material" flag for medium m *
*                        0: no special treatment                       *
*                        1: implicit delta production down to Avionr   *
*                           inside ..dedxf.. routines with recording   *
*                           of the selected values activated           *
*----------------------------------------------------------------------*
*
      PARAMETER ( MNDPDX = 50 )
      PARAMETER ( RMDPDX = 1.15D+00 )
      PARAMETER ( DPDXR1 = 0.15D+00 )
      PARAMETER ( DPDXR2 = 0.70D+00 )
      PARAMETER ( ERDEDX = 0.15D+00 * 0.15D+00 )
      PARAMETER ( MDPDXH =  4 )
*  Toln10 = 2 x log (10)
      PARAMETER ( TOLN10 = 4.605170185988091 D+00 )
*
      LOGICAL LDELTA, LPDETB, LETFUN
      COMMON / DPDXCM / P0DPDX (MPDPDX,MXXMDF), P1DPDX (MPDPDX,MXXMDF),
     &                  TMDPDX (MXXMDF), T0DPDX (MXXMDF),
     &                  TEDPDX (MXXMDF), D0DPDX (MXXMDF),
     &                  AVIONP (MXXMDF), RHORFL (MXXMDF),
     &                  GASPFL (MXXMDF), CCSTER (MXXMDF),
     &                  AMSTER (MXXMDF), XOSTER (MXXMDF),
     &                  XESTER (MXXMDF), AASTER (MXXMDF),
     &                  D0STER (MXXMDF), AVIONT (MXXMDF),
     &                  ETDPDX (MXXMDF), ALMASS (MPDPDX), PTHRMX,
     &                  FRSTIP (MXXMDF), ANPICM (MXXMDF),
     &                  FALTMT (MXXMDF), MALTMT (MXXMDF),
     &                  MSDPDX (MXXMDF), NBDPDX (MXXMDF),
     &                  KDPDXT (MPDPDX,MXXMDF),
     &                  LDELTA (MXXMDF), LPDETB (MXXMDF),
     &                  IJDPDX (-6:NALLWP), LETFUN
      SAVE / DPDXCM /
*/

    const Int_t     mndpdx = 50;
    const Double_t  rmdpdx = 1.15e0;
    const Double_t  dpdxr1 = 0.15e0;
    const Double_t  dpdxr2 = 0.70e0;
    const Double_t  erdedx = 0.15e0 * 0.15e0;
    const Int_t     mdpdxh = 4;

    typedef struct {
	Double_t p0dpdx [mxxmdf][mpdpdx];
	Double_t p1dpdx [mxxmdf][mpdpdx];
	Double_t tmdpdx [mxxmdf];
	Double_t t0dpdx [mxxmdf];
	Double_t tedpdx [mxxmdf];
	Double_t d0dpdx [mxxmdf];
	Double_t avionp [mxxmdf];
	Double_t rhorfl [mxxmdf];
	Double_t gaspfl [mxxmdf];
	Double_t ccster [mxxmdf];
	Double_t amster [mxxmdf];
	Double_t xoster [mxxmdf];
	Double_t xester [mxxmdf];
	Double_t aaster [mxxmdf];
	Double_t d0ster [mxxmdf];
	Double_t aviont [mxxmdf];
	Double_t etdpdx [mxxmdf];
	Double_t almass [mpdpdx];
	Double_t pthrmx;
	Double_t frstip [mxxmdf];
	Double_t anpicm [mxxmdf];
	Double_t faltmt [mxxmdf];
	Int_t    maltmt [mxxmdf];
	Int_t    msdpdx [mxxmdf];
	Int_t    nbdpdx [mxxmdf];
	Int_t    kdpdxt [mxxmdf][mpdpdx];
	Int_t    ldelta [mxxmdf];
	Int_t    lpdetb [mxxmdf];
	Int_t    ijdpdx [nallwp + 7];
	Int_t    letfun;
    } dpdxcmCommon;
#define DPDXCM COMMON_BLOCK(DPDXCM,dpdxcm)
COMMON_BLOCK_DEF(dpdxcmCommon, DPDXCM);
}
