#ifndef ROOT_ECommon
#define ROOT_ECommon

#ifndef __CFORTRAN_LOADED
#include "cfortran.h"
#endif

extern "C" {
// COMMON /EMCALJETS/ NJET, ETJ(100), ETAJ(100,2), PHIJ(100,2), NCELLJ(100)

    typedef struct {
	Int_t   njet;
	Float_t etj[100];
	Float_t etaj[2][100];
	Float_t phij[2][100];
	Int_t   ncellj[100];
    } EmcalJetsCommon;

#define EMCALJETS COMMON_BLOCK(EMCALJETS,emcaljets)
COMMON_BLOCK_DEF(EmcalJetsCommon,EMCALJETS);

// COMMON /EMCALGEO/ etaCellSize, phiCellSize

    typedef struct {
	Float_t etaCellSize;
	Float_t phiCellSize;
    } EmcalCellGeoCommon;

#define EMCALCELLGEO COMMON_BLOCK(EMCALCELLGEO,emcalcellgeo)
COMMON_BLOCK_DEF(EmcalCellGeoCommon,EMCALCELLGEO);

// COMMON /EMCALJETPARAM/ cone_rad, et_seed, ej_min, et_min
    typedef struct {
	Float_t coneRad;
	Float_t etSeed;
	Float_t ejMin;
	Float_t etMin;
    } EmcalJetParamCommon;

#define EMCALJETPARAM COMMON_BLOCK(EMCALJETPARAM,emcaljetparam)
COMMON_BLOCK_DEF(EmcalJetParamCommon,EMCALJETPARAM);
    
}
#endif



