#ifndef ROOT_UA1Common
#define ROOT_UA1Common

#ifndef __CFORTRAN_LOADED
#include "cfortran.h"
#endif

extern "C" {
// COMMON /UA1JETS/ NJET, ETJ(100), ETAJ(100,2), PHIJ(100,2), NCELLJ(100)

    typedef struct {
	Int_t   njet;
	Float_t etj[100];
	Float_t etaj[2][100];
	Float_t phij[2][100];
	Int_t   ncellj[100];
    } UA1JetsCommon;

#define UA1JETS COMMON_BLOCK(UA1JETS,ua1jets)
COMMON_BLOCK_DEF(UA1JetsCommon,UA1JETS);

// COMMON /UA1CELL/ etaCellSize, phiCellSize

    typedef struct {
	Float_t etaCellSize;
	Float_t phiCellSize;
    } UA1CellCommon;

#define UA1CELL COMMON_BLOCK(UA1CELL,ua1cell)
COMMON_BLOCK_DEF(UA1CellCommon,UA1CELL);

// COMMON /UA1PARA/ cone_rad, et_seed, ej_min, et_min
    typedef struct {
	Float_t coneRad;
	Float_t etSeed;
	Float_t ejMin;
	Float_t etMin;
    } UA1ParaCommon;

#define UA1PARA COMMON_BLOCK(UA1PARA,ua1para)
COMMON_BLOCK_DEF(UA1ParaCommon,UA1PARA);
    
}
#endif



