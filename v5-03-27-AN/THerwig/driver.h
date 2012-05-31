#include "Rtypes.h"
#include "cfortran.h"
//
// COMMON/VVJIN/QQIN
//
extern "C" {

    typedef struct {
        char QQIN[50];
} vvjinCommon;
#define VVJIN COMMON_BLOCK(VVJIN,vvjin)
COMMON_BLOCK_DEF(vvjinCommon,VVJIN);
}
