#include "Fdimpar.h"
extern "C" {
/*
*$ CREATE STEPSZ.ADD
*COPY STEPSZ
*
*=== stepsz ===========================================================*
*
*----------------------------------------------------------------------*
*                                                                      *
*   Common stepsz for setting the minimum and maximum step sizes on a  *
*                 a region by region basis: very useful for vacuum re- *
*                 gions with magnetic filed and for saving time ( and  *
*                 accuracy ) with the new plc and lca algorithm in     *
*                 Emf and Fluka                                        *
*                                                                      *
*          W A R N I N G !!!!! At the moment implemented only for      *
*          electron and positron transport in Emf and for charged      *
*          particles transport in Fluka with the new multiple scat-    *
*          tering module!!!!!!                                         *
*                                                                      *
*                  created by A. Ferrari & P. Sala on 14-jan-1990      *
*                                                                      *
*          included in:                                                *
*                        fiprou                                        *
*                        flukam (main)                                 *
*                        kashea                                        *
*                        kaskad                                        *
*                        electr (new version)                          *
*                        mageas                                        *
*                        magnew                                        *
*                        zeroin                                        *
*                                                                      *
*                        Stepmn  = minimum step size (cm)              *
*                        Stepmx  = maximum step size (cm)              *
*                        Mxxrgn = maximum number of regions            *
*                                                                      *
*----------------------------------------------------------------------*
*/

//      COMMON / STEPSZ / STEPMN ( MXXRGN ), STEPMX ( MXXRGN )
typedef struct {
    Double_t stepmn[mxxrgn];
    Double_t stepmx[mxxrgn];
} stepszCommon;

#define STEPSZ COMMON_BLOCK(STEPSZ,stepsz)
COMMON_BLOCK_DEF(stepszCommon,STEPSZ);
 
}

