#ifndef FLTCLCM_H
#define FLTCLCM_H 1

#include "cfortran.h"
#include "Rtypes.h"

extern "C" {
//*$ create ltclcm.add
//*copy ltclcm
//*
//*=== ltclcm ===========================================================*
//*
//*----------------------------------------------------------------------*
//*                                                                      *
//*     lattice cell common:                                             *
//*                                                                      *
//*     created on 09 december 1993  by    alfredo ferrari & paola sala  *
//*                                                   infn - milan       *
//*                                                                      *
//*     last change on 10-dec-93     by    alfredo ferrari               *
//*                                                                      *
//*                                                                      *
//*----------------------------------------------------------------------*
//*

typedef struct {
   Int_t    mlattc;
   Int_t    newlat;
   Int_t    mlatld;
   Int_t    mlatm1;
   Int_t    mltsen;
   Int_t    mltsm1;
} ltclcmCommon;
#define LTCLCM COMMON_BLOCK(LTCLCM,ltclcm)
COMMON_BLOCK_DEF(ltclcmCommon,LTCLCM);
}

#endif
