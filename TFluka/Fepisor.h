#ifndef FEPISOR_H
#define FEPISOR_H 1

#include "Rtypes.h"
#include "cfortran.h"
extern "C" {
//*$ create episor.add
//*copy episor
//*
//*=== episor ===========================================================*
//*
//*----------------------------------------------------------------------*
//*                                                                      *
//*     Include file: episor                                             *
//*                                                                      *
//*     version       march 1996     by   Alfredo Ferrari, INFN - Milan  *
//*                                                                      *
//*     Last change on 09-mar-02     by    Alfredo Ferrari               *
//*                                                                      *
//*     Included in the following subroutines or functions:              *
//*                                                                      *
//*            BDNOPT                                                    *
//*            EPILOG                                                    *
//*            EVTDAT                                                    *
//*            FEEDER                                                    *
//*            FLUKAM                                                    *
//*            MGDRAW                                                    *
//*            SOURCE                                                    *
//*                                                                      *
//*     Description of the common block(s) and variable(s)               *
//*                                                                      *
//*      Whasou(1-18) = user variables                                   *
//*            Tkesum = total kinetic energy of the primaries of the     *
//*                     user written source                              *
//*            Lussrc = flag to inform that the user written source was  *
//*                     used                                             *
//*            Sdusou = user character variable                          *
//*            Lsouit = source is called iteratively until it is .true.  *
//*                                                                      *
//*----------------------------------------------------------------------*
//*

typedef struct {
   Double_t whasou[18];
   Double_t tkesum;
   Int_t    lussrc;
   Int_t    lsouit;
} episorCommon;
#define EPISOR COMMON_BLOCK(EPISOR,episor)
COMMON_BLOCK_DEF(episorCommon,EPISOR);

typedef struct {
   Char_t   sdusou[8];
} chepsrCommon;
#define CHEPSR COMMON_BLOCK(CHEPSR,chepsr)
COMMON_BLOCK_DEF(chepsrCommon,CHEPSR);

typedef struct {
   Double_t xn[3];
   Double_t wn[3];
   Double_t distn;
   Int_t    irltno;
   Int_t    irrnor;
} norlatCommon;
#define NORLAT COMMON_BLOCK(NORLAT,norlat)
COMMON_BLOCK_DEF(norlatCommon,NORLAT);   
}

#endif
