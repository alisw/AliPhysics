#ifndef FSOURCM_H
#define FSOURCM_H 1

#include "Rtypes.h"
#include "cfortran.h"
extern "C" {
//*$ create sourcm.add
//*copy sourcm
//*
//*=== Sourcm ===========================================================*
//*
//*----------------------------------------------------------------------*
//*                                                                      *
//*     Include file: sourcm                                             *
//*                                                                      *
//*     version       march 1996     by   Alfredo Ferrari, INFN - Milan  *
//*                                                                      *
//*     Last change on 09-mar-02     by    Alfredo Ferrari               *
//*                                                                      *
//*     Included in the following subroutines or functions:              *
//*                                                                      *
//*            BDNOPT                                                    *
//*            FLKEND                                                    *
//*            EVTDAT                                                    *
//*            FEEDER                                                    *
//*            FLUKAM                                                    *
//*            MGDRAW                                                    *
//*            SOURCE                                                    *
//*                                                                      *
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
} sourcmCommon;
#define SOURCM COMMON_BLOCK(SOURCM,sourcm)
COMMON_BLOCK_DEF(sourcmCommon,SOURCM);

typedef struct {
   Char_t   sdusou[8];
} chsocmCommon;
#define CHSOCM COMMON_BLOCK(CHSOCM,chsocm)
COMMON_BLOCK_DEF(chsocmCommon,CHSOCM);

typedef struct {
   Double_t xn[3];
   Double_t wn[3];
   Double_t distn;
   Int_t    irltno;
   Int_t    irrnor;
} norlatCommon;
#define NORLAT COMMON_BLOCK(NORLAT,norlat)
COMMON_BLOCK_DEF(norlatCommon,NORLAT);   

typedef struct {
   Double_t xb[3];
   Double_t wb[3];
   Double_t wp[3];
   Double_t xp[3];
   Double_t rin;
   Double_t rout;
   Double_t pinf;
   Double_t dist;
   Int_t    ir;
   Int_t    idbg;
   Int_t    irprim;
   Int_t    nasc;
   Int_t    lsurf;
   Int_t    nbo;
   Int_t    lri;
   Int_t    lro;
   Int_t    kloop;
   Int_t    loop;
   Int_t    itype;
   Int_t    n0a;
} paremCommon;
#define PAREM COMMON_BLOCK(PAREM,parem)
COMMON_BLOCK_DEF(paremCommon,PAREM);   
}

#endif
