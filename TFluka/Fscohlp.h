extern "C" {
//*$ create scohlp.add
//*copy scohlp 
//*
//*=== scohlp ===========================================================*
//*
//*----------------------------------------------------------------------*
//*                                                                      *
//*     Copyright (C) 1991-2005      by    Alfredo Ferrari & Paola Sala  *
//*     All Rights Reserved.                                             *
//*                                                                      *
//*                                                                      *
//*     Created on  05  august 1991  by    Alfredo Ferrari & Paola Sala  *
//*                                                   Infn - Milan       *
//*                                                                      *
//*     Last change on 29-oct-94     by    Alfredo Ferrari               *
//*                                                                      *
//*     Energy/Star binnings/scorings (Comscw):                          *
//*          ISCRNG = 1 --> Energy density  binning                      *
//*          ISCRNG = 2 --> Star   density  binning                      *
//*          ISCRNG = 3 --> Residual nuclei scoring                      *
//*          JSCRNG = # of the binning                                   *
//*     Flux like binnings/estimators (Fluscw):                          *
//*          ISCRNG = 1 --> Boundary crossing estimator                  *
//*          ISCRNG = 2 --> Track  length     binning                    *
//*          ISCRNG = 3 --> Track  length     estimator                  *
//*          ISCRNG = 4 --> Collision density estimator                  *
//*          ISCRNG = 5 --> Yield             estimator                  *
//*          JSCRNG = # of the binning/estimator                         *
//*                                                                      *
//*----------------------------------------------------------------------*

typedef struct {
   Int_t    iscrng;
   Int_t    jscrng;
   Int_t    lsczer;
} scohlpCommon;
#define SCOHLP COMMON_BLOCK(SCOHLP,scohlp)
COMMON_BLOCK_DEF(scohlpCommon,SCOHLP);
}
