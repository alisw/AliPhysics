#ifndef ALI_DECAYER__H
#define ALI_DECAYER__H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "RVersion.h"
#include "TVirtualMCDecayer.h"

typedef TVirtualMCDecayer AliDecayer;

#if ROOT_VERSION_CODE >= ROOT_VERSION(3,04,1)
typedef enum 
{
    kSemiElectronic, kDiElectron, kSemiMuonic, kDiMuon,
    kBJpsiDiMuon, kBJpsiDiElectron, 
    kBPsiPrimeDiMuon, kBPsiPrimeDiElectron, kPiToMu,
    kKaToMu, kNoDecay,
    kHadronicD, kOmega, kPhiKK, kAll, kNoDecayHeavy
} Decay_t;
#endif

#endif //ALI_DECAYER__H
