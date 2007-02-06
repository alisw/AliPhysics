#ifndef ALI_DECAYER__H
#define ALI_DECAYER__H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "RVersion.h"
#include "TVirtualMCDecayer.h"

typedef TVirtualMCDecayer AliDecayer;

#if ROOT_VERSION_CODE >= 197633  //Corresponds to Root v3-04-01
typedef enum 
{
    kBSemiElectronic, kSemiElectronic, kDiElectron, kBSemiMuonic, kSemiMuonic, kDiMuon,
    kBJpsiDiMuon, kBJpsiDiElectron, 
    kBPsiPrimeDiMuon, kBPsiPrimeDiElectron, kPiToMu, kKaToMu, 
    kNoDecay, kHadronicD, kOmega, kPhiKK, 
    kAll, kNoDecayHeavy, kHardMuons, kBJpsi,
    kWToMuon,kWToCharm, kWToCharmToMuon, kZDiMuon
} Decay_t;
#endif

#endif //ALI_DECAYER__H

