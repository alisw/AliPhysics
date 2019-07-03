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
  kBSemiElectronic, kSemiElectronic, kDiElectron, kBSemiMuonic, kDSemiMuonic, kSemiMuonic, kDiMuon, kJpsiDiMuon,
    kBJpsiDiMuon, kBJpsiDiElectron, 
    kBPsiPrimeDiMuon, kBPsiPrimeDiElectron, kPiToMu, kKaToMu, 
    kNoDecay, kHadronicD, kHadronicDWithout4Bodies, kOmega, kLambda, kPhiKK, 
    kAll, kNoDecayHeavy, kHardMuons, kBJpsi,  kBJpsiUndecayed,
    kWToMuon,kWToCharm, kWToCharmToMuon, kZDiMuon, kZDiElectron, kNeutralPion, kAllMuonic,
    kChiToJpsiGammaToMuonMuon, kChiToJpsiGammaToElectronElectron, kNoDecayBeauty, kPsiPrimeJpsiDiElectron,
    kElectronEM, kGammaEM, kDiElectronEM, kBeautyUpgrade,
    kHadronicDWithV0,kHadronicDWithout4BodiesWithV0,kAllEM,
    kLcpKpi, kLcpK0S, kHFYellowReport, kHadronicDPionicD0, kHadronicDWithV0PionicD0, kHadronicDWithout4BodiesPionicD0,
    kHadronicDWithout4BodiesWithV0PionicD0,kHadronicDPionicD0pure,kHadronicDPionicD0K,kHadronicDPionicD0pi, kLcpK0SBDTsig, kEtaPrime
} Decay_t;
#endif

#endif //ALI_DECAYER__H

