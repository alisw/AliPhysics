/* This file is property of and copyright                                 *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/// @file   AliAnaConvCorrPhoton.h
/// @author Svein Lindal
/// @brief  Class used to find correlations between photons and charged tracks
 
#ifndef ALIANACONVCORRPHOTON_CXX
#define ALIANACONVCORRPHOTON_CXX

#include "AliAnaConvCorrBase.h"

class AliAODConversionPhoton;
class TClonesArray;

class AliAnaConvCorrPhoton : public AliAnaConvCorrBase {

public:

  AliAnaConvCorrPhoton(); 
  AliAnaConvCorrPhoton(TString name); 
  virtual ~AliAnaConvCorrPhoton();

  ///Correlation photon with tracks
  virtual void CorrelateWithHadrons(const AliAODConversionPhoton * const photon, const TClonesArray * const tracks, const Bool_t isolated, const Bool_t decayParticle);

  //Process particles identified as pion / eta decay 
  void SkipDecayParticles() { fSkipDecayParticles = kTRUE; }
  void DoDecayParticles() { fSkipDecayParticles = kFALSE; }
  void DoDecayOnly() { fSkipDecayParticles = kFALSE; fDecayOnly = kTRUE; }

 private:

  
  Bool_t fSkipDecayParticles; //Process particles identified as pion / eta decay particles
  Bool_t fDecayOnly;

  AliAnaConvCorrPhoton(const AliAnaConvCorrPhoton&); // not implemented
  AliAnaConvCorrPhoton& operator=(const AliAnaConvCorrPhoton&); // not implemented
  ClassDef(AliAnaConvCorrPhoton, 1); // example of analysis

};

#endif
