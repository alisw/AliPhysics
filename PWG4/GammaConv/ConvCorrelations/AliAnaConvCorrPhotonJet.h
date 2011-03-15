/* This file is property of and copyright                                 *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/// @file   AliAnaConvCorrPhoton.h
/// @author Svein Lindal
/// @brief  Class used to find correlations between photons and jets

#ifndef ALIANACONVCORRPHOTONJET_CXX
#define ALIANACONVCORRPHOTONJET_CXX

#include "AliAnaConvCorrBase.h"

class AliAODConversionParticle;
class TClonesArray;

class AliAnaConvCorrPhotonJet : public AliAnaConvCorrBase {

public:

  AliAnaConvCorrPhotonJet(); 
  AliAnaConvCorrPhotonJet(TString name); 
  virtual ~AliAnaConvCorrPhotonJet();
  
  //Correlate photon with jets
  virtual void CorrelateWithHadrons(const AliAODConversionParticle * const photon, const TClonesArray * const jets, const Bool_t isolated);
  
 private:

  AliAnaConvCorrPhotonJet(const AliAnaConvCorrPhotonJet&); // not implemented
  AliAnaConvCorrPhotonJet& operator=(const AliAnaConvCorrPhotonJet&); // not implemented
  ClassDef(AliAnaConvCorrPhotonJet, 1); // 

};

#endif
