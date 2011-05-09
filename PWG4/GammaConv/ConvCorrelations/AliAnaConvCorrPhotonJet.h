/* This file is property of and copyright                                 *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/// @file   AliAnaConvCorrPhoton.h
/// @author Svein Lindal
/// @brief  Class used to find correlations between photons and jets

#ifndef ALIANACONVCORRPHOTONJET_CXX
#define ALIANACONVCORRPHOTONJET_CXX

#include "AliAnaConvCorrBase.h"

class TH1F;

class AliAODConversionParticle;
class TClonesArray;
class AliAODJet;
class AliAnaConvCorrPhotonJet : public AliAnaConvCorrBase {

public:

  AliAnaConvCorrPhotonJet(); 
  AliAnaConvCorrPhotonJet(TString name); 
  virtual ~AliAnaConvCorrPhotonJet();
  
  //Correlate photon with jets
  virtual void CorrelateWithHadrons(const AliAODConversionParticle * const photon, const TClonesArray * const jets, const Bool_t isolated);
  Double_t ExtractFromJet(AliAODJet * jet, const AliAODConversionParticle * const particle)  const;
  Bool_t IsParticleInJet(AliAODJet * jet, Int_t nTracks, Int_t * trackIds) const;
  void DoJetAnalysisGamma(AliAODJet * jet, const TClonesArray * const photons, const  TClonesArray *const pions ) const;
  void CreateHistograms();

 private:

  AliAnaConvCorrPhotonJet(const AliAnaConvCorrPhotonJet&); // not implemented
  AliAnaConvCorrPhotonJet& operator=(const AliAnaConvCorrPhotonJet&); // not implemented
  ClassDef(AliAnaConvCorrPhotonJet, 1); // 

  TH1F * fhPtFracGamma;// = new TH1F("fhPtFracGamma", "fhPtFracGamma", 100, 0, 10);
  TH1F * fhPtFracPion;// = new TH1F("fhPtFracPion", "fhPtFracPion", 100, 0, 10);

};

#endif
