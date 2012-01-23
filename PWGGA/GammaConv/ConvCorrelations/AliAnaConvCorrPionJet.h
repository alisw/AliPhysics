/* This file is property of and copyright                                 *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/// @file   AliAnaConvCorrPion.h
/// @author Svein Lindal
/// @brief  Class used to find correlations between photons and jets

#ifndef ALIANACONVCORRPIONJET_CXX
#define ALIANACONVCORRPIONJET_CXX

#include "AliAnaConvCorrBase.h"

class AliAODConversionParticle;
class TClonesArray;

class AliAnaConvCorrPionJet : public AliAnaConvCorrBase {

public:

  AliAnaConvCorrPionJet(); 
  AliAnaConvCorrPionJet(TString name); 
  virtual ~AliAnaConvCorrPionJet();
  
  //Correlate photon with jets
  virtual void CorrelateWithHadrons(const AliAODConversionParticle * const photon, const TClonesArray * const jets, const Bool_t isolated);
  
 private:

  AliAnaConvCorrPionJet(const AliAnaConvCorrPionJet&); // not implemented
  AliAnaConvCorrPionJet& operator=(const AliAnaConvCorrPionJet&); // not implemented
  ClassDef(AliAnaConvCorrPionJet, 1); // 

};

#endif
