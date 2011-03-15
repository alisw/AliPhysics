/* This file is property of and copyright                                 *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/// @file   AliAnaConvCorrPion.h
/// @author Svein Lindal
/// @brief  Class used to find correlations between pions and charged tracks


#ifndef ALIANACONVCORRPION_CXX
#define ALIANACONVCORRPION_CXX

#include "AliAnaConvCorrBase.h"

class AliAODConversionParticle;
class TClonesArray;

class AliAnaConvCorrPion : public AliAnaConvCorrBase {

public:

  AliAnaConvCorrPion(); 
  AliAnaConvCorrPion(TString name);
  virtual ~AliAnaConvCorrPion();

  //Correlate pions with charged tracks
  virtual void CorrelateWithHadrons(AliAODConversionParticle * pion, const TClonesArray * tracks, const Bool_t isolated, const Int_t nSpawn, const Int_t * const spawn );
  
 private:

  //Get array of track labels of the 4 decay electrons (2gamma * 2 electrons)
  void GetTrackLabels(const AliAODConversionParticle * pion, const TClonesArray * photons, Int_t* trackLabels);

  AliAnaConvCorrPion(const AliAnaConvCorrPion&); // not implemented
  AliAnaConvCorrPion& operator=(const AliAnaConvCorrPion&); // not implemented
  ClassDef(AliAnaConvCorrPion, 1); // example of analysis

};

#endif
