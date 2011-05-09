/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Author: Svein Lindal <slindal@fys.uio.no>                      *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/// @file   AliAnaConvCorrPion.cxx
/// @author Svein Lindal
/// @brief  Class used to fill calculate correlation between photons - jets

#include "AliAnaConvCorrPionJet.h"
#include "AliAODTrack.h"
#include "TClonesArray.h"
#include "AliAODConversionParticle.h"
#include "AliAODJet.h"

#include <iostream>
// Gamma - jet correlation analysis task
// Authors: Svein Lindal


using namespace std;
ClassImp(AliAnaConvCorrPionJet)

//________________________________________________________________________________
AliAnaConvCorrPionJet::AliAnaConvCorrPionJet() :
AliAnaConvCorrBase("photonJet") {
  //consctructor
}
//________________________________________________________________________________
AliAnaConvCorrPionJet::AliAnaConvCorrPionJet(TString name) :
AliAnaConvCorrBase(name) {
  //consctructor
}


//________________________________________________________________________________
AliAnaConvCorrPionJet::~AliAnaConvCorrPionJet() {
  //destructor
}

///_______________________________________________________________________________
void AliAnaConvCorrPionJet::CorrelateWithHadrons(const AliAODConversionParticle * const pion, const TClonesArray * const jets, const Bool_t isolated) {

  FillTriggerCounters(pion->Pt(), isolated);

  //See header file for documentation
  if (jets) {
      
    for(int ij = 0; ij < jets->GetEntriesFast(); ij++) {
      AliAODJet * jet = dynamic_cast<AliAODJet*>(jets->At(ij));
      if(jet) {
	FillHistograms(pion->Pt(), jet->Pt(), GetDPhi(pion->Phi() - jet->Phi()), pion->Eta() - jet->Eta(), isolated);
      }
    }
  }
}
