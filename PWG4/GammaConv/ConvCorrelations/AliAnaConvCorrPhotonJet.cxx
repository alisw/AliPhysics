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

/// @file   AliAnaConvCorrPhoton.cxx
/// @author Svein Lindal
/// @brief  Class used to fill calculate correlation between photons - jets

#include "AliAnaConvCorrPhotonJet.h"
#include "AliAODTrack.h"
#include "TClonesArray.h"
#include "AliAODConversionParticle.h"
#include "AliAODJet.h"

#include <iostream>
// Gamma - jet correlation analysis task
// Authors: Svein Lindal


using namespace std;
ClassImp(AliAnaConvCorrPhotonJet)

//________________________________________________________________________________
AliAnaConvCorrPhotonJet::AliAnaConvCorrPhotonJet() :
AliAnaConvCorrBase("photonJet") {
  //consctructor
}
//________________________________________________________________________________
AliAnaConvCorrPhotonJet::AliAnaConvCorrPhotonJet(TString name) :
AliAnaConvCorrBase(name) {
  //consctructor
}


//________________________________________________________________________________
AliAnaConvCorrPhotonJet::~AliAnaConvCorrPhotonJet() {
  //destructor
}

///_______________________________________________________________________________
void AliAnaConvCorrPhotonJet::CorrelateWithHadrons(const AliAODConversionParticle * const photon, const TClonesArray * const tracks, const Bool_t isolated) {

  //See header file for documentation
  if (tracks) {
      
    for(int ij = 0; ij < tracks->GetEntriesFast(); ij++) {
      AliAODJet * jet = dynamic_cast<AliAODJet*>(tracks->At(ij));
      if(jet) {
	FillHistograms(photon->Pt(), jet->Pt(), GetDPhi(photon->Phi() - jet->Phi()), photon->Eta() - jet->Eta(), isolated);
      }
    }
  }
}
