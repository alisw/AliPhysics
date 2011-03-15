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
/// @brief  Class used to run conversion gamma/pion - hadron/jet analysis



#include "AliAnaConvCorrPion.h"
#include "AliAODTrack.h"
#include "TClonesArray.h"
#include "AliAODConversionParticle.h"

#include <iostream>


using namespace std;
ClassImp(AliAnaConvCorrPion)

//________________________________________________________________________
AliAnaConvCorrPion::AliAnaConvCorrPion() :
AliAnaConvCorrBase("pion_hadron_corr") {
  //consctructor
}
//________________________________________________________________________
AliAnaConvCorrPion::AliAnaConvCorrPion(TString name) :
AliAnaConvCorrBase(name) {
  //consctructor
}


//________________________________________________________________________________
AliAnaConvCorrPion::~AliAnaConvCorrPion() {
  //destructor
}


///_______________________________________________________________________________
void AliAnaConvCorrPion::GetTrackLabels(const AliAODConversionParticle * pion, const TClonesArray * photons, Int_t* trackLabels) {
  ///Get the track labels of the electrons reconstructed as gamma forming the pion

  for(Int_t i = 0; i< 2; i++) {
    AliAODConversionParticle * gamma = dynamic_cast<AliAODConversionParticle*>(photons->At(pion->GetTrackLabel(i)));

    if(gamma) { 
      for(Int_t j = 0; j< 2; j++) {
	cout << "index " << i + j + ((i>=j)?1:0) << " " << gamma->GetTrackLabel(j) << endl;
	trackLabels[ i*2+ j] = gamma->GetTrackLabel(j);
      }

    } else {
      AliGammaConversionAODObject * aodO = dynamic_cast<AliGammaConversionAODObject*>(photons->At(pion->GetTrackLabel(i)));
      if(aodO) {
	trackLabels[i*2] = aodO->GetLabel1();
	trackLabels[i*2 + 1] = aodO->GetLabel2();
      } else {
	cout << "AliAnaConvCorrPion::GetTrackLabels() :: Not good!!!"<<endl;
      }
    }
  }
}
 

///__________________________________________________________________________
void AliAnaConvCorrPion::CorrelateWithHadrons(AliAODConversionParticle * pion, const TClonesArray * tracks,  const Bool_t isolated, const Int_t nSpawn, const Int_t * const spawn) {
  //See header file for documentation

  FillTriggerCounters(pion->Pt(), isolated);


  if (tracks) {
    for(int ij = 0; ij < tracks->GetEntriesFast(); ij++) {
      AliAODTrack * track = dynamic_cast<AliAODTrack*>(tracks->At(ij));
      if(track) {
	if(pion->IsMySpawn(track->GetID(), nSpawn, spawn)) continue;
	
	if (track->Pt() < GetCorrelatedPt() ) continue;

	FillHistograms(pion->Pt(), track->Pt(), GetDPhi(pion->Phi() - track->Phi()), pion->Eta() - track->Eta(), isolated);

      }
    }
  }
}
