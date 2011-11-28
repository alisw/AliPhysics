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

#include "TRefArray.h"
#include "TH1F.h"
#include <iostream>
// Gamma - jet correlation analysis task
// Authors: Svein Lindal


using namespace std;
ClassImp(AliAnaConvCorrPhotonJet)

//________________________________________________________________________________
AliAnaConvCorrPhotonJet::AliAnaConvCorrPhotonJet() :
AliAnaConvCorrBase("photonJet"),
  fhPtFracGamma(NULL), 
  fhPtFracPion(NULL)
{
  //consctructor
}
//________________________________________________________________________________
AliAnaConvCorrPhotonJet::AliAnaConvCorrPhotonJet(TString name) :
  AliAnaConvCorrBase(name), 
  fhPtFracGamma(NULL), 
  fhPtFracPion(NULL)
{
  //consctructor
}


//________________________________________________________________________________
AliAnaConvCorrPhotonJet::~AliAnaConvCorrPhotonJet() {
  //destructor
}


void AliAnaConvCorrPhotonJet::CreateHistograms() {
  
  CreateBaseHistograms();
  fhPtFracGamma = new TH1F("fhPtFracGamma", "fhPtFracGamma", 100, 0, 10);
  GetHistograms()->Add(fhPtFracGamma);
  fhPtFracPion = new TH1F("fhPtFracPion", "fhPtFracPion", 100, 0, 10);
  GetHistograms()->Add(fhPtFracPion);

}

void AliAnaConvCorrPhotonJet::DoJetAnalysisGamma(AliAODJet * jet, const TClonesArray * const photons, const  TClonesArray *const pions ) const{
  
  Int_t trackIDs[4] = {-9, -9, -9, -9};

  for(Int_t i = 0; i < photons->GetEntriesFast(); i++) {
    AliAODConversionParticle * photon = dynamic_cast<AliAODConversionParticle*>(photons->At(i));
    if(photon) {
      trackIDs[0] = photon->GetLabel1();
      trackIDs[1] = photon->GetLabel2();
      if(IsParticleInJet(jet, 2, trackIDs)){
	fhPtFracGamma->Fill(photon->Pt()/jet->Pt());
      }
    }
  }



  for(Int_t i = 0; i < pions->GetEntriesFast(); i++) {
    AliAODConversionParticle * pion = dynamic_cast<AliAODConversionParticle*>(pions->At(i));
    if(pion) {
      //pion->GetGrandChildren(photons, trackIDs);
      if(IsParticleInJet(jet, 4, trackIDs)){
	fhPtFracPion->Fill(pion->Pt()/jet->Pt());
      }
    }
  }
  
  
  

  
}

//________________________________________________________________________________
Bool_t AliAnaConvCorrPhotonJet::IsParticleInJet(AliAODJet * jet, Int_t nTracks, Int_t * trackIds) const {

  Int_t mTracks = 0;
  TRefArray * refTracks = jet->GetRefTracks();
  for(Int_t jRef = 0; jRef < refTracks->GetEntriesFast(); jRef++) {
    AliAODTrack * track = dynamic_cast<AliAODTrack*>(refTracks->At(jRef));
    if(track) {
      for(Int_t it = 0; it < nTracks; it++) {
	if (track->GetID() == trackIds[it]) {
	  mTracks++;
	}
      }
    }
  }

  //cout <<mTracks << " " << (mTracks > 1) << endl;
  return (mTracks > 1);
}


//________________________________________________________________________________
Double_t AliAnaConvCorrPhotonJet::ExtractFromJet(AliAODJet * jet, const AliAODConversionParticle * const particle) const {
  
  Float_t jetPt = jet->Pt();
  cout << "Jet pt before and after: " << jetPt << "    ";

  TRefArray * refTracks = jet->GetRefTracks();
  for(Int_t jRef = 0; jRef < refTracks->GetEntriesFast(); jRef++) {
    AliAODTrack * track = dynamic_cast<AliAODTrack*>(refTracks->At(jRef));
    if(track) {
      if (track->GetID() == particle->GetLabel1() || track->GetID() == particle->GetLabel2()) {
	cout << " - " << track->Pt() << "  ";
	jetPt = jetPt - track->Pt();
      } else {
	//cout << track->Pt() << endl;
      }
    } else {
      cout <<"FUUUUUUUUUUUUUUUUUCJK"<<endl;
    }
  }
  
  cout << jetPt << endl;
  return jetPt;
}



///_______________________________________________________________________________
void AliAnaConvCorrPhotonJet::CorrelateWithHadrons(const AliAODConversionParticle * const photon, const TClonesArray * const jets, const Bool_t isolated) {
  FillTriggerCounters(photon->Pt(), isolated);
  //See header file for documentation
  if (jets) {
    for(int ij = 0; ij < jets->GetEntriesFast(); ij++) {
      AliAODJet * jet = dynamic_cast<AliAODJet*>(jets->At(ij));
      if(jet) {
	Double_t jetPt = ExtractFromJet(jet, photon);
	FillHistograms(photon->Pt(), jetPt, GetDPhi(photon->Phi() - jet->Phi()), photon->Eta() - jet->Eta(), isolated);
      }
    }
  }
}
