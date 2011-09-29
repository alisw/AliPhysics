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
#include "AliAODConversionPhoton.h"
#include "THnSparse.h"
#include "TH2F.h"

#include <iostream>


using namespace std;
ClassImp(AliAnaConvCorrPion)

//________________________________________________________________________________
AliAnaConvCorrPion::AliAnaConvCorrPion() :
AliAnaConvCorrBase("pion_hadron_corr"), 
  fhdPhiVsInvMassPi0(NULL), 
  fhdPhiVsInvMassEta(NULL), 
  fhPtVsInvMass(NULL)
{
  //consctructor
}
//________________________________________________________________________________
AliAnaConvCorrPion::AliAnaConvCorrPion(TString name) :
  AliAnaConvCorrBase(name),
  fhdPhiVsInvMassPi0(NULL), 
  fhdPhiVsInvMassEta(NULL), 
  fhPtVsInvMass(NULL)
{
  //consctructor
}


//________________________________________________________________________________
AliAnaConvCorrPion::~AliAnaConvCorrPion() {
  //destructor
}

///________________________________________________________________________________
void AliAnaConvCorrPion::CreateHistograms() {
  //Create histograms
  CreateBaseHistograms();

  const Int_t dim = 4;
  Int_t bins[dim] = {200, 200, 32, 14}; 
  Double_t min[dim] = {0, 0, -TMath::PiOver2(), 0.1};
  Double_t max[dim] = {100, 100, 3*TMath::PiOver2(), 0.17};

  fhdPhiVsInvMassPi0 = new THnSparseF("fhdPhiVsInvMassPi0", "fhdPhiVsInvMassPi0", dim, bins, min, max);

  min[3] = 450;
  max[3] = 650;
  bins[3] = 20;
  fhdPhiVsInvMassEta = new THnSparseF("fhdPhiVsInvMassEta", "fhdPhiVsInvMassEta", dim, bins, min, max);
  fhPtVsInvMass = new TH2F("fhPtVsInvMass", "Pt Vs inv mass", GetTriggerBins()->GetSize() -1, GetTriggerBins()->GetArray(), 400, 0, 1);

  GetHistograms()->Add(fhPtVsInvMass);
  GetHistograms()->Add(fhdPhiVsInvMassPi0);
  GetHistograms()->Add(fhdPhiVsInvMassEta);

}

///________________________________________________________________________________
void AliAnaConvCorrPion::GetTrackLabels(const AliAODConversionPhoton * pion, const TClonesArray * photons, Int_t* trackLabels) {
  ///Get the track labels of the electrons reconstructed as gamma forming the pion

  for(Int_t i = 0; i< 2; i++) {
    AliAODConversionPhoton * gamma = dynamic_cast<AliAODConversionPhoton*>(photons->At(pion->GetTrackLabel(i)));

    if(gamma) { 
      for(Int_t j = 0; j< 2; j++) {
	cout << "index " << i + j + ((i>=j)?1:0) << " " << gamma->GetTrackLabel(j) << endl;
	trackLabels[ i*2+ j] = gamma->GetTrackLabel(j);
      }
    }
  }
}
 

///________________________________________________________________________________
void AliAnaConvCorrPion::CorrelateWithHadrons(AliAODConversionPhoton * pion, const TClonesArray * tracks,  const Bool_t isolated, const Int_t nSpawn, const Int_t * const spawn) {
  //See header file for documentation

  fhPtVsInvMass->Fill(pion->Pt(), pion->M());
  FillTriggerCounters(pion->Pt(), isolated);
  
  if (tracks) {
    for(int ij = 0; ij < tracks->GetEntriesFast(); ij++) {
      AliAODTrack * track = dynamic_cast<AliAODTrack*>(tracks->At(ij));
      if(track) {
	//if(pion->IsMySpawn(track->GetID(), nSpawn, spawn)) continue;
	
	//	if (track->Pt() < GetCorrelatedPt() ) continue;
	Double_t x[4] = {pion->Pt(), track->Pt(), GetDPhi(pion->Phi() - track->Phi()), TMath::Abs(pion->M()) };
	if(  (pion->M() > 0.1) &&  (pion->M() < 0.17) ){
	  fhdPhiVsInvMassPi0->Fill(x); 
	} else if ((pion->M() > 0.5) &&  (pion->M() < 0.6 ) ) {
	  fhdPhiVsInvMassEta->Fill(x); 
	}
	FillHistograms(pion->Pt(), track->Pt(), GetDPhi(pion->Phi() - track->Phi()), pion->Eta() - track->Eta(), isolated);
      }
    }
  }
}
