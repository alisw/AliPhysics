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

/// @file   AliAnaConvCorrBase.cxx
/// @author Svein Lindal
/// @brief  Base class for analysation of conversion particle - track correlations


#include "AliAnaConvCorrBase.h"
#include "AliAODTrack.h"

#include "TClonesArray.h"
#include "TH1F.h"
#include "TH3.h"
#include "TList.h"
#include "AliAODConversionParticle.h"

#include <iostream>

// Gamma - jet correlation analysis task
// Authors: Svein Lindal


using namespace std;
ClassImp(AliAnaConvCorrBase)

//________________________________________________________________________
AliAnaConvCorrBase::AliAnaConvCorrBase(TString name, TString title = "title") : TNamed(name, title),
  fHistograms(NULL),
  fAxesList(), 
  fTrigAxisList(), 
  fTrackAxisList(),
  fAxistPt(),
  fAxiscPt(), 
  fAxisdEta(), 
  fAxisdPhi(),
  fAxisIso(), 
  fAxisMEEta(), 
  fAxisMEPhi(),
  fCorrSparse(NULL),
  fTrigSparse(NULL),
  fTrackSparse(NULL)
{
  //Constructor
  fAxesList.SetOwner(kFALSE);
  fTrackAxisList.SetOwner(kFALSE);
  fTrigAxisList.SetOwner(kFALSE);

  SetUpDefaultBins();
 
}


//________________________________________________________________________________
AliAnaConvCorrBase::~AliAnaConvCorrBase() {
  ///destructor
}

//________________________________________________________________________________
void AliAnaConvCorrBase::CreateHistograms() {
  CreateBaseHistograms();
}

///________________________________________________________________________________
void AliAnaConvCorrBase::SetUpDefaultBins() {
  //Set up default bins
  fAxisdEta.Set(40, -1.6, 1.6);
  fAxisdEta.SetNameTitle("dEta", "delta eta");

  fAxisdPhi.Set(32, -TMath::PiOver2(), 3*TMath::PiOver2());
  fAxisdPhi.SetNameTitle("dPhi", "delta Phi");

  Double_t tptbins[14] = {3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 12.5, 15, 20, 25, 30, 50, 100};
  fAxistPt.Set(13, tptbins);
  fAxistPt.SetNameTitle("tPt", "trigger Pt");

  Double_t cptbins[18] = {1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 12.5, 15, 20, 25, 30, 50, 100};
  fAxiscPt.Set(17, cptbins);
  fAxiscPt.SetNameTitle("cPt", "track Pt");

  fAxisIso.Set(3, -0.5, 2.5);
  fAxisIso.SetNameTitle("iso", "isolation");

  fAxesList.AddAt(&fAxisdEta, 0);
  fAxesList.AddAt(&fAxisdPhi, 1);
  fAxesList.AddAt(&fAxistPt, 2);
  fAxesList.AddAt(&fAxiscPt, 3);
  fAxesList.AddAt(&fAxisIso, 4);

  fAxisMEEta.Set(320, -0.8, 0.8);
  fAxisMEEta.SetNameTitle("eta", "eta");
  
  fAxisMEPhi.Set(256, 0, TMath::TwoPi());
  fAxisMEPhi.SetNameTitle("phi", "phi");

  fTrackAxisList.AddAt(&fAxisMEEta, 0);
  fTrackAxisList.AddAt(&fAxisMEPhi, 1);
  fTrackAxisList.AddAt(&fAxistPt, 2);
  fTrackAxisList.AddAt(&fAxiscPt, 3);
  fTrackAxisList.AddAt(&fAxisIso, 4);

  fTrigAxisList.AddAt(&fAxisMEEta, 0);
  fTrigAxisList.AddAt(&fAxisMEPhi, 1);
  fTrigAxisList.AddAt(&fAxistPt, 2);
  fTrigAxisList.AddAt(&fAxisIso, 3);


}


//________________________________________________________________________
void AliAnaConvCorrBase::CreateBaseHistograms() {
  //Create histograms add, to outputlis

  //cout << "Creating histograms for "<< GetName() << endl;

  fHistograms = new TList();
  fHistograms->SetOwner(kTRUE);
  fHistograms->SetName(fName);



  fCorrSparse = CreateSparse(GetName(), GetTitle(), &fAxesList);
  fHistograms->Add(fCorrSparse);

  fTrackSparse = CreateSparse(Form("%s_%s", GetName(), "METrack"), Form("%s %s", GetTitle(), "ME Tracks"), &fTrackAxisList);
  fHistograms->Add(fTrackSparse);

  fTrigSparse = CreateSparse(Form("%s_%s", GetName(), "METrig"), Form("%s %s", GetTitle(), "ME Triggers"), &fTrigAxisList);
  fHistograms->Add(fTrigSparse);

}

///________________________________________________________________________
THnSparseF * AliAnaConvCorrBase::CreateSparse(TString nameString, TString titleString, TList * axesList) {
  //Create sparse
  const Int_t dim = axesList->GetSize();
  
  //cout << nameString << " " << titleString << " " <<   "    dimesion: " << dim << endl;

  TAxis * axes[dim];
  Int_t   bins[dim];
  Double_t min[dim];
  Double_t max[dim];

  for(Int_t i = 0; i<dim; i++) {
	TAxis * axis = dynamic_cast<TAxis*>(axesList->At(i));
	if(axis) axes[i] = axis;
	else {
	  cout << "AliAnalysisTaskdPhi::CreateSparse: Error error, all the axes are not present in axis list" << endl;
	  return NULL;
	}
  }

  for(Int_t i = 0; i<dim; i++) {
	//cout << axes[i]->GetTitle() << endl;
 	bins[i] = axes[i]->GetNbins(); 
	min[i] = axes[i]->GetBinLowEdge(1);
	max[i] = axes[i]->GetBinUpEdge(axes[i]->GetNbins());
  }

  THnSparseF * sparse = new THnSparseF(Form("%s", nameString.Data()), 
									   Form("%s", titleString.Data()), 
									   dim, bins, min, max);
  
  for(Int_t i = 0; i<dim; i++) {
	sparse->GetAxis(i)->SetNameTitle(axes[i]->GetName(), axes[i]->GetTitle() );
	if(axes[i]->GetXbins()->GetSize() > 0) {
	  sparse->SetBinEdges(i, axes[i]->GetXbins()->GetArray() );
	}
  }
  return sparse;
}


///____________________________________________________________________________
// void AliAnaConvCorrBase::FillTriggerCounters(Float_t tPt, Bool_t isolated){ 
//   //Fill histogram with trigger counters

//   fHNTriggers[0]->Fill(tPt);
  
//   if(isolated) {
//     fHNTriggers[isolated]->Fill(tPt);
    
//   }
// }

// ///_____________________________________________________________________________
// void AliAnaConvCorrBase::FillHistograms(Float_t tPt, Float_t cPt, Float_t dPhi, Float_t dEta, Bool_t isolated) {
//   //Fill histograms

//   if(dEta) { ;}
//   //fHdPhi[0]->Fill(tPt, cPt, dPhi);
//   if(isolated) {
//     //fHdPhi[isolated]->Fill(tPt, cPt, dPhi);
//   }
// }

//_______________________________________________________________________________

void AliAnaConvCorrBase::PrintStatistics()  { 
  
  // }
}


// //_______________________________________________________________________________
// void AliAnaConvCorrBase::FillTriggerCounters(const AliAODConversionParticle * particle, Int_t leading) {

// }


//________________________________________________________________
void AliAnaConvCorrBase::CorrelateWithTracks(AliAODConversionParticle * particle, TObjArray * tracks, Int_t const tIDs[4], Int_t isolated = 0) {
  //Correlate particle with tracks


   //FillTriggerCounters(particle, isolated);

  Int_t nDim = fAxesList.GetSize();
  Double_t dphivalues[nDim];
  Double_t trackValues[nDim];

  Double_t trigValues[nDim - 1];
  trigValues[0] = particle->Eta();
  trigValues[1] = particle->Phi();
  trigValues[2] = particle->Pt();
  trigValues[3] = isolated;
  
  fTrigSparse->Fill(trigValues);

  if(nDim > 4) {
	dphivalues[5] = particle->M();
	trackValues[5] = particle->M();
	trigValues[4] = particle->M();
  }

  for(int ij = 0; ij < tracks->GetEntriesFast(); ij++) {
	AliVTrack * track = static_cast<AliVTrack*>(tracks->UncheckedAt(ij));
	Int_t tid = track->GetID();

	if((tid > 0) && (tid == tIDs[0] || tid == tIDs[1] || tid == tIDs[2] || tid == tIDs[3]) ) {
	  continue;
	}
	
	dphivalues[0] = particle->Eta() - track->Eta();
	dphivalues[1] = GetDPhi(particle->Phi() - track->Phi());
	dphivalues[2] = particle->Pt();
	dphivalues[3] = track->Pt();
	dphivalues[4] = isolated;

	trackValues[0] = track->Eta();
	trackValues[1] = track->Phi();
	trackValues[2] = particle->Pt();
	trackValues[3] = track->Pt();
	trackValues[4] = isolated;


	
	fCorrSparse->Fill(dphivalues);
	fTrackSparse->Fill(trackValues);
  }
}

