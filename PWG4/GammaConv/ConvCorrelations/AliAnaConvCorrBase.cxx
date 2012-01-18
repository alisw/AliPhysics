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
  fAxesList(NULL), 
  fAxistPt(),
  fAxiscPt(), 
  fAxisdEta(), 
  fAxisdPhi(),
  fSparse(NULL)
{
  //Constructor

  SetUpDefaultBins();
  
  fAxesList.SetOwner(kTRUE);


}


//________________________________________________________________________________
AliAnaConvCorrBase::~AliAnaConvCorrBase() {
  ///destructor
}


void AliAnaConvCorrBase::CreateHistograms() {
  CreateBaseHistograms();
}

///________________________________________________________________________________
void AliAnaConvCorrBase::SetUpDefaultBins() {
  //Set up default bins
  Double_t ptbins[19] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 12.5, 15, 20, 25, 30, 50, 100};
  fAxisdEta.Set(160, -1.6, 1.6);
  fAxisdEta.SetNameTitle("dEta", "delta eta");
  fAxesList.AddAt(&fAxisdEta, 0);

  fAxisdPhi.Set(64, -TMath::PiOver2(), 3*TMath::PiOver2());
  fAxisdPhi.SetNameTitle("dPhi", "delta Phi");
  fAxesList.AddAt(&fAxisdPhi, 1);

  fAxistPt.Set(18, ptbins);
  fAxistPt.SetNameTitle("tPt", "trigger Pt");
  fAxesList.AddAt(&fAxistPt, 2);

  fAxiscPt.Set(18, ptbins);
  fAxiscPt.SetNameTitle("cPt", "track Pt");
  fAxesList.AddAt(&fAxiscPt, 3);

  for(int iIso = 0; iIso < 2; iIso++) {
    fHNTriggers[iIso] = NULL;
  }
}

//________________________________________________________________________
void AliAnaConvCorrBase::CreateBaseHistograms() {
  //Create histograms add, to outputlis

  cout << "Creating histograms for "<< GetName() << endl;

  fHistograms = new TList();
  fHistograms->SetOwner(kFALSE);
  fHistograms->SetName(fName);

  for(int iIso = 0; iIso < 2; iIso++) {

    fHNTriggers[iIso] = new TH1F(Form("%s_%s_fNTriggers", fName.Data(), (iIso==0)?"nonIso":"isolated"), 
								 Form("%s_%s_fNTriggers", fName.Data(), (iIso==0)?"nonIso":"isolated"), 
								 fAxistPt.GetNbins(), fAxistPt.GetXbins()->GetArray());
    fHNTriggers[iIso]->Sumw2();
    fHistograms->Add(fHNTriggers[iIso]);
  
  }

  fSparse = CreateSparse(GetName(), GetTitle(), &fAxesList);
  fHistograms->Add(fSparse);

}

///________________________________________________________________________
THnSparseF * AliAnaConvCorrBase::CreateSparse(TString nameString, TString titleString, TList * axesList) {
  //Create sparse
  const Int_t dim = axesList->GetSize();
  
  cout << "dimesion: " << dim << endl;

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
	cout << axes[i]->GetTitle() << endl;
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
  //Print some statistics between each file
  for(Int_t i = 1; i <= fHNTriggers[0]->GetNbinsX(); i++) {
    Int_t nTrig = (Int_t) fHNTriggers[0]->GetBinContent(i+1);
    cout << "triggers: " << nTrig << endl;

  }
}


//_______________________________________________________________________________
void AliAnaConvCorrBase::FillTriggerCounters(const AliAODConversionParticle * particle, Bool_t leading) {
  fHNTriggers[leading]->Fill(particle->Pt());
}


//________________________________________________________________
void AliAnaConvCorrBase::CorrelateWithTracks(AliAODConversionParticle * particle, TObjArray * tracks, Int_t const tIDs[4], Bool_t isolated /*= kFALSE*/) {
  //Correlate particle with tracks

  FillTriggerCounters(particle, isolated);

 Int_t nDim = fAxesList.GetSize();
  Double_t dphivalues[nDim];
  
  for(int ij = 0; ij < tracks->GetEntriesFast(); ij++) {
	AliVTrack * track = static_cast<AliVTrack*>(tracks->UncheckedAt(ij));
	Int_t tid = track->GetID();

	if((tid > 0) && (tid == tIDs[0] || tid == tIDs[1] || tid == tIDs[2] || tid == tIDs[3]) ) {
	  continue;
	}
	
	dphivalues[1] = GetDPhi(particle->Phi() - track->Phi());
	dphivalues[0] = particle->Eta() - track->Eta();
	dphivalues[2] = particle->Pt();
	dphivalues[3] = track->Pt();
	if(nDim > 4) dphivalues[4] = particle->M();
	fSparse->Fill(dphivalues);
  }
}

