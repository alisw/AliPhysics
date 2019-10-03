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
  fAxisCent(), 
  fAxisZ(),
  fAxisTrigEta(), 
  fAxisAssEta(), 
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
  fAxisdEta.Set(32, -1.6, 1.6);
  fAxisdEta.SetNameTitle("dEta", "delta eta");

  fAxisTrigEta.SetNameTitle("tEta", "Eta");
  fAxisTrigEta.Set(320, -0.8, 0.8);

  fAxisAssEta.SetNameTitle("aEta", "Eta");
  fAxisAssEta.Set(360, -0.9, 0.9);

  fAxisdPhi.Set(32, -TMath::PiOver2(), 3*TMath::PiOver2());
  fAxisdPhi.SetNameTitle("dPhi", "delta Phi");

  Double_t tptbins[15] = {2.0, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 12.5, 15, 20, 25, 30, 50, 100};
  fAxistPt.Set(13, tptbins);
  fAxistPt.SetNameTitle("tPt", "trigger Pt");

  Double_t cptbins[19] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 12.5, 15, 20, 25, 30, 50, 100};
  fAxiscPt.Set(18, cptbins);
  fAxiscPt.SetNameTitle("cPt", "track Pt");

  fAxisCent.SetNameTitle("centrality", "centrality");
  fAxisCent.Set(1, -999, 999);
  fAxisZ.SetNameTitle("vtxz", "vtxz");
  Double_t zbins[6] = {-10, -5, -1.5, 1.5, 5, 10 }; 
  fAxisZ.Set(5, zbins);

  fAxisIso.Set(1, -0.5, 2.5);
  fAxisIso.SetNameTitle("iso", "isolation");

  fAxesList.AddAt(&fAxisdEta, 0);
  fAxesList.AddAt(&fAxisdPhi, 1);
  fAxesList.AddAt(&fAxistPt, 2);
  fAxesList.AddAt(&fAxiscPt, 3);
  fAxesList.AddAt(&fAxisCent, 4);
  fAxesList.AddAt(&fAxisZ, 5);

  fTrackAxisList.AddAt(&fAxisAssEta, 0);
  fTrackAxisList.AddAt(&fAxistPt, 1);
  fTrackAxisList.AddAt(&fAxiscPt, 2);
  fTrackAxisList.AddAt(&fAxisCent, 3);
  fTrackAxisList.AddAt(&fAxisZ, 4);

  fTrigAxisList.AddAt(&fAxisTrigEta, 0);
  fTrigAxisList.AddAt(&fAxistPt, 1);
  fTrigAxisList.AddAt(&fAxisCent, 2);
  fTrigAxisList.AddAt(&fAxisZ, 3);


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



//_______________________________________________________________________________
void AliAnaConvCorrBase::FillCounters(TObjArray * particles, TObjArray * tracks, Float_t cent, Float_t vtxz) {
  //Fill ME Counters
  const Int_t nbins = fAxistPt.GetNbins();
  Bool_t tmap[nbins];
  for(Int_t ptbin = 0; ptbin < nbins; ptbin++){
    tmap[ptbin] = kFALSE;
  }


  Double_t trackValues[fTrackAxisList.GetSize()];
  trackValues[3] = cent;
  trackValues[4] = vtxz;

  for(Int_t ip = 0; ip < particles->GetEntriesFast(); ip++){
    AliAODConversionParticle * particle = static_cast<AliAODConversionParticle*>(particles->At(ip));

    Int_t tbin = fAxistPt.FindFixBin(particle->Pt());
    if (tbin > 0 && tbin < nbins + 1) {
      if(tmap[tbin - 1] == kTRUE) {
	continue;
      } else {
	tmap[tbin -1 ] = kTRUE;

	if( fTrackAxisList.GetSize() > 5){
	  trackValues[5] = particle->M();
	}

	for(int ij = 0; ij < tracks->GetEntriesFast(); ij++) {
	  AliVTrack * track = static_cast<AliVTrack*>(tracks->UncheckedAt(ij));
	  trackValues[0] = track->Eta();
	  trackValues[1] = particle->Pt();
	  trackValues[2] = track->Pt();
	  fTrackSparse->Fill(trackValues);	
	}
      }
    }
  }
}

//________________________________________________________________
void AliAnaConvCorrBase::CorrelateWithTracks(AliAODConversionParticle * particle, TObjArray * tracks, Int_t const tIDs[4], Float_t cent, Float_t vtxz) {
  //Correlate particle with tracks


  const Int_t nDim = fAxesList.GetSize();
  Double_t dphivalues[nDim];
  dphivalues[4] = cent;
  dphivalues[5] = vtxz;


  Double_t trigValues[fTrigAxisList.GetSize()];
  trigValues[0] = particle->Eta();
  trigValues[1] = particle->Pt();
  trigValues[2] = cent;
  trigValues[3] = vtxz;
  
  if(nDim > 6) {
    dphivalues[6] = particle->M();
    trigValues[4] = particle->M();
  }

  fTrigSparse->Fill(trigValues);

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
    fCorrSparse->Fill(dphivalues);
  }
}

