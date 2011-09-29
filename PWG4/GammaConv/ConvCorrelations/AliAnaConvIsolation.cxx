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

/// @file   AliAnalysisTaskGammaJet.cxx
/// @author Svein Lindal
/// @brief  Class used to run isolation studies of conversion gamma/pions


#include "AliAnaConvIsolation.h"
#include "AliAODTrack.h"

#include "TObject.h"
#include "TClonesArray.h"
#include "TH2F.h"
#include "TList.h"
#include "AliAODConversionPhoton.h"


#include <iostream>

// Gamma - jet correlation analysis task
// Author: Svein Lindal


using namespace std;
ClassImp(AliAnaConvIsolation)

//________________________________________________________________________
AliAnaConvIsolation::AliAnaConvIsolation () : TObject(),
  fIsoCurve(NULL),
  fCurveFunction("0.1*x"),
  fConeSize(0), 
  fMinPt(0.1),
  fMaxPtThreshold(0),
  fSumPtThreshold(0),
  fMaxPtFraction(0),
  fSumPtFraction(0),
  fHistograms(NULL),
  fHistogramMaxPt(50)  
{
  //Constructor
  for(Int_t i = 0; i < 2; i++){
    fhMaxPtInCone[i] = NULL;
    fhSumPtInCone[i] = NULL;
    fhSumPtVsMaxPt[i] = NULL;
    fhPtCandidates[i] = NULL;
    fhTrackMult[i] = NULL;
  }
  
  
}

//________________________________________________________________________
AliAnaConvIsolation::AliAnaConvIsolation(Float_t coneSize, Float_t maxPtThreshold, Float_t sumPtThreshold, Float_t maxPtFraction, Float_t sumPtFraction) :
  TObject(), 
  fIsoCurve(NULL),
  fCurveFunction("0.1*x"),
  fConeSize(coneSize), 
  fMinPt(0.1),
  fMaxPtThreshold(maxPtThreshold),
  fSumPtThreshold(sumPtThreshold),
  fMaxPtFraction(maxPtFraction),
  fSumPtFraction(sumPtFraction),
  fHistograms(NULL),
  fHistogramMaxPt(50)
{
  //Constructor
  for(Int_t i = 0; i < 2; i++){
    fhMaxPtInCone[i] = NULL;
    fhSumPtInCone[i] = NULL;
    fhSumPtVsMaxPt[i] = NULL;
    fhPtCandidates[i] = NULL;
    fhTrackMult[i] = NULL;
  }
}


//________________________________________________________________________________
AliAnaConvIsolation::~AliAnaConvIsolation() {
  //Destructor

}


//________________________________________________________________________
void AliAnaConvIsolation::CreateHistograms()
{
	//Create histograms
  if(!fHistograms) fHistograms = new TList();
  fHistograms->SetName(Form("Isolation_histo_cone_%f_maxPt_%f_sumPt_%f", fConeSize, fSumPtThreshold, fMaxPtThreshold));

  fIsoCurve = new TF1("Isolation_curve", fCurveFunction.Data(), 0, 100);
  fHistograms->Add(fIsoCurve);

  cout << "Creatin isolation histograms conesize :" << fConeSize << endl;

  //Create histograms add, to outputlis
  for(Int_t i = 0; i < 2; i++){
    fhMaxPtInCone[i] = new TH2F(Form("fhMaxPtInCone_%s_%f", (i==0)? "nonIso" : "isolated", fConeSize), 
				Form("Max pt nonIso particle in cone %f vs candidate pt", fConeSize), 
				200, 0, fHistogramMaxPt, 200, 0, fHistogramMaxPt);
    fHistograms->Add(fhMaxPtInCone[i]);
    
    fhSumPtInCone[i] = new TH2F(Form("fhSumPtInCone_%s_%f",  (i==0)? "nonIso" : "isolated", fConeSize), 
			     Form("Sum pt in cone %f vs candidate pt %s", fConeSize,  (i==0)? "nonIsoay" : "isolated"), 
			     200, 0, fHistogramMaxPt, 200, 0, fHistogramMaxPt);
    fHistograms->Add(fhSumPtInCone[i]);

    fhSumPtVsMaxPt[i] = new TH2F(Form("fhSumPtVsMaxPt_%s_%f",  (i==0)? "nonIso" : "isolated", fConeSize), 
				 Form("fhSumPtVsMaxPt_%s_%f",  (i==0)? "nonIso" : "isolated", fConeSize), 
				 200, 0, fHistogramMaxPt, 200, 0, fHistogramMaxPt);
    fHistograms->Add(fhSumPtVsMaxPt[i]);
  }

  for(Int_t iIso = 0; iIso < 2; iIso++){
    fhPtCandidates[iIso] = new TH1F(Form("fhPtCandidates_%s_%f", (iIso==0)? "nonIso" : "isolated", fConeSize),
				    Form("Pt of %s candidates, cone size %f", (iIso==0)? "nonIsolated" : "isolated", fConeSize),
				 200, 0 , fHistogramMaxPt);
    fhPtCandidates[iIso]->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    fhPtCandidates[iIso]->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
    fhPtCandidates[iIso]->SetMarkerStyle(kFullCircle);
    fHistograms->Add(fhPtCandidates[iIso]);
  }

  // for(int iDec = 0; iDec < 2; iDec++) {
  //   for(int iIso = 0; iIso < 2; iIso++) {

  //     fHistSumPt[iIso][iDec] = new TH1F(Form("fHistSumPt_%f_%s_%s", fConeSize, (iIso==0)?"nonIso":"isolated" , (iDec==0)?"noDec":"decay" ),  
  // 					Form("P_{T} distribution cone %f %s %s", fConeSize, (iIso==0)?"nonIso":"isolated" , (iDec==0)?"noDec":"decay" ), 
  // 				     150, 0.1, 50);
  //     fHistSumPt[iIso][iDec]->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  //     fHistSumPt[iIso][iDec]->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  //     fHistSumPt[iIso][iDec]->SetMarkerStyle(kFullCircle);
  //     fHistograms->Add(fHistSumPt[iIso][iDec]);
     
      
  //     fHistMaxPt[iIso][iDec] = new TH1F(Form("fHistMaxPt_%f_%s_%s", fConeSize, (iIso==0)?"nonIso":"isolated" , (iDec==0)?"noDec":"decay" ),  
  // 					Form("P_{T} distribution cone %f %s %s", fConeSize, (iIso==0)?"nonIso":"isolated" , (iDec==0)?"noDec":"decay" ), 
  // 				     150, 0.1, 50);
  //     fHistMaxPt[iIso][iDec]->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  //     fHistMaxPt[iIso][iDec]->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  //     fHistMaxPt[iIso][iDec]->SetMarkerStyle(kFullCircle);
  //     fHistograms->Add(fHistMaxPt[iIso][iDec]);
     
 
  //   }
  // }


  for(Int_t iIso = 0; iIso < 2; iIso++){
    fhTrackMult[iIso] = new TH1F(Form("fhTrackMult_%s_%f", (iIso==0)? "nonIso" : "isolated", fConeSize),
				    Form("Pt of %s candidates, cone size %f", (iIso==0)? "nonIsolated" : "isolated", fConeSize),
				 150, 0 , 150);
    fhTrackMult[iIso]->GetXaxis()->SetTitle("n tracks in event");
    fhTrackMult[iIso]->GetYaxis()->SetTitle("dN/dNTracks");
    fhTrackMult[iIso]->SetMarkerStyle(kFullCircle);
    fHistograms->Add(fhTrackMult[iIso]);
  }

}


//_________________________________________________________________________
Bool_t AliAnaConvIsolation::IsIsolated(AliAODConversionPhoton * particle, const TClonesArray * const tracks, const Int_t nSpawn, const Int_t * const spawn, Bool_t &leading) {
  //See header file for documentation

  leading = kTRUE;

  Float_t ptSum = 0.0;
  Float_t ptMax = 0.0;

  for(int it = 0; it < tracks->GetEntriesFast(); it++) {
    AliAODTrack * track = dynamic_cast<AliAODTrack*>(tracks->At(it));

    if (track) {

      if(track->Pt() < fMinPt) continue;

			if(nSpawn && spawn) { ;}

      ///Ignore tracks that are grandchildren of pion
      // if ( particle->IsMySpawn(track->GetID(), nSpawn, spawn)) 
      // 	continue;
      
      
      if ( IsInCone(particle->Eta() - track->Eta(), particle->Phi() - track->Phi(), fConeSize) ) {
	ptSum += track->Pt();
        if(track->Pt() > ptMax) ptMax = track->Pt();
	if(track->Pt() > particle->Pt()) leading = kFALSE;
      }
    } else {
      cout << "Bad track"<<endl;
    }
  }
  

  Bool_t isolated = EvaluateIsolationCriteria( ptSum, particle->Pt());

  FillHistograms(particle->Pt(), ptMax, ptSum, isolated, tracks->GetEntriesFast());

  return isolated;

}


//_________________________________________________________________________
Bool_t AliAnaConvIsolation::IsIsolated(const AliAODConversionPhoton * const particle, const TClonesArray * const tracks, Bool_t &leading ) {
  //See header file for documentation

  leading = kTRUE;

  Float_t ptSum = 0.0;
  Float_t ptMax = 0.0;

  for(int it = 0; it < tracks->GetEntriesFast(); it++) {
  
    AliAODTrack * track = dynamic_cast<AliAODTrack*>(tracks->At(it));
   
    if (track) {
    
      if(track->Pt() < fMinPt) continue;
      
      if ( (track->GetID() == particle->GetTrackLabel(0)) || track->GetID() == particle->GetTrackLabel(1) )  
	continue;
      
      if ( IsInCone(particle->Eta() - track->Eta(), particle->Phi() - track->Phi(), fConeSize) ) {
	ptSum += track->Pt();
	if(track->Pt() > ptMax) ptMax = track->Pt();
	if(track->Pt() > particle->Pt()) leading = kFALSE;
      }
    } else {
      cout << "Bad track"<<endl;
    }
  }
  

  Bool_t isolated = EvaluateIsolationCriteria( ptSum, particle->Pt());

  FillHistograms(particle->Pt(), ptMax, ptSum, isolated, tracks->GetEntriesFast());

  
  return isolated;

}

///___________________________________________________________________________
void AliAnaConvIsolation::FillHistograms(Float_t pt, Float_t ptMax, Float_t ptSum, Bool_t isolated, Int_t nTracks) {
  //Fill histograms
  fhMaxPtInCone[isolated]->Fill(pt, ptMax);
  fhSumPtInCone[isolated]->Fill(pt, ptSum);
  fhSumPtVsMaxPt[isolated]->Fill(ptMax, ptSum);
  fhPtCandidates[isolated]->Fill(pt);
  fhTrackMult[isolated]->Fill(nTracks);
}


///_____________________________________________________________________________
Bool_t AliAnaConvIsolation::EvaluateIsolationCriteria(Float_t ptSum, Float_t pt) const {
  //Evaluate isolation criteria
  return (ptSum < fIsoCurve->Eval(pt));

}
