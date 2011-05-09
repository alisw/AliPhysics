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
AliAnaConvCorrBase::AliAnaConvCorrBase(TString name) : TObject(),
  fName(name),
  fHistograms(NULL),
  fNPhiBins(32),
  fdPhiBins(NULL),
  fPtBins(NULL)
{
  //Constructor
  fPtBins = new TArrayD(41);
  for(Int_t i = 0; i < 10; i++) {
    fPtBins->SetAt(i*0.5, i);
    fPtBins->SetAt(5 + i*0.5, i + 10);
    
    fPtBins->SetAt(10. + i, i+20);
    fPtBins->SetAt(20. + 2*i, i+30);
  }

  fPtBins->SetAt(50., fPtBins->GetSize() -1);

  fdPhiBins = new TArrayD(fNPhiBins + 1);
  for(Int_t i = 0; i < fNPhiBins+1; i++) {
    fdPhiBins->SetAt(-TMath::PiOver2() + i*TMath::TwoPi()/fNPhiBins, i);
  }
  
  for(int iIso = 0; iIso < 2; iIso++) {
    fHdPhi[iIso] = NULL;
    fHNTriggers[iIso] = NULL;
  }



}


//________________________________________________________________________________
AliAnaConvCorrBase::~AliAnaConvCorrBase() {
  ///destructor
  if(fPtBins)
    delete fPtBins;
  fPtBins = NULL;
  
  if(fdPhiBins)
    delete fdPhiBins;
  fdPhiBins = NULL;
  

}


void AliAnaConvCorrBase::CreateHistograms() {
  CreateBaseHistograms();
}

//________________________________________________________________________
void AliAnaConvCorrBase::CreateBaseHistograms() {
  //Create histograms add, to outputlis

  cout << "Createing histograms for "<< fName.Data() << endl;

  fHistograms = new TList();
  fHistograms->SetOwner(kFALSE);
  fHistograms->SetName(fName);

  for(int iIso = 0; iIso < 2; iIso++) {

    fHdPhi[iIso] = new TH3F(Form("%s_%s_dPhi", fName.Data(),  (iIso==0)?"nonIso":"isolated"),
			    Form("%s_%s_dPhi", fName.Data(),  (iIso==0)?"nonIso":"isolated"),
			    fPtBins->GetSize() -1, fPtBins->GetArray(),
			    fPtBins->GetSize() - 1, fPtBins->GetArray(),
			    fdPhiBins->GetSize() - 1, fdPhiBins->GetArray());
    fHdPhi[iIso]->Sumw2();
    
    fHistograms->Add(fHdPhi[iIso]);

    fHNTriggers[iIso] = new TH1F(Form("%s_%s_fNTriggers", fName.Data(), (iIso==0)?"nonIso":"isolated"), 
				 Form("%s_%s_fNTriggers", fName.Data(), (iIso==0)?"nonIso":"isolated"), 
				 fPtBins->GetSize() - 1, fPtBins->GetArray());
    fHNTriggers[iIso]->Sumw2();
    fHistograms->Add(fHNTriggers[iIso]);
  
  }

}


///____________________________________________________________________________
void AliAnaConvCorrBase::FillTriggerCounters(Float_t tPt, Bool_t isolated){ 
  //Fill histogram with trigger counters

  fHNTriggers[0]->Fill(tPt);
  
  if(isolated) {
    fHNTriggers[isolated]->Fill(tPt);
    
  }
}

///_____________________________________________________________________________
void AliAnaConvCorrBase::FillHistograms(Float_t tPt, Float_t cPt, Float_t dPhi, Float_t dEta, Bool_t isolated) {
  //Fill histograms

  if(dEta) { ;}
  fHdPhi[0]->Fill(tPt, cPt, dPhi);
  if(isolated) {
    fHdPhi[isolated]->Fill(tPt, cPt, dPhi);
  }
}

//_______________________________________________________________________________

void AliAnaConvCorrBase::PrintStatistics()  { 
  //Print some statistics between each file
  for(Int_t i = 1; i <= fHNTriggers[0]->GetNbinsX(); i++) {
    Int_t nTrig = (Int_t) fHNTriggers[0]->GetBinContent(i+1);
    cout << "triggers: " << nTrig << endl;

  }
}
