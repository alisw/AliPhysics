/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Satyajit Jena      .                                           *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*-------------------------------------------------------------------------
 *
 *             AliEbyEChargeFluctuationAnalysis base class
 *           This class deals with the Charge   fluctuation 
 *                Origin: Satyajit Jena <sjena@cern.ch>
 *
 ------------------------------------------------------------------------*/

#include <Riostream.h>
#include <TFile.h>
#include <TSystem.h>
#include <TF1.h>
#include <TH3F.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TH1I.h>
#include <TParticle.h>
#include <TList.h>

#include <AliExternalTrackParam.h>
#include <AliAODEvent.h>
#include <AliESDEvent.h>
#include <AliESDtrack.h>
#include <AliAODTrack.h>
#include <AliMCParticle.h>
#include <AliPID.h>
#include <AliStack.h>
#include <AliCFContainer.h>
#include <AliCFEffGrid.h>
#include <AliCFDataGrid.h>
#include <AliTPCPIDResponse.h>
#include <AliESDpid.h>
#include "AliEbyEEventBase.h"
#include "AliEbyEChargeFluctuationAnalysis.h"


ClassImp(AliEbyEChargeFluctuationAnalysis)

//____________________________________________________________________//
AliEbyEChargeFluctuationAnalysis::AliEbyEChargeFluctuationAnalysis() : 
  TObject(), 
  fListCFQA(0),
  fListMeasureCF(0),
  fEbyEBase(0)
{
  InitHistos();
  //Default constructor
}


//____________________________________________________________________//
AliEbyEChargeFluctuationAnalysis::~AliEbyEChargeFluctuationAnalysis() {
  //Default destructor
  if(fEbyEBase) delete fEbyEBase;
  if(fListCFQA) delete fListCFQA;
  if(fListMeasureCF) delete fListMeasureCF;

}

//____________________________________________________________________//
void AliEbyEChargeFluctuationAnalysis::InitHistos() {

  fListCFQA = new TList();
  fListCFQA->SetName("MFQaList");

  fListMeasureCF = new TList();
  fListMeasureCF->SetName("MFMeasureList");

  TH1F *fEvtStat = new TH1F("hEventCounterMH"," Event Counters for MF Analysis", 20,0,20);
  fListCFQA->Add(fEvtStat); // --:0


  // Int_t kCentBin = fEbyEBase->GetCentralityBin();
  Int_t kCentBin = 50;  
Int_t kMultBin = 2000;
  Double_t kMultMin = 0.0;
  Double_t kMultMax = 2000;  

  TH2F *fhCentralityV0 = new TH2F("hCentralityV0", "Total Charge Fuctuation", kCentBin,0, (Double_t)kCentBin, 3000,0,30000);
  fhCentralityV0->GetXaxis()->SetTitle(" #-");
  fhCentralityV0->GetYaxis()->SetTitle("Frequency");
  
  fListCFQA->Add(fhCentralityV0); // --:1

  TH2F *fNChargeFluc = new TH2F("hNChargeFluc", "Total Charge Fuctuation",kCentBin,0, (Double_t)kCentBin, kMultBin,kMultMin, kMultMax);
  fNChargeFluc->GetXaxis()->SetTitle(" #-");
  fNChargeFluc->GetYaxis()->SetTitle("Frequency");
  fListMeasureCF->Add(fNChargeFluc); // -->:0

  TH2F *fNpChargeFluc = new TH2F("hNpChargeFluc", "Possitive Charge Fuctuation",kCentBin,0, (Double_t)kCentBin, kMultBin,kMultMin, kMultMax);
  fNpChargeFluc->GetXaxis()->SetTitle(" #-");
  fNpChargeFluc->GetYaxis()->SetTitle("Frequency");
  fListMeasureCF->Add(fNpChargeFluc); // -->:1
  
  TH2F *fNnChargeFluc = new TH2F("hNnChargeFluc", "Negative Charge Fuctuation", kCentBin,0, (Double_t)kCentBin, kMultBin,kMultMin, kMultMax);
  fNnChargeFluc->GetXaxis()->SetTitle(" #-");
  fNnChargeFluc->GetYaxis()->SetTitle("Frequency");
  fListMeasureCF->Add(fNpChargeFluc); // -->:2
  
  TH2F *fNnbpFluc = new TH2F("hNnbpFluc", "[-/+] Fuctuation",kCentBin,0, (Double_t)kCentBin, 100,0,10); 
  fNnbpFluc->GetXaxis()->SetTitle(" #-");
  fNnbpFluc->GetYaxis()->SetTitle("Frequency");
  fListMeasureCF->Add(fNnbpFluc); // -->:3
  
  TH2F *fNpbnFluc = new TH2F("hNpbnFluc", "[+/-] Fuctuation",kCentBin,0, (Double_t)kCentBin, 100,0,10);
  fNpbnFluc->GetXaxis()->SetTitle(" #-");
  fNpbnFluc->GetYaxis()->SetTitle("Frequency");
  fListMeasureCF->Add(fNpbnFluc); // -->:4
  
  TH2F *fNFluc = new TH2F("hNFluc", "N [(+/-) - (-/+)] Fuctuation",kCentBin,0, (Double_t)kCentBin,100,-5,5);
  fNFluc->GetXaxis()->SetTitle(" #-");
  fNFluc->GetYaxis()->SetTitle("Frequency");
  fListMeasureCF->Add(fNFluc); // -->:5
  
}


//____________________________________________________________________//
void AliEbyEChargeFluctuationAnalysis::Analyze(AliESDEvent* esd) {

  // Analysis from ESD
  ((TH1F*)(fListCFQA->At(0)))->Fill(1.);  

  Int_t kCent = 1;

  TObjArray *array = new TObjArray();
  Int_t itr = esd->GetNumberOfTracks();
  for (Int_t iTracks = 0; iTracks < itr; iTracks++) 
    {
      AliESDtrack* track = esd->GetTrack(iTracks);
      if (!track) continue;
      array->Add(track);
    }

  Calculate(array,kCent);

  delete array;
}


//____________________________________________________________________//
void AliEbyEChargeFluctuationAnalysis::Analyze(AliAODEvent* aod) {
  // Analysis from AOD
  ((TH1F*)(fListCFQA->At(0)))->Fill(1.);

  printf("%d\n",aod->GetNTracks());
  

}


//____________________________________________________________________//
void AliEbyEChargeFluctuationAnalysis::Analyze(AliStack* stack) {
  // Analysis from MC stack
  ((TH1F*)(fListCFQA->At(0)))->Fill(1.);  

  printf("%d \n",stack->GetNtrack());

}

//____________________________________________________________________________//
void AliEbyEChargeFluctuationAnalysis::Calculate(TObjArray *gTrackArray, Int_t cent){

  Int_t i = 0;
  AliVParticle* track = 0;
  Int_t pCharge = 0;
  Int_t nCharge = 0;
  Int_t gNtrack = gTrackArray->GetEntries();
  TString gAnalysisLevel = fEbyEBase->GetAnalysisLevel();  
  Double_t par[6] = {0.0};

  for(i = 0; i < gNtrack; i++) {
    if(gAnalysisLevel == "ESD")
      track = dynamic_cast<AliESDtrack*>(gTrackArray->At(i));
    else if(gAnalysisLevel == "AOD")
     track = dynamic_cast<AliAODTrack *>(gTrackArray->At(i));
    else if(gAnalysisLevel == "MC")
      track = dynamic_cast<AliMCParticle *>(gTrackArray->At(i));
    Short_t charge = track->Charge();
    if(charge > 0) pCharge += 1;
    if(charge < 0) nCharge += 1;
  }
  
  par[0] = (Double_t)pCharge;
  par[1] = (Double_t)nCharge;
  par[2] = (Double_t)pCharge + (Double_t)nCharge;
  
  if(pCharge != 0 ) par[3]   = (Double_t)pCharge/(Double_t)nCharge;
  if(nCharge != 0 ) par[4]  = (Double_t)nCharge/(Double_t)pCharge;
  
  par[5] = (Double_t)pCharge - (Double_t)nCharge;
    
  ((TH2F*)(fListMeasureCF->At(0)))->Fill(cent, par[2]); 
  ((TH2F*)(fListMeasureCF->At(1)))->Fill(cent, par[0]); 
  ((TH2F*)(fListMeasureCF->At(2)))->Fill(cent, par[1]); 
  ((TH2F*)(fListMeasureCF->At(3)))->Fill(cent, par[4]); 
  ((TH2F*)(fListMeasureCF->At(4)))->Fill(cent, par[3]); 
  ((TH2F*)(fListMeasureCF->At(5)))->Fill(cent, par[5]); 

}
