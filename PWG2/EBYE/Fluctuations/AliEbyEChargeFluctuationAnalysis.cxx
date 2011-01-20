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
  fEbyEBase(0),
  fhNpAverage(0),   
  fhNnAverage(0),   
  fhNCAverage(0),   
  fhNetAverage(0),  
  fhNnSNpAverage(0)
{
  //InitHistos();
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

void AliEbyEChargeFluctuationAnalysis::ReadFromFile() {
  //To read Average values from input histograms

  TFile *file = TFile::Open("InputFluctuationAnalysis.root");
  fhNpAverage    = (TH1F*)file->Get("hNpAverage");	
  fhNnAverage    = (TH1F*)file->Get("hNnAverage");	
  fhNCAverage    = (TH1F*)file->Get("hNCAverage");
  fhNetAverage   = (TH1F*)file->Get("hNetAverage");
  fhNnSNpAverage = (TH1F*)file->Get("hNnSNpAverage");	

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

  TH2F *fNChargeFluc = new TH2F("hNpCharge", "Possitive Charge",kCentBin,0, (Double_t)kCentBin, kMultBin,kMultMin, kMultMax);
  fNChargeFluc->GetXaxis()->SetTitle(" #-");
  fNChargeFluc->GetYaxis()->SetTitle("Frequency");
  fListMeasureCF->Add(fNChargeFluc); // -->:0

  TH2F *fNpChargeFluc = new TH2F("hNnCharge", "Negative Charge",kCentBin,0, (Double_t)kCentBin, kMultBin,kMultMin, kMultMax);
  fNpChargeFluc->GetXaxis()->SetTitle(" #-");
  fNpChargeFluc->GetYaxis()->SetTitle("Frequency");
  fListMeasureCF->Add(fNpChargeFluc); // -->:1
  
  TH2F *fNnChargeFluc = new TH2F("hNCharge", "[ (+) + (-)]", kCentBin,0, (Double_t)kCentBin, kMultBin,kMultMin, kMultMax);
  fNnChargeFluc->GetXaxis()->SetTitle(" #-");
  fNnChargeFluc->GetYaxis()->SetTitle("Frequency");
  fListMeasureCF->Add(fNpChargeFluc); // -->:2
  
  TH2F *fNnbpFluc = new TH2F("hNPmN", " [ (+) - (-)] ",kCentBin,0, (Double_t)kCentBin, 100,0,10); 
  fNnbpFluc->GetXaxis()->SetTitle(" #-");
  fNnbpFluc->GetYaxis()->SetTitle("Frequency");
  fListMeasureCF->Add(fNnbpFluc); // -->:3
  
  TH2F *fNpbnFluc = new TH2F("hNpByNpa", "[ +/<+> ] ",kCentBin,0, (Double_t)kCentBin, 100,0,10);
  fNpbnFluc->GetXaxis()->SetTitle(" #-");
  fNpbnFluc->GetYaxis()->SetTitle("Frequency");
  fListMeasureCF->Add(fNpbnFluc); // -->:4
  
  TH2F *fNFluc = new TH2F("hNnByNna", "[ -/<-> ]",kCentBin,0, (Double_t)kCentBin,100,-5,5);
  fNFluc->GetXaxis()->SetTitle(" #-");
  fNFluc->GetYaxis()->SetTitle("Frequency");
  fListMeasureCF->Add(fNFluc); // -->:5

  TH2F *fNFluc1 = new TH2F("hNpByPmNa", "[ +/<+ - -> ]",kCentBin,0, (Double_t)kCentBin,100,-5,5);
  fNFluc1->GetXaxis()->SetTitle(" #-");
  fNFluc1->GetYaxis()->SetTitle("Frequency");
  fListMeasureCF->Add(fNFluc1); // -->:6


  TH2F *fNFluc2 = new TH2F("hNnByPmNa", "[ -/<+ - -> ]",kCentBin,0, (Double_t)kCentBin,100,-5,5);
  fNFluc2->GetXaxis()->SetTitle(" #-");
  fNFluc2->GetYaxis()->SetTitle("Frequency");
  fListMeasureCF->Add(fNFluc2); // -->:7

  TH2F *fNFluc3 = new TH2F("hNpByPpNa", "[ +/<+ + -> ]",kCentBin,0, (Double_t)kCentBin,100,-5,5);
  fNFluc3->GetXaxis()->SetTitle(" #-");
  fNFluc3->GetYaxis()->SetTitle("Frequency");
  fListMeasureCF->Add(fNFluc3); // -->:8

  TH2F *fNFluc4 = new TH2F("hNnByPpNa", "[ -/<+ + -> ]",kCentBin,0, (Double_t)kCentBin,100,-5,5);
  fNFluc4->GetXaxis()->SetTitle(" #-");
  fNFluc4->GetYaxis()->SetTitle("Frequency");
  fListMeasureCF->Add(fNFluc4); // -->:9

  TH2F *fNFluc5 = new TH2F("hNpNnByNpNaa", "[ (+)(-)/<(-)(+)> ]",kCentBin,0, (Double_t)kCentBin,100,-5,5);
  fNFluc5->GetXaxis()->SetTitle(" #-");
  fNFluc5->GetYaxis()->SetTitle("Frequency");
  fListMeasureCF->Add(fNFluc5); // -->:10
  
}


//____________________________________________________________________//
void AliEbyEChargeFluctuationAnalysis::Analyze(AliESDEvent* esd) {

  // Analysis from ESD
  ((TH1F*)(fListCFQA->At(0)))->Fill(1.);  

  Int_t kCent =  fEbyEBase->GetCentrality(esd);
 
  
  TObjArray *array = new TObjArray();
  Int_t itr = esd->GetNumberOfTracks();

  // printf(" Centrality bin %d and Number of Good Track %d\n",itr,kCent);

  for (Int_t iTracks = 0; iTracks < itr; iTracks++) 
    {
      AliESDtrack* track = esd->GetTrack(iTracks);
      if (!track) continue;
      if(fEbyEBase->IsInPhaseSpace(track))
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


  ((TH1F*)(fListCFQA->At(0)))->Fill(1,1); 
  if(cent < 0) return;
  ((TH1F*)(fListCFQA->At(0)))->Fill(10,1); 

  Int_t gNtrack = gTrackArray->GetEntries();
  if(gNtrack < 1) return;

  ((TH1F*)(fListCFQA->At(0)))->Fill(20,1); 

  AliVParticle* track = 0;
  Int_t pCharge = 0;
  Int_t nCharge = 0;


  TString gAnalysisLevel = fEbyEBase->GetAnalysisLevel();  
  Double_t par[16];
  for(Int_t i = 0; i < 8; i++) par[i] = 0.0;
  Double_t pr[5];
  for(Int_t i = 0; i < 5; i++) pr[i] = 0.0;

   

  for(Int_t i = 0; i < gNtrack; i++) {
    if(gAnalysisLevel == "ESD")
      track = dynamic_cast<AliESDtrack*>(gTrackArray->At(i));
    Short_t charge = track->Charge();
    if(charge > 0) {
      pCharge += 1; 
      ((TH1F*)(fListCFQA->At(3)))->Fill(cent,1);  
      ((TH1F*)(fListCFQA->At(5)))->Fill(cent,1);  
    }
    if(charge < 0) {
      nCharge += 1; 
      ((TH1F*)(fListCFQA->At(4)))->Fill(cent,1); 
      ((TH1F*)(fListCFQA->At(5)))->Fill(cent,1);   
    }
  }
  
  ((TH1F*)(fListCFQA->At(0)))->Fill(131+cent,pCharge);  
  ((TH1F*)(fListCFQA->At(0)))->Fill(5,nCharge);  
  ((TH1F*)(fListCFQA->At(0)))->Fill(6,pCharge+nCharge);  
  ((TH2F*)(fListCFQA->At(1)))->Fill(cent,pCharge+nCharge);  

  pr[0] = fhNpAverage->GetBinContent(cent+1);   
  pr[1] = fhNnAverage->GetBinContent(cent+1);      
  pr[2] = fhNCAverage->GetBinContent(cent+1);      
  pr[3] = fhNetAverage->GetBinContent(cent+1);     
  pr[4] = fhNnSNpAverage->GetBinContent(cent+1);   


  par[0] = (Double_t)pCharge;
  par[1] = (Double_t)nCharge;
  par[2] = (Double_t)pCharge + (Double_t)nCharge;
  par[3] = (Double_t)pCharge - (Double_t)nCharge;
  par[4] = (Double_t)pCharge*((Double_t)pCharge - 1);
  par[5] = (Double_t)nCharge*((Double_t)nCharge - 1);
  par[6] = (Double_t)pCharge*(Double_t)nCharge;

  if(pr[0] != 0) par[7]  = par[0]/pr[0]; // +/<+>
  if(pr[1] != 0) par[8]  = par[0]/pr[1]; // -/<->
  if(pr[2] != 0) par[9]  = par[0]/pr[2]; // ++-/<++->
  if(pr[2] != 0) par[10] = par[0]/pr[2]; // +--/<++->
  if(pr[3] != 0) par[11] = par[0]/pr[3]; // ++-/<+-->
  if(pr[3] != 0) par[12] = par[0]/pr[3]; // +--/<+-->
  if(pr[4] != 0) par[13] = par[0]/pr[4]; //+*-/<+*->
  
  par[14] = (par[7] + par[8])*(par[7] + par[8]); // eqn - 1
  //  par[15] = (par[4]/pr[0]*[pr[0] + par[8])*(par[7] + par[8]); // eqn - 1

    
  Float_t value  = 0;
  if( pr[0] != 0 || pr[1] != 0)
    value = par[0]/pr[0] + par[1]/pr[1]; 

  printf("Centrality: %2d  NTrack: %8d  +ve : %8d  -ve : %8d  %10.6f %10.6f %10.6f\n",cent, gNtrack, pCharge, nCharge, pr[0], pr[1],value);


  ((TH2F*)(fListMeasureCF->At(0)))->Fill(cent, par[0]); //+
  ((TH2F*)(fListMeasureCF->At(1)))->Fill(cent, par[1]); //-
  ((TH2F*)(fListMeasureCF->At(2)))->Fill(cent, par[2]); //++-
  ((TH2F*)(fListMeasureCF->At(3)))->Fill(cent, par[3]); //+--
  ((TH2F*)(fListMeasureCF->At(4)))->Fill(cent, par[4]); //+/-
  ((TH2F*)(fListMeasureCF->At(5)))->Fill(cent, value); //-/+
  ((TH2F*)(fListMeasureCF->At(6)))->Fill(cent, value*value); //+/-
  ((TH2F*)(fListMeasureCF->At(7)))->Fill(cent, par[7]);
  ((TH2F*)(fListMeasureCF->At(8)))->Fill(pCharge,nCharge); 
  ((TH3F*)(fListMeasureCF->At(9)))->Fill(cent, pCharge,nCharge); 
  
}
