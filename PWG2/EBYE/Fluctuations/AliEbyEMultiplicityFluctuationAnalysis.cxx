
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
 *             AliEbyEMultiplicityFluctuationAnalysis base class
 *           This class deals with the Multiplicity   fluctuation 
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

#include "AliESDtrackCuts.h"
#include "AliEbyEEventBase.h"
#include "AliEbyEMultiplicityFluctuationAnalysis.h"

ClassImp(AliEbyEMultiplicityFluctuationAnalysis)

//____________________________________________________________________//
AliEbyEMultiplicityFluctuationAnalysis::AliEbyEMultiplicityFluctuationAnalysis() : 
  TObject(), 
  fListMFQA(0),
  fListMeasureMF(0),
  fEbyEBase(0),
  fESDtrackCuts(0),
  fIsTreeMode(kFALSE),
  iNpossitive(0),
  iNnegative(0),
  iCentrality(0)
{
  if(fIsTreeMode) fFluctuationTree = NULL;
  InitHistos();
  //Default constructor
}


//____________________________________________________________________//
AliEbyEMultiplicityFluctuationAnalysis::~AliEbyEMultiplicityFluctuationAnalysis() {
  //Default destructor
 
  if(fListMFQA) delete fListMFQA;
  if(fListMeasureMF) delete fListMeasureMF;
  if(fEbyEBase) delete fEbyEBase;
  if(fESDtrackCuts) delete fESDtrackCuts;
  if(fFluctuationTree) delete fFluctuationTree;
  /*  if(fIsTreeMode) delete fIsTreeMode;
  if(iNpossitive) delete iNpossitive;
  if(iNnegative) delete iNnegative;
  if(iCentrality) delete iCentrality;
  */
}

//____________________________________________________________________//
void AliEbyEMultiplicityFluctuationAnalysis::InitHistos() {

  if(fIsTreeMode) { 
    fFluctuationTree = new TTree("EbyET","FluctaionContainer");
    fFluctuationTree->Branch("Positive",&iNpossitive,"iNpossitive/I");
    fFluctuationTree->Branch("Negative",&iNnegative,"iNnegative/I");
    fFluctuationTree->Branch("Centrality",&iCentrality,"iCentrality/I");
  }

  fListMFQA = new TList();
  fListMFQA->SetName("MFQaList");

  fListMeasureMF = new TList();
  fListMeasureMF->SetName("MFMeasureList");

  // Int_t kCentBin = fEbyEBase->GetCentralityBin();
  Int_t kCentBin = 50;  
  Int_t kMultBin = 6000;
  Double_t kMultMin = 0.0;
  Double_t kMultMax = 6000;  
  
  
  TH1F *fEvtStat = new TH1F("hEventCounterMH"," Event Counters for MF Analysis", 20,0,20);
  fListMFQA->Add(fEvtStat); // --:0

  TH2F *fhCentralityV0 = new TH2F("hCentralityV0", "Total Charge Fuctuation", kCentBin,0, (Double_t)kCentBin, 3000,0,30000);
  fListMFQA->Add(fhCentralityV0); // --:1

  TH1F *fhCentDist= new TH1F("hCentDist", "Total Charge Fuctuation", kCentBin,0, (Double_t)kCentBin);
  fListMFQA->Add(fhCentDist); // --:2


  TH1F *fEvtStat1 = new TH1F("hEventNp"," Event Counters for +", 60,0,60);
  fListMFQA->Add(fEvtStat1); // --:3
  
  TH1F *fEvtStat2 = new TH1F("hEventNn"," Event Counters for +", 60,0,60);
  fListMFQA->Add(fEvtStat2); // --:4

  TH1F *fEvtStat3 = new TH1F("hEventNppNn"," Event Counters for + + -", 60,0,60);
  fListMFQA->Add(fEvtStat3); // --:5

  
  
  TH2F *fNpCharge = new TH2F("hNpCharge", "Possitive Charge ",kCentBin,0, (Double_t)kCentBin, kMultBin,kMultMin, kMultMax);
  fListMeasureMF->Add(fNpCharge); // -->:0
  
  TH2F *fNnCharge = new TH2F("hNnCharge", "Negative Charge ", kCentBin,0, (Double_t)kCentBin, kMultBin,kMultMin, kMultMax);
  fListMeasureMF->Add(fNnCharge); // -->:1
  
  TH2F *fNCharge = new TH2F("hNCharge", "Total Charge [+ + -]",kCentBin,0, (Double_t)kCentBin, 2*kMultBin,kMultMin, 2*kMultMax);
  fListMeasureMF->Add(fNCharge); // -->:2

  TH2F *fNetCharge = new TH2F("hNetCharge", "Net Charge [+ - -]",kCentBin,0, (Double_t)kCentBin, kMultBin,kMultMin, kMultMax);
  fListMeasureMF->Add(fNetCharge); // -->:3

  TH2F *fNnbp = new TH2F("hNnbp", "[-/+] ",kCentBin,0, (Double_t)kCentBin, 1000,0,20); 
  fListMeasureMF->Add(fNnbp); // -->:4
  
  TH2F *fNpbn = new TH2F("hNpbn", "[+/-] ",kCentBin,0, (Double_t)kCentBin, 1000,0,20);
  fListMeasureMF->Add(fNpbn); // -->:5
  
  TH2F *fNpbyNnmNnbyNp = new TH2F("hNpbyNnmNnbyNp", "[(+/-) - (-/+)] ",kCentBin,0, (Double_t)kCentBin,100,-5,5);
  fListMeasureMF->Add(fNpbyNnmNnbyNp); // -->:6

  TH2F *fNpStarNn = new TH2F("hNpbyNnmNnbyNp", "[(+)(-)] ",kCentBin,0, (Double_t)kCentBin,kMultBin,kMultMin, kMultMax*kMultMax);
  fListMeasureMF->Add(fNpStarNn); // -->:7

  TH2F *fCorrPvsN = new TH2F("hCorrPvsN", " + ~ -",200,0,6000,200,0,6000);
  fListMeasureMF->Add(fCorrPvsN); // -->:8

  TH3F *fCorCent = new TH3F("hCorrCentrality", "Number of +ve ~ Number of -ve",kCentBin,0, (Double_t)kCentBin,200,0,6000,200,0,6000);
  fListMeasureMF->Add(fCorCent); // -->:9

  
  
}


//____________________________________________________________________//
void AliEbyEMultiplicityFluctuationAnalysis::Analyze(AliESDEvent* esd) {

  // Analysis from ESD
  ((TH1F*)(fListMFQA->At(0)))->Fill(1.);  

  Int_t kCent = 1;
   kCent =  fEbyEBase->GetCentrality(esd);
 
  ((TH1F*)(fListMFQA->At(2)))->Fill(kCent);  
  
  TObjArray *array = new TObjArray();
  Int_t itr = esd->GetNumberOfTracks();
  
  for (Int_t iTracks = 0; iTracks < itr; iTracks++) 
    {
      AliESDtrack* track = esd->GetTrack(iTracks);
      if (!track) continue;
      if(!fESDtrackCuts->AcceptTrack(track)) continue;
      array->Add(track);
    }

  Calculate(array,kCent);

  delete array;
}


//____________________________________________________________________//
void AliEbyEMultiplicityFluctuationAnalysis::Analyze(AliAODEvent* aod) {
  // Analysis from AOD
  ((TH1F*)(fListMFQA->At(0)))->Fill(1.);

  printf("%d\n",aod->GetNTracks());
  

}


//____________________________________________________________________//
void AliEbyEMultiplicityFluctuationAnalysis::Analyze(AliStack* stack) {
  // Analysis from MC stack
  ((TH1F*)(fListMFQA->At(0)))->Fill(1.);  

  printf("%d \n",stack->GetNtrack());

}

//____________________________________________________________________________//
void AliEbyEMultiplicityFluctuationAnalysis::Calculate(TObjArray *gTrackArray, Int_t cent){

  ((TH1F*)(fListMFQA->At(0)))->Fill(1,1); 
  if(cent < 0) return;
  ((TH1F*)(fListMFQA->At(0)))->Fill(10,1); 

  Int_t gNtrack = gTrackArray->GetEntries();
  if(gNtrack < 1) return;

  ((TH1F*)(fListMFQA->At(0)))->Fill(20,1); 

  iNpossitive = 0;
  iNnegative = 0;
  iCentrality = 0;

  AliVParticle* track = 0;
  Int_t pCharge = 0;
  Int_t nCharge = 0;


  TString gAnalysisLevel = fEbyEBase->GetAnalysisLevel();  
  Double_t par[8];
  for(Int_t i = 0; i < 8; i++) par[i] = 0.0;
  Double_t pr[8];
  for(Int_t i = 0; i < 8; i++) pr[i] = 0.0;

   

  for(Int_t i = 0; i < gNtrack; i++) {
    if(gAnalysisLevel == "ESD")
      track = dynamic_cast<AliESDtrack*>(gTrackArray->At(i));
    Short_t charge = track->Charge();
    if(charge > 0) {
      pCharge += 1; 
      ((TH1F*)(fListMFQA->At(3)))->Fill(cent,1);  
      ((TH1F*)(fListMFQA->At(5)))->Fill(cent,1);  
    }
    if(charge < 0) {
      nCharge += 1; 
      ((TH1F*)(fListMFQA->At(4)))->Fill(cent,1); 
      ((TH1F*)(fListMFQA->At(5)))->Fill(cent,1);   
    }
  }
  
  ((TH1F*)(fListMFQA->At(0)))->Fill(131+cent,pCharge);  
  ((TH1F*)(fListMFQA->At(0)))->Fill(5,nCharge);  
  ((TH1F*)(fListMFQA->At(0)))->Fill(6,pCharge+nCharge);  
  ((TH2F*)(fListMFQA->At(1)))->Fill(cent,pCharge+nCharge);  


  par[0] = (Double_t)pCharge;
  par[1] = (Double_t)nCharge;
  par[2] = (Double_t)pCharge + (Double_t)nCharge;
  par[3] = (Double_t)pCharge - (Double_t)nCharge;
  
  if(pCharge != 0 ) par[4]   = (Double_t)pCharge/(Double_t)nCharge;
  if(nCharge != 0 ) par[5]  = (Double_t)nCharge/(Double_t)pCharge;
  
  par[6] = (Double_t)pCharge - (Double_t)nCharge;
  par[7] = (Double_t)pCharge*(Double_t)nCharge;

  //  printf("Centrality: %2d  NTrack: %8d  +ve : %8d  -ve : %8d  %10.6f %10.6f %10.6f\n",
  //cent, gNtrack, pCharge, nCharge, pr[0], pr[1],value);


  ((TH2F*)(fListMeasureMF->At(0)))->Fill(cent, par[0]); //+
  ((TH2F*)(fListMeasureMF->At(1)))->Fill(cent, par[1]); //-
  ((TH2F*)(fListMeasureMF->At(2)))->Fill(cent, par[2]); //++-
  ((TH2F*)(fListMeasureMF->At(3)))->Fill(cent, par[3]); //+--
  ((TH2F*)(fListMeasureMF->At(4)))->Fill(cent, par[4]); //+/-
  ((TH2F*)(fListMeasureMF->At(5)))->Fill(cent, par[5]); //-/+
  ((TH2F*)(fListMeasureMF->At(6)))->Fill(cent, par[6]); //+/-
  ((TH2F*)(fListMeasureMF->At(7)))->Fill(cent, par[7]);
  ((TH2F*)(fListMeasureMF->At(8)))->Fill(pCharge,nCharge); 
  ((TH3F*)(fListMeasureMF->At(9)))->Fill(cent, pCharge,nCharge); 
  
    if(fIsTreeMode) {
      iNpossitive = pCharge;
      iNnegative = nCharge;
      iCentrality = cent;
      //   fFluctuationTree->Fill();
    
    }
  
}
