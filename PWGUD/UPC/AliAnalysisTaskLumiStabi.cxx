/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
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

// c++ headers
#include <iostream>
#include <string.h>

// root headers
#include "TROOT.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TTree.h"
#include "TFile.h"
#include "TList.h"

// aliroot headers
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliVEvent.h"
#include "AliVCuts.h"
#include "AliMultSelection.h"
#include "AliAODZDC.h"

// my headers
#include "AliAnalysisTaskLumiStabi.h"

ClassImp(AliAnalysisTaskLumiStabi);

using std::cout;
using std::endl;


//_____________________________________________________________________________
AliAnalysisTaskLumiStabi::AliAnalysisTaskLumiStabi()
  : AliAnalysisTaskSE(),
    fOutputList(0),
    tOutput(0),
    fRunNumber(0),
    fL0inputs(0),
    fSelectPhysics(0),
    fIsSatellite(0),
//    fCentrality(0),
    fV0McentPercentile(300),
    fTrgClassCINTZAC(0),
    fTrgInputV0M(0),
    hCentralityV0M(0),
    hCentralityV0MandPS(0),
    hCentralityV0MandSat(0),
    hDummyCounter(0),
    hTriggerInputsCounter(0),
    hTriggerClassesCounter(0)
{

}//AliAnalysisTaskLumiStabi


//_____________________________________________________________________________
AliAnalysisTaskLumiStabi::AliAnalysisTaskLumiStabi(const char *name)
  : AliAnalysisTaskSE(name),
    fOutputList(0),
    tOutput(0),
    fRunNumber(0),
    fL0inputs(0),
    fSelectPhysics(0),
    fIsSatellite(0),
//    fCentrality(0),
    fV0McentPercentile(300),
    fTrgClassCINTZAC(0),
    fTrgInputV0M(0),
    hCentralityV0M(0),
    hCentralityV0MandPS(0),
    hCentralityV0MandSat(0),
    hDummyCounter(0),
    hTriggerInputsCounter(0),
    hTriggerClassesCounter(0)
{

  for (Int_t i=0;i<4;i++){
    fZNATDCm[i] = 0.;
    fZNCTDCm[i] = 0.;
  }

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());

}//AliAnalysisTaskLumiStabi

//_____________________________________________________________________________
AliAnalysisTaskLumiStabi::~AliAnalysisTaskLumiStabi()
{
  // Destructor
  if (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() != AliAnalysisManager::kProofAnalysis){
     delete fOutputList;     fOutputList = 0x0;
     delete tOutput;     tOutput = 0x0;
  }

}//~AliAnalysisTaskLumiStabi


//_____________________________________________________________________________
void AliAnalysisTaskLumiStabi::UserCreateOutputObjects()
{

  fOutputList = new TList();
  fOutputList ->SetOwner();

  tOutput = new TTree("tOutput", "tOutput");
  tOutput ->Branch("fRunNumber", &fRunNumber);
  tOutput ->Branch("fL0inputs",&fL0inputs);
  tOutput ->Branch("fIsSatellite", &fIsSatellite);
  tOutput ->Branch("fSelectPhysics", &fSelectPhysics);
  tOutput ->Branch("fTrgClassCINTZAC", &fTrgClassCINTZAC);
  tOutput ->Branch("fTrgInputV0M", &fTrgInputV0M);
  tOutput ->Branch("fV0McentPercentile", &fV0McentPercentile);
  tOutput ->Branch("fZNATDCm", &fZNATDCm,"fZNATDCm[4]/F");
  tOutput ->Branch("fZNCTDCm", &fZNCTDCm,"fZNCTDCm[4]/F");
//  tOutput ->Branch("fCentrality", &fCentrality);
//  fOutputList->Add(tOutput);

  hDummyCounter = new TH1I("hDummyCounter","Number of events per run",ENDRUN-STARTRUN,STARTRUN,ENDRUN);
  fOutputList->Add(hDummyCounter);

  hTriggerClassesCounter = new TH1I("hTriggerClassesCounter","Number of analyzed triggers per run",ENDRUN-STARTRUN,STARTRUN,ENDRUN);
  fOutputList->Add(hTriggerClassesCounter);

  hTriggerInputsCounter = new TH1I("hTriggerInputsCounter","Number of analyzed trigger inputs per run",ENDRUN-STARTRUN,STARTRUN,ENDRUN);
  fOutputList->Add(hTriggerInputsCounter);

  hCentralityV0M = new TH1D("hCentralityV0M",";V0M centrality [%];Counts per mil [-]",1000,0,100);
  fOutputList->Add(hCentralityV0M);

  hCentralityV0MandPS = new TH1D("hCentralityV0MandPS",";V0M centrality [%];Counts per mil [-]",1000,0,100);
  fOutputList->Add(hCentralityV0MandPS);

  hCentralityV0MandSat = new TH1D("hCentralityV0MandSat",";V0M centrality [%];Counts per mil [-]",1000,0,100);
  fOutputList->Add(hCentralityV0MandSat);

  PostData(1, fOutputList);
  PostData(2, tOutput);

}//UserCreateOutputObjects

//_____________________________________________________________________________
void AliAnalysisTaskLumiStabi::UserExec(Option_t *)
{
  //Get analysis manager
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  if( !man) {
     Printf("man object not found!");
     return;
  }

  //Get input handler
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  if( !inputHandler) {
     Printf("AliInputEventHandler object not found!");
     return;
  }

  //Get event
  AliVEvent *event = dynamic_cast<AliVEvent*>(InputEvent());
  if( !event) {
     Printf("AliVEvent object not found!");
     return;
  }
//  Printf("Event %s was loaded",event->GetName());\

  hDummyCounter->Fill(fRunNumber);  // simple counter for basic information

  // ZDC timing decision
  AliAODZDC *ZDCdata = (AliAODZDC*) event->GetZDCData();
  for (Int_t i=0;i<4;i++){
    fZNATDCm[i] = ZDCdata->GetZNATDCm(i);
    fZNCTDCm[i] = ZDCdata->GetZNCTDCm(i);
  }
  fIsSatellite = IsSatellite(ZDCdata);

  //Pick only trigger
  fTrgClassCINTZAC = event->GetFiredTriggerClasses().Contains("CINT7ZAC-B-NOPF-CENT");
  if (!fTrgClassCINTZAC) return;

  fRunNumber = event->GetRunNumber();
  fL0inputs = event->GetHeader()->GetL0TriggerInputs();
  Int_t inputV0M = 7; //V0M in Pb-Pb
  if (fRunNumber == 280234 || fRunNumber == 280235) inputV0M = 13; //V0M in Xe-Xe
//   fTrgClassCINTZAC = event->GetFiredTriggerClasses().Contains("CINT7ZAC-B-NOPF-CENT");
  fTrgInputV0M =  fL0inputs & (1 << (inputV0M-1));

  //Check physics selection
  fSelectPhysics = kTRUE;
  if(inputHandler->GetEventSelection()){
    fSelectPhysics = inputHandler->IsEventSelected();
  }
  else{
    Printf("inputHandler->GetEventSelection() object not found!");
  }

  AliMultSelection *centrality = (AliMultSelection * ) event->FindListObject("MultSelection");
  if(!centrality){
    Printf("AliMultSelection object not found!");
  }
  else{
    fV0McentPercentile = centrality->GetMultiplicityPercentile("V0M",kTRUE);
  }

  if (fTrgClassCINTZAC) hTriggerClassesCounter->Fill(fRunNumber);
  if (fTrgInputV0M) hTriggerInputsCounter->Fill(fRunNumber);
  if (fTrgClassCINTZAC && fTrgInputV0M) hCentralityV0M->Fill(fV0McentPercentile);
  if (fTrgClassCINTZAC && fTrgInputV0M && fSelectPhysics) hCentralityV0MandPS->Fill(fV0McentPercentile);
  if (fTrgClassCINTZAC && fTrgInputV0M && !fIsSatellite) hCentralityV0MandSat->Fill(fV0McentPercentile);

  Printf("centrality->GetCentralityPercentile(V0M): %.f",fV0McentPercentile);
  Printf("fIsSatellite %i",fIsSatellite);
  Printf("fSelectPhysics %i",fSelectPhysics);
  Printf("fTrgClassCINTZAC %i",fTrgClassCINTZAC);
  Printf("fTrgInputV0M %i",fTrgInputV0M);

  tOutput->Fill();

  PostData(1, fOutputList);
  PostData(2, tOutput);

}//UserExec

Bool_t AliAnalysisTaskLumiStabi::IsSatellite(AliAODZDC *data)
{
  for (Int_t i = 0; i < 4; i++){
    if (TMath::Abs(data->GetZNATDCm(i))>2.5 && TMath::Abs(data->GetZNATDCm(i))<25) return kTRUE;
    if (TMath::Abs(data->GetZNCTDCm(i))>2.5 && TMath::Abs(data->GetZNCTDCm(i))<25) return kTRUE;
  }
  return kFALSE;
}//IsSatellite

//_____________________________________________________________________________
void AliAnalysisTaskLumiStabi::Terminate(Option_t *)
{
  cout<<"Analysis complete."<<endl;
}//Terminate

