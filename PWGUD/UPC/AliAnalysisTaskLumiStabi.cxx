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
#include "TH1I.h"
#include "TTree.h"
#include "TFile.h"

// aliroot headers
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliOADBContainer.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"

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
    fCentrality(0),
    fCentralityPercentile(300),
    fTrgClassCINTZAC(0),
    fTrgInputV0M(0),
    hCentralityV0M(0),
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
    fCentrality(0),
    fCentralityPercentile(300),
    fTrgClassCINTZAC(0),
    fTrgInputV0M(0),
    hCentralityV0M(0),
    hTriggerInputsCounter(0),
    hTriggerClassesCounter(0)
{

//  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());

}//AliAnalysisTaskLumiStabi

//_____________________________________________________________________________
AliAnalysisTaskLumiStabi::~AliAnalysisTaskLumiStabi()
{
  // Destructor
  
  // Destructor
  if (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() != AliAnalysisManager::kProofAnalysis){
     delete fOutputList;
     fOutputList = 0x0;
  }

}//~AliAnalysisTaskLumiStabi


//_____________________________________________________________________________
void AliAnalysisTaskLumiStabi::UserCreateOutputObjects()
{
  
//  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
//  AliInputEventHandler *inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());

  fOutputList = new TList();
  fOutputList ->SetOwner();

  tOutput = new TTree("tOutput", "tOutput");
  tOutput ->Branch("fRunNumber", &fRunNumber, "fRunNumber/I");
  tOutput ->Branch("fTrgClassCINTZAC", &fTrgClassCINTZAC, "fTrgClassCINTZAC/O");
  tOutput ->Branch("fTrgInputV0M", &fTrgInputV0M, "fTrgInputV0M/O");
  tOutput ->Branch("fCentralityPercentile", &fCentralityPercentile);
  tOutput ->Branch("fCentrality", &fCentrality);
  fOutputList->Add(tOutput);

  hTriggerClassesCounter = new TH1I("hTriggerClassesCounter","Number of analyzed triggers per run",ENDRUN-STARTRUN,STARTRUN,ENDRUN);
  fOutputList->Add(hTriggerClassesCounter);

  hTriggerInputsCounter = new TH1I("hTriggerInputsCounter","Number of analyzed trigger inputs per run",ENDRUN-STARTRUN,STARTRUN,ENDRUN);
  fOutputList->Add(hTriggerInputsCounter);

  hCentralityV0M = new TH1D("hCentralityV0M",";V0M centrality [%];Counts per mil [-]",1000,0,100);
  fOutputList->Add(hCentralityV0M);

  PostData(1, fOutputList);

}//UserCreateOutputObjects


//_____________________________________________________________________________
void AliAnalysisTaskLumiStabi::UserExec(Option_t *)
{

  AliAODEvent *event = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!event) return;
  Printf("Event %s was loaded",event->GetName());

  fRunNumber = event->GetRunNumber();
  fTrgClassCINTZAC = event->GetFiredTriggerClasses().Contains("CINT7ZAC-B-NOPF-CENT");
//  fTrgInputV0M = event->GetHeader()->IsTriggerInputFired("0V0M");
  UInt_t fL0inputs = event->GetHeader()->GetL0TriggerInputs();
//  Int_t inputV0M = 7; //V0M in Pb-Pb
  Int_t inputV0M = 13; //V0M in Xe-Xe
  fTrgInputV0M =  fL0inputs & (1 << (inputV0M-1));

  fCentrality = (AliMultSelection * ) event->FindListObject("MultSelection");
  if( !fCentrality) {
     //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
     Printf("AliMultSelection object not found!");
  }else{
    fCentralityPercentile = fCentrality->GetMultiplicityPercentile("V0M");
  }

  if (fTrgClassCINTZAC) hTriggerClassesCounter->Fill(fRunNumber);
  if (fTrgInputV0M) hTriggerInputsCounter->Fill(fRunNumber);
  if (fTrgClassCINTZAC && fTrgInputV0M) hCentralityV0M->Fill(fCentralityPercentile);

  Printf("fCentrality->GetCentralityPercentile(V0M): %.f",fCentralityPercentile);
  Printf("fTrgClassCINTZAC %i",fTrgClassCINTZAC);
  Printf("fTrgInputV0M %i",fTrgInputV0M);

  tOutput->Fill();

  PostData(1, fOutputList);

}//UserExec


//_____________________________________________________________________________
void AliAnalysisTaskLumiStabi::Terminate(Option_t *)
{
  cout<<"Analysis complete."<<endl;
}//Terminate

