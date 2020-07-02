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
#include "AliAnalysisTaskSatellites.h"

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
#include "AliESDZDC.h"

// my headers

ClassImp(AliAnalysisTaskSatellites);

using std::cout;
using std::endl;


//_____________________________________________________________________________
AliAnalysisTaskSatellites::AliAnalysisTaskSatellites()
  : AliAnalysisTaskSE(),
    fOutputList(0),
    tOutput(0),
    fRunNumber(0),
    fL0inputs(0),
    fIsSatellite(0),
    fTrgClassCINTZAC(0),
    fTrgInputV0M(0),
    hDummyCounter(0),
    hSatellitesCounter(0),
    hTriggerInputsCounter(0),
    hTriggerClassesCounter(0)
{

}//AliAnalysisTaskSatellites


//_____________________________________________________________________________
AliAnalysisTaskSatellites::AliAnalysisTaskSatellites(const char *name)
  : AliAnalysisTaskSE(name),
    fOutputList(0),
    tOutput(0),
    fRunNumber(0),
    fL0inputs(0),
    fIsSatellite(0),
    fTrgClassCINTZAC(0),
    fTrgInputV0M(0),
    hDummyCounter(0),
    hSatellitesCounter(0),
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

}//AliAnalysisTaskSatellites

//_____________________________________________________________________________
AliAnalysisTaskSatellites::~AliAnalysisTaskSatellites()
{
  // Destructor
  if (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() != AliAnalysisManager::kProofAnalysis){
     delete fOutputList;     fOutputList = 0x0;
     delete tOutput;     tOutput = 0x0;
  }

}//~AliAnalysisTaskSatellites


//_____________________________________________________________________________
void AliAnalysisTaskSatellites::UserCreateOutputObjects()
{

  fOutputList = new TList();
  fOutputList ->SetOwner();

  tOutput = new TTree("tOutput", "tOutput");
  tOutput ->Branch("fRunNumber", &fRunNumber);
  tOutput ->Branch("fL0inputs",&fL0inputs);
  tOutput ->Branch("fIsSatellite", &fIsSatellite);
  tOutput ->Branch("fTrgClassCINTZAC", &fTrgClassCINTZAC);
  tOutput ->Branch("fTrgInputV0M", &fTrgInputV0M);
  tOutput ->Branch("fZNATDCm", &fZNATDCm,"fZNATDCm[4]/F");
  tOutput ->Branch("fZNCTDCm", &fZNCTDCm,"fZNCTDCm[4]/F");
  
  const Int_t STARTRUN = 240000;
  const Int_t ENDRUN = 300000;

  hDummyCounter = new TH1I("hDummyCounter","Number of events per run",ENDRUN-STARTRUN,STARTRUN,ENDRUN);
  fOutputList->Add(hDummyCounter);

  hTriggerClassesCounter = new TH1I("hTriggerClassesCounter","Number of analyzed triggers per run",ENDRUN-STARTRUN,STARTRUN,ENDRUN);
  fOutputList->Add(hTriggerClassesCounter);

  hTriggerInputsCounter = new TH1I("hTriggerInputsCounter","Number of analyzed trigger inputs per run",ENDRUN-STARTRUN,STARTRUN,ENDRUN);
  fOutputList->Add(hTriggerInputsCounter);

  hSatellitesCounter = new TH1I("hSatellitesCounter","Number of satellites in V0M per run",ENDRUN-STARTRUN,STARTRUN,ENDRUN);
  fOutputList->Add(hSatellitesCounter);

  PostData(1, fOutputList);
  PostData(2, tOutput);

}//UserCreateOutputObjects

//_____________________________________________________________________________
void AliAnalysisTaskSatellites::UserExec(Option_t *)
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
//  Printf("Event %s was loaded",event->GetName());

  hDummyCounter->Fill(fRunNumber);  // simple counter for basic information

  // ZDC timing decision
  AliESDZDC *ZDCdata = (AliESDZDC*) event->GetZDCData();
  Int_t detChZNA  = ZDCdata->GetZNATDCChannel();
  Int_t detChZNC  = ZDCdata->GetZNCTDCChannel();
  if (event->GetRunNumber()>=245726 && event->GetRunNumber()<=245793) detChZNA = 10;
  for (Int_t i=0;i<4;i++){
    fZNATDCm[i] = ZDCdata->GetZDCTDCCorrected(detChZNA,i);
    fZNCTDCm[i] = ZDCdata->GetZDCTDCCorrected(detChZNC,i);
  }
  fIsSatellite = IsSatellite(ZDCdata);

  //Pick only trigger
//  fTrgClassCINTZAC = event->GetFiredTriggerClasses().Contains("CINT7ZAC-B-NOPF-CENT");
//  if (!fTrgClassCINTZAC) return;

  fRunNumber = event->GetRunNumber();
  fL0inputs = event->GetHeader()->GetL0TriggerInputs();
  Int_t inputV0M = 7; //V0M in Pb-Pb
  if (fRunNumber == 280234 || fRunNumber == 280235) inputV0M = 13; //V0M in Xe-Xe
  fTrgInputV0M =  fL0inputs & (1 << (inputV0M-1));
  if (!fTrgInputV0M) return;

//  if (fTrgClassCINTZAC) hTriggerClassesCounter->Fill(fRunNumber);
  if (fTrgInputV0M) hTriggerInputsCounter->Fill(fRunNumber);
  if (fTrgInputV0M && fIsSatellite) hSatellitesCounter->Fill(fRunNumber);

  Printf("fIsSatellite %i",fIsSatellite);
//  Printf("fTrgClassCINTZAC %i",fTrgClassCINTZAC);
  Printf("fTrgInputV0M %i",fTrgInputV0M);

  tOutput->Fill();

  PostData(1, fOutputList);
  PostData(2, tOutput);

}//UserExec

Bool_t AliAnalysisTaskSatellites::IsSatellite(AliESDZDC *data)
{
  Int_t detChZNA  = data->GetZNATDCChannel();
  Int_t detChZNC  = data->GetZNCTDCChannel();
//  if (event->GetRunNumber()>=245726 && event->GetRunNumber()<=245793) detChZNA = 10;
  for (Int_t i=0;i<4;i++){
    fZNATDCm[i] = data->GetZDCTDCCorrected(detChZNA,i);
    fZNCTDCm[i] = data->GetZDCTDCCorrected(detChZNC,i);
  }
  for (Int_t i = 0; i < 4; i++){
    if (TMath::Abs(data->GetZDCTDCCorrected(detChZNA,i))>2.5 && TMath::Abs(data->GetZDCTDCCorrected(detChZNA,i))<25) return kTRUE;
    if (TMath::Abs(data->GetZDCTDCCorrected(detChZNC,i))>2.5 && TMath::Abs(data->GetZDCTDCCorrected(detChZNC,i))<25) return kTRUE;
  }
  return kFALSE;
}//IsSatellite

//_____________________________________________________________________________
void AliAnalysisTaskSatellites::Terminate(Option_t *)
{
  cout<<"Analysis complete."<<endl;
}//Terminate

