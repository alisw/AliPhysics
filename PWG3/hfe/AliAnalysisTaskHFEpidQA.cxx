/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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
//
// Task for PID QA
// Using AliHFEpidQA and AliHFEMCpidQA
// 
// Authors
//   Matus Kalisky <matus.kalisky@cern.ch>
//   Markus Heide <mheide@uni-muenster.de>
//   Markus Fasel <M.Fasel@gsi.de>
//
#include <TH1I.h>
#include <TList.h>
#include <TFile.h>

#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliHFEpidQA.h"
#include "AliHFEtools.h"
#include "AliESDInputHandler.h"

#include "AliHFEtrdPIDqa.h"

#include "AliAnalysisTaskHFEpidQA.h"

ClassImp(AliAnalysisTaskHFEpidQA)

AliAnalysisTaskHFEpidQA::AliAnalysisTaskHFEpidQA():
  AliAnalysisTaskSE("pidQAtask")
  , fPIDqa(NULL)
  , fOutput(NULL)
  , fEvents(NULL)
  , fNNref(NULL)
  , fTRDTotalChargeInSlice0(kFALSE)
{
  //
  // Default Constructor
  //
}

AliAnalysisTaskHFEpidQA::AliAnalysisTaskHFEpidQA(const Char_t *name):
  AliAnalysisTaskSE(name)
  , fPIDqa(NULL)
  , fOutput(NULL)
  , fEvents(NULL)
  , fNNref(NULL)
  , fTRDTotalChargeInSlice0(kFALSE)
{
  //
  // Default Constructor
  //
  DefineOutput(1, TList::Class());

}

AliAnalysisTaskHFEpidQA::~AliAnalysisTaskHFEpidQA(){
  //
  // Destructor
  //
  if(fPIDqa) delete fPIDqa;
  if(fOutput) delete fOutput;
}

void AliAnalysisTaskHFEpidQA::UserCreateOutputObjects(){
  //
  // Create the output
  // Initialize PID QA
  //
  fOutput = new TList;
  fOutput->SetOwner();

  // Counter for number of events
  fOutput->Add((fEvents = new TH1I("nEvents", "NumberOfEvents", 1, 1, 2)));

  fPIDqa = new AliHFEpidQA;
  if(fTRDTotalChargeInSlice0) fPIDqa->SetTRDTotalChargeInSlice0();
  if(HasV0pidQA()) fPIDqa->SetV0pidQA();
  if(HasRecalculateTRDpid()) fPIDqa->SetRecalculateTRDpid();
  if(fNNref) fPIDqa->SetNNref(fNNref);
  fPIDqa->Init();

  TList *tmp = fPIDqa->GetOutput();
  tmp->SetName("PIDqa");
  fOutput->Add(tmp);
  if(HasV0pidQA()){
    tmp = fPIDqa->GetV0pidQA();
    tmp->SetName("V0pidQA");
    fOutput->Add(tmp);
  }
  tmp = 0x0;
  tmp = fPIDqa->GetV0pidMC();
  if(tmp){
    tmp->SetName("V0pidMC");
    fOutput->Add(tmp);
  }

  // Add TRD PID QA object to the output
  fOutput->Add(fPIDqa->GetTRDQA());
}

Bool_t AliAnalysisTaskHFEpidQA::UserNotify(){
  // DEBUG
  //printf("*****\n");
  //printf(" -D Current File Name: %s \n", CurrentFileName());
  return AliAnalysisTask::Notify();
}

void AliAnalysisTaskHFEpidQA::UserExec(Option_t *){
  //
  // Event Loop
  // 
  AliMCEventHandler* mcHandler = (dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()));
  AliESDInputHandler *inh = dynamic_cast<AliESDInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  AliPIDResponse *workingPID = NULL;
  if(inh && (workingPID = inh->GetPIDResponse()))
    fPIDqa->SetPIDResponse(workingPID);
  else fPIDqa->SetPIDResponse(AliHFEtools::GetDefaultPID(mcHandler ? kTRUE : kFALSE, kFALSE));
  
  // check the MC data
  if(fMCEvent && !mcHandler ) return;
  if(fMCEvent &&  !mcHandler->InitOk() ) return;
  if(fMCEvent &&  !mcHandler->TreeK() ) return;
  if(fMCEvent &&  !mcHandler->TreeTR() ) return;
  if(fMCEvent) fPIDqa->SetMCEvent(fMCEvent);
  
  fPIDqa->SetEvent(fInputEvent);
  fPIDqa->Process();
  fEvents->Fill(1.1);
  PostData(1, fOutput);
}

void AliAnalysisTaskHFEpidQA::Terminate(Option_t *){
  //
  // Do Post Processing
  //
}


