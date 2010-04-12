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

#include "AliAnalysisTaskHFEpidQA.h"
#include "AliHFEpidQA.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"

ClassImp(AliAnalysisTaskHFEpidQA)

AliAnalysisTaskHFEpidQA::AliAnalysisTaskHFEpidQA():
    AliAnalysisTaskSE("pidQAtask")
  , fPIDqa(NULL)
  , fOutput(NULL)
  , fEvents(NULL)
{
  //
  // Default Constructor
  //
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskHFEpidQA::AliAnalysisTaskHFEpidQA(const Char_t *name):
    AliAnalysisTaskSE(name)
  , fPIDqa(NULL)
  , fOutput(NULL)
  , fEvents(NULL)
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

  // Counter for number of events
  fOutput->Add((fEvents = new TH1I("nEvents", "NumberOfEvents", 1, 1, 2)));

  fPIDqa = new AliHFEpidQA;
  if(HasV0pidQA()) fPIDqa->SetV0pidQA();
  if(HasRecalculateTRDpid()) fPIDqa->SetRecalculateTRDpid();
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
  
}

void AliAnalysisTaskHFEpidQA::UserExec(Option_t *){
  //
  // Event Loop
  //
  if(fMCEvent) fPIDqa->SetMCEvent(fMCEvent);
  fPIDqa->SetRun((dynamic_cast<AliESDEvent*>(fInputEvent))->GetRunNumber());
  fPIDqa->SetT0((dynamic_cast<AliESDEvent*>(fInputEvent))->GetT0());
  fPIDqa->Process(fInputEvent);
  fEvents->Fill(1.1);
  PostData(1, fOutput);
}

void AliAnalysisTaskHFEpidQA::Terminate(Option_t *){
  //
  // Do Post Processing
  //
}

