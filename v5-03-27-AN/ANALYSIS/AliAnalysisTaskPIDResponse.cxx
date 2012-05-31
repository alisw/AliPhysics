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

/* $Id: AliAnalysisTaskPIDResponse.cxx 43811 2010-09-23 14:13:31Z wiechula $ */
#include <TFile.h>
#include <TChain.h>

#include <AliAnalysisManager.h>
#include <AliInputEventHandler.h>
#include <AliVEventHandler.h>
#include <AliVEvent.h>
#include <AliVParticle.h>
#include <AliVTrack.h>
#include <AliLog.h>
#include <AliPIDResponse.h>

#include "AliAnalysisTaskPIDResponse.h"

ClassImp(AliAnalysisTaskPIDResponse)

//______________________________________________________________________________
AliAnalysisTaskPIDResponse::AliAnalysisTaskPIDResponse():
AliAnalysisTaskSE(),
fIsMC(kFALSE),
fOADBPath(),
fPIDResponse(0x0),
fRun(0),
fOldRun(0),
fRecoPass(0)
{
  //
  // Dummy constructor
  //
}

//______________________________________________________________________________
AliAnalysisTaskPIDResponse::AliAnalysisTaskPIDResponse(const char* name):
AliAnalysisTaskSE(name),
fIsMC(kFALSE),
fOADBPath(),
fPIDResponse(0x0),
fRun(0),
fOldRun(0),
fRecoPass(0)
{
  //
  // Default constructor
  //
  DefineInput(0,TChain::Class());
}

//______________________________________________________________________________
AliAnalysisTaskPIDResponse::~AliAnalysisTaskPIDResponse()
{
  //
  // Destructor
  //
}

//______________________________________________________________________________
void AliAnalysisTaskPIDResponse::UserCreateOutputObjects()
{
  //
  // Create the output QA objects
  //
  
  AliLog::SetClassDebugLevel("AliAnalysisTaskPIDResponse",10);
  
  //input hander
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler=dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
  if (!inputHandler) AliFatal("Input handler needed");
  
  //pid response object
  inputHandler->CreatePIDResponse(fIsMC);
  fPIDResponse=inputHandler->GetPIDResponse();
  if (!fPIDResponse) AliFatal("PIDResponse object was not created");

  fPIDResponse->SetOADBPath(AliAnalysisManager::GetOADBPath());
  if (!fOADBPath.IsNull()) fPIDResponse->SetOADBPath(fOADBPath.Data());
}

//______________________________________________________________________________
void AliAnalysisTaskPIDResponse::UserExec(Option_t */*option*/)
{
  // Setup the PID response functions and fill the QA histograms
  //
  AliVEvent *event=InputEvent();
  if (!event) return;
  fRun=event->GetRunNumber();
  
  if (fRun!=fOldRun){
    SetRecoInfo();
    fOldRun=fRun;
  }

  fPIDResponse->InitialiseEvent(event,fRecoPass);
}

//______________________________________________________________________________
void AliAnalysisTaskPIDResponse::SetRecoInfo()
{
  //
  // Set reconstruction information
  //
  
  //reset information
  fRecoPass=0;
  
  //Get the current file to check the reconstruction pass (UGLY, but not stored in ESD... )
  AliAnalysisManager *mgr=AliAnalysisManager::GetAnalysisManager();
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
  if (!inputHandler) return;
  
  TTree *tree= (TTree*)inputHandler->GetTree();
  TFile *file= (TFile*)tree->GetCurrentFile();
  
  if (!file) {
    AliError("Current file not found, cannot set reconstruction information");
    return;
  }
  
  //find pass from file name (UGLY, but not stored in ESD... )
  TString fileName(file->GetName());
  if (fileName.Contains("pass1") ) {
    fRecoPass=1;
  } else if (fileName.Contains("pass2") ) {
    fRecoPass=2;
  } else if (fileName.Contains("pass3") ) {
    fRecoPass=3;
  }

  fPIDResponse->SetCurrentFile(fileName.Data());
}
