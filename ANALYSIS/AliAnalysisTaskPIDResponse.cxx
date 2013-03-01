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
#include <AliESDpid.h>
#include <AliProdInfo.h>

#include "AliAnalysisTaskPIDResponse.h"

ClassImp(AliAnalysisTaskPIDResponse)

//______________________________________________________________________________
AliAnalysisTaskPIDResponse::AliAnalysisTaskPIDResponse():
AliAnalysisTaskSE(),
fIsMC(kFALSE),
fCachePID(kTRUE),
fOADBPath(),
fSpecialDetResponse(),
fPIDResponse(0x0),
fRun(0),
fOldRun(0),
fRecoPass(0),
fIsTunedOnData(kFALSE),
fRecoPassTuned(0),
fUseTPCEtaCorrection(kFALSE)//TODO: In future, default kTRUE  
{
  //
  // Dummy constructor
  //
}

//______________________________________________________________________________
AliAnalysisTaskPIDResponse::AliAnalysisTaskPIDResponse(const char* name):
AliAnalysisTaskSE(name),
fIsMC(kFALSE),
fCachePID(kTRUE),
fOADBPath(),
fSpecialDetResponse(),
fPIDResponse(0x0),
fRun(0),
fOldRun(0),
fRecoPass(0),
fIsTunedOnData(kFALSE),
fRecoPassTuned(0),
fUseTPCEtaCorrection(kFALSE)//TODO: In future, default kTRUE
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
  fPIDResponse->SetCachePID(fCachePID);
  if (!fOADBPath.IsNull()) fPIDResponse->SetOADBPath(fOADBPath.Data());

  if(fIsTunedOnData) fPIDResponse->SetTunedOnData(kTRUE,fRecoPassTuned);

  if (!fSpecialDetResponse.IsNull()){
    TObjArray *arr=fSpecialDetResponse.Tokenize("; ");
    for (Int_t i=0; i<arr->GetEntriesFast();++i){
      TString resp(arr->At(i)->GetName());
      if (resp.BeginsWith("TPC:")){
        resp.ReplaceAll("TPC:","");
        fPIDResponse->SetCustomTPCpidResponse(resp.Data());
        AliInfo(Form("Setting custom TPC response file: '%s'",resp.Data()));
      }
    }
    delete arr;
  }
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

  fPIDResponse->SetUseTPCEtaCorrection(fUseTPCEtaCorrection);
  fPIDResponse->InitialiseEvent(event,fRecoPass);
  AliESDpid *pidresp = dynamic_cast<AliESDpid*>(fPIDResponse);
  if(pidresp && AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()){
      pidresp->SetEventHandler(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  }
  //create and attach transient PID object
//   if (fCachePID) {
//     fPIDResponse->FillTrackDetectorPID();
//   }
}

//______________________________________________________________________________
void AliAnalysisTaskPIDResponse::SetRecoInfo()
{
  //
  // Set reconstruction information
  //
  
  //reset information
  fRecoPass=0;
  
  //Get UserInfo from the current tree 
  AliAnalysisManager *mgr=AliAnalysisManager::GetAnalysisManager();
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
  if (!inputHandler) return;

  TList *uiList = inputHandler->GetUserInfo();
  AliProdInfo prodInfo(uiList);
  prodInfo.List();

  TTree *tree= (TTree*)inputHandler->GetTree();
  TFile *file= (TFile*)tree->GetCurrentFile();
  if (!file) {
    AliError("Current file not found, cannot set reconstruction information");
    return;
  } else {
    TString fileName(file->GetName());
    fPIDResponse->SetCurrentFile(fileName.Data());
  }

  if (prodInfo.IsMC() == kTRUE) fIsMC=kTRUE;         // protection if user didn't use macro switch
  if ( (prodInfo.IsMC() == kFALSE) && (fIsMC == kFALSE) ) {      // reco pass is needed only for data
    fRecoPass = prodInfo.GetRecoPass();
    if (fRecoPass < 0) {   // as last resort we find pass from file name (UGLY, but not stored in ESDs/AODs before LHC12d )
      TString fileName(file->GetName());
      if (fileName.Contains("pass1") ) {
	fRecoPass=1;
      } else if (fileName.Contains("pass2") ) {
	fRecoPass=2;
      } else if (fileName.Contains("pass3") ) {
	fRecoPass=3;
      } else if (fileName.Contains("pass4") ) {
	fRecoPass=4;
      } else if (fileName.Contains("pass5") ) {
	fRecoPass=5;
      }
    } 
    if (fRecoPass <= 0) {
      AliError(" ******** Failed to find reconstruction pass number *********");
      AliError(" ******** PID information loaded for 'pass 0': parameters unreliable ******");
      AliError("      --> If these are MC data: please set kTRUE first argument of AddTaskPIDResponse");
      AliError("      --> If these are real data: please insert pass number inside the path of your local file ******");
    }
  }

}
