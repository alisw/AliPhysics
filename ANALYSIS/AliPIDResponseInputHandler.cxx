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

/* $Id: AliPIDResponseInputHandler 46193 2010-12-21 09:00:14Z wiechula $ */

//-----------------------------------------------------------------
//        Handler to set up the PID response object and
//         initialise it correctly for each event
//
//      Origin:
//        Jens Wiechula (jens.wiechula@cern.ch)
//        Martin Vala (martin.vala@cern.ch)

//-----------------------------------------------------------------


#include <TFile.h>
#include <TPRegexp.h>

#include <AliLog.h>
#include <AliVEvent.h>
#include "AliAnalysisManager.h"
#include "AliMultiInputEventHandler.h"
#include "AliPIDResponse.h"

#include "AliPIDResponseInputHandler.h"


ClassImp(AliPIDResponseInputHandler)

//_____________________________________________________________________________
AliPIDResponseInputHandler::AliPIDResponseInputHandler(const char *name) :
  AliInputEventHandler(name, name),
  fIsMC(kFALSE),
  fPIDResponse(0x0),
  fRun(0),
  fOldRun(0),
  fRecoPass(0),
  fMCurrentMutliIH(0)
{
//
// Default constructor.
//
   AliDebug(AliLog::kDebug + 10, "<-");
   AliDebug(AliLog::kDebug + 10, "->");
}

//_____________________________________________________________________________
AliPIDResponseInputHandler::~AliPIDResponseInputHandler()
{
//
// Destructor
//
   AliDebug(AliLog::kDebug + 10, "<-");
// 	delete fArrPidResponseMaster;
   AliDebug(AliLog::kDebug + 10, "->");
}

//_____________________________________________________________________________
Bool_t AliPIDResponseInputHandler::Init(Option_t *opt)
{
//
// Init() is called for all mix input handlers.
//
   AliDebug(AliLog::kDebug + 5, Form("<- opt=%s", opt));

   AliDebug(AliLog::kDebug + 5, Form("->"));
   return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliPIDResponseInputHandler::Init(TTree *tree, Option_t *opt)
{
//
// Init(const char*path) is called for all mix input handlers.
// Create event pool if needed
//
   AliDebug(AliLog::kDebug + 5, Form("<- %p %s opt=%s", (void *) tree, tree->GetName(), opt));

   if (fParentHandler) {
      TString tmp = "";
      AliInputEventHandler *ih = 0;
      fMCurrentMutliIH = dynamic_cast<AliMultiInputEventHandler*>(fParentHandler);
      if (fMCurrentMutliIH) {
         ih = fMCurrentMutliIH->GetFirstInputEventHandler();
         if (ih) {
            //pid response object
            ih->CreatePIDResponse(fIsMC);
            fPIDResponse = ih->GetPIDResponse();
            if (!fPIDResponse) AliFatal("PIDResponse object was not created");
         }
      }
   }

   AliDebug(AliLog::kDebug + 5, Form("->"));
   return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliPIDResponseInputHandler::Notify()
{
//
// Notify() is called for all mix input handlers
//
   AliDebug(AliLog::kDebug + 5, Form("<-"));
   AliDebug(AliLog::kDebug + 5, Form("->"));
   return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliPIDResponseInputHandler::Notify(const char *path)
{
//
// Notify(const char*path) is called for all mix input handlers
//
   AliDebug(AliLog::kDebug + 5, Form("<- %s", path));
   AliDebug(AliLog::kDebug + 5, "->");
   return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliPIDResponseInputHandler::BeginEvent(Long64_t entry)
{
//
// BeginEvent(Long64_t entry) is called for all mix input handlers
//
   AliDebug(AliLog::kDebug + 5, Form("<- %lld", entry));

   if (fParentHandler) {
      TString tmp = "";
      AliInputEventHandler *ih = 0;
      fMCurrentMutliIH = dynamic_cast<AliMultiInputEventHandler*>(fParentHandler);
      if (fMCurrentMutliIH) {
         ih = fMCurrentMutliIH->GetFirstInputEventHandler();
         if (ih) {
            //pid response object
            ih->CreatePIDResponse(fIsMC);
            fPIDResponse = ih->GetPIDResponse();
            if (!fPIDResponse) AliFatal("PIDResponse object was not created");

            AliVEvent *event = ih->GetEvent();
            if (!event) return kFALSE;
            fRun = event->GetRunNumber();

            if (fRun != fOldRun) {
              SetRecoInfo();
              fOldRun = fRun;
            }
           fPIDResponse->SetOADBPath(AliAnalysisManager::GetOADBPath());
           fPIDResponse->InitialiseEvent(event,fRecoPass);
         }
      }
   }
   AliDebug(AliLog::kDebug + 5, "->");
   return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliPIDResponseInputHandler::GetEntry()
{
   AliDebug(AliLog::kDebug + 5, "<-");
   AliDebug(AliLog::kDebug + 5, "->");
   return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliPIDResponseInputHandler::FinishEvent()
{
   //
   // FinishEvent() is called for all mix input handlers
   //
   AliDebug(AliLog::kDebug + 5, Form("<-"));
   AliDebug(AliLog::kDebug + 5, Form("->"));
   return kTRUE;
}

//_____________________________________________________________________________
void AliPIDResponseInputHandler::SetRecoInfo()
{
  //
  // Set reconstruction information
  //
  
  //reset information
  fRecoPass=0;
  
  //Get the current file to check the reconstruction pass (UGLY, but not stored in ESD... )
//   AliAnalysisManager *mgr=AliAnalysisManager::GetAnalysisManager();
  AliVEventHandler *inputHandler=fMCurrentMutliIH->GetFirstInputEventHandler();
  if (!inputHandler) return;
  
  TTree *tree= (TTree*)inputHandler->GetTree();
  TFile *file= (TFile*)tree->GetCurrentFile();
  
  if (!file) {
    AliError("Current file not found, cannot set reconstruction information");
    return;
  }
  
  //find pass from file name (UGLY, but not stored in ESD... )
  TString fileName(file->GetName());
  if (fileName.Contains("/pass1")) {
    fRecoPass=1;
  } else if (fileName.Contains("/pass2")) {
    fRecoPass=2;
  } else if (fileName.Contains("/pass3")) {
    fRecoPass=3;
  }
}
