#include "AliAnalysisTaskNucleiSkimAOD.h"

#include <TChain.h>

#include <AliAnalysisManager.h>
#include <AliAODEvent.h>
#include <AliAODHandler.h>
#include <AliAODHeader.h>
#include <AliAODTrack.h>
#include <AliInputEventHandler.h>
#include <AliPIDResponse.h>
#include <AliTOFPIDResponse.h>

AliAnalysisTaskNucleiSkimAOD::AliAnalysisTaskNucleiSkimAOD(const char* name) :
  AliAnalysisTaskSE(name),
  mEventCuts{},
  mPIDresponse{nullptr},
  mOutputList{nullptr}
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}


AliAnalysisTaskNucleiSkimAOD::~AliAnalysisTaskNucleiSkimAOD() {
  delete mPIDresponse;
  delete mOutputList;
}

void AliAnalysisTaskNucleiSkimAOD::UserCreateOutputObjects() {
  // Creation of the histograms, this is called once
  //

  if (!mOutputList) mOutputList = new TList();
  mOutputList->SetOwner();

  //
  // Get PID response object
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  if(!man)
    AliFatal("Could not find manager");
  AliInputEventHandler* handler = dynamic_cast<AliInputEventHandler*> (man->GetInputEventHandler());
  if(!handler)
    AliFatal("No input event handler");
  mPIDresponse = dynamic_cast<AliPIDResponse*>(handler->GetPIDResponse());
  if (!mPIDresponse)
    AliFatal("PIDResponse object was not created"); // Escalated to fatal. Task is unusable without PID response.

  //
  //Define output histograms
  //

  mEventCuts.AddQAplotsToList(mOutputList);
  PostData(1,mOutputList);

}

void AliAnalysisTaskNucleiSkimAOD::UserExec(Option_t *){

  AliAODEvent *event = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!event) return;
  if (!mEventCuts.AcceptEvent(event)) {
    PostData(1, mOutputList);
    return;
  }

  bool selectedEvent = false;

  for (int iEv = 0;iEv < event->GetNumberOfTracks(); ++iEv) {
    AliAODTrack *track = dynamic_cast<AliAODTrack*>(event->GetTrack(iEv));
    selectedEvent = mPIDresponse->NumberOfSigmasTPC(track,AliPID::kHe3) > -5;
    if (selectedEvent) break;
    selectedEvent = mPIDresponse->NumberOfSigmasTPC(track,AliPID::kDeuteron) > -5 && track->Pt() < 1.4 && track->P() < 1.4;
    if (selectedEvent) break;
    selectedEvent = mPIDresponse->NumberOfSigmasTOF(track,AliPID::kDeuteron) > -10;
    if (selectedEvent) break;
  }

  AliAODHandler *oh = (AliAODHandler *)AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
  if (oh)
    oh->SetFillAOD(kFALSE);

  if (selectedEvent && oh)
  {
    oh->SetFillAOD(kTRUE);
    AliAODEvent *eout = dynamic_cast<AliAODEvent *>(oh->GetAOD());
    AliAODEvent *evin = dynamic_cast<AliAODEvent *>(InputEvent());
    TTree *tout = oh->GetTree();
    if (tout)
    {
      TList *lout = tout->GetUserInfo();
      if (lout->FindObject("alirootVersion") == 0)
      {
        TList *lin = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetUserInfo();
        for (Int_t jj = 0; jj < lin->GetEntries() - 1; ++jj)
        {
          lout->Add(lin->At(jj)->Clone(lin->At(jj)->GetName()));
        }
      }
    }

    AliAODHeader *outH = (AliAODHeader *)eout->GetHeader();
    AliAODHeader *inH = (AliAODHeader *)evin->GetHeader();
    *outH = *inH;

    AliTOFHeader *outTOFh = const_cast<AliTOFHeader *>(eout->GetTOFHeader());
    *outTOFh = *evin->GetTOFHeader();
    *eout->GetVZEROData() = *evin->GetVZEROData();
    *eout->GetTZEROData() = *evin->GetTZEROData();
    new (eout->GetTracks()) TClonesArray(*evin->GetTracks());
    new (eout->GetVertices()) TClonesArray(*evin->GetVertices());
    new (eout->GetV0s()) TClonesArray(*evin->GetV0s());
  }

  PostData(1, mOutputList);
}
