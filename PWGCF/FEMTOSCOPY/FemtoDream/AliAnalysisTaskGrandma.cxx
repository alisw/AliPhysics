/*
 * AliAnalysisTaskGrandma.cxx
 *
 *  Created on: Sep 11, 2018
 *      Author: hohlweger
 */

#include "AliAnalysisTaskGrandma.h"
#include "AliLog.h"
ClassImp(AliAnalysisTaskGrandma)

AliAnalysisTaskGrandma::AliAnalysisTaskGrandma()
    : AliAnalysisTaskSE(),
      fIsMC(false),
      fQA(nullptr),
      fEvent(nullptr),
      fEvtCuts(nullptr),
      fEvtHistList(nullptr),
      fFemtoTrack(nullptr),
      fTrackCuts(nullptr),
      fTrackCutHistList(nullptr),
      fTrackCutHistMCList(nullptr),
      fAntiTrackCuts(nullptr),
      fAntiTrackCutHistList(nullptr),
      fAntiTrackCutHistMCList(nullptr) {
}

AliAnalysisTaskGrandma::AliAnalysisTaskGrandma(const char* name, bool isMC)
    : AliAnalysisTaskSE(name),
      fIsMC(isMC),
      fQA(nullptr),
      fEvent(nullptr),
      fEvtCuts(nullptr),
      fEvtHistList(nullptr),
      fFemtoTrack(nullptr),
      fTrackCuts(nullptr),
      fTrackCutHistList(nullptr),
      fTrackCutHistMCList(nullptr),
      fAntiTrackCuts(nullptr),
      fAntiTrackCutHistList(nullptr),
      fAntiTrackCutHistMCList(nullptr) {
  DefineOutput(1, TList::Class());  //Output for the Event Class and Pair Cleaner
  DefineOutput(2, TList::Class());  //Output for the Event Cuts
//  DefineOutput(3, TList::Class());  //Output for the Event Cuts
//  DefineOutput(4, TList::Class());  //Output for the Event Cuts
//  DefineOutput(5, TList::Class());  //Output for the Event Cuts
//  DefineOutput(6, TList::Class());  //Output for the Event Cuts
}

AliAnalysisTaskGrandma::~AliAnalysisTaskGrandma() {

}

void AliAnalysisTaskGrandma::UserCreateOutputObjects() {
  //first flag: Pile up rejection based on multiple vertices
  // (not used anymore since we event selection is done by AliEventCuts)
  //second flag: for the QA output of the AliEventCuts
  // might want to turn this off for systematics
  fEvent = new AliFemtoDreamEvent(true, true, GetCollisionCandidates());
  fFemtoTrack = new AliFemtoDreamTrack();
  fFemtoTrack->SetUseMCInfo(fIsMC);

  fQA = new TList();
  fQA->SetOwner();
  fQA->SetName("QA");
  fQA->Add(fEvent->GetEvtCutList());

  if (fEvtCuts) {
    fEvtCuts->InitQA();
    if (fEvtCuts->GetHistList()) {
      fEvtHistList = fEvtCuts->GetHistList();
    }
  } else {
    AliWarning("Event cuts are missing! \n");
  }
  PostData(1, fQA);
  PostData(2, fEvtHistList);
}
void AliAnalysisTaskGrandma::UserExec(Option_t *) {
  AliAODEvent *Event = static_cast<AliAODEvent*>(fInputEvent);
  if (!Event) {
    AliWarning("No Input Event");
  } else {
    fEvent->SetEvent(Event);
    if (fEvtCuts->isSelected(fEvent)) {
    }
  }
  PostData(1, fQA);
  PostData(2, fEvtHistList);
  return;
}
