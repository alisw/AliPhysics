/*
 * AliAnalysisTaskFemtoPlotration.cxx
 *
 *  Created on: 24 Nov 2017
 *      Author: bernhardhohlweger
 */

#include "AliAnalysisTaskFemtoDream.h"
#include "AliLog.h"
ClassImp(AliAnalysisTaskFemtoDream)

AliAnalysisTaskFemtoDream::AliAnalysisTaskFemtoDream()
:AliAnalysisTaskSE()
,fTrackBufferSize(0)
,fESDAnalysis(false)
,fMinBookingME(false)
,fMinBookingSample(false)
,fMVPileUp(false)
,fEvtCutQA(false)
,fIsMC(false)
,fAnalysis()
,fQA(0)
,fEvtCuts()
,fEvtHistList(0)
,fTrackCuts()
,fTrackCutHistList(0)
,fTrackCutHistMCList(0)
,fAntiTrackCuts()
,fAntiTrackCutHistList(0)
,fAntiTrackCutHistMCList(0)
,fv0Cuts()
,fv0CutHistList(0)
,fv0CutHistMCList(0)
,fAntiv0Cuts()
,fAntiv0CutHistList(0)
,fAntiv0CutHistMCList(0)
,fCascCuts()
,fCascCutList(0)
,fCascCutMCList(0)
,fAntiCascCuts()
,fAntiCascCutList(0)
,fAntiCascCutMCList(0)
,fConfig()
,fResults(0)
,fResultQA(0)
,fResultsSample(0)
,fResultQASample(0)
{}

AliAnalysisTaskFemtoDream::AliAnalysisTaskFemtoDream(const char *name,bool isESD,bool isMC)
:AliAnalysisTaskSE(name)
,fTrackBufferSize(0)
,fESDAnalysis(isESD)
,fMinBookingME(false)
,fMinBookingSample(false)
,fMVPileUp(false)
,fEvtCutQA(false)
,fIsMC(isMC)
,fAnalysis()
,fQA(0)
,fEvtCuts()
,fEvtHistList(0)
,fTrackCuts()
,fTrackCutHistList(0)
,fTrackCutHistMCList(0)
,fAntiTrackCuts()
,fAntiTrackCutHistList(0)
,fAntiTrackCutHistMCList(0)
,fv0Cuts()
,fv0CutHistList(0)
,fv0CutHistMCList(0)
,fAntiv0Cuts()
,fAntiv0CutHistList(0)
,fAntiv0CutHistMCList(0)
,fCascCuts()
,fCascCutList(0)
,fCascCutMCList(0)
,fAntiCascCuts()
,fAntiCascCutList(0)
,fAntiCascCutMCList(0)
,fConfig()
,fResults(0)
,fResultQA(0)
,fResultsSample(0)
,fResultQASample(0)
{
  DefineOutput(1, TList::Class());  //Output for the Event Class and Pair Cleaner
  DefineOutput(2, TList::Class());  //Output for the Event Cuts
  DefineOutput(3, TList::Class());  //Output for the Track Cuts
  DefineOutput(4, TList::Class());  //Output for the Antitrack Cuts
  DefineOutput(5, TList::Class());  //Output for the v0 Cuts
  DefineOutput(6, TList::Class());  //Output for the Antiv0 Cuts
  DefineOutput(7, TList::Class());  //Output for the Cascade Cuts
  DefineOutput(8, TList::Class());  //Output for the Anti Cascade
  DefineOutput(9, TList::Class());  //Output for the Results
  DefineOutput(10, TList::Class());  //Output for the Results QA
  DefineOutput(11, TList::Class());  //Output for the Sample
  DefineOutput(12, TList::Class());  //Output for the Sample QA
  if (fIsMC) {
    DefineOutput(13, TList::Class());  //Output for the Track Cut MC Info
    DefineOutput(14, TList::Class()); //Output for the AntiTrack Cut MC Info
    DefineOutput(15, TList::Class()); //Output for the v0 Cut MC Info
    DefineOutput(16, TList::Class()); //Output for the Antiv0 Cut MC Info
    DefineOutput(17, TList::Class()); //Output for the Xi Cut MC Info
    DefineOutput(18, TList::Class()); //Output for the AntiXi Cut MC Info
  }
}

AliAnalysisTaskFemtoDream::~AliAnalysisTaskFemtoDream() {
  if (fAnalysis) {
    delete fAnalysis;
  }
}

void AliAnalysisTaskFemtoDream::UserCreateOutputObjects() {
  fAnalysis=new AliFemtoDreamAnalysis();
  fAnalysis->SetTrackBufferSize(fTrackBufferSize);
  fAnalysis->SetMVPileUp(fMVPileUp);
  fAnalysis->SetEvtCutQA(fEvtCutQA);

  if (fEvtCuts) {
    fAnalysis->SetEventCuts(fEvtCuts);
  } else {
    AliFatal("Event cuts missing!");
  }
  if (fTrackCuts) {
    fAnalysis->SetTrackCuts(fTrackCuts);
  } else {
    AliFatal("TrackCuts missing!");
  }
  if (fAntiTrackCuts) {
    fAnalysis->SetAntiTrackCuts(fAntiTrackCuts);
  } else {
    AliFatal("Antitrack cuts missing");
  }
  if (fv0Cuts) {
    fAnalysis->Setv0Cuts(fv0Cuts);
  } else {
    AliFatal("v0 cuts missing");
  }
  if (fAntiv0Cuts) {
    fAnalysis->SetAntiv0Cuts(fAntiv0Cuts);
  } else {
    AliFatal("Antiv0 cuts missing");
  }
  if (fCascCuts) {
    fAnalysis->SetCascadeCuts(fCascCuts);
  } else {
    AliFatal("Cascade cuts missing");
  }
  if (fAntiCascCuts) {
    fAnalysis->SetAntiCascadeCuts(fAntiCascCuts);
  } else {
    AliFatal("AntiCascade cuts missing");
  }
  if (fConfig) {
    fAnalysis->SetCollectionConfig(fConfig);
    fMinBookingME=fConfig->GetMinimalBookingME();
    fMinBookingSample=fConfig->GetMinimalBookingSample();
  } else {
    AliFatal("Event Collection Config missing");
  }

  fAnalysis->Init(fIsMC,GetCollisionCandidates());
  if ((!fMinBookingME)||(!fMinBookingSample)) {
    if (fAnalysis->GetQAList()) {
      fQA=fAnalysis->GetQAList();
    } else {
      fQA=new TList();
      fQA->SetName("QA");
      fQA->SetOwner();
    }
  } else {
    fQA=new TList();
    fQA->SetName("QA");
    fQA->SetOwner();
  }

  if (!fMinBookingME) {
    if (fAnalysis->GetResultQAList()) {
      fResultQA=fAnalysis->GetResultQAList();
    } else {
      fResultQA=new TList();
      fResultQA->SetName("ResultQA");
      fResultQA->SetOwner();
    }
  } else {
    fResultQA=new TList();
    fResultQA->SetName("ResultQA");
    fResultQA->SetOwner();
  }

  if (!fMinBookingSample) {
    if (fAnalysis->GetResultSampleQAList()) {
      fResultQASample=fAnalysis->GetResultSampleQAList();
    }else {
      fResultQASample=new TList();
      fResultQASample->SetOwner();
      fResultQASample->SetName("ResultsQASample");
    }
  } else {
    fResultQASample=new TList();
    fResultQASample->SetOwner();
    fResultQASample->SetName("ResultsQASample");
  }

  if (!fEvtCuts->GetMinimalBooking()) {
    if (fAnalysis->GetEventCutHists()) {
      fEvtHistList=fAnalysis->GetEventCutHists();
    }
  } else {
    fEvtHistList=new TList();
    fEvtHistList->SetName("EventCuts");
    fEvtHistList->SetOwner();
  }
  if (fAnalysis->GetTrackCutHists()) {
    fTrackCutHistList=fAnalysis->GetTrackCutHists();
  } else {
    fTrackCutHistList=new TList();
    fTrackCutHistList->SetName("TrackCuts");
    fTrackCutHistList->SetOwner();
  }
  if (fAnalysis->GetAntitrackCutHists()) {
    fAntiTrackCutHistList=fAnalysis->GetAntitrackCutHists();
  } else {
    fAntiTrackCutHistList=new TList();
    fAntiTrackCutHistList->SetName("AntiTrackCuts");
    fAntiTrackCutHistList->SetOwner();
  }
  if (fAnalysis->Getv0CutHist()) {
    fv0CutHistList=fAnalysis->Getv0CutHist();
  } else {
    fv0CutHistList=new TList();
    fv0CutHistList->SetName("v0Cuts");
    fv0CutHistList->SetOwner();
  }
  if (fAnalysis->GetAntiv0CutHist()) {
    fAntiv0CutHistList=fAnalysis->GetAntiv0CutHist();
  } else {
    fAntiv0CutHistList=new TList();
    fAntiv0CutHistList->SetName("Antiv0Cuts");
    fAntiv0CutHistList->SetOwner();
  }
  if (fAnalysis->GetCascadeCutHist()) {
    fCascCutList=fAnalysis->GetCascadeCutHist();
  } else {
    fCascCutList=new TList();
    fCascCutList->SetName("CascadeCuts");
    fCascCutList->SetOwner();
  }
  if (fAnalysis->GetAntiCascadeCutHist()) {
    fAntiCascCutList=fAnalysis->GetAntiCascadeCutHist();
  } else {
    fAntiCascCutList=new TList();
    fAntiCascCutList->SetName("AntiCascadeCuts");
    fAntiCascCutList->SetOwner();
  }
  //Results we always post
  if (fAnalysis->GetResultList()) {
    fResults=fAnalysis->GetResultList();
  } else {
    AliWarning("Results List not Available");
    fResults=new TList();
    fResults->SetOwner();
    fResults->SetName("ResultsSample");
  }
  if (fAnalysis->GetResultSampleList()) {
    fResultsSample=fAnalysis->GetResultSampleList();
  } else {
    AliWarning("Results Sample List not Available");
    fResultsSample=new TList();
    fResultsSample->SetOwner();
    fResultsSample->SetName("ResultsSample");
  }
  PostData(1,fQA);
  PostData(2,fEvtHistList);
  PostData(3,fTrackCutHistList);
  PostData(4,fAntiTrackCutHistList);
  PostData(5,fv0CutHistList);
  PostData(6,fAntiv0CutHistList);
  PostData(7,fCascCutList);
  PostData(8,fAntiCascCutList);
  PostData(9,fResults);
  PostData(10,fResultQA);
  PostData(11,fResultsSample);
  PostData(12,fResultQASample);
  if (fIsMC) {
    if (!fTrackCuts->GetMinimalBooking()) {
      if (fTrackCuts->GetIsMonteCarlo()) {
        fTrackCutHistMCList=fTrackCuts->GetMCQAHists();
      }
    } else {
      fTrackCutHistMCList=new TList();
      fTrackCutHistMCList->SetName("MCTrkCuts");
      fTrackCutHistMCList->SetOwner();
    }
    if (!fAntiTrackCuts->GetMinimalBooking()) {
      if (fAntiTrackCuts->GetIsMonteCarlo()) {
        fAntiTrackCutHistMCList=fAntiTrackCuts->GetMCQAHists();
      }
    } else {
      fAntiTrackCutHistMCList=new TList();
      fAntiTrackCutHistMCList->SetName("MCAntiTrkCuts");
      fAntiTrackCutHistMCList->SetOwner();
    }
    if (!fv0Cuts->GetMinimalBooking()) {
      if (fv0Cuts->GetIsMonteCarlo()) {
        fv0CutHistMCList=fv0Cuts->GetMCQAHists();
      }
    } else {
      fv0CutHistMCList=new TList();
      fv0CutHistMCList->SetName("MCv0Cuts");
      fv0CutHistMCList->SetOwner();
    }
    if (!fAntiv0Cuts->GetMinimalBooking()) {
      if (fAntiv0Cuts->GetIsMonteCarlo()) {
        fAntiv0CutHistMCList=fAntiv0Cuts->GetMCQAHists();
      }
    } else {
      fAntiv0CutHistMCList=new TList();
      fAntiv0CutHistMCList->SetName("MCAntiv0Cuts");
      fAntiv0CutHistMCList->SetOwner();
    }
    if (!fCascCuts->GetMinimalBooking()) {
      if (fCascCuts->GetIsMonteCarlo()) {
        fCascCutMCList=fCascCuts->GetMCQAHists();
      }
    } else {
      fCascCutMCList=new TList();
      fCascCutMCList->SetName("MCCascCuts");
      fCascCutMCList->SetOwner();
    }
    if (!fAntiCascCuts->GetMinimalBooking()) {
      if (fAntiCascCuts->GetIsMonteCarlo()) {
        fAntiCascCutMCList=fAntiCascCuts->GetMCQAHists();
      }
    } else {
      fAntiCascCutMCList=new TList();
      fAntiCascCutMCList->SetName("MCAntiCascCuts");
      fAntiCascCutMCList->SetOwner();
    }
    PostData(13,fTrackCutHistMCList);
    PostData(14,fAntiTrackCutHistMCList);
    PostData(15,fv0CutHistMCList);
    PostData(16,fAntiv0CutHistMCList);
    PostData(17,fCascCutMCList);
    PostData(18,fAntiCascCutMCList);
  }
}

void AliAnalysisTaskFemtoDream::UserExec(Option_t *) {
  if (fESDAnalysis) {
    AliESDEvent *Event=static_cast<AliESDEvent*>(fInputEvent);
    if (!Event) {
      AliWarning("No Input Event");
    } else {
      fAnalysis->Make(Event);
    }
  } else {
    AliAODEvent *Event=static_cast<AliAODEvent*>(fInputEvent);
    if (!Event) {
      AliWarning("No Input Event");
    } else {
      fAnalysis->Make(Event);
    }
  }
  PostData(1,fQA);
  PostData(2,fEvtHistList);
  PostData(3,fTrackCutHistList);
  PostData(4,fAntiTrackCutHistList);
  PostData(5,fv0CutHistList);
  PostData(6,fAntiv0CutHistList);
  PostData(7,fCascCutList);
  PostData(8,fAntiCascCutList);
  PostData(9,fResults);
  PostData(10,fResultQA);
  PostData(11,fResultsSample);
  PostData(12,fResultQASample);
  if (fIsMC) {
    PostData(13,fTrackCutHistMCList);
    PostData(14,fAntiTrackCutHistMCList);
    PostData(15,fv0CutHistMCList);
    PostData(16,fAntiv0CutHistMCList);
    PostData(17,fCascCutMCList);
    PostData(18,fAntiCascCutMCList);
  }
}
