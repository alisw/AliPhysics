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
,fMVPileUp(false)
,fEvtCutQA(false)
,fIsMC(false)
,fAnalysis()
,fQA()
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
,fAntiCascCuts()
,fAntiCascCutList(0)
,fConfig()
,fResults()
,fResultQA()
{}

AliAnalysisTaskFemtoDream::AliAnalysisTaskFemtoDream(const char *name,bool isMC)
:AliAnalysisTaskSE(name)
,fTrackBufferSize(0)
,fMVPileUp(false)
,fEvtCutQA(false)
,fIsMC(isMC)
,fAnalysis()
,fQA()
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
,fAntiCascCuts()
,fAntiCascCutList(0)
,fConfig()
,fResults()
,fResultQA()
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
  if (fIsMC) {
    DefineOutput(11, TList::Class());  //Output for the Track Cut MC Info
    DefineOutput(12, TList::Class()); //Output for the AntiTrack Cut MC Info
    DefineOutput(13, TList::Class()); //Output for the v0 Cut MC Info
    DefineOutput(14, TList::Class()); //Output for the Antiv0 Cut MC Info
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
  } else {
    AliFatal("Event Collection Config missing");
  }
  fAnalysis->Init();

  if (fAnalysis->GetQAList()) {
    fQA=fAnalysis->GetQAList();
  } else {
    AliFatal("QA Histograms not available");
  }
  if (fAnalysis->GetEventCutHists()) {
    fEvtHistList=fAnalysis->GetEventCutHists();
  } else {
    AliFatal("Event Cut Histograms not available");
  }
  if (fAnalysis->GetTrackCutHists()) {
    fTrackCutHistList=fAnalysis->GetTrackCutHists();
  } else {
    AliFatal("Event Cut Histograms not available");
  }
  if (fAnalysis->GetAntitrackCutHists()) {
    fAntiTrackCutHistList=fAnalysis->GetAntitrackCutHists();
  } else {
    AliFatal("Event Cut Histograms not available");
  }
  if (fAnalysis->Getv0CutHist()) {
//      &&
//      fAnalysis->Getv0PosDaugHist()&&fAnalysis->Getv0NegDaugHist()) {
    fv0CutHistList=fAnalysis->Getv0CutHist();
//    TList *PosDaug=fAnalysis->Getv0PosDaugHist();
//    PosDaug->SetName("ProtonDaughter");
//    fv0CutHistList->Add(PosDaug);
//    TList *NegDaug=fAnalysis->Getv0NegDaugHist();
//    NegDaug->SetName("PionDaughter");
//    fv0CutHistList->Add(NegDaug);
  } else {
    AliFatal("v0 Cut Histograms not available");
  }
  if (fAnalysis->GetAntiv0CutHist()) {
//  &&
//      fAnalysis->GetAntiv0PosDaugHist()&&fAnalysis->GetAntiv0NegDaugHist()) {
    fAntiv0CutHistList=fAnalysis->GetAntiv0CutHist();
//    TList *PosDaug=fAnalysis->GetAntiv0PosDaugHist();
//    PosDaug->SetName("PionDaughter");
//    fAntiv0CutHistList->Add(PosDaug);
//    TList *NegDaug=fAnalysis->GetAntiv0NegDaugHist();
//    NegDaug->SetName("AntiProtonDaughter");
//    fAntiv0CutHistList->Add(NegDaug);
  } else {
    AliFatal("Antiv0 Cut Histograms not available");
  }
  if (fAnalysis->GetCascadeCutHist()) {
    fCascCutList=fAnalysis->GetCascadeCutHist();
  } else {
    AliFatal("Cascade Cut Histograms not availabl");
  }
  if (fAnalysis->GetAntiCascadeCutHist()) {
    fAntiCascCutList=fAnalysis->GetAntiCascadeCutHist();
  } else {
    AliFatal("AntiCascade Cut Histograms not available");
  }
  if (fAnalysis->GetResultList()) {
    fResults=fAnalysis->GetResultList();
  } else {
    AliFatal("Results List not Available");
  }
  if (fAnalysis->GetResultQAList()) {
    fResultQA=fAnalysis->GetResultQAList();
  } else {
    AliFatal("Results QA List not Available");
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

  if (fIsMC) {
    if (fTrackCuts->GetMCQAHists()) {
      fTrackCutHistMCList=fTrackCuts->GetMCQAHists();
    } else {
      AliFatal("No Track Cut MC Histograms!");
    }
    if (fAntiTrackCuts->GetMCQAHists()) {
      fAntiTrackCutHistMCList=fAntiTrackCuts->GetMCQAHists();
    } else {
      AliFatal("No Antitrack Cut MC Histograms!");
    }
    if (fv0Cuts->GetMCQAHists()) {
      fv0CutHistMCList=fv0Cuts->GetMCQAHists();
    } else {
      AliFatal("No v0 cut MC Histograms");
    }
    if (fAntiv0Cuts->GetMCQAHists()) {
      fAntiv0CutHistMCList=fAntiv0Cuts->GetMCQAHists();
    } else {
      AliFatal("No Antiv0 cut MC Histograms");
    }
    PostData(11,fTrackCutHistMCList);
    PostData(12,fAntiTrackCutHistMCList);
    PostData(13,fv0CutHistMCList);
    PostData(14,fAntiv0CutHistMCList);
  }
}

void AliAnalysisTaskFemtoDream::UserExec(Option_t *) {
  AliAODEvent *Event=static_cast<AliAODEvent*>(fInputEvent);
  if (!Event) {
    AliFatal("No Input Event");
  } else {
    fAnalysis->Make(Event);
    if (fAnalysis->GetQAList()) {
       fQA=fAnalysis->GetQAList();
     } else {
       AliFatal("QA Histograms not available");
     }
    if (fAnalysis->GetEventCutHists()) {
      fEvtHistList=fAnalysis->GetEventCutHists();
    } else {
      AliFatal("Event Cut Histograms not available");
    }
    if (fAnalysis->GetTrackCutHists()) {
      fTrackCutHistList=fAnalysis->GetTrackCutHists();
    } else {
      AliFatal("Event Cut Histograms not available");
    }
    if (fAnalysis->GetAntitrackCutHists()) {
      fAntiTrackCutHistList=fAnalysis->GetAntitrackCutHists();
    } else {
      AliFatal("Event Cut Histograms not available");
    }
    if (fAnalysis->Getv0CutHist()) {//&&
//        fAnalysis->Getv0PosDaugHist()&&fAnalysis->Getv0NegDaugHist()) {
      fv0CutHistList=fAnalysis->Getv0CutHist();
//      TList *PosDaug=fAnalysis->Getv0PosDaugHist();
//      PosDaug->SetName("ProtonDaughter");
//      fv0CutHistList->Add(PosDaug);
//      TList *NegDaug=fAnalysis->Getv0NegDaugHist();
//      NegDaug->SetName("PionDaughter");
//      fv0CutHistList->Add(NegDaug);
    } else {
      AliFatal("v0 Cut Histograms not available");
    }
    if (fAnalysis->GetAntiv0CutHist()) {//&&
//        fAnalysis->GetAntiv0PosDaugHist()&&fAnalysis->GetAntiv0NegDaugHist()) {
      fAntiv0CutHistList=fAnalysis->GetAntiv0CutHist();
//      TList *PosDaug=fAnalysis->GetAntiv0PosDaugHist();
//      PosDaug->SetName("PionDaughter");
//      fAntiv0CutHistList->Add(PosDaug);
//      TList *NegDaug=fAnalysis->GetAntiv0NegDaugHist();
//      NegDaug->SetName("AntiProtonDaughter");
//      fAntiv0CutHistList->Add(NegDaug);
    } else {
      AliFatal("Antiv0 Cut Histograms not available");
    }
    if (fAnalysis->GetCascadeCutHist()) {
      fCascCutList=fAnalysis->GetCascadeCutHist();
    } else {
      AliFatal("Cascade Cut Histograms not availabl");
    }
    if (fAnalysis->GetAntiCascadeCutHist()) {
      fAntiCascCutList=fAnalysis->GetAntiCascadeCutHist();
    } else {
      AliFatal("AntiCascade Cut Histograms not availabl");
    }
    if (fAnalysis->GetResultList()) {
      fResults=fAnalysis->GetResultList();
    } else {
      AliFatal("Results List not Available");
    }
    if (fAnalysis->GetResultQAList()) {
      fResultQA=fAnalysis->GetResultQAList();
    } else {
      AliFatal("Results QA List not Available");
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
    if (fIsMC) {
      if (fTrackCuts->GetMCQAHists()) {
        fTrackCutHistMCList=fTrackCuts->GetMCQAHists();
      } else {
        AliFatal("No Track Cut MC Histograms!");
      }
      if (fAntiTrackCuts->GetMCQAHists()) {
        fAntiTrackCutHistMCList=fAntiTrackCuts->GetMCQAHists();
      } else {
        AliFatal("No Antitrack Cut MC Histograms!");
      }
      if (fv0Cuts->GetMCQAHists()) {
        fv0CutHistMCList=fv0Cuts->GetMCQAHists();
      } else {
        AliFatal("No v0 cut MC Histograms");
      }
      if (fAntiv0Cuts->GetMCQAHists()) {
        fAntiv0CutHistMCList=fAntiv0Cuts->GetMCQAHists();
      } else {
        AliFatal("No Antiv0 cut MC Histograms");
      }
      PostData(11,fTrackCutHistMCList);
      PostData(12,fAntiTrackCutHistMCList);
      PostData(13,fv0CutHistMCList);
      PostData(14,fAntiv0CutHistMCList);
    }
  }
}
