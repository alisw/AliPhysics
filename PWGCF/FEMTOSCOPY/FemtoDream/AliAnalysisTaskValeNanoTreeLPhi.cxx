#include "AliAnalysisTaskValeNanoTreeLPhi.h"
#include "AliFemtoDreamBasePart.h"
#include "AliLog.h"
#include "AliNanoAODTrack.h"
#include "AliVEvent.h"
#include "TH1F.h"
#include "TList.h"
#include "AliAODMCParticle.h"

ClassImp(AliAnalysisTaskValeNanoTreeLPhi)
    AliAnalysisTaskValeNanoTreeLPhi::AliAnalysisTaskValeNanoTreeLPhi()
    : AliAnalysisTaskSE(),
      fIsMC(false),
      fUseOMixing(false),
      fInvMassCutSBdown(0.0),
      fInvMassCutSBup(0.0),
      fTrigger(AliVEvent::kINT7),
      fQA(nullptr),
      fEvtList(nullptr),
      fLambdaList(nullptr),
      fAntiLambdaList(nullptr),
      fKaonPlusList(nullptr),
      fKaonMinusList(nullptr),
      fPhiList(nullptr),
      fResults(nullptr),
      fResultsQA(nullptr),
      fInputEvent(nullptr),
      fEvent(nullptr),
      fTrack(nullptr),
      fLambda(nullptr),
      fPhiParticle(nullptr),
      fEventCuts(nullptr),
      fPosKaonCuts(nullptr),
      fNegKaonCuts(nullptr),
      fPhiCuts(nullptr),
      fLambdaCuts(nullptr),
      fAntiLambdaCuts(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fSample(nullptr),
      fGTI(nullptr),
      fTrackBufferSize()
{
}

AliAnalysisTaskValeNanoTreeLPhi::AliAnalysisTaskValeNanoTreeLPhi(
    const char *name, bool isMC)
    : AliAnalysisTaskSE(name),
      fIsMC(isMC),
      fUseOMixing(false),
      fInvMassCutSBdown(0.0),
      fInvMassCutSBup(0.0),
      fTrigger(AliVEvent::kINT7),
      fQA(nullptr),
      fEvtList(nullptr),
      fLambdaList(nullptr),
      fAntiLambdaList(nullptr),
      fKaonPlusList(nullptr),
      fKaonMinusList(nullptr),
      fPhiList(nullptr),
      fResults(nullptr),
      fResultsQA(nullptr),
      fInputEvent(nullptr),
      fEvent(nullptr),
      fTrack(nullptr),
      fLambda(nullptr),
      fPhiParticle(nullptr),
      fEventCuts(nullptr),
      fPosKaonCuts(nullptr),
      fNegKaonCuts(nullptr),
      fPhiCuts(nullptr),
      fLambdaCuts(nullptr),
      fAntiLambdaCuts(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fSample(nullptr),
      fGTI(nullptr),
      fTrackBufferSize(2000)
{
  DefineOutput(1, TList::Class()); //Output for the Event Class and Pair Cleaner
  DefineOutput(2, TList::Class()); //Output for the Event Cuts
  DefineOutput(3, TList::Class()); //Output for the Lambda Cuts
  DefineOutput(4, TList::Class()); //Output for the AntiLambda Cuts
  DefineOutput(5, TList::Class()); //Output for the KaonPlus Cuts
  DefineOutput(6, TList::Class()); //Output for the KaonMinus Cuts
  DefineOutput(7, TList::Class()); //Output for the Phi Cuts
  DefineOutput(8, TList::Class()); //Output for the Results
  DefineOutput(9, TList::Class()); //Output for the Results QA
  }

  AliAnalysisTaskValeNanoTreeLPhi::~AliAnalysisTaskValeNanoTreeLPhi() {}
