#include "AliAnalysisTaskNanoAODFemtoDreamLambdaPhi.h"
#include "AliFemtoDreamBasePart.h"
#include "AliLog.h"
#include "AliNanoAODTrack.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "TH1F.h"
#include "TList.h"
#include "AliAODMCParticle.h"

ClassImp(AliAnalysisTaskNanoAODFemtoDreamLambdaPhi)
    AliAnalysisTaskNanoAODFemtoDreamLambdaPhi::AliAnalysisTaskNanoAODFemtoDreamLambdaPhi()
    : AliAnalysisTaskSE(),
      fIsMC(false),
      fUseOMixing(false),
      fInvMassCutSBdown(0.0),
      fInvMassCutSBup(0.0),
      fTrigger(AliVEvent::kINT7),
      fQA(nullptr),
      fEvtList(nullptr),
      fLambdaList(nullptr),
      fLambdaMCList(nullptr),
      fAntiLambdaList(nullptr),
      fAntiLambdaMCList(nullptr),
      fKaonPlusList(nullptr),
      fKaonPlusMCList(nullptr),
      fKaonMinusList(nullptr),
      fKaonMinusMCList(nullptr),
      fPhiList(nullptr),
      fPhiMCList(nullptr),
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

AliAnalysisTaskNanoAODFemtoDreamLambdaPhi::AliAnalysisTaskNanoAODFemtoDreamLambdaPhi(
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
      fLambdaMCList(nullptr),
      fAntiLambdaList(nullptr),
      fAntiLambdaMCList(nullptr),
      fKaonPlusList(nullptr),
      fKaonPlusMCList(nullptr),
      fKaonMinusList(nullptr),
      fKaonMinusMCList(nullptr),
      fPhiList(nullptr),
      fPhiMCList(nullptr),
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

  if(isMC) {
    DefineOutput(10, TList::Class()); //Output for the Lambda MC
    DefineOutput(11, TList::Class()); //Output for the AntiLambda MC
    DefineOutput(12, TList::Class()); //Output for the KaonPlus MC
    DefineOutput(13, TList::Class()); //Output for the KaonMinus MC
    DefineOutput(14, TList::Class()); //Output for the Phi MC
  }
  }

  AliAnalysisTaskNanoAODFemtoDreamLambdaPhi::~AliAnalysisTaskNanoAODFemtoDreamLambdaPhi() {}

  void AliAnalysisTaskNanoAODFemtoDreamLambdaPhi::UserCreateOutputObjects()
  {

    fEvent = new AliFemtoDreamEvent(true, true, fTrigger);

    fTrack = new AliFemtoDreamTrack();
    fTrack->SetUseMCInfo(
        fPosKaonCuts->GetIsMonteCarlo() || fNegKaonCuts->GetIsMonteCarlo());

    fLambda = new AliFemtoDreamv0();
    fLambda->SetUseMCInfo(
        fLambdaCuts->GetIsMonteCarlo() || fAntiLambdaCuts->GetIsMonteCarlo());
    fLambda->SetPDGCode(3122);
    fLambda->SetPDGDaughterPos(2212);
    fLambda->GetPosDaughter()->SetUseMCInfo(
        fLambdaCuts->GetIsMonteCarlo() || fAntiLambdaCuts->GetIsMonteCarlo());
    fLambda->SetPDGDaughterNeg(211);
    fLambda->GetNegDaughter()->SetUseMCInfo(
        fLambdaCuts->GetIsMonteCarlo() || fAntiLambdaCuts->GetIsMonteCarlo());

    fPhiParticle = new AliFemtoDreamv0();
    fPhiParticle->SetPDGCode(fPhiCuts->GetPDGv0());
    fPhiParticle->SetUseMCInfo(fIsMC);
    fPhiParticle->SetPDGDaughterPos(
        fPhiCuts->GetPDGPosDaug()); // order +sign doesnt play a role
    fPhiParticle->GetPosDaughter()->SetUseMCInfo(fIsMC);
    fPhiParticle->SetPDGDaughterNeg(
        fPhiCuts->GetPDGNegDaug()); // only used for MC Matching
    fPhiParticle->GetNegDaughter()->SetUseMCInfo(fIsMC);

    fGTI = new AliVTrack *[fTrackBufferSize];

    if (!fEventCuts)
    {
      AliFatal("Event Cuts not set!");
    }
    fEventCuts->InitQA();

    fQA = new TList();
    fQA->SetOwner();
    fQA->SetName("QA");
    fQA->Add(fEvent->GetEvtCutList());

    if (!fEventCuts->GetMinimalBooking())
    {
      fEvtList = fEventCuts->GetHistList();
    }
    else
    {
      fEvtList = new TList();
      fEvtList->SetName("EventCuts");
      fEvtList->SetOwner();
    }

    if (!fPosKaonCuts)
    {
      AliFatal("Track Cuts for Particle One not set!");
    }
    fPosKaonCuts->Init();
    fPosKaonCuts->SetName("KaonPlus");

    // if (fPosKaonCuts->GetIsMonteCarlo())
    // {
    //   fPosKaonCuts->SetMCName("MCParticle1");
    //   fKaonPlusList->Add(fPosKaonCuts->GetMCQAHists());
    // }

    if (!fNegKaonCuts)
    {
      AliFatal("Track Cuts for Particle One not set!");
    }
    fNegKaonCuts->Init();
    fNegKaonCuts->SetName("KaonMinus");

    // if (fNegKaonCuts->GetIsMonteCarlo())
    // {
    //   fNegKaonCuts->SetMCName("MCParticle2");
    //   fKaonMinusList->Add(fNegKaonCuts->GetMCQAHists());
    // }

    if (!fPhiCuts)
    {
      AliFatal("Cuts for the phi not set!");
    }
    fPhiCuts->Init();
    fPhiCuts->SetName("Phi");

    // if (fPhiCuts->GetIsMonteCarlo())
    // {
    //   fPhiCuts->SetMCName("MCPhi");
    //   fPhiList->Add(fPhiCuts->GetMCQAHists());
    // }

    if (!fLambdaCuts)
    {
      AliFatal("Track Cuts for Lambda not set!");
    }
    fLambdaCuts->Init();
    fLambdaCuts->SetName("Lambda");

    // if (fLambdaCuts->GetIsMonteCarlo())
    // {
    //   fLambdaCuts->SetMCName("MCLambda");
    //   fLambdaList->Add(fLambdaCuts->GetMCQAHists());
    // }

    if (!fAntiLambdaCuts)
    {
      AliFatal("Track Cuts for AntiLambda not set!");
    }
    fAntiLambdaCuts->Init();
    fAntiLambdaCuts->SetName("AntiLambda");

    // if (fAntiLambdaCuts->GetIsMonteCarlo())
    // {
    //   fAntiLambdaCuts->SetMCName("MCAntiProton");
    //   fAntiLambdaList->Add(fAntiLambdaCuts->GetMCQAHists());
    // }

    fKaonPlusList = fPosKaonCuts->GetQAHists();
    fKaonMinusList = fNegKaonCuts->GetQAHists();
    fLambdaList = fLambdaCuts->GetQAHists();
    fAntiLambdaList = fAntiLambdaCuts->GetQAHists();
    fPhiList = fPhiCuts->GetQAHists();

    fPairCleaner =
        new AliFemtoDreamPairCleaner(0, 2, fConfig->GetMinimalBookingME());
    fPartColl =
        new AliFemtoDreamPartCollection(fConfig, fConfig->GetMinimalBookingME());

    fResultsQA = new TList();
    fResultsQA->SetOwner();
    fResultsQA->SetName("ResultsQA");
    if (fConfig->GetUseEventMixing())
    {
      fResults = fPartColl->GetHistList();
      if (!fConfig->GetMinimalBookingME())
      {
        fResultsQA->Add(fPartColl->GetQAList());
        fResultsQA->Add(fPairCleaner->GetHistList());
      }
    }
    else
    {
      fResults = new TList();
      fResults->SetOwner();
      fResults->SetName("Results");
    }

    PostData(1, fQA);
    PostData(2, fEvtList);
    PostData(3, fLambdaList);
    PostData(4, fAntiLambdaList);
    PostData(5, fKaonPlusList);
    PostData(6, fKaonMinusList);
    PostData(7, fPhiList);
    PostData(8, fResults);
    PostData(9, fResultsQA);

    if (fLambdaCuts->GetIsMonteCarlo())
    {
      if (!fLambdaCuts->GetMinimalBooking())
      {
        fLambdaMCList = fLambdaCuts->GetMCQAHists();
      }
      else
      {
        fLambdaMCList = new TList();
        fLambdaMCList->SetName("MCv0Cuts");
        fLambdaMCList->SetOwner();
      }
      PostData(10, fLambdaMCList);
    }
    if (fAntiLambdaCuts->GetIsMonteCarlo())
    {
      if (!fAntiLambdaCuts->GetMinimalBooking())
      {
        fAntiLambdaMCList = fAntiLambdaCuts->GetMCQAHists();
      }
      else
      {
        fAntiLambdaMCList = new TList();
        fAntiLambdaMCList->SetName("MCAntiv0Cuts");
        fAntiLambdaMCList->SetOwner();
      }
      PostData(11, fAntiLambdaMCList);
    }

    if (fPosKaonCuts->GetIsMonteCarlo())
    {
      if (!fPosKaonCuts->GetMinimalBooking())
      {
        fKaonPlusMCList = fPosKaonCuts->GetMCQAHists();
      }
      else
      {
        fKaonPlusMCList = new TList();
        fKaonPlusMCList->SetName("MCPosKaonCuts");
        fKaonPlusMCList->SetOwner();
      }
      PostData(12, fKaonPlusMCList);
    }

    if (fNegKaonCuts->GetIsMonteCarlo())
    {
      if (!fNegKaonCuts->GetMinimalBooking())
      {
        fKaonMinusMCList = fNegKaonCuts->GetMCQAHists();
      }
      else
      {
        fKaonMinusMCList = new TList();
        fKaonMinusMCList->SetName("MCNegKaonCuts");
        fKaonMinusMCList->SetOwner();
      }
      PostData(13, fKaonMinusMCList);
    }

    if (fPhiCuts->GetIsMonteCarlo())
    {
      if (!fPhiCuts->GetMinimalBooking())
      {
        fPhiMCList = fPhiCuts->GetMCQAHists();
      }
      else
      {
        fPhiMCList = new TList();
        fPhiMCList->SetName("MCPhiCuts");
        fPhiMCList->SetOwner();
      }
      PostData(14, fPhiMCList);
    }
  }

void AliAnalysisTaskNanoAODFemtoDreamLambdaPhi::UserExec(Option_t *) {
  AliVEvent *fInputEvent = InputEvent();

  // PREAMBLE - CHECK EVERYTHING IS THERE
  if (!fInputEvent) {
    AliError("No Input event");
    return;
  }

  // EVENT SELECTION
  fEvent->SetEvent(fInputEvent);
  if (!fEventCuts->isSelected(fEvent)) return;

  ResetGlobalTrackReference();
  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
    AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
    if (!track) {
      AliFatal("No Standard AOD");
      return;
    }
    StoreGlobalTrackReference(track);
  }
  fTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);

  std::vector<AliFemtoDreamBasePart> KaonPlus;
  std::vector<AliFemtoDreamBasePart> KaonMinus;
  std::vector<AliFemtoDreamBasePart> PhiParticles;
  std::vector<AliFemtoDreamBasePart> Lambdas;
  std::vector<AliFemtoDreamBasePart> AntiLambdas;
  AliAODEvent *aodEvt = dynamic_cast<AliAODEvent *>(fInputEvent);
  fLambda->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  for (int iv0 = 0;
       iv0 < static_cast<TClonesArray *>(aodEvt->GetV0s())->GetEntriesFast();
       ++iv0)
  {
    AliAODv0 *casc = aodEvt->GetV0(iv0);
    fLambda->Setv0(fInputEvent, casc);
    if (fLambdaCuts->isSelected(fLambda))
    {
      Lambdas.push_back(*fLambda);
    }
    if (fAntiLambdaCuts->isSelected(fLambda))
    {
      AntiLambdas.push_back(*fLambda);
    }
  }

  static float massKaon =
      TDatabasePDG::Instance()->GetParticle(fPosKaonCuts->GetPDGCode())->Mass();

  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
    AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));

    if (!track) continue;
    fTrack->SetTrack(track, fInputEvent);
    fTrack->SetInvMass(massKaon);
    if (fPosKaonCuts->isSelected(fTrack)) {
      KaonPlus.push_back(*fTrack);
    }
    if (fNegKaonCuts->isSelected(fTrack)) {
      KaonMinus.push_back(*fTrack);
    }
  }

  fPhiParticle->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  for (const auto &posK : KaonPlus) {
    for (const auto &negK : KaonMinus) {
      fPhiParticle->Setv0(posK, negK, fInputEvent, false, false, true);
      fPhiParticle->SetParticleOrigin(AliFemtoDreamBasePart::kPhysPrimary);

      if (fPhiCuts->isSelected(fPhiParticle)) {
        fPhiParticle->SetCPA(
            gRandom->Uniform());  // cpacode needed for CleanDecay v0;
        PhiParticles.push_back(*fPhiParticle);
      }
    }
  }

  if (fIsMC)
  {
    AliAODInputHandler *eventHandler =
        dynamic_cast<AliAODInputHandler *>(AliAnalysisManager::GetAnalysisManager()
                                               ->GetInputEventHandler());
    AliMCEvent *fMC = eventHandler->MCEvent();

    for (int iPart = 0; iPart < (fMC->GetNumberOfTracks()); iPart++)
    {
      AliAODMCParticle *mcPart = (AliAODMCParticle *)fMC->GetTrack(iPart);
      if (TMath::Abs(mcPart->Eta()) < 0.8 && mcPart->IsPhysicalPrimary())
      {
        if (mcPart->GetPdgCode() == fPosKaonCuts->GetPDGCode())
        {
          fPosKaonCuts->FillGenerated(mcPart->Pt());
        }
        else if (mcPart->GetPdgCode() == fNegKaonCuts->GetPDGCode())
        {
          fNegKaonCuts->FillGenerated(mcPart->Pt());
        }
        else if (mcPart->GetPdgCode() == fLambdaCuts->GetPDGv0())
        {
          fLambdaCuts->FillGenerated(mcPart->Pt());
        }
        else if (mcPart->GetPdgCode() == fAntiLambdaCuts->GetPDGv0())
        {
          fAntiLambdaCuts->FillGenerated(mcPart->Pt());
        }
        else if (mcPart->GetPdgCode() == fPhiCuts->GetPDGv0())
        {
          fPhiCuts->FillGenerated(mcPart->Pt());
        }
      }
    }
  }

  fPairCleaner->CleanDecayAndDecay(&Lambdas, &PhiParticles, 0);
  fPairCleaner->CleanDecayAndDecay(&AntiLambdas, &PhiParticles, 1);
  fPairCleaner->ResetArray();
  fPairCleaner->StoreParticle(Lambdas);         // 0
  fPairCleaner->StoreParticle(AntiLambdas);     // 1
  fPairCleaner->StoreParticle(PhiParticles);     // 2
  fPairCleaner->StoreParticle(KaonPlus);  // 3
  fPairCleaner->StoreParticle(KaonMinus);  // 4

  if (fPairCleaner->GetCounter() > 0) {
        fPartColl->SetEvent(fPairCleaner->GetCleanParticles(),
                            fEvent->GetZVertex(), fEvent->GetRefMult08(),
                            fEvent->GetV0MCentrality());
    }

    PostData(1, fQA);
    PostData(2, fEvtList);
    PostData(3, fLambdaList);
    PostData(4, fAntiLambdaList);
    PostData(5, fKaonPlusList);
    PostData(6, fKaonMinusList);
    PostData(7, fPhiList);
    PostData(8, fResults);
    PostData(9, fResultsQA);

    if (fLambdaCuts->GetIsMonteCarlo())
    {
      PostData(10, fLambdaMCList);
    }
    if (fAntiLambdaCuts->GetIsMonteCarlo())
    {
      PostData(11, fAntiLambdaMCList);
    }
    if (fPosKaonCuts->GetIsMonteCarlo())
    {
      PostData(12, fKaonPlusMCList);
    }
    if (fNegKaonCuts->GetIsMonteCarlo())
    {
      PostData(13, fKaonMinusMCList);
    }
    if (fPhiCuts->GetIsMonteCarlo())
    {
      PostData(14, fPhiMCList);
    }
}

void AliAnalysisTaskNanoAODFemtoDreamLambdaPhi::ResetGlobalTrackReference() {
  for (UShort_t i = 0; i < fTrackBufferSize; i++) {
    fGTI[i] = 0;
  }
}

void AliAnalysisTaskNanoAODFemtoDreamLambdaPhi::StoreGlobalTrackReference(
    AliVTrack *track) {
  // see AliFemtoDreamAnalysis for details
  AliNanoAODTrack *nanoTrack = dynamic_cast<AliNanoAODTrack *>(track);
  const int trackID = track->GetID();
  if (trackID < 0) {
    return;
  }
  if (trackID >= fTrackBufferSize) {
    printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n", trackID,
           fTrackBufferSize);
    return;
  }

  if (fGTI[trackID]) {
    if ((!nanoTrack->GetFilterMap()) && (!track->GetTPCNcls())) {
      return;
    }
    if (dynamic_cast<AliNanoAODTrack *>(fGTI[trackID])->GetFilterMap() ||
        fGTI[trackID]->GetTPCNcls()) {
      printf("Warning! global track info already there!");
      printf("         TPCNcls track1 %u track2 %u",
             (fGTI[trackID])->GetTPCNcls(), track->GetTPCNcls());
      printf("         FilterMap track1 %u track2 %u\n",
             dynamic_cast<AliNanoAODTrack *>(fGTI[trackID])->GetFilterMap(),
             nanoTrack->GetFilterMap());
    }
  }
  (fGTI[trackID]) = track;
}
