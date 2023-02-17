#include "AliAnalysisTaskNanoXiPi.h"
#include "AliNanoAODTrack.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"

ClassImp(AliAnalysisTaskNanoXiPi)
    AliAnalysisTaskNanoXiPi::AliAnalysisTaskNanoXiPi()
    : AliAnalysisTaskSE(),
      fisLightWeight(false),
      fIsMC(false),
      fQA(nullptr),
      fEvent(nullptr),
      fEventCuts(nullptr),
      fEvtList(nullptr),
      fTrack(nullptr),
      fPion(nullptr),
      fPionList(nullptr),
      fPionMCList(nullptr),
      fAntiPion(nullptr),
      fAntiPionList(nullptr),
      fAntiPionMCList(nullptr),
      fCascade(nullptr),
      fXi(nullptr),
      fXiList(nullptr),
      fXiMCList(nullptr),
      fAntiXi(nullptr),
      fAntiXiList(nullptr),
      fAntiXiMCList(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fResults(nullptr),
      fResultsQA(nullptr),
      fSample(nullptr),
      fResultsSample(nullptr),
      fResultsSampleQA(nullptr),
      fTrackBufferSize(2000),
      fGTI(nullptr)
{
}

AliAnalysisTaskNanoXiPi::AliAnalysisTaskNanoXiPi(const char *name, bool isMC)
    : AliAnalysisTaskSE(name),
      fisLightWeight(false),
      fIsMC(false),
      fQA(nullptr),
      fEvent(nullptr),
      fEventCuts(nullptr),
      fEvtList(nullptr),
      fTrack(nullptr),
      fPion(nullptr),
      fPionList(nullptr),
      fPionMCList(nullptr),
      fAntiPion(nullptr),
      fAntiPionList(nullptr),
      fAntiPionMCList(nullptr),
      fCascade(nullptr),
      fXi(nullptr),
      fXiList(nullptr),
      fXiMCList(nullptr),
      fAntiXi(nullptr),
      fAntiXiList(nullptr),
      fAntiXiMCList(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fResults(nullptr),
      fResultsQA(nullptr),
      fSample(nullptr),
      fResultsSample(nullptr),
      fResultsSampleQA(nullptr),
      fTrackBufferSize(2000),
      fGTI(nullptr)
{
  DefineOutput(1, TList::Class());  // Output for the Event Class and Pair Cleaner
  DefineOutput(2, TList::Class());  // Output for the Event Cuts
  DefineOutput(3, TList::Class());  // Output for the Pion Cuts
  DefineOutput(4, TList::Class());  // Output for the AntiPion Cuts
  DefineOutput(5, TList::Class());  // Output for the Cascade Cuts
  DefineOutput(6, TList::Class());  // Output for the AntiCascade Cuts
  DefineOutput(7, TList::Class());  // Output for the Results
  DefineOutput(8, TList::Class());  // Output for the Results QA
  DefineOutput(9, TList::Class());  // Output for the Results Sample
  DefineOutput(10, TList::Class()); // Output for the Results Sample QA
  if (isMC)
  {
    DefineOutput(11, TList::Class()); // Output for the Track MC
    DefineOutput(12, TList::Class()); // Output for the Anti Track MC
    DefineOutput(13, TList::Class()); // Output for the Xi MC
    DefineOutput(14, TList::Class()); // Output for the Anti Xi MC
  }
}

AliAnalysisTaskNanoXiPi::~AliAnalysisTaskNanoXiPi()
{
  if (fEvent)
  {
    delete fEvent;
  }
  if (fEventCuts)
  {
    delete fEventCuts;
  }
  if (fTrack)
  {
    delete fTrack;
  }
  if (fPion)
  {
    delete fPion;
  }
  if (fAntiPion)
  {
    delete fAntiPion;
  }
  if (fCascade)
  {
    delete fCascade;
  }
  if (fXi)
  {
    delete fXi;
  }
  if (fAntiXi)
  {
    delete fAntiXi;
  }
  if (fPairCleaner)
  {
    delete fPairCleaner;
  }
  if (fPartColl)
  {
    delete fPartColl;
  }
  if (fSample)
  {
    delete fSample;
  }
}

void AliAnalysisTaskNanoXiPi::UserCreateOutputObjects()
{
  fGTI = new AliVTrack *[fTrackBufferSize];

  if (!fEventCuts)
  {
    AliError("No Event cuts \n");
  }
  else
  {
    fEventCuts->InitQA();
  }
  if (!fPion)
  {
    AliError("No Pion cuts \n");
  }
  else
  {
    fPion->Init();
  }
  if (!fAntiPion)
  {
    AliError("No AntiPion cuts \n");
  }
  else
  {
    fAntiPion->Init();
  }
  if (!fXi)
  {
    AliError("No Xi cuts \n");
  }
  else
  {
    fXi->Init();
  }
  if (!fAntiXi)
  {
    AliError("No AntiXi cuts \n");
  }
  else
  {
    fAntiXi->Init();
  }
  if (!fConfig)
  {
    AliError("No Correlation Config \n");
  }
  else
  {
    fPartColl = new AliFemtoDreamPartCollection(fConfig,
                                                fConfig->GetMinimalBookingME());
    // Xi+ - π- , Xi-, π+
    fPairCleaner = new AliFemtoDreamPairCleaner(4, 3,
                                                fConfig->GetMinimalBookingME());
  }

  // Event
  fEvent = new AliFemtoDreamEvent(true, !fisLightWeight,
                                  GetCollisionCandidates(), true);
  fEvent->SetMultiplicityEstimator(fConfig->GetMultiplicityEstimator());

  // Tracks
  fTrack = new AliFemtoDreamTrack();
  fTrack->SetUseMCInfo(
      fPion->GetIsMonteCarlo() || fAntiPion->GetIsMonteCarlo());

  // Xi
  fCascade = new AliFemtoDreamCascade();
  fCascade->SetUseMCInfo(fXi->GetIsMonteCarlo() || fAntiXi->GetIsMonteCarlo());
  // PDG Codes should be set assuming Xi- to also work for Xi+
  fCascade->SetPDGCode(3312);
  fCascade->SetPDGDaugPos(2212);
  fCascade->GetPosDaug()->SetUseMCInfo(
      fXi->GetIsMonteCarlo() || fAntiXi->GetIsMonteCarlo());
  fCascade->SetPDGDaugNeg(211);
  fCascade->GetNegDaug()->SetUseMCInfo(
      fXi->GetIsMonteCarlo() || fAntiXi->GetIsMonteCarlo());
  fCascade->SetPDGDaugBach(211);
  fCascade->GetBach()->SetUseMCInfo(
      fXi->GetIsMonteCarlo() || fAntiXi->GetIsMonteCarlo());
  fCascade->Setv0PDGCode(3122);

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

  fPionList = fPion->GetQAHists();
  fAntiPionList = fAntiPion->GetQAHists();

  fXiList = fXi->GetQAHists();
  fAntiXiList = fAntiXi->GetQAHists();

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

  fResultsSampleQA = new TList();
  fResultsSampleQA->SetOwner();
  fResultsSampleQA->SetName("ResultsSampleQA");

  fResultsSample = new TList();
  fResultsSample->SetOwner();
  fResultsSample->SetName("ResultsSample");

  PostData(1, fQA);
  PostData(2, fEvtList);
  PostData(3, fPionList);
  PostData(4, fAntiPionList);
  PostData(5, fXiList);
  PostData(6, fAntiXiList);
  PostData(7, fResults);
  PostData(8, fResultsQA);
  PostData(9, fResultsSample);
  PostData(10, fResultsSampleQA);

  if (fPion->GetIsMonteCarlo())
  {
    if (!fPion->GetMinimalBooking())
    {
      fPionMCList = fPion->GetMCQAHists();
    }
    else
    {
      fPionMCList = new TList();
      fPionMCList->SetName("MCTrkCuts");
      fPionMCList->SetOwner();
    }
    PostData(11, fPionMCList);
  }
  if (fAntiPion->GetIsMonteCarlo())
  {
    if (!fAntiPion->GetMinimalBooking())
    {
      fAntiPionMCList = fAntiPion->GetMCQAHists();
    }
    else
    {
      fAntiPionMCList = new TList();
      fAntiPionMCList->SetName("MCAntiTrkCuts");
      fAntiPionMCList->SetOwner();
    }
    PostData(12, fAntiPionMCList);
  }
  if (fXi->GetIsMonteCarlo())
  {
    if (!fXi->GetMinimalBooking())
    {
      fXiMCList = fXi->GetMCQAHists();
    }
    else
    {
      fXiMCList = new TList();
      fXiMCList->SetName("MCXiCuts");
      fXiMCList->SetOwner();
    }
    PostData(13, fXiMCList);
  }
  if (fAntiXi->GetIsMonteCarlo())
  {
    if (!fAntiXi->GetMinimalBooking())
    {
      fAntiXiMCList = fAntiXi->GetMCQAHists();
    }
    else
    {
      fAntiXiMCList = new TList();
      fAntiXiMCList->SetName("MCAntiv0Cuts");
      fAntiXiMCList->SetOwner();
    }
    PostData(14, fAntiXiMCList);
  }
}

void AliAnalysisTaskNanoXiPi::UserExec(Option_t *option)
{
  //  AliVEvent *fInputEvent = InputEvent();
  if (!fInputEvent)
  {
    AliError("No input event");
    return;
  }
  fEvent->SetEvent(fInputEvent);
  if (!fEventCuts->isSelected(fEvent))
  {
    return;
  }

  // PION SELECTION
  ResetGlobalTrackReference();
  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack)
  {
    AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
    if (!track)
    {
      AliFatal("No Standard AOD");
      return;
    }
    StoreGlobalTrackReference(track);
  }
  std::vector<AliFemtoDreamBasePart> Pions;     //π+
  std::vector<AliFemtoDreamBasePart> AntiPions; //π-
  fTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack)
  {
    AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
    fTrack->SetTrack(track, fInputEvent);
    if (fPion->isSelected(fTrack))
    {
      Pions.push_back(*fTrack);
    }
    if (fAntiPion->isSelected(fTrack))
    {
      AntiPions.push_back(*fTrack);
    }
  }

  std::vector<AliFemtoDreamBasePart> Xis;
  std::vector<AliFemtoDreamBasePart> AntiXis;
  AliAODEvent *aodEvt = dynamic_cast<AliAODEvent *>(fInputEvent);

  for (int iCasc = 0;
       iCasc < static_cast<TClonesArray *>(aodEvt->GetCascades())->GetEntriesFast();
       ++iCasc)
  {
    AliAODcascade *casc = aodEvt->GetCascade(iCasc);
    fCascade->SetCascade(fInputEvent, casc);
    if (fXi->isSelected(fCascade))
    {
      Xis.push_back(*fCascade);
    }
    if (fAntiXi->isSelected(fCascade))
    {
      AntiXis.push_back(*fCascade);
    }
  }

  // loop once over the MC stack to calculate Efficiency/Purity
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
        if (mcPart->GetPdgCode() == fPion->GetPDGCode())
        {
          fPion->FillGenerated(mcPart->Pt());
        }
        else if (mcPart->GetPdgCode() == fAntiPion->GetPDGCode())
        {
          fAntiPion->FillGenerated(mcPart->Pt());
        }
        else if (mcPart->GetPdgCode() == fXi->GetPDGv0())
        {
          fXi->FillGenerated(mcPart->Pt());
        }
        else if (mcPart->GetPdgCode() == fAntiXi->GetPDGv0())
        {
          fAntiXi->FillGenerated(mcPart->Pt());
        }
      }
    }
  }

  fPairCleaner->ResetArray();
  fPairCleaner->CleanTrackAndDecay(&Pions, &AntiXis, 0, true);
  fPairCleaner->CleanTrackAndDecay(&AntiPions, &Xis, 1, true);
  fPairCleaner->CleanTrackAndDecay(&Pions, &Xis, 2, true);
  fPairCleaner->CleanTrackAndDecay(&AntiPions, &AntiXis, 3, true);

  fPairCleaner->CleanDecay(&Xis, 0);
  fPairCleaner->CleanDecay(&AntiXis, 1);
  fPairCleaner->CleanDecayAndDecay(&Xis, &AntiXis, 2);

  fPairCleaner->StoreParticle(Pions);
  fPairCleaner->StoreParticle(AntiPions);
  fPairCleaner->StoreParticle(Xis);
  fPairCleaner->StoreParticle(AntiXis);
  if (fPairCleaner->GetCounter() > 0)
  {
    if (fConfig->GetUseEventMixing())
    {
      fPartColl->SetEvent(fPairCleaner->GetCleanParticles(),
                          fEvent->GetZVertex(), fEvent->GetMultiplicity(),
                          fEvent->GetV0MCentrality());
    }
  }

  PostData(1, fQA);
  PostData(2, fEvtList);
  PostData(3, fPionList);
  PostData(4, fAntiPionList);
  PostData(5, fXiList);
  PostData(6, fAntiXiList);
  PostData(7, fResults);
  PostData(8, fResultsQA);
  PostData(9, fResultsSample);
  PostData(10, fResultsSampleQA);

  if (fPion->GetIsMonteCarlo())
  {
    PostData(11, fPionMCList);
  }
  if (fAntiPion->GetIsMonteCarlo())
  {
    PostData(12, fAntiPionMCList);
  }
  if (fXi->GetIsMonteCarlo())
  {
    PostData(13, fXiMCList);
  }
  if (fAntiXi->GetIsMonteCarlo())
  {
    PostData(14, fAntiXiMCList);
  }
}

//____________________________________________________________________________________________________
void AliAnalysisTaskNanoXiPi::ResetGlobalTrackReference()
{
  // see AliFemtoDreamAnalysis for details
  for (int i = 0; i < fTrackBufferSize; i++)
  {
    fGTI[i] = 0;
  }
}

//____________________________________________________________________________________________________
void AliAnalysisTaskNanoXiPi::StoreGlobalTrackReference(AliVTrack *track)
{
  // see AliFemtoDreamAnalysis for details
  AliNanoAODTrack *nanoTrack = dynamic_cast<AliNanoAODTrack *>(track);
  const int trackID = track->GetID();
  if (trackID < 0)
  {
    return;
  }
  if (trackID >= fTrackBufferSize)
  {
    printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n", trackID,
           fTrackBufferSize);
    return;
  }

  if (fGTI[trackID])
  {
    if ((!nanoTrack->GetFilterMap()) && (!track->GetTPCNcls()))
    {
      return;
    }
    if (dynamic_cast<AliNanoAODTrack *>(fGTI[trackID])->GetFilterMap() || fGTI[trackID]->GetTPCNcls())
    {
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
