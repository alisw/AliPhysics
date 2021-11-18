/*
 * AliAnalysisTaskNanoLK.cxx
 *
 *  Created on: May 13, 2019
 *      Author: schmollweger
 */
#include "AliAnalysisTaskNanoLK.h"
#include "AliNanoAODTrack.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"

ClassImp(AliAnalysisTaskNanoLK)
    AliAnalysisTaskNanoLK::AliAnalysisTaskNanoLK()
    : AliAnalysisTaskSE(),
      fisLightWeight(false),
      fIsMC(false),
      fQA(nullptr),
      fEvent(nullptr),
      fEventCuts(nullptr),
      fEvtList(nullptr),
      fTrack(nullptr),
      fPosKaonCuts(nullptr),
      fPosKaonList(nullptr),
      fPosKaonMCList(nullptr),
      fNegKaonCuts(nullptr),
      fNegKaonList(nullptr),
      fNegKaonMCList(nullptr),
      fv0(nullptr),
      fLambdaCuts(nullptr),
      fLambdaList(nullptr),
      fLambdaMCList(nullptr),
      fAntiLambdaCuts(nullptr),
      fAntiLambdaList(nullptr),
      fAntiLambdaMCList(nullptr),
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

AliAnalysisTaskNanoLK::AliAnalysisTaskNanoLK(const char *name, bool isMC)
    : AliAnalysisTaskSE(name),
      fisLightWeight(false),
      fIsMC(false),
      fQA(nullptr),
      fEvent(nullptr),
      fEventCuts(nullptr),
      fEvtList(nullptr),
      fTrack(nullptr),
      fPosKaonCuts(nullptr),
      fPosKaonList(nullptr),
      fPosKaonMCList(nullptr),
      fNegKaonCuts(nullptr),
      fNegKaonList(nullptr),
      fNegKaonMCList(nullptr),
      fv0(nullptr),
      fLambdaCuts(nullptr),
      fLambdaList(nullptr),
      fLambdaMCList(nullptr),
      fAntiLambdaCuts(nullptr),
      fAntiLambdaList(nullptr),
      fAntiLambdaMCList(nullptr),
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
  DefineOutput(1, TList::Class());  //Output for QA of event and pair cleaner
  DefineOutput(2, TList::Class());  //Output for the Event Cuts
  DefineOutput(3, TList::Class());  //Output for the Positive Kaon Cuts
  DefineOutput(4, TList::Class());  //Output for the Negative Kaon Cuts
  DefineOutput(5, TList::Class());  //Output for the Lambda Cuts
  DefineOutput(6, TList::Class());  //Output for the AntiLambda Cuts
  DefineOutput(7, TList::Class());  //Output for the Results
  DefineOutput(8, TList::Class());  //Output for the Results QA
  DefineOutput(9, TList::Class());  //Output for the Results Sample
  DefineOutput(10, TList::Class()); //Output for the Results Sample QA
  if (isMC)
  {
    DefineOutput(11, TList::Class()); //Output for the Positive Kaon MC
    DefineOutput(12, TList::Class()); //Output for the Negative Kaon MC
    DefineOutput(13, TList::Class()); //Output for the Lambda MC
    DefineOutput(14, TList::Class()); //Output for the Anti Lambda MC
  }
}

///Destructor
AliAnalysisTaskNanoLK::~AliAnalysisTaskNanoLK()
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
  if (fPosKaonCuts)
  {
    delete fPosKaonCuts;
  }
  if (fNegKaonCuts)
  {
    delete fNegKaonCuts;
  }
  if (fv0)
  {
    delete fv0;
  }
  if (fLambdaCuts)
  {
    delete fLambdaCuts;
  }
  if (fAntiLambdaCuts)
  {
    delete fAntiLambdaCuts;
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

///Define the container outputs

void AliAnalysisTaskNanoLK::UserCreateOutputObjects()
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
  if (!fPosKaonCuts)
  {
    AliError("No Positive Kaon cuts \n");
  }
  else
  {
    fPosKaonCuts->Init();
  }
  if (!fNegKaonCuts)
  {
    AliError("No Negative Kaon cuts \n");
  }
  else
  {
    fNegKaonCuts->Init();
  }
  if (!fLambdaCuts)
  {
    AliError("No Lambda cuts \n");
  }
  else
  {
    fLambdaCuts->Init();
  }
  if (!fAntiLambdaCuts)
  {
    AliError("No AntiLambda cuts \n");
  }
  else
  {
    fAntiLambdaCuts->Init();
  }
  if (!fConfig)
  {
    AliError("No Correlation Config \n");
  }
  else
  {
    fPartColl = new AliFemtoDreamPartCollection(fConfig,
                                                fConfig->GetMinimalBookingME());
    //Arguments for the pair cleaner as follows:
    //1. How many pairs of Tracks + Decays do you want to clean?
    //(for the purpose of this tutorial, we are going to treat the
    //second track as a decay ;)
    //2. How many decays and decays do you want to clean
    //3. Minimal booking == true means no histograms are created and filled
    //might be handy for systematic checks, in order to reduce the memory
    //usage
    fPairCleaner = new AliFemtoDreamPairCleaner(4, 3,
                                                fConfig->GetMinimalBookingME());
    if (fConfig->GetUsePhiSpinning())
    {
      fSample = new AliFemtoDreamControlSample(fConfig);
    }
  }
  fEvent = new AliFemtoDreamEvent(true, !fisLightWeight,
                                  GetCollisionCandidates(), true);
  fEvent->SetMultiplicityEstimator(fConfig->GetMultiplicityEstimator());

  ///Setting up the Track object
  fTrack = new AliFemtoDreamTrack();
  fTrack->SetUseMCInfo(
      fPosKaonCuts->GetIsMonteCarlo() || fNegKaonCuts->GetIsMonteCarlo());
  ///Setting up the V0 object
  fv0 = new AliFemtoDreamv0();
  fv0->SetUseMCInfo(
      fLambdaCuts->GetIsMonteCarlo() || fAntiLambdaCuts->GetIsMonteCarlo());
  //PDG Codes should be set assuming Lambda to also work for AntiLambda
  fv0->SetPDGCode(3122);
  fv0->SetPDGDaughterPos(2212); //Λ->p π-
  fv0->GetPosDaughter()->SetUseMCInfo(
      fLambdaCuts->GetIsMonteCarlo() || fAntiLambdaCuts->GetIsMonteCarlo());
  fv0->SetPDGDaughterNeg(211);
  fv0->GetNegDaughter()->SetUseMCInfo(
      fLambdaCuts->GetIsMonteCarlo() || fAntiLambdaCuts->GetIsMonteCarlo());

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

  fPosKaonList = fPosKaonCuts->GetQAHists();
  fNegKaonList = fNegKaonCuts->GetQAHists();
  fLambdaList = fLambdaCuts->GetQAHists();
  fAntiLambdaList = fAntiLambdaCuts->GetQAHists();

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

  if (fConfig->GetUsePhiSpinning())
  {
    fResultsSample = fSample->GetHistList();

    if (!fConfig->GetMinimalBookingSample())
    {
      fResultsSampleQA->Add(fSample->GetQAList());
      fResultsQA->Add(fPairCleaner->GetHistList());
    }
  }
  else
  {
    fResultsSample = new TList();
    fResultsSample->SetOwner();
    fResultsSample->SetName("ResultsSample");
  }

  PostData(1, fQA);
  PostData(2, fEvtList);
  PostData(3, fPosKaonList);
  PostData(4, fNegKaonList);
  PostData(5, fLambdaList);
  PostData(6, fAntiLambdaList);
  PostData(7, fResults);
  PostData(8, fResultsQA);
  PostData(9, fResultsSample);
  PostData(10, fResultsSampleQA);

  if (fPosKaonCuts->GetIsMonteCarlo())
  {
    if (!fPosKaonCuts->GetMinimalBooking())
    {
      fPosKaonMCList = fPosKaonCuts->GetMCQAHists();
    }
    else
    {
      fPosKaonMCList = new TList();
      fPosKaonMCList->SetName("MCTrkCuts");
      fPosKaonMCList->SetOwner();
    }
    PostData(11, fPosKaonMCList);
  }
  if (fNegKaonCuts->GetIsMonteCarlo())
  {
    if (!fNegKaonCuts->GetMinimalBooking())
    {
      fNegKaonMCList = fNegKaonCuts->GetMCQAHists();
    }
    else
    {
      fNegKaonMCList = new TList();
      fNegKaonMCList->SetName("MCAntiTrkCuts");
      fNegKaonMCList->SetOwner();
    }
    PostData(12, fNegKaonMCList);
  }

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
    PostData(13, fLambdaMCList);
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
    PostData(14, fAntiLambdaMCList);
  }
}

void AliAnalysisTaskNanoLK::UserExec(Option_t *option)
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
  std::vector<AliFemtoDreamBasePart> PosKaons;
  std::vector<AliFemtoDreamBasePart> NegKaons;
  fTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack)
  {
    AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
    fTrack->SetTrack(track, fInputEvent);
    if (fPosKaonCuts->isSelected(fTrack))
    {
      PosKaons.push_back(*fTrack);
    }
    if (fNegKaonCuts->isSelected(fTrack))
    {
      NegKaons.push_back(*fTrack);
    }
  }

  std::vector<AliFemtoDreamBasePart> Lambdas;
  std::vector<AliFemtoDreamBasePart> AntiLambdas;
  AliAODEvent *aodEvt = dynamic_cast<AliAODEvent *>(fInputEvent);
  fv0->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  for (int iv0 = 0;
       iv0 < static_cast<TClonesArray *>(aodEvt->GetV0s())->GetEntriesFast();
       ++iv0)
  {
    AliAODv0 *casc = aodEvt->GetV0(iv0);
    fv0->Setv0(fInputEvent, casc);
    if (fLambdaCuts->isSelected(fv0))
    {
      Lambdas.push_back(*fv0);
    }
    if (fAntiLambdaCuts->isSelected(fv0))
    {
      AntiLambdas.push_back(*fv0);
    }
  }

  //loop once over the MC stack to calculate Efficiency/Purity
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
      }
    }
  }

  fPairCleaner->ResetArray();
  fPairCleaner->CleanTrackAndDecay(&PosKaons, &Lambdas, 0);
  fPairCleaner->CleanTrackAndDecay(&NegKaons, &AntiLambdas, 1);
  fPairCleaner->CleanTrackAndDecay(&PosKaons, &AntiLambdas, 2);
  fPairCleaner->CleanTrackAndDecay(&NegKaons, &Lambdas, 3);

  fPairCleaner->CleanDecay(&Lambdas, 0);
  fPairCleaner->CleanDecay(&AntiLambdas, 1);
  fPairCleaner->CleanDecayAndDecay(&Lambdas, &AntiLambdas, 2);

  fPairCleaner->StoreParticle(PosKaons);
  fPairCleaner->StoreParticle(NegKaons);
  fPairCleaner->StoreParticle(Lambdas);
  fPairCleaner->StoreParticle(AntiLambdas);

  if (fPairCleaner->GetCounter() > 0)
  {
    if (fConfig->GetUseEventMixing())
    {
      fPartColl->SetEvent(fPairCleaner->GetCleanParticles(),
                          fEvent->GetZVertex(), fEvent->GetMultiplicity(),
                          fEvent->GetV0MCentrality());
    }
    if (fConfig->GetUsePhiSpinning())
    {
      fSample->SetEvent(fPairCleaner->GetCleanParticles(), fEvent);
    }
  }

  PostData(1, fQA);
  PostData(2, fEvtList);
  PostData(3, fPosKaonList);
  PostData(4, fNegKaonList);
  PostData(5, fLambdaList);
  PostData(6, fAntiLambdaList);
  PostData(7, fResults);
  PostData(8, fResultsQA);
  PostData(9, fResultsSample);
  PostData(10, fResultsSampleQA);

  if (fPosKaonCuts->GetIsMonteCarlo())
  {
    PostData(11, fPosKaonMCList);
  }
  if (fNegKaonCuts->GetIsMonteCarlo())
  {
    PostData(12, fNegKaonMCList);
  }
  if (fLambdaCuts->GetIsMonteCarlo())
  {
    PostData(13, fLambdaMCList);
  }
  if (fAntiLambdaCuts->GetIsMonteCarlo())
  {
    PostData(14, fAntiLambdaMCList);
  }
}

//____________________________________________________________________________________________________
void AliAnalysisTaskNanoLK::ResetGlobalTrackReference()
{
  // see AliFemtoDreamAnalysis for details
  for (int i = 0; i < fTrackBufferSize; i++)
  {
    fGTI[i] = 0;
  }
}

//____________________________________________________________________________________________________
void AliAnalysisTaskNanoLK::StoreGlobalTrackReference(AliVTrack *track)
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
