/*
 * AliAnalysisTaskNanoLKr.cxx
 *
 *  Created on: 29 November, 2021
 *      Author: rossana
 */
#include "AliAnalysisTaskNanoLKr.h"

#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliNanoAODTrack.h"

/// constructor 1. We use ClassImp to "generate properly documentation"
ClassImp(AliAnalysisTaskNanoLKr) AliAnalysisTaskNanoLKr::AliAnalysisTaskNanoLKr()
    : AliAnalysisTaskSE(),
      fisLightWeight(false),
      fIsMC(false),
      fTrigger(AliVEvent::kINT7),
      fQA(nullptr),
      fEvent(nullptr),
      fEventCuts(nullptr),
      fEvtList(nullptr),
      fTrack(nullptr),
      fPosKaon(nullptr),
      fPosKaonList(nullptr),
      fPosKaonMCList(nullptr),
      fNegKaon(nullptr),
      fNegKaonList(nullptr),
      fNegKaonMCList(nullptr),
      fv0(nullptr),
      fLambda(nullptr),
      fLambdaList(nullptr),
      fLambdaMCList(nullptr),
      fAntiLambda(nullptr),
      fAntiLambdaList(nullptr),
      fAntiLambdaMCList(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fResults(nullptr),
      fResultsQA(nullptr),
      fSample(nullptr),
      fTrackBufferSize(2000),
      fGTI(nullptr)  /// technicality
{}

/// second constructor
AliAnalysisTaskNanoLKr::AliAnalysisTaskNanoLKr(const char *name, bool isMC)
    : AliAnalysisTaskSE(name),
      fisLightWeight(false),
      fIsMC(isMC),
      fTrigger(AliVEvent::kINT7),
      fQA(nullptr),
      fEvent(nullptr),
      fEventCuts(nullptr),
      fEvtList(nullptr),
      fTrack(nullptr),
      fPosKaon(nullptr),
      fPosKaonList(nullptr),
      fPosKaonMCList(nullptr),
      fNegKaon(nullptr),
      fNegKaonList(nullptr),
      fNegKaonMCList(nullptr),
      fv0(nullptr),
      fLambda(nullptr),
      fLambdaList(nullptr),
      fLambdaMCList(nullptr),
      fAntiLambda(nullptr),
      fAntiLambdaList(nullptr),
      fAntiLambdaMCList(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fResults(nullptr),
      fResultsQA(nullptr),
      fSample(nullptr),
      fTrackBufferSize(2000),
      fGTI(nullptr) {
    DefineOutput(1, TList::Class());  // Output for the Event Class and Pair Cleaner
    DefineOutput(2, TList::Class());  // Output for the Event Cuts
    DefineOutput(3, TList::Class());  // Output for the Positive Kaon Cuts
    DefineOutput(4, TList::Class());  // Output for the Negative Kaon Cuts
    DefineOutput(5, TList::Class());  // Output for the Lambda Cuts
    DefineOutput(6, TList::Class());  // Output for the AntiLambda Cuts
    DefineOutput(7, TList::Class());  // Output for the Results
    DefineOutput(8, TList::Class());  // Output for the Results QA

    if (isMC) {
        DefineOutput(9, TList::Class());   // Output for the Track MC
        DefineOutput(10, TList::Class());  // Output for the Anti Track MC
        DefineOutput(11, TList::Class());  // Output for the V0 MC
        DefineOutput(12, TList::Class());  // Output for the Anti V0 MC
    }
}

/// destructor (not necessary? vale didn t write it)
AliAnalysisTaskNanoLKr::~AliAnalysisTaskNanoLKr() {
    if (fEvent) {
        delete fEvent;
    }
    if (fEventCuts) {
        delete fEventCuts;
    }
    if (fTrack) {
        delete fTrack;
    }
    if (fPosKaon) {
        delete fPosKaon;
    }
    if (fNegKaon) {
        delete fNegKaon;
    }
    if (fv0) {
        delete fv0;
    }
    if (fLambda) {
        delete fLambda;
    }
    if (fAntiLambda) {
        delete fAntiLambda;
    }
    if (fPairCleaner) {
        delete fPairCleaner;
    }
    if (fPartColl) {
        delete fPartColl;
    }
    if (fSample) {
        delete fSample;
    }
}

/// AliAnalysisTaskNanoLambdaKaon::~AliAnalysisTaskNanoLambdaKaon() {} <- Vale put this

void AliAnalysisTaskNanoLKr::UserCreateOutputObjects() {
    fGTI = new AliVTrack *[fTrackBufferSize];  /// go through detail of this

    if (!fEventCuts) {
        AliError("No Event cuts \n");
    } else {
        fEventCuts->InitQA();  /// more details: go to header alifemtodreameventcut
    }
    if (!fPosKaon) {
        AliError("No PosKaon cuts \n");
    } else {
        fPosKaon->Init();
        fPosKaon->SetName("KaonPlus");
    }
    if (!fNegKaon) {
        AliError("No NegKaon cuts \n");
    } else {
        fNegKaon->Init();
        fNegKaon->SetName("KaonMinus");
    }
    if (!fLambda) {
        AliError("No Lambda cuts \n");
    } else {
        fLambda->Init();
        fLambda->SetName("Lambda");
    }
    if (!fAntiLambda) {
        AliError("No AntiLambda cuts \n");
    } else {
        fAntiLambda->Init();
        fAntiLambda->SetName("AntiLambda");
    }

    if (!fConfig) {
        AliError("No Correlation Config \n");
    } else {
        fPartColl = new AliFemtoDreamPartCollection(
            fConfig,
            fConfig->GetMinimalBookingME());  /// new: create a new object, allocating memory

        fPairCleaner = new AliFemtoDreamPairCleaner(4, 3, fConfig->GetMinimalBookingME());
    }

    fEvent = new AliFemtoDreamEvent(true, true, fTrigger);

    fTrack = new AliFemtoDreamTrack();

    fTrack->SetUseMCInfo(fPosKaon->GetIsMonteCarlo() || fNegKaon->GetIsMonteCarlo());

    fv0 = new AliFemtoDreamv0();
    fv0->SetUseMCInfo(fLambda->GetIsMonteCarlo() || fAntiLambda->GetIsMonteCarlo());
    // PDG Codes should be set assuming Lambda0 to also work for AntiLambda
    fv0->SetPDGCode(3122);  /// lambda
    fv0->SetPDGDaughterPos(2212);

    fv0->GetPosDaughter()->SetUseMCInfo(fLambda->GetIsMonteCarlo() || fAntiLambda->GetIsMonteCarlo());
    fv0->SetPDGDaughterNeg(211);
    fv0->GetNegDaughter()->SetUseMCInfo(fLambda->GetIsMonteCarlo() || fAntiLambda->GetIsMonteCarlo());

    fQA = new TList();
    fQA->SetOwner();
    fQA->SetName("QA");
    if (!fConfig->GetMinimalBookingME() && fEvent && fEvent->GetEvtCutList()) {
        fQA->Add(fEvent->GetEvtCutList());
    }

    if (!fEventCuts->GetMinimalBooking()) {
        fEvtList = fEventCuts->GetHistList();
    } else {
        fEvtList = new TList();
        fEvtList->SetName("EventCuts");
        fEvtList->SetOwner();
    }

    fPosKaonList = fPosKaon->GetQAHists();
    fNegKaonList = fNegKaon->GetQAHists();
    fLambdaList = fLambda->GetQAHists();
    fAntiLambdaList = fAntiLambda->GetQAHists();

    if (fPartColl && fPartColl->GetHistList()) {
        fResults = fPartColl->GetHistList();
    }

    if (!fConfig->GetMinimalBookingME() && fPartColl && fPartColl->GetQAList()) {
        fResultsQA = fPartColl->GetQAList();
    } else {
        fResultsQA = new TList();
        fResultsQA->SetName("ResultsQA");
        fResultsQA->SetOwner(true);
    }

    if (!fConfig->GetMinimalBookingME() && fPairCleaner && fPairCleaner->GetHistList()) {
        fQA->Add(fPairCleaner->GetHistList());
    }

    PostData(1, fQA);  /// Postdata: to fill outputs already defined
    PostData(2, fEvtList);
    PostData(3, fPosKaonList);
    PostData(4, fNegKaonList);
    PostData(5, fLambdaList);
    PostData(6, fAntiLambdaList);

    PostData(7, fResults);
    PostData(8, fResultsQA);

    if (fPosKaon->GetIsMonteCarlo()) {
        if (!fPosKaon->GetMinimalBooking()) {
            fPosKaonMCList = fPosKaon->GetMCQAHists();
        } else {
            fPosKaonMCList = new TList();
            fPosKaonMCList->SetName("MCPosKaonCuts");
            fPosKaonMCList->SetOwner();
        }
        PostData(9, fPosKaonMCList);
    }
    if (fNegKaon->GetIsMonteCarlo()) {
        if (!fNegKaon->GetMinimalBooking()) {
            fNegKaonMCList = fNegKaon->GetMCQAHists();
        } else {
            fNegKaonMCList = new TList();
            fNegKaonMCList->SetName("MCNegKaonCuts");
            fNegKaonMCList->SetOwner();
        }
        PostData(10, fNegKaonMCList);
    }

    if (fLambda->GetIsMonteCarlo()) {
        if (!fLambda->GetMinimalBooking()) {
            fLambdaMCList = fLambda->GetMCQAHists();
        } else {
            fLambdaMCList = new TList();
            fLambdaMCList->SetName("MCv0Cuts");
            fLambdaMCList->SetOwner();
        }
        PostData(11, fLambdaMCList);
    }
    if (fAntiLambda->GetIsMonteCarlo()) {
        if (!fAntiLambda->GetMinimalBooking()) {
            fAntiLambdaMCList = fAntiLambda->GetMCQAHists();
        } else {
            fAntiLambdaMCList = new TList();
            fAntiLambdaMCList->SetName("MCAntiv0Cuts");
            fAntiLambdaMCList->SetOwner();
        }
        PostData(12, fAntiLambdaMCList);
    }
}

void AliAnalysisTaskNanoLKr::UserExec(Option_t *option) {
    AliVEvent *fInputEvent = InputEvent();  // maybe take off this comment
    if (!fInputEvent) {
        AliError("No input event");
        return;
    }
    fEvent->SetEvent(fInputEvent);
    if (!fEventCuts->isSelected(fEvent)) {
        return;
    }

    // KAON SELECTION
    ResetGlobalTrackReference();
    for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
        AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
        if (!track) {
            AliFatal("No Standard AOD");
            return;
        }
        StoreGlobalTrackReference(track);
    }
    std::vector<AliFemtoDreamBasePart> PosKaons;
    std::vector<AliFemtoDreamBasePart> NegKaons;
    fTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
    for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
        AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
        fTrack->SetTrack(track, fInputEvent);
        if (fPosKaon->isSelected(fTrack)) {
            PosKaons.push_back(*fTrack);
        }
        if (fNegKaon->isSelected(fTrack)) {
            NegKaons.push_back(*fTrack);
        }
    }

    std::vector<AliFemtoDreamBasePart> Lambdas;
    std::vector<AliFemtoDreamBasePart> AntiLambdas;
    AliAODEvent *aodEvt = dynamic_cast<AliAODEvent *>(fInputEvent);
    fv0->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
    for (int iv0 = 0; iv0 < static_cast<TClonesArray *>(aodEvt->GetV0s())->GetEntriesFast(); ++iv0) {
        AliAODv0 *casc = aodEvt->GetV0(iv0);
        fv0->Setv0(fInputEvent, casc);
        if (fLambda->isSelected(fv0)) {
            Lambdas.push_back(*fv0);
        }
        if (fAntiLambda->isSelected(fv0)) {
            AntiLambdas.push_back(*fv0);
        }
    }

    // loop once over the MC stack to calculate Efficiency/Purity
    if (fIsMC) {
        AliAODInputHandler *eventHandler =
            dynamic_cast<AliAODInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
        AliMCEvent *fMC = eventHandler->MCEvent();

        for (int iPart = 0; iPart < (fMC->GetNumberOfTracks()); iPart++) {
            AliAODMCParticle *mcPart = (AliAODMCParticle *)fMC->GetTrack(iPart);
            if (TMath::Abs(mcPart->Eta()) < 0.8 && mcPart->IsPhysicalPrimary()) {
                if (mcPart->GetPdgCode() == fPosKaon->GetPDGCode()) {
                    fPosKaon->FillGenerated(mcPart->Pt());
                } else if (mcPart->GetPdgCode() == fNegKaon->GetPDGCode()) {
                    fNegKaon->FillGenerated(mcPart->Pt());
                } else if (mcPart->GetPdgCode() == fLambda->GetPDGv0()) {
                    fLambda->FillGenerated(mcPart->Pt());
                } else if (mcPart->GetPdgCode() == fAntiLambda->GetPDGv0()) {
                    fAntiLambda->FillGenerated(mcPart->Pt());
                }
            }
        }
    }

    fPairCleaner->ResetArray();
    fPairCleaner->CleanTrackAndDecay(&PosKaons, &Lambdas, 0);
    fPairCleaner->CleanTrackAndDecay(&PosKaons, &AntiLambdas, 1);
    fPairCleaner->CleanTrackAndDecay(&NegKaons, &Lambdas, 2);
    fPairCleaner->CleanTrackAndDecay(&NegKaons, &AntiLambdas, 3);

    fPairCleaner->CleanDecay(&Lambdas, 0);
    fPairCleaner->CleanDecay(&AntiLambdas, 1);
    fPairCleaner->CleanDecayAndDecay(&Lambdas, &AntiLambdas, 2);

    fPairCleaner->StoreParticle(PosKaons);
    fPairCleaner->StoreParticle(NegKaons);
    fPairCleaner->StoreParticle(Lambdas);
    fPairCleaner->StoreParticle(AntiLambdas);

    if (fPairCleaner->GetCounter() > 0) {
        if (fConfig->GetUseEventMixing()) {
            fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(), fEvent->GetMultiplicity(),
                                fEvent->GetV0MCentrality());
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

    if (fPosKaon->GetIsMonteCarlo()) {
        PostData(9, fPosKaonMCList);
    }
    if (fNegKaon->GetIsMonteCarlo()) {
        PostData(10, fNegKaonMCList);
    }
    if (fLambda->GetIsMonteCarlo()) {
        PostData(11, fLambdaMCList);
    }
    if (fAntiLambda->GetIsMonteCarlo()) {
        PostData(12, fAntiLambdaMCList);
    }
}

//____________________________________________________________________________________________________
void AliAnalysisTaskNanoLKr::ResetGlobalTrackReference() {
    // see AliFemtoDreamAnalysis for details
    for (int i = 0; i < fTrackBufferSize; i++) {
        fGTI[i] = 0;
    }
}

//____________________________________________________________________________________________________
void AliAnalysisTaskNanoLKr::StoreGlobalTrackReference(AliVTrack *track) {
    // see AliFemtoDreamAnalysis for details
    AliNanoAODTrack *nanoTrack = dynamic_cast<AliNanoAODTrack *>(track);
    const int trackID = track->GetID();
    if (trackID < 0) {
        return;
    }
    if (trackID >= fTrackBufferSize) {
        printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n", trackID, fTrackBufferSize);
        return;
    }

    if (fGTI[trackID]) {
        if ((!nanoTrack->GetFilterMap()) && (!track->GetTPCNcls())) {
            return;
        }
        if (dynamic_cast<AliNanoAODTrack *>(fGTI[trackID])->GetFilterMap() || fGTI[trackID]->GetTPCNcls()) {
            printf("Warning! global track info already there!");
            printf("         TPCNcls track1 %u track2 %u", (fGTI[trackID])->GetTPCNcls(), track->GetTPCNcls());
            printf("         FilterMap track1 %u track2 %u\n",
                   dynamic_cast<AliNanoAODTrack *>(fGTI[trackID])->GetFilterMap(), nanoTrack->GetFilterMap());
        }
    }
    (fGTI[trackID]) = track;
}
