#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEvent.h"
#include "AliMCSpectraWeights.h"
#include "AliMultSelection.h"
#include "AliPhysicsSelection.h"
#include "AliStack.h"
#include "AliVEvent.h"
#include "AliVMultiplicity.h"
#include "TArrayD.h"
#include "TChain.h"
#include "TList.h"
#include "TString.h"
#include <iostream>

#include "AliMCSpectraWeightsAnalysisTask.h"

class AliMCSpectraWeightsAnalysisTask;

ClassImp(AliMCSpectraWeightsAnalysisTask);

AliMCSpectraWeightsAnalysisTask::AliMCSpectraWeightsAnalysisTask()
    : AliAnalysisTaskSE(), // default constructor
      fDebugLevel(0), fOutputList(0), fEvent(0), fMCEvent(0), fMCStack(0),
      fIsESD(1), fIsMC(1), fUseMultiplicity(kTRUE),
      fTriggerMask(AliVEvent::kMB), fMCSpectraWeights(0),
      // fstCollisionSystem("pp"),
      // fstMCTrainOutput(""),
      fHistMCPartCorr(0), fHistMCGenPrimTrack(0), fHistMCFractions(0),
      fHistDataFractions(0), fHistMCWeights(0), fBinsMultCent(0), fBinsPt(0),
      fBinsEta(0), fBinsZv(0) {
    // this is used by root for IO purposes, it needs to remain empty
}

AliMCSpectraWeightsAnalysisTask::AliMCSpectraWeightsAnalysisTask(
    const char* name)
    : AliAnalysisTaskSE(name), fDebugLevel(0), fOutputList(0), fEvent(0),
      fMCEvent(0), fMCStack(0), fIsESD(1), fIsMC(1), fUseMultiplicity(kTRUE),
      fTriggerMask(AliVEvent::kMB), fMCSpectraWeights(0),
      // fstCollisionSystem("pp"),
      // fstMCTrainOutput(""),
      fHistMCPartCorr(0), fHistMCGenPrimTrack(0), fHistMCFractions(0),
      fHistDataFractions(0), fHistMCWeights(0), fBinsMultCent(0), fBinsPt(0),
      fBinsEta(0), fBinsZv(0) {
    DefineInput(0, TChain::Class());
    // Set default binning
    Double_t binsMultCentDefault[2] = {0, 10000};
    Double_t binsPtDefault[69] = {
        0.,   0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45, 0.5,  0.55,
        0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95, 1.0,  1.1,  1.2,  1.3,
        1.4,  1.5,  1.6,  1.7,  1.8,  1.9,  2.0,  2.2,  2.4,  2.6,  2.8,  3.0,
        3.2,  3.4,  3.6,  3.8,  4.0,  4.5,  5.0,  5.5,  6.0,  6.5,  7.0,  8.0,
        9.0,  10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0,
        26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0, 45.0, 50.0};
    Double_t binsEtaDefault[31] = {
        -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5,
        -0.4, -0.3, -0.2, -0.1, 0.,   0.1,  0.2,  0.3,  0.4,  0.5,  0.6,
        0.7,  0.8,  0.9,  1.0,  1.1,  1.2,  1.3,  1.4,  1.5};
    Double_t binsZvDefault[13] = {-30., -25., -20., -15., -10., -5., 0.,
                                  5.,   10.,  15.,  20.,  25.,  30.};
    SetBinsPt(68, binsPtDefault);
    SetBinsEta(30, binsEtaDefault);
    SetBinsMultCent(1, binsMultCentDefault);
    SetBinsZv(12, binsZvDefault);

    DefineOutput(1, TList::Class());
}

AliMCSpectraWeightsAnalysisTask::~AliMCSpectraWeightsAnalysisTask() {
    if (fOutputList)
        delete fOutputList;
}

void AliMCSpectraWeightsAnalysisTask::UserCreateOutputObjects() {
    OpenFile(1, "recreate");
    fOutputList = new TList();
    fOutputList->SetOwner();
    if (fIsMC && fMCSpectraWeights) {
        // printf("AliMCSpectraWeightsAnalysisTask:: Having non zero
        // AliMCSpectraWeights obj\n");
        if (fDebugLevel > 0)
            printf("AliMCSpectraWeightsAnalysisTask:: obj status: %d\n",
                   fMCSpectraWeights->GetTaskStatus());
        /// Standard track histogram pt:eta:zV:multcent
        Int_t nBinsTrack[4] = {fBinsPt->GetSize() - 1, fBinsEta->GetSize() - 1,
                               fBinsZv->GetSize() - 1,
                               fBinsMultCent->GetSize() - 1};
        Double_t minTrack[4] = {fBinsPt->GetAt(0), fBinsEta->GetAt(0),
                                fBinsZv->GetAt(0), fBinsMultCent->GetAt(0)};
        Double_t maxTrack[4] = {
            fBinsPt->GetAt(fBinsPt->GetSize() - 1),
            fBinsEta->GetAt(fBinsEta->GetSize() - 1),
            fBinsZv->GetAt(fBinsZv->GetSize() - 1),
            fBinsMultCent->GetAt(fBinsMultCent->GetSize() - 1)};

        fHistMCGenPrimTrack = new THnF("fHistMCGenPrimTrackWithWeights",
                                       "Histogram for generated MC Tracks", 4,
                                       nBinsTrack, minTrack, maxTrack);
        fHistMCGenPrimTrack->SetBinEdges(0, fBinsPt->GetArray());
        fHistMCGenPrimTrack->SetBinEdges(1, fBinsEta->GetArray());
        fHistMCGenPrimTrack->SetBinEdges(2, fBinsZv->GetArray());
        fHistMCGenPrimTrack->SetBinEdges(3, fBinsMultCent->GetArray());
        fHistMCGenPrimTrack->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
        fHistMCGenPrimTrack->GetAxis(1)->SetTitle("#eta");
        fHistMCGenPrimTrack->GetAxis(2)->SetTitle("Zv (cm)");
        fHistMCGenPrimTrack->GetAxis(3)->SetTitle("true multiplicity (MC)");
        fHistMCGenPrimTrack->Sumw2();

        fHistMCPartCorr = fMCSpectraWeights->GetHistMCGenPrimTrackParticles();
        fHistDataFractions = fMCSpectraWeights->GetHistDataFraction();
        fHistMCFractions = fMCSpectraWeights->GetHistMCFraction();
        fHistMCWeights = fMCSpectraWeights->GetHistMCWeights();

        fOutputList->Add(fHistMCGenPrimTrack);
        fOutputList->Add(fHistMCPartCorr);
        fOutputList->Add(fHistDataFractions);
        fOutputList->Add(fHistMCFractions);
        fOutputList->Add(fHistMCWeights);
        // fOutputList->Add(fMCSpectraWeights);
    }
    // else printf("AliMCSpectraWeightsAnalysisTask:: Either running not MC or
    // object of AliMCSpectraWeights is null pointer\n");

    PostData(1, fOutputList);
}

void AliMCSpectraWeightsAnalysisTask::UserExec(Option_t* option) {
    if (!fIsMC) {
        printf("AliMCSpectraWeightsAnalysisTask:: WARNING task only works for "
               "MC!!\n");
        return;
    }
    AliInputEventHandler* inputHandler =
        (AliInputEventHandler*)AliAnalysisManager::GetAnalysisManager()
            ->GetInputEventHandler();
    if (!inputHandler) {
        Printf("ERROR: Could not receive inputHandler");
        return;
    }
    //check if trigger mask is fine
    if(!inputHandler->IsEventSelected() & fTriggerMask) return;
    AliPhysicsSelection* physicsSelection =
        static_cast<AliPhysicsSelection*>(inputHandler->GetEventSelection());
    if (!physicsSelection) {
        Printf("ERROR: Could not receive physicsSelection");
        return;
    }
    fEvent = dynamic_cast<AliVEvent*>(InputEvent());
    if (!fEvent) {
        printf("ERROR: fEvent not available\n");
        return;
    }

    fMCEvent = dynamic_cast<AliMCEvent*>(MCEvent());
    if (!fMCEvent) {
        printf("ERROR: fMCEvent not available\n");
        return;
    }
    fMCStack = fMCEvent->Stack();
    if (!fMCStack) {
        printf("ERROR: fMCStack not available\n");
        return;
    }

    /// event rejection/selection
    AliESDEvent* ESDevent = dynamic_cast<AliESDEvent*>(fEvent);
    if (ESDevent->IsIncompleteDAQ())
        return;
    // if (fUtils->IsSPDClusterVsTrackletBG(event)) return;
    if (ESDevent->IsPileupFromSPD(5, 0.8))
        return;

    // Bool_t isEventTriggered = inputHandler->IsEventSelected() &
    // AliMCSpectraWeightsAnalysisTask::GetTriggerMask();
    Double_t zVertEvent = fEvent->GetPrimaryVertex()->GetZ();
    if (TMath::Abs(zVertEvent) > 10.)
        return;

    // Double_t multEvent =
    // AliMCSpectraWeightsAnalysisTask::GetEventMultCent(fEvent);

    // if(fDebugLevel>0) printf("AliMCSpectraWeightsAnalysisTask:: Event got
    // accepted; now analysing\n");

    ///------------------- Loop over Generated Tracks (True
    ///MC)------------------------------
    if (fMCSpectraWeights->GetTaskStatus() <
        AliMCSpectraWeights::TaskState::kMCSpectraObtained)
        fMCSpectraWeights->FillMCSpectra(fMCEvent);
    else
        for (Int_t iParticle = 0; iParticle < fMCStack->GetNtrack();
             iParticle++) {
            if (!fMCStack->IsPhysicalPrimary(iParticle))
                continue; // reject non physical primaries;
            TParticle* mcGenParticle = fMCStack->Particle(iParticle);
            if (!mcGenParticle) {
                printf("ERROR: mcGenParticle  not available\n");
                continue;
            }

            if ((TMath::Abs(mcGenParticle->GetPDG()->Charge()) < 0.01))
                continue; // neutral particle rejection

            /// \li Acceptance cuts for generated particles
            Float_t eta = mcGenParticle->Eta();
            Float_t pt = mcGenParticle->Pt();

            // TPC acceptance
            if (eta < -0.8)
                continue;
            if (eta > 0.8)
                continue;
            // pt range of published spectra
            if (pt < 0.15)
                continue;
            if (pt > 50)
                continue;
            Double_t dWeight =
                fMCSpectraWeights->GetMCSpectraWeight(mcGenParticle, fMCEvent);
            if (fDebugLevel > 0)
                printf("AliMCSpectraWeightsAnalysisTask:: got weight factor "
                       "%lf for pid: %d\n",
                       dWeight, mcGenParticle->GetPdgCode());
            Double_t mcGenPrimTrackValue[4] = {
                mcGenParticle->Pt(), mcGenParticle->Eta(), zVertEvent,
                fMCSpectraWeights->GetMultOrCent()};
            fHistMCGenPrimTrack->Fill(mcGenPrimTrackValue, dWeight);
        }
    PostData(1, fOutputList);
}

/// Function to either fill Multiplicity or centrality
///
/// \param AliVEvent event to be analised
///
/// \return Double_t with centrality or multiplicity
Double_t AliMCSpectraWeightsAnalysisTask::GetEventMultCent(AliVEvent* event) {
    if (fUseMultiplicity) {
        AliVMultiplicity* multiplicity = event->GetMultiplicity();
        if (!multiplicity) {
            printf("ERROR: multiplicity not available\n");
            return 999;
        }
        Int_t mult = multiplicity->GetNumberOfTracklets();
        return mult;
    } else {
        Float_t centralityF = -1;
        AliMultSelection* MultSelection =
            (AliMultSelection*)fEvent->FindListObject("MultSelection");
        if (MultSelection) {
            centralityF = MultSelection->GetMultiplicityPercentile("V0M");
            if (centralityF > 100)
                return 999;
            return centralityF;
        } else {
            AliInfo("Didn't find MultSelection!");
            return 999;
        }
    }
}
