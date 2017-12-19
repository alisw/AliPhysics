#include <iostream>
#include <string>
#include <functional>
#include <numeric>
#include <random>
#include <algorithm>

#include "TH1F.h"
#include "TH3F.h"
#include "TList.h"
#include "TClonesArray.h"

#include "AliAODEvent.h"
#include "AliAODForwardMult.h"
#include "AliVVZERO.h"
#include "AliMultSelection.h"
#include "AliAODCentralMult.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisDataSlot.h"
#include "AliInputEventHandler.h"


#include "AliAnalysisC2Utils.h"
#include "AliAnalysisTaskValidation.h"

using std::cout;
using std::endl;


AliAnalysisTaskValidation::AliAnalysisTaskValidation()
  : AliAnalysisTaskSE(),
    fIsValidEvent(false),
    fEventValidators(),
    fOutputList(0),
    fQA_event_discard_flow(0),
    fQA_track_discard_flow(0),
    fEventCuts(0),
    fUtils(),
    fFMDV0(0),
    fFMDV0_post(0),
    fFMDV0A(0),
    fFMDV0A_post(0),
    fFMDV0C(0),
    fFMDV0C_post(0)
{
}

AliAnalysisTaskValidation::AliAnalysisTaskValidation(const char *name)
  : AliAnalysisTaskSE(name),
    fIsValidEvent(false),
    fEventValidators(),
    fOutputList(0),
    fQA_event_discard_flow(0),
    fQA_track_discard_flow(0),
    fEventCuts(0),
    fUtils(),
    fFMDV0(0),
    fFMDV0_post(0),
    fFMDV0A(0),
    fFMDV0A_post(0),
    fFMDV0C(0),
    fFMDV0C_post(0)
{
  // Apply all cuts by default
  fEventValidators.push_back(EventValidation::kNoEventCut);
  fEventValidators.push_back(EventValidation::kIsAODEvent);
  fEventValidators.push_back(EventValidation::kPassesAliEventCuts);
  fEventValidators.push_back(EventValidation::kHasFMD);
  fEventValidators.push_back(EventValidation::kHasEntriesFMD);
  fEventValidators.push_back(EventValidation::kHasEntriesV0);
  fEventValidators.push_back(EventValidation::kHasValidVertex);
  fEventValidators.push_back(EventValidation::kHasMultSelection);
  // fEventValidators.push_back(EventValidation::kNotOutOfBunchPU);
  // fEventValidators.push_back(EventValidation::kNotMultiVertexPU);
  fEventValidators.push_back(EventValidation::kNotSPDPU);
  fEventValidators.push_back(EventValidation::kNotSPDClusterVsTrackletBG);
  fEventValidators.push_back(EventValidation::kPassesFMD_V0CorrelatioCut);

  // Default track cuts
  fTrackValidators.push_back(TrackValidation::kNoTrackCut);
  fTrackValidators.push_back(TrackValidation::kTPCOnly);
  fTrackValidators.push_back(TrackValidation::kEtaCut);
  fTrackValidators.push_back(TrackValidation::kPtCut);

  // Define output slot
  DefineOutput(1, TList::Class());
  DefineOutput(2, this->Class());

  // Enable mulivertex pileup cuts
  fEventCuts.fPileUpCutMV = true;
}

AliAnalysisTaskValidation* AliAnalysisTaskValidation::ConnectTask(const char *suffix) {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskValidation", "No analysis manager to connect to.");
    return NULL;
  }
  
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer(Form("event_selection_%s", suffix),
			 TList::Class(),
			 AliAnalysisManager::kOutputContainer,
			 Form("%s", mgr->GetCommonFileName()));

  AliAnalysisDataContainer *cExchange =
    mgr->CreateContainer("event_selection_xchange",
			 AliAnalysisTaskValidation::Class(),
			 AliAnalysisManager::kExchangeContainer,
			 Form("%s", mgr->GetCommonFileName()));
  
  AliAnalysisTaskValidation *taskValidation = new AliAnalysisTaskValidation("TaskValidation");
  if (!taskValidation) {
    ::Error("CreateTasks", "Failed to add task!");
    return NULL;
  }
  mgr->AddTask(taskValidation);
  AliAnalysisDataContainer *inputContainer = mgr->GetCommonInputContainer();
  if(!inputContainer) {
    ::Error("CreateTasks", "No input container available. Failed to add task!");
    return NULL;
  }
  mgr->ConnectInput(taskValidation, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskValidation, 1, coutput1);
  mgr->ConnectOutput(taskValidation, 2, cExchange);
  return taskValidation;
}

AliAnalysisDataContainer* AliAnalysisTaskValidation::GetExchangeContainter() {
  // This is the container defined in the ConnectTask function
  return this->GetOutputSlot(2)->GetContainer();
}

void AliAnalysisTaskValidation::CreateQAHistograms(TList* outlist) {
  /// Event discard flow histogram
  this->fQA_event_discard_flow = new TH1F("qa_discard_flow",
					  "QA event discard flow",
					  this->fEventValidators.size(),
					  0,
					  this->fEventValidators.size());
  TAxis *discardedEvtsAx = this->fQA_event_discard_flow->GetXaxis();

  for (UInt_t idx = 0; idx < this->fEventValidators.size(); idx++) {
    switch (this->fEventValidators[idx]) {
    case EventValidation::kNoEventCut:
      discardedEvtsAx->SetBinLabel(idx + 1, "No cuts"); break;
    case EventValidation::kIsAODEvent:
      discardedEvtsAx->SetBinLabel(idx + 1, "AOD event"); break;
    case EventValidation::kHasFMD:
      discardedEvtsAx->SetBinLabel(idx + 1, "Has FMD"); break;
    case EventValidation::kHasEntriesFMD:
      discardedEvtsAx->SetBinLabel(idx + 1, "Has entries FMD"); break;
    case EventValidation::kHasEntriesV0:
      discardedEvtsAx->SetBinLabel(idx + 1, "Has entries V0"); break;
    case EventValidation::kPassesAliEventCuts:
      discardedEvtsAx->SetBinLabel(idx + 1, "AliEventCuts"); break;
    case EventValidation::kPassesFMD_V0CorrelatioCut:
      discardedEvtsAx->SetBinLabel(idx + 1, "FMD V0 correlation"); break;
    case EventValidation::kHasValidVertex:
      discardedEvtsAx->SetBinLabel(idx + 1, "Valid vertex"); break;
    case EventValidation::kHasMultSelection:
      discardedEvtsAx->SetBinLabel(idx + 1, "Has MultSelection"); break;
    case EventValidation::kNotOutOfBunchPU:
      discardedEvtsAx->SetBinLabel(idx + 1, "Not out-of-bunch PU"); break;
    case EventValidation::kNotMultiVertexPU:
      discardedEvtsAx->SetBinLabel(idx + 1, "Not multi-vertex PU"); break;
    case EventValidation::kNotSPDPU:
      discardedEvtsAx->SetBinLabel(idx + 1, "Not SPD PU"); break;
    case EventValidation::kNotSPDClusterVsTrackletBG:
      discardedEvtsAx->SetBinLabel(idx + 1, "SPD clstrs vs BG cut"); break;
    }
  }
  outlist->Add(this->fQA_event_discard_flow);

  /// Track discard flow
  this->fQA_track_discard_flow = new TH1F("qa_tack_discard_flow",
					  "QA track discard flow",
					  this->fTrackValidators.size(),
					  0,
					  this->fTrackValidators.size());
  TAxis *discardedTrksAx = this->fQA_track_discard_flow->GetXaxis();

  for (UInt_t idx = 0; idx < this->fTrackValidators.size(); idx++) {
    switch (this->fTrackValidators[idx]) {
    case TrackValidation::kNoTrackCut:
      discardedTrksAx->SetBinLabel(idx + 1, "No cuts"); break;
    case TrackValidation::kTPCOnly:
      discardedTrksAx->SetBinLabel(idx + 1, "!TPC-Only"); break;
    case TrackValidation::kEtaCut:
      discardedTrksAx->SetBinLabel(idx + 1, "!Eta cut"); break;
    case TrackValidation::kPtCut:
      discardedTrksAx->SetBinLabel(idx + 1, "!Pt cut"); break;
    }
  }
  outlist->Add(this->fQA_track_discard_flow);  
}

void AliAnalysisTaskValidation::UserCreateOutputObjects() {
  // Stop right here if there are no Validators to work with
  if (this->fEventValidators.size() == 0) {
    AliFatal("No event validators specified!");
  }

  // Setup output list
  this->fOutputList = new TList();
  this->fOutputList->SetOwner();

  // Create QA histograms in Event selection
  fEventCuts.AddQAplotsToList(this->fOutputList);
  this->CreateQAHistograms(this->fOutputList);

  // FMD V0 QA histograms
  this->fFMDV0 = new TH2F("FMDV0", "FMD vs V0 pre cut;FMD;V0;",
			  2000, 0, 2000, 2000, 0, 2000);
  this->fOutputList->Add(this->fFMDV0);

  this->fFMDV0_post = new TH2F("FMDV0_post", "FMD vs V0 post cut;FMD;V0;",
			  2000, 0, 2000, 2000, 0, 2000);
  this->fOutputList->Add(this->fFMDV0_post);

  this->fFMDV0A = new TH2F("FMDV0A", "FMD vs V0A;FMD;V0A;",
			   1000, 0, 1000, 1000, 0, 1000);
  this->fOutputList->Add(this->fFMDV0A);
  this->fFMDV0A_post = new TH2F("FMDV0A_post", "FMD vs V0A post cut;FMD;V0A;",
				1000, 0, 1000, 1000, 0, 1000);
  this->fOutputList->Add(this->fFMDV0A_post);

  this->fFMDV0C = new TH2F("FMDV0C", "FMD vs V0C;FMD;V0C;",
			   1000, 0, 1000, 1000, 0, 1000);
  this->fOutputList->Add(this->fFMDV0C);
  this->fFMDV0C_post = new TH2F("FMDV0C_post", "FMD vs V0C post cut;FMD;V0C;",
				1000, 0, 1000, 1000, 0, 1000);
  this->fOutputList->Add(this->fFMDV0C_post);

  // Slot 0 is reserved; 1 needs to be called here to get at least empty histograms
  PostData(1, fOutputList);
  // Do NOT post the exchange container here!
  // I.E.: Not PostData(2, this);
}

void AliAnalysisTaskValidation::UserExec(Option_t *)
{
  this->fIsValidEvent = true;
  for (UInt_t idx = 0; idx < this->fEventValidators.size(); idx++) {
    switch (this->fEventValidators[idx]) {
    case EventValidation::kNoEventCut:
      this->fIsValidEvent = this->NoCut(); break;
    case EventValidation::kIsAODEvent:
      this->fIsValidEvent = this->IsAODEvent(); break;
    case EventValidation::kHasFMD:
      this->fIsValidEvent = this->HasFMD(); break;
    case EventValidation::kHasEntriesFMD:
      this->fIsValidEvent = this->HasEntriesFMD(); break;
    case EventValidation::kHasEntriesV0:
      this->fIsValidEvent = this->HasEntriesV0(); break;
    case EventValidation::kPassesAliEventCuts:
      this->fIsValidEvent = this->PassesAliEventCuts(); break;
    case EventValidation::kPassesFMD_V0CorrelatioCut:
      this->fIsValidEvent = this->PassesFMDV0CorrelatioCut(true); break;
    case EventValidation::kHasValidVertex:
      this->fIsValidEvent = this->HasValidVertex(); break;
    case EventValidation::kHasMultSelection:
      this->fIsValidEvent = this->HasMultSelection(); break;
    case EventValidation::kNotOutOfBunchPU:
      this->fIsValidEvent = this->NotOutOfBunchPU(); break;
    case EventValidation::kNotMultiVertexPU:
      this->fIsValidEvent = this->NotMultiVertexPU(); break;
    case EventValidation::kNotSPDPU:
      this->fIsValidEvent = this->NotSPDPU(); break;
    case EventValidation::kNotSPDClusterVsTrackletBG:
      this->fIsValidEvent = this->NotSPDClusterVsTrackletBG(); break;
    }
    if (this->fIsValidEvent) {
      this->fQA_event_discard_flow->Fill(idx);
    } else {
      // Stop checking once this event has been flaged as invalid
      break;
    }
  }
  PostData(1, fOutputList);
  // Make the current pointer to this task avaialble in the exchange container
  PostData(2, this);
}

Bool_t AliAnalysisTaskValidation::IsAODEvent() {
  return dynamic_cast<AliAODEvent*>(this->InputEvent()) ? true : false;
}

Bool_t AliAnalysisTaskValidation::HasFMD() {
  AliAODEvent* aod = dynamic_cast<AliAODEvent*>(this->InputEvent());
  AliAODForwardMult* aodForward =
    dynamic_cast<AliAODForwardMult*>(aod->FindListObject("Forward"));
  if (!aodForward) {
    return false;
  }
  return true;
}

Bool_t AliAnalysisTaskValidation::HasEntriesFMD() {
  if (this->HasFMD() && this->GetFMDhits().size() > 0) {
    return true;
  }
  return false;
}

Bool_t AliAnalysisTaskValidation::HasEntriesV0() {
  if (this->GetV0hits().size() > 0) {
    return true;
  }
  return false;
}

Bool_t AliAnalysisTaskValidation::PassesAliEventCuts() {
  return this->fEventCuts.AcceptEvent(this->InputEvent());
}

Bool_t AliAnalysisTaskValidation::PassesFMDV0CorrelatioCut(Bool_t fill_qa) {
  // Overlap regions between the two detectors
  // Float_t fmd_v0a_overlap[2] = {2.8, 5.03};
  // Float_t fmd_v0c_overlap[2] = {-3.4, -2.01};
  Tracks v0hits = this->GetV0hits();
  Tracks fmdhits = this->GetFMDhits();

  // Calculate hits on each detector in overlap region
  Float_t nV0A_hits =
    std::accumulate(v0hits.begin(), v0hits.end(), 0,
		    [](Float_t a, AliAnalysisTaskValidation::Track t) {
		      return a + ((2.8 < t.eta && t.eta < 5.03) ? t.weight : 0.0f);
		    });
  Float_t nFMD_fwd_hits =
    std::accumulate(fmdhits.begin(), fmdhits.end(), 0,
		    [](Float_t a, AliAnalysisTaskValidation::Track t) {
		      return a + ((2.8 < t.eta && t.eta < 5.03) ? t.weight : 0.0f);
		    });
  Float_t nV0C_hits =
    std::accumulate(v0hits.begin(), v0hits.end(), 0,
		    [](Float_t a, AliAnalysisTaskValidation::Track t) {
		      return a + ((-3.4 < t.eta && t.eta < 2.01) ? t.weight : 0.0f);
		    });
  Float_t nFMD_bwd_hits =
    std::accumulate(fmdhits.begin(), fmdhits.end(), 0,
		    [](Float_t a, AliAnalysisTaskValidation::Track t) {
		      return a + ((-3.4 < t.eta && t.eta < 2.01) ? t.weight : 0.0f);
		    });
  if (fill_qa) {
    this->fFMDV0->Fill(nFMD_bwd_hits + nFMD_fwd_hits, nV0C_hits + nV0A_hits);
    this->fFMDV0A->Fill(nFMD_fwd_hits, nV0A_hits);
    this->fFMDV0C->Fill(nFMD_bwd_hits, nV0C_hits);
  }

  // Cut on V0 - FMD outliers outliers
  //  if (nV0A_hits + nV0C_hits < (nFMD_fwd_hits + nFMD_bwd_hits - 40)) {
  if (nV0A_hits + nV0C_hits < 1.5*(nFMD_fwd_hits + nFMD_bwd_hits) - 20) {
    return false;
  }
  if (fill_qa) {
    this->fFMDV0_post->Fill(nFMD_bwd_hits + nFMD_fwd_hits, nV0C_hits + nV0A_hits);
    this->fFMDV0A_post->Fill(nFMD_fwd_hits, nV0A_hits);
    this->fFMDV0C_post->Fill(nFMD_bwd_hits, nV0C_hits);
  }

  return true;
}

Bool_t AliAnalysisTaskValidation::HasMultSelection() {
  return
    dynamic_cast< AliMultSelection* >(InputEvent()->FindListObject("MultSelection")) ?
    true : false;
}

Bool_t AliAnalysisTaskValidation::HasValidVertex() {
  if (!this->InputEvent()->GetPrimaryVertex()
      || !this->InputEvent()->GetPrimaryVertex()->GetZ()) {
    return false;
  } else {
    return true;
  }
}

AliAnalysisTaskValidation::Tracks AliAnalysisTaskValidation::GetFMDhits() const {
    // Relies on the event being vaild (no extra checks if object exists done here)
  AliAODForwardMult* aodForward =
    static_cast<AliAODForwardMult*>(fInputEvent->FindListObject("Forward"));
  // Shape of d2Ndetadphi: 200, -4, 6, 20, 0, 2pi
  const TH2D& d2Ndetadphi = aodForward->GetHistogram();
  Int_t nEta = d2Ndetadphi.GetXaxis()->GetNbins();
  Int_t nPhi = d2Ndetadphi.GetYaxis()->GetNbins();
  Tracks ret_vector;
  // FMD has no pt resolution!
  Double_t pt = 0;
  for (Int_t iEta = 1; iEta <= nEta; iEta++) {
    Int_t valid = Int_t(d2Ndetadphi.GetBinContent(iEta, 0));
    if (!valid) {
       // No data expected for this eta
      continue;
    }
    Float_t eta = d2Ndetadphi.GetXaxis()->GetBinCenter(iEta);
    for (Int_t iPhi = 1; iPhi <= nPhi; iPhi++) {
      // Bin content is most likely number of particles!
      Float_t mostProbableN = d2Ndetadphi.GetBinContent(iEta, iPhi);
      if (mostProbableN > 0) {
	Float_t phi = d2Ndetadphi.GetYaxis()->GetBinCenter(iPhi);
	ret_vector.push_back(AliAnalysisTaskValidation::Track(eta, phi, pt, mostProbableN));
      }
    }
  }
  // Above, we looped through the histogram in an orderly
  // fashion. This means that the eta and phi values are in increasing
  // orders. For eta this is not an issue. We will just always have
  // the trigger (or associated) with a larger eta. In a eta1, eta2
  // plots, this will show by all values being in one triangle. That
  // is ok. However, for phi this is not ok if we require
  // eta1==eta2. In that case, phi2 > phi1. This leads to a very
  // shifted correlation function. Especially if the eta bin width is
  // larger than the resolution. Bottom line, I suffle this vector
  // here to avoid this mess.

  // #times I thought this was unnecessary but turned it back on later: 3

  std::random_device rd;
  std::default_random_engine engine{rd()};
  std::shuffle(std::begin(ret_vector), std::end(ret_vector), engine);
  return ret_vector;
}

AliAnalysisTaskValidation::Tracks AliAnalysisTaskValidation::GetV0hits() const {
  // Relies on the event being vaild (no extra checks if object exists done here)
  Tracks ret_vector;
  AliVVZERO* vzero = this->InputEvent()->GetVZEROData();
  const Int_t nV0_channels = 64;
  // V0 has no pt resolution!
  Double_t pt = 0;
  for (Int_t ich = 0; ich < nV0_channels; ich++) {
    //vzero->GetMultiplicity(ich);
    Float_t amp = this->InputEvent()->GetVZEROEqMultiplicity(ich);
    if (amp > 0) {
      Float_t eta = 0.5 * (vzero->GetVZEROEtaMin(ich) + vzero->GetVZEROEtaMax(ich));
      Float_t phi = vzero->GetVZEROAvgPhi(ich);
      ret_vector.push_back(
	AliAnalysisTaskValidation::Track(eta, AliAnalysisC2Utils::Wrap02pi(phi), pt, amp));
    }
  }
  return ret_vector;
}

AliAnalysisTaskValidation::Tracks AliAnalysisTaskValidation::GetTracks() {
  auto in_tracks = this->GetAllCentralBarrelTracks();
  Tracks out_tracks;
  for (auto obj: *in_tracks) {
    auto tr = static_cast<AliAODTrack*>(obj);
    auto track_valid = true;
    for (UInt_t idx = 0; idx < this->fEventValidators.size(); idx++) {
      switch (this->fTrackValidators[idx]) {
      case TrackValidation::kNoTrackCut:
	break;
      case TrackValidation::kTPCOnly:
	track_valid = tr->TestFilterBit(128); break;
      case TrackValidation::kEtaCut:
	track_valid = (TMath::Abs(tr->Eta()) < 0.9); break;
      case TrackValidation::kPtCut:
	track_valid = (tr->Pt() < 0.3); break;
      }
      if (track_valid) {
	this->fQA_event_discard_flow->Fill(idx);
      } else {
	// Stop checking once this event has been flaged as invalid
	break;
      }
    }
    if (track_valid) {
      Double_t weight = 1.0;
      out_tracks
	.push_back(AliAnalysisTaskValidation::Track(tr->Eta(),
						    AliAnalysisC2Utils::Wrap02pi(tr->Phi()),
						    tr->Pt(),
						    weight));
    }
  }
  PostData(1, fOutputList);
  return out_tracks;
}

AliAnalysisTaskValidation::Tracks AliAnalysisTaskValidation::GetTracklets() const {
  const Double_t dummy_pt = 0;
  AliAnalysisTaskValidation::Tracks ret_vector;

  // Using the aptly named parent class AliVMultiplicity
  AliVMultiplicity* mult = this->fInputEvent->GetMultiplicity();
  Int_t nTracklets = mult->GetNumberOfTracklets();
  for (Int_t i=0; i < nTracklets; i++){
    // Using a dphi cut in units of mrad; This cut is motivated in
    // https://aliceinfo.cern.ch/Notes/sites/aliceinfo.cern.ch.Notes/files/notes/analysis/lmilano/2017-Aug-11-analysis_note-note.pdf
    auto dphi  = mult->GetDeltaPhi(i);
    if (TMath::Abs(dphi) * 1000 > 5) {
      continue;
    }
    auto phi   = mult->GetPhi(i);
    auto eta   = -TMath::Log(TMath::Tan(mult->GetTheta(i)/2));
    ret_vector.push_back(AliAnalysisTaskValidation::Track(eta, phi, dummy_pt, 1));
  }
  return ret_vector;
}

AliAnalysisTaskValidation::Tracks AliAnalysisTaskValidation::GetSPDclusters() const {
  // Relies on the event being vaild (no extra checks if object exists done here)
  AliAODCentralMult* aodcmult = static_cast<AliAODCentralMult*>
    (fInputEvent->FindListObject("CentralClusters")); // Shape of d2Ndetadphi: 200, -4, 6, 20, 0, 2pi
  const TH2D& d2Ndetadphi = aodcmult->GetHistogram();
  const Int_t nEta = d2Ndetadphi.GetXaxis()->GetNbins();
  const Int_t nPhi = d2Ndetadphi.GetYaxis()->GetNbins();
  const Double_t pt = 0;
  AliAnalysisTaskValidation::Tracks ret_vector;

  for (Int_t iEta = 1; iEta <= nEta; iEta++) {
    Int_t valid = Int_t(d2Ndetadphi.GetBinContent(iEta, 0));
    if (!valid) {
      // No data expected for this eta
      continue;
    }
    Float_t eta = d2Ndetadphi.GetXaxis()->GetBinCenter(iEta);
    // FIXME: For now, drop everything outside of -1.7 < eta 1.7
    if (eta < -1.7 || eta > 1.7) {
      continue;
    }
    for (Int_t iPhi = 1; iPhi <= nPhi; iPhi++) {
      // Bin content is most likely number of particles!
      Float_t mostProbableN = d2Ndetadphi.GetBinContent(iEta, iPhi);
      if (mostProbableN > 0) {
	Float_t phi = d2Ndetadphi.GetYaxis()->GetBinCenter(iPhi);
	ret_vector.push_back(AliAnalysisTaskValidation::Track(eta, phi, pt, mostProbableN));
      }
    }
  }
  // See the reasoning for the following suffle in the code that gets the hits on the FMD
  std::random_device rd;
  std::default_random_engine engine{rd()};
  std::shuffle(std::begin(ret_vector), std::end(ret_vector), engine);
  return ret_vector;
}

TClonesArray* AliAnalysisTaskValidation::GetAllCentralBarrelTracks() {
  // If we are dealing with an ESD event, we have to have an AOD handler as well!
  // We get all the particles/tracks from this AOD handler.
  AliAODEvent* aodEvent = dynamic_cast< AliAODEvent* >(this->InputEvent());
  return aodEvent->GetTracks();
}

TClonesArray* AliAnalysisTaskValidation::GetAllMCTruthTracks() {
  // If we are dealing with an ESD event, we have to have an AOD handler as well!
  // We get all the particles/tracks from this AOD handler.
  AliAODEvent* aodEvent = dynamic_cast< AliAODEvent* >(this->InputEvent());
  // Yes, the following is realy "aodEvent" not mcEvent :P
  return
    dynamic_cast<TClonesArray*>(aodEvent->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
}

Bool_t AliAnalysisTaskValidation::UserNotify() {
  // Turn off all branches
  // return true;
  this->fInputHandler->GetTree()->SetBranchStatus("*", false, 0);
  // Turn back on all branches which got dumped into the top-level.
  // These are ~250 branches - hurray!
  this->fInputHandler->GetTree()->SetBranchStatus("f*", true);
  // Multiplicity framework; very expensive to read!
  // this->fInputHandler->GetTree()->SetBranchStatus("MultSelection", true);
  // Turn on headers;
  this->fInputHandler->GetTree()->SetBranchStatus("header", true);
  // Somehow vertices have to be turned on as glob...
  this->fInputHandler->GetTree()->SetBranchStatus("vertices.*", true);
  // ... but tracks individually for sub branches
  this->fInputHandler->GetTree()->SetBranchStatus("tracks", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fUniqueID", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fBits", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fMomentum[3]", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fPosition[3]", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fMomentumAtDCA[3]", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fPositionAtDCA[2]", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fRAtAbsorberEnd", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fChi2perNDF", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fChi2MatchTrigger", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fITSchi2", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fFlags", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fLabel", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fTOFLabel[3]", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fTrackLength", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fITSMuonClusterMap", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fMUONtrigHitsMapTrg", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fMUONtrigHitsMapTrk", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fFilterMap", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fTPCFitMap.fUniqueID", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fTPCFitMap.fBits", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fTPCFitMap.fNbits", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fTPCFitMap.fNbytes", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fTPCFitMap.fAllBits", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fTPCClusterMap.fUniqueID", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fTPCClusterMap.fBits", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fTPCClusterMap.fNbits", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fTPCClusterMap.fNbytes", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fTPCClusterMap.fAllBits", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fTPCSharedMap.fUniqueID", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fTPCSharedMap.fBits", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fTPCSharedMap.fNbits", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fTPCSharedMap.fNbytes", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fTPCSharedMap.fAllBits", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fTPCnclsF", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fTPCNCrossedRows", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fID", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fCharge", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fType", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fPIDForTracking", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fCaloIndex", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fProdVertex", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fTrackPhiOnEMCal", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fTrackEtaOnEMCal", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fTrackPtOnEMCal", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fIsMuonGlobalTrack", true);
  this->fInputHandler->GetTree()->SetBranchStatus("tracks.fMFTClusterPattern", true);


  
  // this->fInputHandler->GetTree()->SetBranchStatus("tracks.*", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("v0s.*", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("cascades.*", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("tracklets", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("jets.*", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("emcalCells", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("phosCells", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("caloClusters.*", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("emcalTrigger", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("phosTrigger", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fmdClusters.*", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("pmdClusters.*", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("hmpidRings.*", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("dimuons.*", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("AliAODTZERO", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("AliAODVZERO", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("AliAODZDC", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("AliAODAD", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("AliTOFHeader", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("trdTracks.*", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("Forward", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("ForwardEP", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("CentralClusters", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("MultSelection", true);
  
  // this->fInputHandler->GetTree()->SetBranchStatus("fUniqueID", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fBits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fName", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fTitle", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fMagneticField", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fMuonMagFieldScale", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fCentrality", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fEventplane", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fEventplaneMag", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fEventplaneQx", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fEventplaneQy", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fZDCN1Energy", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fZDCP1Energy", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fZDCN2Energy", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fZDCP2Energy", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fZDCEMEnergy[2]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fNQTheta", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fQTheta", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fTriggerMask", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fTriggerMaskNext50", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fFiredTriggers", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fRunNumber", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fRefMult", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fRefMultPos", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fRefMultNeg", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fNMuons", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fNDimuons", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fNGlobalMuons", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fNGlobalDimuons", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fDAQAttributes", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fEventType", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fOrbitNumber", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fPeriodNumber", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fBunchCrossNumber", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fRefMultComb05", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fRefMultComb08", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fRefMultComb10", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fTriggerCluster", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fDiamondXY[2]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fDiamondCovXY[3]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fDiamondZ", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fDiamondSig2Z", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fOfflineTrigger", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fESDFileName", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fEventNumberESDFile", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fNumberESDTracks", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fL0TriggerInputs", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fL1TriggerInputs", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fL2TriggerInputs", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fITSClusters[6]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fTPConlyRefMult", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fVZEROEqFactors[64]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fT0spread[4]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fIRInt2InteractionsMap.fUniqueID", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fIRInt2InteractionsMap.fBits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fIRInt2InteractionsMap.fNbits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fIRInt2InteractionsMap.fNbytes", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fIRInt2InteractionsMap.fAllBits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fIRInt1InteractionsMap.fUniqueID", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fIRInt1InteractionsMap.fBits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fIRInt1InteractionsMap.fNbits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fIRInt1InteractionsMap.fNbytes", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fIRInt1InteractionsMap.fAllBits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("vertices", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("vertices.fUniqueID", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("vertices.fBits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("vertices.fName", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("vertices.fTitle", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("vertices.fPosition[3]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("vertices.fChi2perNDF", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("vertices.fID", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("vertices.fBCID", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("vertices.fType", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("vertices.fNprong", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("vertices.fNContributors", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("vertices.fParent", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("vertices.fDaughters", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("v0s", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("v0s.fUniqueID", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("v0s.fBits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("v0s.fSecondaryVtx", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("v0s.fCharge", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("v0s.fNProngs", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("v0s.fNDCA", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("v0s.fNPID", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("v0s.fPx", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("v0s.fPy", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("v0s.fPz", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("v0s.fd0", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("v0s.fDCA", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("v0s.fPID", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("v0s.fDcaV0ToPrimVertex", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("v0s.fOnFlyStatus", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("cascades", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("cascades.fUniqueID", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("cascades.fBits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("cascades.fSecondaryVtx", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("cascades.fCharge", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("cascades.fNProngs", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("cascades.fNDCA", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("cascades.fNPID", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("cascades.fPx", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("cascades.fPy", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("cascades.fPz", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("cascades.fd0", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("cascades.fDCA", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("cascades.fPID", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("cascades.fDcaV0ToPrimVertex", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("cascades.fOnFlyStatus", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("cascades.fDecayVertexXi", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("cascades.fChargeXi", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("cascades.fDcaXiDaughters", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("cascades.fDcaXiToPrimVertex", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("cascades.fDcaBachToPrimVertex", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("cascades.fMomBachX", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("cascades.fMomBachY", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("cascades.fMomBachZ", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fUniqueID", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fBits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fName", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fTitle", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fNTracks", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fTheta", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fPhi", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fDeltaPhi", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fLabels", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fLabelsL2", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fFiredChips[2]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fITSClusters[6]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fFastOrFiredChips.fUniqueID", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fFastOrFiredChips.fBits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fFastOrFiredChips.fNbits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fFastOrFiredChips.fNbytes", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fFastOrFiredChips.fAllBits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fClusterFiredChips.fUniqueID", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fClusterFiredChips.fBits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fClusterFiredChips.fNbits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fClusterFiredChips.fNbytes", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fClusterFiredChips.fAllBits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("jets", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("jets.fUniqueID", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("jets.fBits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("jets.fBackgEnergy[2]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("jets.fEffectiveArea[2]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("jets.fEffectiveAreaError[2]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("jets.fNeutralFraction", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("jets.fPtSubtracted[2]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("jets.fPtLeadingConstituent", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("jets.fTrigger", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fUniqueID", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fBits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fName", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fTitle", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fNCells", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fHGLG", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fCellNumber", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fAmplitude", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fTime", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fEFraction", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fMCLabel", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fType", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fUniqueID", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fBits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fName", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fTitle", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fNCells", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fHGLG", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fCellNumber", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fAmplitude", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fTime", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fEFraction", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fMCLabel", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fType", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("caloClusters", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("caloClusters.fUniqueID", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("caloClusters.fBits", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("caloClusters.fEnergy", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("caloClusters.fPosition[3]", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("caloClusters.fChi2", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("caloClusters.fPID[13]", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("caloClusters.fID", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("caloClusters.fNLabel", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("caloClusters.fLabel", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("caloClusters.fFilterMap", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("caloClusters.fType", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("caloClusters.fClusterMCEdepFraction", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("caloClusters.fDistToBadChannel", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("caloClusters.fDispersion", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("caloClusters.fM20", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("caloClusters.fM02", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("caloClusters.fEmcCpvDistance", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("caloClusters.fTrackDx", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("caloClusters.fTrackDz", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("caloClusters.fNExMax", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("caloClusters.fTOF", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("caloClusters.fCoreEnergy", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("caloClusters.fTracksMatched", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("caloClusters.fNCells", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("caloClusters.fCellsAbsId", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("caloClusters.fCellsAmpFraction", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("caloClusters.fCellsMCEdepFractionMap", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fUniqueID", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fBits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fName", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fTitle", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fNEntries", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fCurrent", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fColumn", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fRow", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fAmplitude", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fTime", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fNL0Times", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fL1TimeSum", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fTriggerBits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fL1Threshold[4]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fL1V0[2]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fL1FrameMask", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fL1DCALThreshold[4]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fL1SubRegion", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fL1DCALFrameMask", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fMedian[2]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fTriggerBitWord", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fL1DCALV0[2]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fUniqueID", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fBits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fName", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fTitle", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fNEntries", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fCurrent", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fColumn", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fRow", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fAmplitude", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fTime", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fNL0Times", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fL1TimeSum", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fTriggerBits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fL1Threshold[4]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fL1V0[2]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fL1FrameMask", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fL1DCALThreshold[4]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fL1SubRegion", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fL1DCALFrameMask", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fMedian[2]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fTriggerBitWord", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fL1DCALV0[2]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fmdClusters", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fmdClusters.fUniqueID", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fmdClusters.fBits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fmdClusters.fEnergy", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fmdClusters.fPosition[3]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fmdClusters.fChi2", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fmdClusters.fPID[13]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fmdClusters.fID", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fmdClusters.fNLabel", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fmdClusters.fLabel", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fmdClusters.fFilterMap", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fmdClusters.fType", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fmdClusters.fClusterMCEdepFraction", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fmdClusters.fProdVertex", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("fmdClusters.fPrimTrack", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("pmdClusters", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("pmdClusters.fUniqueID", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("pmdClusters.fBits", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("pmdClusters.fEnergy", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("pmdClusters.fPosition[3]", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("pmdClusters.fChi2", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("pmdClusters.fPID[13]", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("pmdClusters.fID", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("pmdClusters.fNLabel", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("pmdClusters.fLabel", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("pmdClusters.fFilterMap", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("pmdClusters.fType", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("pmdClusters.fClusterMCEdepFraction", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("pmdClusters.fAssocCluster", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("hmpidRings", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("hmpidRings.fUniqueID", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("hmpidRings.fBits", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("hmpidRings.fHmpidAODtrkId", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("hmpidRings.fHmpidAODqn", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("hmpidRings.fHmpidAODcluIdx", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("hmpidRings.fHmpidAODtrkTheta", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("hmpidRings.fHmpidAODtrkPhi", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("hmpidRings.fHmpidAODsignal", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("hmpidRings.fHmpidAODocc", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("hmpidRings.fHmpidAODchi2", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("hmpidRings.fHmpidAODtrkX", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("hmpidRings.fHmpidAODtrkY", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("hmpidRings.fHmpidAODmipX", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("hmpidRings.fHmpidAODmipY", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("hmpidRings.fHmpidAODpid[5]", true);
  // // this->fInputHandler->GetTree()->SetBranchStatus("hmpidRings.fHMPIDmom[3]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("dimuons", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("dimuons.fUniqueID", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("dimuons.fBits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("dimuons.fMu[2]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fUniqueID", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fBits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fT0TOF[3]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fPileup", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fSattelite", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fBackground", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fT0TOFbest[3]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fT0VertexRaw", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fT0zVertex", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fT0Amp[26]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fPileupBits.fUniqueID", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fPileupBits.fBits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fPileupBits.fNbits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fPileupBits.fNbytes", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fPileupBits.fAllBits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fUniqueID", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fBits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fBBtriggerV0A", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fBGtriggerV0A", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fBBtriggerV0C", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fBGtriggerV0C", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fMultiplicity[64]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fBBFlag[64]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fBGFlag[64]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fV0ATime", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fV0CTime", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fV0ADecision", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fV0CDecision", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fTriggerChargeA", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fTriggerChargeC", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fTriggerBits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fIsBB[64][21]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fIsBG[64][21]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fUniqueID", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fBits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fZEM1Energy", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fZEM2Energy", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fZNCTowerEnergy[5]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fZNATowerEnergy[5]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fZPCTowerEnergy[5]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fZPATowerEnergy[5]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fZNCTowerEnergyLR[5]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fZNATowerEnergyLR[5]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fZPCTowerEnergyLR[5]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fZPATowerEnergyLR[5]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fZDCParticipants", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fZDCPartSideA", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fZDCPartSideC", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fImpactParameter", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fImpactParamSideA", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fImpactParamSideC", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fZDCTDCSum", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fZDCTDCDifference", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fZNCTDC", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fZNATDC", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fZPCTDC", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fZPATDC", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fZNCTDCm[4]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fZNATDCm[4]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fZPCTDCm[4]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fZPATDCm[4]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fIsZNAfired", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fIsZNCfired", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fIsZPAfired", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fIsZPCfired", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fUniqueID", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fBits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fBBtriggerADA", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fBGtriggerADA", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fBBtriggerADC", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fBGtriggerADC", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fMultiplicity[16]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fBBFlag[16]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fBGFlag[16]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fADATime", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fADCTime", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fADADecision", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fADCDecision", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fTriggerChargeA", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fTriggerChargeC", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fTriggerBits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fIsBB[16][21]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fIsBG[16][21]", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fUniqueID", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fBits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fDefaultEventTimeValue", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fDefaultEventTimeRes", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fNbins", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fTOFtimeResolution", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fT0spread", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fNumberOfTOFclusters", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("fNumberOfTOFtrgPads", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("trdTracks", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("trdTracks.fUniqueID", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("trdTracks.fBits", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("trdTracks.fGlobalStack", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("trdTracks.fPID", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("trdTracks.fLayerMask", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("trdTracks.fA", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("trdTracks.fFlagsTiming", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("trdTracks.fTracklets", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("trdTracks.fTrackMatch", true);
  // this->fInputHandler->GetTree()->SetBranchStatus("trdTracks.fLabel", true);

  // this->fInputHandler->GetTree()->SetBranchStatus("*", true);
  return true;
}
