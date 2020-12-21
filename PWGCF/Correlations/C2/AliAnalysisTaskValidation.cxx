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
#include "AliMCEvent.h"
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
    fTrackValidators(),
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

AliAnalysisTaskValidation::AliAnalysisTaskValidation(const char *name, bool is_reconstructed)
  : AliAnalysisTaskSE(name),
    fIsValidEvent(false),
    fEventValidators(),
    fTrackValidators(),
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
  if (is_reconstructed) {
    fEventValidators.push_back(EventValidation::kIsAODEvent);
    fEventValidators.push_back(EventValidation::kPassesAliEventCuts);
    fEventValidators.push_back(EventValidation::kHasFMD);
    fEventValidators.push_back(EventValidation::kHasEntriesFMD);
    fEventValidators.push_back(EventValidation::kHasEntriesV0);
    fEventValidators.push_back(EventValidation::kHasValidVertex);
    fEventValidators.push_back(EventValidation::kHasMultSelection);
    // This one kills another 60% of the events in LHC15o HIR :/
    // fEventValidators.push_back(EventValidation::kNotOutOfBunchPU);
    fEventValidators.push_back(EventValidation::kNotMultiVertexPU);
    fEventValidators.push_back(EventValidation::kNotSPDPU);
    fEventValidators.push_back(EventValidation::kNotSPDClusterVsTrackletBG);
    fEventValidators.push_back(EventValidation::kPassesFMD_V0CorrelatioCut);
  }
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

AliAnalysisTaskValidation* AliAnalysisTaskValidation::ConnectTask(const char *suffix,
								  bool is_reconstructed) {
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
  
  auto *taskValidation = new AliAnalysisTaskValidation("TaskValidation", is_reconstructed);
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
    auto eta   = -TMath::Log(TMath::Tan(mult->GetTheta(i)/2));
    // Drop everything outside of -1.7 < eta 1.7 to avoid overlas with the FMD
    if (eta < -1.7 || eta > 1.7) {
      continue;
    }
    auto phi   = mult->GetPhi(i);
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

TClonesArray* AliAnalysisTaskValidation::GetAllMCTruthTracksAsTClonesArray() {
  // If we are dealing with an ESD event, we have to have an AOD handler as well!
  // We get all the particles/tracks from this AOD handler.
  auto aodEvent = dynamic_cast< AliAODEvent* >(this->InputEvent());
  if (!aodEvent) {
    AliFatal("No AOD event found");
  }
  auto tr_arr = dynamic_cast<TClonesArray*>
    (aodEvent->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
  if(!tr_arr){
    AliFatal("No MC array found in AOD");
  }
  return tr_arr;
}

AliAnalysisTaskValidation::Tracks AliAnalysisTaskValidation::GetMCTruthTracks() {
  AliAnalysisTaskValidation::Tracks ret_vector;
  if (this->IsAODEvent()) {
    auto mc_tracks = this->GetAllMCTruthTracksAsTClonesArray();
    // Avoid reallocation of the vector in the loop
    ret_vector.reserve(mc_tracks->GetEntriesFast());
    TIter next_tr(mc_tracks);
    AliAODMCParticle* mc_tr = 0;
    auto weight = 1;
    while ((mc_tr = static_cast<AliAODMCParticle*>(next_tr()))) {
      if (!mc_tr->IsPrimary()) continue;
      if (mc_tr->Charge() == 0) continue;
      auto tr = AliAnalysisTaskValidation::Track(mc_tr->Eta(), mc_tr->Phi(), mc_tr->Pt(), weight);
      ret_vector.push_back(tr);
    }
  } else { // ESD event
    auto mcEvent = this->MCEvent();
    if (!mcEvent) {
      AliFatal("This is not a Monte Carlo event.");
    }
    Int_t ntracks = mcEvent->GetNumberOfTracks();
    std::cout << "ntracks " << ntracks << std::endl;
    auto valid = 0;
    for (Int_t iTrack=0; iTrack < ntracks; iTrack++) {
      auto mc_p = static_cast<AliMCParticle*>(mcEvent->GetTrack(iTrack));
      if (!mcEvent->Stack()->IsPhysicalPrimary(mc_p->GetLabel())) continue;
      if (mc_p->Charge() == 0) continue;
      auto weight = 1;
      auto tr = AliAnalysisTaskValidation::Track(mc_p->Eta(), mc_p->Phi(), mc_p->Pt(), weight);
      ret_vector.push_back(tr);
      valid++;
    }
    std::cout << "processed " << valid << std::endl;
  }
  return ret_vector;
}

Bool_t AliAnalysisTaskValidation::UserNotify() {
  // If this is MC we have to read all the branches
  // Also, we should check for other tasks here!
  if (!this->MCEvent()) {
    // Turn off all branches
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
  }
  return true;
}
