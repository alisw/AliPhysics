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

#include "AliAnalysisC2Utils.h"
#include "AliAnalysisTaskValidation.h"

using std::cout;
using std::endl;


AliAnalysisTaskValidation::AliAnalysisTaskValidation()
  : AliAnalysisTaskSE(),
    fIsValidEvent(false),
    fValidators(),
    fOutputList(0),
    fQADiscard_flow(0),
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
    fValidators(),
    fOutputList(0),
    fQADiscard_flow(0),
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
  fValidators.push_back(EventValidation::kNoCut);
  fValidators.push_back(EventValidation::kIsAODEvent);
  fValidators.push_back(EventValidation::kPassesAliEventCuts);
  fValidators.push_back(EventValidation::kHasFMD);
  fValidators.push_back(EventValidation::kHasEntriesFMD);
  fValidators.push_back(EventValidation::kHasEntriesV0);
  fValidators.push_back(EventValidation::kHasValidVertex);
  fValidators.push_back(EventValidation::kHasMultSelection);
  // fValidators.push_back(EventValidation::kNotOutOfBunchPU);
  fValidators.push_back(EventValidation::kNotMultiVertexPU);
  fValidators.push_back(EventValidation::kNotSPDPU);
  fValidators.push_back(EventValidation::kNotSPDClusterVsTrackletBG);
  fValidators.push_back(EventValidation::kPassesFMD_V0CorrelatioCut);

  // Define output slot
  DefineOutput(1, TList::Class());
  DefineOutput(2, this->Class());

  // Enable mulivertex pileup cuts
  fEventCuts.fPileUpCutMV = true;
}

void AliAnalysisTaskValidation::CreateQAHistograms(TList* outlist) {
  this->fQADiscard_flow = new TH1F("qa_discard_flow", "QA event discard flow",
				   this->fValidators.size(), 0, this->fValidators.size());
  TAxis *discardedEvtsAx = this->fQADiscard_flow->GetXaxis();

  for (UInt_t idx = 0; idx < this->fValidators.size(); idx++) {
    switch (this->fValidators[idx]) {
    case EventValidation::kNoCut:
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
  outlist->Add(this->fQADiscard_flow);
}

void AliAnalysisTaskValidation::UserCreateOutputObjects() {
  // Stop right here if there are no Validators to work with
  if (this->fValidators.size() == 0) {
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
  for (UInt_t idx = 0; idx < this->fValidators.size(); idx++) {
    switch (this->fValidators[idx]) {
    case EventValidation::kNoCut:
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
      this->fQADiscard_flow->Fill(idx);
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
  if (nV0A_hits + nV0C_hits < (nFMD_fwd_hits + nFMD_bwd_hits - 40)) {
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
  
  AliVVZERO* vzero = static_cast<AliAODEvent*>(this->InputEvent())->GetVZEROData();
  const Int_t nV0_channels = 64;
  // V0 has no pt resolution!
  Double_t pt = 0;
  for (Int_t ich = 0; ich < nV0_channels; ich++) {
    //vzero->GetMultiplicity(ich);
    Float_t amp = static_cast<AliAODEvent*>(this->InputEvent())->GetVZEROEqMultiplicity(ich);
    if (amp > 0) {
      Float_t eta = 0.5 * (vzero->GetVZEROEtaMin(ich) + vzero->GetVZEROEtaMax(ich));
      Float_t phi = vzero->GetVZEROAvgPhi(ich);
      ret_vector.push_back(
	AliAnalysisTaskValidation::Track(eta, AliAnalysisC2Utils::Wrap02pi(phi), pt, amp));
    }
  }
  return ret_vector;
}

AliAnalysisTaskValidation::Tracks AliAnalysisTaskValidation::GetSPDtracklets() const {
  AliAODEvent* aodEvent = dynamic_cast< AliAODEvent* >(this->InputEvent());
  AliAODTracklets* aodTracklets = aodEvent->GetTracklets();

  const Double_t pt = 0;
  AliAnalysisTaskValidation::Tracks ret_vector;
  for (Int_t i = 0; i < aodTracklets->GetNumberOfTracklets(); i++) {
    auto eta = aodTracklets->GetEta(i);
    auto phi = aodTracklets->GetPhi(i);
    ret_vector.push_back(AliAnalysisTaskValidation::Track(eta, phi, pt, 1));
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

