#include <iostream>
#include <vector>

#include "TH1.h"
#include "THn.h"
#include "TList.h"
#include "TMath.h"

#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliMCEvent.h"
#include "AliMultSelection.h"
#include "AliStack.h"
#include "AliVEvent.h"

#include "AliAnalysisTaskC2.h"
#include "AliAnalysisC2Utils.h"

using std::cout;
using std::endl;
using std::vector;

//________________________________________________________________________
AliAnalysisTaskC2::AliAnalysisTaskC2()
  : AliAnalysisTaskSE(),
    fSettings(),
    fOutputList(0),
    fEventCounter(0),
    fEtaPhiZvtx_max_res(0),
    fsingleHists(0),
    fpairHists(0)
{
}

//________________________________________________________________________
AliAnalysisTaskC2::AliAnalysisTaskC2(const char *name)
  : AliAnalysisTaskSE(name),
    fSettings(),
    fOutputList(0),
    fEventCounter(0),
    fEtaPhiZvtx_max_res(0),
    fsingleHists(0),
    fpairHists(0)
{
  // Rely on validation task for event and track selection
  DefineInput(1, AliAnalysisTaskValidation::Class());
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskC2::UserCreateOutputObjects()
{
  // Setup output list if it was not done already
  if (!this->fOutputList){
    this->fOutputList = new TList();
    this->fOutputList->SetOwner();
  }

  // Some values which are used throughout the hist setup
  const Int_t nMult = this->fSettings.fMultBinEdges.size() - 1;
  const Int_t nPtbins = this->fSettings.fPtBinEdges.size() - 1;
  const Int_t nptPairBins = AliAnalysisC2Utils::ComputePtPairBin(nPtbins, nPtbins) + 1;
  const Int_t nZvtxBins = this->fSettings.fNZvtxBins;

  // Setup up 2-particle density histogram eta1, eta2, phi1, phi2, mult, zvtx
  auto setup_pair_hist =
    [this, nMult, nptPairBins, nZvtxBins] (const char* name, vector<Double_t> eta1Edges, vector<Double_t> eta2Edges)
    {
      Int_t nbins_pairs[cPairsDims::kNdimensions] = {0};
      nbins_pairs[cPairsDims::kEta1] = Int_t(eta1Edges.size() - 1);
      nbins_pairs[cPairsDims::kEta2] = Int_t(eta2Edges.size() - 1);
      nbins_pairs[cPairsDims::kPhi1] = this->fSettings.fNPhiBins;
      nbins_pairs[cPairsDims::kPhi2] = this->fSettings.fNPhiBins;
      nbins_pairs[cPairsDims::kPtPair] = nptPairBins;
      nbins_pairs[cPairsDims::kMult] = Int_t(this->fSettings.fMultBinEdges.size() - 1);
      nbins_pairs[cPairsDims::kZvtx] = nZvtxBins;

      Double_t xmin_pairs[cPairsDims::kNdimensions] = {0};
      xmin_pairs[cPairsDims::kEta1]   = 0; // Dummy
      xmin_pairs[cPairsDims::kEta2]   = 0; // Dummy
      xmin_pairs[cPairsDims::kPhi1]   = this->fSettings.fPhiAcceptanceLowEdge;
      xmin_pairs[cPairsDims::kPhi2]   = this->fSettings.fPhiAcceptanceLowEdge;
      xmin_pairs[cPairsDims::kPtPair] = 0;  // ptpairbin starts at 0 (OBS: this it not GeV!)
      xmin_pairs[cPairsDims::kMult]   = 0;    // mult dummy
      xmin_pairs[cPairsDims::kZvtx]   = this->fSettings.fZVtxAcceptanceLowEdge;
  
      Double_t xmax_pairs[cPairsDims::kNdimensions] = {0};
      xmax_pairs[cPairsDims::kEta1]   =  1;  // Dummy
      xmax_pairs[cPairsDims::kEta2]   =  1;  // Dummy
      xmax_pairs[cPairsDims::kPhi1]   =  this->fSettings.fPhiAcceptanceUpEdge;
      xmax_pairs[cPairsDims::kPhi2]   =  this->fSettings.fPhiAcceptanceUpEdge;
      xmax_pairs[cPairsDims::kPtPair] =  Double_t(nptPairBins);  // ptpairbin max
      xmax_pairs[cPairsDims::kMult]   =  1;  // Dummy
      xmax_pairs[cPairsDims::kZvtx]   =  this->fSettings.fZVtxAcceptanceUpEdge;

      THn* pair = new THnF(name,
			   "<N_{1}N_{2}>;#eta_{1};#eta_{2};#phi_{1};#phi_{2};p_{pair};mult;z_{vtx};",
			   cPairsDims::kNdimensions, nbins_pairs, xmin_pairs, xmax_pairs);
      pair->GetAxis(cPairsDims::kMult)->Set(nMult, &this->fSettings.fMultBinEdges[0]);
      pair->GetAxis(cPairsDims::kEta1)->Set(eta1Edges.size() - 1, &eta1Edges[0]);
      pair->GetAxis(cPairsDims::kEta2)->Set(eta2Edges.size() - 1, &eta2Edges[0]);

      this->fpairHists.push_back(pair);
      this->fOutputList->Add(pair);
    };
  setup_pair_hist("pairs_bwdbwd", this->fSettings.fEtaEdgesBwd, this->fSettings.fEtaEdgesBwd);
  setup_pair_hist("pairs_bwdits", this->fSettings.fEtaEdgesBwd, this->fSettings.fEtaEdgesIts);
  setup_pair_hist("pairs_bwdfwd", this->fSettings.fEtaEdgesBwd, this->fSettings.fEtaEdgesFwd);
  setup_pair_hist("pairs_itsits", this->fSettings.fEtaEdgesIts, this->fSettings.fEtaEdgesIts);
  setup_pair_hist("pairs_itsfwd", this->fSettings.fEtaEdgesIts, this->fSettings.fEtaEdgesFwd);
  setup_pair_hist("pairs_fwdfwd", this->fSettings.fEtaEdgesFwd, this->fSettings.fEtaEdgesFwd);

  // Single histogram: eta, phi, pt, mult, zvtx
  auto setup_single_hist =
    [this, nMult, nPtbins, nZvtxBins] (const char* name, vector<Double_t> etaEdges)
    {
      Int_t nbins_singles[cSinglesDims::kNdimensions] = {
	Int_t(etaEdges.size() - 1),
	this->fSettings.fNPhiBins,
	nPtbins,
        nMult,
	nZvtxBins
      };
      Double_t xmin_singles[cSinglesDims::kNdimensions] = {
	0,  // Dummy; eta min
	this->fSettings.fPhiAcceptanceLowEdge,
	0,  // Dummy; pt min
	0,  // Dummy; mult min
        this->fSettings.fZVtxAcceptanceLowEdge
      };
      Double_t xmax_singles[cSinglesDims::kNdimensions] = {
	1, // Dummy eta max
	this->fSettings.fPhiAcceptanceUpEdge,
	1,  // Dummy pt max
	1,  // Dummy mult max
	this->fSettings.fZVtxAcceptanceUpEdge
      };
      THn* single = new THnF(name, "<N>;#eta;#phi;p_{T};mult;z_{vtx};",
			     cSinglesDims::kNdimensions, nbins_singles, xmin_singles, xmax_singles);
      single->GetAxis(cSinglesDims::kEta)->Set(etaEdges.size() - 1, &etaEdges[0]);
      single->GetAxis(cSinglesDims::kPt)->Set(nPtbins, &this->fSettings.fPtBinEdges[0]);
      single->GetAxis(cSinglesDims::kMult)->Set(nMult, &this->fSettings.fMultBinEdges[0]);

      this->fsingleHists.push_back(single);
      this->fOutputList->Add(single);
    };
  setup_single_hist("singles_bwd", this->fSettings.fEtaEdgesBwd);
  setup_single_hist("singles_its", this->fSettings.fEtaEdgesIts);
  setup_single_hist("singles_fwd", this->fSettings.fEtaEdgesFwd);
  
  // Event counter
  Int_t nbins_evt_counter[cEventCounterDims::kNdimensions] = {
    nMult,
    nZvtxBins
  };
  Double_t xmin_evt_counter[cEventCounterDims::kNdimensions] = {
    0, // Dummy mult min
    this->fSettings.fZVtxAcceptanceLowEdge
  };
  Double_t xmax_evt_counter[cEventCounterDims::kNdimensions] = {
    0, // Dummy mult max
    this->fSettings.fZVtxAcceptanceUpEdge
  };
  this->fEventCounter = new THnF("eventCounter", "Event counter;mult;z_{vtx};",
				 cEventCounterDims::kNdimensions,
				 nbins_evt_counter, xmin_evt_counter, xmax_evt_counter);
  this->fEventCounter->GetAxis(cEventCounterDims::kMult)->Set(nMult, &this->fSettings.fMultBinEdges[0]);
  this->fOutputList->Add(this->fEventCounter);

  {
    const Int_t ndims = 3;
    const Int_t nbins[ndims] = {200, 20, 100};
    const Double_t lowerBounds[ndims] = {-4.0, 0, -10.0};
    const Double_t upperBounds[ndims] = {6.0, 2*TMath::Pi(), 10.0};
    this->fEtaPhiZvtx_max_res = new THnF("etaPhiZvtx_max_res",
					 "etaPhiZvtx_max_res;#eta;#phi;#zvtx;",
					 ndims,
					 nbins,
					 lowerBounds,
					 upperBounds
					 );
    this->fOutputList->Add(this->fEtaPhiZvtx_max_res);
  }
  
  AliLog::SetGlobalLogLevel(AliLog::kError);
  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskC2::UserExec(Option_t *)
{
  // Get the event validation object
  AliAnalysisTaskValidation* ev_val = dynamic_cast<AliAnalysisTaskValidation*>(this->GetInputData(1));

  AliMCEvent*   mcEvent = this->MCEvent();

  if (!ev_val->IsValidEvent()){
    PostData(1, this->fOutputList);
    return;
  }
  const Float_t multiplicity = this->GetEventClassifierValue();
  const Double_t zvtx = (this->InputEvent()->GetPrimaryVertex())
    ? this->InputEvent()->GetPrimaryVertex()->GetZ()
    : -999;
  // Check that the given event is not over/under in z_vtx or mult
  {
    TAxis* multAxis = fEventCounter->GetAxis(cEventCounterDims::kMult);
    if (multAxis->FindBin(multiplicity) == 0
      || multAxis->FindBin(multiplicity) == multAxis->GetNbins() + 1) {
      return;
    }
  }
  {
    TAxis* zvtxAxis = fEventCounter->GetAxis(cEventCounterDims::kZvtx);
    if (zvtxAxis->FindBin(zvtx) == 0
      || zvtxAxis->FindBin(zvtx) == zvtxAxis->GetNbins() + 1) {
      return;
    }
  }
  const Double_t evWeight = (this->fSettings.kMCTRUTH == this->fSettings.fDataType)
    ? mcEvent->GenEventHeader()->EventWeight()
    : 1;
  // Load all valid tracks/hits used in the following
  const auto tracks = this->GetValidTracks();

  // fill the valid tracks into the appropriate QA histograms
  for (auto &t: tracks) {
    const Double_t stuffing[] = {t.eta, t.phi, zvtx};
    this->fEtaPhiZvtx_max_res->Fill(stuffing, t.weight);
  }

  // Fill the event counter histogram
  {
    Double_t stuffing[3] = {multiplicity, zvtx};
    this->fEventCounter->Fill(stuffing, evWeight);
  }
  
  // Fill the single particle histograms
  for (auto const &track: tracks) {
    // Fill the single track histograms
    Double_t stuffing[5] = {track.eta,
			    AliAnalysisC2Utils::Wrap02pi(track.phi),
			    track.pt,
			    multiplicity,
			    zvtx};
    for (auto single: this->fsingleHists) {
      single->Fill(stuffing, evWeight * track.weight);
    }
  }
  struct Bin_indices{
    Int_t i_eta;
    Int_t i_phi;
    Int_t i_pt;
    Int_t i_mult;
    Int_t i_zvtx;
    AliAnalysisTaskValidation::Track track;
  };
  // Find the bin indices here instead of doing it (implicity in the loop) which is super slow
  vector<Bin_indices> track_bin_idxs;
  for (auto t: tracks) {
    THn *hist;
    if (t.eta < this->fSettings.fEtaEdgesIts[0]) {
      hist = this->fsingleHists[0];
    } else if (t.eta >= this->fSettings.fEtaEdgesIts[0] && t.eta < this->fSettings.fEtaEdgesFwd[0]) {
      hist = this->fsingleHists[1];
    } else {
      hist = this->fsingleHists[2];
    }
    track_bin_idxs.push_back(Bin_indices {
	hist->GetAxis(cSinglesDims::kEta)->FindFixBin(t.eta),
	hist->GetAxis(cSinglesDims::kPhi)->FindFixBin(AliAnalysisC2Utils::Wrap02pi(t.phi)),
	hist->GetAxis(cSinglesDims::kPt)->FindFixBin(t.pt),
	hist->GetAxis(cSinglesDims::kMult)->FindFixBin(multiplicity),
	hist->GetAxis(cSinglesDims::kZvtx)->FindFixBin(zvtx),
	t});
  }
  // Fill the pair particle histograms
  Bin_indices *trigger, *assoc;
  for (UInt_t i = 0; i < track_bin_idxs.size(); i++) {
    // Do not pair with itself and drop mirrored pairs
    for (UInt_t j = i + 1; j < track_bin_idxs.size(); j++) {
      // assign trigger/assoc
      // If we bin in pt, use pt as distinction, else use eta
      if (this->fSettings.fPtBinEdges.size() > 2) {
	if (track_bin_idxs[i].track.pt <= track_bin_idxs[j].track.pt) {
	  trigger = &track_bin_idxs[j];
	  assoc = &track_bin_idxs[i];
	}
	else {
	  trigger = &track_bin_idxs[i];
	  assoc = &track_bin_idxs[j];
	}
      }
      else {
	if (track_bin_idxs[i].track.eta <= track_bin_idxs[j].track.eta) {
	  trigger = &track_bin_idxs[j];
	  assoc = &track_bin_idxs[i];
	}
	else {
	  trigger = &track_bin_idxs[i];
	  assoc = &track_bin_idxs[j];
	}
      }
      // FIXME: This will break if I bin ITS region in pt but not the FMD part!
      Int_t pt1Bin = 1; // this->fsingleHists[0]->GetAxis(cSinglesDims::kPt)->FindFixBin(assoc->pt);
      Int_t pt2Bin = 1; // this->fsingleHists[0]->GetAxis(cSinglesDims::kPt)->FindFixBin(trigger->pt);
      Int_t stuffing[7] =
	{assoc->i_eta,
	 trigger->i_eta,
	 assoc->i_phi,
	 trigger->i_phi,
	 AliAnalysisC2Utils::ComputePtPairBin(pt1Bin, pt2Bin) + 1, // bin indices start at 1...
	 assoc->i_mult,
	 assoc->i_zvtx};
      this->fpairHists[this->GetPairHistIndex(trigger->track.eta, assoc->track.eta)]
	->AddBinContent(stuffing, assoc->track.weight * trigger->track.weight * evWeight);
    }
  }
  PostData(1, this->fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskC2::Terminate(Option_t *)
{
  // PostData(1, this->fOutputList);
}

UInt_t AliAnalysisTaskC2::GetPairHistIndex(Float_t trigger, Float_t assoc) {
  auto is_bwd =
    [this] (Float_t eta)
    {return eta < this->fSettings.fEtaEdgesIts[0];};
  auto is_its =
    [this] (Float_t eta)
    {return eta >= this->fSettings.fEtaEdgesIts[0] && eta < this->fSettings.fEtaEdgesFwd[0];};
  auto is_fwd =
    [this] (Float_t eta)
    {return eta >= this->fSettings.fEtaEdgesFwd[0];};
  if (is_bwd(assoc) && is_bwd(trigger)) return 0; // bwd bwd
  if (is_bwd(assoc) && is_its(trigger)) return 1; // bwd its
  if (is_bwd(assoc) && is_fwd(trigger)) return 2; // bwd fwd
  if (is_its(assoc) && is_its(trigger)) return 3; // its its
  if (is_its(assoc) && is_fwd(trigger)) return 4; // its fwd
  if (is_fwd(assoc) && is_fwd(trigger)) return 5; // fwd fwd
  return -1;
}


AliAnalysisTaskValidation::Tracks AliAnalysisTaskC2::GetValidTracks() {
  // Get the event validation object
  AliAnalysisTaskValidation* ev_val = dynamic_cast<AliAnalysisTaskValidation*>(this->GetInputData(1));

  AliAnalysisTaskValidation::Tracks ret_vector;
  // Append central tracklets
  if (this->fSettings.kRECON == this->fSettings.fDataType) {
    // Are we running on SPD clusters? If so add them to our track vector
    if (this->fSettings.fUseSPDclusters) {
      AliError("SPD clusters not yet implemented");
    }
    else if (this->fSettings.fUseSPDtracklets) {
      auto spdhits = ev_val->GetSPDtracklets();
      ret_vector.insert(ret_vector.end(), spdhits.begin(), spdhits.end());
    }
    // Append the fmd hits to this vector if we are looking at reconstructed data,
    // All hits on the FMD (above the internally used threshold) are "valid"
    if (this->fSettings.fUseFMD) {
      auto fmdhits = ev_val->GetFMDhits();
      ret_vector.insert(ret_vector.end(), fmdhits.begin(), fmdhits.end());
    } else if (this->fSettings.fUseV0){
      auto v0amps = ev_val->GetV0hits();
      ret_vector.insert(ret_vector.end(), v0amps.begin(), v0amps.end());
    }
  }
  // MC truth case:
  else {
    AliError("MC truth is not yet implemented");
    // auto allMCtracks = this->GetValidCentralTracks();
    // this->fValidTracks.insert(this->fValidTracks.end(), allMCtracks.begin(), allMCtracks.end());
  }
  return ret_vector;
}

Float_t AliAnalysisTaskC2::GetEventClassifierValue() {
  if (this->fSettings.fMultEstimator == this->fSettings.fMultEstimatorValidTracks){
    return this->GetValidTracks().size();
  }
  else {
    // Event is invalid if no multselection is present; ie. tested in IsValidEvent() already
    AliMultEstimator *multEstimator =
      (dynamic_cast< AliMultSelection* >(this->InputEvent()->FindListObject("MultSelection")))
      ->GetEstimator(this->fSettings.fMultEstimator);
    // const Float_t multiplicity = ((Float_t)multEstimator->GetValue()) / multEstimator->GetMean();
    const Float_t multiplicity = multEstimator->GetPercentile();
    return multiplicity;
  }
}
