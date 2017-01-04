#include <iostream>
#include <vector>

#include "TH1.h"
#include "THn.h"
#include "TList.h"
#include "TMath.h"
#include "TProfile2D.h"

#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliLog.h"
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
  : AliAnalysisTaskC2Base(),
    fEventCounter(0),
    fsingleHists(0),
    fpairHists(0)
{
}

//________________________________________________________________________
AliAnalysisTaskC2::AliAnalysisTaskC2(const char *name)
  : AliAnalysisTaskC2Base(name),
    fEventCounter(0),
    fsingleHists(0),
    fpairHists(0)
{
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskC2::UserCreateOutputObjects()
{
  // Call the base class to set up histograms there
  AliAnalysisTaskC2Base::UserCreateOutputObjects();
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

  AliLog::SetGlobalLogLevel(AliLog::kError);
  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskC2::UserExec(Option_t *)
{
  // Run setup for base class
  this->SetupEventForBase();

  // AODEvent() returns a null pointer for whatever reason...
  AliMCEvent*   mcEvent = this->MCEvent();
  AliAODEvent* aodEvent = dynamic_cast< AliAODEvent* >(this->InputEvent());

  if (!this->IsValidEvent()){
    PostData(1, this->fOutputList);
    return;
  }
  const Float_t multiplicity = this->GetEventClassifierValue();
  this->fmultDistribution->Fill(this->GetEventClassifierValue());

  const Double_t evWeight = (this->fSettings.kMCTRUTH == this->fSettings.fDataType)
    ? mcEvent->GenEventHeader()->EventWeight()
    : 1;
  const Double_t zvtx = (this->InputEvent()->GetPrimaryVertex())
    ? this->InputEvent()->GetPrimaryVertex()->GetZ()
    : -999;

  for (auto const &track: this->GetValidTracks()) {
    // Fill the single track histograms
    Double_t stuffing[5] = {track.eta,
			    AliAnalysisC2Utils::WrapAngle(track.phi, fsingleHists[0]->GetAxis(cSinglesDims::kPhi)),
			    track.pt,
			    multiplicity,
			    zvtx};
    for (auto single: this->fsingleHists) {
      single->Fill(stuffing, evWeight * track.weight);
    }
  }
  {
    Double_t stuffing[3] = {multiplicity, zvtx};
    this->fEventCounter->Fill(stuffing, evWeight);
  }
  // n: number of valid tracks
  const vector< AliAnalysisC2NanoTrack > tracks = this->GetValidTracks();
  const AliAnalysisC2NanoTrack *trigger, *assoc;
  for (Int_t i = 0; i < tracks.size(); i++) {
    // Do not pair with itself and drop mirrored pairs
    for (Int_t j = i + 1; j < tracks.size(); j++) {
      // assign trigger/assoc
      // If we bin in pt, use pt as distinction, else use eta
      if (this->fSettings.fPtBinEdges.size() > 2) {
	if (tracks[i].pt <= tracks[j].pt) {
	  trigger = &tracks[j];
	  assoc = &tracks[i];
	}
	else {
	  trigger = &tracks[i];
	  assoc = &tracks[j];
	}
      }
      else {
	if (tracks[i].eta <= tracks[j].eta) {
	  trigger = &tracks[j];
	  assoc = &tracks[i];
	}
	else {
	  trigger = &tracks[i];
	  assoc = &tracks[j];
	}
      }
      // FIXME: This will break if I bin ITS region in pt but not the FMD part!
      Int_t pt1Bin = this->fsingleHists[0]->GetAxis(cSinglesDims::kPt)->FindFixBin(assoc->pt);
      Int_t pt2Bin = this->fsingleHists[0]->GetAxis(cSinglesDims::kPt)->FindFixBin(trigger->pt);
      Double_t stuffing[7] =
	{assoc->eta,
	 trigger->eta,
	 AliAnalysisC2Utils::WrapAngle(assoc->phi, fpairHists[0]->GetAxis(cPairsDims::kPhi1)),
	 AliAnalysisC2Utils::WrapAngle(trigger->phi, fpairHists[0]->GetAxis(cPairsDims::kPhi2)),
	 Double_t(AliAnalysisC2Utils::ComputePtPairBin(pt1Bin, pt2Bin) + 0.5),  // +.5 to hit the bin center
	 multiplicity,
	 zvtx};
      for (auto pair: this->fpairHists) {
	// Only envoc fill if we are in the right eta region. If we
	// filled, move to the next track Ie. Asume that histograms do
	// not overlap in phace-space so each pair can only go into one
	// histogram.  This is so ugly because ROOT does not report
	// back if the fill was successful or in an overflow...
	if (// Check bounderies for eta1
	    pair->GetAxis(cPairsDims::kEta1)->GetBinLowEdge(1) <= stuffing[cPairsDims::kEta1]
	    && stuffing[cPairsDims::kEta1] < pair->GetAxis(cPairsDims::kEta1)->
	    GetBinUpEdge(pair->GetAxis(cPairsDims::kEta1)->GetNbins())
	    // Check bounderies for eta1
	    && pair->GetAxis(cPairsDims::kEta2)->GetBinLowEdge(1) <= stuffing[cPairsDims::kEta2]
	    && stuffing[cPairsDims::kEta2] < pair->GetAxis(cPairsDims::kEta2)->
	    GetBinUpEdge(pair->GetAxis(cPairsDims::kEta2)->GetNbins())) {
	  pair->Fill(stuffing, assoc->weight * trigger->weight * evWeight);
	  continue;
	}
      }
    }
    PostData(1, this->fOutputList);
  }
}
//________________________________________________________________________
void AliAnalysisTaskC2::Terminate(Option_t *)
{
  // PostData(1, this->fOutputList);
}
