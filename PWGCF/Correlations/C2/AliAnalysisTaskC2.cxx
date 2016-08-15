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
  
//________________________________________________________________________
AliAnalysisTaskC2::AliAnalysisTaskC2()
  : AliAnalysisTaskC2Base(),
    fEventCounter(0),
    fSingles(0),
    fPairs(0)
{
}

//________________________________________________________________________
AliAnalysisTaskC2::AliAnalysisTaskC2(const char *name)
  : AliAnalysisTaskC2Base(name),
    fEventCounter(0),
    fSingles(0),
    fPairs(0)
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
  // Setup up 2-particle density histogram
  // eta1, eta2, phi1, phi2, mult, zvtx
  Int_t nptPairBins =
    AliAnalysisC2Utils::ComputePtPairBin(this->fSettings.fPtBinEdges.size() - 1,
					 this->fSettings.fPtBinEdges.size() - 1)
    + 1;
  Int_t nMult = this->fSettings.fMultBinEdges.size() - 1;
  Int_t nEta = this->fSettings.fEtaBinEdges.size() - 1;

  Int_t nbins_pairs[cPairsDims::kNdimensions] = {0};
  nbins_pairs[cPairsDims::kEta1] = nEta;
  nbins_pairs[cPairsDims::kEta2] = nEta;
  nbins_pairs[cPairsDims::kPhi1] = this->fSettings.fNPhiBins;
  nbins_pairs[cPairsDims::kPhi2] = this->fSettings.fNPhiBins;
  nbins_pairs[cPairsDims::kPtPair] = nptPairBins;
  nbins_pairs[cPairsDims::kMult] = Int_t(this->fSettings.fMultBinEdges.size() - 1);
  nbins_pairs[cPairsDims::kZvtx] = this->fSettings.fNZvtxBins;

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
  
  this->fPairs = new THnS("pairs",
			  "<N_{1}N_{2}>;#eta_{1};#eta_{2};#phi_{1};#phi_{2};p_{pair};mult;z_{vtx};",
			  cPairsDims::kNdimensions, nbins_pairs, xmin_pairs, xmax_pairs);
  this->fPairs->GetAxis(cPairsDims::kMult)->Set(nMult, &this->fSettings.fMultBinEdges[0]);
  this->fPairs->GetAxis(cPairsDims::kEta1)->Set(nEta, &this->fSettings.fEtaBinEdges[0]);
  this->fPairs->GetAxis(cPairsDims::kEta2)->Set(nEta, &this->fSettings.fEtaBinEdges[0]);
  // exclude_over_under_flow(this->fPairs);
  this->fOutputList->Add(fPairs);

  Int_t nPtbins = this->fSettings.fPtBinEdges.size() - 1;
  // eta, phi, pt, mult, zvtx
  // Reuse values from pairs histogram where appropriate
  Int_t nbins_singles[cSinglesDims::kNdimensions] = {
    nbins_pairs[cPairsDims::kEta1],
    nbins_pairs[cPairsDims::kPhi1],
    nPtbins,
    nbins_pairs[cPairsDims::kMult],
    nbins_pairs[cPairsDims::kZvtx]};
  Double_t xmin_singles[cSinglesDims::kNdimensions] = {
    xmin_pairs[cPairsDims::kEta1],
    xmin_pairs[cPairsDims::kPhi1],
    0,  // Dummy
    xmin_pairs[cPairsDims::kMult],
    xmin_pairs[cPairsDims::kZvtx]};
  Double_t xmax_singles[cSinglesDims::kNdimensions] = {
    xmax_pairs[cPairsDims::kEta1],
    xmax_pairs[cPairsDims::kPhi1],
    1,  // Dummy
    xmax_pairs[cPairsDims::kMult],
    xmax_pairs[cPairsDims::kZvtx]};
  this->fSingles = new THnF("singles", "<N>;#eta;#phi;p_{T};mult;z_{vtx};",
			    cSinglesDims::kNdimensions, nbins_singles, xmin_singles, xmax_singles);
  this->fSingles->GetAxis(cSinglesDims::kEta)->Set(nEta, &this->fSettings.fEtaBinEdges[0]);
  this->fSingles->GetAxis(cSinglesDims::kPt)->Set(nPtbins, &this->fSettings.fPtBinEdges[0]);
  this->fSingles->GetAxis(cSinglesDims::kMult)->Set(nMult, &this->fSettings.fMultBinEdges[0]);
  this->fOutputList->Add(this->fSingles);

  // Event counter
  Int_t nbins_evt_counter[cEventCounterDims::kNdimensions] = {
    nbins_pairs[cPairsDims::kMult],
    nbins_pairs[cPairsDims::kZvtx]
  };
  Double_t xmin_evt_counter[cEventCounterDims::kNdimensions] = {
    xmin_pairs[cPairsDims::kMult],
    xmin_pairs[cPairsDims::kZvtx]
  };
  Double_t xmax_evt_counter[cEventCounterDims::kNdimensions] = {
    xmax_pairs[cPairsDims::kMult],
    xmax_pairs[cPairsDims::kZvtx]
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
  // this->fmultDistribution->Fill(this->GetEventClassifierValue());

  const Double_t weight = (this->fSettings.kMCTRUTH == this->fSettings.fDataType)
    ? mcEvent->GenEventHeader()->EventWeight()
    : 1;
  const Double_t zvtx = (this->InputEvent()->GetPrimaryVertex())
    ? this->InputEvent()->GetPrimaryVertex()->GetZ()
    : -999;

  for (auto const &track: this->GetValidTracks()) {
    // Fill the single track histogram
    Double_t stuffing[5] = {track.eta,
			    AliAnalysisC2Utils::WrapAngle(track.phi,
							  this->fSingles->GetAxis(cSinglesDims::kPhi)),
			    track.pt,
			    multiplicity,
			    zvtx};
    this->fSingles->Fill(stuffing, weight);
  }
  {
    Double_t stuffing[3] = {multiplicity, zvtx};
    this->fEventCounter->Fill(stuffing, weight);
  }
  // n: number of valid tracks
  const std::vector< AliAnalysisC2NanoTrack > tracks = this->GetValidTracks();
  const AliAnalysisC2NanoTrack *trigger, *assoc;
  for (Int_t i = 0; i < tracks.size(); i++) {
    // Do not pair with itself and drop mirrored pairs
    for (Int_t j = i + 1; j < tracks.size(); j++) {
      // assign trigger/assoc depending on pt
      if (tracks[i].pt <= tracks[j].pt) {
	trigger = &tracks[j];
	assoc = &tracks[i];
      }
      else {
	trigger = &tracks[i];
	assoc = &tracks[j];
      }
      Int_t pt1Bin = this->fSingles->GetAxis(cSinglesDims::kPt)->FindFixBin(assoc->pt);
      Int_t pt2Bin = this->fSingles->GetAxis(cSinglesDims::kPt)->FindFixBin(trigger->pt);
      Double_t stuffing[8] =
	{assoc->eta,
	 trigger->eta,
	 AliAnalysisC2Utils::WrapAngle(assoc->phi, this->fPairs->GetAxis(cPairsDims::kPhi1)),
	 AliAnalysisC2Utils::WrapAngle(trigger->phi, this->fPairs->GetAxis(cPairsDims::kPhi2)),
	 Double_t(AliAnalysisC2Utils::ComputePtPairBin(pt1Bin, pt2Bin) + 0.5),  // +.5 to hit the bin center
	 multiplicity,
	 zvtx};
      this->fPairs->Fill(stuffing, weight * weight);
    }
  }
  PostData(1, this->fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskC2::Terminate(Option_t *)
{
  // PostData(1, this->fOutputList);
}
