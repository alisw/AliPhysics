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
    fPairs(0),
    fmultDistribution(0),
    fNchDistribution(0)
{
}

//________________________________________________________________________
AliAnalysisTaskC2::AliAnalysisTaskC2(const char *name)
  : AliAnalysisTaskC2Base(name),
    fEventCounter(0),
    fSingles(0),
    fPairs(0),
    fmultDistribution(0),
    fNchDistribution(0)
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
  Int_t nbins_pairs[cPairsDims::kNdimensions] = {
    this->fSettings.fNEtaBins, this->fSettings.fNEtaBins,
    this->fSettings.fNPhiBins, this->fSettings.fNPhiBins,
    nptPairBins,
    Int_t(this->fSettings.fMultBinEdges.size() - 1),
    this->fSettings.fNZvtxBins};
  Double_t xmin_pairs[cPairsDims::kNdimensions] = {
    this->fSettings.fEtaAcceptanceLowEdge, this->fSettings.fEtaAcceptanceLowEdge,
    this->fSettings.fPhiAcceptanceLowEdge, this->fSettings.fPhiAcceptanceLowEdge,
    0,  // ptpairbin starts at 0 (OBS: this it not GeV!)
    0,  // mult dummy
    this->fSettings.fZVtxAcceptanceLowEdge
    };
  Double_t xmax_pairs[cPairsDims::kNdimensions] = {
    this->fSettings.fEtaAcceptanceUpEdge, this->fSettings.fEtaAcceptanceUpEdge,
    this->fSettings.fPhiAcceptanceUpEdge, this->fSettings.fPhiAcceptanceUpEdge,
    Double_t(nptPairBins),  // ptpairbin max
    1,  // mult dummy
    this->fSettings.fZVtxAcceptanceUpEdge};
  this->fPairs = new THnS("pairs",
			  "<N_{1}N_{2}>;#eta_{1};#eta_{2};#phi_{1};#phi_{2};p_{pair};mult;z_{vtx};",
			  cPairsDims::kNdimensions, nbins_pairs, xmin_pairs, xmax_pairs);
  this->fPairs->GetAxis(cPairsDims::kMult)
    ->Set(nMult, &this->fSettings.fMultBinEdges[0]);
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

  this->fNchDistribution = new TProfile2D("NchDistribution", "Mean N_{ch};cent;z_{vtx};",
					  nbins_pairs[cPairsDims::kMult],
					  xmin_pairs[cPairsDims::kMult],
					  xmax_pairs[cPairsDims::kMult],
					  nbins_pairs[cPairsDims::kZvtx],
					  xmin_pairs[cPairsDims::kZvtx],
					  xmax_pairs[cPairsDims::kZvtx]);
  this->fNchDistribution->GetXaxis()->Set(nbins_pairs[cPairsDims::kMult], &this->fSettings.fMultBinEdges[0]);
  this->fOutputList->Add(this->fNchDistribution);

  // Debug hists
  this->fmultDistribution = new TH1F("multDistribution", "multDistribution", 210, 0, 210);
  this->fOutputList->Add(this->fmultDistribution);

  // this->fRndmGenerator = new TRandom3();
  
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
  // Event is invalid if no multselection is present; ie. tested in base class
  AliMultEstimator *multEstimator =
    (dynamic_cast< AliMultSelection* >(this->InputEvent()->FindListObject("MultSelection")))
    ->GetEstimator(this->fSettings.fMultEstimator);
  const Float_t multiplicity = multEstimator->GetPercentile();
  this->fmultDistribution->Fill(multiplicity);

  const Double_t weight = (this->fSettings.kMCTRUTH == this->fSettings.fDataType)
    ? mcEvent->GenEventHeader()->EventWeight()
    : 1;
  const Double_t zvtx = (this->InputEvent()->GetPrimaryVertex())
    ? this->InputEvent()->GetPrimaryVertex()->GetZ()
    : -999;

  this->fNchDistribution->Fill(multiplicity, zvtx, multEstimator->GetValue());
  // Yes, the following is realy "aodEvent" not mcEvent :P
  TClonesArray* tracksArray = (this->fSettings.kMCTRUTH == this->fSettings.fDataType)
    ? dynamic_cast<TClonesArray*>(aodEvent->GetList()->FindObject(AliAODMCParticle::StdBranchName()))
    : aodEvent->GetTracks();

  if (!tracksArray || tracksArray->GetSize()==0){
    AliWarning("Could not retrieve track array");
    PostData(1, this->fOutputList);
    return;
  }

  std::vector< cNano_track > tracks;
  TIter nextTrack(tracksArray);
  while (TObject* obj = nextTrack()){
    // The naming around VTrack, VParticle, AODTrack mcParticle is a mess!
    // Take-away message: They all derive from AliVParticle one way or another.
    AliVParticle* particle = dynamic_cast< AliVParticle* >(obj);
    if (!this->IsValidParticle(particle)){
      continue;
    }
    // eta, phi, pt, q are pure virtual in AliVParticle, so it should be safe to just call them on AliVParticle
    cNano_track tmp_track = {
      particle->Eta(), //fRndmGenerator->Rndm() * 1.6 - 0.8;
      particle->Phi(), //fRndmGenerator->Rndm() * 2*TMath::Pi();
      particle->Pt()
    };

    // check if this is a under/overflow bin. If so, discard the track
    Int_t ptBinIndex = this->fSingles->GetAxis(cSinglesDims::kPt)->FindFixBin(tmp_track.pt);
    if (ptBinIndex == 0 || ptBinIndex == this->fSingles->GetAxis(cSinglesDims::kPt)->GetNbins() + 1)
      continue;
    
    // Fill the single track histogram
    Double_t stuffing[5] = {tmp_track.eta,
			    AliAnalysisC2Utils::WrapAngle(tmp_track.phi,
							  this->fSingles->GetAxis(cSinglesDims::kPhi)),
			    tmp_track.pt,
			    multiplicity,
			    zvtx};
    this->fSingles->Fill(stuffing, weight);
    // Save the track for the pair histogram
    tracks.push_back(tmp_track);
  }
  // If we have at least one valid track, it is a valid event
  // If there are no tracks, do not count the event since this would be an undefined correlation
  if (tracks.size() == 0){
    // This might be slightly inconsistent and misses the PostData!
    this->fDiscardedEvents->Fill(cDiscardEventReasons::noTracksInPtRegion);
    PostData(1, this->fOutputList);
    return;
  }
  {
    Double_t stuffing[3] = {multiplicity, zvtx};
    this->fEventCounter->Fill(stuffing, weight);
  }
  // n: number of valid tracks
  for (std::vector<cNano_track>::size_type iTrack = 0; iTrack < tracks.size(); iTrack++) {
      for (std::vector<cNano_track>::size_type jTrack = 0; jTrack < tracks.size(); jTrack++) {
	if (iTrack == jTrack)
	  continue;
	// gurantee that pt1 is always smaller than pt2. This takes care of mirrored pairs
	if (tracks[iTrack].pt > tracks[jTrack].pt)
	  continue;
	Int_t pt1Bin = this->fSingles->GetAxis(cSinglesDims::kPt)->FindFixBin(tracks[iTrack].pt);
	Int_t pt2Bin = this->fSingles->GetAxis(cSinglesDims::kPt)->FindFixBin(tracks[jTrack].pt);
	Double_t stuffing[8] =
	  {tracks[iTrack].eta,
	   tracks[jTrack].eta,
	   AliAnalysisC2Utils::WrapAngle(tracks[iTrack].phi, this->fPairs->GetAxis(cPairsDims::kPhi1)),
	   AliAnalysisC2Utils::WrapAngle(tracks[jTrack].phi, this->fPairs->GetAxis(cPairsDims::kPhi2)),
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
