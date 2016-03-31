#include <iostream>
#include <vector>

#include "THnSparse.h"
#include "THn.h"
#include "TList.h"
#include "TMath.h"
#include "TVectorF.h"
#include "TMatrix.h"
#include "TRandom3.h"

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
    fcents(0),
    fRndmGenerator(0)
{
}

//________________________________________________________________________
AliAnalysisTaskC2::AliAnalysisTaskC2(const char *name, Int_t mode)
  : AliAnalysisTaskC2Base(name, mode),
    fEventCounter(0),
    fSingles(0),
    fPairs(0),
    fcents(0),
    fRndmGenerator(0)
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
  // eta1, eta2, phi1, phi2, cent, zvtx
  const Int_t etaNbins = 15;
  // phiNbins must be divisable by 2, but not by 4; so that we can later shift it by pi/2 (2pi is total int.)
  // The idea is to have the deltaPhi histogram with a bin centered arround 0
  const Int_t phiNbins = 26;
  const Double_t pt_bin_edges[] = {3.0, 4.0, 6.0, 8.0, 15.0};
  const Int_t npTbins = sizeof(pt_bin_edges) / sizeof(pt_bin_edges[0]) - 1;
  const Int_t nptPairBins = AliAnalysisC2Utils::ComputePtPairBin(npTbins, npTbins) + 1;
  const Double_t cent_bin_edges[] = {0, 20, 40, 90};
  const Int_t nCent    = sizeof(cent_bin_edges)/sizeof(cent_bin_edges[0]) - 1;
  const Int_t nZvtx    = 20;
  const Int_t nbins[cPairsDims::kNdimensions] = {etaNbins, etaNbins,
						 phiNbins, phiNbins,
						 nptPairBins,
						 nCent, nZvtx};
  const Double_t xmin[cPairsDims::kNdimensions] = {-0.8, -0.8,
						   0, 0,
						   0,
						   0, -10};
  const Double_t xmax[cPairsDims::kNdimensions] = {0.8, 0.8,
						   2*TMath::Pi(), 2*TMath::Pi(),
						   Double_t(nptPairBins),
						   100, 10};
  this->fPairs = new THnC("pairs",
			  "<N_{1}N_{2}>;#eta_{1};#eta_{2};#phi_{1};#phi_{2};p_{pair};cent;z_{vtx};",
			  cPairsDims::kNdimensions, nbins, xmin, xmax);
  this->fPairs->GetAxis(cPairsDims::kCent)->Set(nCent, cent_bin_edges);
  // exclude_over_under_flow(this->fPairs);
  this->fOutputList->Add(fPairs);

  {
    // eta, phi, pt, cent, zvtx
    // Reuse values from pairs histogram where appropriate
    Int_t nbins_singles[cSinglesDims::kNdimensions] = {nbins[cPairsDims::kEta1],
						       nbins[cPairsDims::kPhi1],
						       npTbins,
						       nbins[cPairsDims::kCent],
						       nbins[cPairsDims::kZvtx]};
    Double_t xmin_singles[cSinglesDims::kNdimensions] = {xmin[cPairsDims::kEta1],
							 xmin[cPairsDims::kPhi1],
							 0,  // Dummy
							 xmin[cPairsDims::kCent],
							 xmin[cPairsDims::kZvtx]};
    Double_t xmax_singles[cSinglesDims::kNdimensions] = {xmax[cPairsDims::kEta1],
							 xmax[cPairsDims::kPhi1],
							 1,  // Dummy
							 xmax[cPairsDims::kCent],
							 xmax[cPairsDims::kZvtx]};
    this->fSingles = new THnS("singles", "<N_{1}>;#eta_{1};#phi_{1};p_{T};cent;z_{vtx};",
			      cSinglesDims::kNdimensions, nbins_singles, xmin_singles, xmax_singles);
    this->fSingles->GetAxis(cSinglesDims::kPt)->Set(npTbins, pt_bin_edges);
    this->fSingles->GetAxis(cSinglesDims::kCent)->Set(nCent, cent_bin_edges);
    this->fOutputList->Add(this->fSingles);
  }
  {
    // Event counter
    Int_t nbins_evt_counter[cEventCounterDims::kNdimensions] =
      {nbins[cEventCounterDims::kCent], nbins[cPairsDims::kZvtx]};
    Double_t xmin_evt_counter[cEventCounterDims::kNdimensions] =
      {xmin[cEventCounterDims::kCent], xmin[cPairsDims::kZvtx]};
    Double_t xmax_evt_counter[cEventCounterDims::kNdimensions] =
      {xmax[cEventCounterDims::kCent], xmax[cPairsDims::kZvtx]};
    this->fEventCounter = new THnF("eventCounter", "Event counter;cent;z_{vtx};",
				   cEventCounterDims::kNdimensions,
				   nbins_evt_counter, xmin_evt_counter, xmax_evt_counter);
    this->fEventCounter->GetAxis(cEventCounterDims::kCent)->Set(nCent, cent_bin_edges);
    this->fOutputList->Add(this->fEventCounter);
  }

  // Debug hists
  this->fcents = new TH1F("cents", "cents", 210, 0, 210);
  this->fOutputList->Add(this->fcents);

  // this->fRndmGenerator = new TRandom3();
  
  AliLog::SetGlobalLogLevel(AliLog::kError);
  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskC2::UserExec(Option_t *)
{
  // AODEvent() returns a null pointer for whatever reason...
  AliMCEvent*   mcEvent = this->MCEvent();
  AliAODEvent* aodEvent = dynamic_cast< AliAODEvent* >(this->InputEvent());
  if (!this->IsValidEvent())
    return;
  AliMultSelection *multSelection =
    dynamic_cast< AliMultSelection* >(this->InputEvent()->FindListObject("MultSelection"));
  if( !multSelection){
    AliWarning("AliMultSelection object not found!");
    return;
  }
  
  const Float_t centrality = multSelection->GetMultiplicityPercentile("V0M");
  this->fcents->Fill(centrality);

  const Double_t weight = (kMCTRUTH == fMode)
    ? 1 //mcEvent->GenEventHeader()->EventWeight()
    : 1;
  const Double_t zvtx = (this->InputEvent()->GetPrimaryVertex())
    ? this->InputEvent()->GetPrimaryVertex()->GetZ()
    : -999;
  
  // Yes, the following is realy "aodEvent" not mcEvent :P
  TClonesArray* tracksArray = (kMCTRUTH == fMode)
    ? dynamic_cast<TClonesArray*>(aodEvent->GetList()->FindObject(AliAODMCParticle::StdBranchName()))
    : aodEvent->GetTracks();
  if (!tracksArray || tracksArray->GetSize()==0){
    AliWarning("Could not retrieve track array");
    return;
  }

  std::vector< cNano_track > tracks;
  TIter nextTrack(tracksArray);
  while (TObject* obj = nextTrack()){
    // The naming around VTrack, VParticle, AODTrack mcParticle is a mess!
    // Take-away message: They all derive from AliVParticle one way or another.
    AliVParticle* particle = dynamic_cast< AliVParticle* >(obj);
    if (!this->IsValidParticle(particle))
      continue;
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
			    centrality,
			    zvtx};
    this->fSingles->Fill(stuffing, weight);
    // Save the track for the pair histogram
    tracks.push_back(tmp_track);
  }
  // If we have at least one valid track, it is a valid event
  // If there are no tracks, do not count the event since this would be an undefined correlation
  if (tracks.size() == 0){
    // This might be slightly inconsistent and misses the PostData!
    //this->fDiscardedEvents->Fill(cDiscardEventReasons::noTracksInPtRegion);
    //return;
  }
  {
    Double_t stuffing[3] = {centrality, zvtx};
    this->fEventCounter->Fill(stuffing, weight);
  }
  // Fill the single particle density
  for (std::vector<cNano_track>::size_type iTrack = 0; iTrack < tracks.size(); iTrack++) {

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
	   centrality,
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
