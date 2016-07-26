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
    fSingleParticleForEfficiency(0)
{
}

//________________________________________________________________________
AliAnalysisTaskC2::AliAnalysisTaskC2(const char *name, Int_t mode)
  : AliAnalysisTaskC2Base(name, mode),
    fEventCounter(0),
    fSingleParticleForEfficiency(0)
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
  const Int_t ndim_pairs   = 6;
  const Int_t ndim_singles = 4;
  const Int_t etaNbins = 15;
  // phiNbins must be divisable by 2, but not by 4; so that we can later shift it by pi/2 (2pi is total int.)
  // The idea is to have the deltaPhi histogram with a bin centered arround 0
  const Int_t phiNbins = 26;
  const Double_t mult_bin_edges[] = {0, 1, 5, 50, 80, 101};
  const Int_t nMult    = sizeof(mult_bin_edges)/sizeof(mult_bin_edges[0]) - 1;
  const Int_t nZvtx    = 20;
  enum                            {kEta1,    kEta2,    kPhi1,    kPhi2,    kMult, kZvtx};
  const Int_t nbins[ndim_pairs] = {etaNbins, etaNbins, phiNbins, phiNbins, nMult, nZvtx};
  const Double_t xmin[ndim_pairs] = {-0.8, -0.8, 0, 0, 0, -10};
  const Double_t xmax[ndim_pairs] = {0.8, 0.8, 2*TMath::Pi(), 2*TMath::Pi(), 101, 10};
  
  this->fPairs = new THnS("pairs", "<N_{1}N_{2}>;#eta_{1};#eta_{2};#phi_{1};#phi_{2};mult;z_{vtx};",
			  ndim_pairs, nbins, xmin, xmax);
  this->fPairs->GetAxis(kMult)->Set(nMult, mult_bin_edges);
  // exclude_over_under_flow(this->fPairs);
  this->fOutputList->Add(fPairs);

  {
    const Int_t ndims_single1 = 4;
    Int_t nbins_single1[ndims_single1] = {nbins[kEta1], nbins[kPhi1], nbins[kMult], nbins[kZvtx]};
    Double_t xmin_single1[ndims_single1] = {xmin[kEta1], xmin[kPhi1], xmin[kMult], xmin[kZvtx]};
    Double_t xmax_single1[ndims_single1] = {xmax[kEta1], xmax[kPhi1], xmax[kMult], xmax[kZvtx]};
    this->fSingle1 = new THnF("single1", "<N_{1}>;#eta_{1};#phi_{1};mult;z_{vtx};",
			      ndims_single1, nbins_single1, xmin_single1, xmax_single1);
    this->fSingle1->GetAxis(2)->Set(nMult, mult_bin_edges);
    this->fOutputList->Add(this->fSingle1);
  }
  {
    const Int_t ndims_single2 = 4;
    Int_t nbins_single2[ndims_single2] = {nbins[kEta2], nbins[kPhi2], nbins[kMult], nbins[kZvtx]};
    Double_t xmin_single2[ndims_single2] = {xmin[kEta2], xmin[kPhi2], xmin[kMult], xmin[kZvtx]};
    Double_t xmax_single2[ndims_single2] = {xmax[kEta2], xmax[kPhi2], xmax[kMult], xmax[kZvtx]};
    this->fSingle2 = new THnF("single2", "<N_{2}>;#eta_{2};#phi_{2};mult;z_{vtx};",
				    ndims_single2, nbins_single2, xmin_single2, xmax_single2);
    this->fSingle2->GetAxis(2)->Set(nMult, mult_bin_edges);
    this->fOutputList->Add(this->fSingle2);
  }
  {
    // Event counter
    const Int_t ndims_evt_counter = 2;
    Int_t nbins_evt_counter[ndims_evt_counter] = {nbins[kMult], nbins[kZvtx]};
    Double_t xmin_evt_counter[ndims_evt_counter] = {xmin[kMult], xmin[kZvtx]};
    Double_t xmax_evt_counter[ndims_evt_counter] = {xmax[kMult], xmax[kZvtx]};
    this->fEventCounter = new THnF("eventConter", "Event counter;mult;z_{vtx};",
					 ndims_evt_counter, nbins_evt_counter, xmin_evt_counter, xmax_evt_counter);
    this->fEventCounter->GetAxis(0)->Set(nMult, mult_bin_edges);
    this->fOutputList->Add(this->fEventCounter);
  }

  // The strange binning is deliberate! Think about it!  Thought about
  // it! The delta hist can be filled with values of [-2, 2] if the
  // data is straight from the tree. But only [-2 + binwidth/2, 2 -
  // bindwidth] if the data is prebinned!
  const Float_t etaBinWidth = this->fPairs->GetAxis(0)->GetBinWidth(1);
  const Float_t phiBinWidth = this->fPairs->GetAxis(2)->GetBinWidth(1);

  this->fPairsConventional = AliAnalysisC2Utils::CreateTransformedHist(this->fPairs);
  this->fOutputList->Add(fPairsConventional);

  // Debug hists
  this->fdNdeta = new TH1F("dNdeta", "dNdeta", etaNbins, -1, 1);
  this->fOutputList->Add(this->fdNdeta);

  this->fmultDistribution = new TH1F("multDistribution", "multDistribution", 210, 0, 210);
  this->fOutputList->Add(this->fmultDistribution);

  // Histogram to construct efficiency in eta, phi, pt, mult, z_vtx
  const Int_t ndims_eff = 5;
  const Double_t pt_bin_edges[] = {0., 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
				   1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6,
				   3.0, 3.4, 3.8, 4.2, 4.6, 5.0, 5.5, 6.0, 7.0,
				   8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 100.0};
  const Int_t npTbins = sizeof(pt_bin_edges)/sizeof(pt_bin_edges[0]) - 1;
  Int_t nbins_eff[ndims_eff] = {etaNbins, phiNbins, npTbins, nbins[kMult], nbins[kZvtx]};
  // Dummy values for pT; the axis is later changed to variable bin sizes!
  Double_t xmin_eff[ndims_eff] = {xmin[kEta2], xmin[kPhi2], 0.0, xmin[kMult], xmin[kZvtx]};
  Double_t xmax_eff[ndims_eff] = {xmax[kEta2], xmax[kPhi2], 100.0, xmax[kMult], xmax[kZvtx]};
  this->fSingleParticleForEfficiency = new THnF("singleParticleForEfficiency",
						"Single particle density for efficiency evaluation;"
						"#eta;#phi;p_{T};mult;z_{vtx};",
						ndims_eff, nbins_eff, xmin_eff, xmax_eff);
  this->fSingleParticleForEfficiency->GetAxis(2)->Set(npTbins, pt_bin_edges);
  this->fSingleParticleForEfficiency->GetAxis(3)->Set(nMult, mult_bin_edges);
  this->fOutputList->Add(this->fSingleParticleForEfficiency);

  this->fRndmGenerator = new TRandom3();
  
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
  
  const Float_t multiplicity = multSelection->GetMultiplicityPercentile("V0M");
  this->fmultDistribution->Fill(multiplicity);

  const Double_t weight = (kMCTRUTH == fMode)
    ? 1 //mcEvent->GenEventHeader()->EventWeight()
    : 1;
  const Double_t zvtx = (this->InputEvent()->GetPrimaryVertex())
    ? this->InputEvent()->GetPrimaryVertex()->GetZ()
    : 0;
  
  std::vector< cNano_track > tracks1, tracks2;
  // Yes, the following is realy "aodEvent" not mcEvent :P
  TClonesArray* tracks = (kMCTRUTH == fMode)
    ? dynamic_cast<TClonesArray*>(aodEvent->GetList()->FindObject(AliAODMCParticle::StdBranchName()))
    : aodEvent->GetTracks();
  if (!tracks || tracks->GetSize()==0){
    AliWarning("Could not retrieve track array");
    return;
  }

  TIter nextTrack(tracks);
  while (TObject* obj = nextTrack()){
    // The naming around VTrack, VParticle, AODTrack mcParticle is a mess!
    // Take-away message: They all derive from AliVParticle one way or another.
    AliVParticle* particle = dynamic_cast< AliVParticle* >(obj);
    if (!this->IsValidParticle(particle))
      continue;
    // eta, phi, pt, q are pure virtual in AliVParticle, so it should be safe to just call them on AliVParticle
    cNano_track tmp_track = {
      .eta = particle->Eta(), //fRndmGenerator->Rndm() * 1.6 - 0.8;
      .phi = particle->Phi(), //fRndmGenerator->Rndm() * 2*TMath::Pi();
      .pt = particle->Pt()};
    this->fdNdeta->Fill(tmp_track.eta, weight);
    {
      Double_t stuffing[5] = {tmp_track.eta, tmp_track.phi, tmp_track.pt, multiplicity, zvtx};
      this->fSingleParticleForEfficiency->Fill(stuffing);
    }

    if ((tmp_track.pt > 1) && (tmp_track.pt < 2))
      tracks1.push_back(tmp_track);
    if ((tmp_track.pt > 2) && (tmp_track.pt < 4))
      tracks2.push_back(tmp_track);
  }
  // If we do not have tracks we have no pairs
  // in this scenario, the value of a "correlation" is undefined!?
  // Only fill histogram if we have at least one pair (ie. one trigger and one assocs)
  if (tracks1.size() == 0 || tracks2.size() == 0)
    return;
  // increment event counter and keep the stuffing in scope
  {
    Double_t stuffing[3] = {multiplicity, zvtx};
    this->fEventCounter->Fill(stuffing, weight);
  }

  for (std::vector<cNano_track>::size_type iType1 = 0; iType1 < tracks1.size(); iType1++) {
    Double_t stuffing[4] = {tracks1[iType1].eta,
			    AliAnalysisC2Utils::WrapAngle(tracks1[iType1].phi, this->fSingle1->GetAxis(1)),
			    multiplicity,
			    zvtx};
    this->fSingle1->Fill(stuffing, weight);
  }
  for (std::vector<cNano_track>::size_type iType2 = 0; iType2 < tracks2.size(); iType2++) {
    Double_t stuffing[4] = {tracks2[iType2].eta,
			    AliAnalysisC2Utils::WrapAngle(tracks2[iType2].phi, this->fSingle2->GetAxis(1)),
			    multiplicity,
			    zvtx};
    this->fSingle2->Fill(stuffing, weight);
  }
  // n: number of valid tracks
  // What is the number of valid pairs? Untriggered: n^2-n? This does contain mirrored pairs!
  for (std::vector<cNano_track>::size_type iType1 = 0; iType1 < tracks1.size(); iType1++) {
    for (std::vector<cNano_track>::size_type iType2 = 0; iType2 < tracks2.size(); iType2++) {
      {
	Double_t stuffing[6] =
	  {tracks1[iType1].eta,
	   tracks2[iType2].eta,
	   AliAnalysisC2Utils::WrapAngle(tracks1[iType1].phi, this->fPairs->GetAxis(2)),
	   AliAnalysisC2Utils::WrapAngle(tracks2[iType2].phi, this->fPairs->GetAxis(3)),
	   multiplicity,
	   zvtx};
	this->fPairs->Fill(stuffing, weight * weight);
      	Double_t stuffing_trans[6];
	AliAnalysisC2Utils::TransformPoints(stuffing, stuffing_trans);
	stuffing_trans[2] = AliAnalysisC2Utils::WrapAngle(stuffing_trans[2], this->fPairsConventional->GetAxis(2));
	this->fPairsConventional->Fill(stuffing_trans, weight * weight);
      }
    }
  }
  PostData(1, this->fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskC2::Terminate(Option_t *)
{
  // PostData(1, this->fOutputList);
}
