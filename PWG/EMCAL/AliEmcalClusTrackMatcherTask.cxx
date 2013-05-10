// $Id$
//
// Track/cluster matcher
// 
// Author: C.Loizides, S.Aiola

#include "AliEmcalClusTrackMatcherTask.h"

#include <TClonesArray.h>

#include "AliEmcalParticle.h"
#include "AliLog.h"
#include "AliPicoTrack.h"
#include "AliVCluster.h"
#include "AliVTrack.h"

ClassImp(AliEmcalClusTrackMatcherTask)

//________________________________________________________________________
AliEmcalClusTrackMatcherTask::AliEmcalClusTrackMatcherTask() : 
  AliAnalysisTaskEmcal("AliEmcalClusTrackMatcherTask",kFALSE),
  fMaxDistance(0.06)
{
  // Constructor.

  for(Int_t icent=0; icent<8; ++icent) {
    for(Int_t ipt=0; ipt<9; ++ipt) {
      for(Int_t ieta=0; ieta<2; ++ieta) {
	fHistMatchEta[icent][ipt][ieta] = 0;
	fHistMatchPhi[icent][ipt][ieta] = 0;
      }
    }
  }
}

//________________________________________________________________________
AliEmcalClusTrackMatcherTask::AliEmcalClusTrackMatcherTask(const char *name, Bool_t histo) : 
  AliAnalysisTaskEmcal(name,histo),
  fMaxDistance(0.06)
{
  // Standard constructor.

  for(Int_t icent=0; icent<8; ++icent) {
    for(Int_t ipt=0; ipt<9; ++ipt) {
      for(Int_t ieta=0; ieta<2; ++ieta) {
	fHistMatchEta[icent][ipt][ieta] = 0;
	fHistMatchPhi[icent][ipt][ieta] = 0;
      }
    }
  }
}

//________________________________________________________________________
AliEmcalClusTrackMatcherTask::~AliEmcalClusTrackMatcherTask()
{
  // Destructor.
}

//________________________________________________________________________
void AliEmcalClusTrackMatcherTask::UserCreateOutputObjects()
{
  // Create my user objects.

  if (!fCreateHisto)
    return;

  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  const Int_t nCentChBins = fNcentBins * 2;

  for(Int_t icent=0; icent<nCentChBins; ++icent) {
    for(Int_t ipt=0; ipt<9; ++ipt) {
      for(Int_t ieta=0; ieta<2; ++ieta) {
	TString nameEta(Form("fHistMatchEta_%i_%i_%i",icent,ipt,ieta));
	fHistMatchEta[icent][ipt][ieta] = new TH1F(nameEta, nameEta, 400, -0.2, 0.2);
	TString namePhi(Form("fHistMatchPhi_%i_%i_%i",icent,ipt,ieta));
	fHistMatchPhi[icent][ipt][ieta] = new TH1F(namePhi, namePhi, 400, -0.2, 0.2);
	fOutput->Add(fHistMatchEta[icent][ipt][ieta]);
	fOutput->Add(fHistMatchPhi[icent][ipt][ieta]);
      }
    }
  }

  PostData(1, fOutput);
}

//________________________________________________________________________
Int_t AliEmcalClusTrackMatcherTask::GetMomBin(Double_t p) const
{
  // Get momenum bin.

  Int_t pbin=-1;
  if (p<0.5) 
    pbin=0;
  else if (p>=0.5 && p<1.0) 
    pbin=1;
  else if (p>=1.0 && p<1.5) 
    pbin=2;
  else if (p>=1.5 && p<2.) 
    pbin=3;
  else if (p>=2. && p<3.) 
    pbin=4;
  else if (p>=3. && p<4.) 
    pbin=5;
  else if (p>=4. && p<5.) 
    pbin=6;
  else if (p>=5. && p<8.) 
    pbin=7;
  else if (p>=8.) 
    pbin=8;

  return pbin;
}

//________________________________________________________________________
Bool_t AliEmcalClusTrackMatcherTask::Run() 
{
  // Run the matching for the selected options.
  
  const Double_t maxd2 = fMaxDistance*fMaxDistance;

  const Int_t nC = fCaloClusters->GetEntries();
  const Int_t nT = fTracks->GetEntries();

  // set the links between tracks and clusters
  for (Int_t c = 0; c < nC; ++c) {
    AliEmcalParticle *partC = static_cast<AliEmcalParticle*>(fCaloClusters->At(c));
    if (!partC)
      continue;
    if (!AcceptEmcalPart(partC))
      continue;
    for (Int_t t = 0; t < nT; ++t) {
      AliEmcalParticle *partT = static_cast<AliEmcalParticle*>(fTracks->At(t));
      if (!partT)
	continue;
      if (!AcceptEmcalPart(partT))
        continue;
      AliVCluster *clust = partC->GetCluster();
      AliVTrack   *track = partT->GetTrack()  ;
      Double_t deta = 999;
      Double_t dphi = 999;
      AliPicoTrack::GetEtaPhiDiff(track, clust, dphi, deta);
      Double_t d2 = deta * deta + dphi * dphi;
      if (d2 > maxd2)
        continue;
      Double_t d = TMath::Sqrt(d2);
      partC->AddMatchedObj(t, d);
      partT->AddMatchedObj(c, d);

      if (fCreateHisto) {
	Int_t mombin = GetMomBin(track->P());
	Int_t centbinch = fCentBin;
	if (track->Charge()<0) 
	  centbinch += fNcentBins;
	Int_t etabin = 0;
	if(track->Eta() > 0) 
	  etabin = 1;
	    
	fHistMatchEta[centbinch][mombin][etabin]->Fill(deta);
	fHistMatchPhi[centbinch][mombin][etabin]->Fill(dphi);
      }
    }
  }

  for (Int_t c = 0; c < nC; ++c) {
    AliEmcalParticle *partC = static_cast<AliEmcalParticle*>(fCaloClusters->At(c));
    if (!partC)
      continue;
    if (partC->GetNumberOfMatchedObj() <= 0)
      continue;
    const UInt_t matchedId = partC->GetMatchedObjId();
    AliEmcalParticle *partT = static_cast<AliEmcalParticle*>(fTracks->At(matchedId));
    AliVCluster *clust = partC->GetCluster();
    AliVTrack   *track = partT->GetTrack()  ;
    Double_t deta = 999;
    Double_t dphi = 999;
    AliPicoTrack::GetEtaPhiDiff(track, clust, dphi, deta);
    clust->SetEmcCpvDistance(matchedId);
    clust->SetTrackDistance(deta, dphi);
  }

  for (Int_t t = 0; t < nT; ++t) {
    AliEmcalParticle *partT = static_cast<AliEmcalParticle*>(fTracks->At(t));
    if (!partT)
      continue;
    if (partT->GetNumberOfMatchedObj() <= 0)
      continue;
    AliVTrack *track = partT->GetTrack();
    track->SetEMCALcluster(partT->GetMatchedObjId());
  }

  return kTRUE;
}
