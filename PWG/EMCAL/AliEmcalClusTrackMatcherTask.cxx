// $Id$
//
// Track/cluster matcher
// 
// Author: C.Loizides, S.Aiola

#include "AliEmcalClusTrackMatcherTask.h"

#include <TClonesArray.h>

#include "AliParticleContainer.h"
#include "AliEmcalParticle.h"
#include "AliLog.h"
#include "AliPicoTrack.h"
#include "AliVCluster.h"
#include "AliVTrack.h"

ClassImp(AliEmcalClusTrackMatcherTask)

//________________________________________________________________________
AliEmcalClusTrackMatcherTask::AliEmcalClusTrackMatcherTask() : 
  AliAnalysisTaskEmcal("AliEmcalClusTrackMatcherTask",kFALSE),
  fMaxDistance(0.06),
  fHistMatchEtaAll(0),
  fHistMatchPhiAll(0)
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
  fMaxDistance(0.06),
  fHistMatchEtaAll(0),
  fHistMatchPhiAll(0)
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
void AliEmcalClusTrackMatcherTask::ExecOnce()
{
  // Initialize the analysis.

  if (fParticleCollArray.GetEntriesFast()<2) {
    AliError(Form("Wrong number of particle collections (%d), required 2",fParticleCollArray.GetEntriesFast()));
    return;
  }

  for (Int_t i = 0; i < 2; i++) {
    AliParticleContainer *cont = static_cast<AliParticleContainer*>(fParticleCollArray.At(i));
    cont->SetClassName("AliEmcalParticle");
  }

  AliAnalysisTaskEmcal::ExecOnce();
}

//________________________________________________________________________
void AliEmcalClusTrackMatcherTask::UserCreateOutputObjects()
{
  // Create my user objects.

  if (!fCreateHisto)
    return;

  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  const Int_t nCentChBins = fNcentBins * 2;

  fHistMatchEtaAll = new TH1F("fHistMatchEtaAll", "fHistMatchEtaAll", 400, -0.2, 0.2);
  fHistMatchPhiAll = new TH1F("fHistMatchPhiAll", "fHistMatchPhiAll", 400, -0.2, 0.2);
  fOutput->Add(fHistMatchEtaAll);
  fOutput->Add(fHistMatchPhiAll);

  for(Int_t icent=0; icent<nCentChBins; ++icent) {
    for(Int_t ipt=0; ipt<9; ++ipt) {
      for(Int_t ieta=0; ieta<2; ++ieta) {
	TString nameEta(Form("fHistMatchEta_%i_%i_%i",icent,ipt,ieta));
	fHistMatchEta[icent][ipt][ieta] = new TH1F(nameEta, nameEta, 400, -0.2, 0.2);
	fHistMatchEta[icent][ipt][ieta]->SetXTitle("#Delta#eta");
	TString namePhi(Form("fHistMatchPhi_%i_%i_%i",icent,ipt,ieta));
	fHistMatchPhi[icent][ipt][ieta] = new TH1F(namePhi, namePhi, 400, -0.2, 0.2);
	fHistMatchPhi[icent][ipt][ieta]->SetXTitle("#Delta#phi");
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

  AliParticleContainer *tracks = static_cast<AliParticleContainer*>(fParticleCollArray.At(0));
  AliParticleContainer *clusters = static_cast<AliParticleContainer*>(fParticleCollArray.At(1));

  AliEmcalParticle *partC = 0;
  AliEmcalParticle *partT = 0;

  const Double_t maxd2 = fMaxDistance*fMaxDistance;

  // set the links between tracks and clusters
  clusters->ResetCurrentID();
  while ((partC = static_cast<AliEmcalParticle*>(clusters->GetNextAcceptParticle()))) {
    AliVCluster *clust = partC->GetCluster();

    tracks->ResetCurrentID();
    while ((partT = static_cast<AliEmcalParticle*>(tracks->GetNextAcceptParticle()))) {
      AliVTrack *track = partT->GetTrack();
      Double_t deta = 999;
      Double_t dphi = 999;
      AliPicoTrack::GetEtaPhiDiff(track, clust, dphi, deta);
      Double_t d2 = deta * deta + dphi * dphi;
      if (d2 > maxd2)
        continue;

      Double_t d = TMath::Sqrt(d2);
      partC->AddMatchedObj(tracks->GetCurrentID(), d);
      partT->AddMatchedObj(clusters->GetCurrentID(), d);

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
        fHistMatchEtaAll->Fill(deta);
        fHistMatchPhiAll->Fill(dphi);
      }
    }
  }

  clusters->ResetCurrentID();
  while ((partC = static_cast<AliEmcalParticle*>(clusters->GetNextAcceptParticle()))) {
    if (partC->GetNumberOfMatchedObj() <= 0)
      continue;
    const UInt_t matchedId = partC->GetMatchedObjId();
    partT = static_cast<AliEmcalParticle*>(tracks->GetParticle(matchedId));
    AliVCluster *clust = partC->GetCluster();
    AliVTrack   *track = partT->GetTrack()  ;
    Double_t deta = 999;
    Double_t dphi = 999;
    AliPicoTrack::GetEtaPhiDiff(track, clust, dphi, deta);
    clust->SetEmcCpvDistance(matchedId);
    clust->SetTrackDistance(deta, dphi);
  }
  
  tracks->ResetCurrentID();
  while ((partT = static_cast<AliEmcalParticle*>(tracks->GetNextAcceptParticle()))) {
    if (partT->GetNumberOfMatchedObj() <= 0)
      continue;
    AliVTrack *track = partT->GetTrack();
    track->SetEMCALcluster(partT->GetMatchedObjId());
  }

  return kTRUE;
}
