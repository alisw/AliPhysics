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
  AliAnalysisTaskEmcal("AliEmcalClusTrackMatcherTask"),
  fMaxDistance(0.1)
{
  // Constructor.
}

//________________________________________________________________________
AliEmcalClusTrackMatcherTask::AliEmcalClusTrackMatcherTask(const char *name) : 
  AliAnalysisTaskEmcal(name),
  fMaxDistance(0.1)
{
  // Standard constructor.
}

//________________________________________________________________________
AliEmcalClusTrackMatcherTask::~AliEmcalClusTrackMatcherTask()
{
  // Destructor.
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
}

