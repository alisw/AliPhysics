// $Id$
//
// Track/cluster matcher
// 
// Author: C.Loizides, S.Aiola

#include <TClonesArray.h>

#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliPicoTrack.h"
#include "AliEmcalParticle.h"
#include "AliLog.h"

#include "AliEmcalClusTrackMatcherTask.h"

ClassImp(AliEmcalClusTrackMatcherTask)

//________________________________________________________________________
AliEmcalClusTrackMatcherTask::AliEmcalClusTrackMatcherTask() : 
  AliAnalysisTaskEmcal("AliEmcalClusTrackMatcherTask"),
  fDoClusTrack(1),
  fDoTrackClus(0),
  fMaxDistance(0.1)
{
  // Constructor.
}

//________________________________________________________________________
AliEmcalClusTrackMatcherTask::AliEmcalClusTrackMatcherTask(const char *name) : 
  AliAnalysisTaskEmcal(name),
  fDoClusTrack(1),
  fDoTrackClus(0),
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
  
  if (fDoClusTrack) 
    DoMatching(fCaloClusters, fTracks);

  if (fDoTrackClus) 
    DoMatching(fTracks, fCaloClusters);
  
  return kTRUE;
}

//________________________________________________________________________
void AliEmcalClusTrackMatcherTask::DoMatching(TClonesArray *array1, TClonesArray *array2) 
{
  // Do the actual matching.

  if (!array1 || !array2)
    return;

  const Int_t n1 = array1->GetEntries();
  const Int_t n2 = array2->GetEntries();
  for (Int_t i = 0; i < n1; ++i) {
    AliEmcalParticle *part1 = static_cast<AliEmcalParticle*>(array1->At(i));
    if (!part1)
      continue;
      
    AliVCluster *cluster1 = part1->GetCluster();
    AliVTrack   *track1   = part1->GetTrack()  ;

    if ((!cluster1 || !AcceptCluster(cluster1)) && 
        (!track1 || !AcceptTrack(track1)))
      continue;
     
    part1->ResetMatchedObjects();

    for (Int_t j = 0; j < n2; ++j) {

      AliEmcalParticle *part2 = static_cast<AliEmcalParticle*>(array2->At(j));
      if (!part2)
	continue;

      AliVCluster *cluster2 = part2->GetCluster();
      AliVTrack   *track2   = part2->GetTrack()  ;

      if ((!cluster2 || !AcceptCluster(cluster2)) && 
          (!track2 || !AcceptTrack(track2)))
	continue;

      Double_t deta = 999;
      Double_t dphi = 999;
      if (track1 && cluster2)
	AliPicoTrack::GetEtaPhiDiff(track1, cluster2, dphi, deta);
      else if (track2 && cluster1)
	AliPicoTrack::GetEtaPhiDiff(track2, cluster1, dphi, deta);
      else
	continue;
	   
      Double_t d = TMath::Sqrt(deta * deta + dphi * dphi);
      if(d < fMaxDistance) 
	part1->AddMatchedObj(j, d);
    }
    
    if (part1->GetNumberOfMatchedObj() > 0) {
      if (track1) {
	track1->SetEMCALcluster(part1->GetMatchedObjId());
      } else if (cluster1) {
	const UInt_t matchedId = part1->GetMatchedObjId();
	Double_t deta = 999;
	Double_t dphi = 999;
	AliEmcalParticle *part2 = static_cast<AliEmcalParticle*>(array2->At(matchedId));
	AliVTrack *track2 = 0;
        if (part2)
          track2 = part2->GetTrack();
	if (track2) {
	  AliPicoTrack::GetEtaPhiDiff(track2, cluster1, dphi, deta);
	  cluster1->SetEmcCpvDistance(matchedId);
	  cluster1->SetTrackDistance(deta, dphi);
	}
      }
    }
  }
}
