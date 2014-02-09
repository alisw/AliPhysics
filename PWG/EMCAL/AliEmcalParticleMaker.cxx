// $Id$
//
// Class to make emcal particles in AOD/ESD events.
//
// Author: S.Aiola

#include <TClonesArray.h>

#include "AliLog.h"
#include "AliVTrack.h"
#include "AliVCluster.h"
#include "AliEmcalParticle.h"

#include "AliEmcalParticleMaker.h"

ClassImp(AliEmcalParticleMaker)

//________________________________________________________________________
AliEmcalParticleMaker::AliEmcalParticleMaker() : 
  AliAnalysisTaskEmcal("AliEmcalParticleMaker",kFALSE),
  fTracksOutName("EmcalTracks"),
  fCaloOutName("EmcalClusters"),
  fTracksOut(0),
  fCaloClustersOut(0)
{
  // Constructor.
}

//________________________________________________________________________
AliEmcalParticleMaker::AliEmcalParticleMaker(const char *name) : 
  AliAnalysisTaskEmcal(name,kFALSE),
  fTracksOutName("EmcalTracks"),
  fCaloOutName("EmcalClusters"),
  fTracksOut(0),
  fCaloClustersOut(0)
{
  // Constructor.
}

//________________________________________________________________________
AliEmcalParticleMaker::~AliEmcalParticleMaker()
{
  // Destructor.
}

//________________________________________________________________________
void AliEmcalParticleMaker::ExecOnce()
{
  // Init the analysis.

  AliAnalysisTaskEmcal::ExecOnce();

  if (!fInitialized)
    return;

  if (!fTracksOutName.IsNull()) {
    fTracksOut = new TClonesArray("AliEmcalParticle");
    fTracksOut->SetName(fTracksOutName);
    AddObjectToEvent(fTracksOut);
  }

  if (!fCaloOutName.IsNull()) {
    fCaloClustersOut = new TClonesArray("AliEmcalParticle");
    fCaloClustersOut->SetName(fCaloOutName);
    AddObjectToEvent(fCaloClustersOut);    
  }
}

//________________________________________________________________________
Bool_t AliEmcalParticleMaker::Run() 
{
  // Create the emcal particles

  if (fTracks && fTracksOut) {
    // clear container (normally a null operation as the event should clean it already)
    fTracksOut->Delete();
    const Int_t Ntracks = fTracks->GetEntries();
    for (Int_t iTracks = 0; iTracks < Ntracks; ++iTracks) {
      AliVTrack *track = static_cast<AliVTrack*>(fTracks->At(iTracks));
      AliEmcalParticle *ep = new ((*fTracksOut)[iTracks]) AliEmcalParticle(track, iTracks);
      if (0&&fCaloClusters)
	ep->SetMatchedPtr(fCaloClusters);
    }
  }

  if (fCaloClusters && fCaloClustersOut) {
    // clear container (normally a null operation as the event should clean it already)
    fCaloClustersOut->Delete();
    const Int_t Nclusters = fCaloClusters->GetEntries();
    for (Int_t iClusters = 0, iN=0; iClusters < Nclusters; ++iClusters) {
      AliVCluster *cluster = static_cast<AliVCluster*>(fCaloClusters->At(iClusters));
      /* Commented because for simplicity prefer to keep indices aligned with clusters (CL)
        if (!cluster->IsEMCAL()) continue;
      */
      AliEmcalParticle *ep = new ((*fCaloClustersOut)[iN++]) AliEmcalParticle(cluster, iClusters, fVertex[0], fVertex[1], fVertex[2]);
      if (0&&fTracks)
	ep->SetMatchedPtr(fTracks);
    }
  }
  return kTRUE;
}
