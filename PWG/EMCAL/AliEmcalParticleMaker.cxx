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

    if (!(InputEvent()->FindListObject(fTracksOutName))) {
      InputEvent()->AddObject(fTracksOut);
    }
    else {
      fInitialized = kFALSE;
      AliFatal(Form("%s: Container with same name %s already present. Aborting", GetName(), fTracksOutName.Data()));
      return;
    }
  }

  if (!fCaloOutName.IsNull()) {
    fCaloClustersOut = new TClonesArray("AliEmcalParticle");
    fCaloClustersOut->SetName(fCaloOutName);
    
    // post output in event if not yet present
    if (!(InputEvent()->FindListObject(fCaloOutName))) {
      InputEvent()->AddObject(fCaloClustersOut);
    }
    else {
      fInitialized = kFALSE;
      AliFatal(Form("%s: Container with same name %s already present. Aborting", GetName(), fCaloOutName.Data()));
      return;
    }
  }
}

//________________________________________________________________________
Bool_t AliEmcalParticleMaker::Run() 
{
  // Create the emcal particles

  if (fTracks && fTracksOut) {
    // clear container (normally a null operation as the event should clean it already)
    fTracksOut->Delete();

    // loop over tracks
    const Int_t Ntracks = fTracks->GetEntries();
    for (Int_t iTracks = 0; iTracks < Ntracks; ++iTracks) {
      
      AliVTrack *track = dynamic_cast<AliVTrack*>(fTracks->At(iTracks));
      new ((*fTracksOut)[iTracks]) AliEmcalParticle(track, iTracks);
    }
  }

  if (fCaloClusters && fCaloClustersOut) {
    // clear container (normally a null operation as the event should clean it already)
    fCaloClustersOut->Delete();

    // loop over clusters
    const Int_t Nclusters = fCaloClusters->GetEntries();
    for (Int_t iClusters = 0; iClusters < Nclusters; ++iClusters) {
      AliVCluster *cluster = dynamic_cast<AliVCluster*>(fCaloClusters->At(iClusters));
      new ((*fCaloClustersOut)[iClusters]) AliEmcalParticle(cluster, iClusters, fVertex[0], fVertex[1], fVertex[2]);
    }
  }
  return kTRUE;
}
