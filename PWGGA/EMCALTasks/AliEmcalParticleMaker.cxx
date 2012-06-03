// $Id: AliEmcalParticleMaker.cxx 56635 2012-05-22 22:49:56Z loizides $
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
  AliAnalysisTaskEmcal("AliEmcalParticleMaker"),
  fTracksOutName("EmcalTracks"),
  fCaloOutName("EmcalClusters"),
  fTracksOut(0),
  fCaloClustersOut(0)
{
  // Constructor.
}

//________________________________________________________________________
AliEmcalParticleMaker::AliEmcalParticleMaker(const char *name) : 
  AliAnalysisTaskEmcal(name),
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
void AliEmcalParticleMaker::UserCreateOutputObjects()
{
  // Create my user objects.

  fTracksOut = new TClonesArray("AliEmcalParticle");
  fTracksOut->SetName(fTracksOutName);

  fCaloClustersOut = new TClonesArray("AliEmcalParticle");
  fCaloClustersOut->SetName(fCaloOutName);
}

//________________________________________________________________________
Bool_t AliEmcalParticleMaker::Run() 
{
  // Create th emcal particles

  // add tracks to event if not yet there
  if (!(InputEvent()->FindListObject(fTracksOutName))) {
    InputEvent()->AddObject(fTracksOut);
  }
  if (!(InputEvent()->FindListObject(fCaloOutName))) {
    InputEvent()->AddObject(fCaloClustersOut);
  }

  // clear container (normally a null operation as the event should clean it already)
  fTracksOut->Delete();
  fCaloClustersOut->Delete();

  // loop over tracks
  const Int_t Ntracks = fTracks->GetEntries();
  for (Int_t iTracks = 0; iTracks < Ntracks; ++iTracks) {

    AliVTrack *track = dynamic_cast<AliVTrack*>(fTracks->At(iTracks));
    new ((*fTracksOut)[iTracks]) AliEmcalParticle(track, iTracks);
  }

  // loop over clusters
  const Int_t Nclusters = fCaloClusters->GetEntries();
  for (Int_t iClusters = 0; iClusters < Nclusters; ++iClusters) {

    AliVCluster *cluster = dynamic_cast<AliVCluster*>(fCaloClusters->At(iClusters));
    new ((*fCaloClustersOut)[iClusters]) AliEmcalParticle(cluster, iClusters);
  }

  return kTRUE;
}
