// $Id$
//
// Jet randomizer task.
//
// Author: S.Aiola, C.Loizides

#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TRandom3.h>

#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliVCluster.h"
#include "AliEMCALDigit.h"
#include "AliEMCALRecPoint.h"
#include "AliPicoTrack.h"
#include "AliEMCALGeometry.h"
#include "AliLog.h"

#include "AliJetRandomizerTask.h"

ClassImp(AliJetRandomizerTask)

//________________________________________________________________________
AliJetRandomizerTask::AliJetRandomizerTask() : 
  AliJetModelBaseTask("AliJetRandomizerTask")
{
  // Default constructor.

  SetSuffix("Randomized");
}

//________________________________________________________________________
AliJetRandomizerTask::AliJetRandomizerTask(const char *name) : 
  AliJetModelBaseTask(name)
{
  // Standard constructor.

  SetSuffix("Randomized");
}

//________________________________________________________________________
AliJetRandomizerTask::~AliJetRandomizerTask()
{
  // Destructor
}

//________________________________________________________________________
void AliJetRandomizerTask::Run() 
{
  // Randomize particles.

  if (fNClusters > 0 && fOutClusters) {
    const Int_t nClusters = fClusters->GetEntriesFast();
    for (Int_t i = 0; i < nClusters; ++i) {
      AliVCluster *cluster = dynamic_cast<AliVCluster*>(fClusters->At(i));
      if (!cluster)
	continue;
      if (!cluster->IsEMCAL())
	continue;

      Float_t pos[3];
      cluster->GetPosition(pos);
      TVector3 clusVec(pos);

      AddCluster(cluster->E(), clusVec.Eta());
    }
  }
 
  if (fNTracks > 0 && fOutTracks) {
    const Int_t nTracks = fTracks->GetEntriesFast();
    for (Int_t i = 0; i < nTracks; ++i) {
      AliPicoTrack *track = dynamic_cast<AliPicoTrack*>(fTracks->At(i));
      if (!track)
	continue;

      AddTrack(track->Pt(), track->Eta());
    }
  }
}

//________________________________________________________________________
void AliJetRandomizerTask::UserExec(Option_t *) 
{
  // Execute per event.

  if (!fCopyArray) {
    AliWarning("fCopyArray == kFALSE not allowed for AliJetRandomizerTask, will set kTRUE");
    fCopyArray = kTRUE;
  }

  Init();

  Run();
}
