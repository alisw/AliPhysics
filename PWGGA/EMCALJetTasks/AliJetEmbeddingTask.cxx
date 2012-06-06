// $Id$
//
// Jet embedding task.
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

#include "AliJetEmbeddingTask.h"

ClassImp(AliJetEmbeddingTask)

//________________________________________________________________________
AliJetEmbeddingTask::AliJetEmbeddingTask() : 
  AliJetModelBaseTask("AliJetEmbeddingTask")
{
  // Default constructor.
  SetSuffix("Embedded");
}

//________________________________________________________________________
AliJetEmbeddingTask::AliJetEmbeddingTask(const char *name) : 
  AliJetModelBaseTask(name)
{
  // Standard constructor.
  SetSuffix("Embedded");
}

//________________________________________________________________________
AliJetEmbeddingTask::~AliJetEmbeddingTask()
{
  // Destructor
}


//________________________________________________________________________
void AliJetEmbeddingTask::Run() 
{
  // Embed particles.
  
  if (fNClusters > 0 && fOutClusters) {
    for (Int_t i = 0; i < fNClusters; ++i) {
      AddCluster();
    }
  }
 
  if (fNTracks > 0 && fOutTracks) {
    for (Int_t i = 0; i < fNTracks; ++i) {
      AddTrack();
    }
  }
}

//________________________________________________________________________
void AliJetEmbeddingTask::UserExec(Option_t *) 
{
  // Execute per event.

  Init();

  Run();
}
