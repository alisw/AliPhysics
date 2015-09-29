// $Id$
//
// Jet embedding task.
//
// Author: S.Aiola, C.Loizides

#include <TRandom3.h>

#include "AliJetEmbeddingTask.h"

ClassImp(AliJetEmbeddingTask)

//________________________________________________________________________
AliJetEmbeddingTask::AliJetEmbeddingTask() : 
  AliJetModelBaseTask("AliJetEmbeddingTask"),
  fMassless(kFALSE),
  fNeutralFraction(0),
  fNeutralMass(0.135),
  fMass(0.1396)
{
  // Default constructor.
  SetSuffix("Embedded");
}

//________________________________________________________________________
AliJetEmbeddingTask::AliJetEmbeddingTask(const char *name) : 
  AliJetModelBaseTask(name),
  fMassless(kFALSE),
  fNeutralFraction(0),
  fNeutralMass(0.135),
  fMass(0.1396)
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
    if (fCopyArray) 
      CopyClusters();
    for (Int_t i = 0; i < fNClusters; ++i) {
      AddCluster();
    }
  }
 
  if (fNTracks > 0 && fOutTracks) {
    if (fCopyArray) 
      CopyTracks();
    for (Int_t i = 0; i < fNTracks; ++i) {
      Double_t mass = fMass;
      Short_t charge = 1;
      if(fNeutralFraction>0.) {
	Double_t rnd = gRandom->Rndm();
	if(rnd<fNeutralFraction) {
	  charge = 0;
	  mass = fNeutralMass;
	}
      }
      if(fMassless) mass = 0.;
      AddTrack(-1,-999,-1,0,0,0,0,kFALSE,0,charge,mass);
    }
  }
}
