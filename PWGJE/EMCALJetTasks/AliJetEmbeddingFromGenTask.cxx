// $Id: AliJetEmbeddingFromGenTask.cxx 57324 2012-06-21 04:33:52Z loizides $
//
// Jet embedding task.
//
// Author: S.Aiola, C.Loizides

#include "AliJetEmbeddingFromGenTask.h"

#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TRandom3.h>

#include "AliAnalysisManager.h"
#include "AliEMCALDigit.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecPoint.h"
#include "AliGenerator.h"
#include "AliLog.h"
#include "AliPicoTrack.h"
#include "AliStack.h"
#include "AliVCluster.h"
#include "AliVEvent.h"

ClassImp(AliJetEmbeddingFromGenTask)

//________________________________________________________________________
AliJetEmbeddingFromGenTask::AliJetEmbeddingFromGenTask() : 
  AliJetModelBaseTask("AliJetEmbeddingFromGenTask"),
  fGen(0)
{
  // Default constructor.
  SetSuffix("EmbeddedFromGen");
}

//________________________________________________________________________
AliJetEmbeddingFromGenTask::AliJetEmbeddingFromGenTask(const char *name) : 
  AliJetModelBaseTask(name),
  fGen(0)
{
  // Standard constructor.
  SetSuffix("EmbeddedFromGen");
}

//________________________________________________________________________
AliJetEmbeddingFromGenTask::~AliJetEmbeddingFromGenTask()
{
  // Destructor
}

//________________________________________________________________________
void AliJetEmbeddingFromGenTask::UserExec(Option_t *) 
{
  // Execute per event.

  if (!fIsInit) {
    ExecOnce();
    fIsInit = 1;
  }
  Run();
}

//________________________________________________________________________
void AliJetEmbeddingFromGenTask::Run() 
{
  // Embed particles.

  fGen->Generate();

  AliStack *stack = fGen->GetStack();
  const Int_t nprim = stack->GetNprimary();
  for (Int_t i=0;i<nprim;++i) {
    TParticle *part = stack->Particle(i);
  }
}
