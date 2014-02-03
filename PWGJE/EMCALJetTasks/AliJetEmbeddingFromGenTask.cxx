// $Id$
//
// Jet embedding task.
//
// Author: S.Aiola, C.Loizides

#include "AliJetEmbeddingFromGenTask.h"

#include <TClonesArray.h>
#include <TFolder.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TRandom3.h>
#include "AliAnalysisManager.h"
#include "AliEMCALDigit.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecPoint.h"
#include "AliGenerator.h"
#include "AliHeader.h"
#include "AliLog.h"
#include "AliPicoTrack.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliStack.h"
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
Bool_t AliJetEmbeddingFromGenTask::ExecOnce() 
{
  // Exec only once.

  if (!gAlice) {
    new AliRun("gAlice","The ALICE Off-line Simulation Framework");
    delete gRandom;
    gRandom = new TRandom3(0);
  }

  TFolder *folder = new TFolder(GetName(),GetName());
  AliRunLoader *rl = new AliRunLoader(folder);
  rl->MakeHeader();
  rl->MakeStack();
  AliStack *stack = rl->Stack();
  fGen->SetStack(stack);
  fGen->Init();

  if (!(InputEvent()->FindListObject(fTracksName))) {
    fOutTracks = new TClonesArray("AliPicoTrack", 1000);
    fOutTracks->SetName(fTracksName);
    InputEvent()->AddObject(fOutTracks);
    fNTracks = 0;
  }
  
  return kTRUE;
}

//________________________________________________________________________
void AliJetEmbeddingFromGenTask::Run() 
{
  // Embed particles.

  if(fOutTracks)
   fOutTracks->Delete();

  if (fCopyArray) 
    CopyTracks();

  AliStack *stack = fGen->GetStack();
  stack->Reset();
  fGen->Generate();
  const Int_t nprim = stack->GetNprimary();
  for (Int_t i=0;i<nprim;++i) {
    if (!stack->IsPhysicalPrimary(i))
      continue;
    TParticle *part = stack->Particle(i);
    TParticlePDG *pdg = part->GetPDG(1);
    if (!pdg) 
      continue;
    Int_t c = (Int_t)(TMath::Abs(pdg->Charge()));
    if (c==0) 
      continue;
    Double_t pt = part->Pt();
    Double_t eta = part->Eta();
    Double_t phi = part->Phi();
    if (eta<fEtaMin)
      continue;
    if (eta>fEtaMax)
      continue;
    if (phi<fPhiMin)
      continue;
    if (phi>fPhiMax)
      continue;
    if (pt<fPtMin)
      continue;
    if (pt>fPtMax)
      continue;
    AddTrack(pt, eta, phi,0,0,0,0,0,0,part->GetMass());
  }
}
