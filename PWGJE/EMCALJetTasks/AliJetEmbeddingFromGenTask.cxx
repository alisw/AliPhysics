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
#include <TProfile.h>
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
#include "AliGenPythiaEventHeader.h"

ClassImp(AliJetEmbeddingFromGenTask)

//________________________________________________________________________
AliJetEmbeddingFromGenTask::AliJetEmbeddingFromGenTask() : 
  AliJetModelBaseTask("AliJetEmbeddingFromGenTask"),
  fGen(0),
  fHistTrials(0),
  fHistXsection(0),
  fHistPtHard(0)
{
  // Default constructor.
  SetSuffix("EmbeddedFromGen");
}

//________________________________________________________________________
AliJetEmbeddingFromGenTask::AliJetEmbeddingFromGenTask(const char *name, Bool_t drawqa) :
  AliJetModelBaseTask(name,drawqa),
  fGen(0),
  fHistTrials(0),
  fHistXsection(0),
  fHistPtHard(0)
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
void AliJetEmbeddingFromGenTask::UserCreateOutputObjects()
{
  // Create user output.

  if (!fQAhistos)
    return;

  AliJetModelBaseTask::UserCreateOutputObjects();

  fHistTrials = new TH1F("fHistTrials", "fHistTrials", 1, 0, 1);
  fHistTrials->GetYaxis()->SetTitle("trials");
  fOutput->Add(fHistTrials);

  fHistXsection = new TProfile("fHistXsection", "fHistXsection", 1, 0, 1);
  fHistXsection->GetYaxis()->SetTitle("xsection");
  fOutput->Add(fHistXsection);

  fHistPtHard = new TH1F("fHistPtHard", "fHistPtHard", 500, 0., 500.);
  fHistPtHard->GetXaxis()->SetTitle("p_{T,hard} (GeV/c)");
  fHistPtHard->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistPtHard);

  PostData(1, fOutput);
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

  FillPythiaHistograms();
}

//________________________________________________________________________
void AliJetEmbeddingFromGenTask::FillPythiaHistograms() {
  //Get PYTHIA info: pt-hard, x-section, trials

  if (!fQAhistos)
    return;

  AliRunLoader *rl = AliRunLoader::Instance();
  AliGenPythiaEventHeader *genPH = dynamic_cast<AliGenPythiaEventHeader*>(rl->GetHeader()->GenEventHeader());
  if(genPH) {
    Printf("found pythia event header. pThard: %f Trials: %d xsec: %f",genPH->GetPtHard(),genPH->Trials(),genPH->GetXsection());

    Float_t xsec = genPH->GetXsection();
    Int_t trials = genPH->Trials();
    Float_t pthard = genPH->GetPtHard();

    fHistXsection->Fill(0.5,xsec);
    fHistTrials->Fill(0.5,trials);
    fHistPtHard->Fill(pthard);
  }
}
