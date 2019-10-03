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
#include "AliVCluster.h"
#include "AliVEvent.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHerwigEventHeader.h"
#include "AliEmcalPythiaInfo.h"
#include "AliPythiaRndm.h"
#include "AliHerwigRndm.h"
ClassImp(AliJetEmbeddingFromGenTask)

//________________________________________________________________________
AliJetEmbeddingFromGenTask::AliJetEmbeddingFromGenTask() : 
  AliJetModelBaseTask("AliJetEmbeddingFromGenTask"),
  fGen(0),
  fMassless(kFALSE),
  fChargedOnly(kFALSE),
  fToyModelFragmentation(kFALSE),
  fToyModelFraction(0.9),
  fHistPt(0),
  fHistEtaPhi(0),
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
  fMassless(kFALSE),
  fChargedOnly(kFALSE),
  fToyModelFragmentation(kFALSE),
  fToyModelFraction(0.9),
  fHistPt(0),
  fHistEtaPhi(0),
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

  fHistPt = new TH1F("fHistpt","fHistPt;#it{p}_{T};N",100,0.,100.);
  fOutput->Add(fHistPt);

  fHistEtaPhi = new TH2F("fHistEtapHI","fHistEtaPhi;#eta;#varphi",100,-3.,3.,100.,0.,TMath::TwoPi());
  fOutput->Add(fHistEtaPhi);

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
  gAlice->SetRunLoader(rl);
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

  if(!fPythiaInfoName.IsNull()) {
    if (!(InputEvent()->FindListObject(fPythiaInfoName))) {
      fPythiaInfo = new AliEmcalPythiaInfo("PythiaInfo");
      fPythiaInfo->SetName(fPythiaInfoName);
      InputEvent()->AddObject(fPythiaInfo);
    }
  }
  return kTRUE;
}

//________________________________________________________________________
void AliJetEmbeddingFromGenTask::Run() 
{
  // Embed particles.

  if (fCopyArray) 
    CopyTracks();
  if(fGenType<3){
  AliPythiaRndm::SetPythiaRandom(new TRandom3());
  AliPythiaRndm::GetPythiaRandom()->SetSeed(clock()+gSystem->GetPid());}
  if(fGenType==3){
  AliHerwigRndm::SetHerwigRandom(new TRandom3());
  AliHerwigRndm::GetHerwigRandom()->SetSeed(clock()+gSystem->GetPid());}
 
  AliStack *stack = fGen->GetStack();
  stack->Reset();
  fGen->Generate();
  const Int_t nprim = stack->GetNprimary();
  // reject if partons are missing from stack for some reason
  if(nprim < 8) return;
  if(fPythiaInfo) {
    TParticle *part6 = stack->Particle(6);
    TParticle *part7 = stack->Particle(7);
    
    fPythiaInfo->SetPartonFlag6(TMath::Abs(part6->GetPdgCode()));
    fPythiaInfo->SetParton6(part6->Pt(), part6->Eta(), part6->Phi(), part6->GetMass());
    
    fPythiaInfo->SetPartonFlag7(TMath::Abs(part7->GetPdgCode()));
    fPythiaInfo->SetParton7(part7->Pt(), part7->Eta(), part7->Phi(), part7->GetMass());
  }

  for (Int_t i=0;i<nprim;++i) {
    if (!stack->IsPhysicalPrimary(i))
      continue;
    TParticle *part = stack->Particle(i);
    TParticlePDG *pdg = part->GetPDG(1);
    if (!pdg) 
      continue;
    Int_t c = (Int_t)(TMath::Abs(pdg->Charge()));
    if (fChargedOnly && c==0)  continue;
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
    Double_t mass = part->GetMass();
    if(fMassless) mass = 0.;
    if(fToyModelFragmentation) pt=pt*fToyModelFraction;
    fHistPt->Fill(pt);
    fHistEtaPhi->Fill(eta,phi);
    AddTrack(pt, eta, phi,0,0,0,0,0,0,c,mass);
  }

  FillPythiaHistograms();
}

//________________________________________________________________________
void AliJetEmbeddingFromGenTask::FillPythiaHistograms() {
  //Get PYTHIA info: pt-hard, x-section, trials

  if (!fQAhistos)
    return;

  AliRunLoader *rl = AliRunLoader::Instance();
  AliGenPythiaEventHeader *genPP=0x0;
  AliGenHerwigEventHeader *genPH=0x0;
  if(fGenType<3)  genPP = dynamic_cast<AliGenPythiaEventHeader*>(rl->GetHeader()->GenEventHeader());
  if(fGenType==3)  genPH = dynamic_cast<AliGenHerwigEventHeader*>(rl->GetHeader()->GenEventHeader());
  if(fGenType<3 && genPP) {
    Float_t xsec = genPP->GetXsection();
    Int_t trials = genPP->Trials();
    Float_t pthard = genPP->GetPtHard();
    Float_t ptWeight=genPP->EventWeight(); 
    
    if(fPythiaInfo) fPythiaInfo->SetPythiaEventWeight(ptWeight);
    fHistXsection->Fill(0.5,xsec);
    fHistTrials->Fill(0.5,trials);
    fHistPtHard->Fill(pthard);
  }

   if(fGenType==3 && genPH) {
    Float_t xsec = genPH->Weight();
    Int_t trials = genPH->Trials();
    Float_t ptWeight=genPH->Weight();  
    if(fPythiaInfo) fPythiaInfo->SetPythiaEventWeight(ptWeight);
    fHistXsection->Fill(0.5,xsec);
    fHistTrials->Fill(0.5,trials);

  }


}
