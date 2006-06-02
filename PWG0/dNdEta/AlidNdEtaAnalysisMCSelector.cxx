/* $Id$ */

#include "AlidNdEtaAnalysisMCSelector.h"

#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TVector3.h>
#include <TH1F.h>
#include <TH3F.h>
#include <TTree.h>

#include <AliLog.h>
#include <AliGenEventHeader.h>
#include <AliHeader.h>

#include "dNdEtaAnalysis.h"


ClassImp(AlidNdEtaAnalysisMCSelector)

AlidNdEtaAnalysisMCSelector::AlidNdEtaAnalysisMCSelector() :
  AlidNdEtaAnalysisSelector(),
  fVertex(0),
  fPartEta(0),
  fEvents(0)
{
  //
  // Constructor. Initialization of pointers
  //
}

AlidNdEtaAnalysisMCSelector::~AlidNdEtaAnalysisMCSelector()
{
  //
  // Destructor
  //
}

void AlidNdEtaAnalysisMCSelector::Init(TTree *tree)
{
   AlidNdEtaAnalysisSelector::Init(tree);

  tree->SetBranchStatus("ESD", 0);

  fVertex = new TH3F("vertex_check", "vertex_check", 50, -50, 50, 50, -50, 50, 50, -50, 50);
  fPartEta = new TH1F("dndeta_check", "dndeta_check", 120, -6, 6);
  fPartEta->Sumw2();
}

Bool_t AlidNdEtaAnalysisMCSelector::Process(Long64_t entry)
{
  //

  if (AliSelector::Process(entry) == kFALSE)
    return kFALSE;

  TTree* particleTree = GetKinematics();
  if (!particleTree)
  {
    AliDebug(AliLog::kError, "Kinematics not available");
    return kFALSE;
  }

  AliHeader* header = GetHeader();
  if (!header)
  {
    AliDebug(AliLog::kError, "Header not available");
    return kFALSE;
  }

  // get the MC vertex
  AliGenEventHeader* genHeader = header->GenEventHeader();

  TArrayF vtxMC(3);
  genHeader->PrimaryVertex(vtxMC);

  particleTree->SetBranchStatus("*", 0);
  particleTree->SetBranchStatus("fDaughter[2]", 1);
  particleTree->SetBranchStatus("fPdgCode", 1);
  particleTree->SetBranchStatus("fPx", 1);
  particleTree->SetBranchStatus("fPy", 1);
  particleTree->SetBranchStatus("fPz", 1);
  particleTree->SetBranchStatus("fVx", 1);
  particleTree->SetBranchStatus("fVy", 1);
  particleTree->SetBranchStatus("fVz", 1);

  TParticle* particle = 0;
  particleTree->SetBranchAddress("Particles", &particle);

  Int_t nPrim  = header->GetNprimary();
  Int_t nTotal = header->GetNtrack();

  for (Int_t i_mc = nTotal - nPrim; i_mc < nTotal; ++i_mc)
  {
    particleTree->GetEntry(i_mc);

    if (!particle)
      continue;

    if (IsPrimaryCharged(particle, nPrim) == kFALSE)
      continue;

    AliDebug(AliLog::kDebug+1, Form("Accepted primary %d, unique ID: %d", i_mc, particle->GetUniqueID()));

    fdNdEtaAnalysis->FillTrack(vtxMC[2], particle->Eta());
    fVertex->Fill(particle->Vx(), particle->Vy(), particle->Vz());

    fPartEta->Fill(particle->Eta());
  }
  fdNdEtaAnalysis->FillEvent(vtxMC[2]);

  ++fEvents;

  return kTRUE;
}

void AlidNdEtaAnalysisMCSelector::Terminate()
{
  AlidNdEtaAnalysisSelector::Terminate();

  fPartEta->Scale(1.0/fEvents);
  fPartEta->Scale(1.0/fPartEta->GetBinWidth(1));

  TCanvas* canvas = new TCanvas("control", "control", 900, 450);
  canvas->Divide(2, 1);

  canvas->cd(1);
  fVertex->Draw();

  canvas->cd(2);
  fPartEta->Draw();
}
