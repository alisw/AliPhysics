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
#include <TFile.h>

#include <AliLog.h>
#include <AliGenEventHeader.h>
#include <AliHeader.h>

#include "dNdEta/dNdEtaAnalysis.h"
#include "AliPWG0Helper.h"


ClassImp(AlidNdEtaAnalysisMCSelector)

AlidNdEtaAnalysisMCSelector::AlidNdEtaAnalysisMCSelector() :
  AliSelectorRL(),
  fdNdEtaAnalysis(0),
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

void AlidNdEtaAnalysisMCSelector::SlaveBegin(TTree * tree)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  AliSelectorRL::SlaveBegin(tree);

  fdNdEtaAnalysis = new dNdEtaAnalysis("dndeta", "dndeta");
}

void AlidNdEtaAnalysisMCSelector::Init(TTree *tree)
{
  AliSelectorRL::Init(tree);

  tree->SetBranchStatus("ESD", 0);

  fVertex = new TH3F("vertex_check", "vertex_check", 50, -50, 50, 50, -50, 50, 50, -50, 50);
  fPartEta = new TH1F("dndeta_check", "dndeta_check", 120, -6, 6);
  fPartEta->Sumw2();
}

Bool_t AlidNdEtaAnalysisMCSelector::Process(Long64_t entry)
{
  // fill the dNdEtaAnalysis class from the monte carlo

  if (AliSelectorRL::Process(entry) == kFALSE)
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

    if (AliPWG0Helper::IsPrimaryCharged(particle, nPrim) == kFALSE)
      continue;

    AliDebug(AliLog::kDebug+1, Form("Accepted primary %d, unique ID: %d", i_mc, particle->GetUniqueID()));

    fdNdEtaAnalysis->FillTrack(vtxMC[2], particle->Eta(), particle->Pt(), 1);
    fVertex->Fill(particle->Vx(), particle->Vy(), particle->Vz());

    fPartEta->Fill(particle->Eta());
  }
  fdNdEtaAnalysis->FillEvent(vtxMC[2], 1);

  ++fEvents;

  return kTRUE;
}

void AlidNdEtaAnalysisMCSelector::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  AliSelectorRL::SlaveTerminate();

  // Add the histograms to the output on each slave server
  if (!fOutput)
  {
    AliDebug(AliLog::kError, Form("ERROR: Output list not initialized."));
    return;
  }

  fOutput->Add(fdNdEtaAnalysis);
}

void AlidNdEtaAnalysisMCSelector::Terminate()
{
  //

  AliSelectorRL::Terminate();

  fdNdEtaAnalysis = dynamic_cast<dNdEtaAnalysis*> (fOutput->FindObject("dndeta"));

  if (!fdNdEtaAnalysis)
  {
    AliDebug(AliLog::kError, Form("ERROR: Histograms not available %p", (void*) fdNdEtaAnalysis));
    return;
  }

  TFile* fout = new TFile("analysis_mc.root","RECREATE");

  fdNdEtaAnalysis->SaveHistograms();

  fout->Write();
  fout->Close();

  fPartEta->Scale(1.0/fEvents);
  fPartEta->Scale(1.0/fPartEta->GetBinWidth(1));

  TCanvas* canvas = new TCanvas("control", "control", 900, 450);
  canvas->Divide(2, 1);

  canvas->cd(1);
  fVertex->Draw();

  canvas->cd(2);
  fPartEta->Draw();
}
