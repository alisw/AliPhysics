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
#include <AliStack.h>

#include "dNdEta/dNdEtaAnalysis.h"
#include "AliPWG0Helper.h"


ClassImp(AlidNdEtaAnalysisMCSelector)

AlidNdEtaAnalysisMCSelector::AlidNdEtaAnalysisMCSelector() :
  AliSelectorRL(),
  fdNdEtaAnalysis(0),
  fVertex(0),
  fPartEta(0),
  fPartPt(0),
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
  fVertex = new TH3F("vertex_check", "vertex_check", 50, -50, 50, 50, -50, 50, 50, -50, 50);
  fPartEta = new TH1F("dndeta_check", "dndeta_check", 120, -6, 6);
  fPartPt =  new TH1F("dndeta_check_pt", "dndeta_check_pt", 1000, 0, 10);
  fPartEta->Sumw2();
}

void AlidNdEtaAnalysisMCSelector::Init(TTree *tree)
{
  AliSelectorRL::Init(tree);

  tree->SetBranchStatus("*", 0);
}

Bool_t AlidNdEtaAnalysisMCSelector::Process(Long64_t entry)
{
  // fill the dNdEtaAnalysis class from the monte carlo

  if (AliSelectorRL::Process(entry) == kFALSE)
    return kFALSE;

  AliStack* stack = GetStack();
  if (!stack)
  {
    AliDebug(AliLog::kError, "Stack not available");
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

  // loop over mc particles
  Int_t nPrim  = stack->GetNprimary();

  for (Int_t iMc = 0; iMc < nPrim; ++iMc)
  {
    TParticle* particle = stack->Particle(iMc);

    if (!particle)
      continue;

    if (AliPWG0Helper::IsPrimaryCharged(particle, nPrim) == kFALSE)
      continue;

    AliDebug(AliLog::kDebug+1, Form("Accepted primary %d, unique ID: %d", iMc, particle->GetUniqueID()));

    fdNdEtaAnalysis->FillTrack(vtxMC[2], particle->Eta(), particle->Pt(), 1);
    fVertex->Fill(particle->Vx(), particle->Vy(), particle->Vz());

    fPartEta->Fill(particle->Eta());

    if (TMath::Abs(particle->Eta()) < 0.8)
      fPartPt->Fill(particle->Pt());
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
  fOutput->Add(fPartPt);
}

void AlidNdEtaAnalysisMCSelector::Terminate()
{
  //

  AliSelectorRL::Terminate();

  fdNdEtaAnalysis = dynamic_cast<dNdEtaAnalysis*> (fOutput->FindObject("dndeta"));
  fPartPt = dynamic_cast<TH1F*> (fOutput->FindObject("dndeta_check_pt"));

  if (!fdNdEtaAnalysis || !fPartPt)
  {
    AliDebug(AliLog::kError, Form("ERROR: Histograms not available %p %p", (void*) fdNdEtaAnalysis, (void*) fPartPt));
    return;
  }

  fdNdEtaAnalysis->Finish(0, -1);

  TFile* fout = new TFile("analysis_mc.root","RECREATE");

  fdNdEtaAnalysis->SaveHistograms();

  fout->Write();
  fout->Close();

  if (fPartEta)
  {
    fPartEta->Scale(1.0/fEvents);
    fPartEta->Scale(1.0/fPartEta->GetBinWidth(1));

    TCanvas* canvas = new TCanvas("control", "control", 900, 450);
    canvas->Divide(2, 1);

    canvas->cd(1);
    fVertex->Draw();

    canvas->cd(2);
    fPartEta->Draw();
  }

  if (fPartPt)
  {
    fPartPt->Scale(1.0/fEvents);
    fPartPt->Scale(1.0/fPartPt->GetBinWidth(1));

    new TCanvas("control2", "control2", 500, 500);
    fPartPt->Draw();
  }
}
