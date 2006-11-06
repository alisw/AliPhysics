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
  fdNdEtaAnalysisTr(0),
  fdNdEtaAnalysisTrVtx(0),
  fVertex(0),
  fPartPt(0)
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
  fdNdEtaAnalysisTr = new dNdEtaAnalysis("dndetaTr", "dndetaTr");
  fdNdEtaAnalysisTrVtx = new dNdEtaAnalysis("dndetaTrVtx", "dndetaTrVtx");
  fVertex = new TH3F("vertex_check", "vertex_check", 50, -50, 50, 50, -50, 50, 50, -50, 50);
  for (Int_t i=0; i<3; ++i)
  {
    fPartEta[i] = new TH1F(Form("dndeta_check_%d", i), Form("dndeta_check_%d", i), 120, -6, 6);
    fPartEta[i]->Sumw2();
  }
  fPartPt =  new TH1F("dndeta_check_pt", "dndeta_check_pt", 1000, 0, 10);
  fPartPt->Sumw2();
}

void AlidNdEtaAnalysisMCSelector::Init(TTree *tree)
{
  AliSelectorRL::Init(tree);

  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("fTriggerMask", 1);
  tree->SetBranchStatus("fSPDVertex*", 1);
}

Bool_t AlidNdEtaAnalysisMCSelector::Process(Long64_t entry)
{
  // fill the dNdEtaAnalysis class from the monte carlo

  if (AliSelectorRL::Process(entry) == kFALSE)
    return kFALSE;

  // Check prerequisites
  if (!fESD)
  {
    AliDebug(AliLog::kError, "ESD branch not available");
    return kFALSE;
  }

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

  Bool_t vertexReconstructed = AliPWG0Helper::IsVertexReconstructed(fESD);
  Bool_t eventTriggered = AliPWG0Helper::IsEventTriggered(fESD);

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

    if (eventTriggered)
    {
      fdNdEtaAnalysisTr->FillTrack(vtxMC[2], particle->Eta(), particle->Pt(), 1);
      if (vertexReconstructed)
        fdNdEtaAnalysisTrVtx->FillTrack(vtxMC[2], particle->Eta(), particle->Pt(), 1);
    }

    if (TMath::Abs(vtxMC[2]) < 10)
    {
      fPartEta[0]->Fill(particle->Eta());

      if (vtxMC[2] < 0)
        fPartEta[1]->Fill(particle->Eta());
      else
        fPartEta[2]->Fill(particle->Eta());
    }

    if (TMath::Abs(particle->Eta()) < 0.8)
      fPartPt->Fill(particle->Pt());
  }
  fdNdEtaAnalysis->FillEvent(vtxMC[2], 1);
  if (eventTriggered)
  {
    fdNdEtaAnalysisTr->FillEvent(vtxMC[2], 1);
    if (vertexReconstructed)
      fdNdEtaAnalysisTrVtx->FillEvent(vtxMC[2], 1);
  }

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
  fOutput->Add(fdNdEtaAnalysisTr);
  fOutput->Add(fdNdEtaAnalysisTrVtx);
  fOutput->Add(fPartPt);
  for (Int_t i=0; i<3; ++i)
    fOutput->Add(fPartEta[i]);
}

void AlidNdEtaAnalysisMCSelector::Terminate()
{
  //

  AliSelectorRL::Terminate();

  fdNdEtaAnalysis = dynamic_cast<dNdEtaAnalysis*> (fOutput->FindObject("dndeta"));
  fdNdEtaAnalysisTr = dynamic_cast<dNdEtaAnalysis*> (fOutput->FindObject("dndetaTr"));
  fdNdEtaAnalysisTrVtx = dynamic_cast<dNdEtaAnalysis*> (fOutput->FindObject("dndetaTrVtx"));
  fPartPt = dynamic_cast<TH1F*> (fOutput->FindObject("dndeta_check_pt"));
  for (Int_t i=0; i<3; ++i)
    fPartEta[i] = dynamic_cast<TH1F*> (fOutput->FindObject(Form("dndeta_check_%d", i)));

  if (!fdNdEtaAnalysis || !fdNdEtaAnalysisTr || !fdNdEtaAnalysisTrVtx || !fPartPt)
  {
    AliDebug(AliLog::kError, Form("ERROR: Histograms not available %p %p", (void*) fdNdEtaAnalysis, (void*) fPartPt));
    return;
  }

  fdNdEtaAnalysis->Finish(0, -1);
  fdNdEtaAnalysisTr->Finish(0, -1);
  fdNdEtaAnalysisTrVtx->Finish(0, -1);

  Int_t events = (Int_t) fdNdEtaAnalysis->GetVtxHistogram()->Integral();
  fPartPt->Scale(1.0/events);
  fPartPt->Scale(1.0/fPartPt->GetBinWidth(1));

  if (fPartEta[0])
  {
    TH1* vtxHist = fdNdEtaAnalysis->GetVtxHistogram();
    Int_t events1 = (Int_t) vtxHist->Integral(vtxHist->GetXaxis()->FindBin(-9.9), vtxHist->GetXaxis()->FindBin(-0.1));
    Int_t events2 = (Int_t) vtxHist->Integral(vtxHist->GetXaxis()->FindBin(0.1), vtxHist->GetXaxis()->FindBin(9.9));

    fPartEta[0]->Scale(1.0 / (events1 + events2));
    fPartEta[1]->Scale(1.0 / events1);
    fPartEta[2]->Scale(1.0 / events2);

    for (Int_t i=0; i<3; ++i)
      fPartEta[i]->Scale(1.0 / fPartEta[i]->GetBinWidth(1));

    new TCanvas("control", "control", 500, 500);
    for (Int_t i=0; i<3; ++i)
    {
      fPartEta[i]->SetLineColor(i+1);
      fPartEta[i]->Draw((i==0) ? "" : "SAME");
    }
  }

  TFile* fout = new TFile("analysis_mc.root","RECREATE");

  fdNdEtaAnalysis->SaveHistograms();
  fdNdEtaAnalysisTr->SaveHistograms();
  fdNdEtaAnalysisTrVtx->SaveHistograms();
  fPartPt->Write();

  fout->Write();
  fout->Close();

  if (fPartPt)
  {
    new TCanvas("control2", "control2", 500, 500);
    fPartPt->DrawCopy();
  }
}
