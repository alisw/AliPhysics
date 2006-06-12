/* $Id$ */

#include "AlidNdEtaVertexRecEffSelector.h"

#include <TCanvas.h>
#include <TH1F.h>
#include <TTree.h>
#include <TParticle.h>
#include <TFile.h>

#include <AliLog.h>
#include <AliGenEventHeader.h>
#include <AliHeader.h>
#include <AliESD.h>
#include <AliESDVertex.h>

//
// This class plots the vertex reconstruction efficiency
// If a vertex was reconstructed is decided by the function CheckVertex()
// In any case the *generated* multiplicity is filled into the histogram
//

ClassImp(AlidNdEtaVertexRecEffSelector)

const Float_t AlidNdEtaVertexRecEffSelector::fkEtaRange = 0.9;

AlidNdEtaVertexRecEffSelector::AlidNdEtaVertexRecEffSelector() :
  AliSelectorRL(),
  fdNGen(0),
  fdNRec(0),
  fVtxGen(0),
  fVtxRec(0)
{
  //
  // Constructor. Initialization of pointers
  //
}

AlidNdEtaVertexRecEffSelector::~AlidNdEtaVertexRecEffSelector()
{
  //
  // Destructor
  //

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
}

void AlidNdEtaVertexRecEffSelector::SlaveBegin(TTree * tree)
{
  // initializes the histograms

  AliSelectorRL::SlaveBegin(tree);

  fdNGen = new TH1F("dNGen", "dNGen", 90, 0, 50);
  fdNRec = dynamic_cast<TH1F*>(fdNGen->Clone("dNRec"));

  fdNGen->SetTitle("Generated Events;dN_{Gen};Count");
  fdNRec->SetTitle("Events with reconstructed vertex;dN_{Gen};Count");
  
  fVtxGen = new TH1F("VtxGen", "VtxGen", 200, -20, 20);
  fVtxRec = dynamic_cast<TH1F*>(fVtxGen->Clone("VtxRec"));
}

Bool_t AlidNdEtaVertexRecEffSelector::CheckVertex()
{
  //
  // check if the vertex has been reconstructed well enough
  //  
  if (!fESD)
    return kFALSE;

  const AliESDVertex* vtxESD = fESD->GetVertex();

  // the vertex should be reconstructed
  if (strcmp(vtxESD->GetName(),"default")==0)
    return kFALSE;

  Double_t vtxRes[3];
  vtxRes[0] = vtxESD->GetXRes();
  vtxRes[1] = vtxESD->GetYRes();
  vtxRes[2] = vtxESD->GetZRes();

  // the resolution should be reasonable???
  if (vtxRes[2]==0 || vtxRes[2]>0.1)
    return kFALSE;

  return kTRUE;
}

Bool_t AlidNdEtaVertexRecEffSelector::Process(Long64_t entry)
{
  //
  // fills fdNGen and fdNRec
  //

  if (AliSelectorRL::Process(entry) == kFALSE)
    return kFALSE;

  // check prerequisites
  if (!fESD)
  {
    AliDebug(AliLog::kError, "ESD branch not available");
    return kFALSE;
  }

  AliHeader* header = GetHeader();
  if (!header)
  {
    AliDebug(AliLog::kError, "Header not available");
    return kFALSE;
  }


  // loop over mc particles
  TTree* particleTree = GetKinematics();
  TParticle* particle = 0;
  particleTree->SetBranchAddress("Particles", &particle);

  Int_t nPrim  = header->GetNprimary();
  Int_t nTotal = header->GetNtrack();

  Int_t n = 0;

  for (Int_t iMc = nTotal - nPrim; iMc < nTotal; ++iMc)
  {
    particleTree->GetEntry(iMc);

    if (!particle)
      continue;

    if (IsPrimaryCharged(particle, nPrim) == kFALSE)
      continue;

    if (TMath::Abs(particle->Eta()) < fkEtaRange)
      ++n;
  }// end of mc particle

  Float_t dN = (Float_t) n / (fkEtaRange*2);

  fdNGen->Fill(dN);

  // check vertex reconstruction
  if (CheckVertex() != kFALSE)
    fdNRec->Fill(dN);

  AliGenEventHeader* genHeader = header->GenEventHeader();

  TArrayF vtxMC(3);
  genHeader->PrimaryVertex(vtxMC);
  
  fVtxGen->Fill(vtxMC[2]);

  if (CheckVertex() != kFALSE)
    fVtxRec->Fill(vtxMC[2]);

  return kTRUE;
}

void AlidNdEtaVertexRecEffSelector::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  AliSelectorRL::SlaveTerminate();

  // Add the histograms to the output on each slave server
  if (!fOutput)
  {
    AliDebug(AliLog::kError, "ERROR: Output list not initialized");
    return;
  }

  fOutput->Add(fdNGen);
  fOutput->Add(fdNRec);
  fOutput->Add(fVtxGen);
  fOutput->Add(fVtxRec);
}

void AlidNdEtaVertexRecEffSelector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  AliSelectorRL::Terminate();

  if (!fOutput)
  {
    AliDebug(AliLog::kError, "ERROR: Output list not initialized");
    return;
  }

  fdNGen = dynamic_cast<TH1F*> (fOutput->FindObject("dNGen"));
  fdNRec = dynamic_cast<TH1F*> (fOutput->FindObject("dNRec"));
  fVtxGen = dynamic_cast<TH1F*> (fOutput->FindObject("VtxGen"));
  fVtxRec = dynamic_cast<TH1F*> (fOutput->FindObject("VtxRec"));
  if (!fdNGen || !fdNRec || !fVtxGen || !fVtxRec)
  {
    AliDebug(AliLog::kError, Form("ERROR: Histograms not available %p %p %p %p", (void*) fdNGen, (void*) fdNRec, (void*) fVtxGen, (void*) fVtxRec));
    return;
  }

  TFile* fout = new TFile("vertexRecEff.root","RECREATE");

  fdNGen->Write();
  fdNRec->Write();

  fVtxGen->Write();
  fVtxRec->Write();

  fout->Write();
  fout->Close();

  TCanvas* canvas = new TCanvas("dN", "dN", 900, 450);
  canvas->Divide(2, 1);

  canvas->cd(1);
  fdNGen->Draw();
  fdNRec->SetLineColor(kRed);
  fdNRec->Draw("SAME");

  canvas->cd(2);
  fVtxGen->Draw();
  fVtxRec->SetLineColor(kRed);
  fVtxRec->Draw("SAME");
}
