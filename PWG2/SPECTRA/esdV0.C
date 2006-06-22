#define esdV0_cxx
// The class definition in esdV0.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called everytime a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("esdV0.C")
// Root > T->Process("esdV0.C","some options")
// Root > T->Process("esdV0.C+")
//

#include "esdV0.h"
#include <TPDGCode.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "AliPID.h"

esdV0::esdV0(TTree *) :
  TSelector(),
  fChain(0),
  fESD(0),
  fHistMassK0(0),
  fHistMassLambda(0),
  fHistMassAntiLambda(0),
  fHistMassLambdaVsProb(0),
  fHistMassLambdaCut(0),
  fHistMassAntiLambdaCut(0)
{
  // Constructor. Initialization of pointers
}

esdV0::~esdV0() {
  // Remove all pointers

  delete fESD;

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
}

void esdV0::Begin(TTree *)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
}

void esdV0::SlaveBegin(TTree * tree)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  Init(tree);
  
  TString option = GetOption();

  // Create histograms on each slave server
  fHistMassK0 = new TH1F("hMassK0", "K^{0} candidates", 100, 0.4, 0.6);
  fHistMassK0->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  fHistMassK0->SetOption("E");
  fHistMassK0->SetMarkerStyle(kFullCircle);
  fHistMassK0->SetStats(kFALSE);

  fHistMassLambda = new TH1F("hMassLambda", "#Lambda^{0} candidates", 75, 1.05, 1.2);
  fHistMassLambda->GetXaxis()->SetTitle("M(p#pi^{-}) [GeV/c^{2}]");
  fHistMassLambda->SetOption("E");
  fHistMassLambda->SetMarkerStyle(kFullCircle);
  fHistMassLambda->SetMarkerColor(kRed);
  fHistMassLambda->SetFillColor(4);
  fHistMassLambda->SetStats(kFALSE);

  fHistMassAntiLambda = new TH1F("hMassAntiLambda", "#bar{#Lambda}^{0} candidates", 75, 1.05, 1.2);
  fHistMassAntiLambda->GetXaxis()->SetTitle("M(#bar{p}#pi^{+}) [GeV/c^{2}]");
  fHistMassAntiLambda->SetOption("E");
  fHistMassAntiLambda->SetMarkerStyle(kFullCircle);
  fHistMassAntiLambda->SetMarkerColor(kRed);
  fHistMassAntiLambda->SetFillColor(4);
  fHistMassAntiLambda->SetStats(kFALSE);

  fHistMassLambdaVsProb = new TH2F("hMassLambdaVsProb", "#Lambda^{0} and #bar{#Lambda}^{0} candidates", 21, -2.5, 102.5, 75, 1.05, 1.2);
  fHistMassLambdaVsProb->GetXaxis()->SetTitle("prob(p) [%]");
  fHistMassLambdaVsProb->GetYaxis()->SetTitle("M(p#pi) [GeV/c^{2}]");
  fHistMassLambdaVsProb->SetOption("BOX");
  fHistMassLambdaVsProb->SetFillStyle(0);
  fHistMassLambdaVsProb->SetStats(kFALSE);

  fHistMassLambdaCut = new TH1F("hMassLambdaCut", "#Lambda^{0} candidates", 75, 1.05, 1.2);
  fHistMassLambdaCut->GetXaxis()->SetTitle("M(p#pi^{-}) [GeV/c^{2}]");
  fHistMassLambdaCut->SetOption("E");
  fHistMassLambdaCut->SetMarkerStyle(kFullCircle);
  fHistMassLambdaCut->SetMarkerColor(kRed);
  fHistMassLambdaCut->SetStats(kFALSE);

  fHistMassAntiLambdaCut = new TH1F("hMassAntiLambdaCut", "#bar{#Lambda}^{0} candidates", 75, 1.05, 1.2);
  fHistMassAntiLambdaCut->GetXaxis()->SetTitle("M(#bar{p}#pi^{+}) [GeV/c^{2}]");
  fHistMassAntiLambdaCut->SetOption("E");
  fHistMassAntiLambdaCut->SetMarkerStyle(kFullCircle);
  fHistMassAntiLambdaCut->SetMarkerColor(kRed);
  fHistMassAntiLambdaCut->SetStats(kFALSE);
}

void esdV0::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses of the tree
  // will be set. It is normaly not necessary to make changes to the
  // generated code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running with PROOF.

  // Set branch addresses
  if (tree == 0) return;
  fChain = tree;

  fChain->SetBranchAddress("ESD",&fESD);
}

Bool_t esdV0::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. Typically here the branch pointers
  // will be retrieved. It is normaly not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed.
  
  return kTRUE;
}


Bool_t esdV0::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either TTree::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing
  // keyed objects with PROOF, the object is already loaded and is available
  // via the fObject pointer.
  //
  // This function should contain the "body" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.

  // WARNING when a selector is used with a TChain, you must use
  //  the pointer to the current TTree to call GetEntry(entry).
  //  The entry is always the local entry number in the current tree.
  //  Assuming that fChain is the pointer to the TChain being processed,
  //  use fChain->GetEntry(entry).

  fChain->GetEntry(entry);

  if (!fESD) return kFALSE;

  Int_t nv0s = fESD->GetNumberOfV0s();

  Double_t mass;
  Double_t pPrior[5] = {0, 0, 0.85, 0.1, 0.05};
  Double_t p[5];
  Double_t pSum, pNorm;
  for (Int_t iV0 = 0; iV0 < fESD->GetNumberOfV0s(); iV0++) {
    AliESDv0* v0 = fESD->GetV0(iV0);
    if (!v0) continue;

    fHistMassK0->Fill(v0->GetEffMass());

    v0->ChangeMassHypothesis(kLambda0);
    mass = v0->GetEffMass();
    fESD->GetTrack(v0->GetPindex())->GetESDpid(p);
    pSum = 0;
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) pSum += p[i] * pPrior[i];
    if (pSum <= 0) pSum = 1.;
    pNorm = p[AliPID::kProton] * pPrior[AliPID::kProton] / pSum;
    fHistMassLambdaVsProb->Fill(100.*pNorm, mass); 
    if (pNorm > 0.1) fHistMassLambda->Fill(mass);
    
    v0->ChangeMassHypothesis(kLambda0Bar);
    mass = v0->GetEffMass();
    fESD->GetTrack(v0->GetNindex())->GetESDpid(p);
    pSum = 0;
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) pSum += p[i] * pPrior[i];
    if (pSum <= 0) pSum = 1.;
    pNorm = p[AliPID::kProton] * pPrior[AliPID::kProton] / pSum;
    fHistMassLambdaVsProb->Fill(100.*pNorm, mass); 
    if (pNorm > 0.1) fHistMassAntiLambda->Fill(mass);
   
  }//loop v0s
  
  return kTRUE;
}



void esdV0::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  // Add the histograms to the output on each slave server
  fOutput->Add(fHistMassK0);
  fOutput->Add(fHistMassLambda);
  fOutput->Add(fHistMassAntiLambda);
  fOutput->Add(fHistMassLambdaVsProb);
  fOutput->Add(fHistMassLambdaCut);
  fOutput->Add(fHistMassAntiLambdaCut);
}

void esdV0::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  fHistMassK0 = dynamic_cast<TH1F*>(fOutput->FindObject("hMassK0"));
  fHistMassLambda = dynamic_cast<TH1F*>(fOutput->FindObject("hMassLambda"));
  fHistMassAntiLambda = dynamic_cast<TH1F*>(fOutput->FindObject("hMassAntiLambda"));
  fHistMassLambdaVsProb = dynamic_cast<TH2F*>(fOutput->FindObject("hMassLambdaVsProb"));
  fHistMassLambdaCut = dynamic_cast<TH1F*>(fOutput->FindObject("hMassLambdaCut"));
  fHistMassAntiLambdaCut = dynamic_cast<TH1F*>(fOutput->FindObject("hMassAntiLambdaCut"));

  TFile* file = TFile::Open("V0.root", "RECREATE");
  fHistMassK0->Write();
  fHistMassLambda->Write();
  fHistMassAntiLambda->Write();
  fHistMassLambdaVsProb->Write();
  fHistMassLambdaCut->Write();
  fHistMassAntiLambdaCut->Write();
  file->Close();
  delete file;

  if (!gROOT->IsBatch()) 
    {
      TCanvas *c1 = new TCanvas("c1","V0s",10,10,610,610);
      c1->SetFillColor(10);
      c1->SetHighLightColor(10);
      c1->Divide(2,1);
 
      c1->cd(1);
      fHistMassK0->DrawCopy("E");
      fHistMassK0->Fit("gaus","q","",0.49,0.51);
      TVirtualPad* pad = (TVirtualPad*)c1->cd(2);
      pad->Divide(1,2);
      pad->cd(1);
      fHistMassLambda->DrawCopy("E");
      pad->cd(2);
      fHistMassAntiLambda->DrawCopy("E");
    }
}
