#define AliAnalysisTaskRLPt_cxx

#include "TChain.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TParticle.h"

#include "AliESD.h"
#include "AliLog.h"
#include "AliStack.h"

#include "AliAnalysisTaskRL.h"
#include "AliAnalysisTaskRLPt.h"

ClassImp(AliAnalysisTaskRLPt)

//________________________________________________________________________
AliAnalysisTaskRLPt::AliAnalysisTaskRLPt(const char *name) :AliAnalysisTaskRL(name,""), fESD(0), fHistPt(0) {
  // Constructor.
  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0, TH1F::Class());
}

//___________________________________________________________________________
void AliAnalysisTaskRLPt::ConnectInputData(Option_t *) {
  // Initialize branches.
  printf("   ConnectInputData of task %s\n", GetName());
  if (!fESD) {
    char ** address = (char **)GetBranchAddress(0, "ESD");
    if (address) fESD = (AliESD*)(*address);
    if (!fESD) {
      fESD = new AliESD();
      SetBranchAddress(0, "ESD", &fESD);
    }
  }
}

//___________________________________________________________________________
void AliAnalysisTaskRLPt::CreateOutputObjects() {
  printf("   CreateOutputObjects of task %s\n", GetName());
  if (!fHistPt) {
    fHistPt = new TH1F("fHistPt","This is the Pt distribution",15,0.1,3.1);
    fHistPt->SetStats(kTRUE);
    fHistPt->GetXaxis()->SetTitle("P_{T} [GeV]");
    fHistPt->GetYaxis()->SetTitle("#frac{dN}{dP_{T}}");
    fHistPt->GetXaxis()->SetTitleColor(1);
    fHistPt->SetMarkerStyle(kFullCircle);
  }
}

//________________________________________________________________________
void AliAnalysisTaskRLPt::Exec(Option_t *) {
  // Task making a pt distribution.
  // Get input data
  TTree *tinput = (TTree*)GetInputData(0);
  Long64_t ientry = tinput->GetReadEntry();
  if (AliAnalysisTaskRL::GetEntry(ientry) == kFALSE) {
    printf("Couldn't get event from the runLoader\n");
    return;
  }  
  if (!fESD) return;

  AliStack* stack = GetStack();
  if (!stack) {
  AliDebug(AliLog::kError, "Stack not available");
    //return kFALSE;
  }
  // loop over mc particles
  Int_t nPrim = stack->GetNprimary();
  printf("Particles: %d - Tracks: %d \n",nPrim,fESD->GetNumberOfTracks());
  for(Int_t i = 0; i < nPrim; i++) {
    TParticle * particle = stack->Particle(i); 
    if(TMath::Abs(particle->Eta()) > 1.0) continue;
    fHistPt->Fill(particle->Pt());
  }

  // Post final data. It will be written to a file with option "RECREATE"
  PostData(0, fHistPt);
}      

//________________________________________________________________________
void AliAnalysisTaskRLPt::Terminate(Option_t *) {
  // Draw some histogram at the end.
  if (!gROOT->IsBatch()) {
    TCanvas *c1 = new TCanvas("c1","Pt",10,10,310,310);
    c1->SetFillColor(10); c1->SetHighLightColor(10);
    c1->cd(1)->SetLeftMargin(0.15); c1->cd(1)->SetBottomMargin(0.15);  
    c1->cd(1)->SetLogy();
    fHistPt->DrawCopy("E");
  }
}
