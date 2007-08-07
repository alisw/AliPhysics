#define AliAnalysisTaskPt_cxx
#include "TROOT.h"
#include "TChain.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "Riostream.h"

#include "AliAnalysisTask.h"

#include "AliESD.h"

#include "AliAnalysisTaskPt.h"

ClassImp(AliAnalysisTaskPt)

//________________________________________________________________________
AliAnalysisTaskPt::AliAnalysisTaskPt(const char *name) :AliAnalysisTask(name,""), fESD(0), fHistPt(0) {
  // Constructor.
  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0, TH1F::Class());
}

//________________________________________________________________________
void AliAnalysisTaskPt::ConnectInputData(Option_t *) {
  printf("   ConnectInputData %s\n", GetName());

  char ** address = (char **)GetBranchAddress(0, "ESD");
  if (address) {
    fESD = (AliESD*)(*address);
  }
  else  {
    fESD = new AliESD();
    SetBranchAddress(0, "ESD", &fESD);
  }
}

//________________________________________________________________________
void AliAnalysisTaskPt::CreateOutputObjects() {
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
void AliAnalysisTaskPt::Exec(Option_t *) {
  // Task making a pt distribution.
  // Get input data
  TChain *chain = (TChain*)GetInputData(0);
  Long64_t ientry = chain->GetReadEntry();
  if (!fESD) return;

  cout<<"Entry: "<<ientry<<" - Tracks: "<<fESD->GetNumberOfTracks()<<endl;
  
  for(Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
    AliESDtrack * track = fESD->GetTrack(iTracks);
    //UInt_t status = track->GetStatus();
    Double_t momentum[3];
    track->GetPxPyPz(momentum);
    Double_t Pt = sqrt(pow(momentum[0],2) + pow(momentum[1],2));
    fHistPt->Fill(Pt);
  }//track loop 
  // Post final data. It will be written to a file with option "RECREATE"
  PostData(0, fHistPt);
}      

//________________________________________________________________________
void AliAnalysisTaskPt::Terminate(Option_t *) {
  // Draw some histogram at the end.
  if (!gROOT->IsBatch()) {
    TCanvas *c1 = new TCanvas("c1","Pt",10,10,510,510);
    c1->SetFillColor(10); c1->SetHighLightColor(10);
    c1->cd(1)->SetLeftMargin(0.15); c1->cd(1)->SetBottomMargin(0.15);  
    c1->cd(1)->SetLogy();
    fHistPt = (TH1F*)GetOutputData(0);
    if (fHistPt) fHistPt->DrawCopy("E");
  }
}
