#define AliAnalysisTaskFemto_cxx
#include "TROOT.h"
#include "TChain.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TFile.h"

#include "AliAnalysisTask.h"

#include "AliESDEvent.h"

#include "AliAnalysisTaskFemto.h"

ClassImp(AliAnalysisTaskFemto)

//________________________________________________________________________
  AliAnalysisTaskFemto::AliAnalysisTaskFemto(const char *name) :AliAnalysisTask(name,""), fESD(0), fHistPt(0) {
  // Constructor.
  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
  //  DefineInput(1, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0, TH1F::Class());
}

//________________________________________________________________________
void AliAnalysisTaskFemto::ConnectInputData(Option_t *) {
  printf("   ConnectInputData %s\n", GetName());

//   char ** address = (char **)GetBranchAddress(0, "ESD");
//   cout << "Got tree address " << address << endl;
//   if (address) {
//     fESD = (AliESDEvent*)(*address);
//   }
//   else  {
//     fESD = new AliESDEvent();
//     fESD->
//     //    SetBranchAddress(0, "ESD", &fESD);
//   }
  TChain *inpchain = (TChain *) GetInputData(0);
  fESD = new AliESDEvent();
  fESD->ReadFromTree(inpchain);
  
  

//   char ** addressfriend = (char **)GetBranchAddress(1, "ESDfriend");
//   cout << "Got friend address " << addressfriend << endl;
//   if (addressfriend) {
//     fESDfriend = (AliESDfriend*)(*address);
//   }
//   else  {
//     fESDfriend = new AliESDfriend();
//     SetBranchAddress(1, "ESDfriend", &fESDfriend);
//   }


//   TString fileName(chain->GetCurrentFile()->GetName());
//   fileName.ReplaceAll("AliESDs", "AliESDfriends");
//   cout << "Reading friend " << fileName.Data() << endl;;
//   chain->AddFriend("esdFriendTree",fileName.Data());
//   SetBranchAddress(0, "ESDfriend",&fESDfriend);

}

//________________________________________________________________________
void AliAnalysisTaskFemto::CreateOutputObjects() {
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
void AliAnalysisTaskFemto::Exec(Option_t *) {
  // Task making a pt distribution.
  // Get input data
  TChain *chain = (TChain*)GetInputData(0);
  Long64_t ientry = chain->GetReadEntry();
  if (!fESD) return;
  //  if (!chain->GetCurrentFile())
    //    return 0;
    
  printf("Tracks: %d \n",fESD->GetNumberOfTracks());
  
  if (fESD->GetNumberOfTracks() > 0) {
    cout << "Chain has " << chain->GetListOfFriends() << " friends" << endl;
    cout << "First track is " << fESD->GetTrack(0)->Get1P() << " momentum" << endl;
    
    fReader->SetESDSource(fESD);
    //    fReader->SetESDfriendSource(*fESDfriend);
    
    fManager->ProcessEvent();
  
    // Post final data. It will be written to a file with option "RECREATE"
    PostData(0, fHistPt);
  }
}      

//________________________________________________________________________
void AliAnalysisTaskFemto::Terminate(Option_t *) {
  // Draw some histogram at the end.
//   if (!gROOT->IsBatch()) {
//     TCanvas *c1 = new TCanvas("c1","Pt",10,10,310,310);
//     c1->SetFillColor(10);
//     c1->SetHighLightColor(10);
    
//     c1->cd(1)->SetLeftMargin(0.15);
//     c1->cd(1)->SetBottomMargin(0.15);  
//     c1->cd(1)->SetLogy();
//     //fHistPt = (TH1F*)GetOutputData(0);
//     if (fHistPt) fHistPt->DrawCopy("E");
//   }
  
}
//________________________________________________________________________
void AliAnalysisTaskFemto::SetFemtoReader(AliFemtoEventReaderESDChain *aReader)
{
  fReader = aReader;
}
//________________________________________________________________________
void AliAnalysisTaskFemto::SetFemtoManager(AliFemtoManager *aManager)
{
  fManager = aManager;
}
//________________________________________________________________________
// void AliAnalysisTaskFemto::SetFriendAddress(AliESDfriend **aFriendAddress)
// {
//   fESDfriend = aFriendAddress;
// }

