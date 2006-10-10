// The following methods are defined in this file:
//    Begin():        called everytime a loop on the tree starts, a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the slave servers.
//    Process():      called for each event, in this function you decide what to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree, a convenient place to draw/fit your histograms.

#include "EsdQa.h"    //class header
#include <TCanvas.h>  //Terminate()
#include <TChain.h>
#include <TF1.h>
#include <TBenchmark.h>
#include <TH2F.h>
#include <fstream>    //caf()      
#include <TProof.h>   //caf()
#include <TDSet.h>    //caf()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void EsdQa::Begin(TTree *)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void EsdQa::SlaveBegin(TTree *tree)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   Init(tree);

   TString option = GetOption();

   // create histograms on each slave server
   fCkovP    = new TH2F("CkovP" , "#theta_{c}, [rad];P, [GeV]", 150,   0,  7  ,100, 0, 1); 
   fSigP     = new TH2F("SigP"  ,"#sigma_{#theta_c}"          , 150,   0,  7  ,100, 0, 1e20);
   fMipXY    = new TH2F("MipXY" ,"mip position"               , 260,   0,130  ,252,0,126); 
   fDifXY    = new TH2F("DifXY" ,"diff"                       , 260, -10, 10  ,252,-10,10); 

   fProb[0] = new TH1F("PidE" ,"PID: e yellow #mu magenta"  ,100,0,1); fProb[0]->SetLineColor(kYellow);
   fProb[1] = new TH1F("PidMu","pid of #mu"                 ,100,0,1); fProb[1]->SetLineColor(kMagenta);
   fProb[2] = new TH1F("PidPi","PID: #pi red K green p blue",100,0,1); fProb[2]->SetLineColor(kRed);
   fProb[3] = new TH1F("PidK" ,"pid of K"                   ,100,0,1); fProb[3]->SetLineColor(kGreen);
   fProb[4] = new TH1F("PidP" ,"pid of p"                   ,100,0,1); fProb[4]->SetLineColor(kBlue);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void EsdQa::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses of the tree
   // will be set. It is normaly not necessary to make changes to the
   // generated code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running with PROOF.

   // Set branch addresses
   if ( !tree ) 
     return ;
   fChain = tree ;
   fChain->SetBranchAddress("ESD", &fEsd) ;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t EsdQa::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. Typically here the branch pointers
   // will be retrieved. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed.

   // Get branch pointers

   return kTRUE;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t EsdQa::Process(Long64_t entry)
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
   //  use fChain->GetTree()->GetEntry(entry).

  fChain->GetTree()->GetEntry(entry);

  for(Int_t iTrk=0;iTrk<fEsd->GetNumberOfTracks();iTrk++){
     AliESDtrack *pTrk=fEsd->GetTrack(iTrk);
     
     if(pTrk->GetRICHsignal()<0) continue;
     
     fCkovP->Fill(pTrk->GetP(),pTrk->GetRICHsignal()) ; 
     fSigP ->Fill(pTrk->GetP(),TMath::Sqrt(pTrk->GetRICHchi2()));
     
     Float_t xm,ym;  pTrk->GetRICHmipXY(xm,ym);  fMipXY->Fill(xm,ym);
     Float_t xd,yd;  pTrk->GetRICHdxdy(xd,yd);   fDifXY->Fill(xd,yd);
     
     Double_t pid[5];  pTrk->GetRICHpid(pid); for(Int_t i =0;i<5;i++) fProb[i]->Fill(pid[i]);
  }//tracks loop 
     
  return kTRUE;
}//Process()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void EsdQa::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
  
  // Add the histograms to the output on each slave server
  
  fOutput->Add(fCkovP);
  fOutput->Add(fSigP); 
  fOutput->Add(fMipXY);
  fOutput->Add(fDifXY);
  
  for(Int_t i=0;i<5;i++) fOutput->Add(fProb[i]);
}//SlaveTerminate()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void EsdQa::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
  
  fCkovP   = dynamic_cast<TH2F*>(fOutput->FindObject("CkovP")) ;
  fSigP    = dynamic_cast<TH2F*>(fOutput->FindObject("SigP")) ; 
  fMipXY   = dynamic_cast<TH2F*>(fOutput->FindObject("MipXY")) ;
  fDifXY   = dynamic_cast<TH2F*>(fOutput->FindObject("DifXY")) ;
  
  fProb[0] = dynamic_cast<TH1F*>(fOutput->FindObject("PidE")) ;
  fProb[1] = dynamic_cast<TH1F*>(fOutput->FindObject("PidMu")) ;
  fProb[2] = dynamic_cast<TH1F*>(fOutput->FindObject("PidPi")) ;
  fProb[3] = dynamic_cast<TH1F*>(fOutput->FindObject("PidK")) ;
  fProb[4] = dynamic_cast<TH1F*>(fOutput->FindObject("PidP")) ;

  Float_t n=1.292; //mean freon ref idx 
  TF1 *pPi=new TF1("RiPiTheo","acos(sqrt(x*x+[0]*[0])/(x*[1]))",1.2,7); pPi->SetLineWidth(1); pPi->SetParameter(1,n); 
  AliPID ppp;                 pPi->SetLineColor(kRed);   pPi->SetParameter(0,AliPID::ParticleMass(AliPID::kPion));    //mass
  TF1 *pK=(TF1*)pPi->Clone(); pK ->SetLineColor(kGreen); pK ->SetParameter(0,AliPID::ParticleMass(AliPID::kKaon)); 
  TF1 *pP=(TF1*)pPi->Clone(); pP ->SetLineColor(kBlue);  pP ->SetParameter(0,AliPID::ParticleMass(AliPID::kProton)); 

  TCanvas *pC=new TCanvas("c1","ESD QA");pC->SetFillColor(10); pC->SetHighLightColor(10); pC->Divide(3,2);
  pC->cd(1); fCkovP->Draw(); pPi->Draw("same"); pK->Draw("same"); pP->Draw("same");   pC->cd(2); fMipXY->Draw();   pC->cd(3); fProb[0]->Draw(); fProb[1]->Draw("same"); 
  pC->cd(4); fSigP ->Draw();                                                          pC->cd(5); fDifXY->Draw();   pC->cd(6); fProb[2]->Draw(); fProb[3]->Draw("same"); fProb[4]->Draw("same"); 
  
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void loc()
{
  TChain* pChain =new TChain("esdTree");
  pChain->Add("AliESDs.root");

  pChain->Process("EsdQa.C+");	
}

void caf()
{
  gBenchmark->Start("PRooF exec");
  TChain* pChain =new TChain("esdTree");
  
  ifstream list; list.open("list.txt");

  TString file;
  while(list.good()) {
    list>>file;
    if (!file.Contains("root")) continue; //it's wrong file name
    pChain->Add(file.Data());
  }
  list.close();
  
  pChain->GetListOfFiles()->Print();
  
  TVirtualProof *pProof=TProof::Open("kir@lxb6046.cern.ch");	
  pProof->UploadPackage("ESD.par");
  pProof->EnablePackage("ESD");
  
  pChain->MakeTDSet()->Process("EsdQa.C+");
  
  gBenchmark->Show("PRooF exec");
}

