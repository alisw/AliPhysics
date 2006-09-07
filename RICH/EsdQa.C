// The class definition in esdAna.h has been generated automatically
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
// Root > T->Process("esdAna.C")
// Root > T->Process("esdAna.C","some options")
// Root > T->Process("esdAna.C+")
//

#include "EsdQa.h"
#include <AliESD.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TGrid.h>
#include <TAlienCollection.h>
#include <TMath.h>

void EsdQa::Begin(TTree *)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

void EsdQa::SlaveBegin(TTree *tree)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   Init(tree);

   TString option = GetOption();

   // create histograms on each slave server
   fCkovMomH = new TH2F("CkovMomH", "Ckov angel,[rad];P [GeV]", 150, 0,15, 100 , 0, 1); 
   fsigma2   = new TH1F("sigma2","#sigma_{#theta_C}^2",2000,1e9,2e10);
   fX        = new TH1F("X","mipX[cm]",260,0,130);
   fY        = new TH1F("Y","mipY[cm]",252,0,126);
   fdist     = new TH2F("dist","dist mip-track[cm]",150,0,15,100,0,10);

   fProbH[0] = new TH1F("pidE" ,"pid of e"  ,100,0,1);
   fProbH[1] = new TH1F("pidMu","pid of #mu",100,0,1);
   fProbH[2] = new TH1F("pidPi","pid of #pi",100,0,1);
   fProbH[3] = new TH1F("pidK" ,"pid of K"  ,100,0,1);
   fProbH[4] = new TH1F("pidP" ,"pid of p"  ,100,0,1);
}

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
     
     Float_t dx, dy;
     pTrk->GetRICHdxdy(dx,dy);
     Float_t r;
         if(dx<0 && dy<0) r=999.;
        else r = TMath::Sqrt(dx*dx+dy*dy);
     fdist->Fill(pTrk->GetP(),r);
     if(pTrk->GetRICHsignal()<0) continue;
     fCkovMomH->Fill(pTrk->GetP(),pTrk->GetRICHsignal()) ; 
     Double_t pid[5];
     pTrk->GetRICHpid(pid);
     fProbH[0]->Fill(pid[0]);    fProbH[1]->Fill(pid[1]);    fProbH[2]->Fill(pid[2]);    fProbH[3]->Fill(pid[3]);    fProbH[4]->Fill(pid[4]);
     fsigma2->Fill(pTrk->GetRICHchi2());
     Float_t xmip, ymip; 
     pTrk->GetRICHmipXY(xmip,ymip);
     fX->Fill(xmip); fY->Fill(ymip);
    } 
     
  return kTRUE;
}

void EsdQa::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
  
  // Add the histograms to the output on each slave server
  
  fOutput->Add(fCkovMomH) ;
  fOutput->Add(fsigma2) ; 
  fOutput->Add(fX) ;
  fOutput->Add(fY) ;
  fOutput->Add(fdist) ;
  fOutput->Add(fProbH[0]);  fOutput->Add(fProbH[1]);
  fOutput->Add(fProbH[2]);  fOutput->Add(fProbH[3]);
  fOutput->Add(fProbH[4]);
}

void EsdQa::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
  
  fCkovMomH = dynamic_cast<TH2F*>(fOutput->FindObject("CkovMomH")) ;
  fsigma2 = dynamic_cast<TH1F*>(fOutput->FindObject("sigma2")) ; 
  fX = dynamic_cast<TH1F*>(fOutput->FindObject("X")) ;
  fY = dynamic_cast<TH1F*>(fOutput->FindObject("Y")) ;
  fdist = dynamic_cast<TH2F*>(fOutput->FindObject("dist")) ;
  fProbH[0] = dynamic_cast<TH1F*>(fOutput->FindObject("pidE")) ;
  fProbH[1] = dynamic_cast<TH1F*>(fOutput->FindObject("pidMu")) ;
  fProbH[2] = dynamic_cast<TH1F*>(fOutput->FindObject("pidPi")) ;
  fProbH[3] = dynamic_cast<TH1F*>(fOutput->FindObject("pidK")) ;
  fProbH[4] = dynamic_cast<TH1F*>(fOutput->FindObject("pidP")) ;


  TFile * file = TFile::Open("esdAnaM.root", "RECREATE");
  fCkovMomH->Write() ; 
  fX->Write();
  fY->Write();
  fdist->Write();
  for(Int_t i=0;i<5;i++) fProbH[i] -> Write();
  fsigma2->Write();
//  fHistNeutralMul->Write() ; 
  file->Close() ;
  delete file ;
  if (!gROOT->IsBatch()) {
    TCanvas *c1 = new TCanvas("c1","Annalisa's best trial");c1->SetFillColor(10); c1->SetHighLightColor(10);
    fCkovMomH->DrawCopy() ;
   }
}


void run()
{
  TGrid::Connect("alien://",0,0,"t");
  TAlienCollection *pAlienColl = new TAlienCollection("newIIcento.xml") ;

  TChain* pChain = new TChain("esdTree");

  pAlienColl->Reset() ;
  while (pAlienColl->Next()) {
    char esdFile[255] ;
    sprintf(esdFile, "%s", pAlienColl->GetTURL("")) ;
    Printf("Adding alien file %s", esdFile);
    pChain->Add(esdFile) ;
  }

  EsdQa *pSel =new EsdQa;
  pChain->Process(pSel);	
}
