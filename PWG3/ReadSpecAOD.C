#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include "TObject.h"
#include "TStyle.h"
#include "TROOT.h"

#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODEventInfo.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODCluster.h"
#include "AliAODDimuon.h"

#endif

void ReadSpecAOD(const char *fileName = "AliMuonAOD.root") {

  gSystem->Load("libTree.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("STEERBase/libSTEERBase");
  gSystem->Load("AOD/libAOD");
  gROOT->LoadMacro("AliAODEventInfo.cxx++");
  gROOT->LoadMacro("AliAODDimuon.cxx++");
   

  gStyle->SetOptStat(111111);
  gStyle->SetFrameFillColor(10);
  gStyle->SetCanvasColor(10);
  gStyle->SetOptStat(0);
  TH1F *mass=new TH1F("mass","Invariant mass",100,0.,20.);
  TH1F *rap=new TH1F("rap","Rapidity",100,-5.,0.);
  TH1F *cost=new TH1F("cost","Cost_CS",100,-1.,1.);
  TH1F *pt=new TH1F("pt","Pt",100,0.,50.);
  TH1F *ptmuon=new TH1F("ptmuon","single muon Pt",100,0.,50.);

  // open input file and get the TTree
  TFile inFile(fileName, "READ");
  TTree *aodTree = (TTree*)inFile.Get("AOD");

  AliAODEvent *aod = (AliAODEvent*)aodTree->GetUserInfo()->FindObject("AliAODEvent");

  TClonesArray *Dimuons;
  TClonesArray *tracks;
  TClonesArray *vertices;
  AliAODEventInfo *MuonInfos;
  
  aodTree->SetBranchAddress("Dimuons",&Dimuons);
  aodTree->SetBranchAddress("tracks",&tracks);
  aodTree->SetBranchAddress("vertices",&vertices);
  aodTree->SetBranchAddress("MuonInfos",&MuonInfos);

  // loop over events
  Int_t nEvents = aodTree->GetEntries();
  for (Int_t nEv = 0; nEv < nEvents; nEv++) {
    cout << "Event: " << nEv+1 << "/" << nEvents << endl;
    aodTree->GetEntry(nEv);
    // loop over tracks
    Int_t nTracks = tracks->GetEntries();
    for (Int_t nTr = 0; nTr < nTracks; nTr++) {
      AliAODTrack *tr = (AliAODTrack *)tracks->At(nTr);
      ptmuon->Fill(tr->Pt());
      // print track info
      cout << nTr+1 << "/" << nTracks << ": track pt: " << tr->Pt();
      if (tr->GetProdVertex()) {
	cout << ", vertex z of this track: " << tr->GetProdVertex()->GetZ();
      }
      cout << endl;
    }
    // loop over dimuons
    
    Int_t nDimuons = Dimuons->GetEntries();
    cout << nDimuons << " dimuon(s)" << endl;
    for(Int_t nDi=0; nDi < nDimuons; nDi++){
      AliAODDimuon *di=(AliAODDimuon*)Dimuons->At(nDi);
       if((MuonInfos->MUON_Unlike_HPt_L0())){
         mass->Fill(di->M());
         pt->Fill(di->Pt());
         rap->Fill(di->Y());
         cost->Fill(di->CostCS());
         cout << "Dimuon: " << nDi << " q: " << di->Charge() 
	      << " m: " << di->M() << " Y: " << di->Y() << " Pt: " << di->Pt()<< " CostCS: " << di->CostCS() << endl  ;
       }
    }
//     // loop over vertices
//     Int_t nVtxs = vertices->GetEntries();
//     for (Int_t nVtx = 0; nVtx < nVtxs; nVtx++) {
//       
//       // print track info
//       cout << nVtx+1 << "/" << nVtxs << ": vertex z position: " <<vertices->At(nVtx)->GetZ() << endl;
//     }
  }
  inFile.Close();
  TCanvas *c1=new TCanvas();
  c1->Show();
  c1->Divide(2,2);
  c1->ForceUpdate();
  c1->cd(1);
  mass->DrawClone();
  c1->cd(2);
  rap->DrawClone();
  c1->cd(3);
  pt->DrawClone();
  c1->cd(4);
  cost->DrawClone();
  
  TCanvas *c2 = new TCanvas();
  c2->cd();
  ptmuon->DrawClone();
  return;
}
