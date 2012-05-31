//////////////////////////////////////////////////////////////////
// Macro to check ESD info provided by the 2 SPD layers         //
// Provides:                                                    //
//  1 canvas with                                               //
//      - distribution of the SPD vertex Z-coord                //
//      - SPD vertex Z-coord and sigma for each event           //
//  1 canvas for tracklets with                                 //
//      - #tracklets                                            //
//      - #tracklets vs #clusters (inner layer)                 //
//      - deltaPhi                                              //
//      - Phi                                                   //
//  1 canvas for clusters (inner layer) not associated with     //
//      - #single clusters                                      // 
//      - Phi                                                   //
//      - Theta                                                 //
//                                                              //
//  Maria.Nicassio@ba.infn.it                                   //
//  Domenico.Elia@ba.infn.it                                    //
//////////////////////////////////////////////////////////////////
                                                                                                                         
#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"                                                                                                               
#include "AliESD.h"                                                                                                               
#include "AliESDEvent.h"
#include "AliESDHeader.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliMultiplicity.h"
#endif

/* $Id$ */
void ShowSPDESDs (Int_t RunStart, Int_t RunStop) {

  Char_t fileName[256];

  Char_t* dir = "/home/elia/alice/pp/data/first";            

  TH1F* hSPDVertex = new TH1F("SPDVertex","",80,-20.,20.);
  TH1F* hSPDVertexPerEvent = new TH1F("SPDVertexPerEvent","",100,0.,100.);

  TH1F* hnTracklets = new TH1F("nTracklets","",200,0.,200.);
  TH2F* hnTracklets_nCl1 = new TH2F("nTracklets_nCl1","",200,0.,200.,200,0.,200.);
  TH1F* hDePhiTracklets = new TH1F("DePhiTracklets","",200,-0.2,0.2);
  TH1F* hPhiTracklets = new TH1F("PhiTracklets","",600,0.,2*TMath::Pi());
  TH1F* hThetaTracklets = new TH1F("ThetaTracklets","",300,0.,TMath::Pi());
  TH1F* hnSingleClustersLay1 = new TH1F("nSingleClustersLay1","",200,0.,200.);
  TH1F* hPhiSingleClustersLay1 = new TH1F("PhiSingleClustersLay1","",600,0.,2*TMath::Pi());
  TH1F* hThetaSingleClustersLay1 = new TH1F("ThetaSingleClustersLay1","",600,0.,TMath::Pi());

  // Loop over runs ...
  for (Int_t run=RunStart; run<RunStop+1; run++) {

    sprintf(fileName,"%s/AliESDsRAW_%04d.root",dir,run);

    // Open input file and get the TTree
    TFile inFile(fileName, "READ");
  
    TTree *esdTree = (TTree*)inFile.Get("esdTree");
    AliESDEvent *esd = new AliESDEvent();
    esd->ReadFromTree(esdTree);
  
    // Loop over events
    Int_t nESDEvents = esdTree->GetEntries();
    for (Int_t iEv = 0; iEv < nESDEvents; iEv++) {
//    cout << "Event: " << iEv+1 << "/" << nESDEvents << endl;

      // Read events
      esdTree->GetEvent(iEv);      

      // Get the ESD vertex
      const AliESDVertex* vtxESD = esd->GetVertex();

      Double_t ESDvtx[3];
      Double_t sigmaESDvtx[3];
      vtxESD->GetXYZ(ESDvtx);
      vtxESD->GetSigmaXYZ(sigmaESDvtx);
//    cout<<"SPD vtx: x="<<ESDvtx[0]<<" y="<<ESDvtx[1]<<" z="<<ESDvtx[2]<<endl;
//    cout<<"sigma SPD vtx: x="<<sigmaESDvtx[0]<<" y="<<sigmaESDvtx[1]<<" z="<<sigmaESDvtx[2]<<endl;
      hSPDVertex->Fill(ESDvtx[2]);          
      if (ESDvtx[2]!=0.) hSPDVertexPerEvent->SetBinError(iEv,sigmaESDvtx[2]);         
      hSPDVertexPerEvent->Fill(iEv,ESDvtx[2]);                                         

      // Get the SPD reconstructed multiplicity
      const AliMultiplicity* SPDMult = esd->GetMultiplicity();

      Int_t nTracklets = SPDMult->GetNumberOfTracklets();
      // Loop over tracklets
      for (Int_t itr=0; itr<nTracklets; ++itr) {
        Float_t thetaTr = SPDMult->GetTheta(itr);
        Float_t phiTr = SPDMult->GetPhi(itr);
        Float_t dePhiTr = SPDMult->GetDeltaPhi(itr);
        hPhiTracklets->Fill(phiTr);
        hDePhiTracklets->Fill(dePhiTr);
        hThetaTracklets->Fill(thetaTr);
      }

      Int_t nSingleCl1 = SPDMult->GetNumberOfSingleClusters();
      Int_t nCl1 = nSingleCl1 + nTracklets;
      // Loop over unassociated clusters
      for (Int_t icl=0; icl<nSingleCl1; ++icl) {
        Float_t phiSing = SPDMult->GetPhiSingle(icl);
        Float_t thetaSing = SPDMult->GetThetaSingle(icl);
        hPhiSingleClustersLay1->Fill(phiSing);
        hThetaSingleClustersLay1->Fill(thetaSing);
      }   
      hnTracklets_nCl1->Fill(nCl1,nTracklets);
      hnTracklets->Fill(nTracklets);
      hnSingleClustersLay1->Fill(nSingleCl1);

    } // end loop over events
  } // end loop over runs
  
  TFile* fout = new TFile("out_ShowSPDESDs.root","RECREATE");
  
  TCanvas *cVertex = new TCanvas("cVertex","SPD vertex");
  cVertex->Divide(1,2);
  cVertex->cd(1);
  hSPDVertex->GetYaxis()->SetTitle("Entries");
  hSPDVertex->GetXaxis()->SetTitle("Reconstructed vertex [cm]");
  hSPDVertex->Draw();
  cVertex->cd(2);
  hSPDVertexPerEvent->SetMarkerStyle(21);
  hSPDVertexPerEvent->SetMarkerColor(2);
  hSPDVertexPerEvent->SetMarkerSize(0.4);
  hSPDVertexPerEvent->GetXaxis()->SetTitle("Event number");
  hSPDVertexPerEvent->GetYaxis()->SetTitle("Reconstructed vertex [cm]");
  hSPDVertexPerEvent->Draw("p");                                                                                            
  TCanvas *cTracklets = new TCanvas("cTracklets","SPD Tracklets");
  cTracklets->Divide(2,2);
  cTracklets->cd(1);
  hnTracklets->GetXaxis()->SetTitle("Number of tracklets");
  hnTracklets->GetYaxis()->SetTitle("Entries");
  hnTracklets->Draw();
  cTracklets->cd(2);
  hnTracklets_nCl1->GetXaxis()->SetTitle("# clusters (inner layer)");
  hnTracklets_nCl1->GetYaxis()->SetTitle("# tracklets");
  hnTracklets_nCl1->Draw();
  cTracklets->cd(3);
  hDePhiTracklets->GetXaxis()->SetTitle("#Delta#phi [rad]");
  hDePhiTracklets->GetYaxis()->SetTitle("Entries");
  hDePhiTracklets->Draw(); 
  cTracklets->cd(4);
  hPhiTracklets->GetXaxis()->SetTitle("#phi [rad]");
  hPhiTracklets->GetYaxis()->SetTitle("Entries");
  hPhiTracklets->Draw();

  TCanvas *cSingleClusters = new TCanvas("cSingleClusters","Unassociated clusters");
  cSingleClusters->Divide(2,2);
  cSingleClusters->cd(1);
  hnSingleClustersLay1->GetXaxis()->SetTitle("# clusters (inner layer)");
  hnSingleClustersLay1->GetYaxis()->SetTitle("Entries");
  hnSingleClustersLay1->Draw();
  cSingleClusters->cd(2);
  hPhiSingleClustersLay1->GetXaxis()->SetTitle("#phi [rad]");
  hPhiSingleClustersLay1->GetYaxis()->SetTitle("Entries");
  hPhiSingleClustersLay1->Draw();
  cSingleClusters->cd(3);
  hThetaSingleClustersLay1->GetXaxis()->SetTitle("#theta [rad]");
  hThetaSingleClustersLay1->GetYaxis()->SetTitle("Entries");
  hThetaSingleClustersLay1->Draw();

  hSPDVertex->Write();
  hSPDVertexPerEvent->Write();
  hnTracklets_nCl1->Write();
  hnTracklets->Write();

  hDePhiTracklets->Write();
  hPhiTracklets->Write();
  hThetaTracklets->Write();
  hnSingleClustersLay1->Write(); 
  hPhiSingleClustersLay1->Write(); 
  hThetaSingleClustersLay1->Write();

  fout->Close(); 
  return;

}
