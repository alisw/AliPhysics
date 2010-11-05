/*************************************************************************
* Macro correctSPDdNdEta                                                 *
* To read data and correction histograms produced with                   *
* AliAnalysisTaskSPDdNdEta and apply corrections to data                 *    
*                                                                        *
* Author:  M. Nicassio (INFN Bari)                                       *
* Contact: Maria.Nicassio@ba.infn.it, Domenico.Elia@ba.infn.it           *
**************************************************************************/
/*
To do 
1. test prim
2. check errors on corrections
3. trigger efficiency correction 
*/


#include "Riostream.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TList.h" 
#include "TTree.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TCanvas.h"

void correctSPDdNdEta (Char_t* fileRaw, Char_t* fileRawBkgCorr, Char_t* fileAlphaCorr, Char_t* fileMCBkgCorr, Char_t* fileMC, 
                       Bool_t kPrimariesTest=kFALSE,
                       Int_t binVtxStart=5, Int_t binVtxStop=15,
                       Double_t cutCorr=0., Double_t cutBkg=1., Double_t cutBkgMC=1., Float_t scaleBkg=0.255, Float_t scaleBkgMC=.255)   { 

//lim vtx: 5-15 ->[-10cm,10cm[; 0-20 ->[-20cm,20cm[  

  const Int_t nBinMultCorr = 200; 
  Float_t nMaxMult = 20000.;
  Double_t binLimMultCorr[nBinMultCorr+1];
  binLimMultCorr[0] = 0.;
  binLimMultCorr[1] = 1.;
  for (Int_t i = 2; i<=nBinMultCorr;++i) {
    binLimMultCorr[i] = (i-1)*nMaxMult/nBinMultCorr;
  }

  const Int_t nBinEtaCorr = 60; 
  Float_t etaMax = 3.;
  Float_t etaMin = -3.;
  Double_t *binLimEtaCorr = new Double_t[nBinEtaCorr+1];
  for (Int_t i = 0; i<nBinEtaCorr+1; ++i) {
    binLimEtaCorr[i] = (Double_t) etaMin+i*(etaMax*2.)/nBinEtaCorr;
  }
  
  Float_t nBinsPerPseudorapidityUnit = nBinEtaCorr/(binLimEtaCorr[nBinEtaCorr]-binLimEtaCorr[0]);
  cout<<"Number of bins per pseudorapidity unit"<<nBinsPerPseudorapidityUnit<<endl;

  gStyle->SetOptLogy(kFALSE);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetFrameFillColor(0); 
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetHistLineWidth(2);
  gStyle->SetLabelSize(0.05, "x");
  gStyle->SetLabelSize(0.05, "y");
  gStyle->SetTitleSize(0.05, "x");
  gStyle->SetTitleSize(0.05, "y");
  gStyle->SetTitleOffset(1.1, "x");
  gStyle->SetTitleOffset(0.9, "y"); 
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPalette(1,0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);

  //Input files
  TFile *f_MC = new TFile(fileMC);
  TFile *f_acorr = new TFile(fileAlphaCorr);
  TFile *f_mcbkgcorr = new TFile(fileMCBkgCorr);
  TFile *f_rawbkg = new TFile(fileRawBkgCorr);
  TFile *f_raw = new TFile(fileRaw);

  TList *l_MC = (TList*)f_MC->Get("cOutput");
  TList *l_acorr = (TList*)f_acorr->Get("cOutput");
  TList *l_mcbkgcorr = (TList*)f_mcbkgcorr->Get("cOutput");
  TList *l_rawbkg = (TList*)f_rawbkg->Get("cOutput");
  TList *l_raw = (TList*)f_raw->Get("cOutput");

  //Histogram to be corrected at tracklet level
  TH2F *hSPDEta_2Draw = new TH2F(*((TH2F*)l_raw->FindObject("fHistSPDRAWEtavsZ")));
  hSPDEta_2Draw->SetNameTitle("SPDEta_2Draw","Reconstructed tracklets");

  //Corrections at track level
  // Combinatorial background: beta correction from data
  TH2F *hCombBkg = (TH2F*)l_rawbkg->FindObject("fHistSPDRAWEtavsZ");
  hCombBkg->Scale(scaleBkg);

  TH2F *hCombBkgCorrData = new TH2F("backgroundCorrDATA","Combinatorial background correction",nBinEtaCorr,binLimEtaCorr,20,-20.,20.);
  hCombBkgCorrData->GetXaxis()->SetTitle("#eta");
  hCombBkgCorrData->GetYaxis()->SetTitle("z_{SPDvtx} [cm]");
  hCombBkgCorrData->Sumw2();
  hCombBkgCorrData->Divide(hCombBkg,hSPDEta_2Draw,1.,1.);

  // Combinatorial background: beta correction from MC
  TH2F *hCombBkgMC = (TH2F*)l_mcbkgcorr->FindObject("fHistSPDRAWEtavsZ"); 
  hCombBkgMC->Scale(scaleBkgMC);
  TH2F *hCombBkgMCDen = (TH2F*)l_acorr->FindObject("fHistSPDRAWEtavsZ");

  TH2F *hCombBkgCorrMC = new TH2F("backgroundCorrMC","Combinatorial background correction",nBinEtaCorr,binLimEtaCorr,20,-20.,20.);
  hCombBkgCorrMC->GetXaxis()->SetTitle("#eta");
  hCombBkgCorrMC->GetYaxis()->SetTitle("z_{SPDvtx} [cm]");
  hCombBkgCorrMC->Sumw2();
  hCombBkgCorrMC->Divide(hCombBkgMC,hCombBkgMCDen,1.,1.);
  
  // 1-beta in MC
  TH2F *hCombBkgCorr2MC = new TH2F("backgroundCorr2MC","Combinatorial background correction",nBinEtaCorr,binLimEtaCorr,20,-20.,20.);  
  for (Int_t ieta=0; ieta<nBinEtaCorr; ieta++) 
    for (Int_t jz=0; jz<20; jz++) 
      hCombBkgCorr2MC->SetBinContent(ieta+1,jz+1,1.);
  hCombBkgCorr2MC->Add(hCombBkgCorrMC,-1.);

//  new TCanvas();
//  hCombBkgCorr2MC->Draw();

  // Errors on beta
  //..data
  for (Int_t ieta=0; ieta<nBinEtaCorr; ieta++) {
    for (Int_t jz=0; jz<20; jz++) {
      hCombBkgCorrData->SetBinError(ieta+1,jz+1,scaleBkg*TMath::Sqrt(hCombBkg->GetBinContent(ieta+1,jz+1))/hSPDEta_2Draw->GetBinContent(ieta+1,jz+1));
    }
  }

  //..MC--> no, computed in alpha


  // Alpha correction: MC only
  // alpha here is 1/alpha 
  TH2F *hPrimaries_evtTrVtx = new TH2F(*((TH2F*) l_acorr->FindObject("fHistNonDetectableCorrNum"))); // all prim in sel events
  
  TH2F *hAlphaCorr = new TH2F("AlphaCorr","",nBinEtaCorr,binLimEtaCorr,20,-20.,20.);
  hAlphaCorr->GetXaxis()->SetTitle("#eta");
  hAlphaCorr->GetYaxis()->SetTitle("z_{vtx} [cm]");
//  hAlphaCorr->Divide(hPrimaries_evtTrVtx,hCombBkgMCDen,1.,1.);
//  hAlphaCorr->Divide(hCombBkgCorr2MC);
  hAlphaCorr->Multiply(hCombBkgCorr2MC,hCombBkgMCDen);
  hAlphaCorr->Divide(hPrimaries_evtTrVtx);
  
    // Errors on alpha

  for (Int_t ieta=0; ieta<nBinEtaCorr; ieta++) {
    for (Int_t jz=0; jz<20; jz++) {
      hAlphaCorr->SetBinError(ieta+1,jz+1,TMath::Sqrt((hCombBkgCorr2MC->GetBinContent(ieta+1,jz+1))*(hCombBkgMCDen->GetBinContent(ieta+1,jz+1)))/hPrimaries_evtTrVtx->GetBinContent(ieta+1,jz+1));
    }
  }

  

  // Number of events and normalization histogram
  TH1D *hSPDvertex = (TH1D*)l_raw->FindObject("fHistSPDvtxAnalysis");
  Double_t nEvtsRec=hSPDvertex->Integral(binVtxStart+1,binVtxStop);
  cout<<"Number of reconstructed events in the selected vertex range: "<<nEvtsRec <<endl;

  //Cutting corrections at the edges of the acceptance
  for (Int_t jeta=0; jeta<nBinEtaCorr; jeta++) {
    for (Int_t kvtx=0; kvtx<20; kvtx++) {
      if (hAlphaCorr->GetBinContent(jeta+1,kvtx+1)<cutCorr) {
        hAlphaCorr->SetBinContent(jeta+1,kvtx+1,0.);
        hAlphaCorr->SetBinError(jeta+1,kvtx+1,0.);
      }
      if (hCombBkgCorrData->GetBinContent(jeta+1,kvtx+1)>cutBkg) { 
        hCombBkgCorrData->SetBinContent(jeta+1,kvtx+1,0.);
        hCombBkgCorrData->SetBinError(jeta+1,kvtx+1,0.);
      }
      if (hCombBkgCorrMC->GetBinContent(jeta+1,kvtx+1)>cutBkgMC) {
        hCombBkgCorrMC->SetBinContent(jeta+1,kvtx+1,0.);
        hCombBkgCorrMC->SetBinError(jeta+1,kvtx+1,0.);
      }
    }
  }

  new TCanvas();
  hAlphaCorr->Draw();
  new TCanvas();
  hCombBkgCorrMC->Draw();
  new TCanvas();
  hCombBkgCorrData->Draw();


  // MC distribution --> ok if running on samples for centrality bins

  TH2F *hVertexMC2D = (TH2F*)l_MC->FindObject("fHistSelEvts");  // pay att here!! change if comparing to true MC centrality 
  TH1D *hMCvertex = new TH1D(*(hVertexMC2D->ProjectionY("MCvertex",0,-1,"e")));  //MC vertex distrib all events
  TH2F *Eta = (TH2F*)l_MC->FindObject("fHistNonDetectableCorrNum"); // test comparing with MC in same subset (centrality selected)

//  TH2F *Eta = (TH2F*)l_MC->FindObject("fHistAllPrimaries"); // to compare with MC selection!!
  TH1D *hdNdEta = new TH1D("dNdEta","Pseudorapidity ",nBinEtaCorr,binLimEtaCorr);
  hdNdEta->Sumw2();
  hdNdEta = Eta->ProjectionX("dNdEta",0,-1,"e");
  Double_t nEvtsTot=hMCvertex->Integral(0,hMCvertex->GetNbinsX()+1);
  Double_t nEvtsTot10= hMCvertex->Integral(binVtxStart+1,binVtxStop);

  cout<<"Number of generated events in +-20cm: "<<nEvtsTot <<endl;
  cout<<"Number of generated events in the selected vertex range: "<<nEvtsTot10 <<endl;

  hdNdEta->Scale(nBinsPerPseudorapidityUnit/nEvtsTot);
  hdNdEta->GetXaxis()->SetTitle("#eta");
  hdNdEta->GetYaxis()->SetTitle("dN_{ch}/d#eta");
  new TCanvas();
  hdNdEta->Draw();

  //Corrected distributions
  TH2F *hSPDEta_2D_bkgCorr =  new TH2F("SPDEta_2D_bkgCorr", "",nBinEtaCorr,binLimEtaCorr,20,-20.,20.);
  hSPDEta_2D_bkgCorr->Sumw2();

  TH2F *hSPDEta_2D_fullyCorr =  new TH2F("SPDEta_2D_fullyCorr", "",nBinEtaCorr,binLimEtaCorr,20,-20.,20.);
  hSPDEta_2D_fullyCorr->Sumw2();

  TH1F *hnorm_fullyCorr =  new TH1F("norm_fullyCorr", "",nBinEtaCorr,binLimEtaCorr);
  hnorm_fullyCorr->Sumw2();


  //dN/deta
  TH1F *hdNdEta_raw = new TH1F("dNdEta_raw","",nBinEtaCorr,binLimEtaCorr);
  hdNdEta_raw->Sumw2();
  hdNdEta_raw->GetXaxis()->SetTitle("#eta");
  hdNdEta_raw->GetYaxis()->SetTitle("dN_{ch}/#eta");
  hdNdEta_raw->Add(hSPDEta_2Draw->ProjectionX("SPDEta_2Draw_x",binVtxStart+1,binVtxStop,"e"));
  hdNdEta_raw->Scale(nBinsPerPseudorapidityUnit/nEvtsRec);
  new TCanvas();
  hdNdEta_raw->Draw(); 

  TH1F *hdNdEta_bkgCorr = new TH1F("dNdEta_bkgCorr","",nBinEtaCorr,binLimEtaCorr);
  hdNdEta_bkgCorr->Sumw2();
  TH1F *hdNdEta_fullyCorr = new TH1F("dNdEta_fullyCorr","",nBinEtaCorr,binLimEtaCorr);
  hdNdEta_fullyCorr->Sumw2();

  //Apply corrections at
  //...track level
  hSPDEta_2D_bkgCorr->Add(hSPDEta_2Draw,-1);
  hSPDEta_2D_bkgCorr->Multiply(hCombBkgCorrData);
  hSPDEta_2D_bkgCorr->Add(hSPDEta_2Draw); 
  hSPDEta_2D_fullyCorr->Add(hSPDEta_2D_bkgCorr);
  hSPDEta_2D_fullyCorr->Divide(hAlphaCorr); 
  new TCanvas();
  hSPDEta_2D_bkgCorr->Draw();
  new TCanvas();
  hSPDEta_2D_fullyCorr->Draw();

  //Filling normalization histogram
  for (Int_t ivtx=binVtxStart; ivtx<binVtxStop; ivtx++) {
    Double_t nEvtsivtx = hSPDvertex->GetBinContent(ivtx+1);
    for (Int_t jeta=0; jeta<nBinEtaCorr; jeta++) {
      if (hAlphaCorr->GetBinContent(jeta+1,ivtx+1)!=0.) { 
        hnorm_fullyCorr->Fill(hAlphaCorr->GetXaxis()->GetBinCenter(jeta+1),nEvtsivtx);
      }
    }
  }
  
  for (Int_t jeta=0; jeta<nBinEtaCorr; jeta++) {
    hnorm_fullyCorr->SetBinError(jeta+1,0);
  }
  new TCanvas();
  hnorm_fullyCorr->Draw();

  hdNdEta_bkgCorr->Add(hSPDEta_2D_bkgCorr->ProjectionX("SPDEta_2D_bkgCorr_x",binVtxStart+1,binVtxStop,"e"));
  hdNdEta_bkgCorr->Scale(1./nEvtsRec); 
  hdNdEta_fullyCorr->Divide(hSPDEta_2D_fullyCorr->ProjectionX("SPDEta_2D_fullyCorr_x",binVtxStart+1,binVtxStop,"e"),hnorm_fullyCorr);

  hdNdEta_bkgCorr->Scale(nBinsPerPseudorapidityUnit);
  hdNdEta_fullyCorr->Scale(nBinsPerPseudorapidityUnit);

  new TCanvas();
  hdNdEta_bkgCorr->Draw();
  new TCanvas();
  hdNdEta_fullyCorr->Draw();
  hdNdEta->Draw("same");
  // Create mask for MC prim in SPD acceptance
  TH1D *hdNdEta_test = new TH1D("dNdEta_test","Pseudorapidity ",nBinEtaCorr,binLimEtaCorr);  
  TH2F *hMask = new TH2F("Mask","",nBinEtaCorr,binLimEtaCorr,20,-20,20);
  hMask->Sumw2();
  for (Int_t jeta=0; jeta<nBinEtaCorr; jeta++) {
    for (Int_t kvtx=0; kvtx<20; kvtx++) {
       if (hSPDEta_2D_fullyCorr->GetBinContent(jeta+1,kvtx+1)!=0) hMask->SetBinContent(jeta+1,kvtx+1,1.);

    }
  }
  hMask->Multiply(Eta);
  hdNdEta_test->Add(hMask->ProjectionX("Mask_x",binVtxStart+1,binVtxStop,"e"));
  hdNdEta_test->Scale(1./nEvtsTot10);
  hdNdEta_test->Scale(nBinsPerPseudorapidityUnit);

  TH1F* hRatiotest=new TH1F("ratiotest","",nBinEtaCorr,binLimEtaCorr);
  hRatiotest->GetXaxis()->SetTitle("#eta");
  hRatiotest->GetYaxis()->SetTitle("Generated/Corrected");

  for (Int_t i=0;i<nBinEtaCorr;i++) hdNdEta_test->SetBinError(i+1,0.);
  hRatiotest->Divide(hdNdEta_test,hdNdEta_fullyCorr,1.,1.);
  new TCanvas();
  hRatiotest->Draw("p,histo");


  // Generated/corrected ratios
  TH2F* hRatiodNdEta2D=new TH2F("ratiodNdEta2D","",nBinEtaCorr,binLimEtaCorr,20,-20,20);
  hRatiodNdEta2D->Divide(hSPDEta_2D_fullyCorr,Eta,1.,1.);
  new TCanvas();
  hRatiodNdEta2D->Draw();

  TH1F* hRatiodNdEta=new TH1F("ratiodNdEta","",nBinEtaCorr,binLimEtaCorr);  // Add ratio for two corrected dN/dEta 
  hRatiodNdEta->GetXaxis()->SetTitle("#eta");
  hRatiodNdEta->GetYaxis()->SetTitle("Generated/Corrected");

  for (Int_t i=0;i<nBinEtaCorr;i++) hdNdEta->SetBinError(i+1,0.);
  hRatiodNdEta->Divide(hdNdEta,hdNdEta_fullyCorr,1.,1.); 
  new TCanvas();
  hRatiodNdEta->Draw("p,histo"); 
  // Draw dN/dEta
  TCanvas *ccorrected = new TCanvas("c3","Corrected distributions");
  hdNdEta->SetLineWidth(3);
  hdNdEta->SetMinimum(0.);
  hdNdEta->SetMaximum(4000);
  hdNdEta->Draw("histo");

  hdNdEta_raw->SetLineColor(kGreen);
  hdNdEta_raw->SetLineWidth(3);
  hdNdEta_raw->Draw("histo,same");

  hdNdEta_bkgCorr->SetMarkerStyle(21);
  hdNdEta_bkgCorr->SetMarkerColor(kBlue);

  hdNdEta_bkgCorr->Draw("same,p,histo");

  hdNdEta_fullyCorr->SetMarkerStyle(20);
  hdNdEta_fullyCorr->SetMarkerColor(kRed);
  hdNdEta_fullyCorr->Draw("p,same,histo");

  TLegend *leg1 = new TLegend(0.2,0.7,0.7,0.9,NULL,"brNDC");
  leg1->SetFillColor(0);
  leg1->SetTextFont(62);
  leg1->SetTextSize(0.0304569);
  TLegendEntry *entry1=leg1->AddEntry(hdNdEta_raw,"Reconstructed","l");
                entry1=leg1->AddEntry(hdNdEta,"Generated","l");
                entry1=leg1->AddEntry(hdNdEta_bkgCorr,"Bkg corrected","p");
                entry1=leg1->AddEntry(hdNdEta_fullyCorr,"Corrected (triggered events with vertex)","p");

  leg1->Draw();
/*  new TCanvas();
  // plot the relative stat error for this correction
  TH2F* hStatErrTrackToPart = new TH2F("staterrperctracktopart","",nBinEtaCorr,binLimEtaCorr,20,-20.,20.);
  hStatErrTrackToPart->GetXaxis()->SetTitle("#eta");
  hStatErrTrackToPart->GetYaxis()->SetTitle("z_{vtx} [cm]");
  hStatErrTrackToPart->GetZaxis()->SetTitle("Statistical error (%)");
  for (Int_t kvtx=0; kvtx<20; kvtx++) {
    for (Int_t jeta=0; jeta<nBinEtaCorr; jeta++) {
      if (hTrackToParticleCorr->GetBinContent(jeta+1,kvtx+1)) 
      hStatErrTrackToPart->SetBinContent(jeta+1,kvtx+1,(hTrackToParticleCorr->GetBinError(jeta+1,kvtx+1))/(hTrackToParticleCorr->GetBinContent(jeta+1,kvtx+1))*100);
      else hStatErrTrackToPart->SetBinContent(jeta+1,kvtx+1,0);
    }
  }
  new TCanvas();
  hStatErrTrackToPart->Draw();
  new TCanvas();
  hRatiodNdEta_trackToParticle->Draw();
new TCanvas();
hTrackToParticleCorr->Draw();
  // Write histos
*/
  TFile *fout= new TFile("SPDdNdEta.root","RECREATE"); 


  // save more histos
  hCombBkgCorrData->Write();
  hCombBkgCorrMC->Write();
  hCombBkgCorr2MC->Write();
  hAlphaCorr->Write();
  fout->Write();
  fout->Close();

}
