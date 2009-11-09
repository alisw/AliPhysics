/*************************************************************************
* Macro correctSPDdNdEta                                                 *
* To read data and correction histograms produced with                   *
* AliAnalysisTaskSPDdNdEta and apply corrections to data                 *    
*                                                                        *
* Author:  M. Nicassio (INFN Bari)                                       *
* Contact: Maria.Nicassio@ba.infn.it, Domenico.Elia@ba.infn.it           *
**************************************************************************/

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

void correctSPDdNdEtapp (Char_t* fileRaw, Char_t* fileCorr, Char_t* fileMC, Char_t* fileAccPb,
                       Bool_t kPrimariesTest, Bool_t kAccPb,
                       Int_t binVtxStart, Int_t binVtxStop,
                       Double_t cutAccPb, Double_t cutAcc, Double_t cutBkg, 
                       Double_t cutAlgEff, Double_t cutDec, Double_t cutVtx, Double_t cutTrig)  { 

//lim vtx: 5-15 ->[-10cm,10cm[; 0-20 ->[-20cm,20cm[  

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
  TFile *f_accPb = new TFile(fileAccPb);
  TFile *f_corr = new TFile(fileCorr);
  TFile *f_MC = new TFile(fileMC);
  TFile *f_raw = new TFile(fileRaw);

  TList *l_corr = (TList*)f_corr->Get("cOutput");
  TList *l_MC = (TList*)f_MC->Get("cOutput");
  TList *l_raw = (TList*)f_raw->Get("cOutput");

  //Corrections at track level
  //Background
  TH2F *hBackgroundCorr = new TH2F("backgroundCorr","Background correction",60,-3.,3.,20,-20.,20.);
  hBackgroundCorr->Sumw2();
  TH2F *hBackgroundCorrNum = (TH2F*)l_corr->At(32);

  Int_t corrHistoName;
  if (kPrimariesTest) {
    corrHistoName = 32;
  } else {
    corrHistoName = 33;
  }

  TH2F *hBackgroundCorrDen = (TH2F*)l_corr->At(corrHistoName);
  hBackgroundCorr->Divide(hBackgroundCorrNum,hBackgroundCorrDen,1.,1.); 
  
  for (Int_t ieta=0; ieta<60; ieta++) {
    for (Int_t jz=0; jz<20; jz++) {
      Float_t den = hBackgroundCorrNum->GetBinContent(ieta+1,jz+1);
      if (den!=0) {
        Float_t num = hBackgroundCorrDen->GetBinContent(ieta+1,jz+1)-den;
        hBackgroundCorr->SetBinError(ieta+1,jz+1,(TMath::Sqrt(num)/den)/TMath::Power(hBackgroundCorr->GetBinContent(ieta+1,jz+1),2));
      } else {
        hBackgroundCorr->SetBinError(ieta+1,jz+1,0);
      }
    }
  }

  //Algorithm efficiency
  TH2F *hTrackletRecEff = new TH2F("efficiencyCorr","Tracklet algorithm + SPD efficiency",60,-3,3,20,-20,20);
  hTrackletRecEff->Sumw2();
  TH2F *hTrackletEffCorr_num = new TH2F(*((TH2F*)l_corr->At(34)));  
  hTrackletEffCorr_num->Rebin2D(2,2,"");
  hTrackletRecEff->Divide(hBackgroundCorrNum,hTrackletEffCorr_num,1.,1.,"b");   // 1/Efficiency=correction_weight 

  //Correction for particles that do not reach the sensitive layer
  TH2F *hDecPartLay2Corr = new TH2F("decPartLay2Corr","Correction for particles not reaching the detector",60,-3.,3.,20,-20.,20.);
  hDecPartLay2Corr->Sumw2();
  TH2F *hPrimaries_evtTrVtx = new TH2F(*((TH2F*) l_corr->At(35)));
  TH2F *hDetectablePrimariesLay2_evtTrVtx = new TH2F(*((TH2F*) l_corr->At(36)));
  hDetectablePrimariesLay2_evtTrVtx->Rebin2D(2,2,""); 
  hDecPartLay2Corr->Divide(hDetectablePrimariesLay2_evtTrVtx,hPrimaries_evtTrVtx,1.,1.,"b");

  //Acceptance
  TH2F *hAccCorr = new TH2F("accCorr_test","SPD acceptance correction",120,-3,3,40,-20,20);
  hAccCorr->Sumw2();
  hAccCorr->Rebin2D();
  TH2F *hAccNum = new TH2F(*((TH2F*)l_corr->At(34)));
  hAccNum->Sumw2();
  TH2F *hAccDen = new TH2F(*((TH2F*)l_corr->At(36)));
  hAccDen->Sumw2();
  hAccNum->Rebin2D();
  hAccDen->Rebin2D();
  hAccCorr->Divide(hAccNum,hAccDen,1.,1.,"b");
  TH2F *hAccCorrRebin = new TH2F(*hAccCorr);
  hAccCorrRebin->Sumw2();
//  hAccCorrRebin->Rebin2D();

  TH2F *hAccCorrPb = (TH2F*) f_accPb->Get("AccTrac");

  //Vertex-trigger correction
  TH2F *hVertexTrigCorrEtaNum = (TH2F*)l_corr->At(37);
  TH2F *hTriggerCorrEta_norm  = (TH2F*)l_corr->At(38);
  TH2F *hVertexCorrEta_norm   = (TH2F*)l_corr->At(35);   
  TH2F *hTrigCorrEtaNumNSD = (TH2F*)l_corr->At(39);

  TH2F *hVertexCorrEta = new TH2F("vertexCorrEta","",60,-3.,3.,20,-20.,20.);
  hVertexCorrEta->Sumw2();
  hVertexCorrEta->Divide(hVertexCorrEta_norm,hTriggerCorrEta_norm,1.,1.,"b");

  TH2F *hTriggerCorrEta = new TH2F("triggerCorrEta","",60,-3.,3.,20,-20.,20.);
  hTriggerCorrEta->Sumw2();
  hTriggerCorrEta->Divide(hTriggerCorrEta_norm,hVertexTrigCorrEtaNum,1.,1.,"b");

  TH2F *hVertexTrigCorrEta = new TH2F("vertexTrigEtaCorr","Correction for MB trigger and vertexer inefficiency",60,-3.,3.,20,-20.,20.);
  hVertexTrigCorrEta->Sumw2();
  hVertexTrigCorrEta->Divide(hVertexCorrEta_norm,hVertexTrigCorrEtaNum,1.,1.,"b");

  TH2F *hVertexTrigCorrEtaNSD = new TH2F("vertexTrigEtaCorrNSD","",60,-3.,3.,20,-20.,20.);
  hVertexTrigCorrEtaNSD->Sumw2();
  hVertexTrigCorrEtaNSD->Divide(hTrigCorrEtaNumNSD,hVertexCorrEta_norm,1.,1.);

  TH2F *hTrigCorrEtaNSD = new TH2F("triggerCorrEtaNSD","",60,-3.,3.,20,-20.,20.);
  hTrigCorrEtaNSD->Divide(hTriggerCorrEta_norm,hTrigCorrEtaNumNSD,1.,1.);
  TH2F *hTrigCorrEtaNSDInv = new TH2F("triggerCorrEtaNSDInv","",60,-3.,3.,20,-20.,20.);
  hTrigCorrEtaNSDInv->Divide(hTrigCorrEtaNumNSD,hTriggerCorrEta_norm,1.,1.);

  TH2F *hTracksinTriggeredEvts_NSD = new TH2F(*(TH2F*)l_corr->At(40)); 
  for (Int_t ieta=0; ieta<60; ieta++) {
    for (Int_t jz=0; jz<20; jz++) {
      Float_t den = hTrigCorrEtaNumNSD->GetBinContent(ieta+1,jz+1);
      if (den!=0) {
        Float_t acc=hTracksinTriggeredEvts_NSD->GetBinContent(ieta+1,jz+1)/hTrigCorrEtaNumNSD->GetBinContent(ieta+1,jz+1);
        Float_t binerr = TMath::Sqrt(acc*(1-acc)/den);
        Float_t nonbinerr = TMath::Sqrt(hTriggerCorrEta_norm->GetBinContent(ieta+1,jz+1)-hTracksinTriggeredEvts_NSD->GetBinContent(ieta+1,jz+1))/den;
        hTrigCorrEtaNSD->SetBinError(ieta+1,jz+1,TMath::Sqrt(TMath::Power(binerr,2)+TMath::Power(nonbinerr,2)));
      } else {
        hTrigCorrEtaNSD->SetBinError(ieta+1,jz+1,0);
      }
    }
  }

  //Corrections at event level
  TH2F *hVertexTrigCorrNum = (TH2F*)l_corr->At(41);
  TH2F *hTriggerCorr_norm  = (TH2F*)l_corr->At(43);
  TH2F *hVertexCorr_norm   = (TH2F*)l_corr->At(42);
  TH2F *hTriggerCorrNumNSD = (TH2F*)l_corr->At(44);

  Int_t nBinMultCorr = 28;
  Double_t binLimMultCorr[nBinMultCorr+1];
  for (Int_t i = 0; i<nBinMultCorr+1; ++i) {
    if (i<21) binLimMultCorr[i] = (Double_t) i;
    else if (i<27) binLimMultCorr[i] = (Double_t) 20 +5*(i-20);
    else if (i==27) binLimMultCorr[i] = 100;
    else if (i==28) binLimMultCorr[i] = 200;
  }


  TH2F *hVertexCorr = new TH2F("vertexCorr","Vertex inefficiency correction",nBinMultCorr,binLimMultCorr,20,-20.,20.); 
  hVertexCorr->Sumw2();
  hVertexCorr->Divide(hVertexCorr_norm,hTriggerCorr_norm,1.,1.,"b");

  TH2F *hTriggerCorr = new TH2F("triggerCorr","MB trigger correction - inelastic events",nBinMultCorr,binLimMultCorr,20,-20.,20.);
  hTriggerCorr->Sumw2();
  hTriggerCorr->Divide(hTriggerCorr_norm,hVertexTrigCorrNum,1.,1.,"b");
 
  TH2F *hTriggerCorrNSD = new TH2F("triggerCorrNSD","MB trigger correction - NSD events",nBinMultCorr,binLimMultCorr,20,-20.,20.);
  hTriggerCorrNSD->Sumw2();
  hTriggerCorrNSD->Divide(hTriggerCorr_norm,hTriggerCorrNumNSD,1.,1.);
  TH2F *hTriggerCorrNSD2 = new TH2F("triggerCorrNSD2","",nBinMultCorr,binLimMultCorr,20,-20.,20.);
  hTriggerCorrNSD2->Sumw2();
  hTriggerCorrNSD2->Divide(hTriggerCorrNumNSD,hTriggerCorr_norm,1.,1.);


  TH2F *hTriggeredEvts_NSD = new TH2F(*(TH2F*)l_corr->At(45)); 
  for (Int_t im=0; im<28; im++) {
    for (Int_t jz=0; jz<20; jz++) {
      Float_t den = hTriggerCorrNumNSD->GetBinContent(im+1,jz+1);
      if (den!=0) {
        Float_t acc=hTriggeredEvts_NSD->GetBinContent(im+1,jz+1)/hTriggerCorrNumNSD->GetBinContent(im+1,jz+1);
        Float_t binerr = TMath::Sqrt(acc*(1-acc)/den);
        Float_t nonbinerr = TMath::Sqrt(hTriggerCorr_norm->GetBinContent(im+1,jz+1)-hTriggeredEvts_NSD->GetBinContent(im+1,jz+1))/den;
        hTriggerCorrNSD->SetBinError(im+1,jz+1,TMath::Sqrt(TMath::Power(binerr,2)+TMath::Power(nonbinerr,2)));
      } else {
        hTriggerCorrNSD->SetBinError(im+1,jz+1,0);
      }
    }
  }

  TH2F *hVertexMC2D = (TH2F*)l_MC->At(41);   
  TH1D *hMCvertex = new TH1D(*(hVertexMC2D->ProjectionY("MCvertex",0,-1,"e")));  //MC vertex distrib all events
//  hMCvertex->Draw(); 

  //Raw distributions

  //...track level
  TH2F *hSPDEta_2Draw_accBin = new TH2F(*((TH2F*)l_raw->At(2)));  
  hSPDEta_2Draw_accBin->SetNameTitle("SPDEta_2Draw_accBin","");
  TH2F *hSPDEta_2Draw = new TH2F(*((TH2F*)l_raw->At(2)));
  hSPDEta_2Draw->SetNameTitle("SPDEta_2Draw","Reconstructed tracklets");

  //...event level
  TH2F *hESDVertexvsNtracklets = new TH2F(*((TH2F*)l_raw->At(0)));
  hESDVertexvsNtracklets->SetNameTitle("ESDVertexvsNtracklets","Triggered events with vertex");
  TH2F *hESDVertexvsNtracklets_trigEvts = new TH2F(*((TH2F*)l_raw->At(1)));
  hESDVertexvsNtracklets_trigEvts->SetNameTitle("ESDVertexvsNtracklets_trigEvts","Triggered events");
  
  TH2F *hESDVertexvsNtracklets_Corr = new TH2F("ESDVertexvsNtracklets_Corr","",nBinMultCorr,binLimMultCorr,20,-20.,20.);
  TH2F *hESDVertexvsNtracklets_CorrNSD = new TH2F("ESDVertexvsNtracklets_CorrNSD","",nBinMultCorr,binLimMultCorr,20,-20.,20.);
 
  TH1D *hSPDvertex = (TH1D*)l_raw->At(19);
  if (kPrimariesTest) {
    hSPDEta_2Draw = new TH2F(*((TH2F*)l_MC->At(32)));
    hSPDEta_2Draw->SetNameTitle("SPDEta_2Draw","");
    hESDVertexvsNtracklets = new TH2F(*((TH2F*)l_MC->At(42)));
    hESDVertexvsNtracklets_trigEvts = new TH2F(*((TH2F*)l_MC->At(43)));

    hSPDvertex = hVertexCorr_norm->ProjectionY("MCvertex_RV",0,-1,"e");  
  }

  if (!kPrimariesTest) hSPDEta_2Draw->Rebin2D(2,2,"");

  //Vertex correction

  //Computing the percentage of triggered events with vertex and tracklet mult!=0 for each rec vertex bin
  //to distribute triggered events with tracklet mult = 0 in vertex bins.
  TH1F* fHistZdistribTrigVtxEvt = new TH1F("fHistZdistribTrigVtxEvt", "",20,-20.,20.);
  fHistZdistribTrigVtxEvt->Sumw2();
  for (Int_t ivtx=0; ivtx<20; ivtx++) {
    Float_t fracZ = (hESDVertexvsNtracklets->ProjectionY("_y",2,hESDVertexvsNtracklets->GetNbinsX()+1,"e"))->GetBinContent(ivtx+1)/
                    Float_t (hESDVertexvsNtracklets->ProjectionY("_y",2,hESDVertexvsNtracklets->GetNbinsX()+1,"e")->Integral(0,hESDVertexvsNtracklets->GetNbinsY()+1));
  
    fHistZdistribTrigVtxEvt->SetBinContent(ivtx+1,fracZ);
    //cout<<" Fractions of rec events: Num "<<hESDVertexvsNtracklets->ProjectionY("_y",2,hESDVertexvsNtracklets->GetNbinsX()+1,"e")->GetBinContent(ivtx+1)<<
    //      " Den  "<<hESDVertexvsNtracklets->ProjectionY("_y",2,hESDVertexvsNtracklets->GetNbinsX()+1,"e")->Integral(0,hESDVertexvsNtracklets->GetNbinsY()+1)<<
    //      " Frac: "<<fracZ<<endl;
  }

  //Computing the MC correction for the previous % 
  TH1F* fHistEvtCorrFactor = new TH1F("fHistEvtCorrFactor", "",20,-20.,20.);
  fHistEvtCorrFactor->Sumw2();
  for (Int_t ivtx=0; ivtx<20; ivtx++) {

    Float_t fracZnum = (hTriggerCorr_norm->GetBinContent(1,ivtx+1))/Float_t(hTriggerCorr_norm->ProjectionY("_y",1,1,"e")->Integral(0,hTriggerCorr_norm->GetNbinsY()+1));

    //cout<<"Vertex bin "<<ivtx+1<<" Content per bin "<<hTriggerCorr_norm->GetBinContent(1,ivtx+1)<<" Integral "<<hTriggerCorr_norm->ProjectionY("_y",1,1,"e")->Integral(0,hTriggerCorr_norm->GetNbinsY()+1)<<" MC fraction of triggered events in 0-multiplicity bin "<<fracZnum<<endl;

    Float_t fracZden = ((hVertexCorr_norm->ProjectionY("_y",2,hVertexCorr_norm->GetNbinsX()+1,"e"))->GetBinContent(ivtx+1))/Float_t(hVertexCorr_norm->ProjectionY("_y",2,hVertexCorr_norm->GetNbinsX()+1,"e")->Integral(0,hVertexCorr_norm->GetNbinsY()+1));  

    //cout<<" MC fraction of rec events (multSPD>0). Num "<<(hVertexCorr_norm->ProjectionY("_y",2,hVertexCorr_norm->GetNbinsX()+1,"e"))->GetBinContent(ivtx+1)<<" Den "<<hVertexCorr_norm->ProjectionY("_y",2,hVertexCorr_norm->GetNbinsX()+1,"e")->Integral(0,hVertexCorr_norm->GetNbinsY()+1)<<" fractions "<<fracZden<<endl;

    Float_t ratio =0.;
    Float_t ratioerr =0.;
    if (fracZden!=0.) {
      ratio = fracZnum/fracZden;
    } 
    if(ratio!=0.) {
      fHistEvtCorrFactor->SetBinContent(ivtx+1,ratio);
      fHistEvtCorrFactor->SetBinError(ivtx+1,ratioerr);
    }
  }

  TH1F* hprod = new TH1F("prod", "",20,-20.,20.);
  hprod->Sumw2();

  hprod->Multiply(fHistEvtCorrFactor,fHistZdistribTrigVtxEvt,1.,1.);

  //MC distributions

  Int_t histoNameRV;
  if (kPrimariesTest) {
    histoNameRV = 46; //MC pseudorapidity, MC vertex  
  } else {
    histoNameRV = 47; //MC pseudorapidity, rec vertex  
  }

  TH2F *Eta_RV = (TH2F*)l_MC->At(histoNameRV);
  TH2F *Eta = (TH2F*)l_MC->At(48);
//  Eta->Draw();

  TH1D *hdNdEta_RV = new TH1D("dNdEta_RV","Pseudorapidity ",60,-3,3);
  hdNdEta_RV=Eta_RV->ProjectionX("dNdEta_RV",binVtxStart+1,binVtxStop,"e");
  Double_t nEvtsRec=hSPDvertex->Integral(binVtxStart+1,binVtxStop);
  cout<<"Number of reconstructed events: "<<nEvtsRec <<endl;
  hdNdEta_RV->Scale(10./nEvtsRec);

  TH1D *hdNdEta_Trigg = new TH1D("dNdEta_Trigg","Pseudorapidity ",60,-3,3);
  TH2F *hEta_Trigg = new TH2F(*((TH2F*)l_MC->At(38))); 
  hdNdEta_Trigg=hEta_Trigg->ProjectionX("dNdEta_Trigg",0,-1,"e");
  Double_t nEvtsTrigg=0;
  nEvtsTrigg =hTriggerCorr_norm->Integral();
  hdNdEta_Trigg->Scale(10./nEvtsTrigg);
  cout<<"Number of triggered events: "<<nEvtsTrigg <<endl;
  TH1D *hdNdEta = new TH1D("dNdEta","Pseudorapidity ",60,-3.,3.);
  TH1D *hdNdEta10 = new TH1D("dNdEta10","Pseudorapidity ",60,-3.,3.);
  hdNdEta->Sumw2();
  hdNdEta = Eta->ProjectionX("dNdEta",0,-1,"e");
  hdNdEta10 = Eta->ProjectionX("dNdEta10",binVtxStart+1,binVtxStop,"e");
  Double_t nEvtsTot=hMCvertex->Integral(0,hMCvertex->GetNbinsX()+1);
  Double_t nEvtsTot10= hMCvertex->Integral(binVtxStart+1,binVtxStop);

  cout<<"Number of generated events in +-20cm: "<<nEvtsTot <<endl;
  cout<<"Number of generated events in +-10cm: "<<nEvtsTot10 <<endl;

  hdNdEta->Scale(10./nEvtsTot);
  hdNdEta10->Scale(10./nEvtsTot10);
  hdNdEta->GetXaxis()->SetTitle("Pseudorapidity #eta");
  hdNdEta->GetYaxis()->SetTitle("dN_{ch}/d#eta");
//  hdNdEta->Draw();

  TH1F *hdNdEtaNSD = (TH1F*)l_MC->At(51);
  hdNdEtaNSD->Sumw2();
  Double_t nEvtsNSD=0;
  TH1F *hproctype = (TH1F*)l_MC->At(52); 
  nEvtsNSD =hproctype->GetBinContent(2);
//   cout<<"NSD"<<endl;
  cout<<"Number of generated NSD events in +-20cm: "<<nEvtsNSD <<endl;
  hdNdEtaNSD->Scale(10./nEvtsNSD);

  for (Int_t ibin=0; ibin<60; ibin++) hdNdEtaNSD->SetBinError(ibin+1,0);

  TH1F *hdNdEtaMCInel = (TH1F*)l_MC->At(49);
  TH1F *hdNdEtaMCNSD = (TH1F*)l_MC->At(69);
  Double_t ntracksInel = hdNdEtaMCInel->GetBinContent(4);
  Double_t ntracksNSD = hdNdEtaMCNSD->GetBinContent(4); 
  cout<<"dNdEtaMCInel: "<<ntracksInel/nEvtsTot<<endl;
  cout<<"dNdEtaMCNSD: "<<ntracksNSD/nEvtsNSD<<endl;

  //Corrected distributions
  TH2F *hSPDEta_2D_bkgCorr =  new TH2F("SPDEta_2D_bkgCorr", "",60,-3.,3.,20,-20.,20.);
  hSPDEta_2D_bkgCorr->Sumw2();
  TH2F *hSPDEta_2D_bkgEffCorr =  new TH2F("SPDEta_2D_bkgEffCorr", "",60,-3.,3.,20,-20.,20.);
  hSPDEta_2D_bkgEffCorr->Sumw2();
  TH2F *hSPDEta_2D_bkgEffAccCorr =  new TH2F("SPDEta_2D_bkgEffAccCorr", "",120,-3.,3.,40,-20.,20.);
  hSPDEta_2D_bkgEffAccCorr->Sumw2();
  TH2F *hSPDEta_2D_fullyCorr =  new TH2F("SPDEta_2D_fullyCorr", "",60,-3.,3.,20,-20.,20.);
  hSPDEta_2D_fullyCorr->Sumw2();
  TH2F *hSPDEta_2D_triggeredEvts =  new TH2F("SPDEta_2D_triggeredEvts","",60,-3.,3.,20,-20.,20.);
  hSPDEta_2D_triggeredEvts->Sumw2();
  TH2F *hSPDEta_2D_allEvts =  new TH2F("SPDEta_2D_allEvts","",60,-3.,3.,20,-20.,20.);
  hSPDEta_2D_allEvts->Sumw2();
  TH2F *hSPDEta_2D_NSDEvts =  new TH2F("SPDEta_2D_NSDEvts","",60,-3.,3.,20,-20.,20.);
  hSPDEta_2D_NSDEvts->Sumw2();


  TH1F *hnorm_bkgEffAccCorr =  new TH1F("norm_bkgEffAccCorr", "",60,-3.,3.);
  hnorm_bkgEffAccCorr->Sumw2();
  TH1F *hnorm_fullyCorr =  new TH1F("norm_fullyCorr", "",60,-3.,3.);
  hnorm_fullyCorr->Sumw2();
  TH1F *hnorm_triggered =  new TH1F("norm_triggered", "",60,-3.,3.);
  hnorm_triggered->Sumw2();
  TH1F *hnorm =  new TH1F("norm", "",60,-3.,3.);
  hnorm->Sumw2();
  TH1F *hnormNSD =  new TH1F("normNSD", "",60,-3.,3.); 
  hnormNSD->Sumw2();

  TH1F *hnorm_gen =  new TH1F("norm_gen", "",60,-3.,3.);
  hnorm_gen->Sumw2();

  //dN/deta
  TH1F *hdNdEta_raw = new TH1F("dNdEta_raw","", 60,-3.,3.);
  hdNdEta_raw->Sumw2();
  hdNdEta_raw->Add(hSPDEta_2Draw->ProjectionX("SPDEta_2Draw_x",binVtxStart+1,binVtxStop,"e"));
  hdNdEta_raw->Scale(10./nEvtsRec);
  
  TH1F *hdNdEta_bkgCorr = new TH1F("dNdEta_bkgCorr","", 60,-3.,3.);
  hdNdEta_bkgCorr->Sumw2();
  TH1F *hdNdEta_bkgEffCorr = new TH1F("dNdEta_bkgEffCorr","", 60,-3.,3.);
  hdNdEta_bkgEffCorr->Sumw2();
  TH1F *hdNdEta_bkgEffAccCorr = new TH1F("dNdEta_bkgEffAccCorr","", 60,-3.,3.);
  hdNdEta_bkgEffAccCorr->Sumw2();
  TH1F *hdNdEta_fullyCorr = new TH1F("dNdEta_fullyCorr","", 60,-3.,3.);
  hdNdEta_fullyCorr->Sumw2();
  TH1F *hdNdEta_triggeredEvts = new TH1F("dNdEta_triggeredEvts","", 60,-3.,3.);
  hdNdEta_triggeredEvts->Sumw2();
  TH1F *hdNdEta_allEvts = new TH1F("dNdEta_allEvts","", 60,-3.,3.);
  hdNdEta_allEvts->Sumw2();
  TH1F *hdNdEta_NSDEvts = new TH1F("dNdEta_NSDEvts","", 60,-3.,3.);
  hdNdEta_NSDEvts->Sumw2();

  //Cutting either on weights or on relative errors of weights
  // Acceptance Pb-Pb
  for (Int_t jeta=0; jeta<120; jeta++) {
    for (Int_t kvtx=0; kvtx<40; kvtx++) {
      if (hAccCorrPb->GetBinContent(jeta+1,kvtx+1) < cutAccPb) {
        hAccCorrPb->SetBinContent(jeta+1,kvtx+1,0.);
        hAccCorrPb->SetBinError(jeta+1,kvtx+1,0.);
      }
    }
  }
  for (Int_t jeta=0; jeta<60; jeta++) {
    for (Int_t kvtx=0; kvtx<20; kvtx++) {
      //Acceptance p-p
      Float_t relErr = hAccCorr->GetBinError(jeta+1,kvtx+1)/hAccCorr->GetBinContent(jeta+1,kvtx+1);
      if (relErr>cutAcc) {
        hAccCorr->SetBinContent(jeta+1,kvtx+1,0.);
        hAccCorr->SetBinError(jeta+1,kvtx+1,0.);
      }
      // Background
      Float_t relBinErr = hBackgroundCorr->GetBinError(jeta+1,kvtx+1)/hBackgroundCorr->GetBinContent(jeta+1,kvtx+1);
      if (relBinErr>cutBkg) {
        hBackgroundCorr->SetBinContent(jeta+1,kvtx+1,0.);
        hBackgroundCorr->SetBinError(jeta+1,kvtx+1,0.);
      }
      //Alg eff
      relBinErr = hTrackletRecEff->GetBinError(jeta+1,kvtx+1)/hTrackletRecEff->GetBinContent(jeta+1,kvtx+1);
      if (relBinErr>cutAlgEff) {
        hTrackletRecEff->SetBinContent(jeta+1,kvtx+1,0.);
        hTrackletRecEff->SetBinError(jeta+1,kvtx+1,0.);
      }
      //Dec particles
      relBinErr = hDecPartLay2Corr->GetBinError(jeta+1,kvtx+1)/hDecPartLay2Corr->GetBinContent(jeta+1,kvtx+1);
      if (relBinErr>cutDec) {
        hDecPartLay2Corr->SetBinContent(jeta+1,kvtx+1,0.);
        hDecPartLay2Corr->SetBinError(jeta+1,kvtx+1,0.);
      }
      //Vtx Corr
      relBinErr = hVertexCorrEta->GetBinError(jeta+1,kvtx+1)/hVertexCorrEta->GetBinContent(jeta+1,kvtx+1);
      if (relBinErr>cutVtx) {
        hVertexCorrEta->SetBinContent(jeta+1,kvtx+1,0.);
        hVertexCorrEta->SetBinError(jeta+1,kvtx+1,0.);
      }
      relBinErr = hTriggerCorrEta->GetBinError(jeta+1,kvtx+1)/hTriggerCorrEta->GetBinContent(jeta+1,kvtx+1);
      if (relBinErr>cutTrig) {
        hTriggerCorrEta->SetBinContent(jeta+1,kvtx+1,0.);
        hTriggerCorrEta->SetBinError(jeta+1,kvtx+1,0.);
      }
      relBinErr = hTrigCorrEtaNSD->GetBinError(jeta+1,kvtx+1)/hTrigCorrEtaNSD->GetBinContent(jeta+1,kvtx+1);
      if (relBinErr>cutTrig) {
        hTrigCorrEtaNSD->SetBinContent(jeta+1,kvtx+1,0.);
        hTrigCorrEtaNSD->SetBinError(jeta+1,kvtx+1,0.);
      }
      if (hAccCorr->GetBinContent(jeta+1,kvtx+1)==0) {
        hDecPartLay2Corr->SetBinContent(jeta+1,kvtx+1,0.);  //show these correction only in bins inside the SPD acceptance
        hDecPartLay2Corr->SetBinError(jeta+1,kvtx+1,0.);
        hVertexTrigCorrEta->SetBinContent(jeta+1,kvtx+1,0.);
        hVertexTrigCorrEta->SetBinError(jeta+1,kvtx+1,0.);
        hVertexTrigCorrEtaNSD->SetBinContent(jeta+1,kvtx+1,0.);
        hVertexTrigCorrEtaNSD->SetBinError(jeta+1,kvtx+1,0.);
        hTrigCorrEtaNSDInv->SetBinContent(jeta+1,kvtx+1,0.);
        hTrigCorrEtaNSDInv->SetBinError(jeta+1,kvtx+1,0.);

      }
    } 
  }

  //Apply corrections at
  //...track level
  hSPDEta_2D_bkgCorr->Add(hSPDEta_2Draw);
  hSPDEta_2D_bkgCorr->Multiply(hBackgroundCorr);
  hSPDEta_2D_bkgEffCorr->Divide(hSPDEta_2D_bkgCorr,hTrackletRecEff,1.,1.);
  if (kAccPb) {
    hSPDEta_2D_bkgEffAccCorr->Divide(hSPDEta_2Draw_accBin,hAccCorrPb,1.,1.);
    hSPDEta_2D_bkgEffAccCorr->Rebin2D(2,2,"");
  } else {
    hSPDEta_2D_bkgEffAccCorr->Rebin2D(2,2,"");
    hSPDEta_2D_bkgEffAccCorr->Divide(hSPDEta_2Draw,hAccCorr,1.,1.);
//    hSPDEta_2D_bkgEffAccCorr->Divide(hSPDEta_2Draw_accBin,hAccCorr,1.,1.);
  }
//  hSPDEta_2D_bkgEffAccCorr->Rebin2D(2,2,"");
  hSPDEta_2D_bkgEffAccCorr->Multiply(hBackgroundCorr);
  hSPDEta_2D_bkgEffAccCorr->Divide(hTrackletRecEff);
  hSPDEta_2D_fullyCorr->Add(hSPDEta_2D_bkgEffAccCorr);
  hSPDEta_2D_fullyCorr->Divide(hDecPartLay2Corr);
  hSPDEta_2D_triggeredEvts->Add(hSPDEta_2D_fullyCorr);
  hSPDEta_2D_triggeredEvts->Divide(hVertexCorrEta);
  hSPDEta_2D_allEvts->Add(hSPDEta_2D_triggeredEvts);
  hSPDEta_2D_allEvts->Divide(hTriggerCorrEta);
  hSPDEta_2D_NSDEvts->Add(hSPDEta_2D_triggeredEvts);
  hSPDEta_2D_NSDEvts->Divide(hTrigCorrEtaNSD);

  //...event level
  Double_t evtTrigNoTracklets = hESDVertexvsNtracklets_trigEvts->ProjectionY("_y",1,1,"e")->Integral(0,hESDVertexvsNtracklets_trigEvts->GetNbinsX()+1); 
  cout<<"Number of events triggered with tracklet multiplicity = 0: "<<evtTrigNoTracklets<<endl;
  for (Int_t ivtx=0;ivtx<20;ivtx++) { 
    hESDVertexvsNtracklets_trigEvts->SetBinContent(1,ivtx+1,(hprod->GetBinContent(ivtx+1))*(evtTrigNoTracklets));
    hESDVertexvsNtracklets_trigEvts->SetBinError(1,ivtx+1,TMath::Sqrt((TMath::Power(hprod->GetBinContent(ivtx+1),2)*evtTrigNoTracklets)+(TMath::Power(evtTrigNoTracklets,2)*TMath::Power(hprod->GetBinError(ivtx+1),2))));
//    cout<<"bincont "<<hESDVertexvsNtracklets_trigEvts->GetBinContent(1,ivtx+1)<<" binerr "<<hESDVertexvsNtracklets_trigEvts->GetBinError(1,ivtx+1)<<endl;
  }

  hESDVertexvsNtracklets_Corr->Add(hESDVertexvsNtracklets_trigEvts);
  hESDVertexvsNtracklets_Corr->Divide(hTriggerCorr);
  hESDVertexvsNtracklets_CorrNSD->Add(hESDVertexvsNtracklets_trigEvts);
  hESDVertexvsNtracklets_CorrNSD->Divide(hTriggerCorrNSD);

  //Filling normalization histograms

  TH1D *hSPDvertexCorr_y = new TH1D("SPDvertexCorr_y","",20,-20,20);
  hSPDvertexCorr_y = hESDVertexvsNtracklets_Corr->ProjectionY("SPDvertexCorr_y",0,-1,"e");
  cout<<"Number of events (tracklets corrected)"<<hSPDvertexCorr_y->Integral(binVtxStart+1,binVtxStop)<<endl;  
  for (Int_t ivtx=binVtxStart; ivtx<binVtxStop; ivtx++) {
    Double_t nEvtsivtx = hSPDvertex->GetBinContent(ivtx+1);
    Double_t nEvtsivtxCorr = hSPDvertexCorr_y->GetBinContent(ivtx+1);
    Double_t nEvtsivtxCorr_triggered = hESDVertexvsNtracklets_trigEvts->ProjectionY("_y",0,-1,"e")->GetBinContent(ivtx+1);

    Double_t nEvtsivtxCorr_NSD = hESDVertexvsNtracklets_CorrNSD->ProjectionY("_y",0,-1,"e")->GetBinContent(ivtx+1);

    //Double_t nEvtGen = hMCvertex->GetBinContent(ivtx+1);
    //cout<<"iBinVtx: "<<ivtx<<" nEvtsGen-> "<<nEvtGen<<" nEvtsCorr-> "<<nEvtsivtxCorr<<endl;
    for (Int_t jeta=0; jeta<60; jeta++) {
      if (hSPDEta_2D_bkgEffAccCorr->GetBinContent(jeta+1,ivtx+1)!=0.)
        hnorm_bkgEffAccCorr->Fill(hBackgroundCorr->GetXaxis()->GetBinCenter(jeta+1),nEvtsivtx);
      if (hSPDEta_2D_fullyCorr->GetBinContent(jeta+1,ivtx+1)!=0.)
        hnorm_fullyCorr->Fill(hBackgroundCorr->GetXaxis()->GetBinCenter(jeta+1),nEvtsivtx);
      if (hSPDEta_2D_triggeredEvts->GetBinContent(jeta+1,ivtx+1)!=0.)
        hnorm_triggered->Fill(hBackgroundCorr->GetXaxis()->GetBinCenter(jeta+1),nEvtsivtxCorr_triggered);
      if (hSPDEta_2D_allEvts->GetBinContent(jeta+1,ivtx+1)!=0.)
        hnorm->Fill(hBackgroundCorr->GetXaxis()->GetBinCenter(jeta+1),nEvtsivtxCorr);
      if (hSPDEta_2D_NSDEvts->GetBinContent(jeta+1,ivtx+1)!=0.)
        hnormNSD->Fill(hBackgroundCorr->GetXaxis()->GetBinCenter(jeta+1),nEvtsivtxCorr_NSD);
    }
  }
  
  for (Int_t jeta=0; jeta<60; jeta++) {
    hnorm_bkgEffAccCorr->SetBinError(jeta+1,TMath::Sqrt(hnorm_bkgEffAccCorr->GetBinContent(jeta+1)));
    hnorm_fullyCorr->SetBinError(jeta+1,TMath::Sqrt(hnorm_fullyCorr->GetBinContent(jeta+1)));
    hnorm_triggered->SetBinError(jeta+1,TMath::Sqrt(hnorm_triggered->GetBinContent(jeta+1)));
    hnorm->SetBinError(jeta+1,TMath::Sqrt(hnorm->GetBinContent(jeta+1))); 
    hnormNSD->SetBinError(jeta+1,TMath::Sqrt(hnorm->GetBinContent(jeta+1)));  
  }

  hdNdEta_bkgCorr->Add(hSPDEta_2D_bkgCorr->ProjectionX("SPDEta_2D_bkgCorr_x",binVtxStart+1,binVtxStop,"e"));
  hdNdEta_bkgCorr->Scale(1./nEvtsRec);
  hdNdEta_bkgEffCorr->Add(hSPDEta_2D_bkgEffCorr->ProjectionX("SPDEta_2D_bkgEffCorr_x",binVtxStart+1,binVtxStop,"e"));
  hdNdEta_bkgEffCorr->Scale(1./nEvtsRec);
  hdNdEta_bkgEffAccCorr->Divide(hSPDEta_2D_bkgEffAccCorr->ProjectionX("SPDEta_2D_bkgEffAccCorr_x",binVtxStart+1,binVtxStop,"e"),hnorm_bkgEffAccCorr);
  hdNdEta_fullyCorr->Divide(hSPDEta_2D_fullyCorr->ProjectionX("SPDEta_2D_fullyCorr_x",binVtxStart+1,binVtxStop,"e"),hnorm_fullyCorr);
  hdNdEta_triggeredEvts->Divide(hSPDEta_2D_triggeredEvts->ProjectionX("SPDEta_2D_triggeredEvts_x",binVtxStart+1,binVtxStop,"e"),hnorm_triggered);
  hdNdEta_allEvts->Divide(hSPDEta_2D_allEvts->ProjectionX("SPDEta_2D_allEvts_x",binVtxStart+1,binVtxStop,"e"),hnorm);
  hdNdEta_NSDEvts->Divide(hSPDEta_2D_NSDEvts->ProjectionX("SPDEta_2D_NSDEvts_x",binVtxStart+1,binVtxStop,"e"),hnormNSD);

  hdNdEta_bkgCorr->Scale(10.);
  hdNdEta_bkgEffCorr->Scale(10.);
  hdNdEta_bkgEffAccCorr->Scale(10.);
  hdNdEta_fullyCorr->Scale(10.);
  hdNdEta_triggeredEvts->Scale(10.);
  hdNdEta_allEvts->Scale(10.);
  hdNdEta_NSDEvts->Scale(10.);

  TH1F* hRationorm=new TH1F("rationorm","",60,-3,3);
  hRationorm->Add(hnorm);
  hRationorm->Scale(1./nEvtsTot10);

  TH1D *hdNdEta_test = new TH1D("dNdEta_test","Pseudorapidity ",60,-3,3);  
  TH2F *hMask_allEvts = new TH2F("Mask_allEvts","",60,-3,3,20,-20,20);
  hMask_allEvts->Sumw2();

  TH2F *Eta_test = new TH2F("Eta_test","",60,-3.,3.,20,-20.,20.);
  if (kPrimariesTest) {
    for (Int_t kvtx=0; kvtx<20; kvtx++) {
      for (Int_t jeta=0; jeta<60; jeta++) {
        if (hSPDEta_2D_allEvts->GetBinContent(jeta+1,kvtx+1)!=0) {
          hMask_allEvts->SetBinContent(jeta+1,kvtx+1,1);
        }
      }
    }
    //dN/dEta gen test
    Eta_test->Multiply(Eta,hMask_allEvts,1.,1.);
    hdNdEta_test->Divide(Eta_test->ProjectionX("eta_x",binVtxStart+1,binVtxStop,"e"),hnorm,1.,1.);
    hdNdEta_test->Scale(10.);
  }

  // Generated/corrected ratios
  TH2F* hRatiodNdEta2D=new TH2F("ratiodNdEta2D","",60,-3,3,20,-20,20);
//  hRatiodNdEta2D->Divide(hSPDEta_2D_fullyCorr,Eta_RV,1.,1.);
  if (kPrimariesTest) hRatiodNdEta2D->Divide(hSPDEta_2D_allEvts,Eta_test,1.,1.);
  else hRatiodNdEta2D->Divide(hSPDEta_2D_allEvts,Eta,1.,1.);

  TH1F* hRatiodNdEta=new TH1F("ratiodNdEta","",60,-3,3);
  hRatiodNdEta->GetXaxis()->SetTitle("Pseudorapidity #eta");
  hRatiodNdEta->GetYaxis()->SetTitle("Generated/Corrected");

  TH1F* hRatiodNdEta10=new TH1F("ratiodNdEta10","",60,-3,3);
  hRatiodNdEta10->Divide(hdNdEta,hdNdEta10);
//  new TCanvas();
//  hRatiodNdEta10->Draw();

  TH1F* hRatiodNdEtaNSD=new TH1F("ratiodNdEtaNSD","",60,-3,3);
  hRatiodNdEtaNSD->GetXaxis()->SetTitle("Pseudorapidity #eta");
  hRatiodNdEtaNSD->GetYaxis()->SetTitle("Generated/Corrected");

  TH1F* hRatiodNdEta_vertex=new TH1F("ratiodNdEta_vertex","",60,-3,3);
  TH1F* hRatiodNdEta_triggered=new TH1F("ratiodNdEta_triggered","",60,-3,3);

  for (Int_t i=0;i<60;i++) hdNdEta->SetBinError(i+1,0.);
  for (Int_t i=0;i<60;i++) hdNdEta_RV->SetBinError(i+1,0.);

  if (kPrimariesTest) hRatiodNdEta->Divide(hdNdEta_test,hdNdEta_allEvts,1.,1.); 
  else hRatiodNdEta->Divide(hdNdEta,hdNdEta_allEvts,1.,1.);
  hRatiodNdEtaNSD->Divide(hdNdEtaNSD,hdNdEta_NSDEvts,1.,1.);
  hRatiodNdEta_vertex->Divide(hdNdEta_RV,hdNdEta_fullyCorr,1.,1.);
  hRatiodNdEta_triggered->Divide(hdNdEta_Trigg,hdNdEta_triggeredEvts,1.,1.);
 
  Float_t max = 14 ;

  //Draw a canvas with correction maps at track level
  TCanvas *ccorr = new TCanvas("c1","Tracklet level corrections and data");
  ccorr->Divide(3,2);
  ccorr->cd(1);
  hBackgroundCorr->GetXaxis()->SetTitle("Pseudorapidity #eta");
  hBackgroundCorr->GetYaxis()->SetTitle("z_{vertex} [cm]");
  hBackgroundCorr->GetXaxis()->SetRangeUser(-2,2);
  hBackgroundCorr->GetYaxis()->SetRangeUser(-10.,8.);
  hBackgroundCorr->SetMinimum(0.);
  hBackgroundCorr->SetMaximum(1.);
  hBackgroundCorr->Draw("colz");
  ccorr->cd(2);
  hTrackletRecEff->GetXaxis()->SetTitle("Pseudorapidity #eta");
  hTrackletRecEff->GetYaxis()->SetTitle("z_{vertex} [cm]");
  hTrackletRecEff->GetXaxis()->SetRangeUser(-2,2);
  hTrackletRecEff->GetYaxis()->SetRangeUser(-10.,8.);
  hTrackletRecEff->Draw("colz");
  ccorr->cd(3);
  hAccCorr->GetXaxis()->SetTitle("Pseudorapidity #eta");
  hAccCorr->GetYaxis()->SetTitle("z_{vertex} [cm]");
  hAccCorr->GetXaxis()->SetRangeUser(-2,2);
  hAccCorr->GetYaxis()->SetRangeUser(-10.,8.);
  hAccCorr->Draw("colz");
  ccorr->cd(4);
  hDecPartLay2Corr->GetXaxis()->SetTitle("Pseudorapidity #eta");
  hDecPartLay2Corr->GetYaxis()->SetTitle("z_{vertex} [cm]");
  hDecPartLay2Corr->GetXaxis()->SetRangeUser(-2,2);
  hDecPartLay2Corr->GetYaxis()->SetRangeUser(-10.,8.);
  hDecPartLay2Corr->SetMinimum(0.9);
  hDecPartLay2Corr->SetMaximum(1.);
  hDecPartLay2Corr->Draw("colz");
  ccorr->cd(5);
  hVertexTrigCorrEta->GetXaxis()->SetTitle("Pseudorapidity #eta");
  hVertexTrigCorrEta->GetYaxis()->SetTitle("z_{vertex} [cm]");
  hVertexTrigCorrEta->GetXaxis()->SetRangeUser(-2,2);
  hVertexTrigCorrEta->GetYaxis()->SetRangeUser(-10.,8.);
  hVertexTrigCorrEta->SetMinimum(0.99);
  hVertexTrigCorrEta->SetMaximum(1.01);
  hVertexTrigCorrEta->Draw("colz");
  ccorr->cd(6);
  hSPDEta_2Draw->GetXaxis()->SetTitle("Pseudorapidity #eta");
  hSPDEta_2Draw->GetYaxis()->SetTitle("z_{vertex} [cm]");
  hSPDEta_2Draw->GetXaxis()->SetRangeUser(-2,2);
  hSPDEta_2Draw->GetYaxis()->SetRangeUser(-10.,8.);
  hSPDEta_2Draw->Draw("colz");

  //Draw a canvas with correction maps at event level
  TCanvas *ccorrev = new TCanvas("c2","Event level corrections and data");
  ccorrev->Divide(3,3);
  ccorrev->cd(1);
  hVertexCorr->GetXaxis()->SetTitle("Tracklet multiplicity");
  hVertexCorr->GetYaxis()->SetTitle("z_{MCvtx} [cm]");
  hVertexCorr->GetXaxis()->SetRangeUser(0.,20.);
  hVertexCorr->Draw("colz");
  ccorrev->cd(2);
  hTriggerCorr->GetXaxis()->SetTitle("Tracklet multiplicity");
  hTriggerCorr->GetYaxis()->SetTitle("z_{MCvtx} [cm]");
  hTriggerCorr->GetXaxis()->SetRangeUser(0.,20.);
  hTriggerCorr->Draw("colz");
  ccorrev->cd(3);
  hTriggerCorrNSD->GetXaxis()->SetTitle("Tracklet multiplicity");
  hTriggerCorrNSD->GetYaxis()->SetTitle("z_{MCvtx} [cm]");
  hTriggerCorrNSD->GetXaxis()->SetRangeUser(0.,20.);
  hTriggerCorrNSD->Draw("colz");
  ccorrev->cd(4);
  fHistZdistribTrigVtxEvt->GetYaxis()->SetTitle("Fraction of rec events per vertex bin");
  fHistZdistribTrigVtxEvt->GetXaxis()->SetTitle("z_{SPDvtx} [cm]");
  fHistZdistribTrigVtxEvt->Draw();
  ccorrev->cd(5);
  fHistEvtCorrFactor->GetYaxis()->SetTitle("Correction weight");
  fHistEvtCorrFactor->GetXaxis()->SetTitle("z_{MCvtx} [cm]");
  fHistEvtCorrFactor->Draw();   
  ccorrev->cd(6);
  hESDVertexvsNtracklets->GetYaxis()->SetTitle("z_{SPDvtx} [cm]");
  hESDVertexvsNtracklets->GetXaxis()->SetTitle("Tracklet mulitplicity");
  hESDVertexvsNtracklets->GetXaxis()->SetRangeUser(0.,20.);
  hESDVertexvsNtracklets->Draw("colz");
  ccorrev->cd(7);
  hESDVertexvsNtracklets_trigEvts->GetYaxis()->SetTitle("z_{SPDvtx} [cm]");
  hESDVertexvsNtracklets_trigEvts->GetXaxis()->SetTitle("Tracklet mulitplicity");
  hESDVertexvsNtracklets_trigEvts->GetXaxis()->SetRangeUser(0.,20.);
  hESDVertexvsNtracklets_trigEvts->Draw("colz"); //bin mult=0 already distributed in z

  // Draw dN/dEta
  hdNdEta_raw->SetLineColor(kGreen);
  hdNdEta_raw->SetLineWidth(3);

  hdNdEta_RV->SetLineWidth(3);
  hdNdEta_RV->SetMinimum(0.);
  hdNdEta_RV->SetMaximum(max);
  hdNdEta_RV->SetLineColor(kRed);
//  hdNdEta_RV->Draw("histo");
  hdNdEta->SetLineWidth(3);
  hdNdEta->SetMinimum(0.);
  hdNdEta->SetMaximum(max);

  new TCanvas();
  hdNdEta->Draw("histo");
//  hdNdEta_RV->Draw("histo,same");

  hdNdEta_raw->Draw("histo,same");

  hdNdEta_bkgCorr->SetMarkerStyle(21);
  hdNdEta_bkgCorr->SetMarkerColor(kBlue);

  hdNdEta_bkgCorr->Draw("same,p,histo");

  hdNdEta_bkgEffCorr->SetMarkerStyle(22);
  hdNdEta_bkgEffCorr->SetMarkerColor(kGreen);

  hdNdEta_bkgEffAccCorr->SetMarkerStyle(23);
  hdNdEta_bkgEffAccCorr->SetMarkerColor(kRed);

  hdNdEta_bkgEffCorr->Draw("p,same,e");
  hdNdEta_bkgEffAccCorr->Draw("p,same,e");
  hdNdEta_fullyCorr->SetMarkerStyle(20);
  hdNdEta_fullyCorr->SetMarkerColor(kBlue);
  hdNdEta_fullyCorr->Draw("p,same,e");

  hdNdEta_triggeredEvts->SetMarkerStyle(21);
  hdNdEta_triggeredEvts->SetMarkerColor(kOrange);
  hdNdEta_triggeredEvts->Draw("p,same,e");

  hdNdEta_allEvts->SetMarkerStyle(21);
  hdNdEta_allEvts->SetMarkerColor(kRed);
  hdNdEta_allEvts->Draw("p,same,e");

  TLegend *leg1 = new TLegend(0.2,0.7,0.7,0.9,NULL,"brNDC");
  leg1->SetFillColor(0);
  leg1->SetTextFont(62);
  leg1->SetTextSize(0.0304569);
  TLegendEntry *entry1=leg1->AddEntry(hdNdEta_raw,"Reconstructed","l");
                entry1=leg1->AddEntry(hdNdEta,"Generated","l");
                entry1=leg1->AddEntry(hdNdEta_bkgCorr,"Bkg corrected","p");
                entry1=leg1->AddEntry(hdNdEta_bkgEffCorr,"Bkg+AlgEff+SPDEff corrected","p");
                entry1=leg1->AddEntry(hdNdEta_bkgEffAccCorr,"Bkg+AlgEff+SPDAccEff corrected","p");
                entry1=leg1->AddEntry(hdNdEta_fullyCorr,"Corrected (triggered events with vertex)","p");
                entry1=leg1->AddEntry(hdNdEta_triggeredEvts,"Corrected (triggered events)","p");
                entry1=leg1->AddEntry(hdNdEta_allEvts,"Corrected (all events)","p");

  leg1->Draw();

  new TCanvas();

  hdNdEta->Draw("histo");
  hdNdEta_NSDEvts->SetMarkerStyle(21);
  hdNdEtaNSD->SetLineColor(kGreen);
  hdNdEta_allEvts->Draw("p,same,e");
  hdNdEta_NSDEvts->Draw("same,p,e");
  hdNdEtaNSD->Draw("histo,same");
  
  TLegend *leg2 = new TLegend(0.3,0.7,0.7,0.9,NULL,"brNDC");
  leg2->SetFillColor(0);
  leg2->SetTextFont(62);

  TLegendEntry *entry2=leg2->AddEntry(hdNdEta_allEvts,"Corrected INEL","p");
                entry2=leg2->AddEntry(hdNdEta_NSDEvts,"Corrected NSD","p");
                entry2=leg2->AddEntry(hdNdEta,"Generated INEL","l");
                entry2=leg2->AddEntry(hdNdEtaNSD,"Generated NSD","l");

  leg2->Draw();

  //Trigger-vertex efficiency 
  TH1F* fHistMultAllNonDiff = (TH1F*)l_MC->At(63);
  TH1F* fHistMultAllSingleDiff = (TH1F*)l_MC->At(64);
  TH1F* fHistMultAllDoubleDiff = (TH1F*)l_MC->At(65);
  TH1F* fHistMultTrVtxNonDiff = (TH1F*)l_MC->At(66);
  TH1F* fHistMultTrVtxSingleDiff = (TH1F*)l_MC->At(67);
  TH1F* fHistMultTrVtxDoubleDiff = (TH1F*)l_MC->At(68);

  TH1F* fHistTrigVtxEffvsMultEtacutNonDiff = new TH1F("fHistTrigVtxEffvsMultEtacutNonDiff","",20,-0.5,19.5);
  fHistTrigVtxEffvsMultEtacutNonDiff->GetXaxis()->SetTitle("Multiplicity in |#eta|<1.4");
  fHistTrigVtxEffvsMultEtacutNonDiff->GetYaxis()->SetTitle("Efficiency");
  fHistTrigVtxEffvsMultEtacutNonDiff->Sumw2();
  TH1F* fHistTrigVtxEffvsMultEtacutSingleDiff = new TH1F("fHistTrigVtxEffvsMultEtacutSingleDiff","",20,-0.5,19.5);
  fHistTrigVtxEffvsMultEtacutSingleDiff->Sumw2();
  TH1F* fHistTrigVtxEffvsMultEtacutDoubleDiff = new TH1F("fHistTrigVtxEffvsMultEtacutDoubleDiff","",20,-0.5,19.5);
  fHistTrigVtxEffvsMultEtacutDoubleDiff->Sumw2();
  fHistTrigVtxEffvsMultEtacutNonDiff->Divide(fHistMultTrVtxNonDiff,fHistMultAllNonDiff,1.,1.);
  fHistTrigVtxEffvsMultEtacutSingleDiff->Divide(fHistMultTrVtxSingleDiff,fHistMultAllSingleDiff,1.,1.);
  fHistTrigVtxEffvsMultEtacutDoubleDiff->Divide(fHistMultTrVtxDoubleDiff,fHistMultAllDoubleDiff,1.,1.);

  new TCanvas();

  fHistTrigVtxEffvsMultEtacutNonDiff->SetMarkerStyle(2);
  fHistTrigVtxEffvsMultEtacutSingleDiff->SetMarkerStyle(25);
  fHistTrigVtxEffvsMultEtacutDoubleDiff->SetMarkerStyle(24);
  fHistTrigVtxEffvsMultEtacutNonDiff->GetYaxis()->SetRangeUser(0.,1.1);
  fHistTrigVtxEffvsMultEtacutNonDiff->GetXaxis()->SetRangeUser(0.,10);
  fHistTrigVtxEffvsMultEtacutNonDiff->Draw("p,histo");
  fHistTrigVtxEffvsMultEtacutSingleDiff->Draw("same,p,histo");
  fHistTrigVtxEffvsMultEtacutDoubleDiff->Draw("same,p,histo");
 
  TLegend *leg3 = new TLegend(0.3,0.2,0.7,0.4444,NULL,"brNDC");
  leg3->SetTextFont(62);
  leg3->SetFillColor(0);
  TLegendEntry *entry3=leg3->AddEntry(fHistTrigVtxEffvsMultEtacutNonDiff,"Non-diffractive","p");
                entry3=leg3->AddEntry(fHistTrigVtxEffvsMultEtacutSingleDiff,"Single-diffractive","p");
                entry3=leg3->AddEntry(fHistTrigVtxEffvsMultEtacutDoubleDiff,"Double-diffractive","p");
  
  leg3->Draw();

  TH1F* hProcessTypeTriggered = (TH1F*)l_MC->At(53); 
  Float_t nInelTriggered = hProcessTypeTriggered->GetBinContent(3);
//  Float_t nNSDTriggered = hProcessTypeTriggered->GetBinContent(2);
  Float_t nNDTriggered = hProcessTypeTriggered->GetBinContent(1);
  Float_t nSDTriggered = hProcessTypeTriggered->GetBinContent(4);
  Float_t nDDTriggered = hProcessTypeTriggered->GetBinContent(5);
  Float_t nND = hproctype->GetBinContent(1);
  Float_t nSD = hproctype->GetBinContent(4);
  Float_t nDD = hproctype->GetBinContent(5);
 
  cout<<"Trigger efficiencies: Inel-->"<<nInelTriggered/nEvtsTot<<" ND-->"<<nNDTriggered/nND<<" SD-->"<<nSDTriggered/nSD<<" DD-->"<<nDDTriggered/nDD<<endl;

  new TCanvas();
  hRatiodNdEtaNSD->GetYaxis()->SetRangeUser(0.9,1.1);
  hRatiodNdEtaNSD->GetXaxis()->SetRangeUser(-2.,1.9);
  hRatiodNdEtaNSD->Draw();
  new TCanvas();
  hRatiodNdEta->GetYaxis()->SetRangeUser(0.9,1.1);
  hRatiodNdEta->GetXaxis()->SetRangeUser(-2.,1.9);
  hRatiodNdEta->Draw();


  // Write histos

  TFile *fout= new TFile("SPDdNdEta.root","RECREATE"); 

  hAccCorr->Write();
  hAccCorrPb->Write();

  hBackgroundCorr->Write();
  hTrackletRecEff->Write();
  hDecPartLay2Corr->Write();

  hSPDEta_2Draw_accBin->Write();
  hSPDEta_2Draw->Write();

  hSPDvertex->Write();
  hSPDvertexCorr_y->Write();
  hSPDEta_2D_bkgCorr->Write();
  hSPDEta_2D_bkgEffCorr->Write();
  hSPDEta_2D_bkgEffAccCorr->Write();
  hSPDEta_2D_fullyCorr->Write();
  hSPDEta_2D_triggeredEvts->Write();
  hSPDEta_2D_allEvts->Write();
  hSPDEta_2D_NSDEvts->Write();

  hnorm_bkgEffAccCorr->Write();
  hnorm_triggered->Write();
  hnorm->Write();
  hnormNSD->Write();
  hnorm_gen->Write();

  hVertexCorr->Write();
  hTriggerCorr->Write();
  hTriggerCorrNSD->Write();
  hTriggerCorrNSD2->Write();

  hVertexCorrEta->Write();
  hTriggerCorrEta->Write();
  hVertexTrigCorrEta->Write();
  hVertexTrigCorrEtaNSD->Write();
  hTrigCorrEtaNSD->Write();
  hTrigCorrEtaNSDInv->Write();

  fHistZdistribTrigVtxEvt->Write();
  fHistEvtCorrFactor->Write();

  hESDVertexvsNtracklets->Write();
  hESDVertexvsNtracklets_trigEvts->Write();
  hESDVertexvsNtracklets_Corr->Write();
  hESDVertexvsNtracklets_CorrNSD->Write();
 
  hdNdEta->Write();
  hdNdEta_RV->Write();
  hdNdEta_bkgCorr->Write();
  hdNdEta_bkgEffCorr->Write();
  hdNdEta_bkgEffAccCorr->Write();
  hdNdEta_fullyCorr->Write();
  hdNdEta_triggeredEvts->Write();
  hdNdEta_allEvts->Write();
  hdNdEta_NSDEvts->Write();

  hdNdEta_raw->Write();

  hRatiodNdEta->Write();
  hRatiodNdEtaNSD->Write();
  hRatiodNdEta_triggered->Write();
  hRatiodNdEta_vertex->Write();
  hRatiodNdEta2D->Write();
  hRationorm->Write();

  Eta_test->Write();
  Eta->Write(); 
  Eta_RV->Write();

  hprod->Write();
  hdNdEta10->Write();
  hRatiodNdEta10->Write();

  fHistTrigVtxEffvsMultEtacutNonDiff->Write(); 
  fHistTrigVtxEffvsMultEtacutSingleDiff->Write();
  fHistTrigVtxEffvsMultEtacutDoubleDiff->Write();
  
  fout->Write();
  fout->Close();

}
