/*************************************************************************
* Macro correctSPDdNdEtapp                                               *
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
                         Double_t cutAlgEff, Double_t cutDec, Double_t cutVtx)  {

//lim vtx: 5-15 ->[-10cm,10cm[; 0-20 ->[-20cm,20cm[  

  gStyle->SetOptLogy(kFALSE);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetHistLineWidth(2);
  gStyle->SetLabelSize(0.05, "x");
  gStyle->SetLabelSize(0.05, "y");
  gStyle->SetTitleSize(0.05, "x");
  gStyle->SetTitleSize(0.05, "y");
  gStyle->SetTitleOffset(1.1, "x");
  gStyle->SetPadBottomMargin(0.14);
//  gStyle->SetPalette(1,0);
  gStyle->SetFillColor(0);


  //Input files
  TFile *f_accPb = new TFile(fileAccPb);
  TFile *f_corr = new TFile(fileCorr);
  TFile *f_MC = new TFile(fileMC);
  TFile *f_raw = new TFile(fileRaw);

  TList *l_corr = (TList*)f_corr->Get("cOutput");
  TList *l_MC = (TList*)f_MC->Get("cOutput");
  TList *l_raw = (TList*)f_raw->Get("cOutput");

  //Corrections
  //Background
  TH2F *hBackgroundCorr = new TH2F("backgroundCorr","Background correction",60,-3.,3.,20,-20.,20.);
  hBackgroundCorr->Sumw2();
  TH2F *hBackgroundCorrNum = (TH2F*)l_corr->At(31);

  Int_t corrHistoName;
  if (kPrimariesTest) {
    corrHistoName = 31;
  } else {
    corrHistoName = 32;
  }

  TH2F *hBackgroundCorrDen = (TH2F*)l_corr->At(corrHistoName);
  hBackgroundCorr->Divide(hBackgroundCorrNum,hBackgroundCorrDen,1.,1.); 

  //Algorithm efficiency
  TH2F *hTrackletRecEff = new TH2F("trackletRecEffCorr","Tracklet algorithm efficiency",60,-3,3,20,-20,20);
  hTrackletRecEff->Sumw2();
  TH2F *hTrackletEffCorr_num = (TH2F*)l_corr->At(33);  
  hTrackletRecEff->Divide(hBackgroundCorrNum,hTrackletEffCorr_num,1.,1.);   // 1/Efficiency=correction_weight

  //Correction for particles that do not reach the sensitive layer
  TH2F *hDecPartLay2Corr = new TH2F("decPartLay2Corr","",60,-3.,3.,20,-20.,20.);
  hDecPartLay2Corr->Sumw2();
  TH2F *hPrimaries_evtTrVtx = (TH2F*) l_corr->At(34);
  TH2F *hDetectablePrimariesLay2_evtTrVtx = (TH2F*) l_corr->At(35);
  hDecPartLay2Corr->Divide(hPrimaries_evtTrVtx,hDetectablePrimariesLay2_evtTrVtx,1.,1.);

  //Acceptance
  TH2F *hAccCorr = new TH2F("accEffCorr_test","",60,-3,3,20,-20,20);
  hAccCorr->Sumw2();
  hAccCorr->Divide(hTrackletEffCorr_num,hDetectablePrimariesLay2_evtTrVtx,1.,1.);

  TH2F *hAccCorrPb = (TH2F*) f_accPb->Get("AccTrac");

  //Vertex-trigger correction
  TH2F *hVertexTrigCorrEtaNum = (TH2F*)l_corr->At(36);
  TH2F *hTriggerCorrEta_norm  = (TH2F*)l_corr->At(37);
  TH2F *hVertexCorrEta_norm   = (TH2F*)l_corr->At(34);   

  TH2F *hVertexCorrEta = new TH2F("vertexCorrEta","",60,-3.,3.,20,-20.,20.);
  hVertexCorrEta->Sumw2();
  hVertexCorrEta->Divide(hTriggerCorrEta_norm,hVertexCorrEta_norm,1.,1.);

  TH2F *hTriggerCorrEta = new TH2F("triggerEtaCorr","",60,-3.,3.,20,-20.,20.);
  hTriggerCorrEta->Sumw2();
  hTriggerCorrEta->Divide(hVertexTrigCorrEtaNum,hTriggerCorrEta_norm,1.,1.);

  TH2F *hVertexTrigCorrEta = new TH2F("vertexTrigEtaCorr","",60,-3.,3.,20,-20.,20.);
  hVertexTrigCorrEta->Sumw2();
  hVertexTrigCorrEta->Divide(hVertexTrigCorrEtaNum,hVertexCorrEta_norm,1.,1.);

  TH2F *hVertexTrigCorrNum = (TH2F*)l_corr->At(38);
  TH2F *hTriggerCorr_norm  = (TH2F*)l_corr->At(40);
  TH2F *hVertexCorr_norm   = (TH2F*)l_corr->At(39);

  TH2F *hVertexCorr = new TH2F("vertexCorr","",200,0.,200.,20,-20.,20.); 
  hVertexCorr->Sumw2();
  hVertexCorr->Divide(hTriggerCorr_norm,hVertexCorr_norm,1.,1.);

  TH2F *hTriggerCorr = new TH2F("triggerCorr","",200,0.,200.,20,-20.,20.);
  hTriggerCorr->Sumw2();
  hTriggerCorr->Divide(hVertexTrigCorrNum,hTriggerCorr_norm,1.,1.);

  TH2F *hVertexTrigCorr = new TH2F("vertexTrigCorr","",200,0.,200.,20,-20.,20.);
  hVertexTrigCorr->Sumw2();
  hVertexTrigCorr->Divide(hVertexTrigCorrNum,hVertexCorr_norm,1.,1.);

  TH2F *hVertexMC2D = (TH2F*)l_MC->At(38);   
  TH1D *hMCvertex = new TH1D(*(hVertexMC2D->ProjectionY("MCvertex",0,-1,"e")));  //MC vertex distrib all events
//  hMCvertex->Draw(); 

  //Raw distributions

  //...track level
  TH2F *hSPDEta_2Draw_accBin = new TH2F(*((TH2F*)l_raw->At(2)));  
  hSPDEta_2Draw_accBin->SetNameTitle("SPDEta_2Draw_accBin","");
  TH2F *hSPDEta_2Draw = new TH2F(*((TH2F*)l_raw->At(2)));
  hSPDEta_2Draw->SetNameTitle("SPDEta_2Draw","");

  //...event level
  TH2F *hESDVertexvsNtracklets = (TH2F*)l_raw->At(0);
  TH2F *hESDVertexvsNtracklets_trigEvts = (TH2F*)l_raw->At(1);
  
  TH2F *hESDVertexvsNtracklets_Corr = new TH2F("ESDVertexvsNtracklets_Corr","",200,0.,200.,20,-20.,20.);
  
  TH1D *hSPDvertex = (TH1D*)l_raw->At(18);
  if (kPrimariesTest) {
    hSPDEta_2Draw = new TH2F(*((TH2F*)l_MC->At(31)));
    hSPDEta_2Draw->SetNameTitle("SPDEta_2Draw","");
    hESDVertexvsNtracklets = new TH2F(*((TH2F*)l_MC->At(39)));
    hESDVertexvsNtracklets_trigEvts = new TH2F(*((TH2F*)l_MC->At(40)));

    hSPDvertex = hVertexCorr_norm->ProjectionY("MCvertex_RV",0,-1,"e");  
  }

  if (!kPrimariesTest) hSPDEta_2Draw->Rebin2D(2,2,"");

  //Vertex correction

  //Computing the percentage of triggered events with vertex and tracklet mult!=0 for each rec vertex bin
  TH1F* fHistZdistribTrigVtxEvt = new TH1F("fHistZdistribTrigVtxEvt", "",20,-20.,20.);
  fHistZdistribTrigVtxEvt->Sumw2();
  for (Int_t ivtx=binVtxStart; ivtx<binVtxStop; ivtx++) {
     
    Float_t fracZ = (hESDVertexvsNtracklets->ProjectionY("_y",2,200,"e")->GetBinContent(ivtx+1))/Float_t (hESDVertexvsNtracklets->ProjectionY("_y",2,200,"e")->Integral(binVtxStart+1,binVtxStop)); 
    fHistZdistribTrigVtxEvt->SetBinContent(ivtx+1,fracZ); 
  }

  //Computing the factor to correct the previous % 
  TH1F* fHistEvtCorrFactor = new TH1F("fHistEvtCorrFactor", "",20,-20.,20.);
  fHistEvtCorrFactor->Sumw2();
  for (Int_t ivtx=0; ivtx<20; ivtx++) {
    Float_t fracZnum = (hTriggerCorr_norm->GetBinContent(1,ivtx+1))/Float_t(hTriggerCorr_norm->ProjectionY("_y",1,1,"e")->Integral());
    Float_t fracZden = (hVertexCorr_norm->ProjectionY("_y",2,200,"e")->GetBinContent(ivtx+1))/Float_t(hVertexCorr_norm->ProjectionY("_y",2,200,"e")->Integral()); 
    Float_t ratio =0.;
    if (fracZden!=0.) ratio = fracZnum/fracZden;
    if(ratio!=0.) fHistEvtCorrFactor->SetBinContent(ivtx+1,ratio);
  }

  TH1F* hprod = new TH1F("prod", "",20,-20.,20.);
  hprod->Sumw2();

  hprod->Multiply(fHistEvtCorrFactor,fHistZdistribTrigVtxEvt,1.,1.);

  //MC distributions

  Int_t histoNameRV;
  if (kPrimariesTest) {
    histoNameRV = 41;
  } else {
    histoNameRV = 42;   
  }

  TH2F *Eta_RV = (TH2F*)l_MC->At(histoNameRV);
  TH2F *Eta = (TH2F*)l_MC->At(43);
//  Eta->Draw();

  TH1D *hdNdEta_RV = new TH1D("dNdEta_RV","Pseudorapidity ",60,-3,3);
  hdNdEta_RV=Eta_RV->ProjectionX("dNdEta_RV",binVtxStart+1,binVtxStop,"e");
  Double_t nEvtsRec=0;
  for (Int_t ivtx=binVtxStart; ivtx<binVtxStop; ivtx++) {
    nEvtsRec += hSPDvertex->GetBinContent(ivtx+1);
  }
  hdNdEta_RV->Scale(10./nEvtsRec);

  TH1D *hdNdEta_Trigg = new TH1D("dNdEta_Trigg","Pseudorapidity ",60,-3,3);
  TH2F *hEta_Trigg = new TH2F(*((TH2F*)l_corr->At(37))); 
  hdNdEta_Trigg=hEta_Trigg->ProjectionX("dNdEta_Trigg",0,-1,"e");
  Double_t nEvtsTrigg=0;
  nEvtsTrigg =hTriggerCorr_norm->Integral();
  hdNdEta_Trigg->Scale(10./nEvtsTrigg);

  TH1D *hdNdEta = new TH1D("dNdEta","Pseudorapidity ",60,-3.,3.);
  hdNdEta->Sumw2();
  hdNdEta = Eta->ProjectionX("dNdEta",0,-1,"e");
  Double_t nEvtsTot=0;
  for (Int_t ivtx=0; ivtx<20; ivtx++) {
    nEvtsTot += hMCvertex->GetBinContent(ivtx+1);
  }
  cout<<"Number of generated events: "<<nEvtsTot <<endl;

  hdNdEta->Scale(10./nEvtsTot);
  hdNdEta->GetXaxis()->SetTitle("Pseudorapidity #eta");
  hdNdEta->GetYaxis()->SetTitle("dN_{ch}/d#eta");
//  hdNdEta->Draw();

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

  TH1F *hnorm_bkgCorr =  new TH1F("norm_bkgCorr", "",60,-3.,3.);
  hnorm_bkgCorr->Sumw2();
  TH1F *hnorm_bkgEffCorr =  new TH1F("norm_bkgEffCorr", "",60,-3.,3.);
  hnorm_bkgEffCorr->Sumw2();
  TH1F *hnorm_bkgEffAccCorr =  new TH1F("norm_bkgEffAccCorr", "",60,-3.,3.);
  hnorm_bkgEffAccCorr->Sumw2();
  TH1F *hnorm_fullyCorr =  new TH1F("norm_fullyCorr", "",60,-3.,3.);
  hnorm_fullyCorr->Sumw2();
  TH1F *hnorm_triggered =  new TH1F("norm_triggered", "",60,-3.,3.);
  hnorm_triggered->Sumw2();
  TH1F *hnorm =  new TH1F("norm", "",60,-3.,3.);
  hnorm->Sumw2();
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

  //Cutting either on weights or on relative errors of weights
  // Acceptance Pb-Pb
  TH2F *hAccErrorsvsAccTr = new TH2F("AccErrorsvsAccTr","",20000,0.,10000.,50,0.,1.);
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
      //Acceptance p-p
      relBinErr = hAccCorr->GetBinError(jeta+1,kvtx+1)/hAccCorr->GetBinContent(jeta+1,kvtx+1);
      if (relBinErr>cutAcc) {
//      if (hAccCorr->GetBinContent(jeta+1,kvtx+1)<cutAcc) {
        hAccCorr->SetBinContent(jeta+1,kvtx+1,0.);
        hAccCorr->SetBinError(jeta+1,kvtx+1,0.);
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
    }
  }

  for (Int_t ieta=0; ieta<60; ieta++) {
    for (Int_t izeta=0; izeta<20; izeta++) {
      if (hAccCorr->GetBinContent(ieta+1,izeta+1)==0)
        hDecPartLay2Corr->SetBinContent(ieta+1,izeta+1,0.);
    }
  }

  for (Int_t ieta=0; ieta<60; ieta++) {
    for (Int_t izeta=0; izeta<20; izeta++) {
      if (hAccCorr->GetBinContent(ieta+1,izeta+1)==0)
        hVertexTrigCorrEta->SetBinContent(ieta+1,izeta+1,0.);
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
    hSPDEta_2D_bkgEffAccCorr->Add(hSPDEta_2Draw);
    hSPDEta_2D_bkgEffAccCorr->Divide(hAccCorr);
  }
  hSPDEta_2D_bkgEffAccCorr->Multiply(hBackgroundCorr);
  hSPDEta_2D_bkgEffAccCorr->Divide(hTrackletRecEff);
  hSPDEta_2D_fullyCorr->Add(hSPDEta_2D_bkgEffAccCorr);
  hSPDEta_2D_fullyCorr->Multiply(hDecPartLay2Corr);
  hSPDEta_2D_triggeredEvts->Add(hSPDEta_2D_fullyCorr);
  hSPDEta_2D_triggeredEvts->Multiply(hVertexCorrEta);
  hSPDEta_2D_allEvts->Add(hSPDEta_2D_triggeredEvts);
  hSPDEta_2D_allEvts->Multiply(hTriggerCorrEta);

  //...event level
  Double_t evtTrigNoTracklets = hESDVertexvsNtracklets_trigEvts->ProjectionY("_y",1,1,"e")->Integral(); 
  for (Int_t ivtx=0;ivtx<20;ivtx++) { 
    hESDVertexvsNtracklets_trigEvts->SetBinContent(1,ivtx+1,(hprod->GetBinContent(ivtx+1))*(evtTrigNoTracklets));
  }

  hESDVertexvsNtracklets_Corr->Add(hESDVertexvsNtracklets_trigEvts);
  hESDVertexvsNtracklets_Corr->Multiply(hTriggerCorr);

  //Filling normalization histograms

  TH1D *hSPDvertexCorr_y = new TH1D("SPDvertexCorr_y","",20,-20,20);
  hSPDvertexCorr_y = hESDVertexvsNtracklets_Corr->ProjectionY("SPDvertexCorr_y",0,-1,"e");
  cout<<"Number of events (tracklets corrected)"<<hSPDvertexCorr_y->Integral(binVtxStart+1,binVtxStop)<<endl;  
  for (Int_t ivtx=binVtxStart; ivtx<binVtxStop; ivtx++) {
    Double_t nEvtsivtx = hSPDvertex->GetBinContent(ivtx+1);
    Double_t nEvtsivtxCorr = hSPDvertexCorr_y->GetBinContent(ivtx+1);
    Double_t nEvtsivtxCorr_triggered = hESDVertexvsNtracklets_trigEvts->ProjectionY("_y",0,-1,"e")->GetBinContent(ivtx+1);
    Double_t nEvtGen = hMCvertex->GetBinContent(ivtx+1);
    cout<<"iBinVtx: "<<ivtx<<" nEvtsGen-> "<<nEvtGen<<" nEvtsCorr-> "<<nEvtsivtxCorr<<endl;
    for (Int_t jeta=0; jeta<60; jeta++) {
      if (hSPDEta_2D_bkgCorr->GetBinContent(jeta+1,ivtx+1)!=0.)
        hnorm_bkgCorr->Fill(hBackgroundCorr->GetXaxis()->GetBinCenter(jeta+1),nEvtsivtx);
      if (hSPDEta_2D_bkgEffCorr->GetBinContent(jeta+1,ivtx+1)!=0.)
        hnorm_bkgEffCorr->Fill(hTrackletRecEff->GetXaxis()->GetBinCenter(jeta+1),nEvtsivtx);
      if (hSPDEta_2D_bkgEffAccCorr->GetBinContent(jeta+1,ivtx+1)!=0.)
        hnorm_bkgEffAccCorr->Fill(hBackgroundCorr->GetXaxis()->GetBinCenter(jeta+1),nEvtsivtx);
      if (hSPDEta_2D_fullyCorr->GetBinContent(jeta+1,ivtx+1)!=0.)
        hnorm_fullyCorr->Fill(hBackgroundCorr->GetXaxis()->GetBinCenter(jeta+1),nEvtsivtx);
      if (hSPDEta_2D_triggeredEvts->GetBinContent(jeta+1,ivtx+1)!=0.)
        hnorm_triggered->Fill(hBackgroundCorr->GetXaxis()->GetBinCenter(jeta+1),nEvtsivtxCorr_triggered);
      if (hSPDEta_2D_allEvts->GetBinContent(jeta+1,ivtx+1)!=0.)
        hnorm->Fill(hBackgroundCorr->GetXaxis()->GetBinCenter(jeta+1),nEvtsivtxCorr);
    }
  }
  
  for (Int_t jeta=0; jeta<60; jeta++) {
    hnorm_bkgCorr->SetBinError(jeta+1,TMath::Sqrt(hnorm_bkgCorr->GetBinContent(jeta+1)));
    hnorm_bkgEffCorr->SetBinError(jeta+1,TMath::Sqrt(hnorm_bkgEffCorr->GetBinContent(jeta+1)));
    hnorm_bkgEffAccCorr->SetBinError(jeta+1,TMath::Sqrt(hnorm_bkgEffAccCorr->GetBinContent(jeta+1)));
    hnorm_fullyCorr->SetBinError(jeta+1,TMath::Sqrt(hnorm_fullyCorr->GetBinContent(jeta+1)));
    hnorm_triggered->SetBinError(jeta+1,TMath::Sqrt(hnorm_triggered->GetBinContent(jeta+1)));
    hnorm->SetBinError(jeta+1,TMath::Sqrt(hnorm->GetBinContent(jeta+1)));
  }

  hdNdEta_bkgCorr->Divide(hSPDEta_2D_bkgCorr->ProjectionX("SPDEta_2D_bkgCorr_x",binVtxStart+1,binVtxStop,"e"),hnorm_bkgCorr);
  hdNdEta_bkgEffCorr->Divide(hSPDEta_2D_bkgEffCorr->ProjectionX("SPDEta_2D_bkgEffCorr_x",binVtxStart+1,binVtxStop,"e"),hnorm_bkgEffCorr);
  hdNdEta_bkgEffAccCorr->Divide(hSPDEta_2D_bkgEffAccCorr->ProjectionX("SPDEta_2D_bkgEffAccCorr_x",binVtxStart+1,binVtxStop,"e"),hnorm_bkgEffAccCorr);
  hdNdEta_fullyCorr->Divide(hSPDEta_2D_fullyCorr->ProjectionX("SPDEta_2D_fullyCorr_x",binVtxStart+1,binVtxStop,"e"),hnorm_fullyCorr);
  hdNdEta_triggeredEvts->Divide(hSPDEta_2D_triggeredEvts->ProjectionX("SPDEta_2D_triggeredEvts_x",binVtxStart+1,binVtxStop,"e"),hnorm_triggered);
  hdNdEta_allEvts->Divide(hSPDEta_2D_allEvts->ProjectionX("SPDEta_2D_allEvts_x",binVtxStart+1,binVtxStop,"e"),hnorm);

  hdNdEta_bkgCorr->Scale(10.);
  hdNdEta_bkgEffCorr->Scale(10.);
  hdNdEta_bkgEffAccCorr->Scale(10.);
  hdNdEta_fullyCorr->Scale(10.);
  hdNdEta_triggeredEvts->Scale(10.);
  hdNdEta_allEvts->Scale(10.);

  TH2F* hRatiodNdEta2D=new TH2F("ratiodNdEta2D","",60,-3,3,20,-20,20);
//  hRatiodNdEta2D->Divide(hSPDEta_2D_fullyCorr,Eta_RV,1.,1.);
  hRatiodNdEta2D->Divide(hSPDEta_2D_allEvts,Eta,1.,1.);

  TH1F* hRationorm=new TH1F("rationorm","",60,-3,3);
  hRationorm->Add(hnorm);
  hRationorm->Scale(1./nEvtsTot);

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
    hdNdEta_test->Divide(Eta_test->ProjectionX("eta_x",0,-1,"e"),hnorm,1.,1.);
    hdNdEta_test->Scale(10.);
  }

  //Ratios generated/corrected distributions
  TH1F* hRatiodNdEta=new TH1F("ratiodNdEta","",60,-3,3);
  hRatiodNdEta->Sumw2(); 

  hRatiodNdEta->SetMarkerColor(2);
  hRatiodNdEta->SetMarkerStyle(22);
  hRatiodNdEta->GetXaxis()->SetTitle("Pseudorapidity #eta");
  hRatiodNdEta->GetYaxis()->SetTitle("Generated/Corrected");

  TH1F* hRatiodNdEta_vertex=new TH1F("ratiodNdEta_vertex","",60,-3,3);
  TH1F* hRatiodNdEta_triggered=new TH1F("ratiodNdEta_triggered","",60,-3,3);


  for (Int_t i=0;i<60;i++) hdNdEta->SetBinError(i+1,0.);
  for (Int_t i=0;i<60;i++) hdNdEta_RV->SetBinError(i+1,0.);

  if (kPrimariesTest) {
    hRatiodNdEta->Divide(hdNdEta_test,hdNdEta_allEvts,1.,1.); 
    
  } else {
    hRatiodNdEta->Divide(hdNdEta,hdNdEta_allEvts,1.,1.);
    hRatiodNdEta_vertex->Divide(hdNdEta_RV,hdNdEta_fullyCorr,1.,1.);
    hRatiodNdEta_triggered->Divide(hdNdEta_Trigg,hdNdEta_triggeredEvts,1.,1.);
  }

  // Draw dN/dEta
  hdNdEta_raw->SetLineColor(kGreen);
  hdNdEta_raw->SetLineWidth(3);

  hdNdEta_RV->SetLineWidth(3);
  hdNdEta_RV->SetMinimum(0.);
  hdNdEta_RV->SetMaximum(15.);
  hdNdEta_RV->SetLineColor(kRed);
//  hdNdEta_RV->Draw("histo");
  hdNdEta->SetLineWidth(3);
  hdNdEta->SetMinimum(0.);
  hdNdEta->SetMaximum(15.);
  hdNdEta->Draw("histo");
//  hdNdEta_RV->Draw("histo,same");

  hdNdEta_raw->Draw("histo,same");

  hdNdEta_bkgCorr->SetMarkerStyle(21);
  hdNdEta_bkgCorr->SetMarkerColor(kBlue);

  hdNdEta_bkgCorr->Draw("same,e");

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

  TLegend *leg1 = new TLegend(0.3,0.7,0.7,0.9,NULL,"brNDC");
  leg1->SetTextFont(62);
  leg1->SetTextSize(0.0304569);
  TLegendEntry *entry1=leg1->AddEntry(hdNdEta_raw,"Reconstructed","l");
//                entry1=leg1->AddEntry(hdNdEta_RV,"Generated","l");
                entry1=leg1->AddEntry(hdNdEta,"Generated","l");
                entry1=leg1->AddEntry(hdNdEta_bkgCorr,"Bkg corrected","p");
                entry1=leg1->AddEntry(hdNdEta_bkgEffCorr,"Bkg+AlgEff+SPDEff corrected","p");
                entry1=leg1->AddEntry(hdNdEta_bkgEffAccCorr,"Bkg+AlgEff+SPDAccEff corrected","p");
                entry1=leg1->AddEntry(hdNdEta_fullyCorr,"Corrected (triggered events with vertex)","p");
                entry1=leg1->AddEntry(hdNdEta_triggeredEvts,"Corrected (triggered events)","p");
                entry1=leg1->AddEntry(hdNdEta_allEvts,"Corrected (all events)","p");

  leg1->Draw();

//  hRatiodNdEta->Draw();

  // Write histos

  TFile *fout= new TFile("SPDdNdEtaCorr.root","RECREATE"); 

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

  hnorm_bkgCorr->Write();
  hnorm_bkgEffCorr->Write();
  hnorm_bkgEffAccCorr->Write();
  hnorm_triggered->Write();
  hnorm->Write();
  hnorm_gen->Write();

  hVertexCorr->Write();
  hVertexTrigCorr->Write();
  hTriggerCorr->Write();
  hVertexCorrEta->Write();
  hTriggerCorrEta->Write();
  hVertexTrigCorrEta->Write();

  fHistZdistribTrigVtxEvt->Write();
  fHistEvtCorrFactor->Write();

  hESDVertexvsNtracklets->Write();
  hESDVertexvsNtracklets_trigEvts->Write();
  hESDVertexvsNtracklets_Corr->Write();
 
  hdNdEta->Write();
  hdNdEta_RV->Write();
  hdNdEta_bkgCorr->Write();
  hdNdEta_bkgEffCorr->Write();
  hdNdEta_bkgEffAccCorr->Write();
  hdNdEta_fullyCorr->Write();
  hdNdEta_triggeredEvts->Write();
  hdNdEta_allEvts->Write();
  hdNdEta_raw->Write();

  hAccErrorsvsAccTr->Write();

  hRatiodNdEta->Write();
  hRatiodNdEta_triggered->Write();
  hRatiodNdEta_vertex->Write();
  hRatiodNdEta2D->Write();
  hRationorm->Write();

  Eta_test->Write();
  Eta->Write(); 

  hprod->Write();
  
  fout->Write();
  fout->Close();

}
