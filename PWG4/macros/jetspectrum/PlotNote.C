#include "TList.h"
#include "TColor.h"
#include "TArrow.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "THnSparse.h"
#include "TDirectoryFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TMath.h"
#include "TF1.h"
#include "TASImage.h"
#include "TPaveText.h"
#include "TParameter.h"
#include "TPaveStats.h"


#include "AliAnalysisHelperJetTasks.h"
#include "AliAnalysisTaskJetServices.h"
#include "AliAnalysisTaskJetSpectrum2.h"
#include <fstream>
#include <iostream>

#include "Normalize2d.C"
#include "ALICEWorkInProgress.C"

using namespace std;

Int_t PlotSpectraPbPb();
Int_t PlotSpectraPbPb2();
Int_t PlotJetBFluctuations();
Int_t PlotJetBFluctuations2(UInt_t iPlotFlag = 0xFFFFF,UInt_t iPlotType = 0xFFFFF,Float_t inPtLo = 50,Float_t inPtUp = 100);
Int_t PlotSubtract();
Int_t PlotJetQA();

// helpers

Int_t fDebug = 0;

enum {kMaxJets = 4};
void HistsFromSingleJets(char* cMask,TList *list,Int_t nJets,TH1F** h,Int_t iRebin = 1,int iFirst = 0);
void HistsFromSingleJets(char* cMask,TList *list,Int_t nJets,TH2F** h,Int_t iRebin = 1,int iFirst = 0);
void ScaleNevent(TH1 *h,TFile *fIn,Int_t iCond = 5,Int_t iMC = 0,Int_t iCen = 0);
void ScaleH1(TH1* h,char* cX = "",char* cY = "",Float_t fScale = 1,Bool_t bWidth = true);
void SetStyleH1(TH1 *h,Int_t color = 0,Int_t fillStyle = 0,Int_t markerStyle = 0);
void ScaleH2(TH2* h,char* cX = "",char* cY = "",char* cZ = "",Float_t fScale = 1,Bool_t bWidth = true);
TH1F* GetPythia8Spectrum(const char *cHist,const char *cPythia,Float_t deta,Int_t iRebin = 1,Int_t iMask = 0);
Double_t GetPoissonFluctuation(TH1 *h1,Double_t areaIn,Double_t areaJet);
TObject* GetObjectFromFile(const char *cName,const char *cFile,const char* cDir,const char* cRep = "pwg4");
TF1* FitLHSgaus(TH1D *hDeltaPt, double minPt = -60., double maxPt = 5., int minIterations = 1, int maxIterations = 1, double maxDeltaMean = 0.01, int mode=1, double minBoundSigma = 3., double maxBoundSigma = 0.5 );
Double_t Gamma(Double_t *x,Double_t *par);
Double_t Gamma0(Double_t *x,Double_t *par);
TH2F* GetTH2PlotB(const char *cPath,Int_t embType=0, Int_t classType=0, Int_t cl=-1, Int_t rp=-1);
void SetHistoAttributes(TH1* h1,Int_t iMarker = kFullCircle,Int_t iColor = kBlack);
void SetGraphAttributes(TGraph* gr,Int_t iMarker = kFullCircle,Int_t iColor = kBlack);

void CloseFiles();

const char *cPrintMask = "110116_%s.eps";
void set_plot_style();
TCanvas *cTmp = 0;

TString picPrefix = "110504_";
TString picSuffix = "png";
TString picName = "";

TH1F* GetPythia8Spectrum(const char *cHist,const char *cPythia,Float_t deta,Int_t iRebin,Int_t iMask){
  TDirectory *opwd = gDirectory;
  TFile *fPythia = (TFile*)gROOT->GetListOfFiles()->FindObject(cPythia);
  
  if(!fPythia) fPythia  = TFile::Open(cPythia);
  opwd->cd();

  static int iCount  = 0;
  TH1F* hist = 0;
  if(iMask<=0){
    hist = (TH1F*)fPythia->Get(cHist);      
    if(hist)hist = (TH1F*)hist->Clone(Form("%s_%d",hist->GetName(),iCount));
  }
  else{
    for(int i = 0;i < iMask;i++){
      TH1F *hTmp = (TH1F*)fPythia->Get(Form(cHist,i));
      if(!hTmp){
	Printf("%s not found",Form(cHist,i));
	return 0;
      }   
      if(!hist)hist = (TH1F*)hTmp->Clone(Form("%s_%d_%d",hTmp->GetName(),iMask,iCount));
      else hist->Add(hTmp);
    }
  }

  if(!hist){
    Printf("%s not found",cHist);
    return 0;
  }  
  // fetch the cross section
  TH1F *hxsec = (TH1F*)fPythia->Get("h1Xsec");

  Float_t xsec =  hxsec->GetBinContent(1);
  Printf("%d xsec = %1.3E",__LINE__,xsec);
  //  if(xsec==0)xsec = 40.79; // tmp hack
  hist->Rebin(iRebin);
  hist->Scale(1./hist->GetBinWidth(1)); 
  //  hist->Scale(1./xsec);
  hist->Scale(1./deta);

  //  Double_t xtotal = 71.39; // xtotal for 7 TeV!!
  Double_t xtotal = 62.01; // xtotal for 7 TeV!!

  // scale with fraction of total NSD cross section
  hist->Scale(1/xtotal);
  Printf("%d fraction = %1.3E of total cross section %1.3E",__LINE__,xsec/xtotal,xtotal);
  hist->SetXTitle("p_{T} (GeV/c)");
  hist->SetYTitle("1/N  dN/dp_{T}dy");
  hist->SetDirectory(opwd);
  iCount++;
  return hist;

}




void PlotNote(){
  set_plot_style();


  UInt_t iPlotFlag = 0xFFF;  
  // plot only RC
  //  iPlotFlag = (1<<0)|(1<<1)|(1<<2);//RC
  //  iPlotFlag = (1<<0);//RC
  //  iPlotFlag = (1<<3)|(1<<5)|(1<<6);//jets
  //  iPlotFlag = (1<<4);//tracks
  //  iPlotFlag = (1<<7);//tracks and poission
  //  iPlotFlag = (1<<0)|(1<<1)|(1<<3)|(1<<4);//similars
  //  PlotJetBFluctuations2(iPlotFlag,1<<0);
  // plot only jets
  // iPlotFlag = (1<<3)|(1<<5)|(1<<6);
  // plot only tracks
  //

  // QM Plots

  // All Random Cones NO RP
  UInt_t iEst = 1<<7|1<<8;
  /*
  iPlotFlag = (1<<0)|(1<<1)|(1<<2)|iEst;//RC
  PlotJetBFluctuations2(iPlotFlag,1<<0);
  iPlotFlag = (1<<3)|(1<<5)|(1<<6)|iEst;//jets
  PlotJetBFluctuations2(iPlotFlag,1<<0);
  iPlotFlag = (1<<4)|iEst;//single tracks
  PlotJetBFluctuations2(iPlotFlag,1<<0);
  // the essentials 
  
  iPlotFlag = (1<<0)|(1<<3)|(1<<4)|iEst;//jets
  PlotJetBFluctuations2(iPlotFlag,1<<0);


  // vs mult
  
  iPlotFlag = (1<<0)|(1<<1)|(1<<4)|iEst;//jets
  PlotJetBFluctuations2(iPlotFlag,1<<1);
  */



  // vs cent wit RP
  iPlotFlag = (1<<4)|iEst;//jets
  PlotJetBFluctuations2(iPlotFlag,1<<0|1<<2);

  iPlotFlag = (1<<0)|iEst;//jets
  PlotJetBFluctuations2(iPlotFlag,1<<0|1<<2);

  return;
  // vs mult with RP
  iPlotFlag = (1<<4)|iEst;//jets
  PlotJetBFluctuations2(iPlotFlag,1<<1|1<<2);

  // vs 

  //  PlotSpectraPbPb2();
  //  PlotSubtract();
}

TList *GetList(const char *cFile,const char *cName){
  TDirectory *opwd = gDirectory;
  TFile *fIn = (TFile*)gROOT->GetListOfFiles()->FindObject(cFile);
  TList *list = 0;
  if(!fIn) fIn  = TFile::Open(cFile);
  if(!fIn)return list;
  opwd->cd();

  TDirectory *dir = (TDirectory*)fIn->Get(Form("PWG4_%s",cName));
  if(!dir){
    Printf("GetList: Directory PWG4_%s not found",cName);
    return list;
  }				  
  list = (TList*)dir->Get(Form("pwg4%s",cName));
  if(!list){
    Printf("GetList: list pwg4%s not found",cName);
    return list;
  }				  
  return list;

}

Int_t PlotJetQA(){

  set_plot_style();
  gStyle->SetTitleYSize(0.7*gStyle->GetTitleYSize());
  gStyle->SetTitleXSize(0.7*gStyle->GetTitleXSize());
  gStyle->SetTitleXOffset(1.2*gStyle->GetTitleXOffset());
  gStyle->SetTitleYOffset(1.5*gStyle->GetTitleYOffset());
  gStyle->SetPadRightMargin(1.1*gStyle->GetPadRightMargin());
  gStyle->SetPadBottomMargin(0.7*gStyle->GetPadBottomMargin());
  

  Bool_t isPP = false;

  TString printType = "eps"; 
  const Int_t nCen = 5;
  TString sCen[nCen] = {" 0-80"," 0-10%%","10-30","30-50","50-80"};
  const Int_t nJets = 3;


  TString cAdd = "Rec"; // decide whcih type to take Rec RecFull Gen GenFull
  TString cAddFull = "RecFull"; // decide whcih type to take Rec RecFull Gen GenFull

  TCanvas *c1 = new TCanvas();

  TString sFile =  "~/alice/data/analysis/train_maf/output_110216.root";
  TString sinputName  = "spec2_clustersAOD_ANTIKT04_B2_Filter00256_Cut00150_Skip00_clustersAOD_ANTIKT04_B1_Filter00256_Cut00150_Skip00_0000000000";
  TString sinputJet = "%% Pb+Pb Anti-k_{T} R = 0.4 (B2)";
  TString sinputPrint = "PbPb_antikt04_B2";
  sinputName  = "spec2_clustersAOD_ANTIKT04_B1_Filter00256_Cut00150_Skip00_clustersAOD_ANTIKT04_B0_Filter00256_Cut00150_Skip00_0000000000";cAdd = "Gen";cAddFull = "GenFull";sinputJet = "%% Pb+Pb Anti-k_{T} R = 0.4 (B0)";sinputPrint = "PbPb_antikt04_B0";

  if(isPP){
    sinputName  = "spec2_clustersAOD_ANTIKT04_B0_Filter00256_Cut00150_Skip00_clustersAOD_ANTIKT04_B0_Filter00256_Cut00150_Skip00_0000000000";sinputJet = "p+p #sqrt{s} = 7 TeV Anti-k_{T} R = 0.4 (B0)";sinputPrint = "pp_antikt04_B0";
    sCen[0] = "";
    sFile = "~/alice/jets/macros/batch/output.root";
  }


  TFile *fIn = TFile::Open(sFile.Data());


  TDirectory *dir = (TDirectory*)fIn->Get("PWG4_services");
  TList *list1 = (TList*)dir->Get("pwg4serv");

  TH2F* hTriggerCount =(TH2F*)list1->FindObject("fh2TriggerCount");

  for(int ic = 0;ic < (isPP?1:nCen);ic++){

    TString cName = Form("PWG4_%s_Class%02d",sinputName.Data(),ic);
    TDirectory *inputDir = (TDirectory*)fIn->Get(cName.Data());

    Float_t fNevents = hTriggerCount->GetBinContent(hTriggerCount->FindBin(ic,AliAnalysisTaskJetServices::kSelected));

    Printf(" %s: %10d events",sCen[ic].Data(),(Int_t)fNevents); 

    if(!inputDir){
      Printf("Dir %s not found",cName.Data());
      continue;
    }
    cName = Form("pwg4%s_Class%02d",sinputName.Data(),ic);
    TList *inputList = (TList*)inputDir->Get(cName.Data());
    if(!inputList){
      Printf("List %s not found",cName.Data());
      continue;
    }

    TH2F *h2EtaPtLead[nJets+2] = {0,};
    TH2F *h2EtaPtLeadFull[nJets+2] = {0,};
    TH2F *h2EtaPtAll[nJets+2] = {0,};
    HistsFromSingleJets(Form("fh2EtaPt%s_j%%d",cAdd.Data()),inputList,1,h2EtaPtLead);
    HistsFromSingleJets(Form("fh2EtaPt%s_j%%d",cAddFull.Data()),inputList,1,h2EtaPtLeadFull);
    HistsFromSingleJets(Form("fh2EtaPt%s_j%%d",cAddFull.Data()),inputList,nJets+1,h2EtaPtAll,1,nJets);

    c1->SetLogz();   
    for(int i = 0; i< nJets+2;i++)Printf("Lead     > %d %p",i,h2EtaPtLead[i]);
    for(int i = 0; i< nJets+2;i++)Printf("LeadFull > %d %p",i,h2EtaPtLeadFull[i]);
    for(int i = 0; i< nJets+2;i++)Printf("All      > %d %p",i,h2EtaPtAll[i]);
    h2EtaPtLead[0]->Draw("colz");    
    c1->Modified();
    c1->Update();
    c1->SaveAs(Form("%s_%s_Cen%02d_EtaPtLead2D.%s",sinputPrint.Data(),cAdd.Data(),ic,printType.Data()));
    if(!gROOT->IsBatch())
      if(getchar()=='q')return 1;

    Float_t fStep = 20;
    Float_t fStepLo = 0;
    Float_t fStepUp = fStepLo+fStep;
    Float_t fMax = 160;
    
    while(fStepUp<fMax){
      Int_t iBinLo = h2EtaPtLead[0]->GetYaxis()->FindBin(fStepLo);
      Int_t iBinUp = h2EtaPtLead[0]->GetYaxis()->FindBin(fStepUp)-1;
      TH1D *hProjLead = h2EtaPtLead[0]->ProjectionX(Form("hEtaPtLead_%d_%d_%d",iBinLo,iBinUp,ic),iBinLo,iBinUp,"e");
      TH1D *hProjLeadFull = h2EtaPtLeadFull[0]->ProjectionX(Form("hEtaPtLeadFull_%d_%d_%d",iBinLo,iBinUp,ic),iBinLo,iBinUp,"e");
      TH1D *hProjAll = h2EtaPtAll[nJets+1]->ProjectionX(Form("hEtaPtAll_%d_%d_%d",iBinLo,iBinUp,ic),iBinLo,iBinUp,"e");

      ScaleH1(hProjLead,"#eta","1/N_{evt} dN_{jet}/d#eta",1./fNevents,true);
      SetStyleH1(hProjLead,kBlue,0,kFullCircle);
      ScaleH1(hProjLeadFull,"#eta","1/N_{evt} dN_{jet}/d#eta",1./fNevents,true);
      SetStyleH1(hProjLeadFull,kRed,0,kFullCircle);
      ScaleH1(hProjAll,"#eta","1/N_{evt} dN_{jet}/d#eta",1./fNevents,true);
      SetStyleH1(hProjAll,kGray,1001,kFullCircle);


      hProjAll->DrawCopy("hist");
      hProjLead->DrawCopy("histsame");
      hProjLeadFull->DrawCopy("histsame");

      Float_t ptUp = fStepUp-h2EtaPtLead[0]->GetYaxis()->GetBinWidth(iBinUp);

      Printf("%3.0f-%3.0f",fStepLo,ptUp);
      TLatex *txt = 0;
      txt = new TLatex(-1.1,1.1*hProjAll->GetMaximum(),Form("%s %s p_{T} = %3.0f-%3.0f GeV",sCen[ic].Data(),sinputJet.Data(),fStepLo,ptUp));
      txt->SetTextSize(gStyle->GetTextSize());
      txt->Draw();
      c1->SaveAs(Form("%s_%s_Cen%02d_%03d-%03d_eta.%s",sinputPrint.Data(),cAdd.Data(),ic,(Int_t)fStepLo,(Int_t)ptUp,printType.Data()));
      c1->Modified();
      c1->Update();
      if(!gROOT->IsBatch())
	if(getchar()=='q')return 1;
      fStepLo += fStep;
      fStepUp += fStep;
      delete hProjAll;
      delete hProjLead;
      delete hProjLeadFull;
    }
    
    // Phi vs. p_T

 
    TH2F *h2PhiPtLead[nJets+2] = {0,};
    TH2F *h2PhiPtLeadFull[nJets+2] = {0,};
    TH2F *h2PhiPtAll[nJets+2] = {0,};
    HistsFromSingleJets(Form("fh2PhiPt%s_j%%d",cAdd.Data()),inputList,1,h2PhiPtLead);
    HistsFromSingleJets(Form("fh2PhiPt%s_j%%d",cAddFull.Data()),inputList,1,h2PhiPtLeadFull);
    HistsFromSingleJets(Form("fh2PhiPt%s_j%%d",cAddFull.Data()),inputList,nJets+1,h2PhiPtAll,1,nJets);

    c1->SetLogz();    
    h2PhiPtLead[0]->Draw("colz");    
    c1->SaveAs(Form("%s_%s_Cen%02d_phiPtLead2D.%s",sinputPrint.Data(),cAdd.Data(),ic,printType.Data()));
    c1->Modified();
    c1->Update();

    if(!gROOT->IsBatch())
      if(getchar()=='q')return 1;

    fStepLo = 0;
    fStepUp = fStepLo+fStep;

    while(fStepUp<fMax){
      Int_t iBinLo = h2PhiPtLead[0]->GetYaxis()->FindBin(fStepLo);
      Int_t iBinUp = h2PhiPtLead[0]->GetYaxis()->FindBin(fStepUp)-1;
      TH1D *hProjLead = h2PhiPtLead[0]->ProjectionX(Form("hPhiPtLead_%d_%d_%d",iBinLo,iBinUp,ic),iBinLo,iBinUp,"e");
      TH1D *hProjLeadFull = h2PhiPtLeadFull[0]->ProjectionX(Form("hPhiPtLeadFull_%d_%d_%d",iBinLo,iBinUp,ic),iBinLo,iBinUp,"e");
      TH1D *hProjAll = h2PhiPtAll[nJets+1]->ProjectionX(Form("hPhiPtAll_%d_%d_%d",iBinLo,iBinUp,ic),iBinLo,iBinUp,"e");

      ScaleH1(hProjLead,"#phi","1/N_{evt} dN_{jet}/d#phi",1./fNevents,true);
      SetStyleH1(hProjLead,kBlue,0,kFullCircle);
      ScaleH1(hProjLeadFull,"#phi","1/N_{evt} dN_{jet}/d#phi",1./fNevents,true);
      SetStyleH1(hProjLeadFull,kRed,0,kFullCircle);
      ScaleH1(hProjAll,"#phi","1/N_{evt} dN_{jet}/d#phi",1./fNevents,true);
      SetStyleH1(hProjAll,kGray,1001,kFullCircle);

      hProjAll->SetMinimum(0);

      hProjAll->DrawCopy("hist");
      hProjLead->DrawCopy("histsame");
      hProjLeadFull->DrawCopy("histsame");

      Float_t ptUp = fStepUp-h2PhiPtLead[0]->GetYaxis()->GetBinWidth(iBinUp);

      Printf("%3.0f-%3.0f",fStepLo,ptUp);
      TLatex *txt = 0;
      txt = new TLatex(0.0,1.1*hProjAll->GetMaximum(),Form("%s %s p_{T} = %3.0f-%3.0f GeV",sCen[ic].Data(),sinputJet.Data(),fStepLo,ptUp));
      txt->SetTextSize(gStyle->GetTextSize());
      txt->Draw();
      c1->SaveAs(Form("%s_%s_Cen%02d_%03d-%03d_phi.%s",sinputPrint.Data(),cAdd.Data(),ic,(Int_t)fStepLo,(Int_t)ptUp,printType.Data()));
      c1->Modified();
      c1->Update();
      if(!gROOT->IsBatch())
	if(getchar()=='q')return 1;
      fStepLo += fStep;
      fStepUp += fStep;
      delete hProjAll;
      delete hProjLead;
      delete hProjLeadFull;
    }
    
    // area vs. p_T


    TH2F *h2AreaPtLead[nJets+2] = {0,};
    TH2F *h2AreaPtLeadFull[nJets+2] = {0,};
    TH2F *h2AreaPtAll[nJets+2] = {0,};

    HistsFromSingleJets(Form("fh2AreaPt%s_j%%d",cAdd.Data()),inputList,1,h2AreaPtLead);
    HistsFromSingleJets(Form("fh2AreaPt%s_j%%d",cAddFull.Data()),inputList,1,h2AreaPtLeadFull);
    HistsFromSingleJets(Form("fh2AreaPt%s_j%%d",cAddFull.Data()),inputList,nJets+1,h2AreaPtAll,1,nJets);

    c1->SetLogz();    
    h2AreaPtLead[0]->Draw("colz");    
    c1->Modified();
    c1->Update();
    c1->SaveAs(Form("%s_%s_Cen%02d_areaPtLead2D.%s",sinputPrint.Data(),cAdd.Data(),ic,printType.Data()));
    if(!gROOT->IsBatch())
      if(getchar()=='q')return 1;

    fStepLo = 0;
    fStepUp = fStepLo+fStep;

    while(fStepUp<fMax){
      Int_t iBinLo = h2AreaPtLead[0]->GetYaxis()->FindBin(fStepLo);
      Int_t iBinUp = h2AreaPtLead[0]->GetYaxis()->FindBin(fStepUp)-1;
      TH1D *hProjLead = h2AreaPtLead[0]->ProjectionX(Form("hAreaPtLead_%d_%d_%d",iBinLo,iBinUp,ic),iBinLo,iBinUp,"e");
      TH1D *hProjLeadFull = h2AreaPtLeadFull[0]->ProjectionX(Form("hAreaPtLeadFull_%d_%d_%d",iBinLo,iBinUp,ic),iBinLo,iBinUp,"e");
      TH1D *hProjAll = h2AreaPtAll[nJets+1]->ProjectionX(Form("hAreaPtAll_%d_%d_%d",iBinLo,iBinUp,ic),iBinLo,iBinUp,"e");

      ScaleH1(hProjLead,"Area A","1/N_{evt} dN_{jet}/dA",1./fNevents,true);
      SetStyleH1(hProjLead,kBlue,0,kFullCircle);
      ScaleH1(hProjLeadFull,"Area A","1/N_{evt} dN_{jet}/dA",1./fNevents,true);
      SetStyleH1(hProjLeadFull,kRed,0,kFullCircle);
      ScaleH1(hProjAll,"Area A","1/N_{evt} dN_{jet}/dA",1./fNevents,true);
      SetStyleH1(hProjAll,kGray,1001,kFullCircle);

      hProjAll->SetMinimum(0);
      hProjAll->DrawCopy("hist");
      hProjLead->DrawCopy("histsame");
      hProjLeadFull->DrawCopy("histsame");

      Float_t ptUp = fStepUp-h2AreaPtLead[0]->GetYaxis()->GetBinWidth(iBinUp);

      Printf("%3.0f-%3.0f",fStepLo,ptUp);
      TLatex *txt = 0;
      txt = new TLatex(0.,1.1*hProjAll->GetMaximum(),Form("%s %s p_{T} = %3.0f-%3.0f GeV",sCen[ic].Data(),sinputJet.Data(),fStepLo,ptUp));
      txt->SetTextSize(gStyle->GetTextSize());
      txt->Draw();
      c1->SaveAs(Form("%s_%s_Cen%02d_%03d-%03d_area.%s",sinputPrint.Data(),cAdd.Data(),ic,(Int_t)fStepLo,(Int_t)ptUp,printType.Data()));
      c1->Modified();
      c1->Update();
      if(!gROOT->IsBatch())
	if(getchar()=='q')return 1;
      fStepLo += fStep;
      fStepUp += fStep;
      delete hProjAll;
      delete hProjLead;
      delete hProjLeadFull;
    }
    
    // Area vs eta:
    TH2F *h2AreaEta[nJets+1] = {0,};
    HistsFromSingleJets(Form("fh2EtaArea%s_j%%d",cAdd.Data()),inputList,nJets,h2AreaEta);
    ScaleH2(h2AreaEta[nJets],"area A","#eta","1/N_{evt} dN_{jets}/dAd#eta",1./fNevents,true);
    c1->SetLogz();     
    h2AreaEta[nJets]->Draw("colz");
    c1->Modified();
    c1->Update();
    c1->SaveAs(Form("%s_%s_Cen%02d_areaEta2D.%s",sinputPrint.Data(),cAdd.Data(),ic,printType.Data()));
    if(!gROOT->IsBatch())
      if(getchar()=='q')return 1;
  }

  return 0;
}

Int_t PlotSpectraPbPb2(){

  // Using now the THNSparse histos
  // 
  const int kMaxFiles = 2;
  TString sinputFile[kMaxFiles];
  TString sinputDir[kMaxFiles];

  TCanvas *c1 = new TCanvas();
  c1->SetLogy();
  THnSparseF *fhnJetPtRec[kMaxFiles] = {0,};
  THnSparseF *fhnTrackPtRec[kMaxFiles] = {0,};
  THnSparseF *fhnEvent[kMaxFiles] = {0,};

  const int nCen = 4;
  const Float_t nColl[nCen] = {1502,742.5,250,46.6};
  const Float_t fCentLo[nCen] = {0,10,30,50};
  const Float_t fCentUp[nCen] = {10,30,50,80};
  const Int_t iColCen[nCen] = {kRed+4,kRed,kBlue,kGray+1};

  Int_t iF = 0;
  sinputFile[iF] = "~/Dropbox/SharedJets/Christian/Files/PWG4_JetTasksOutput_LHC10h_AOD_tmp.root";
  //  sinputDir[iF] = "PWG4_spec2_clustersAOD_ANTIKT04_B2_Filter00128_Cut02000_Skip02_clustersAOD_ANTIKT04_B0_Filter00128_Cut02000_Skip02_0000000000_Class00";
  sinputDir[iF] = "PWG4_spec2_clustersAOD_ANTIKT04_B2_Filter00128_Cut00150_Skip02__0000000000_Class00";
  fhnJetPtRec[iF] =  (THnSparseF*)GetObjectFromFile("fhnJetPtRec",sinputFile[iF].Data(),sinputDir[iF].Data());
  fhnTrackPtRec[iF] =  (THnSparseF*)GetObjectFromFile("fhnTrackPtRec",sinputFile[iF].Data(),sinputDir[iF].Data());


  fhnEvent[iF] = (THnSparseF*) GetObjectFromFile("fhnEvent",sinputFile[iF].Data(),sinputDir[iF].Data());

  Bool_t bFirst1 = true;
  for(int ic = 0;ic <nCen;ic++){
    // project out the jets
    Int_t ibLo =      fhnJetPtRec[iF]->GetAxis(2)->FindBin(fCentLo[ic]+0.001);
    Int_t ibUp =      fhnJetPtRec[iF]->GetAxis(2)->FindBin(fCentUp[ic]-0.001);
    fhnJetPtRec[iF]->GetAxis(2)->SetRange(ibLo,ibUp);
    fhnJetPtRec[iF]->GetAxis(0)->SetRange(3,3); // take all jets
    TH1D *hSpectrumJets = fhnJetPtRec[iF]->Projection(1,"E");
    hSpectrumJets->SetName(Form("hSpectrumJets_C%d",ic));
    hSpectrumJets->SetMarkerColor(iColCen[ic]);
    hSpectrumJets->SetMarkerStyle(kFullCircle);
    hSpectrumJets->Rebin(2);
    hSpectrumJets->SetAxisRange(0,120);
    // project out the tracks
    fhnTrackPtRec[iF]->GetAxis(2)->SetRange(ic+1,ic+1);
    //    fhnTrackPtRec[iF]->GetAxis(0)->SetRange(2,2); // take all tracks
    TH1D *hSpectrumTracks = fhnTrackPtRec[iF]->Projection(1,"E");
    hSpectrumTracks->SetName(Form("hSpectrumTracks_C%d",ic));


    c1->cd();
    if(bFirst1){
      hSpectrumJets->DrawCopy("p");
      // hSpectrumTracks->DrawCopy("p");
      bFirst1 = false;
    }
    else {
      hSpectrumJets->DrawCopy("psame");
      //      hSpectrumTracks->DrawCopy("psame");
    }

    c1->Update();
    if(getchar()=='q')return 1;
  }
  return 0;
}

Int_t PlotSubtract(){

  // Using now the THNSparse histos
  // 

  TCanvas *c1 = new TCanvas(); // create with good aspect ratios...
  c1->SetLogz();

  const int kMaxFiles = 2;
  TString sinputFile[kMaxFiles];
  TString sinputDir[kMaxFiles];
  TString sLegend[kMaxFiles];




  // Parameter
  Float_t fMaxMult[kMaxFiles] = {0,};
  Float_t fMaxRho[kMaxFiles] = {0,};


  TH2F* fh2CentvsRho[kMaxFiles]= {0,};
  TH2F* fh2MultvsRho[kMaxFiles]= {0,};
  


  const Int_t nSub = 3;
  THnSparseF* hnDPtAreaCentMult[kMaxFiles][nSub] = {0,};

  Int_t iF = 0;
  sinputFile[iF] = "~/Dropbox/SharedJets/Christian/Files/PWG4_JetTasksOutput_LHC10h_AOD_tmp.root";
  sinputDir[iF] = "PWG4_JetSubtract_B2";
  sLegend[iF] = "B2 (150 MeV)";
  fMaxMult[iF] = 3500;
  fMaxRho[iF] = 250;
  iF++;
  
  sinputFile[iF] = "~/Dropbox/SharedJets/Christian/Files/PWG4_JetTasksOutput_LHC10h_AOD_tmp.root";
  sinputDir[iF] = "PWG4_JetSubtract_B2_Cut2000";
  sLegend[iF] = "B2 (2000 MeV)";
  fMaxMult[iF] = 300;
  fMaxRho[iF] = 50;
  
  

  for(iF = 0;iF<kMaxFiles;iF++){
    if(sinputFile[iF].Length()==0)continue;
    // fetch all the histos
    fh2CentvsRho[iF]= (TH2F*) GetObjectFromFile("fh2CentvsRho",sinputFile[iF].Data(),sinputDir[iF].Data(),"pwg");
    fh2MultvsRho[iF]= (TH2F*) GetObjectFromFile("fh2MultvsRho",sinputFile[iF].Data(),sinputDir[iF].Data(),"pwg");
    

    for(int iSub = 0;iSub<nSub;iSub++){
      hnDPtAreaCentMult[iF][iSub] = (THnSparseF*) GetObjectFromFile(Form("hnDPtAreaCentMult_%d",iSub),sinputFile[iF].Data(),sinputDir[iF].Data(),"pwg");
    }
    
    TLatex *txt = new TLatex();
    txt->SetNDC();

    fh2CentvsRho[iF]->SetAxisRange(0.,fMaxRho[iF],"Y");
    fh2CentvsRho[iF]->Draw("colz");
    txt->DrawLatex(0.6,0.8,sLegend[iF].Data());
    c1->Update();
    picName = Form("%sCentVsRho_%d.%s",picPrefix.Data(),iF,picSuffix.Data());
    c1->SaveAs(picName.Data());
    if(getchar()=='q')return 1;
    
    fh2MultvsRho[iF]->SetAxisRange(0.,fMaxRho[iF],"Y");
    fh2MultvsRho[iF]->SetAxisRange(0.,fMaxMult[iF]);
    fh2MultvsRho[iF]->Draw("colz");
    txt->DrawLatex(0.9,0.6,sLegend[iF].Data());
    c1->Update();
    picName = Form("%sMultVsRho_%d.%s",picPrefix.Data(),iF,picSuffix.Data());
    c1->SaveAs(picName.Data());
    if(getchar()=='q')return 1;
  }

  return 0;
  Bool_t bFirst1 = true;
  iF = 0;
  for(int ic = 0;ic <10;ic++){  
    Int_t iSub = 2;
    hnDPtAreaCentMult[iF][iSub]->GetAxis(2)->SetRange(ic+1,ic+1);
    TH2D *hDptArea = hnDPtAreaCentMult[iF][iSub]->Projection(0,1,"E");
    hDptArea->SetName(Form("hDptArea_C%d_%d",ic,iSub));
    hDptArea->Draw("colz");
    c1->Update();
    if(getchar()=='q')return 1;
  }

  
  return 0;

}


Int_t PlotSpectraPbPb(){

  // PLot the simple 1D histos from the spectrum task
  //
  const Int_t iRebin = 4;
  const int kMaxFiles = 5;

  Double_t yMin = 0.01;
  Double_t yMax = 1E+07;

  TString sinputFile[kMaxFiles];
  TFile*  inputFile[kMaxFiles] = {kMaxFiles*0};
  TString sinputDir[kMaxFiles];
  TDirectory *inputDir[kMaxFiles] = {kMaxFiles*0};
  TString sinputList[kMaxFiles];  
  TString sCen[kMaxFiles];  
  TList*  inputList[kMaxFiles] = {kMaxFiles*0};
  //  TString sinputLegend[kMaxFiles];
  TString sinputJet[kMaxFiles];
  Int_t kColor[kMaxFiles];
  Float_t nColl[kMaxFiles];

  bool bLogLog = false;
  //  const Int_t kRef = 1; // anti -kT

  // 
  TH1F* hJets[kMaxFiles];

  const int nJets = 2;

  Int_t iJF = 0;
  sinputFile[iJF]  = "~/alice/data/analysis/train/LHC10h_110116/PWG4_JetTasksOutput.root";
  //  sinputFile[iJF]  = "~/alice/data/analysis/train_maf/output_110210.root";
  sinputDir[iJF]  = "PWG4_spec2_clustersAOD_ANTIKT04_B3_Filter00256_Cut00150_Skip00_clustersAOD_ANTIKT04_B1_Filter00256_Cut00150_Skip00_0000000000_Class00";
  sinputList[iJF]  = "pwg4spec2_clustersAOD_ANTIKT04_B3_Filter00256_Cut00150_Skip00_clustersAOD_ANTIKT04_B1_Filter00256_Cut00150_Skip00_0000000000_Class00";
  sinputJet[iJF] = "anti-k_{T} R = 0.4";
  sCen[iJF] = " 0-80%";
  kColor[iJF] = kBlack;
  nColl[iJF] = 453.3; 
  iJF++;

  //  sinputFile[iJF]  = "../output/PWG4_JetTasksOutput_Merge_051205b_Subtract3_UA1.root"; 
  sinputFile[iJF]  = "~/alice/data/analysis/train/LHC10h_110116/PWG4_JetTasksOutput.root";
  //  sinputFile[iJF]  = "~/alice/data/analysis/train_maf/output_110210.root";
  sinputDir[iJF]  = "PWG4_spec2_clustersAOD_ANTIKT04_B3_Filter00256_Cut00150_Skip00_clustersAOD_ANTIKT04_B1_Filter00256_Cut00150_Skip00_0000000000_Class01";
  sinputList[iJF]  = "pwg4spec2_clustersAOD_ANTIKT04_B3_Filter00256_Cut00150_Skip00_clustersAOD_ANTIKT04_B1_Filter00256_Cut00150_Skip00_0000000000_Class01"; 
 sinputJet[iJF] = " 0-10% Anti-k_{T} R = 0.4";
  sCen[iJF] = " 0-10%";
  kColor[iJF] = kRed;
  nColl[iJF] = 1502; 
  iJF++;

  //  sinputFile[iJF]  = "../output/PWG4_JetTasksOutput_Merge_051205b_Subtract3_UA1.root"; 
  sinputFile[iJF]  = "~/alice/data/analysis/train/LHC10h_110116/PWG4_JetTasksOutput.root";
  //  sinputFile[iJF]  = "~/alice/data/analysis/train_maf/output_110210.root";
  sinputDir[iJF]  = "PWG4_spec2_clustersAOD_ANTIKT04_B3_Filter00256_Cut00150_Skip00_clustersAOD_ANTIKT04_B1_Filter00256_Cut00150_Skip00_0000000000_Class02";
  sinputList[iJF]  = "pwg4spec2_clustersAOD_ANTIKT04_B3_Filter00256_Cut00150_Skip00_clustersAOD_ANTIKT04_B1_Filter00256_Cut00150_Skip00_0000000000_Class02"; 
  sinputJet[iJF] = "10-30% Anti-k_{T} R = 0.4";
  sCen[iJF] = "10-30%";
  kColor[iJF] = kRed-4;
  nColl[iJF] = 742.5; 
  iJF++;

  //  sinputFile[iJF]  = "../output/PWG4_JetTasksOutput_Merge_051205b_Subtract3_UA1.root"; 
  sinputFile[iJF]  = "~/alice/data/analysis/train/LHC10h_110116/PWG4_JetTasksOutput.root";
  //  sinputFile[iJF]  = "~/alice/data/analysis/train_maf/output_110210.root";
  sinputDir[iJF]  = "PWG4_spec2_clustersAOD_ANTIKT04_B3_Filter00256_Cut00150_Skip00_clustersAOD_ANTIKT04_B1_Filter00256_Cut00150_Skip00_0000000000_Class03";
  sinputList[iJF]  = "pwg4spec2_clustersAOD_ANTIKT04_B3_Filter00256_Cut00150_Skip00_clustersAOD_ANTIKT04_B1_Filter00256_Cut00150_Skip00_0000000000_Class03"; 
  sinputJet[iJF] = "30-50% Anti-k_{T} R = 0.4";
  sCen[iJF] = "30-50%";
  kColor[iJF] = kRed-7;
  nColl[iJF] = 250; 
  iJF++;

  //  sinputFile[iJF]  = "../output/PWG4_JetTasksOutput_Merge_051205b_Subtract3_UA1.root"; 
  sinputFile[iJF]  = "~/alice/data/analysis/train/LHC10h_110116/PWG4_JetTasksOutput.root";
  //  sinputFile[iJF]  = "~/alice/data/analysis/train_maf/output_110210.root";
  sinputDir[iJF]  = "PWG4_spec2_clustersAOD_ANTIKT04_B3_Filter00256_Cut00150_Skip00_clustersAOD_ANTIKT04_B1_Filter00256_Cut00150_Skip00_0000000000_Class04";
  sinputList[iJF]  = "pwg4spec2_clustersAOD_ANTIKT04_B3_Filter00256_Cut00150_Skip00_clustersAOD_ANTIKT04_B1_Filter00256_Cut00150_Skip00_0000000000_Class04"; 
  sinputJet[iJF] = "Anti-k_{T} R = 0.4";
  sCen[iJF] = "50-80%";
  kColor[iJF] = kRed-9;
  nColl[iJF] = 46.6; 
  iJF++;


  for(int iF = 0; iF < kMaxFiles;iF++){
    if(sinputFile[iF].Length()==0)continue;
    if(gROOT->GetListOfFiles()->FindObject(sinputFile[iF].Data())){
      inputFile[iF] = (TFile*)gROOT->GetListOfFiles()->FindObject(sinputFile[iF].Data());
    }
    else {
      inputFile[iF] = TFile::Open(sinputFile[iF].Data());
    }
    inputDir[iF] = (TDirectory*)inputFile[iF]->Get(sinputDir[iF].Data());
    if(!inputDir[iF]){
      Printf("Dir not found %s",sinputDir[iF].Data());
      continue;
    }
    inputList[iF] = (TList*)inputDir[iF]->Get(sinputList[iF].Data());
    if(!inputList[iF]){
      Printf("List not found %s",sinputList[iF].Data());
      continue;
    }
  }


 

  TCanvas *c1 = new TCanvas("c1","spectra",20,20,800,600);
  TCanvas *c2 = new TCanvas("c2","ratio",820,20,800,600);

  TLegend *leg1 = new TLegend(0.45,0.65,0.7,0.85);
  leg1->SetHeader(Form("%s |#eta| < 0.4",sinputJet[0].Data()));
  leg1->SetFillColor(0);
  leg1->SetTextFont(gStyle->GetTextFont());
  leg1->SetBorderSize(0);

  char *cMaskRec = "fh1PtRecIn_j%d";
  
  c1->SetLogy();
  if(bLogLog)c1->SetLogx();
    
  for(int iF = 0; iF < kMaxFiles;iF++){
    if(sinputFile[iF].Length()==0)continue;
    // Draw the jet spectra
    if(!inputList[iF])continue;
    TH1F *hRec[nJets+1];
    HistsFromSingleJets(cMaskRec,inputList[iF],nJets,hRec,iRebin);
    /// TH1F *hRec[nJets] = (TH1F*)inputList[iF]->FindObject("fh1PtJetsRecIn");    hRec[nJets]->Rebin(iRebin);

    hRec[nJets]->SetMarkerStyle(kFullCircle);
    hRec[nJets]->SetXTitle("p_{T} (GeV/c)");
    hRec[nJets]->SetMarkerColor(kColor[iF]);
    hRec[nJets]->Scale(1./1.6); // deta
    hRec[nJets]->Scale(1./hRec[nJets]->GetBinWidth(1));

    hRec[nJets]->SetYTitle("dN/dp_{T}dy");

    if(true){

      ScaleNevent(hRec[nJets], inputFile[iF],AliAnalysisTaskJetServices::kSelected,0,iF); // <----- careful assumes  iCen == iF      
      hRec[nJets]->SetYTitle("1/N_{evt} dN/dp_{T}");
      //      Printf("%d Ncall %f",iF,nColl[iF]); 
      //      hRec[nJets]->Scale(1./nColl[iF]);
      //      hRec[nJets]->SetYTitle("1/(N_{coll}*(N_{evt}) dN/dp_{T}dy");
      yMin = 1E-8;
      yMax = 1;
    }
    if(sCen[iF].Length()){
      leg1->AddEntry(hRec[nJets],Form("%s",sCen[iF].Data()),"P");
    }

    hJets[iF] = hRec[nJets];
    c1->cd();

  
    if(bLogLog)hRec[nJets]->SetAxisRange(5,100);
    else hRec[nJets]->SetAxisRange(1,100);

    if(iF==0){
      if(yMax>0) hRec[nJets]->SetMaximum(yMax);
      if(yMin>0) hRec[nJets]->SetMinimum(yMin);
      hRec[nJets]->DrawCopy("P");
    }
    else{
     if(sinputJet[iF].Length())hRec[nJets]->DrawCopy("Psame");
    }
    c1->Update();
    if(getchar()=='q')return 1;
  }


  c1->Update();
  picName = "jetspectrumPbPb";
  if(bLogLog)picName += "LogLog";
  leg1->Draw("same");
  TLatex *txt = 0;
  txt = new TLatex(5,0.1,"LHC2010 Pb+Pb #sqrt{s_{NN}} = 2.76 TeV");
  txt->SetTextSize(gStyle->GetTextSize()*0.8);
  txt->Draw();
  ALICEWorkInProgress(c1,"02/15/2011");
  c1->SaveAs(Form(cPrintMask,picName.Data()));
  if(getchar()=='q')return 1;


  c2->cd();
  leg1->Draw("same");




  picName = "jetspectrumLabels";
  c2->SaveAs(Form(cPrintMask,picName.Data()));
  if(getchar()=='q')return 1;

  // scale with nColl
  for(int iF = 0; iF < kMaxFiles;iF++){
    if(sinputFile[iF].Length()==0)continue;
    if(!hJets[iF])continue;
    Printf("%d Ncoll %f",iF,nColl[iF]); 
    hJets[iF]->Scale(1./nColl[iF]);
    hJets[iF]->SetYTitle("1/(N_{coll}*(N_{evt}) dN/dp_{T}dy");
  }

  const Int_t nPSmear = 9;
  TH1F *hPythia[nPSmear];
  // Fetch the smeared pythia spectra

  TString pName;
  c2->SetLogy();
  TH1F* hRatio[kMaxFiles];

  for(int i = 0;i < nPSmear;i++){
    if(i==nPSmear-1){
      pName = "h1JetPtCh";
      hPythia[i] = GetPythia8Spectrum(pName.Data(),"~/alice/sim/pythia8/examples/output/LHCJets_2.75TeV_allpt_R04.root",1.0,iRebin);
    }
    else {
      pName = Form("h1JetPtCh_C%02d_j%%d",i);
      hPythia[i] = GetPythia8Spectrum(pName.Data(),"~/alice/sim/pythia8/examples/output/LHCJets_2.75TeV_allpt_R04.root",1.0,iRebin,nJets);
    } 
    hPythia[i]->SetMarkerStyle(kOpenCircle);
    hPythia[i]->SetMarkerColor(kBlack);
    hPythia[i]->SetYTitle("1/(N_{coll}*N_{evt}) dN/dp_{T}dy");
    hPythia[i]->SetAxisRange(0,100);
    c1->cd();
    hPythia[i]->Draw("CP");
    hPythia[i]->SetMaximum(1E-02);
    hPythia[i]->SetMinimum(1E-09);
    for(int iF = 0; iF < kMaxFiles;iF++){
      if(sinputFile[iF].Length()==0)continue;
      if(!hJets[iF])continue;
      
      c1->cd();
      hJets[iF]->Draw("psame");

      c2->cd();
      hRatio[iF] = (TH1F*)hJets[iF]->Clone(Form("hRatio%d",iF));
      hRatio[iF]->SetMaximum(100);
      hRatio[iF]->SetMinimum(0.01);
      hRatio[iF]->Divide(hPythia[i]);
      hRatio[iF]->SetYTitle("Ratio");
      if(iF==0) hRatio[iF]->DrawCopy("p");
      else hRatio[iF]->DrawCopy("psame");
      c1->Update();
      c2->Update();
    }

    if(getchar()=='q')return 1;
    picName = Form("jetspectrumRatioPbPb_Smear%d",i);
    c2->SaveAs(Form(cPrintMask,picName.Data()));
    picName = Form("jetspectrumPbPb_Smear%d",i);
    c1->SaveAs(Form(cPrintMask,picName.Data()));
    if(getchar()=='q')return 1;

  }
  c1->cd();
  hPythia[nPSmear-1]->DrawCopy("P");
  for(int i = 0;i < nPSmear-1;i++){
    c1->cd();
    hPythia[i]->SetMarkerStyle(kFullCircle);
    hPythia[i]->SetMarkerColor(kRed-9+i);
    hPythia[i]->DrawCopy("psame");

    TH1F *hRatioS = (TH1F*)hPythia[i]->Clone(Form("hPythia_%d",i));
    hRatioS->Divide(hPythia[nPSmear-1]);
    hRatioS->SetMaximum(1E5);
    hRatioS->SetMinimum(0.1);
    hRatioS->SetYTitle("Ratio");
    c2->cd();
    if(i==0)hRatioS->DrawCopy("p");
    else hRatioS->DrawCopy("psame");
    c1->Update();
    c2->Update();
    if(getchar()=='q')return 1;
  }

  if(getchar()=='q')return 1;
  picName = Form("jetSmearRatio");
  c2->SaveAs(Form(cPrintMask,picName.Data()));
  picName = Form("jetSmear");
  c1->SaveAs(Form(cPrintMask,picName.Data()));
  if(getchar()=='q')return 1;

  // fetch the jet spectrum first


  c1->cd();
  TH1F *hPythiaRef = GetPythia8Spectrum(pName.Data(),"~/alice/sim/pythia8/examples/output/LHCJets_2.75TeV_allpt_R04.root",1.0,1);
  hPythiaRef->Draw();
  c1->Update();
  Float_t fMinPt[kMaxFiles];
  for(int iF = 0; iF < kMaxFiles;iF++){
    if(nColl[iF]<=0)continue;
      for(int i = hPythiaRef->GetNbinsX();i >0;i--){
	if(hPythiaRef->GetBinContent(i)>1./nColl[iF]){
	   fMinPt[iF] = hPythiaRef->GetBinCenter(i);
	   break;
	}
      }
      Printf("Ncoll %4.1f Min Pt %3.1f",nColl[iF],fMinPt[iF]);
  }

  char *cFile = "~/alice/sim/pythia8/examples/output/LHCJets_2.75TeV_allpt_R04.root";
  TFile *fIn = (TFile*)gROOT->GetListOfFiles()->FindObject(cFile);
  if(!fIn) fIn  = TFile::Open(cFile);
  if(!fIn)return 0;

  TH1D* hPythia2[nPSmear];
  for(int i = 0;i < nPSmear-1;i++){

    TH2F *hCorrelation = 0;
    for(int ij = 0;ij<2;ij++){
      TH2F* hTmp = (TH2F*)fIn->Get(Form("h2SmearCorrelationCh_C%02d_j%d",i,ij));
      if(hTmp&&!hCorrelation)hCorrelation = (TH2F*)hTmp->Clone(Form("%s_%d",hTmp->GetName(),i));
      else if (hTmp) hCorrelation->Add(hTmp);
    }
    c1->cd();
    
    hPythia2[i] = hCorrelation->ProjectionY(Form("hPythia2_%d",i),10,300);
    hPythia2[i]->SetMarkerStyle(kFullCircle);
    hPythia2[i]->SetMarkerColor(kRed-9+i);
    hPythia2[i]->Rebin(iRebin);
    hPythia2[i]->SetAxisRange(0,100);
    hPythia2[i]->Scale(1./hPythia2[i]->GetBinWidth(1)); 
    //  hist->Scale(1./xsec);
    hPythia2[i]->Scale(1./1.); // deta
    
    //  Double_t xtotal = 71.39; // xtotal for 7 TeV!!
    Double_t xtotal = 62.01; // xtotal for 2.76 TeV!!
    // scale with fraction of total NSD cross section
    hPythia2[i]->Scale(1/xtotal);



    c1->cd();
    if(i==0)hPythia2[i]->DrawCopy("p");
    else hPythia2[i]->DrawCopy("psame");

    TH1F *hRatioS = (TH1F*)hPythia2[i]->Clone(Form("hPythia2_%d",i));
    hRatioS->Divide(hPythia[nPSmear-1]);
    hRatioS->SetMaximum(1E5);
    hRatioS->SetMinimum(0.1);
    hRatioS->SetYTitle("Ratio");
    c2->cd();
    if(i==0)hRatioS->DrawCopy("p");
    else hRatioS->DrawCopy("psame");
    c1->Update();
    c2->Update();
    if(getchar()=='q')return 1;
  }
  return 0;
}

 

void set_plot_style() {
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}


Int_t PlotJetBFluctuations(){
  // plot the diffent background estimates

  const int nCen = 4;
  const Float_t nColl[nCen] = {1502,742.5,250,46.6};
   const Float_t fCentLo[nCen] = {0,10,30,50};
  const Float_t fCentUp[nCen] = {10,30,50,80};
  TH2F *hFrame = new TH2F("hFrame",";#delta p_{T} (GeV/c);Probability/GeV",200,-70,70,100,1E-5,50); 
  hFrame->SetTitleOffset(1.5,"Y");
  hFrame->SetTitleOffset(1.5,"X");
  hFrame->SetLabelSize(hFrame->GetLabelSize("Y")*0.9,"Y");
  hFrame->SetLabelSize(hFrame->GetLabelSize("X")*0.9,"X");
  hFrame->SetTitleSize(hFrame->GetTitleSize("Y")*0.7,"Y");
  hFrame->SetTitleSize(hFrame->GetTitleSize("X")*0.7,"X");


  TString printType = "png"; 
  TString tmpName; 

  TCanvas *c1 = new TCanvas("c11","c1",600,600);
  c1->SetLogy();

  //  TFile::SetCacheFileDir("/tmp/");

  // Martha, single particle jets
  TFile *fM = TFile::Open("~/alice/jets/macros/corrections/tmp/MV_PWG4_JetTasksOutput_AOD_EmbeddingSingleTrack.root");
  TH1D *hDeltaPtM[nCen] = {0};
  TString sDeltaPtM = "";

  // select 
  Float_t fMinPtM = 20;
  Float_t fMaxPtM = 40;
  int iB = 2;

  /*
    0: 0-10%
    1: 10-30%
    2: 30-50%
    3: 50-80%
  */

  for(int ic = 0;ic < nCen;ic++){
    tmpName = Form("PWG4_BkgFluctCent%dB%d",ic,iB);
    sDeltaPtM = Form("anti-k_{T} embedded tracks %2.0f-%2.0f GeV (MV)",fMinPtM,fMaxPtM);
    TDirectory *dir = (TDirectory*)fM->Get(tmpName.Data());
    if(!dir)Printf("Line:%d %s not found",__LINE__,tmpName.Data());
    tmpName = Form("taskBkgFluctCent%dB%d",ic,iB);
    TList *list = (TList*)dir->Get(tmpName.Data());    
    if(!list)Printf("Line:%d %s not found",__LINE__,tmpName.Data());

    TH3F *h3Tmp = (TH3F*)list->FindObject("fh3PtDeltaPtArea");
    hDeltaPtM[ic] = h3Tmp->ProjectionY(Form("hDeltaM%d",ic),h3Tmp->GetXaxis()->FindBin(fMinPtM),h3Tmp->GetXaxis()->FindBin(fMaxPtM),
				0,-1,"E");
    hDeltaPtM[ic]->Rebin(10);
    Float_t fScale = hDeltaPtM[ic]->Integral("width");
    if(fScale)hDeltaPtM[ic]->Scale(1./fScale);
    hDeltaPtM[ic]->SetMarkerStyle(33);
    hDeltaPtM[ic]->SetMarkerColor(kGreen+2);
    hDeltaPtM[ic]->SetLineColor( hDeltaPtM[ic]->GetMarkerColor());
  }    
  
  // fetch the BiAs
  
  TH1D *hBiaL[nCen];

  TFile *fL = TFile::Open("~/alice/jets/macros/corrections/tmp/2011-03-20_lcm_pw4plots.root");

  for(int ic = 0;ic < nCen;ic++){
    if(iB==1)tmpName = "BiA sa antikt centrality";
    else if (iB==2)tmpName = "BiA va antikt centrality";
    TH2F *h2Tmp = (TH2F*)fL->Get(tmpName.Data());
    if(!h2Tmp)Printf("Line:%d %s not found",__LINE__,tmpName.Data());
    Int_t ibLo = h2Tmp->GetXaxis()->FindBin(fCentLo[ic]);
    Int_t ibUp = h2Tmp->GetXaxis()->FindBin(fCentUp[ic])-1;
    Printf("Line:%d bin %d - %d",__LINE__,ibLo,ibUp);
    hBiaL[ic] = h2Tmp->ProjectionY(Form("hBiaL%d",ic),ibLo,ibUp,"E");
    Float_t fScale =  hBiaL[ic]->Integral("width");
    if(fScale)  hBiaL[ic]->Scale(1./fScale);
    hBiaL[ic]->SetMarkerStyle(kFullCircle);
    hBiaL[ic]->SetMarkerColor(kBlue);
    hBiaL[ic]->SetLineColor( hBiaL[ic]->GetMarkerColor());
  }
  
  TH1D *hBiaRC[nCen];
  TString sBiaRC = "";

  for(int ic = 0;ic < nCen;ic++){
    if(iB==1)tmpName = "BiA sa RC skip0 centrality";
    else if (iB==2)tmpName = "BiA va RC skip0 centrality";
    sBiaRC = "BiA RC";
    TH2F *h2Tmp = (TH2F*)fL->Get(tmpName.Data());
    if(!h2Tmp)Printf("Line:%d %s not found",__LINE__,tmpName.Data());
    Int_t ibLo = h2Tmp->GetXaxis()->FindBin(fCentLo[ic]);
    Int_t ibUp = h2Tmp->GetXaxis()->FindBin(fCentUp[ic])-1;
    Printf("Line:%d bin %d - %d",__LINE__,ibLo,ibUp);
    hBiaRC[ic] = h2Tmp->ProjectionY(Form("hBiaRC%d",ic),ibLo,ibUp,"E");
    Float_t fScale =  hBiaRC[ic]->Integral("width");
    if(fScale)  hBiaRC[ic]->Scale(1./fScale);
    hBiaRC[ic]->SetMarkerStyle(kFullCircle);
    hBiaRC[ic]->SetMarkerColor(kRed);
    hBiaRC[ic]->SetLineColor( hBiaRC[ic]->GetMarkerColor());
  }


  TH1D *hBiaRC2[nCen];
  TString sBiaRC2;
  for(int ic = 0;ic < nCen;ic++){
    if(iB==1)tmpName = "BiA sa RC skip2 centrality";
    else if (iB==2)tmpName = "BiA va RC skip2 centrality";
    sBiaRC2 = "BiA RC (excl. 2 leading jets)";
    TH2F *h2Tmp = (TH2F*)fL->Get(tmpName.Data());
    if(!h2Tmp)Printf("Line:%d %s not found",__LINE__,tmpName.Data());
    Int_t ibLo = h2Tmp->GetXaxis()->FindBin(fCentLo[ic]);
    Int_t ibUp = h2Tmp->GetXaxis()->FindBin(fCentUp[ic])-1;
    Printf("Line:%d bin %d - %d",__LINE__,ibLo,ibUp);
    hBiaRC2[ic] = h2Tmp->ProjectionY(Form("hBiaRC2%d",ic),ibLo,ibUp,"E");
    Float_t fScale =  hBiaRC2[ic]->Integral("width");
    if(fScale)  hBiaRC2[ic]->Scale(1./fScale);
    hBiaRC2[ic]->SetMarkerStyle(kOpenCircle);
    hBiaRC2[ic]->SetMarkerColor(kRed);
    hBiaRC2[ic]->SetLineColor( hBiaRC2[ic]->GetMarkerColor());
  }



  c1->SetLogy(0);
  c1->SetLogz();  
  TLatex *txt = new TLatex();
  txt->SetTextFont(gStyle->GetTextFont());
  txt->SetTextSize(gStyle->GetTextSize()*0.6);

  TLatex *txt2 = new TLatex();
  txt2->SetTextFont(gStyle->GetTextFont());
  txt2->SetTextSize(gStyle->GetTextSize()*0.7);
  txt2->SetTextAlign(22);
  txt2->SetTextColor(kRed);


  if(iB==1)tmpName = "background sa vs multiplicity";
  else if(iB==2)tmpName = "background va vs multiplicity";
  TH2F *h2RhoVsMult = (TH2F*)fL->Get(tmpName.Data());
  if(!h2RhoVsMult)Printf("Line:%d %s not found",__LINE__,tmpName.Data());
  c1->SetMargin(0.15,0.18,0.2,0.10);

  h2RhoVsMult->SetTitleOffset(1.5,"Y");
  h2RhoVsMult->SetTitleOffset(1.5,"X");
  h2RhoVsMult->SetLabelSize(h2RhoVsMult->GetLabelSize("Y")*0.9,"Y");
  h2RhoVsMult->SetLabelSize(h2RhoVsMult->GetLabelSize("X")*0.9,"X");
  h2RhoVsMult->SetTitleSize(h2RhoVsMult->GetTitleSize("Y")*0.7,"Y");
  h2RhoVsMult->SetTitleSize(h2RhoVsMult->GetTitleSize("X")*0.7,"X");
  h2RhoVsMult->SetXTitle("input tracks");
  h2RhoVsMult->SetYTitle("#rho (GeV/unit area)");
  h2RhoVsMult->SetAxisRange(0,3000.);
  h2RhoVsMult->Draw("colz");
  txt->Draw();
  txt->DrawLatex(100,180,"LHC 2010 Pb+Pb Run #sqrt{s_{NN}} = 2.76 TeV");
  //  txt2->DrawLatex(800,150,"ALICE Performance");
  //  txt2->DrawLatex(800,140,"01/03/2011");
  c1->Update();
  c1->SaveAs(Form("rhovsmult_B%d.%s",iB,printType.Data()));
  if(getchar()=='q')return 1;

  if(iB==1)tmpName = "background sa vs centrality";
  else if(iB==2)tmpName = "background va vs centrality";
  TH2F *h2RhoVsCent = (TH2F*)fL->Get(tmpName.Data());
  if(!h2RhoVsCent)Printf("Line:%d %s not found",__LINE__,tmpName.Data());
  h2RhoVsCent->SetXTitle("centrality (%)");
  h2RhoVsCent->SetYTitle("#rho (GeV/unit area)");
  
  h2RhoVsCent->SetTitleOffset(1.5,"Y");
  h2RhoVsCent->SetTitleOffset(1.5,"X");
  h2RhoVsCent->SetLabelSize(h2RhoVsCent->GetLabelSize("Y")*0.9,"Y");
  h2RhoVsCent->SetLabelSize(h2RhoVsCent->GetLabelSize("X")*0.9,"X");
  h2RhoVsCent->SetTitleSize(h2RhoVsCent->GetTitleSize("Y")*0.7,"Y");
  h2RhoVsCent->SetTitleSize(h2RhoVsCent->GetTitleSize("X")*0.7,"X");
  h2RhoVsCent->SetAxisRange(3,200,"Y");
  h2RhoVsCent->Draw("colz");
  //  txt->DrawLatex(20,180,"LHC 2010 Pb+Pb Run #sqrt{s_{NN}} = 2.76 TeV");
  //  txt2->DrawLatex(50,150,"ALICE Performance");
  //  txt2->DrawLatex(50,140,"01/03/2011");
  c1->Update();
  c1->SaveAs(Form("rhovscent_B%d.%s",iB,printType.Data()));
  if(getchar()=='q')return 1;

  // fetch the data from bastian...
  Float_t fMinPtB = fMinPtM;
  Float_t fMaxPtB = fMaxPtM;


  TH1D *hDeltaPtB1[nCen] = {0};
  TString sDeltaPtB1 = "";
  for(int ic = 0;ic < nCen;ic++){
    sDeltaPtB1 = Form("anti-k_{T} embedded jet %2.0f-%2.0f GeV",fMinPtB,fMaxPtB);
    TH2F *hTmp = GetTH2PlotB("~/Dropbox/SharedJets/Bastian/Files/110428/",1,0,ic); // emb jets
    if(!hTmp)Printf("%d %s not found",__LINE__,tmpName.Data());
    int ibLo = hTmp->GetYaxis()->FindBin(fMinPtB);
    int ibUp = hTmp->GetYaxis()->FindBin(fMaxPtB)-1;
    hDeltaPtB1[ic] = hTmp->ProjectionX(Form("fHistDeltaPtB1_c%d",ic),ibLo,ibUp,"E");
    hDeltaPtB1[ic]->SetMarkerStyle(kFullSquare);
    hDeltaPtB1[ic]->SetMarkerColor(kBlue);
    hDeltaPtB1[ic]->SetLineColor(hDeltaPtB1[ic]->GetMarkerColor());
    hDeltaPtB1[ic]->Rebin(2);
    Float_t fScale =  hDeltaPtB1[ic]->Integral("width");
    if(fScale)  hDeltaPtB1[ic]->Scale(1./fScale);
  }


  TH1D *hDeltaPtB2[nCen] = {0};
  TString sDeltaPtB2 = "";
  for(int ic = 0;ic < nCen;ic++){
    sDeltaPtB2 = Form("anti-k_{T} embedded track %2.0f-%2.0f GeV BB",fMinPtB,fMaxPtB);
    TH2F *hTmp = GetTH2PlotB("~/Dropbox/SharedJets/Bastian/Files/110430/",0,0,ic); // emb jets
    if(!hTmp)Printf("%d %s not found",__LINE__,tmpName.Data());
    int ibLo = hTmp->GetYaxis()->FindBin(fMinPtB);
    int ibUp = hTmp->GetYaxis()->FindBin(fMaxPtB)-1;
    hDeltaPtB2[ic] = hTmp->ProjectionX(Form("fHistDeltaPtB2_c%d",ic),ibLo,ibUp,"E");
    hDeltaPtB2[ic]->SetMarkerStyle(33);
    hDeltaPtB2[ic]->SetMarkerColor(kBlue+4);
    hDeltaPtB2[ic]->SetLineColor(hDeltaPtB2[ic]->GetMarkerColor());
    hDeltaPtB2[ic]->Rebin(2);
    Float_t fScale =  hDeltaPtB2[ic]->Integral("width");
    if(fScale)  hDeltaPtB2[ic]->Scale(1./fScale);
  }



  c1->SetLogy();
  c1->SetMargin(0.15,0.05,0.2,0.05);

  TF1 *gaus = new TF1("gaus","gaus",-60,2);
  TF1 *gaus2 = new TF1("gaus2","gaus",-60,2);
  Double_t mean = 0;
  Double_t sigma = 0;
  Double_t sigma_err = 0;
  TF1* tmpGaus = 0;
  TF1 *gamma0 = new TF1("gamma0",Gamma0,-40,40,4);
  gamma0->SetParameters(1,1,1,1);

  for(int ic = 0;ic < nCen;ic++){
    TLegend *leg1 = new TLegend(0.2,0.65,0.3,0.93);
    leg1->SetHeader(Form("Pb+Pb %2.0f-%2.0f%% R = 0.4 (B%d)",fCentLo[ic],fCentUp[ic],iB));
    leg1->SetFillColor(0);
    leg1->SetTextFont(gStyle->GetTextFont());
    leg1->SetTextSize(gStyle->GetTextSize()*0.6);
    leg1->SetBorderSize(0);
    hFrame->DrawCopy();

    /*
    hBiaL[ic]->DrawCopy("psame");
    leg1->AddEntry(hBiaL[ic],Form("BiA anti-k_{T}"),"P");
    */

    /*
    hBiaRC2[ic]->DrawCopy("psame");
    tmpGaus = FitLHSgaus(hBiaRC2[ic],sigma,sigma_err,mean);
    tmpGaus->SetRange(-40,40);
    tmpGaus->SetLineColor( hBiaRC2[ic]->GetLineColor());
    tmpGaus->SetLineStyle(kDashed);
    tmpGaus->Draw("same");
    leg1->AddEntry(hBiaRC2[ic],sBiaRC2.Data(),"P");
    leg1->AddEntry(tmpGaus,Form("LHS Gaus fit: #mu = %2.2f, #sigma = %2.2f",mean,sigma),"L");
    */    

    hBiaRC[ic]->DrawCopy("psame");
    hBiaRC[ic]->Fit(gamma0);
    c1->Update();
    if(getchar()=='q')return 1;
    continue;


    tmpGaus = FitLHSgaus(hBiaRC[ic],sigma,sigma_err,mean);
    tmpGaus->SetRange(-40,40);
    tmpGaus->SetLineColor( hBiaRC[ic]->GetLineColor());
    tmpGaus->Draw("same");
    leg1->AddEntry(hBiaRC[ic],sBiaRC.Data(),"P");
    leg1->AddEntry(tmpGaus,Form("LHS Gaus fit: #mu = %2.2f, #sigma = %2.2f",mean,sigma),"L");
    
    hDeltaPtB1[ic]->Draw("psame");
    tmpGaus = FitLHSgaus(hDeltaPtB1[ic],sigma,sigma_err,mean);
    tmpGaus->SetRange(-40,40);
    tmpGaus->SetLineColor( hDeltaPtB1[ic]->GetLineColor());
    tmpGaus->Draw("same");
    leg1->AddEntry(hDeltaPtB1[ic],sDeltaPtB1.Data(),"P");
    leg1->AddEntry(tmpGaus,Form("LHS Gaus fit: #mu = %2.2f, #sigma = %2.2f",mean,sigma),"L");


    hDeltaPtB2[ic]->Draw("psame");
    tmpGaus = FitLHSgaus(hDeltaPtB2[ic],sigma,sigma_err,mean);
    tmpGaus->SetRange(-40,40);
    tmpGaus->SetLineColor( hDeltaPtB2[ic]->GetLineColor());
    tmpGaus->Draw("same");
    leg1->AddEntry(hDeltaPtB2[ic],sDeltaPtB2.Data(),"P");
    leg1->AddEntry(tmpGaus,Form("LHS Gaus fit: #mu = %2.2f, #sigma = %2.2f",mean,sigma),"L");
    
    hDeltaPtM[ic]->DrawCopy("psame");
    tmpGaus = FitLHSgaus(hDeltaPtM[ic],sigma,sigma_err,mean);
    tmpGaus->SetRange(-40,40);
    tmpGaus->SetLineColor( hDeltaPtM[ic]->GetLineColor());
    tmpGaus->Draw("same");
    leg1->AddEntry(hDeltaPtM[ic],sDeltaPtM.Data(),"P");
    leg1->AddEntry(tmpGaus,Form("LHS Gaus fit: #mu = %2.2f, #sigma = %2.2f",mean,sigma),"L");

    leg1->Draw();


    txt2->SetNDC();
    //    txt2->DrawLatex(0.5,0.97,"ALICE Performance 01/03/2011");
    c1->Update();
    c1->SaveAs(Form("deltaPt_pT_%d_%d_B%d_cen%02d.%s",(Int_t)fMinPtM,(Int_t)fMaxPtM,iB,ic,printType.Data()));
    if(getchar()=='q')return 1;


  }



  // now plot different centralities on one plot
  bool bFirst = true;
  hFrame->DrawCopy();
  TLegend *leg2 = new TLegend(0.2,0.65,0.3,0.93);
  leg2->SetHeader(Form("%s R = 0.4 (B%d)",sBiaRC2.Data(),iB));
  TH1D **hTmp1 = hBiaRC2;
  tmpName = "BiARC2";
  leg2->SetFillColor(0);
  leg2->SetTextFont(gStyle->GetTextFont());
  leg2->SetTextSize(gStyle->GetTextSize()*0.6);
  leg2->SetBorderSize(0);
  hFrame->DrawCopy();


  for(int ic = 0;ic<nCen;ic++){
    int iColor = hTmp1[ic]->GetMarkerColor()-4+2*ic;
    hTmp1[ic]->SetLineColor(iColor);
    hTmp1[ic]->SetMarkerColor(iColor);
    if(bFirst){
      hTmp1[ic]->Draw("psame");
      bFirst = false;
    }
    else{
      hTmp1[ic]->Draw("psame");
    }
    leg2->AddEntry(hTmp1[ic],Form("%2.0f-%2.0f%%",fCentLo[ic],fCentUp[ic]),"P");
    
  }
  leg2->Draw("same");
  c1->Update();
  c1->SaveAs(Form("deltaPt_cent_%s_B%d.%s",tmpName.Data(),iB,printType.Data()));
  if(getchar()=='q')return 1;

  bFirst = true;
  for(int ic = nCen-1;ic>=0;ic--){
    hTmp1[ic]->Scale(1./nColl[ic]);
    hTmp1[ic]->SetMinimum(1E-8);
    if(bFirst){
      hTmp1[ic]->Draw("p");
      bFirst = false;
    }
    else{
      hTmp1[ic]->Draw("psame");
    }
    
  }
  leg2->Draw("same");
  c1->Update();
  c1->SaveAs(Form("deltaPt_cent_ncoll_%s_B%d.%s",tmpName.Data(),iB,printType.Data()));



  return 0;


}


Int_t PlotJetBFluctuations2(UInt_t iPlotFlag,UInt_t iPlotType,Float_t inPtLo,Float_t inPtUp){
  // plot the diffent background estimates

  const int nCen = 4;
  const Float_t nColl[nCen] = {1502,742.5,250,46.6};
  const Float_t fCentLo[nCen] = {0,10,30,50};
  const Float_t fCentUp[nCen] = {10,30,50,80};

  Int_t iB = 2;



  // multiplicity (nb. of tracklets) bins
  const Int_t nMult = 17;
  const Int_t multBinWidth = 200;
  Int_t multMin[nMult] = {0,};
  Int_t multMax[nMult] = {0,};
  for(Int_t i=0; i<nMult; ++i){
    multMin[i] = i*multBinWidth;
    multMax[i] = multMin[i]+multBinWidth-1;
  }

  // bins wrt. reactions
  const Int_t nRP = 4; // 0 --> all 1,2,3

  TH2F *hFrame = new TH2F("hFrame",";#delta p_{T} (GeV/c);Probability/GeV",200,-70,100,1000,1E-7,50); 
  hFrame->SetTitleOffset(1.5,"Y");
  hFrame->SetTitleOffset(1.5,"X");
  hFrame->SetLabelSize(hFrame->GetLabelSize("Y")*0.9,"Y");
  hFrame->SetLabelSize(hFrame->GetLabelSize("X")*0.9,"X");
  hFrame->SetTitleSize(hFrame->GetTitleSize("Y")*0.7,"Y");
  hFrame->SetTitleSize(hFrame->GetTitleSize("X")*0.7,"X");

  TCanvas *c1 = new TCanvas("c11","c1",600,600);
  c1->SetLogy();


  // Produce Delta Pt Plots for each centrality and for each type
  const int kDeltaTypes = 9;
  TString sDelta[kDeltaTypes][nRP];
  TH1D  *hDeltaPt[kDeltaTypes][nCen][nRP] = {0,};    
  TGraphErrors *grMeanDeltaPtCent[kDeltaTypes][nRP] = {0,};
  TGraphErrors *grSigmaDeltaPtCent[kDeltaTypes][nRP] = {0,};

  // 
  TH1D  *hDeltaPtMult[kDeltaTypes][nMult][nRP] = {0,};    
  TGraphErrors *grMeanDeltaPtMult[kDeltaTypes][nRP] = {0,};
  TGraphErrors *grSigmaDeltaPtMult[kDeltaTypes][nRP] = {0,};

  // 
  const Int_t kMarkerDelta[kDeltaTypes] = {kFullCircle,kFullCircle,kOpenCircle,kFullSquare,    33,kFullSquare,kFullSquare,          31,        31};
  const Int_t kColorDelta[kDeltaTypes] =         {kRed,      kBlue,       kRed,      kBlue,kBlack,    kBlue+4,     kRed+3,  kOrange+10,kOrange+10};
  const Int_t kFitDelta[kDeltaTypes] =    {          1,          1,          1,          1,     1,          1,          1,           0,         0};


  const Int_t kMarkerRP[nRP] = {kFullSquare,22,29,23}; // first no used
  const Int_t kColorRP[nRP] = {kBlack,kBlue+2,kCyan+2,kGreen+2};

  // 


  const Int_t kNormType = 7; // normaize this spectrum at fixed p_T
  const Float_t kNormPt = 10.; // normaize this spectrum at fixed p_T
  Float_t kNormValue[nMult][nCen] = {0,}; // normaize the spectrum with this value, if negative take from the first delta p_T histo
  Float_t kNormValueMult[nMult][nRP] = {0,}; // normaize the spectrum with this value, if negative take from the first delta p_T histo
  

  // 
  // the cuts for the single centrality
  Float_t fMinPt[kDeltaTypes] = {0};
  Float_t fMaxPt[kDeltaTypes] = {0};

  for(int iD = 0;iD<kDeltaTypes;iD++){
    fMinPt[iD] = inPtLo;
    fMaxPt[iD] = inPtUp;
  }
  
  TString tmpName;



  // 0 Leticias BiA anti-kT
  
  TFile *fL = TFile::Open("~/Dropbox/SharedJets/Leticia/randomcones/pwg4plots.root");
  TFile *fRP = TFile::Open("~/Dropbox/SharedJets/Leticia/randomcones/reactionplane.root");
  if(fDebug)Printf("Line: %d",__LINE__);
  Int_t iDelta = 0;
  Int_t iRP = 0;
  if(fL){
    iDelta = 0;
    // leticia BiA Randome cones
    sDelta[iDelta][0] = "BiA RC";
    if(iB==1)tmpName = "BiA sa RC skip0 centrality";
    else if (iB==2)tmpName = "BiA RC skip0 va centrality";
    TH2F *h2Tmp = (TH2F*)fL->Get(tmpName.Data());
    if(!h2Tmp)Printf("Line:%d %s not found",__LINE__,tmpName.Data());


    for(int ic = 0;ic < nCen;ic++){
      iRP = 0;
      if(!h2Tmp)continue;
      Int_t ibLo = h2Tmp->GetXaxis()->FindBin(fCentLo[ic]);
      Int_t ibUp = h2Tmp->GetXaxis()->FindBin(fCentUp[ic])-1;
      hDeltaPt[iDelta][ic][iRP] = h2Tmp->ProjectionY(Form("hBiaRC%d",ic),ibLo,ibUp,"E");
      Float_t fScale =  hDeltaPt[iDelta][ic][iRP]->Integral("width");
      if(fScale)  hDeltaPt[iDelta][ic][iRP]->Scale(1./fScale);
      if(fRP){
	if(ic==0){tmpName = "bia RC vs RP central";}
	else if(ic == 3){ tmpName = "bia RC vs RP perif";}
	else {continue;}
	TH2F *h2TmpRP = (TH2F*)fRP->Get(tmpName.Data());
	if(!h2TmpRP)Printf("Line:%d %s not found",__LINE__,tmpName.Data());
	for(iRP = 1;iRP<nRP;iRP++){
	  hDeltaPt[iDelta][ic][iRP] = h2TmpRP->ProjectionY(Form("hDeltaPt%d_%d_%d",iDelta,ic,iRP),iRP,iRP);
	  Float_t fScale =  hDeltaPt[iDelta][ic][iRP]->Integral("width");
	  if(fScale)  hDeltaPt[iDelta][ic][iRP]->Scale(1./fScale);
	}
      }
    }
    iRP = 0;


    if(fDebug)Printf("Line: %d",__LINE__);

    if(iB==1)tmpName = "BiA sa RC skip0 multiplicity";
    else if (iB==2)tmpName = "BiA RC skip0 va multiplicity";
    h2Tmp = (TH2F*)fL->Get(tmpName.Data());
    if(!h2Tmp)Printf("Line:%d %s not found",__LINE__,tmpName.Data());
    if(fDebug)Printf("Line: %d",__LINE__);
    for(int im = 0;im < nMult;im++){
      if(!h2Tmp)continue;
      Int_t ibLo = h2Tmp->GetXaxis()->FindBin(multMin[im]);
      Int_t ibUp = h2Tmp->GetXaxis()->FindBin(multMax[im])-1;
      Printf("Line:%d bin %d - %d",__LINE__,ibLo,ibUp);
      hDeltaPtMult[iDelta][im][iRP] = h2Tmp->ProjectionY(Form("hBia%d_M%d",iDelta,im),ibLo,ibUp,"E");
      Float_t fScale =  hDeltaPtMult[iDelta][im][iRP]->Integral("width");
      if(fScale)  hDeltaPtMult[iDelta][im][iRP]->Scale(1./fScale);
    }
    if(fDebug)Printf("Line: %d",__LINE__);
    //------------

    iDelta = 1;
    sDelta[iDelta][0] = "BiA RC (randomized event)";
    if(iB==1)tmpName = "BiA RC random input sa centrality";
    else if (iB==2)tmpName = "BiA RC random input va centrality";
    h2Tmp = (TH2F*)fL->Get(tmpName.Data());
    if(!h2Tmp)Printf("Line:%d %s not found",__LINE__,tmpName.Data());
    if(fDebug)Printf("Line: %d",__LINE__);
    for(int ic = 0;ic < nCen;ic++){
      if(!h2Tmp)continue;
      iRP = 0;
      Int_t ibLo = h2Tmp->GetXaxis()->FindBin(fCentLo[ic]);
      Int_t ibUp = h2Tmp->GetXaxis()->FindBin(fCentUp[ic])-1;
      Printf("Line:%d bin %d - %d",__LINE__,ibLo,ibUp);
      hDeltaPt[iDelta][ic][iRP] = h2Tmp->ProjectionY(Form("hBiaL%d",ic),ibLo,ibUp,"E");
      Float_t fScale =  hDeltaPt[iDelta][ic][iRP]->Integral("width");
      if(fScale)  hDeltaPt[iDelta][ic][iRP]->Scale(1./fScale);
      if(fRP){
	if(ic==0){tmpName = "bia RC randomized input vs RP central";}
	else if(ic == 3){ tmpName = "bia RC randomized input vs RP perif";}
	else {continue;}
	TH2F *h2TmpRP = (TH2F*)fRP->Get(tmpName.Data());
	if(!h2TmpRP)Printf("Line:%d %s not found",__LINE__,tmpName.Data());
	for(iRP = 1;iRP<nRP;iRP++){
	  hDeltaPt[iDelta][ic][iRP] = h2TmpRP->ProjectionY(Form("hDeltaPt%d_%d_%d",iDelta,ic,iRP),iRP,iRP);
	  Float_t fScale =  hDeltaPt[iDelta][ic][iRP]->Integral("width");
	  if(fScale)  hDeltaPt[iDelta][ic][iRP]->Scale(1./fScale);
	}
      }
    }
    iRP = 0;

    if(fDebug)Printf("Line: %d",__LINE__);
    if(iB==1)tmpName = "BiA RC random input sa multiplicity";
    else if (iB==2)tmpName = "BiA RC random input va multiplicity";
    h2Tmp = (TH2F*)fL->Get(tmpName.Data());
    if(!h2Tmp)Printf("Line:%d %s not found",__LINE__,tmpName.Data());
    for(int im = 0;im < nMult;im++){
      if(!h2Tmp)continue;
      Int_t ibLo = h2Tmp->GetXaxis()->FindBin(multMin[im]);
      Int_t ibUp = h2Tmp->GetXaxis()->FindBin(multMax[im])-1;
      Printf("Line:%d bin %d - %d",__LINE__,ibLo,ibUp);
      hDeltaPtMult[iDelta][im][iRP] = h2Tmp->ProjectionY(Form("hBia%d_M%d",iDelta,im),ibLo,ibUp,"E");
      Float_t fScale =  hDeltaPtMult[iDelta][im][iRP]->Integral("width");
      if(fScale)  hDeltaPtMult[iDelta][im][iRP]->Scale(1./fScale);
    }



    // ------------------------------------------
    // leticia random cones skip leading jets
    iDelta = 2;
    sDelta[iDelta][0] = "BiA RC (excl. 2 leading jets)";
    if(fDebug)Printf("Line: %d",__LINE__);
    if(iB==1)tmpName = "BiA sa RC skip2 centrality";
    else if (iB==2)tmpName = "BiA RC skip2 va centrality";
    h2Tmp = (TH2F*)fL->Get(tmpName.Data());
    if(!h2Tmp)Printf("Line:%d %s not found",__LINE__,tmpName.Data());      
    for(int ic = 0;ic < nCen;ic++){
      if(!h2Tmp)continue;

      Int_t ibLo = h2Tmp->GetXaxis()->FindBin(fCentLo[ic]);
      Int_t ibUp = h2Tmp->GetXaxis()->FindBin(fCentUp[ic])-1;
      Printf("Line:%d bin %d - %d",__LINE__,ibLo,ibUp);
      hDeltaPt[iDelta][ic][iRP] = h2Tmp->ProjectionY(Form("hBiaRC2%d",ic),ibLo,ibUp,"E");
      Float_t fScale =  hDeltaPt[iDelta][ic][iRP]->Integral("width");
      if(fScale) hDeltaPt[iDelta][ic][iRP]->Scale(1./fScale);
    }
    if(fDebug)Printf("Line: %d",__LINE__);
    if(iB==1)tmpName = "BiA RC skip2 sa multiplicity";
    else if (iB==2)tmpName = "BiA RC skip2 va multiplicity";
    h2Tmp = (TH2F*)fL->Get(tmpName.Data());
    if(!h2Tmp)Printf("Line:%d %s not found",__LINE__,tmpName.Data());
    for(int im = 0;im < nMult;im++){
      if(!h2Tmp)continue;
      Int_t ibLo = h2Tmp->GetXaxis()->FindBin(multMin[im]);
      Int_t ibUp = h2Tmp->GetXaxis()->FindBin(multMax[im])-1;
      Printf("Line:%d bin %d - %d",__LINE__,ibLo,ibUp);
      hDeltaPtMult[iDelta][im][iRP] = h2Tmp->ProjectionY(Form("hBia%d_M%d",iDelta,im),ibLo,ibUp,"E");
      Float_t fScale =  hDeltaPtMult[iDelta][im][iRP]->Integral("width");
      if(fScale)  hDeltaPtMult[iDelta][im][iRP]->Scale(1./fScale);
    }
  }  

  // fetch the data from bastian...


  // jet embedded
  TString cBB = "/Users/kleinb/Dropbox/SharedJets/Bastian/Files/"; 

  iDelta = 3;
  sDelta[iDelta][0] = Form("anti-k_{T} embedded jet %2.0f-%2.0f GeV",fMinPt[iDelta],fMaxPt[iDelta]);
  for(iRP = 0;iRP<nRP;iRP++){
    for(int ic = 0;ic < nCen;ic++){
      TH2F *hTmp = GetTH2PlotB(cBB.Data(),1,0,ic,iRP-1); // emb jets
      if(!hTmp)Printf("%d %s not found",__LINE__,tmpName.Data());
      int ibLo = hTmp->GetYaxis()->FindBin(fMinPt[iDelta]);
      int ibUp = hTmp->GetYaxis()->FindBin(fMaxPt[iDelta])-1;
      hDeltaPt[iDelta][ic][iRP] = hTmp->ProjectionX(Form("fHistDeltaPtB1_c%d_rp%d",ic,iRP),ibLo,ibUp,"E");
      hDeltaPt[iDelta][ic][iRP]->Rebin(2);
      Float_t fScale = hDeltaPt[iDelta][ic][iRP]->Integral("width");
      if(fScale)  hDeltaPt[iDelta][ic][iRP]->Scale(1./fScale);
    }
  }
  iRP = 0;
  for(iRP = 0;iRP<nRP;iRP++){
    for(int im = 0;im <nMult;im++){
      TH2F *hTmp = GetTH2PlotB(cBB.Data(),1,1,im,iRP-1); // emb jets
      if(!hTmp)Printf("%d %s not found",__LINE__,tmpName.Data());
      int ibLo = hTmp->GetYaxis()->FindBin(fMinPt[iDelta]);
      int ibUp = hTmp->GetYaxis()->FindBin(fMaxPt[iDelta])-1;
      Printf("Line:%d bin %d - %d",__LINE__,ibLo,ibUp);
      hDeltaPtMult[iDelta][im][iRP] = hTmp->ProjectionX(Form("fHistDeltaPtB2_%d_mult%d_rp%d",iDelta,im,iRP),ibLo,ibUp,"E");
      hDeltaPtMult[iDelta][im][iRP]->Rebin(2);
      Float_t fScale =   hDeltaPtMult[iDelta][im][iRP]->Integral("width");
      if(fScale)   hDeltaPtMult[iDelta][im][iRP]->Scale(1./fScale);
    }
  }
  iRP = 0;

  iDelta = 4;
  sDelta[iDelta][0] = Form("anti-k_{T} embedded track %2.0f-%2.0f GeV BB",fMinPt[iDelta],fMaxPt[iDelta]);
  
  for(iRP = 0;iRP<nRP;iRP++){
    for(int ic = 0;ic < nCen;ic++){
      TH2F *hTmp = GetTH2PlotB(cBB.Data(),0,0,ic,iRP-1); // emb jets
      if(!hTmp)Printf("%d %s not found",__LINE__,tmpName.Data());
      int ibLo = hTmp->GetYaxis()->FindBin(fMinPt[iDelta]);
      int ibUp = hTmp->GetYaxis()->FindBin(fMaxPt[iDelta])-1;
      hDeltaPt[iDelta][ic][iRP] = hTmp->ProjectionX(Form("fHistDeltaPtB2_c%d_rp%d",ic,iRP),ibLo,ibUp,"E");
      hDeltaPt[iDelta][ic][iRP]->Rebin(2);
      Float_t fScale =   hDeltaPt[iDelta][ic][iRP]->Integral("width");
      if(fScale)   hDeltaPt[iDelta][ic][iRP]->Scale(1./fScale);
    }
  }
  iRP = 0;
  for(iRP = 0;iRP<nRP;iRP++){
    for(int im = 0;im <nMult;im++){
      TH2F *hTmp = GetTH2PlotB(cBB.Data(),0,1,im,iRP-1); // emb jets
      if(!hTmp)Printf("%d %s not found",__LINE__,tmpName.Data());
      int ibLo = hTmp->GetYaxis()->FindBin(fMinPt[iDelta]);
      int ibUp = hTmp->GetYaxis()->FindBin(fMaxPt[iDelta])-1;
      Printf("Line:%d bin %d - %d",__LINE__,ibLo,ibUp);
      hDeltaPtMult[iDelta][im][iRP] = hTmp->ProjectionX(Form("fHistDeltaPtB2_%d_mult%d_rp%d",iDelta,im,iRP),ibLo,ibUp,"E");
      hDeltaPtMult[iDelta][im][iRP]->Rebin(2);
      Float_t fScale =   hDeltaPtMult[iDelta][im][iRP]->Integral("width");
      if(fScale)   hDeltaPtMult[iDelta][im][iRP]->Scale(1./fScale);
    }
  }
  iRP = 0;
    // Data from marta
    iDelta = 5;
    TString fOTFQoff = "~/Dropbox/SharedJets/OnTheFlyOutput/PWG4_JetTasksOutput_110422_p2_OTFQOff_000139107_Total.root";
    TFile *fQoff = TFile::Open(fOTFQoff.Data());
    sDelta[iDelta][0] = Form("anti-k_{T} embedded unquenched jet %2.0f-%2.0f GeV",fMinPt[iDelta],fMaxPt[iDelta]);
    if(fQoff){
      for(int ic = 0;ic < nCen;ic++){
      tmpName = Form("PWG4_BkgFluct%sCent%dB%d","ANTIKT",ic,iB);
      TDirectory *dir = (TDirectory*)fQoff->Get(tmpName.Data());
      if(!dir)Printf("Line:%d %s not found",__LINE__,tmpName.Data());
      tmpName = Form("taskBkgFluct%sCent%dB%d","ANTIKT",ic,iB);
      TList *list = (TList*)dir->Get(tmpName.Data());
      if(!list)Printf("Line:%d %s not found",__LINE__,tmpName.Data());
      TH3F *h3Tmp = (TH3F*)list->FindObject("fh3PtDeltaPtArea");
      if(!h3Tmp)Printf("Line:%d %s not found",__LINE__,"fh3PtDeltaPtArea");
      hDeltaPt[iDelta][ic][iRP] = h3Tmp->ProjectionY(Form("hDeltaM%d_%d",ic,iDelta),h3Tmp->GetXaxis()->FindBin(fMinPt[iDelta]),h3Tmp->GetXaxis()->FindBin(fMaxPt[iDelta]),
						0,-1,"E");
      if(ic<2)  hDeltaPt[iDelta][ic][iRP]->Rebin(4);
      else if(ic==2)  hDeltaPt[iDelta][ic][iRP]->Rebin(2);
      Float_t fScale =  hDeltaPt[iDelta][ic][iRP]->Integral("width");
      if(fScale) hDeltaPt[iDelta][ic][iRP]->Scale(1./fScale);
      Printf("Line:%d %p",__LINE__,hDeltaPt[iDelta][ic][iRP]);
      }
    }
    else{
      Printf("Could not open %s",fOTFQoff.Data());
    }

    // marta quenched...
    iDelta = 6;
    TString fOTFQon = "~/Dropbox/SharedJets/OnTheFlyOutput/PWG4_JetTasksOutput_110422_p2_OTFQOn_000139107.root";
    TFile *fQon = TFile::Open(fOTFQon.Data());
    sDelta[iDelta][0] = Form("anti-k_{T} embedded quenched jet %2.0f-%2.0f GeV",fMinPt[iDelta],fMaxPt[iDelta]);
    if(fQon){
      for(int ic = 0;ic < nCen;ic++){
	tmpName = Form("PWG4_BkgFluct%sCent%dB%d","ANTIKT",ic,iB);
	TDirectory *dir = (TDirectory*)fQon->Get(tmpName.Data());
	if(!dir)Printf("Line:%d %s not found",__LINE__,tmpName.Data());
	tmpName = Form("taskBkgFluct%sCent%dB%d","ANTIKT",ic,iB);
	TList *list = (TList*)dir->Get(tmpName.Data());
	if(!list)Printf("Line:%d %s not found",__LINE__,tmpName.Data());
	
	TH3F *h3Tmp = (TH3F*)list->FindObject("fh3PtDeltaPtArea");
	hDeltaPt[iDelta][ic][iRP] = h3Tmp->ProjectionY(Form("hDeltaM%d_%d",ic,iDelta),h3Tmp->GetXaxis()->FindBin(fMinPt[iDelta]),h3Tmp->GetXaxis()->FindBin(fMaxPt[iDelta]),
					     0,-1,"E");
	if(ic<2)  hDeltaPt[iDelta][ic][iRP]->Rebin(4);
	else if(ic==2)  hDeltaPt[iDelta][ic][iRP]->Rebin(2);
	Float_t fScale = hDeltaPt[iDelta][ic][iRP]->Integral("width");
	if(fScale)hDeltaPt[iDelta][ic][iRP]->Scale(1./fScale);
      }
    }
    else{
      Printf("Could not open %s",fOTFQon.Data());
    }

    TString sinputFile = "~/Dropbox/SharedJets/Christian/Files/PWG4_JetTasksOutput_LHC10h_AOD_tmp.root";
    TString sinputDir = "PWG4_spec2_clustersAOD_ANTIKT04_B2_Filter00128_Cut00150_Skip02__0000000000_Class00";
    // Plot the real jet spectrum, to be normalized later...
    iDelta = 7;
    sDelta[iDelta][0] = "anti k_{T} scaled jet spectrum"; 
    
    if(true){
      THnSparseF* fhnJetPtRec =  (THnSparseF*)GetObjectFromFile("fhnJetPtRec",sinputFile.Data(),sinputDir.Data());
      fhnJetPtRec->GetAxis(0)->SetRange(3,3); // take all jets

      // in centrality bins....
      for(iRP = 0;iRP<nRP;iRP++){
	for(int ic = 0;ic < nCen;ic++){
	  Int_t icLo = fhnJetPtRec->GetAxis(2)->FindBin(fCentLo[ic]+0.1);
	  Int_t icUp = fhnJetPtRec->GetAxis(2)->FindBin(fCentUp[ic]-0.1);
	  fhnJetPtRec->GetAxis(2)->SetRange(icLo,icUp);

	  // iRP range 0 does rese , this is what we want for all :)
	  fhnJetPtRec->GetAxis(4)->SetRange(iRP,iRP);

	  hDeltaPt[iDelta][ic][iRP] = fhnJetPtRec->Projection(1,"E");
	  hDeltaPt[iDelta][ic][iRP]->SetName(Form("hSpectrumJet_%d_C%d_rp%d",iDelta,ic,iRP));
	}
	fhnJetPtRec->GetAxis(2)->SetRange(0,0);// reset centrality
	for(int im = 0;im < nMult;im++){
	  Int_t ibTLo =   fhnJetPtRec->GetAxis(3)->FindBin(multMin[im]);
	  Int_t ibTUp =   fhnJetPtRec->GetAxis(3)->FindBin(multMax[im]);

	  fhnJetPtRec->GetAxis(3)->SetRange(ibTLo,ibTUp);	  
	   // iRP range 0 does rese , this is what we want for all :)
	  fhnJetPtRec->GetAxis(4)->SetRange(iRP,iRP);

	  hDeltaPtMult[iDelta][im][iRP] = fhnJetPtRec->Projection(1,"E");
	  hDeltaPtMult[iDelta][im][iRP]->SetName(Form("hSpectrumJets_%d_M%d_rp%d",iDelta,im,iRP));
	}
      }
    }
    iRP = 0;


      // fetch the poissonian fluctuations! CKB
      iDelta = 8;
      sDelta[iDelta][0] = "Poissonian limit";
      if(true){
	THnSparseF* fhnTrackPtRec =  (THnSparseF*)GetObjectFromFile("fhnTrackPtRec",sinputFile.Data(),sinputDir.Data());
	fhnTrackPtRec->GetAxis(0)->SetRange(2,2); // take all tracks

      THnSparseF* fhnEvent =  (THnSparseF*)GetObjectFromFile("fhnEvent",sinputFile.Data(),sinputDir.Data());
      TH1F *fh1Centrality = (TH1F*)GetObjectFromFile("fh1Centrality",sinputFile.Data(),sinputDir.Data());
      // Multiplicity
      TH2F *fh2MultRec = (TH2F*)GetObjectFromFile("fh2MultRec",sinputFile.Data(),sinputDir.Data());
      TH1D *h1Mult = (TH1D*)fh2MultRec->ProjectionX("h1Mult");
      

      // in centrality bins....

      for(int ic = 0;ic < nCen;ic++){
	Int_t icLo = fhnTrackPtRec->GetAxis(2)->FindBin(fCentLo[ic]+0.1);
	Int_t icUp = fhnTrackPtRec->GetAxis(2)->FindBin(fCentUp[ic]-0.1);
	fhnTrackPtRec->GetAxis(2)->SetRange(icLo,icUp);

	TH1D *hMult = fhnEvent->Projection(1,"E");
	if(!grSigmaDeltaPtCent[iDelta][iRP]){
	  grSigmaDeltaPtCent[iDelta][iRP] = new TGraphErrors(nCen);
	  SetGraphAttributes(grSigmaDeltaPtCent[iDelta][iRP],kMarkerDelta[iDelta],kColorDelta[iDelta]);
	}

	TH1D *hSpectrumTracks = fhnTrackPtRec->Projection(1,"E");
	hSpectrumTracks->SetName(Form("hSpectrumTracks_C%d",ic));
	
	// scale with bin width?
	
	for(int ib = 1;ib<hSpectrumTracks->GetNbinsX();ib++){
	  Float_t val = hSpectrumTracks->GetBinContent(ib);
	  hSpectrumTracks->SetBinContent(ib,val/hSpectrumTracks->GetBinWidth(ib));
	}
	
	int ibLo = fh1Centrality->FindBin(fCentLo[ic]+0.1);
	int ibUp = fh1Centrality->FindBin(fCentUp[ic]-0.1);
	Float_t fScale = fh1Centrality->Integral(ibLo,ibUp);
	if(fScale>0)hSpectrumTracks->Scale(1./fScale);	
	Double_t sigma = GetPoissonFluctuation(hSpectrumTracks,TMath::Pi()*2.*1.8,
					       TMath::Pi()*0.4*0.4);
	Double_t sigma_error = sigma/1000.;
	Double_t mean = hSpectrumTracks->GetMean();
	Double_t rms = hSpectrumTracks->GetRMS();
	Double_t mult = hSpectrumTracks->Integral("width");
	Double_t cent = (fCentUp[ic]+fCentLo[ic])/2;
	Double_t cent_e = (fCentUp[ic]-fCentLo[ic])/2;
	grSigmaDeltaPtCent[iDelta][iRP]->SetPoint(ic,cent,sigma);
	grSigmaDeltaPtCent[iDelta][iRP]->SetPointError(ic,cent_e,sigma_error);
	delete hSpectrumTracks;
	delete hMult;
      }
      fhnTrackPtRec->GetAxis(2)->SetRange(0,0);// reset
      
      // in multiplicity bins...
      for(iRP = 0;iRP<nRP;iRP++){
	for(int im = 0;im < nMult;im++){
	  Int_t ibELo = fhnEvent->GetAxis(1)->FindBin(multMin[im]);
	  Int_t ibEUp = fhnEvent->GetAxis(1)->FindBin(multMax[im]);
	  fhnEvent->GetAxis(1)->SetRange(ibELo,ibEUp);
	  
	  TH1D *hMult = fhnEvent->Projection(1,"E");
	  
	  if(!grSigmaDeltaPtMult[iDelta][iRP]){
	    grSigmaDeltaPtMult[iDelta][iRP] = new TGraphErrors(nCen);
	    SetGraphAttributes(grSigmaDeltaPtMult[iDelta][iRP],kMarkerDelta[iDelta],kColorDelta[iDelta]);
	    if(iRP>0)SetGraphAttributes(grSigmaDeltaPtMult[iDelta][iRP],kMarkerRP[iRP],kColorDelta[iDelta]); // change only the marke here
	  }

	  // mult range
	  Int_t ibTLo =   fhnTrackPtRec->GetAxis(3)->FindBin(multMin[im]);
	  Int_t ibTUp =   fhnTrackPtRec->GetAxis(3)->FindBin(multMax[im]);
	  fhnTrackPtRec->GetAxis(3)->SetRange(ibTLo,ibTUp);

	  // iRP range 0 does rese , this is what we wantfor all :)
	  fhnTrackPtRec->GetAxis(4)->SetRange(iRP,iRP);
	  // control
	  /*
	  TH2D *h2Tmp = fhnTrackPtRec->Projection(1,4,"E");
	  h2Tmp->Draw("colz");
	  c1->Update();
	  if(getchar()=='q')return 1;
	  */
	  TH1D *hSpectrumTracks = fhnTrackPtRec->Projection(1,"E");
	  hSpectrumTracks->SetName(Form("hSpectrumTracks_M%d_rp%d",im,iRP));
	  
	  // scale with bin width?
	  for(int ib = 1;ib<hSpectrumTracks->GetNbinsX();ib++){
	    Float_t val = hSpectrumTracks->GetBinContent(ib);
	    hSpectrumTracks->SetBinContent(ib,val/hSpectrumTracks->GetBinWidth(ib));
	  }
	  
	  int ibLo = h1Mult->FindBin(multMin[im]);
	  int ibUp = h1Mult->FindBin(multMax[im]);
	  h1Mult->GetXaxis()->SetRange(ibLo,ibUp);
	  Float_t fScale = h1Mult->Integral(ibLo,ibUp);
	  if(fScale>0)hSpectrumTracks->Scale(1./fScale);	

	  Double_t sigma = 0;
	  Double_t sigma_error = 0;
	  Double_t mult = 0;
	  
	  Double_t mult2 = h1Mult->GetMean(1);
	  Double_t mult2_e = h1Mult->GetMean(11); // should be the error on the mult :)

	  if(iRP==0){// all
	   sigma = GetPoissonFluctuation(hSpectrumTracks,TMath::Pi()*2.*1.8,
				  TMath::Pi()*0.4*0.4);
	   sigma_error = sigma/1000.;
	   mult = hSpectrumTracks->Integral("width");
	  }
	  else{
	    // uses only 1/3 of acceptance
	    sigma = GetPoissonFluctuation(hSpectrumTracks,TMath::Pi()*2.*1.8/Float_t(nRP-1),
				  TMath::Pi()*0.4*0.4);
	    sigma_error = sigma/1000.;
	    mult = hSpectrumTracks->Integral("width");
	  }
	  Printf("mult %4.3f mult2 %4.3f",mult,mult2);
	  
	  grSigmaDeltaPtMult[iDelta][iRP]->SetPoint(im,mult2,sigma);
	  grSigmaDeltaPtMult[iDelta][iRP]->SetPointError(im,mult2_e,sigma_error);
	  delete hSpectrumTracks;
	  delete hMult;
	}
      }
      // in multiplicity bins...
	

    }

    for(iDelta = 1;iDelta < kDeltaTypes;iDelta++){
      sDelta[iDelta][1] = sDelta[iDelta][0] + " in plane";
      sDelta[iDelta][2] = sDelta[iDelta][0] + " in between";
      sDelta[iDelta][3] = sDelta[iDelta][0] + " out of plane";
    }

    // THIS IS FOR THE ALICE LABELS AND THE WORK IN PROGRESS PRELIMINARY ETC.
    TLatex *txt2 = new TLatex();
    txt2->SetTextFont(gStyle->GetTextFont());
    txt2->SetTextSize(gStyle->GetTextSize()*0.7);
    txt2->SetTextAlign(22);
    txt2->SetTextColor(kRed);
    

    
    if(iPlotType&(1<<0)){
      for(int ic = 0;ic < nCen;ic++){
	TLegend *leg1 = new TLegend(0.2,0.7,0.3,0.98);
	leg1->SetHeader(Form("Pb+Pb %2.0f-%2.0f%% R = 0.4 (B%d)",fCentLo[ic],fCentUp[ic],iB));
	leg1->SetFillColor(0);
	leg1->SetTextFont(gStyle->GetTextFont());
	leg1->SetTextSize(gStyle->GetTextSize()*0.6);
	leg1->SetBorderSize(0);
	hFrame->DrawCopy();

	for(iRP = 0;iRP<((iPlotType&(1<<2))?nRP:1);iRP++){
	  Double_t mean = 0;
	  Double_t sigma = 0;
	  Double_t sigma_error = 0;
	  Double_t mean_error = 0;
	  TF1* tmpGaus = 0;
	  for(iDelta = 0;iDelta<kDeltaTypes;iDelta++){
	    if(!hDeltaPt[iDelta][ic][iRP]){
	      Printf("%d:%d:%d not found",iDelta,ic,iRP);
	      continue;
	    }
	    SetHistoAttributes(hDeltaPt[iDelta][ic][iRP],kMarkerDelta[iDelta],kColorDelta[iDelta]);
	    if(iRP!=0)SetHistoAttributes(hDeltaPt[iDelta][ic][iRP],kMarkerRP[iRP],kColorRP[iRP]);
	    if(!(iPlotFlag&(1<<iDelta)))continue;

	    if(!grMeanDeltaPtCent[iDelta][iRP]&&kFitDelta[iDelta]){
	      grMeanDeltaPtCent[iDelta][iRP] = new TGraphErrors(nCen);
	      SetGraphAttributes(grMeanDeltaPtCent[iDelta][iRP],kMarkerDelta[iDelta],kColorDelta[iDelta]);
	      if(iRP>0)SetGraphAttributes(grMeanDeltaPtCent[iDelta][iRP],kMarkerRP[iRP],kColorDelta[iRP]);
	    }
	    if(!grSigmaDeltaPtCent[iDelta][iRP]&&kFitDelta[iDelta]){
	      grSigmaDeltaPtCent[iDelta][iRP] = new TGraphErrors(nCen);
	      SetGraphAttributes(grSigmaDeltaPtCent[iDelta][iRP],kMarkerDelta[iDelta],kColorDelta[iDelta]);
	      if(iRP>0)SetGraphAttributes(grSigmaDeltaPtCent[iDelta][iRP],kMarkerRP[iRP],kColorDelta[iRP]);
	    }

	    // take the first iDelta Pt as anchor for norm, worst case does nothin (scale with 1
	    if(kNormValue[ic][iRP]<=0){
	      Int_t ib = hDeltaPt[iDelta][ic][iRP]->FindBin(kNormPt);
	      kNormValue[ic][iRP] = hDeltaPt[iDelta][ic][iRP]->GetBinContent(ib);
	    }
	    if(iDelta==kNormType){
	      Int_t ib = hDeltaPt[iDelta][ic][iRP]->FindBin(kNormPt);
	      Float_t val1 = hDeltaPt[iDelta][ic][iRP]->GetBinContent(ib);
	      if(val1){
		hDeltaPt[iDelta][ic][iRP]->Scale(kNormValue[ic][iRP]/val1);
	      }
	    }


	    hDeltaPt[iDelta][ic][iRP]->DrawCopy("psame");
	    leg1->AddEntry(hDeltaPt[iDelta][ic][iRP],sDelta[iDelta][iRP].Data(),"P");
	    if(kFitDelta[iDelta]){
	      tmpGaus = FitLHSgaus(hDeltaPt[iDelta][ic][iRP]);
	      mean = tmpGaus->GetParameter(1);
	      sigma = tmpGaus->GetParameter(2);
	      mean_error = tmpGaus->GetParError(1);
	      sigma_error = tmpGaus->GetParError(2);	
	      
	      tmpGaus->SetRange(-40,40);
	      tmpGaus->SetLineColor(kColorDelta[iDelta]);
	      tmpGaus->Draw("same");
	      leg1->AddEntry(tmpGaus,Form("LHS Gaus fit: #mu = %2.2f, #sigma = %2.2f",mean,sigma),"L");
	      Double_t cent = (fCentUp[ic]+fCentLo[ic])/2;
	      Double_t cent_e = (fCentUp[ic]-fCentLo[ic])/2;
	      Printf("cent %3.3f +- %3.3f",cent,cent_e);
	      grMeanDeltaPtCent[iDelta][iRP]->SetPoint(ic,cent,mean);
	      grMeanDeltaPtCent[iDelta][iRP]->SetPointError(ic,cent_e,mean_error);
	      
	      grSigmaDeltaPtCent[iDelta][iRP]->SetPoint(ic,cent,sigma);
	      grSigmaDeltaPtCent[iDelta][iRP]->SetPointError(ic,cent_e,sigma_error);
	    }
	  } //deltatypes
	}// RP
	leg1->Draw();
	
	txt2->SetNDC();
	//    txt2->DrawLatex(0.5,0.97,"ALICE Performance 01/03/2011");
	c1->Update();
	c1->SaveAs(Form("%sdeltaPt_pT%s%s_B%d_cen%02d_%03d.%s",
			picPrefix.Data(),
			Form("%03.0f_%03.0f",inPtLo,inPtUp),
			(iPlotType&(1<<2))?"_rp":"",
			iB,ic,iPlotFlag,
			picSuffix.Data()));
	if(!gROOT->IsBatch()){
	  if(getchar()=='q')return 1;
	}
      }

      // Draw the trending plots vs cent...
      c1->SetLogy(0);    
      TH2F* hFrameMeanCent = new TH2F("hMeanCent",";centrality (%);#mu of LHS Gaus",1,0.,100.,1,-10.,10.);
      hFrameMeanCent->SetLabelSize(hFrame->GetLabelSize("Y")*0.9,"Y");
      hFrameMeanCent->SetLabelSize(hFrame->GetLabelSize("X")*0.9,"X");
      hFrameMeanCent->SetTitleSize(hFrame->GetTitleSize("Y")*0.9,"Y");
      hFrameMeanCent->SetTitleSize(hFrame->GetTitleSize("X")*0.9,"X");
      hFrameMeanCent->SetTitleOffset(1.1,"Y");
      hFrameMeanCent->SetTitleOffset(1.1,"X");
      hFrameMeanCent->SetStats(kFALSE);
      hFrameMeanCent->DrawCopy();
      
      TLegend *legSM = new TLegend(0.2,0.7,0.3,0.98);
      legSM->SetHeader(Form("Pb+Pb R = 0.4 (B%d)",iB));
      legSM->SetFillColor(0);
      legSM->SetTextFont(gStyle->GetTextFont());
      legSM->SetTextSize(gStyle->GetTextSize()*0.6);
      legSM->SetBorderSize(0);
      
      
      for(iRP = 0;iRP<((iPlotType&(1<<2))?nRP:1);iRP++){
	for(iDelta = 0;iDelta <kDeltaTypes;iDelta++){
	  if(!(iPlotFlag&(1<<iDelta)))continue;
	  if(!grMeanDeltaPtCent[iDelta][iRP])continue;
	  grMeanDeltaPtCent[iDelta][iRP]->Draw("psame");
	  legSM->AddEntry(grMeanDeltaPtCent[iDelta][iRP],sDelta[iDelta][iRP].Data(),"P");
	}
      }
      legSM->Draw();
      c1->Update();
      c1->SaveAs(Form("%sMeanVsCent_pT%s%s_B%d_%03d.%s",
		      picPrefix.Data(),
		      Form("%03.0f_%03.0f",inPtLo,inPtUp),
		      (iPlotType&(1<<2))?"_rp":"",
		      iB,iPlotFlag,picSuffix.Data()));
	if(!gROOT->IsBatch()){
      if(getchar()=='q')return 1;
	}
      TH2F* hFrameSigmaCent = new TH2F("hSigmaCent",";centrality (%);#sigma of LHS Gaus",1,0.,100.,10,0.,18.);
      hFrameSigmaCent->SetLabelSize(hFrame->GetLabelSize("Y")*0.9,"Y");
      hFrameSigmaCent->SetLabelSize(hFrame->GetLabelSize("X")*0.9,"X");
      hFrameSigmaCent->SetTitleSize(hFrame->GetTitleSize("Y")*0.9,"Y");
      hFrameSigmaCent->SetTitleSize(hFrame->GetTitleSize("X")*0.9,"X");
      hFrameSigmaCent->SetTitleOffset(1.1,"Y");
      hFrameSigmaCent->SetTitleOffset(1.1,"X");
      hFrameSigmaCent->SetStats(kFALSE);
      hFrameSigmaCent->DrawCopy();
      
      TLegend *legSC = new TLegend(0.2,0.7,0.3,0.98);
      legSC->SetHeader(Form("Pb+Pb R = 0.4 (B%d)",iB));
      legSC->SetFillColor(0);
      legSC->SetTextFont(gStyle->GetTextFont());
      legSC->SetTextSize(gStyle->GetTextSize()*0.6);
      legSC->SetBorderSize(0);
      
      for(iRP = 0;iRP<((iPlotType&(1<<2))?nRP:1);iRP++){
	for(iDelta = 0;iDelta <kDeltaTypes;iDelta++){
	  if(!(iPlotFlag&(1<<iDelta)))continue;
	  if(!grSigmaDeltaPtCent[iDelta][iRP])continue;
	  grSigmaDeltaPtCent[iDelta][iRP]->Draw("psame");
	  legSC->AddEntry(grSigmaDeltaPtCent[iDelta][iRP],sDelta[iDelta][iRP].Data(),"P");
	}
      }
      legSC->Draw();
      c1->Update();


      c1->SaveAs(Form("%sSigmaVsCent_pT%s%s_B%d_%03d.%s",
		      picPrefix.Data(),
		      Form("%03.0f_%03.0f",inPtLo,inPtUp),
		      (iPlotType&(1<<2))?"_rp":"",
		      iB,iPlotFlag,picSuffix.Data()));

      if(!gROOT->IsBatch()){
	if(getchar()=='q')return 1;
      }
    }

    // plot the trends vs mult
    if(iPlotType&(1<<1)){
      for(int im = 0;im < nMult;im++){
	TLegend *leg1 = new TLegend(0.2,0.7,0.3,0.98);
	leg1->SetHeader(Form("Pb+Pb mult. %03d-%03d R = 0.4 (B%d)",multMin[im],multMax[im],iB));
	leg1->SetFillColor(0);
	leg1->SetTextFont(gStyle->GetTextFont());
	leg1->SetTextSize(gStyle->GetTextSize()*0.6);
	leg1->SetBorderSize(0);
	hFrame->DrawCopy();
	
	for(iRP = 0;iRP<((iPlotType&(1<<2))?nRP:1);iRP++){
	  Double_t mean = 0;
	  Double_t sigma = 0;
	  Double_t sigma_error = 0;
	  Double_t mean_error = 0;
	  TF1* tmpGaus = 0;
	  for(iDelta = 0;iDelta<kDeltaTypes;iDelta++){
	    if(!hDeltaPtMult[iDelta][im][iRP]){
	      Printf("%d:%d:%d not found",iDelta,im,iRP);
	      continue;
	    }
	    SetHistoAttributes(hDeltaPtMult[iDelta][im][iRP],kMarkerDelta[iDelta],kColorDelta[iDelta]);
	    if(iRP>0)SetHistoAttributes(hDeltaPtMult[iDelta][im][iRP],kMarkerRP[iRP],kColorRP[iRP]);
	    if(!(iPlotFlag&(1<<iDelta)))continue;
	    
	    if(!grMeanDeltaPtMult[iDelta][iRP]&&kFitDelta[iDelta]){
	      grMeanDeltaPtMult[iDelta][iRP] = new TGraphErrors(nMult);
	      SetGraphAttributes(grMeanDeltaPtMult[iDelta][iRP],kMarkerDelta[iDelta],kColorDelta[iDelta]);
	      if(iRP>0)SetGraphAttributes(grMeanDeltaPtMult[iDelta][iRP],kMarkerRP[iRP],kColorRP[iRP]);
	    }
	    if(!grSigmaDeltaPtMult[iDelta][iRP]&&kFitDelta[iDelta]){
	      grSigmaDeltaPtMult[iDelta][iRP] = new TGraphErrors(nMult);
	      SetGraphAttributes(grSigmaDeltaPtMult[iDelta][iRP],kMarkerDelta[iDelta],kColorDelta[iDelta]);
	      if(iRP>0)SetGraphAttributes(grSigmaDeltaPtMult[iDelta][iRP],kMarkerRP[iRP],kColorRP[iRP]);

	    }

	    // take the first iDelta Pt as anchor for norm, worst case does nothin (scale with 1
	    if(kNormValueMult[im][iRP]<=0){
	      Int_t ib = hDeltaPtMult[iDelta][im][iRP]->FindBin(kNormPt);
	      kNormValueMult[im][iRP] = hDeltaPtMult[iDelta][im][iRP]->GetBinContent(ib);
	    }


	    if(iDelta==kNormType){
	      Int_t ib = hDeltaPtMult[iDelta][im][iRP]->FindBin(kNormPt);
	      Float_t val1 = hDeltaPtMult[iDelta][im][iRP]->GetBinContent(ib);
	      if(val1!=0){
		hDeltaPtMult[iDelta][im][iRP]->Scale(kNormValueMult[im][iRP]/val1);
	      }
	    }

	    hDeltaPtMult[iDelta][im][iRP]->DrawCopy("psame");
	    leg1->AddEntry(hDeltaPtMult[iDelta][im][iRP],sDelta[iDelta][iRP].Data(),"P");
	    if(kFitDelta[iDelta]){
	      tmpGaus = FitLHSgaus(hDeltaPtMult[iDelta][im][iRP]);
	      mean = tmpGaus->GetParameter(1);
	      sigma = tmpGaus->GetParameter(2);
	      mean_error = tmpGaus->GetParError(1);
	      sigma_error = tmpGaus->GetParError(2);	
	      
	      tmpGaus->SetRange(-40,40);
	      tmpGaus->SetLineColor(hDeltaPtMult[iDelta][im][iRP]->GetLineColor());
	      tmpGaus->Draw("same");
	      leg1->AddEntry(tmpGaus,Form("LHS Gaus fit: #mu = %2.2f, #sigma = %2.2f",mean,sigma),"L");
	      
	      // to be replaced by actual mean of the mutliplicity
	      Double_t mult = (multMin[im]+multMax[im])/2;
	      Double_t mult_e = (multMax[im]-multMin[im])/2;
	      Printf("mult %3.3f +- %3.3f",mult,mult_e);
	      grMeanDeltaPtMult[iDelta][iRP]->SetPoint(im,mult,mean);
	      grMeanDeltaPtMult[iDelta][iRP]->SetPointError(im,mult_e,mean_error);
	      
	      grSigmaDeltaPtMult[iDelta][iRP]->SetPoint(im,mult,sigma);
	      grSigmaDeltaPtMult[iDelta][iRP]->SetPointError(im,mult_e,sigma_error);
	    }
	  }
	}
	leg1->Draw();
	txt2->SetNDC();
	//    txt2->DrawLatex(0.5,0.97,"ALICE Performance 01/03/2011");
	c1->Update();
	c1->SaveAs(Form("%sdeltaPt_pT%s%s_B%d_mult%02d_%03d.%s",
			picPrefix.Data(),
			Form("%03.0f_%03.0f",inPtLo,inPtUp),
			(iPlotType&(1<<2))?"_rp":"",
			iB,im,iPlotFlag,
			picSuffix.Data()));
	if(!gROOT->IsBatch()){
	  if(getchar()=='q')return 1;	
	}
      }
	// Draw the trending plots vs cent...
	c1->SetLogy(0);    
	TH2F* hFrameMeanMult = new TH2F("hMeanMult",";raw #tracks;#mu of LHS Gaus",500,0.,3200.,12,-10.,10.);
	hFrameMeanMult->SetLabelSize(hFrame->GetLabelSize("Y")*0.9,"Y");
	hFrameMeanMult->SetLabelSize(hFrame->GetLabelSize("X")*0.9,"X");
	hFrameMeanMult->SetTitleSize(hFrame->GetTitleSize("Y")*0.9,"Y");
	hFrameMeanMult->SetTitleSize(hFrame->GetTitleSize("X")*0.9,"X");
	hFrameMeanMult->SetTitleOffset(1.1,"Y");
	hFrameMeanMult->SetTitleOffset(1.1,"X");
	hFrameMeanMult->SetStats(kFALSE);
	hFrameMeanMult->DrawCopy();
	
	TLegend *legMM = new TLegend(0.2,0.7,0.3,0.98);
	legMM->SetHeader(Form("Pb+Pb R = 0.4 (B%d)",iB));
	legMM->SetFillColor(0);
	legMM->SetTextFont(gStyle->GetTextFont());
	legMM->SetTextSize(gStyle->GetTextSize()*0.6);
	legMM->SetBorderSize(0);
	
	
	for(iRP = 0;iRP<((iPlotType&(1<<2))?nRP:1);iRP++){
	  for(iDelta = 0;iDelta <kDeltaTypes;iDelta++){
	    if(!(iPlotFlag&(1<<iDelta)))continue;
	    if(!grMeanDeltaPtMult[iDelta][iRP])continue;
	    grMeanDeltaPtMult[iDelta][iRP]->Draw("psame");
	    legMM->AddEntry(grSigmaDeltaPtMult[iDelta][iRP],sDelta[iDelta][iRP].Data(),"P");
	  }
	}
	legMM->Draw();
	c1->Update();

	c1->SaveAs(Form("%sMeanVsMult_pT%s%s_B%d_%03d.%s",
			picPrefix.Data(),
			Form("%03.0f_%03.0f",inPtLo,inPtUp),
			(iPlotType&(1<<2))?"_rp":"",
			iB,iPlotFlag,picSuffix.Data()));
	if(!gROOT->IsBatch()){
	if(getchar()=='q')return 1;
	}
	TH2F* hFrameSigmaMult = new TH2F("hSigmaMult",";raw #tracks;#sigma of LHS Gaus",500,0.,3200.,10,0.,18.);
	hFrameSigmaMult->SetLabelSize(hFrame->GetLabelSize("Y")*0.9,"Y");
	hFrameSigmaMult->SetLabelSize(hFrame->GetLabelSize("X")*0.9,"X");
	hFrameSigmaMult->SetTitleSize(hFrame->GetTitleSize("Y")*0.9,"Y");
	hFrameSigmaMult->SetTitleSize(hFrame->GetTitleSize("X")*0.9,"X");
	hFrameSigmaMult->SetTitleOffset(1.1,"Y");
	hFrameSigmaMult->SetTitleOffset(1.1,"X");
	hFrameSigmaMult->SetStats(kFALSE);
	hFrameSigmaMult->DrawCopy();
	
	TLegend *legSM = new TLegend(0.2,0.7,0.3,0.98);
	legSM->SetHeader(Form("Pb+Pb R = 0.4 (B%d)",iB));
	legSM->SetFillColor(0);
	legSM->SetTextFont(gStyle->GetTextFont());
	legSM->SetTextSize(gStyle->GetTextSize()*0.6);
	legSM->SetBorderSize(0);
	
	for(iRP = 0;iRP<((iPlotType&(1<<2))?nRP:1);iRP++){
	  for(iDelta = 0;iDelta <kDeltaTypes;iDelta++){
	    if(!(iPlotFlag&(1<<iDelta)))continue;
	    if(!grSigmaDeltaPtMult[iDelta][iRP])continue;
	    grSigmaDeltaPtMult[iDelta][iRP]->Draw("psame");
	    legSM->AddEntry(grSigmaDeltaPtMult[iDelta][iRP],sDelta[iDelta][iRP].Data(),"P");
	  }
	}
	legSM->Draw();
	c1->Update();
	c1->SaveAs(Form("%sSigmaVsMult_pT%s%s_B%d_%03d.%s",
			picPrefix.Data(),
			Form("%03.0f_%03.0f",inPtLo,inPtUp),
			(iPlotType&(1<<2))?"_rp":"",
			iB,iPlotFlag,picSuffix.Data()));
	if(!gROOT->IsBatch()){
	if(getchar()=='q')return 1;
	}
	
    }// trends vs. mult

    // plot the trend with RP and vs. mult
    CloseFiles();
    return 0;


}


    void ScaleH1(TH1* h,char* cX,char* cY,Float_t fScale,Bool_t bWidth){
  if(!h)return;
  if(!fScale)return;

  h->Scale(fScale);
  if(bWidth){
    for(int ib = 1;ib <= h->GetNbinsX();ib++){
      Float_t val = h->GetBinContent(ib);
      Float_t err = h->GetBinError(ib);
      Float_t w = h->GetBinWidth(ib);
      h->SetBinContent(ib,val/w);
      h->SetBinError(ib,err/w);
      // Printf("width %f",w);
    }
  }
  h->SetXTitle(cX);
  h->SetYTitle(cY);

}


void ScaleH2(TH2* h,char* cX,char* cY,char* cZ,Float_t fScale,Bool_t bWidth){
  if(!h)return;
  if(!fScale)return;
  
  h->Scale(fScale);
  if(bWidth){
    for(int ibx = 1;ibx <= h->GetNbinsX();ibx++){
      for(int iby = 1;iby <= h->GetNbinsY();iby++){
	Float_t val = h->GetBinContent(ibx,iby);
	Float_t err = h->GetBinError(ibx,iby);
	Float_t wx = h->GetXaxis()->GetBinWidth(ibx);
	Float_t wy = h->GetYaxis()->GetBinWidth(iby);
	h->SetBinContent(ibx,iby,val/(wx*wy));
	h->SetBinError(ibx,iby,err/(wx*wy));
	// Printf("width %f",w);
      }
    }
  }
  h->SetXTitle(cX);
  h->SetYTitle(cY);
  h->SetZTitle(cZ);

}

void SetStyleH1(TH1 *h,Int_t color,Int_t fillStyle,Int_t markerStyle){
  if(color){
    h->SetLineColor(color);
    h->SetLineWidth(3);
    h->SetMarkerColor(color);
    h->SetFillColor(color);
  }
  if(fillStyle)h->SetFillStyle(fillStyle);
  if(markerStyle)h->SetMarkerStyle(markerStyle);
}

void HistsFromSingleJets(char* cMask,TList *list,Int_t nJets,TH1F** h,Int_t iRebin,int iFirst){

  for(int ij = iFirst;ij < nJets;ij++){
    h[ij] = (TH1F*)list->FindObject(Form(cMask,ij));

    if(h[ij]){
      if(iRebin>1){
        h[ij] = (TH1F*)h[ij]->Clone();
	h[ij]->SetDirectory(gROOT);
        h[ij]->Rebin(iRebin);
      }
      if(ij==iFirst){
        h[nJets] = (TH1F*)h[ij]->Clone(Form(cMask,nJets));
      }
      else{
        h[nJets]->Add(h[ij]);
      }
    }
    else{
      Printf("%s not found",Form(cMask,ij));
    }
  }

}


void HistsFromSingleJets(char* cMask,TList *list,Int_t nJets,TH2F** h,Int_t iRebin,int iFirst){

  for(int ij = iFirst;ij < nJets;ij++){
    h[ij] = (TH2F*)list->FindObject(Form(cMask,ij));
    if(h[ij]){
      if(iRebin>1){
        h[ij] = (TH2F*)h[ij]->Clone();
        h[ij]->SetDirectory(gROOT);
        h[ij]->Rebin(iRebin);
      }
      if(ij==iFirst){
        h[nJets] = (TH2F*)h[ij]->Clone(Form(cMask,nJets));
      }
      else{
        h[nJets]->Add(h[ij]);
      }
    }
    else{
      Printf("%s not found",Form(cMask,ij));
    }
  }

}

void ScaleNevent(TH1* h,TFile *fIn,Int_t iCond,Int_t iMC,Int_t iCen){
  // fetch the trigger histos                                                                                                                                                      

  TDirectory *dir = (TDirectory*)fIn->Get("PWG4_services");
  TList *list1 = (TList*)dir->Get("pwg4serv");

  TH2F* hTriggerCount =(TH2F*)list1->FindObject("fh2TriggerCount");
  Float_t nevents = hTriggerCount->GetBinContent(hTriggerCount->FindBin(0+iCen,iCond)); // we have to scale with number of triggers?                                               

  Float_t xsection = 1;
  if(iMC==7000){
    // 7 TeV sigma ND = 48.85 sigma NEL = 71.39 (Pythia 8... confirm with 6                                                                                                        
    xsection = 71.39;
    Printf("MC Scaling setting nevents=%f to 1, xsection = %E",nevents,xsection);
    nevents = 1;
  }

  Printf("Scaling %s with number of event %f",h->GetName(),nevents);
  if(nevents){
    Float_t scalef = 1./nevents/xsection;
    Printf("scale factor %E",scalef);
    h->Scale(scalef);
  }

}

Double_t Gamma0(Double_t *x,Double_t *par){

  // Fixed multiplicity M


  Double_t p = par[0];
  Double_t b = par[1];
  Double_t A = par[2];
  Double_t xval = x[0] + p/b;
  // xval += par[3];
  if(xval<0)return 0;

  // take log to avoid zeros and NANs
  Double_t f1 = TMath::Log(A*b);
  Double_t f11 = 1;// ROOT::Math::lgamma(p); 
  Double_t f2 = TMath::Log(b * xval)*(p-1);
  Double_t f3 = -1.*b*xval;
  Double_t f = f1-f11+f2+f3;
  f = TMath::Exp(f);
  return f;
}

Double_t Gamma(Double_t *x,Double_t *par){

  // Fixed multiplicity M
  Double_t M = par[0];
  // Paramter p = <p_T>^2/sigma_p_T^2 adn b = <p_T>/sigma_p_T^2
  Double_t b = par[1]/par[2]/par[2];
  Double_t p = par[1] * b;

  Double_t f = M * b / 1; // ROOT::Math::tgamma(M*p) * TMath::Power(M * b * x[0],M*(p-1)) * TMath::Exp(M*b*x[0]); 

  return f;
}


TF1* FitLHSgaus(TH1D *hDeltaPt, double minPt, double maxPt, int minIterations, int maxIterations, double maxDeltaMean, int mode, double minBoundSigma, double maxBoundSigma)
{

  Double_t minPtBound = minPt;
  Double_t maxPtBound = maxPt;

  Int_t nIter = 0;
  Double_t deltaMean = 999.;
  Double_t oldMean   = 0.;

  TF1 *f1 = new TF1("f1","gaus",minPtBound,maxPtBound);
  Double_t sigma = 0.;
  Double_t sigma_error = 0.;
  Double_t mean = 0.;
  Double_t mean_error = 0.;

  while(nIter<minIterations || (nIter<maxIterations && deltaMean>maxDeltaMean)){

     if(nIter>0){ // for initial fit use minPt and maxPt
        if(mode==0){
           maxPtBound = mean+5.;
        }
        if(mode==1){
           minPtBound = mean-(minBoundSigma*sigma);
           maxPtBound = mean+(maxBoundSigma*sigma);
        }
     }

     f1->SetRange(minPtBound, maxPtBound);
     hDeltaPt->Fit(f1,"R0");
     Printf("fit range: %2.2f - %2.2f", minPtBound, maxPtBound);
     mean = f1->GetParameter(1);
     sigma = f1->GetParameter(2);
     mean_error = f1->GetParError(1);
     sigma_error = f1->GetParError(2);

     deltaMean = TMath::Abs(mean-oldMean);
     oldMean = mean;
     Printf("FIT %d:  #mu = %2.2f +- %2.2f (diff %2.2f), #sigma = %2.2f +- %2.2f, #chi_{2}/NDF = %2.2f/%d", nIter, mean, mean_error, deltaMean, sigma, sigma_error, f1->GetChisquare(), f1->GetNDF());
     nIter++;
  }

  f1->SetRange(f1->GetX(1E-5,-70.,0.)-0.75,f1->GetX(1E-5,0.,70.)+0.75);

  return f1;
}



TF1* FitLHSgaus(TH1D *hDeltaPt, double &sigma, double &sigma_error, double &mean)
{

  TF1 *f1 = new TF1("LHSgaus","gaus");
  f1->SetParameters(1,0,11);
  hDeltaPt->Fit(f1,"R0","",-50.,0.);
  //  f1 = hDeltaPt->GetFunction("gaus");
  sigma = f1->GetParameter(2);
  mean = f1->GetParameter(1);
  hDeltaPt->Fit(f1,"R0","",mean-3.*sigma,mean+0.3*sigma);
  mean = f1->GetParameter(1);
  sigma = f1->GetParameter(2);
  sigma_error = f1->GetParError(2);

  return f1;
}




TH2F* GetTH2PlotB(const char *cPath,Int_t embType, Int_t classType, Int_t cl, Int_t rp){
   
   // emb type 0: single tracks, 1: emb jets
   // class type 0: centrality, 1: multiplicity (nb. of input tracks)
   // cl: centrality or multplicity class, -1 all
   // rp -1: all, 0: in plane, 1: in between, 2: out of plane
   
   TString sClType;
   if(classType) sClType = "mult";
   else          sClType = "cent";
   TString sEmbType;
   if(embType==0) sEmbType = "singles";
   if(embType==1) sEmbType = "jets";
   
   TString path(cPath);
   TString file = Form("%s/jetresponse_%s_%s.root",path.Data(),sClType.Data(),sEmbType.Data());
   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(file.Data());
   if(!f)f = TFile::Open(file);

   if(!f){
      Printf("GetTH2PlotB:%d Could not open file %s",__LINE__,file.Data());
      return 0x0;
   }
   
   TString dfName = Form("JetsDeltaPt/rp %d/%s %d",rp,sClType.Data(),cl);
   TDirectory* df = (TDirectory*)f->Get(dfName.Data());
   if(!df){
      Printf("GetTH2PlotB:%d Could not open directory %s",__LINE__,dfName.Data());
      return 0x0;
   }
   
   TString h2Name = "delta pT vs Probe pT";
   TH2F* h2 = (TH2F*)(df->Get(h2Name.Data()))->Clone(Form("h2B_%d",cl));
   if(!h2){
      Printf("GetTH2PlotB:%d Could not get histogram %s",__LINE__,h2Name.Data());
      return 0x0;
   }
   
   return h2;
}

TObject* GetObjectFromFile(const char *cName,const char *cFile,const char* cDir,const char *cRep){
   TDirectory *opwd = gDirectory;
   TFile *fIn = (TFile*)gROOT->GetListOfFiles()->FindObject(cFile);
   if(!fIn)fIn = TFile::Open(cFile);
   opwd->cd();

   if(!fIn){
     Printf("File %s not found",cFile);
     return 0;
   }
   TDirectory *dir = (TDirectory*)fIn->Get(cDir);
   if(!dir){
     Printf("dir %s not found",cDir);
     return 0;
   }
   TString sList(cDir);
   sList.ReplaceAll("PWG4_",cRep);
   TList *list = (TList*)dir->Get(sList.Data());
   if(!list){
     Printf("list %s not found",sList.Data());
     return 0;
   }
   
   TObject *obj = list->FindObject(cName);
   if(!obj){
     Printf("object %s not found",cName);
     return 0;
   }
   return obj;
 }

Double_t GetPoissonFluctuation(TH1 *h1,Double_t areaIn,Double_t areaJet){
  if(!h1)return 0;
  Printf(">>>> %s ",h1->GetName());
  Double_t meanPt = h1->GetMean();
  Double_t rmsPt = h1->GetRMS();
  Double_t mult = h1->Integral("width");

  Double_t multJet = mult/areaIn*areaJet;
  Double_t sigma = TMath::Sqrt(multJet) * TMath::Sqrt(meanPt*meanPt+rmsPt*rmsPt); 
  Printf("MeanPt %6.3f RMS %6.3f Tracks: %d",meanPt,rmsPt,(Int_t)mult);
  Printf("Tracks/jet: %d SigmaJet: %6.3f",(Int_t)multJet,sigma);
  return sigma;
}

void SetHistoAttributes(TH1* h1,Int_t iMarker,Int_t iColor){
  if(!h1)return;
  h1->SetMarkerColor(iColor);
  h1->SetLineColor(iColor);
  h1->SetMarkerStyle(iMarker);
}

void SetGraphAttributes(TGraph* gr,Int_t iMarker,Int_t iColor){
  if(!gr)return;
  gr->SetMarkerColor(iColor);
  gr->SetLineColor(iColor);
  gr->SetMarkerStyle(iMarker);
}


void CloseFiles(){
  TSeqCollection *coll = gROOT->GetListOfFiles();
  for(int i = 0;i<coll->GetEntries();i++){
    TFile *f = (TFile*)coll->At(i);
    Printf("Closing %d %s",i,f->GetName());
    f->Close();
  }
}
