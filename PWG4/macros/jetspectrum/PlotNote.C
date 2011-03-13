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
Int_t PlotJetBFluctuations();
Int_t PlotJetQA();
// helpers


enum {kMaxJets = 4};
void HistsFromSingleJets(char* cMask,TList *list,Int_t nJets,TH1F** h,Int_t iRebin = 1,int iFirst = 0);
void HistsFromSingleJets(char* cMask,TList *list,Int_t nJets,TH2F** h,Int_t iRebin = 1,int iFirst = 0);
void ScaleNevent(TH1 *h,TFile *fIn,Int_t iCond = 5,Int_t iMC = 0,Int_t iCen = 0);
void ScaleH1(TH1* h,char* cX = "",char* cY = "",Float_t fScale = 1,Bool_t bWidth = true);
void SetStyleH1(TH1 *h,Int_t color = 0,Int_t fillStyle = 0,Int_t markerStyle = 0);
void ScaleH2(TH2* h,char* cX = "",char* cY = "",char* cZ = "",Float_t fScale = 1,Bool_t bWidth = true);
TH1F* GetPythia8Spectrum(const char *cHist,const char *cPythia,Float_t deta,Int_t iRebin = 1,Int_t iMask = 0);
const char *cPrintMask = "110116_%s.eps";
void set_plot_style();
TCanvas *cTmp = 0;

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
 
  PlotJetBFluctuations();
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
  TString picName = "jetspectrumPbPb";
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
  const Float_t fCentLo[nCen] = {0,10,30,50};
  const Float_t fCentUp[nCen] = {10,30,50,80};
  TH2F *hFrame = new TH2F("hFrame",";#delta p_{T} (GeV/c);Probability/GeV",200,-70,70,100,1E-5,5); 
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

  TFile::SetCacheFileDir("/tmp/");

  // Martha, single particle jets
  TFile *fM = TFile::Open("http://qgp.uni-muenster.de/~stevero/tmp/jets/PWG4_JetTasksOutput_AOD_EmbeddingSingleTrack.root","CACHEREAD");
  TH1D *hDeltaPtM[nCen] = {0};

  // select 
  Float_t fMinPtM = 40;
  Float_t fMaxPtM = 80;
  int iB = 2;

  /*
    0: 0-10%
    1: 10-30%
    2: 30-50%
    3: 50-80%
  */




  for(int ic = 0;ic < nCen;ic++){
    tmpName = Form("PWG4_BkgFluctCent%dB%d",ic,iB);
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
    hDeltaPtM[ic]->SetMarkerStyle(kFullSquare);
    hDeltaPtM[ic]->SetMarkerColor(kGreen+2);
    hDeltaPtM[ic]->SetLineColor( hDeltaPtM[ic]->GetMarkerColor());
  }    
  
  // fetch the BiAs
  
  TH1D *hBiaL[nCen];

  TFile *fL = TFile::Open("~/alice/jets/macros/corrections/tmp/pwg4plots.root");

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

  for(int ic = 0;ic < nCen;ic++){
    if(iB==1)tmpName = "BiA sa RC centrality";
    else if (iB==2)tmpName = "BiA va RC centrality";
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
  txt2->DrawLatex(800,150,"ALICE Performance");
  txt2->DrawLatex(800,140,"01/03/2011");
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
  txt->DrawLatex(20,180,"LHC 2010 Pb+Pb Run #sqrt{s_{NN}} = 2.76 TeV");
  txt2->DrawLatex(50,150,"ALICE Performance");
  txt2->DrawLatex(50,140,"01/03/2011");
  c1->Update();
  c1->SaveAs(Form("rhovscent_B%d.%s",iB,printType.Data()));
  if(getchar()=='q')return 1;

  // fetch the data from bastian...
  Float_t fMinPtB = fMinPtM;
  Float_t fMaxPtB = fMaxPtM;


  // the embbedded jets, carefull, take only above 60 GeV
  // c = all jets, d = leading jets, e = single tracks
  TFile *fB1 = TFile::Open("~/alice/jets/macros/corrections/tmp/PWG4_JetTasksOutput_e.root");
  TH1D *hDeltaPtB1[nCen] = {0};
  for(int ic = 0;ic < nCen;ic++){
    tmpName = Form("PWG4_JetResponse_%s_%s%02d_B%d_Filter00256_Cut%05d_Skip00","clusters", "ANTIKT", (Int_t)(0.4*10), iB, (Int_t)(0.15*1000));
    TDirectoryFile *df = dynamic_cast<TDirectoryFile*> (fB1->Get(tmpName.Data()));
    if(!df)Printf("%d %s not found",__LINE__,tmpName.Data());
    tmpName.ReplaceAll("PWG4_JetResponse","jetresponse");
    TList *list        = dynamic_cast<TList*> (df->Get(tmpName.Data()));
    if(!list)Printf("%d %s not found",__LINE__,tmpName.Data());
    tmpName = Form("pt_smearing%d",ic+1);
    TH2F *hTmp = (TH2F*)list->FindObject(tmpName.Data());
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


  c1->SetLogy();
  c1->SetMargin(0.15,0.05,0.2,0.05);

  TF1 *gaus = new TF1("gaus","gaus",-60,2);
  TF1 *gaus2 = new TF1("gaus2","gaus",-60,2);
  for(int ic = 0;ic < nCen;ic++){
    TLegend *leg1 = new TLegend(0.2,0.78,0.3,0.93);
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

    hBiaRC[ic]->DrawCopy("psame");
    gaus->SetLineColor(hBiaRC[ic]->GetLineColor());
    hBiaRC[ic]->Fit(gaus,"R0");
    gaus->SetRange(-40.,40.);
    gaus->Draw("same");
    leg1->AddEntry(hBiaRC[ic],Form("BiA random cones (excl. leading jets)"),"P");
    leg1->AddEntry(gaus,Form("LHS Gaus fit: #mu = %2.2f, #sigma = %2.2f",gaus->GetParameter(1),gaus->GetParameter(2)),"L");



    //    hDeltaPtB1[ic]->Draw("psame");

    
    hDeltaPtM[ic]->DrawCopy("psame");
    gaus2->SetLineColor(hDeltaPtM[ic]->GetLineColor());
    hDeltaPtM[ic]->Fit(gaus2,"R0");
    gaus2->SetRange(-40.,40.);
    //    gaus2->Draw("same");
    leg1->AddEntry(hDeltaPtM[ic],Form("anti-k_{T} single tracks %2.0f-%2.0f GeV",fMinPtM,fMaxPtM),"P");
    leg1->Draw();
    txt2->SetNDC();
    txt2->DrawLatex(0.5,0.97,"ALICE Performance 01/03/2011");
    c1->Update();
    c1->SaveAs(Form("deltaPt_B%d_cen%02d.%s",iB,ic,printType.Data()));
    if(getchar()=='q')return 1;


  }



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

