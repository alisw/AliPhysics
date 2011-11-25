#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TCanvas.h>
#include <TGrid.h>
#include <TFile.h>
#include <TList.h>
#include <TGridResult.h>
#include <TPaveStats.h>
#include <TGraph.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLatex.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h>
#endif

void PlotGeneral(TFile* fildat, TCanvas**& clist, Int_t& cnum);
void PlotITSsa(TFile* fildat, TCanvas**& clist, Int_t& cnum);
void PlotSDD(TFile* fildat, TCanvas**& clist, Int_t& cnum);
void GetGainModuleLevelSSD(TFile* fildat, TCanvas**& clist, Int_t& cnum);
void VertexQAMacro(TFile *fildat, TCanvas **&clist, Int_t &cnum);
void PlotSPD(TFile* fildat, TFile* filMC, TCanvas**& clist, Int_t& cnum);
Bool_t PlotITSTPCMatchingEff(TFile *f, TCanvas**& clist,Int_t& cnum);
void SetDrawAtt(Int_t markerstyle,Int_t markercolor,Int_t markersize,Int_t linecolor,Int_t linewidth,TH1 *h1);
Double_t LangausFun(Double_t *x, Double_t *par);
void SaveC(TFile &fout, TCanvas**& clist, Int_t cnum);
TString GetRunNumber();

//  the run number is available to all the functions. Its value is set by AliITSQAchecks
  Int_t gRunNumber = 0;
  Int_t gRunNumberMC = 0;
  TString pdfFileNames="";


//_______________________________________________________________________
void AliITSQAchecks(TString option="grid",
			  Int_t nRun=167713,
			  TString period="LHC11h",
			  TString qaTrain="QA90",
		    TString filenamedata="QAresults.root", TString filenameMC="alien:///alice/data/2011/LHC11h/000167706/ESDs/pass1_HLT/QA90/QAresults.root",Int_t nRunMC=0){
  // THIS MACRO SHOULD BE COMPILED. IT DOES NOT WORK WITH THE INTERPRETER
  // option:  "local" if filenamedata is the name of a local file
  //          "grid" if on alien
  // nRun:    run number
  // period:  LHC period (e.g. LHC11h)
  // qaTrain: QA train specifier  
  //          Empty string if QAresults.root is in the ESDs/pass1_HLT directory 
  // filenamedata: QAresults.root is by default the file name with the results
  // filenameMC: file name for MC comparison. If the names begins with alien:
  //             the file is accessed through alien, otherwise is taken as local
  // nRunMC:  run number for comparison. If filenamMC begings with "alien:" 
  //          the run number is taken from the path. Otherwise, in case of a 
  //          local filenameMC, the run number must be specified here
  // Select here what you want to display
  // the complete selection string is
  // "general ITSSA SPD SDD SSD vertex ITSTPC"
  // Contact:  Stefania Beole': beole@to.infn.it  

/* $Id$ */

  gRunNumber = nRun;
  TString aux(filenameMC);
  if(aux.BeginsWith("alien:")){
    aux=aux.Remove(0,35);
    aux=aux.Remove(6,aux.Length());  
    gRunNumberMC = atoi(aux.Data());
  }
  else {
    gRunNumber = nRunMC;
  }

  TString selection("general ITSSA SPD SDD SSD vertex ITSTPC"); 
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1111);
  TFile *fildat;
  TString path;
  Int_t year=2011;
  if(period.Contains("LHC10")) year=2010;
  else if(period.Contains("LHC09")) year=2009;

  if(option.Contains("local")){
    fildat=new TFile(filenamedata.Data());
    printf("Opened file %s\n",fildat->GetName());
  }else{
    TGrid::Connect("alien:");
    if(qaTrain.Contains("QA")){
      path=Form("/alice/data/%d/%s/%09d/ESDs/pass1_HLT/%s/",year,period.Data(),nRun,qaTrain.Data());
    } else {
      path=Form("/alice/data/%d/%s/%09d/ESDs/pass1_HLT/",year,period.Data(),nRun);
    }
    filenamedata = "alien://"+path+"QAresults.root";
    fildat=TFile::Open(filenamedata.Data());
  }
  if(option.Contains("local") && filenameMC.Contains("alien"))TGrid::Connect("alien:");
  TFile* filMC=TFile::Open(filenameMC.Data());
  TCanvas** clist;
  Int_t cnum;
  char rn[10];
  sprintf(rn,"%d",gRunNumber);
  TString strRN(rn);
  TString founame="Outfil"+strRN+".root";
  TFile fout(founame,"recreate");
  if(selection.Contains("general")){
    PlotGeneral(fildat,clist,cnum); 
    printf("GENERAL - cnum = %d\n",cnum);
    SaveC(fout,clist,cnum);

  }
  if(selection.Contains("ITSSA")){
    PlotITSsa(fildat,clist,cnum); 
    printf("ITSSA - cnum = %d\n",cnum);
    SaveC(fout,clist,cnum);
  }
  if(selection.Contains("SDD")){
    PlotSDD(fildat,clist,cnum); 
    printf("SDD - cnum = %d\n",cnum);
    SaveC(fout,clist,cnum);
  }
  if(selection.Contains("SSD")){
    GetGainModuleLevelSSD(fildat,clist,cnum);
    printf("SSD - cnum = %d\n",cnum);
    SaveC(fout,clist,cnum);
  }
  if(selection.Contains("vertex")){
    VertexQAMacro(fildat,clist,cnum);
    printf("VERTEX - cnum = %d\n",cnum);
    SaveC(fout,clist,cnum);
  }
  if(selection.Contains("SPD")){
    PlotSPD(fildat,filMC,clist,cnum);
    printf("SPD - cnum = %d\n",cnum);
    SaveC(fout,clist,cnum);
  }
  if(selection.Contains("ITSTPC")){
    PlotITSTPCMatchingEff(fildat,clist,cnum);
    printf("ITSTPC - cnum = %d\n",cnum);
    SaveC(fout,clist,cnum);
  }

  fout.Close();

  // merge the pdf files
  TString command("gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=merged");
  command=command+strRN+".pdf "+pdfFileNames;
  gSystem->Exec(command.Data());
  printf(" Merging the pdf file:  %s \n",command.Data());
  
}

//_______________________________________________________________________
void PlotGeneral(TFile* fildat, TCanvas**& clist, Int_t& cnum){
  TDirectoryFile* df=(TDirectoryFile*)fildat->Get("SDD_Performance");
  if(!df){
    printf("SDD_Performance MISSING -> Exit\n");
    return;
  }
  TList* l=(TList*)df->Get("coutputRP");
  if(!df){
    printf("coutputRP TList MISSING -> Exit\n");
    return;
  }
  cnum=1; // number of canvases 
  clist= new TCanvas* [1];//array of pointers to TCanvases
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1111);
  TH1F* hcllay=(TH1F*)l->FindObject("hCluInLay");
  TH1F* hev=(TH1F*)l->FindObject("hNEvents");
  Int_t nTotEvents=hev->GetBinContent(2);
  Int_t nTrigEvents=hev->GetBinContent(3);
  Int_t nEvents=nTotEvents;
  printf("---- Statistics ----\n");
  printf("Number of Events = %d\n",nTotEvents);
  if(nTrigEvents>0){ 
    printf("Number of Triggered Events = %d\n",nTrigEvents);
    nEvents=nTrigEvents;
  }else{
    printf("No request on the trigger done when running the task\n");
  }
  if(hcllay){
    Double_t norm=hcllay->GetBinContent(1);
    if(norm>0.){
      hcllay->Scale(1./norm);
      hcllay->SetTitle("");
      hcllay->GetXaxis()->SetRange(2,7);
      hcllay->SetMinimum(0.);
      hcllay->SetMaximum(1.1);
      hcllay->SetMarkerStyle(23);
      TString ctitle=GetRunNumber()+"General checks: PointPerLayer";
      TCanvas* ceffL=new TCanvas("ceffL",ctitle,1000,800);
      clist[0]=ceffL;
      // ceffL->Divide(1,2);
      // ceffL->cd(1);
      ceffL->SetGridy();
      hcllay->Draw(); 
      TLatex* tg=new TLatex(0.15,0.2,"Fraction of tracks with point in ITS layer");
      tg->SetTextSize(0.04);
      tg->SetNDC();
      tg->SetTextColor(1);
      tg->Draw();
      TString testo="Run "+GetRunNumber();
      TLatex* tg2 = new TLatex(0.15,0.85,testo.Data());
      tg2->SetTextSize(0.04);
      tg2->SetNDC();
      tg2->SetTextColor(2);
      tg2->Draw();
      hcllay->GetXaxis()->SetTitle("Layer");
      hcllay->GetYaxis()->SetTitle("Fraction of tracks with point in layer");
      ceffL->Update();
      ceffL->SaveAs("track_points_per_layer.pdf");
      pdfFileNames+=" track_points_per_layer.pdf";
    }
  }
}


//_______________________________________________________________________
//////////////////////////////////////////////////////////////////////
/// ITSsa ////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
 void PlotITSsa(TFile* fildat, TCanvas**& clist, Int_t& cnum){
    TDirectoryFile* df=(TDirectoryFile*)fildat->Get("TracksITSsa");
    if(!df) df=(TDirectoryFile*)fildat->Get("ITSsaTracks");
    if(!df){
      printf("ITSsa_Performance MISSING -> Exit\n");
      return;
    }
 
    TList* l=(TList*)df->Get("clistITSsaTracks");
    if(!df){
      printf("clistITSsaTracks TList MISSING -> Exit\n");
      return;
    }
    cnum=2; // number of canvases 
    clist= new TCanvas* [2];//array of pointers to TCanvases
     gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1111);


  TH1F* hPtTPCITS=(TH1F*)l->FindObject("hPtTPCITS");
  TH1F* hPtITSsa=(TH1F*)l->FindObject("hPtITSsa");
  TH1F* hPtITSpureSA=(TH1F*)l->FindObject("hPtITSpureSA");

  TH2F* hEtaPhiTPCITS=(TH2F*)l->FindObject("hEtaPhiTPCITS");
  TH2F* hEtaPhiITSsa=(TH2F*)l->FindObject("hEtaPhiITSsa");
  TH2F* hEtaPhiITSpureSA=(TH2F*)l->FindObject("hEtaPhiITSpureSA");
  TH1F* hChi2TPCITS=(TH1F*)l->FindObject("hChi2TPCITS");
  TH1F* hChi2ITSsa=(TH1F*)l->FindObject("hChi2ITSsa");

  TH1F* hRatio=(TH1F*)hPtTPCITS->Clone("hRatio");
  TH1F* hRatio1=(TH1F*)hPtTPCITS->Clone("hRatio1");
   hRatio->Add(hPtITSsa);
  hRatio->Divide(hPtITSpureSA);
  hRatio->SetStats(0);
  hRatio1->Divide(hPtITSsa);
  hRatio1->SetStats(0);

  TString ctitle=GetRunNumber()+"ITS standalone: performance vs Pt";
  TCanvas* cITSsa1=new TCanvas("cITSsa1",ctitle,1200,1200);
  clist[0]=cITSsa1;
  cITSsa1->Divide(1,3);
  cITSsa1->cd(1);
  // hPtITSpureSA->Draw();
  // hPtITSpureSA->GetXaxis()->SetTitle("Pt (GeV/c)");
  // gPad->Update();
  // TPaveStats *st1=(TPaveStats*)hPtITSpureSA->GetListOfFunctions()->FindObject("stats");
  // st1->SetY1NDC(0.71);
  // st1->SetY2NDC(0.9);
  hPtTPCITS->SetLineColor(2);
  hPtTPCITS->GetXaxis()->SetTitle("Pt (GeV/c)");
  //  hPtTPCITS->Draw("sames");
  hPtTPCITS->Draw();
  gPad->Update();
  TPaveStats *st2=(TPaveStats*)hPtTPCITS->GetListOfFunctions()->FindObject("stats");
  st2->SetY1NDC(0.71);
  st2->SetY2NDC(0.9);
  st2->SetTextColor(2);

  hPtITSsa->SetLineColor(4);
  hPtITSsa->Draw("sames");
  gPad->Update();
  TPaveStats *st3=(TPaveStats*)hPtITSsa->GetListOfFunctions()->FindObject("stats");
  st3->SetY1NDC(0.51);
  st3->SetY2NDC(0.7);
  st3->SetTextColor(4);
  TLegend* leg=new TLegend(0.5,0.5,0.69,0.79);
  leg->SetFillColor(0);
  TLegendEntry* ent=leg->AddEntry(hPtTPCITS,"TPC+ITS","L");
  ent->SetTextColor(hPtTPCITS->GetLineColor());
  ent=leg->AddEntry(hPtITSsa,"ITSsa","L");
  ent->SetTextColor(hPtITSsa->GetLineColor());
   // to be used only with pp data (ITS pure SA)  
 // ent=leg->AddEntry(hPtITSpureSA,"ITS pureSA","L");
  //ent->SetTextColor(hPtITSpureSA->GetLineColor());
  leg->Draw();
  cITSsa1->cd(2);
   gPad->SetGridx();
   gPad->SetGridy();
   hRatio1->GetXaxis()->SetTitle("Pt (GeV/c)");
   hRatio1->GetYaxis()->SetTitle("TPCITS/ITSsa");
   //   hRatio->GetYaxis()->SetTitle("(TPCITS+ITSsa)/ITSpureSA");
   hRatio1->DrawCopy();
   TLatex* tratio=new TLatex(0.2,0.75,"TPC+ITS/ITSsa vs Pt");
   tratio->SetNDC();
   tratio->SetTextColor(1);
  tratio->Draw();
  cITSsa1->cd(3);
  hChi2ITSsa->Scale(1./hChi2ITSsa->GetEntries());
  hChi2TPCITS->Scale(1./hChi2TPCITS->GetEntries());
  hChi2TPCITS->SetLineColor(2);
  hChi2TPCITS->Draw("");
   TLatex* tchi=new TLatex(0.25,0.85,"chi2 vs Pt");
   tchi->SetNDC();
   tchi->SetTextColor(1);
  tchi->Draw();
  gPad->Update();
  TPaveStats *stc2=(TPaveStats*)hChi2TPCITS->GetListOfFunctions()->FindObject("stats");
  stc2->SetY1NDC(0.71);
  stc2->SetY2NDC(0.9);
  stc2->SetTextColor(2);
  //  c2->Update();
  hChi2ITSsa->SetLineColor(4);
  hChi2ITSsa->Draw("sames");
  gPad->Update();
  TPaveStats *stc3=(TPaveStats*)hChi2ITSsa->GetListOfFunctions()->FindObject("stats");
  stc3->SetY1NDC(0.51);
  stc3->SetY2NDC(0.7);
  stc3->SetTextColor(4);
  leg->Draw();

  cITSsa1->Update();
  cITSsa1->SaveAs("ITSsa1.pdf");
  pdfFileNames+=" ITSsa1.pdf";
  gStyle->SetPalette(1);
  hEtaPhiITSpureSA->SetStats(0);
  hEtaPhiITSpureSA->SetTitle("ITS pureSA");
  hEtaPhiITSsa->SetStats(0);
  hEtaPhiITSsa->SetTitle("ITSsa tracks");
  hEtaPhiTPCITS->SetStats(0);
  hEtaPhiTPCITS->SetTitle("TPC+ITS tracks");
  ctitle=GetRunNumber()+"Eta-phi distribution for ITSsa and TPC+ITS tracks";
  TCanvas* cITSsa2=new TCanvas("cITSsa2",ctitle,1200,800);
  clist[1]=cITSsa2;
  //  cITSsa2->Divide(3,1); for ITSpuresa
  cITSsa2->Divide(2,1);
  cITSsa2->cd(1);
  // TPad* p1=new TPad("p1","Tracking: tracks distribution in EtaPhi",0,0.5,1.,1.);
  // p1->Divide(3,1);
  // p1->cd (1);
  // hEtaPhiITSpureSA->Draw("colz");
  // hEtaPhiITSpureSA->GetXaxis()->SetTitle("Eta");
  // hEtaPhiITSpureSA->GetYaxis()->SetTitle("Phi");
  //  cITSsa2->cd(2);
  //  p1->cd(2);
  hEtaPhiITSsa->Draw("colz");
  hEtaPhiITSsa->GetXaxis()->SetTitle("Eta");
  hEtaPhiITSsa->GetYaxis()->SetTitle("Phi");
  cITSsa2->cd(2);
  //  p1->cd(3);
  hEtaPhiTPCITS->Draw("colz");
  hEtaPhiTPCITS->GetXaxis()->SetTitle("Eta");
  hEtaPhiTPCITS->GetYaxis()->SetTitle("Phi");
  //  c4->cd(4);
  cITSsa2->SaveAs("ITSsa2.pdf");  
  pdfFileNames+=" ITSsa2.pdf";
}

//-----------------------------------------------------
// ///////////  Plot SDD ////////////
//_______________________________________________________________________
void SetDrawAtt(Int_t markerstyle,Int_t markercolor,Int_t markersize,Int_t linecolor,Int_t linewidth,TH1 *h1){ 


  h1->SetMarkerStyle(markerstyle);
  h1->SetMarkerColor(markercolor);
  h1->SetMarkerSize(markersize);
  h1->SetLineColor(linecolor);
  h1->SetLineWidth(linewidth);
}

//_______________________________________________________________________
Double_t LangausFun(Double_t *x, Double_t *par) {

  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation), 
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location

  // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;


  // MP shift correction
  mpc = par[1] - mpshift * par[0]; 

  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];

  step = (xupp-xlow) / np;

  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);

    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }

  return (par[2] * step * sum * invsq2pi / par[3]);

}
//_________________________________________________________________________
void PlotSDD(TFile* fildat, TCanvas**& clist, Int_t& cnum){
  TDirectoryFile* df=(TDirectoryFile*)fildat->Get("SDD_Performance");
    if(!df){
      printf("SDD_Performance MISSING -> Exit\n");
      return;
    }
    TList* l=(TList*)df->Get("coutputRP");
    if(!df){
      printf("coutputRP TList MISSING -> Exit\n");
      return;
    }
    cnum=1; // number of canvases 
    clist= new TCanvas* [cnum];//array of pointers to TCanvases
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1111);
  TH1F* htimT=(TH1F*)l->FindObject("hDrTimTPAll");
  TH1F* htimTe=(TH1F*)l->FindObject("hDrTimTPExtra");
  TH1F* htimTne=(TH1F*)l->FindObject("hDrTimTPNoExtra");
  htimT->Rebin(4);
  htimTe->Rebin(4);
  htimTne->Rebin(4);
  htimT->SetLineWidth(2);
  htimTe->SetLineWidth(2);
  htimTne->SetLineWidth(2);  
  // TH1F* hev=(TH1F*)l->FindObject("hNEvents");
  // Int_t nTotEvents=hev->GetBinContent(2);
  // Int_t nTrigEvents=hev->GetBinContent(3);
  // Int_t nEvents=nTotEvents;
  // printf("---- Statistics ----\n");
  // printf("Number of Events = %d\n",nTotEvents);
  // if(nTrigEvents>0){ 
  //   printf("Number of Triggered Events = %d\n",nTrigEvents);
  //   nEvents=nTrigEvents;
  // }else{
  //   printf("No request on the trigger done when running the task\n");
  // }
  // if(hcllay){
  //   Double_t norm=hcllay->GetBinContent(1);
  //   if(norm>0.){
  //     hcllay->Scale(1./norm);
  //     hcllay->SetTitle("");
  //     hcllay->GetXaxis()->SetRange(2,7);
  //     hcllay->SetMinimum(0.);
  //     hcllay->SetMaximum(1.1);
  //     hcllay->SetMarkerStyle(23);
  //     TCanvas* ceffL=new TCanvas("ceffL","General: PointPerLayer",800,1000);
  //     clist[0]=ceffL;
  //     ceffL->Divide(1,2);
  //     ceffL->cd(1);
  //     ceffL->SetGridy();
  //     hcllay->Draw();
  //     hcllay->GetXaxis()->SetTitle("Layer");
  //     hcllay->GetYaxis()->SetTitle("Fraction of tracks with point in layer");
  //     ceffL->Update();
  //   }
  // }
  TH1F* hSigTim[8];
  TGraphErrors* gmpv=new TGraphErrors(0);
  TGraphErrors* gsigg=new TGraphErrors(0);
  TGraphErrors* gsigl=new TGraphErrors(0);
  gmpv->SetTitle("");
  gsigg->SetTitle("");
  gsigl->SetTitle("");
  Int_t iPoint=0;
  TF1 *lfun = new TF1("LangausFun",LangausFun,50.,300.,4);
  for(Int_t it=0; it<8; it++){
    hSigTim[it]=(TH1F*)l->FindObject(Form("hSigTimeInt%d",it));
    if(hSigTim[it]->GetEntries()>200){
      lfun->SetLineWidth(2);
      lfun->SetParameter(0,5.);
      lfun->SetParameter(1,80.);
      lfun->SetParameter(2,hSigTim[it]->GetEntries()/10.);
      lfun->SetParameter(3,10.);
      lfun->SetParLimits(3,0.,20);

      //      hSigTim[it]->Fit("LangausFun","QLR");
      hSigTim[it]->Fit("LangausFun","ON");
      hSigTim[it]->GetXaxis()->SetTitle(Form("dE/dx, time interval %d",it+1));
      hSigTim[it]->GetYaxis()->SetTitle("Events");
      Float_t mpv=lfun->GetParameter(1);
      Float_t empv=lfun->GetParError(1);
      Float_t sig=lfun->GetParameter(3);
      Float_t esig=lfun->GetParError(3);
      Float_t sigl=lfun->GetParameter(0);
      Float_t esigl=lfun->GetParError(0);
      gmpv->SetPoint(iPoint,(Float_t)it,mpv);
      gmpv->SetPointError(iPoint,0.,empv);
      gsigg->SetPoint(iPoint,(Float_t)it,sig);
      gsigg->SetPointError(iPoint,0.,esig);
      gsigl->SetPoint(iPoint,(Float_t)it,sigl);
      gsigl->SetPointError(iPoint,0.,esigl);
      ++iPoint;
      printf("Bin %d - MPV=%.3f  \t SigmaLandau=%.3f  \t SigmaGaus=%.3f\n",it,mpv,sigl,sig);
    }
  }
  TString ctitle=GetRunNumber()+"SDD: DriftTime - dE/dx";
  TCanvas* ctim=new TCanvas("ctim",ctitle,800,1000);
  clist[0]=ctim;
  ctim->Divide(1,2);
  ctim->cd(1);
  // htimT->Draw();
  // htimTe->SetLineColor(2);
  // htimTe->Draw("same");
  htimTne->SetLineColor(4);
  htimTne->Draw("");
  htimTne->GetXaxis()->SetTitle("Drift Time (ns)");
  htimTne->GetYaxis()->SetTitle("TrackPoints");
  htimTne->GetYaxis()->SetTitleOffset(1.2);
  // TLatex* ta=new TLatex(0.5,0.85,"All Clusters");
  // ta->SetNDC();
  // ta->SetTextColor(1);
  // ta->Draw();
  // TLatex* te=new TLatex(0.5,0.8,"Extra Clusters");
  // te->SetNDC();
  // te->SetTextColor(2);
  // te->Draw();
  //  TLatex* tn=new TLatex(0.3,0.3,"Non-Extra Clusters");
  TLatex* tn=new TLatex(0.3,0.3,"Clusters on SDD modules");
  tn->SetNDC();
  tn->SetTextColor(4);
  tn->Draw();
  TLine* tlin3=new TLine(450.,0.,450.,htimTne->GetMaximum());
    tlin3->SetLineColor(2);
    tlin3->SetLineWidth(2);
    tlin3->SetLineStyle(7);
    tlin3->Draw("same");
  TLine* tlin4=new TLine(620.,0.,620.,htimTne->GetMaximum());
    tlin4->SetLineColor(2);
    tlin4->SetLineWidth(2);
    tlin4->SetLineStyle(7);
    tlin4->Draw("same");
  TLatex* tlimit1=new TLatex(0.2,0.5,"Range for t0");
  tlimit1->SetNDC();
  tlimit1->SetTextColor(2);
  tlimit1->Draw();
  TLine* tlin5=new TLine(6200.,0.,6200.,htimTne->GetMaximum());
    tlin5->SetLineColor(2);
    tlin5->SetLineStyle(7);
    tlin5->SetLineWidth(2);
    tlin5->Draw("same");
  TLine* tlin6=new TLine(5150.,0.,5150.,htimTne->GetMaximum());
    tlin6->SetLineColor(2);
    tlin6->SetLineWidth(2);
    tlin6->SetLineStyle(7);
    tlin6->Draw("same");
  TLatex* tlimit2=new TLatex(0.6,0.5,"Range for falling edge");
  tlimit2->SetNDC();
  tlimit2->SetTextColor(2);
  tlimit2->Draw();

  //  ctim->Update();
  //  TCanvas* cpars=new TCanvas("cpars","Params",800,600);
  ctim->cd(2);
  gPad->SetLeftMargin(0.14);
  gPad->SetFrameLineWidth(2);
  gPad->SetTickx();
  gPad->SetTicky();
  gmpv->SetMarkerStyle(20);
  gmpv->SetMinimum(75);
  gmpv->SetMaximum(90);
  gmpv->GetXaxis()->SetLimits(-0.2,6.8);
  gmpv->Draw("AP");
  gmpv->GetXaxis()->SetTitle("Drift Time interval number");
  gmpv->GetYaxis()->SetTitle("Landau MPV (keV)");
  gmpv->GetXaxis()->SetTitleSize(0.05);
  gmpv->GetYaxis()->SetTitleSize(0.05);
  gmpv->GetYaxis()->SetTitleOffset(1.2);
  TLatex* tex=new TLatex(0.2,0.75,"dE/dx MPV vs Drift time interval");
  tex->SetNDC();
  tex->SetTextColor(1);
  tex->Draw();
 //  cpars->Update();
  ctim->Update();
  ctim->SaveAs("SDD.pdf");
  pdfFileNames+=" SDD.pdf";
} 

//_______________________________________________________________________
//////////////// SSD ///////////////////////
//_______________________________________________________________________
void GetGainModuleLevelSSD(TFile* fildat, TCanvas**& clist, Int_t& cnum)
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1111);
  gStyle->SetPalette(1,0);
  cnum=1;
  clist=new TCanvas*[1]; 

  TDirectoryFile* df=(TDirectoryFile*)fildat->Get("PWG1dEdxSSDQA");
  TList* listin=(TList*)df->Get("SSDdEdxQA");
  if(!listin) return;
  TH2F* fHistQ=0x0;
  fHistQ=(TH2F*)listin ->FindObject("QACharge");
  fHistQ->SetStats(111);
  fHistQ->SetTitle("SSD Charge vs module number");
if(!fHistQ) return;
  TH2F* fHistCR=(TH2F*)listin ->FindObject("QAChargeRatio");
  fHistCR->SetStats(0);
  fHistCR->SetTitle("SSD Charge Ratio vs module number");

  if(!fHistCR) return;

  TH1F* fHistMPVs=new TH1F("SSD HistMPVS","HistMPVs;MPV;N",75,70,95);
	
  TH1F* fHistCRmean=new TH1F("SSD HistCRmean","HistCRmean;CRmean;N",200,-1,1);
	
  TH1F *fMPVGraph = new TH1F("SSD MPV","MPVgraph;Module number;MPV",1698,-0.5,1697.5);
  fMPVGraph->SetMarkerColor(kRed);
  fMPVGraph->SetMarkerSize(0.5);
  fMPVGraph->SetMarkerStyle(22);
  fMPVGraph->SetStats(111111);
  
  TH1F *fCRmeanGraph = new TH1F("SSD CRmeangraph","CRmeangraph;Module number;MPV",1698,-0.5,1697.5);
  fCRmeanGraph->SetMarkerColor(kBlue);
  fCRmeanGraph->SetMarkerSize(0.5);
  fCRmeanGraph->SetMarkerStyle(23);
  fCRmeanGraph->SetStats(111111);

  Float_t mpv[1698];
  Int_t ntofit=200;
 
  //  ofstream outfiletxtbad;
  //outfiletxtbad.open("QALHC11eCPass1_bis/badModules.txt");	
    for (int i =0;i<1698;i++)
    {
      //      cout<<i<<endl;
      TString tmpQ("Q");
      tmpQ+=i;
      TString tmpCR("CR");
      tmpCR+=i;
      TH1D* fHist1DCR= fHistCR->ProjectionY(tmpCR,i+1,i+1);
      Double_t mean=fHist1DCR->GetMean();
      if(!(TMath::Abs(mean)<1.0)||fHist1DCR->GetEntries()<10)
	  continue;
      fHistCRmean->Fill(mean);
      fCRmeanGraph->SetBinContent(i+1,mean);
      fCRmeanGraph->SetBinError(i+1,fHist1DCR->GetRMS());
      fCRmeanGraph->GetYaxis()->SetTitle("CR");
      TH1D* fHist1DQ=fHistQ->ProjectionY(tmpQ,i+1,i+1);
      //check bad modules
      if(fHist1DQ->GetEntries()<ntofit)
	{
	  //outfiletxtbad<<"Low statistic \t module= "<<i<<" netries="<<fHist1DQ->GetEntries()<<endl;
	  continue;
	}
      else
	{
	  tmpQ+="fit";
	  Float_t range=fHist1DQ->GetBinCenter(fHist1DQ->GetMaximumBin());
	  TF1 *f1 = new TF1(tmpQ,LangausFun,range*0.45,range*3.0,4);
	  f1->SetParameters(7.0,range,1.0,5.5);
	  Float_t normalization=fHist1DQ->GetEntries()*fHist1DQ->GetXaxis()->GetBinWidth(2)/f1->Integral(range*0.45,range*3.0);
	  f1->SetParameters(7.0,range,normalization,5.5);
	  //f1->SetParameters(7.0,range,fHist1DQ->GetMaximum(),5.5);
	  f1->SetParNames("sigma Landau","MPV","N","sigma Gaus");
	  f1->SetParLimits(0,2.0,100.0);
	  f1->SetParLimits(3,0.0,100.0);
	  if(fHist1DQ->Fit(tmpQ,"BRQON")==0)
	    {
	      mpv[i]=f1->GetParameter(1);
	      fHistMPVs->Fill(mpv[i]);	
	      fMPVGraph->SetBinContent(i+1,f1->GetParameter(1));
	      fMPVGraph->SetBinError(i+1,f1->GetParError(1));
	      if(mpv[i]<75.0)
		{
		  //outfiletxtbad<<"MPV lower than 75 \t module="<<i<<endl;
		}	
	      if(mpv[i]>100.0)
		{
		  // outfiletxtbad<<"MPV higher than 100 \t module="<<i<<endl;		  
		}
	      if(f1->GetParError(1)>1.0)
		{
		  //outfiletxtbad<<"MPV high error on MPV  \t module="<<i<<endl;				
		}
	    }
	  else
	    {
	      mpv[i]=1;
	      //outfiletxtbad<<"BAD FIT \t module="<<i<<endl;
	      continue;
	    }	
	}	
    }	
  
    TString ctitle=GetRunNumber()+"SSD Calibration 1";
  TCanvas *c1SSD = new TCanvas("c1SSD",ctitle,1000,1000);
  clist[0]=c1SSD;
  c1SSD->Divide(2,3);
  c1SSD->cd(1);
  fHistQ->DrawCopy("colz");
  c1SSD->cd(2);
  fHistCR->DrawCopy("colz");

  //  TCanvas *c2SSD = new TCanvas("c2SSD","SSD Calibration 2",1000,1000);
  //clist[1]=c2SSD;
  //c2SSD->Divide(2,2);
  c1SSD->cd(3);
  fMPVGraph->DrawCopy();  
     TLine* tlin0=new TLine(0.,80.,1698.,80.);
    tlin0->SetLineColor(2);
    tlin0->SetLineWidth(2);
    tlin0->Draw("same");
     TLine* tlin01=new TLine(0.,90.,1698.,90.);
    tlin01->SetLineColor(2);
    tlin01->SetLineWidth(2);
    tlin01->Draw("same");
  c1SSD->cd(4); 
  fHistMPVs->DrawCopy();
  c1SSD->cd(5);
  fCRmeanGraph->DrawCopy();
    TLine* tlin1=new TLine(0.,0.2,1698.,0.2);
    tlin1->SetLineColor(2);
    tlin1->SetLineWidth(2);
    tlin1->Draw("same");
    TLine* tlin2=new TLine(0.,-0.2,1698.,-0.2);
    tlin2->SetLineColor(2);
    tlin2->SetLineWidth(2);
    tlin2->Draw("same");
  TLatex* ta1=new TLatex(0.2,0.8,"SSD Calibration");
  ta1->SetNDC();
  ta1->SetTextSize(0.05);
  ta1->SetTextColor(2);
  ta1->Draw("same");
  c1SSD->cd(6);
  fHistCRmean->DrawCopy();
  c1SSD->Update();
  c1SSD->SaveAs("SSD.pdf");
  pdfFileNames+=" SSD.pdf";
}

//_______________________________________________________________________
void VertexQAMacro(TFile *fildat, TCanvas **&clist, Int_t &cnum){

	TDirectoryFile *dir = (TDirectoryFile*)fildat->Get("Vertex_Performance");
	if(!dir){
		Printf("Vertex directory not found... check!");
	}
	
	TList *lt = (TList*)dir->Get("cOutputVtxESD");
	
	cnum = 1;
	clist = new TCanvas*[1];

	
	TH1F *xVtxSPD = (TH1F*)lt->FindObject("fhSPDVertexX");
	TH1F *yVtxSPD = (TH1F*)lt->FindObject("fhSPDVertexY");
	TH1F *zVtxSPD = (TH1F*)lt->FindObject("fhSPDVertexZ");
	
	TH1F *zVtxSPD_Zonly = (TH1F*)lt->FindObject("fhSPDVertexZonly");
	
	if(!zVtxSPD_Zonly){
		Printf("using SPD 3D histo, Zonly not available");
	    zVtxSPD_Zonly = (TH1F*)lt->FindObject("fhSPDVertexZ");	
	}
	
	TH1F *xVtxTRK = (TH1F*)lt->FindObject("fhTRKVertexX");
	TH1F *yVtxTRK = (TH1F*)lt->FindObject("fhTRKVertexY");
	TH1F *zVtxTRK = (TH1F*)lt->FindObject("fhTRKVertexZ");
	
	TH2F *hntrksSPDvsSPDcls = (TH2F*)lt->FindObject("fhntrksSPDvsSPDcls");
    TH2F *hntrksZvsSPDcls = (TH2F*)lt->FindObject("fhntrksZvsSPDcls");
	
	Bool_t histoCorelation = kTRUE;
	
	if(!hntrksZvsSPDcls){
		Printf("skipping the second part, no 2D histos available");
		histoCorelation=kFALSE;	
	}

	TString ctitle=GetRunNumber()+"TRKandSPD3DxVtx";
	TCanvas *TRK_SPD3D_Vtx = new TCanvas("TRKandSPD3DVtx",ctitle,1000,1000);
	TRK_SPD3D_Vtx->Divide(3,2);
	clist[0]=TRK_SPD3D_Vtx;
	gStyle->SetOptFit(111);

	TRK_SPD3D_Vtx->cd(1);
	xVtxSPD->SetMarkerStyle(20);
	xVtxSPD->SetLineWidth(3);
	xVtxSPD->SetMarkerColor(kBlue+2);
	TF1 *fx = new TF1("gaus", "gaus", -1, 1);
  	xVtxTRK->SetMarkerStyle(20);
	xVtxTRK->SetLineWidth(4);
	xVtxTRK->SetLineColor(2);
	xVtxTRK->Draw("PE");
	xVtxTRK->Fit("gaus", "M");
	xVtxSPD->Draw("PE SAME");
	xVtxTRK->GetXaxis()->SetRangeUser(-0.05, 0.15);
	xVtxSPD->GetXaxis()->SetRangeUser(-0.05, 0.15);
	
	TLatex* tVTX1=new TLatex(0.15,0.85,"VertexSPD");
    tVTX1->SetNDC();
    tVTX1->SetTextColor(kBlue+2);
    tVTX1->Draw();
	TLatex* tVTX2=new TLatex(0.15,0.8,"VertexTRK");
    tVTX2->SetNDC();
    tVTX2->SetTextColor(2);
    tVTX2->Draw();
	
	TRK_SPD3D_Vtx->cd(2);
	yVtxSPD->SetMarkerStyle(20);
	yVtxSPD->SetLineWidth(3);
	yVtxSPD->SetMarkerColor(kBlue+2);
	TF1 *fy = new TF1("gaus", "gaus", -1, 1);
	yVtxTRK->SetMarkerStyle(20);
	yVtxTRK->SetLineWidth(3);
	yVtxTRK->SetLineColor(2);
	yVtxTRK->Draw("PE");
	yVtxTRK->Fit("gaus", "M");
	yVtxSPD->Draw("PE SAME");
	yVtxTRK->GetXaxis()->SetRangeUser(0.15, 0.4);
	yVtxSPD->GetXaxis()->SetRangeUser(0.15, 0.4);

	TLatex* tVTX3=new TLatex(0.15,0.85,"VertexSPD");
  tVTX3->SetNDC();
  tVTX3->SetTextColor(kBlue+2);
  tVTX3->Draw();
	TLatex* tVTX4=new TLatex(0.15,0.8,"VertexTRK");
  tVTX4->SetNDC();
  tVTX4->SetTextColor(2);
  tVTX4->Draw();
 	
	
	
	TRK_SPD3D_Vtx->cd(3);
	
	TF1 *fz = new TF1("gaus", "gaus", -20, 20);
	zVtxTRK->SetMarkerStyle(20);
	zVtxTRK->SetLineWidth(3);
	zVtxTRK->SetMarkerColor(2);
	zVtxTRK->SetLineColor(2);
	zVtxTRK->Draw("PE");
	zVtxTRK->Fit("gaus", "M");
	zVtxSPD->SetMarkerStyle(20);
	zVtxSPD->SetLineWidth(1);
	zVtxSPD->SetLineColor(kBlue+2);
	zVtxSPD->SetMarkerColor(kBlue+2);
	zVtxSPD->SetMarkerSize(0.8);
	zVtxSPD->Draw("PE SAME");
	TLatex* tVTX5=new TLatex(0.15,0.85,"VertexSPD");
  tVTX5->SetNDC();
  tVTX5->SetTextColor(kBlue+2);
  tVTX5->Draw();
	TLatex* tVTX6=new TLatex(0.15,0.8,"VertexTRK");
  tVTX6->SetNDC();
  tVTX6->SetTextColor(2);
  tVTX6->Draw();

	
	
	//	TCanvas *corrContrSPDClusters = new TCanvas("corrContrSPDClusters", "corrContrSPDClusters");
	// corrContrSPDClusters->Divide(3,1);
	//clist[1]=corrContrSPDClusters;
		TRK_SPD3D_Vtx->cd(4);
		//	corrContrSPDClusters->cd(1);
	zVtxSPD_Zonly->SetLineWidth(3);
	zVtxSPD_Zonly->SetLineColor(kBlue+2);
	zVtxSPD_Zonly->Draw();
	TLatex* tVTX7=new TLatex(0.15,0.8,"Vertex Z only");
  tVTX7->SetNDC();
  tVTX7->SetTextColor(2);
  tVTX7->Draw();
	
	if(histoCorelation){
			TRK_SPD3D_Vtx->cd(5);
			//corrContrSPDClusters->cd(2);
		hntrksSPDvsSPDcls->SetMarkerStyle(20);
		hntrksSPDvsSPDcls->Draw();
	
		TRK_SPD3D_Vtx->cd(6);
		//		corrContrSPDClusters->cd(3);
		hntrksZvsSPDcls->SetMarkerStyle(20);
		hntrksZvsSPDcls->Draw();	
   }	
	TRK_SPD3D_Vtx->SaveAs("vertex.pdf");
	pdfFileNames+=" vertex.pdf";
	delete fx;
	delete fy;
	delete fz;
}

//_______________________________________________________________________
void PlotSPD(TFile *fildat, TFile *filMC, TCanvas **&clist, Int_t &cnum){

  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111);
  cnum=2; // number of canvases 
  clist= new TCanvas* [2];//array of pointers to TCanvases

   TDirectoryFile *spddata = (TDirectoryFile*)fildat->Get("SPD_Performance");
    spddata->cd();
   TList *fListData = (TList*)spddata->Get("coutput1");
    
   TString fTitleData = "Data";
   TString fTitleMc = "MC";


  Double_t nevtsData = ((TH1I*)(fListData->FindObject("hEventsProcessed")))->GetEntries();
  printf("   #events in %s : %f \n",fTitleData.Data(),nevtsData);

  TDirectoryFile *spdmc = (TDirectoryFile*)filMC->Get("SPD_Performance");
  spdmc->cd();
  TList *fListMc = (TList*)spdmc->Get("coutput1");
  Double_t nevtsMc = ((TH1I*)(fListMc->FindObject("hEventsProcessed")))->GetEntries();
  printf("   #events in %s : %f \n",fTitleMc.Data(),nevtsMc);

  TH2F *trackData = (TH2F*)fListData->FindObject("hSPDphivsSPDeta");
  trackData->SetTitle(Form("%s %s",trackData->GetTitle(),fTitleData.Data()));
  TH1D *trackDataPhi = trackData->ProjectionY();
  if(!trackDataPhi) printf("NO 1 \n");

  trackDataPhi->SetTitle("Tracklets vs Phi");
  TH1D *trackDataEta = trackData->ProjectionX();
  if(!trackDataEta) printf("NO 2 \n");
  trackDataEta->SetTitle("Tracklets vs eta");

  TH1F etaData, phiData;
  trackDataEta->Copy(etaData);
  trackDataPhi->Copy(phiData);

  TH1F etaFrac, phiFrac, mcEta, mcPhi;
  trackDataEta->Copy(etaFrac);
  trackDataPhi->Copy(phiFrac);

  TH2F *trackMc = (TH2F*)fListMc->FindObject("hSPDphivsSPDeta");
    trackMc->SetTitle(Form("%s %s",trackMc->GetTitle(),fTitleMc.Data()));

  TString ctitle = GetRunNumber()+"tracklets";
  TCanvas *tracklets = new TCanvas("tracklets",ctitle,1200,600);
  clist[0]=tracklets;
  tracklets->Divide(2,1);
  tracklets->cd(1);
  tracklets->cd(1)->SetRightMargin(0.15);
  trackData->SetTitle(Form("Run %d",gRunNumber));
  trackData->DrawCopy("colz");
  tracklets->cd(2);
  tracklets->cd(2)->SetRightMargin(0.15);
  trackMc->SetTitle(Form("Run %d",gRunNumberMC));
  TH1D *h = (TH1D*)trackMc->DrawCopy("colz");
  tracklets->SaveAs("SPDtracklets.pdf");
  pdfFileNames+=" SPDtracklets.pdf";

  TH1D *trackMcPhi = trackMc->ProjectionY();
  trackMcPhi->SetTitle(Form("%s",h->GetTitle()));
  TH1D *trackMcEta = trackMc->ProjectionX();
  trackMcEta->SetTitle(Form("%s",h->GetTitle()));

  TH1F etaMc, phiMc;
  trackMcEta->Copy(etaMc);
  trackMcPhi->Copy(phiMc);

  trackMcEta->Copy(mcEta);
  trackMcPhi->Copy(mcPhi);

  etaFrac.Scale(1./etaFrac.GetEntries());
  mcEta.Scale(1./mcEta.GetEntries());
  etaFrac.Add(&mcEta,-1);
  etaFrac.Divide(&mcEta);

  phiFrac.Scale(1./phiFrac.GetEntries());
  mcPhi.Scale(1./mcPhi.GetEntries());
  phiFrac.Add(&mcPhi,-1);
  phiFrac.Divide(&mcPhi);

  ctitle = GetRunNumber()+"tracklets and ratios vs eta and phi";
  TCanvas *track = new TCanvas("track",ctitle,1200,1200);
  clist[1]=track;
  track->Divide(2,2);
  track->cd(1);
  phiData.SetLineColor(kRed);
  phiData.SetLineWidth(2);
  phiData.Scale(1./phiData.GetEntries());
  phiData.DrawCopy();
  phiMc.Scale(1./phiMc.GetEntries());
  TLatex* tphi=new TLatex(0.5,0.85,Form("Red = %d; Blue = %d",gRunNumber,gRunNumberMC));
  tphi->SetNDC();
  tphi->SetTextSize(0.04);
   tphi->SetTextColor(1);
  tphi->Draw();  
  phiMc.DrawCopy("same");
  track->cd(2);
  etaData.SetLineColor(kRed);
  etaData.SetLineWidth(2);
  etaData.Scale(1./etaData.GetEntries());
  etaData.DrawCopy();
  etaMc.Scale(1./etaMc.GetEntries());
  tphi->Draw();  
  etaMc.DrawCopy("same");
  track->cd(3);
  phiFrac.SetLineColor(1);
  TLatex* tratio=new TLatex(0.2,0.85,Form("Run %d / Run %d",gRunNumber,gRunNumberMC));
  tratio->SetNDC();
  tratio->SetTextColor(1);
  phiFrac.DrawCopy();
  tratio->Draw();  
  track->cd(4);
  etaFrac.SetLineColor(1);
  etaFrac.DrawCopy();  
  tratio->Draw();  
  track->SaveAs("SPD_eta_phi.pdf");
  pdfFileNames+=" SPD_eta_phi.pdf";
}

//_______________________________________________________________________
Bool_t PlotITSTPCMatchingEff(TFile *f, TCanvas**& clist,Int_t& cnum) {

  cnum=1;
  clist = new TCanvas*[1];

  //  clist = new TCanvas* [1];
  TString ctitle = GetRunNumber()+"ITS-TPC match";
  TCanvas* cITSTPCmatch = new TCanvas("cITSTPCmatch",ctitle,10,10,1200,600);
  clist[0]=cITSTPCmatch;
  cITSTPCmatch->Divide(2,1);
  cITSTPCmatch->cd(1);
  gPad->SetGridy();
  gPad->SetLogx();
  cITSTPCmatch->cd(2);
  gPad->SetGridy();
  gPad->SetLogx();

  //  clist = cITSTPCmatch;

  if(!f) return kFALSE;

  TList *list=0;
  TList *listSPD=0;
  TDirectoryFile *dir=0;

  // count active SPD HSs
  dir=(TDirectoryFile*)f->GetDirectory("SPD_Performance");
  if(dir) listSPD = (TList*)dir->Get("coutput1");
  if(!dir) return kFALSE;

  Float_t spdFrac[2]={0.,0.};
  TH1F *hnHSsSPD=new TH1F("hnHSsSPD","Active HSs in SPD layers 1 and 2; layer; HSs",2,0.5,2.5);
  if(listSPD) {
    //listSPD->Print();
    TH1F *hFiredChip = (TH1F*)listSPD->FindObject("hFiredChip");
    Int_t nHSsInner=0,nHSsOuter=0;
    for(Int_t i=0;i<400;i++) if(hFiredChip->GetBinContent(i)>0) nHSsInner++;
    for(Int_t i=400;i<1200;i++) if(hFiredChip->GetBinContent(i)>0) nHSsOuter++;
    nHSsInner = (Int_t)(nHSsInner/10);
    nHSsOuter = (Int_t)(nHSsOuter/10);
    hnHSsSPD->SetBinContent(1,nHSsInner);
    hnHSsSPD->SetBinContent(2,nHSsOuter);
    spdFrac[0]=(Float_t)nHSsInner/40.;
    spdFrac[1]=(Float_t)nHSsOuter/80.;
  }
  TGraph *spdFrac0=new TGraph(1);
  spdFrac0->SetPoint(0,0.08,spdFrac[0]);
  spdFrac0->SetMarkerColor(1); spdFrac0->SetMarkerStyle(20);
  TGraph *spdFrac1=new TGraph(1);
  spdFrac1->SetPoint(0,0.08,spdFrac[1]);
  spdFrac1->SetMarkerColor(1); spdFrac1->SetMarkerStyle(24);
  TLegend *l2=new TLegend(0.1,0.62,0.5,0.93);
  l2->SetBorderSize(1);
  l2->AddEntry(spdFrac0,"Frac. active SPD0","p");
  l2->AddEntry(spdFrac1,"Frac. active SPD1","p");

  //
  // Efficiencies for CENTRAL
  //
  dir=(TDirectoryFile*)f->GetDirectory("ITS_Performance");
  if(dir) list = (TList*)dir->Get("cOutputITS_3500_10000");
  if(!list) return kFALSE;

  TH1F *fHistPtTPCInAcc = (TH1F*)list->FindObject("fHistPtTPCInAcc");
  TH1F *fHistPtITSMI6InAcc = (TH1F*)list->FindObject("fHistPtITSMI6InAcc");
  TH1F *fHistPtITSMI5InAcc = (TH1F*)list->FindObject("fHistPtITSMI5InAcc");
  TH1F *fHistPtITSMI4InAcc = (TH1F*)list->FindObject("fHistPtITSMI4InAcc");
  TH1F *fHistPtITSMI3InAcc = (TH1F*)list->FindObject("fHistPtITSMI3InAcc");
  TH1F *fHistPtITSMI2InAcc = (TH1F*)list->FindObject("fHistPtITSMI2InAcc");
  TH1F *fHistPtITSMISPDInAcc = (TH1F*)list->FindObject("fHistPtITSMISPDInAcc");
  TH1F *fHistPtITSMIoneSPDInAcc = (TH1F*)list->FindObject("fHistPtITSMIoneSPDInAcc");
  TH1F *fHistPtITSTPCsel = (TH1F*)list->FindObject("fHistPtITSTPCsel");
  TH1F *fHistPtITSMIge2InAcc = (TH1F*)fHistPtITSMI6InAcc->Clone("fHistPtITSMIge2InAcc");
  fHistPtITSMIge2InAcc->Add(fHistPtITSMI5InAcc);
  fHistPtITSMIge2InAcc->Add(fHistPtITSMI4InAcc);
  fHistPtITSMIge2InAcc->Add(fHistPtITSMI3InAcc);
  fHistPtITSMIge2InAcc->Add(fHistPtITSMI2InAcc);


  TLegend *l3=new TLegend(0.5,0.62,0.95,0.93);
  l3->SetBorderSize(1);
  cITSTPCmatch->cd(1);
  fHistPtITSMIge2InAcc->SetTitle("Fraction of prolonged tracks with N ITS points: central");
  fHistPtITSMIge2InAcc->SetYTitle("ITS+TPC / TPC");
  fHistPtITSMIge2InAcc->Divide(fHistPtITSMIge2InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMIge2InAcc->SetMaximum(1.6);
  fHistPtITSMIge2InAcc->SetMinimum(0);
  fHistPtITSMIge2InAcc->GetXaxis()->SetRangeUser(0.1,30);
  fHistPtITSMIge2InAcc->Draw();
  l3->AddEntry(fHistPtITSMIge2InAcc,">=2 cls","l");
  fHistPtITSMI6InAcc->Divide(fHistPtITSMI6InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMI6InAcc->SetLineColor(2);
  l3->AddEntry(fHistPtITSMI6InAcc,"6 cls","l");
  fHistPtITSMI6InAcc->Draw("same");
  fHistPtITSMI5InAcc->Divide(fHistPtITSMI5InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMI5InAcc->SetLineColor(3);
  l3->AddEntry(fHistPtITSMI5InAcc,"5 cls","l");
  fHistPtITSMI5InAcc->Draw("same");
  fHistPtITSMI4InAcc->Divide(fHistPtITSMI4InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMI4InAcc->SetLineColor(4);
  l3->AddEntry(fHistPtITSMI4InAcc,"4 cls","l");
  fHistPtITSMI4InAcc->Draw("same");
  fHistPtITSMI3InAcc->Divide(fHistPtITSMI3InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMI3InAcc->SetLineColor(6);
  l3->AddEntry(fHistPtITSMI3InAcc,"3 cls","l");
  fHistPtITSMI3InAcc->Draw("same");
  fHistPtITSMI2InAcc->Divide(fHistPtITSMI2InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMI2InAcc->SetLineColor(7);
  l3->AddEntry(fHistPtITSMI2InAcc,"2 cls","l");
  fHistPtITSMI2InAcc->Draw("same");
  fHistPtITSMISPDInAcc->Divide(fHistPtITSMISPDInAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMISPDInAcc->SetLineColor(9);
  l3->AddEntry(fHistPtITSMISPDInAcc,"2SPD + any","l");
  fHistPtITSMISPDInAcc->Draw("same");
  fHistPtITSMIoneSPDInAcc->Divide(fHistPtITSMIoneSPDInAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMIoneSPDInAcc->SetLineColor(15);
  l3->AddEntry(fHistPtITSMIoneSPDInAcc,">=1SPD + any","l");
  fHistPtITSMIoneSPDInAcc->Draw("same");
  fHistPtITSTPCsel->Divide(fHistPtITSTPCsel,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSTPCsel->SetLineColor(kAzure+1);
  l3->AddEntry(fHistPtITSTPCsel,">=1SPD + any + d_{0} cut","l");
  fHistPtITSTPCsel->Draw("same");
  fHistPtITSMIge2InAcc->Draw("same");
  l3->Draw();
  l2->Draw();
  spdFrac0->Draw("p");
  spdFrac1->Draw("p");

  //
  // Efficiencies for PERIPHERAL
  //
  dir=(TDirectoryFile*)f->GetDirectory("ITS_Performance");
  if(dir) list = (TList*)dir->Get("cOutputITS_70_310");
  if(!list) return kFALSE;

  fHistPtTPCInAcc = (TH1F*)list->FindObject("fHistPtTPCInAcc");
  fHistPtITSMI6InAcc = (TH1F*)list->FindObject("fHistPtITSMI6InAcc");
  fHistPtITSMI5InAcc = (TH1F*)list->FindObject("fHistPtITSMI5InAcc");
  fHistPtITSMI4InAcc = (TH1F*)list->FindObject("fHistPtITSMI4InAcc");
  fHistPtITSMI3InAcc = (TH1F*)list->FindObject("fHistPtITSMI3InAcc");
  fHistPtITSMI2InAcc = (TH1F*)list->FindObject("fHistPtITSMI2InAcc");
  fHistPtITSMISPDInAcc = (TH1F*)list->FindObject("fHistPtITSMISPDInAcc");
  fHistPtITSMIoneSPDInAcc = (TH1F*)list->FindObject("fHistPtITSMIoneSPDInAcc");
  fHistPtITSTPCsel = (TH1F*)list->FindObject("fHistPtITSTPCsel");
  fHistPtITSMIge2InAcc = (TH1F*)fHistPtITSMI6InAcc->Clone("fHistPtITSMIge2InAcc");
  fHistPtITSMIge2InAcc->Add(fHistPtITSMI5InAcc);
  fHistPtITSMIge2InAcc->Add(fHistPtITSMI4InAcc);
  fHistPtITSMIge2InAcc->Add(fHistPtITSMI3InAcc);
  fHistPtITSMIge2InAcc->Add(fHistPtITSMI2InAcc);


  cITSTPCmatch->cd(2);
  fHistPtITSMIge2InAcc->SetTitle("Fraction of prolonged tracks with N ITS points: peripheral");
  fHistPtITSMIge2InAcc->SetYTitle("ITS+TPC / TPC");
  fHistPtITSMIge2InAcc->Divide(fHistPtITSMIge2InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMIge2InAcc->SetMaximum(1.6);
  fHistPtITSMIge2InAcc->SetMinimum(0);
  fHistPtITSMIge2InAcc->GetXaxis()->SetRangeUser(0.1,30);
  fHistPtITSMIge2InAcc->Draw();
  fHistPtITSMI6InAcc->Divide(fHistPtITSMI6InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMI6InAcc->SetLineColor(2);
  fHistPtITSMI6InAcc->Draw("same");
  fHistPtITSMI5InAcc->Divide(fHistPtITSMI5InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMI5InAcc->SetLineColor(3);
  fHistPtITSMI5InAcc->Draw("same");
  fHistPtITSMI4InAcc->Divide(fHistPtITSMI4InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMI4InAcc->SetLineColor(4);
  fHistPtITSMI4InAcc->Draw("same");
  fHistPtITSMI3InAcc->Divide(fHistPtITSMI3InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMI3InAcc->SetLineColor(6);
  fHistPtITSMI3InAcc->Draw("same");
  fHistPtITSMI2InAcc->Divide(fHistPtITSMI2InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMI2InAcc->SetLineColor(7);
  fHistPtITSMI2InAcc->Draw("same");
  fHistPtITSMISPDInAcc->Divide(fHistPtITSMISPDInAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMISPDInAcc->SetLineColor(9);
  fHistPtITSMISPDInAcc->Draw("same");
  fHistPtITSMIoneSPDInAcc->Divide(fHistPtITSMIoneSPDInAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMIoneSPDInAcc->SetLineColor(15);
  fHistPtITSMIoneSPDInAcc->Draw("same");
  fHistPtITSTPCsel->Divide(fHistPtITSTPCsel,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSTPCsel->SetLineColor(kAzure+1);
  fHistPtITSTPCsel->Draw("same");
  fHistPtITSMIge2InAcc->Draw("same");
  l3->Draw();
  l2->Draw();
  spdFrac0->Draw("p");
  spdFrac1->Draw("p");
  cITSTPCmatch->SaveAs("TPCITSmatching.pdf");
  pdfFileNames+=" TPCITSmatching.pdf";
  return kTRUE;
}

//_______________________________________________________________________
void SaveC(TFile &fout, TCanvas**& clist, Int_t cnum){
  TDirectory *current = gDirectory;
  fout.cd();
  for(Int_t i=0;i<cnum;i++)clist[i]->Write();
  delete[] clist;
  current->cd();
}

//_______________________________________________________________________
TString GetRunNumber(){
  // returns a string with the run number
  char rn[10];
  sprintf(rn,"%d  ",gRunNumber);
  TString str(rn);
  return str;
}

