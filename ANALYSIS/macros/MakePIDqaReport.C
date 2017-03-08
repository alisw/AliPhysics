#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLine.h"
#include "TList.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TObjArray.h"
#include "TString.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TF1.h"
#include "TPDF.h"
#include "TColor.h"
#include "AliPID.h"
#endif

void SetupStyle();
TH2* Get2DHistogramfromList(TList *pidqalist, const char* listname, const char* histoname);
void AddFit(TH2* h2d);
void PublishCanvas(TList *qaList, const char* det, const char* name, TString nadd="");
void SetupPadStyle();

void LoadLibs();
Int_t CheckLoadLibrary(const char* library);

TCanvas *fCanvas=0x0;

/*

Example (require aliroot environment)

root.exe -l -b -q $ALICE_ROOT/ANALYSIS/macros/MakePIDqaReport.C'("PIDqa.root")'

*/

void MakePIDqaReport(const char* inputFile, const char* outputFile="PIDqaReport.pdf", TString dirInFile = "")
{
  //
  // Make a pdf file with the efficiency report
  //

  LoadLibs();
  SetupStyle();

  TFile f(inputFile);
  if (!f.IsOpen()){
    printf("Could not open file '%s'\n",f.GetName());
    return;
  }
  
  TString listName = "PIDqa";
  if (dirInFile != "")
    listName = listName.Prepend(Form("%s/", dirInFile.Data()));
  
  printf("%s", listName.Data());
  TList *qaList = (TList*) f.Get(listName.Data());
  if (!qaList){
    printf("Could not find list '%s' in file '%s'\n",listName.Data(), f.GetName());
    return;
  }

  fCanvas=new TCanvas;

  TPDF p(outputFile);

  //
  // Invariant mass plots
  //


  //
  // Make QA info
  //

  // ITS PID
  PublishCanvas(qaList,"ITS","hNsigmaP_ITS_%s");

  // TPC PID
  TList *qaListTPC = (TList*)qaList->FindObject("TPC");
  if (qaListTPC){
    PublishCanvas(qaListTPC,"TPCBasic","hNsigmaP_TPC_Basic_%s");
    PublishCanvas(qaListTPC,"TPCV0","hNsigmaP_TPC_V0_%s");
    //   if (man->GetCurrentPeriod()=="11h"){
      //     PublishCanvas(qaListTPC,"TPC","hNsigmaP_TPC_Basic_%s_Hybrid","Hybrid");
    //     PublishCanvas(qaListTPC,"TPC","hNsigmaP_TPC_Basic_%s_OROChigh","OROChigh");
    //   }
  }
  else {
    printf("Could not find list '%s/TPC' in file '%s'\n", listName.Data(), f.GetName());
  }

  // TPC PID after 3 sigma TOF cut
  PublishCanvas(qaList,"TPC_TOF","hNsigmaP_TPC_TOF_%s");

  // TOF PID
  PublishCanvas(qaList,"TOF","hNsigmaP_TOF_%s");

  // TRD PID
  TList *qaListTRD = (TList*)qaList->FindObject("TRD");
  if (qaListTRD) PublishCanvas(qaListTRD,"TRDTruncatedMean","hNsigmaP_TOFTPC_TRD_TruncatedMean_%s");
  TList *qaListTRDV0 = (TList*)qaListTRD->FindObject("TRDV0");
  if (qaListTRDV0) PublishCanvas(qaListTRDV0,"TRDTruncatedMeanV0","hNsigmaP_TOFTPC_TRD_TruncatedMeanV0_%s");

  if (qaListTRD){

      
      TH2 *hLikeP_TRD_4tls_electron=Get2DHistogramfromList(qaListTRD,"TRDLikelihood","hLike1DP_TRD_4tls_electron_Likelihood");
      TH2 *hLikeP_TRD_4tls_pion=Get2DHistogramfromList(qaListTRD,"TRDLikelihood","hLike1DP_TRD_4tls_pion_Likelihood");
      TH2 *hLikeP_TRD_5tls_electron=Get2DHistogramfromList(qaListTRD,"TRDLikelihood","hLike1DP_TRD_5tls_electron_Likelihood");
      TH2 *hLikeP_TRD_5tls_pion=Get2DHistogramfromList(qaListTRD,"TRDLikelihood","hLike1DP_TRD_5tls_pion_Likelihood");
      TH2 *hLikeP_TRD_6tls_electron=Get2DHistogramfromList(qaListTRD,"TRDLikelihood","hLike1DP_TRD_6tls_electron_Likelihood");
      TH2 *hLikeP_TRD_6tls_pion=Get2DHistogramfromList(qaListTRD,"TRDLikelihood","hLike1DP_TRD_6tls_pion_Likelihood");

      fCanvas->Divide(2,3);
      fCanvas->cd(1);
      gPad->SetLogz();
      if(hLikeP_TRD_4tls_electron) hLikeP_TRD_4tls_electron->Draw("colz");
      fCanvas->cd(2);
      gPad->SetLogz();
      if(hLikeP_TRD_4tls_pion) hLikeP_TRD_4tls_pion->Draw("colz");
      fCanvas->cd(3);
      gPad->SetLogz();
      if(hLikeP_TRD_5tls_electron) hLikeP_TRD_5tls_electron->Draw("colz");
      fCanvas->cd(4);
      gPad->SetLogz();
      if(hLikeP_TRD_5tls_pion) hLikeP_TRD_5tls_pion->Draw("colz");
      fCanvas->cd(5);
      gPad->SetLogz();
      if(hLikeP_TRD_6tls_electron) hLikeP_TRD_6tls_electron->Draw("colz");
      fCanvas->cd(6);
      gPad->SetLogz();
      if(hLikeP_TRD_6tls_pion) hLikeP_TRD_6tls_pion->Draw("colz");

      fCanvas->Update();
      fCanvas->Clear();

      fCanvas->Divide(2,3);
      TH2 *hTPCnsigmaPnoTRD=Get2DHistogramfromList(qaListTRD,"TRDLikelihood","hTPCnsigmaPnoTRD_Likelihood");
      TH2 *hTPCnsigmaPTRD1D_5tls=Get2DHistogramfromList(qaListTRD,"TRDLikelihood","hTPCnsigmaPLQ1D_5tls_Likelihood");
      TH2 *hTPCnsigmaPTRD1D_6tls=Get2DHistogramfromList(qaListTRD,"TRDLikelihood","hTPCnsigmaPLQ1D_6tls_Likelihood");
      TH2 *hTPCnsigmaPTRD2D_5tls=Get2DHistogramfromList(qaListTRD,"TRDLikelihood","hTPCnsigmaPLQ2D_5tls_Likelihood");
      TH2 *hTPCnsigmaPTRD2D_6tls=Get2DHistogramfromList(qaListTRD,"TRDLikelihood","hTPCnsigmaPLQ2D_6tls_Likelihood");

      if(hTPCnsigmaPnoTRD && hTPCnsigmaPTRD1D_5tls && hTPCnsigmaPTRD1D_6tls && hTPCnsigmaPTRD2D_5tls && hTPCnsigmaPTRD2D_6tls){
	  fCanvas->cd(1);
	  TPaveText pt(.1,.1,.9,.9,"NDC");
	  pt.SetBorderSize(1);
	  pt.SetFillColor(0);
	  pt.SetTextSizePixels(10);
	  pt.AddText("ratio TPC nsigma before/after TRD PID");
          pt.Draw();
	  fCanvas->cd(2);
	  gPad->SetLogz();
	  hTPCnsigmaPTRD1D_5tls->Divide(hTPCnsigmaPnoTRD);
	  hTPCnsigmaPTRD1D_5tls->Draw("colz");
	  fCanvas->cd(3);
	  gPad->SetLogz();
	  hTPCnsigmaPTRD1D_6tls->Divide(hTPCnsigmaPnoTRD);
	  hTPCnsigmaPTRD1D_6tls->Draw("colz");
	  fCanvas->cd(4);
	  gPad->SetLogz();
	  hTPCnsigmaPTRD2D_5tls->Divide(hTPCnsigmaPnoTRD);
	  hTPCnsigmaPTRD2D_5tls->Draw("colz");
	  fCanvas->cd(5);
	  gPad->SetLogz();
	  hTPCnsigmaPTRD2D_6tls->Divide(hTPCnsigmaPnoTRD);
	  hTPCnsigmaPTRD2D_6tls->Draw("colz");
	  fCanvas->Update();
	  fCanvas->Clear();
      }

  }

  if (qaListTRDV0){

      fCanvas->Divide(2,3);
      TH2 *hLikeP_TRDV0_4tls_electron=Get2DHistogramfromList(qaListTRDV0,"TRDLikelihoodV0","hLike1DP_TRD_4tls_electron_LikelihoodV0");
      TH2 *hLikeP_TRDV0_4tls_pion=Get2DHistogramfromList(qaListTRDV0,"TRDLikelihoodV0","hLike1DP_TRD_4tls_pion_LikelihoodV0");
      TH2 *hLikeP_TRDV0_5tls_electron=Get2DHistogramfromList(qaListTRDV0,"TRDLikelihoodV0","hLike1DP_TRD_5tls_electron_LikelihoodV0");
      TH2 *hLikeP_TRDV0_5tls_pion=Get2DHistogramfromList(qaListTRDV0,"TRDLikelihoodV0","hLike1DP_TRD_5tls_pion_LikelihoodV0");
      TH2 *hLikeP_TRDV0_6tls_electron=Get2DHistogramfromList(qaListTRDV0,"TRDLikelihoodV0","hLike1DP_TRD_6tls_electron_LikelihoodV0");
      TH2 *hLikeP_TRDV0_6tls_pion=Get2DHistogramfromList(qaListTRDV0,"TRDLikelihoodV0","hLike1DP_TRD_6tls_pion_LikelihoodV0");

      fCanvas->cd(1);
      gPad->SetLogz();
      if(hLikeP_TRDV0_4tls_electron) hLikeP_TRDV0_4tls_electron->Draw("colz");
      fCanvas->cd(2);
      gPad->SetLogz();
      if(hLikeP_TRDV0_4tls_pion) hLikeP_TRDV0_4tls_pion->Draw("colz");
      fCanvas->cd(3);
      gPad->SetLogz();
      if(hLikeP_TRDV0_5tls_electron) hLikeP_TRDV0_5tls_electron->Draw("colz");
      fCanvas->cd(4);
      gPad->SetLogz();
      if(hLikeP_TRDV0_5tls_pion) hLikeP_TRDV0_5tls_pion->Draw("colz");
      fCanvas->cd(5);
      gPad->SetLogz();
      if(hLikeP_TRDV0_6tls_electron) hLikeP_TRDV0_6tls_electron->Draw("colz");
      fCanvas->cd(6);
      gPad->SetLogz();
      if(hLikeP_TRDV0_6tls_pion) hLikeP_TRDV0_6tls_pion->Draw("colz");

      fCanvas->Update();
      fCanvas->Clear();

      TCanvas* fCanvas=new TCanvas;
      fCanvas->Divide(2,3);
      TH2 *hTPCnsigmaPnoTRDV0=Get2DHistogramfromList(qaListTRDV0,"TRDLikelihoodV0","hTPCnsigmaPnoTRD_LikelihoodV0");
      TH2 *hTPCnsigmaPTRD1DV0_5tls=Get2DHistogramfromList(qaListTRDV0,"TRDLikelihoodV0","hTPCnsigmaPLQ1D_5tls_LikelihoodV0");
      TH2 *hTPCnsigmaPTRD1DV0_6tls=Get2DHistogramfromList(qaListTRDV0,"TRDLikelihoodV0","hTPCnsigmaPLQ1D_6tls_LikelihoodV0");
      TH2 *hTPCnsigmaPTRD2DV0_5tls=Get2DHistogramfromList(qaListTRDV0,"TRDLikelihoodV0","hTPCnsigmaPLQ2D_5tls_LikelihoodV0");
      TH2 *hTPCnsigmaPTRD2DV0_6tls=Get2DHistogramfromList(qaListTRDV0,"TRDLikelihoodV0","hTPCnsigmaPLQ2D_6tls_LikelihoodV0");

      if(hTPCnsigmaPnoTRDV0 && hTPCnsigmaPTRD1DV0_5tls && hTPCnsigmaPTRD1DV0_6tls && hTPCnsigmaPTRD2DV0_5tls && hTPCnsigmaPTRD2DV0_6tls){
	  fCanvas->cd(1);
	  TPaveText pt(.1,.1,.9,.9,"NDC");
	  pt.SetBorderSize(1);
	  pt.SetFillColor(0);
	  pt.SetTextSizePixels(10);
	  pt.AddText("V0 ratio TPC nsigma before/after TRD PID");
          pt.Draw();
	  fCanvas->cd(2);
	  gPad->SetLogz();
	  hTPCnsigmaPTRD1DV0_5tls->Divide(hTPCnsigmaPnoTRDV0);
	  hTPCnsigmaPTRD1DV0_5tls->Draw("colz");
	  fCanvas->cd(3);
	  gPad->SetLogz();
	  hTPCnsigmaPTRD1DV0_6tls->Divide(hTPCnsigmaPnoTRDV0);
	  hTPCnsigmaPTRD1DV0_6tls->Draw("colz");
	  fCanvas->cd(4);
	  gPad->SetLogz();
	  hTPCnsigmaPTRD2DV0_5tls->Divide(hTPCnsigmaPnoTRDV0);
	  hTPCnsigmaPTRD2DV0_5tls->Draw("colz");
	  fCanvas->cd(5);
	  gPad->SetLogz();
	  hTPCnsigmaPTRD2DV0_6tls->Divide(hTPCnsigmaPnoTRDV0);
	  hTPCnsigmaPTRD2DV0_6tls->Draw("colz");
	  fCanvas->Update();
	  fCanvas->Clear();
      }
  }
 


  // TPC Response info
  TObjArray *qaInfo=(TObjArray*)qaList->FindObject("QAinfo");
  TObjArray *tpcInfo=0x0;
  if (qaInfo && (tpcInfo=(TObjArray*)qaInfo->FindObject("TPC_info"))){
    TObjArray *tpcSplineInfo=(TObjArray*)tpcInfo->FindObject("TPC_spline_names");
    TObjArray *tpcConfigInfo=(TObjArray*)tpcInfo->FindObject("TPC_config_info");
    fCanvas->Divide(1,2);
  
    TPaveText pt(.1,.1,.9,.9,"NDC");
    pt.SetBorderSize(1);
    pt.SetFillColor(0);
    pt.SetTextSizePixels(16);

    if (tpcSplineInfo){
      for (Int_t i=0; i<tpcSplineInfo->GetEntriesFast();++i) pt.AddText(tpcSplineInfo->At(i)->GetName());
    }
    
    TPaveText pt2(.1,.1,.9,.9,"NDC");
    pt2.SetBorderSize(1);
    pt2.SetFillColor(0);
    pt2.SetTextSizePixels(16);
    if (tpcConfigInfo){
      for (Int_t i=0; i<tpcConfigInfo->GetEntriesFast();++i) pt2.AddText(tpcConfigInfo->At(i)->GetName());
    }

    fCanvas->cd(1);
    pt.Draw();
    fCanvas->cd(2);
    pt2.Draw();
    
    fCanvas->Update();
    fCanvas->Clear();
  }
  
  delete qaList;

  p.Close();
  delete fCanvas;
}

void SetupStyle()
{
  const Int_t NCont=255;

  TStyle *st = new TStyle("mystyle","mystyle");
  gROOT->GetStyle("Plain")->Copy((*st));
  st->SetTitleX(0.1);
  st->SetTitleW(0.8);
  st->SetTitleH(0.08);
  st->SetStatX(.9);
  st->SetStatY(.9);
  st->SetNumberContours(NCont);
  st->SetPalette(1,0);
  st->SetOptStat("erm");
  st->SetOptFit(0);
  st->SetGridColor(kGray+1);
  st->SetPadGridX(kTRUE);
  st->SetPadGridY(kTRUE);
  st->SetPadTickX(kTRUE);
  st->SetPadTickY(kTRUE);
  st->cd();

  const Int_t NRGBs = 5;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };

  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);

}

TH2* Get2DHistogramfromList(TList *pidqalist, const char* listname, const char* histoname)
{
  TList *histolist = (TList *)pidqalist->FindObject(listname);
  if (!histolist) {printf(" list not found \n");  return 0x0; }
  TH2* histo = (TH2*)histolist->FindObject(histoname);
  //   if (!histo) {printf(" histogram not found \n");  return 0x0; }
  return histo;
}

void AddFit(TH2* h2d)
{
  //
  // Fit in slices and draw mean and sigma
  //
  TF1 *f1 = new TF1("f1", "gaus");
  f1->SetRange(-1.5,1.5);
  TObjArray aSlices;
  h2d->FitSlicesY(f1, 0,-1, 0, "QNR",&aSlices);
  aSlices.SetOwner(1);
  TH1* hMean=(TH1*)aSlices.At(1);
  TH1* hSigma=(TH1*)aSlices.At(2);
  TH1* hChi2=(TH1*)aSlices.At(3);
  hChi2->Scale(1./10.);
  aSlices.AddAt(0x0,1);
  aSlices.AddAt(0x0,2);
  aSlices.AddAt(0x0,3);
  hMean->SetMarkerStyle(20);
  hMean->SetMarkerSize(0.3);
  hMean->SetOption("same");
  h2d->GetListOfFunctions()->Add(hMean);
  hSigma->SetMarkerStyle(20);
  hSigma->SetMarkerSize(0.3);
  hSigma->SetOption("same");
  hSigma->SetMarkerColor(kMagenta);
  h2d->GetListOfFunctions()->Add(hSigma);
  hChi2->SetOption("same");
  hChi2->SetMarkerColor(kMagenta + 2);
  hChi2->SetLineColor(kMagenta + 2);
  h2d->GetListOfFunctions()->Add(hChi2);

  TLine *l=0x0;
  l=new TLine(h2d->GetXaxis()->GetXmin(),0,h2d->GetXaxis()->GetXmax(),0);
  l->SetLineStyle(2);
  h2d->GetListOfFunctions()->Add(l);
  l=new TLine(h2d->GetXaxis()->GetXmin(),1,h2d->GetXaxis()->GetXmax(),1);
  l->SetLineStyle(2);
  h2d->GetListOfFunctions()->Add(l);
}

void PublishCanvas(TList *qaList, const char* det, const char* name, TString nadd)
{
  //
  // draw all nSigma + signal histo
  //


  TObjArray arrHistos;

  TPaveText pt(.1,.1,.9,.9,"NDC");
  pt.SetBorderSize(1);
  pt.SetFillColor(0);
  pt.SetTextSizePixels(16);
  pt.AddText(Form("%s PID QA",det));
  if (!nadd.IsNull()){
    pt.AddText(nadd.Data());
    nadd.Prepend("_");
  }
  arrHistos.Add(&pt);

  TH2 *hSig=Get2DHistogramfromList(qaList,det,Form("hSigP_%s",det));
  if (hSig){
    hSig->SetOption("colz");
    arrHistos.Add(hSig);
  }

  for (Int_t i=0;i<AliPID::kSPECIESC;++i){
    //     for (Int_t i=0;i<AliPID::kSPECIES;++i){
    if (i==(Int_t)AliPID::kMuon) continue;
    TH2 *h=Get2DHistogramfromList(qaList,det,Form(name,AliPID::ParticleName(i)));
    if (!h) continue;
    h->SetOption("colz");
    AddFit(h);
    arrHistos.Add(h);
  }

  Int_t nPads=arrHistos.GetEntriesFast();
  Int_t nCols = (Int_t)TMath::Ceil( TMath::Sqrt(nPads) );
  Int_t nRows = (Int_t)TMath::Ceil( (Double_t)nPads/(Double_t)nCols );

  
  fCanvas->Divide(nCols,nRows);


  for (Int_t i=0; i<nPads;++i) {
    fCanvas->cd(i+1);
    SetupPadStyle();
    arrHistos.At(i)->Draw();
  }

  fCanvas->Update();
  fCanvas->Clear();

}

void SetupPadStyle()
{
  gPad->SetLogx();
  gPad->SetLogz();
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetTickx();
  gPad->SetTicky();
}

void LoadLibs()
{
  CheckLoadLibrary("libCore");
  CheckLoadLibrary("libPhysics");
  CheckLoadLibrary("libMinuit");
  CheckLoadLibrary("libGui");
  CheckLoadLibrary("libXMLParser");
  
  CheckLoadLibrary("libGeom");
  CheckLoadLibrary("libVMC");
  
  CheckLoadLibrary("libNet");
  CheckLoadLibrary("libTree");
  CheckLoadLibrary("libProof");
  
  CheckLoadLibrary("libSTEERBase");
  CheckLoadLibrary("libESD");
  CheckLoadLibrary("libCDB");
  CheckLoadLibrary("libRAWDatabase");
  CheckLoadLibrary("libRAWDatarec");
  CheckLoadLibrary("libANALYSIS");
  CheckLoadLibrary("libSTEER");
  
  CheckLoadLibrary("libSTAT");
  
  CheckLoadLibrary("libAOD");
  CheckLoadLibrary("libOADB");
  CheckLoadLibrary("libANALYSISalice");
  CheckLoadLibrary("libCORRFW");
  
  
  CheckLoadLibrary("libTPCbase");
}

Int_t CheckLoadLibrary(const char* library)
{
  // checks if a library is already loaded, if not loads the library
  
  if (strlen(gSystem->GetLibraries(library, "", kFALSE)) > 0)
    return 1;
  
  return gSystem->Load(library);
}
