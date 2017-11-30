#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TKey.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TList.h"
#include "TMath.h"
#include "TTree.h"
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
void PublishCanvas(TList *qaList, const char* det, const char* name, TString nadd="", TString outputformat="");
void CheckPIDInTracking(TList* qaListTPC,const char* listname);
void SetupPadStyle();
void StoreTrendingVars(TH1* hMean, TH1* hSigma);


void LoadLibs();
Int_t CheckLoadLibrary(const char* library);

TCanvas *fCanvas=0x0;
TTree *fTree=0x0;

const Int_t nDetsForTrending=4; // ITS, TPC, TPC_TOF, TOF
TString nameDetsForTrending[nDetsForTrending]={"ITS","TPC_Basic","TPC_TOF","TOF"};
const Int_t nSpecForTrending=5; // e, pi, K ,p, d
TString nameSpecForTrending[nSpecForTrending]={"electron","pion","kaon","proton","deuteron"};
const Int_t nVarsForTrending=2; // mean, sigma of nsigma
TString nameVarsForTrending[nVarsForTrending]={"meannSigma","signSigma"};
const Int_t nMomsForTrending=2; // low pt, high pt
Float_t momForTrending[nDetsForTrending][2]={0.3,1.2, //ITS
					     0.3,1.5, //TPC
					     0.3,1.5, //TPC_TOF
					     1.,2.}; //TOF

const Int_t totTrending=nSpecForTrending*nDetsForTrending*nMomsForTrending*nVarsForTrending;
Float_t vecForTrend[totTrending];



void MakeTrendingPIDQA(const char* inputFile, TString dirInFile = "PIDqa", TString detList="all", Int_t runNumber=111111, const TString outputformat="png,root",const char* outputFile="PIDqaReport.pdf")
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
  
  printf("%s\n", listName.Data());
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
  detList.ToUpper();

  //
  // Make QA info
  //
  TFile * origf=0x0;
  if(detList.Contains("ALL")){
    fTree=new TTree("trending","tree of trending variables");
    fTree->Branch("nrun",&runNumber,"nrun/I");
    fTree->Fill();
  }else{
    origf=new TFile("trending.root","update");
    if(origf){
      TList* l=(TList*)origf->GetListOfKeys();
      for(Int_t j=0; j<l->GetEntries(); j++){
	TKey* o=(TKey*)l->At(j);
	TString cname=o->GetClassName();
	if(cname.Contains("TTree")){
	  fTree=(TTree*)origf->Get(o->GetName());
	  break;
	}
      }
    }
  }
  

  // ITS PID
  if(detList.Contains("ALL") || detList.Contains("ITS")) PublishCanvas(qaList,"ITS","hNsigmaP_ITS_%s","",outputformat);

  // TPC PID
  if(detList.Contains("ALL") || detList.Contains("TPC")){
    TList *qaListTPC = (TList*)qaList->FindObject("TPC");
    if (qaListTPC){
      CheckPIDInTracking(qaListTPC,"TPCBasic");
      PublishCanvas(qaListTPC,"TPCBasic","hNsigmaP_TPC_Basic_%s","",outputformat);
      PublishCanvas(qaListTPC,"TPCV0","hNsigmaP_TPC_V0_%s","",outputformat);
      //   if (man->GetCurrentPeriod()=="11h"){
      //     PublishCanvas(qaListTPC,"TPC","hNsigmaP_TPC_Basic_%s_Hybrid","Hybrid");
      //     PublishCanvas(qaListTPC,"TPC","hNsigmaP_TPC_Basic_%s_OROChigh","OROChigh");
      //   }
    }
    else {
      printf("Could not find list '%s/TPC' in file '%s'\n", listName.Data(), f.GetName());
    }
  
    // TPC Response info
    TObjArray *qaInfo=(TObjArray*)qaList->FindObject("QAinfo");
    TObjArray *tpcInfo=0x0;
    if (qaInfo && (tpcInfo=(TObjArray*)qaInfo->FindObject("TPC_info"))){
      TH1F* tpcSplineInfo=(TH1F*)tpcInfo->FindObject("TPCsplineNames");
      TH1F* tpcConfigInfo=(TH1F*)tpcInfo->FindObject("TPCconfigInfo");
      fCanvas->Divide(1,2);
      
      TPaveText pt(.1,.1,.9,.9,"NDC");
      pt.SetBorderSize(1);
      pt.SetFillColor(0);
      pt.SetTextSizePixels(16);
      
      if (tpcSplineInfo){
	for (Int_t i=1; i<=tpcSplineInfo->GetNbinsX();++i) pt.AddText(tpcSplineInfo->GetXaxis()->GetBinLabel(i));
      }
      
      TPaveText pt2(.1,.1,.9,.9,"NDC");
      pt2.SetBorderSize(1);
      pt2.SetFillColor(0);
      pt2.SetTextSizePixels(16);
      if (tpcConfigInfo){
       	for (Int_t i=1; i<=tpcConfigInfo->GetNbinsX();++i) pt2.AddText(tpcConfigInfo->GetXaxis()->GetBinLabel(i));
      }
      
      fCanvas->cd(1);
      pt.Draw();
      fCanvas->cd(2);
      pt2.Draw();
      if(outputformat.Contains("png")) fCanvas->SaveAs("TPCSplinesAndConfigInfo.png");
      if(outputformat.Contains("eps")) fCanvas->SaveAs("TPCSplinesAndConfigInfo.eps");
      fCanvas->Update();
      fCanvas->Clear();
    }

    // TPC PID after 3 sigma TOF cut
    PublishCanvas(qaList,"TPC_TOF","hNsigmaP_TPC_TOF_%s","",outputformat);

  }
  
  // TOF PID
  if(detList.Contains("ALL") || detList.Contains("TOF")) PublishCanvas(qaList,"TOF","hNsigmaP_TOF_%s","",outputformat);

  // TRD PID
  if(detList.Contains("ALL") || detList.Contains("TRD")){
    fCanvas->Divide(2,3);
    TH2 *hLikeP_TRD_3tls_electron=Get2DHistogramfromList(qaList,"TRD","hLikeP_TRD_3tls_electron");
    TH2 *hLikeP_TRD_3tls_pion=Get2DHistogramfromList(qaList,"TRD","hLikeP_TRD_3tls_pion");
    TH2 *hLikeP_TRD_4tls_electron=Get2DHistogramfromList(qaList,"TRD","hLikeP_TRD_4tls_electron");
    TH2 *hLikeP_TRD_4tls_pion=Get2DHistogramfromList(qaList,"TRD","hLikeP_TRD_4tls_pion");
    TH2 *hLikeP_TRD_5tls_electron=Get2DHistogramfromList(qaList,"TRD","hLikeP_TRD_5tls_electron");
    TH2 *hLikeP_TRD_5tls_pion=Get2DHistogramfromList(qaList,"TRD","hLikeP_TRD_5tls_pion");
    
    /*
     *  cTRDnsigma[countcanvas]->cd(1);
     *  TPaveText pt3TRD(.02,.02,.49,.52);
     *  pt3TRD.SetTextAlign(11);
     *  pt3TRD.SetTextSizePixels(16);
     *  pt3TRD.AddText(Form(" TRD PID QA %s.%s.%d", first.Data(), man->GetCurrentPeriod().Data(), pass));
     *  pt3TRD.Draw();
     */
    fCanvas->cd(1);
    SetupPadStyle();
    if(hLikeP_TRD_3tls_electron) hLikeP_TRD_3tls_electron->Draw("colz");
    fCanvas->cd(2);
    SetupPadStyle();
    if(hLikeP_TRD_3tls_pion) hLikeP_TRD_3tls_pion->Draw("colz");
    fCanvas->cd(3);
    SetupPadStyle();
    if(hLikeP_TRD_4tls_electron) hLikeP_TRD_4tls_electron->Draw("colz");
    fCanvas->cd(4);
    SetupPadStyle();
    if(hLikeP_TRD_4tls_pion) hLikeP_TRD_4tls_pion->Draw("colz");
    fCanvas->cd(5);
    SetupPadStyle();
    if(hLikeP_TRD_5tls_electron) hLikeP_TRD_5tls_electron->Draw("colz");
    fCanvas->cd(6);
    SetupPadStyle();
    if(hLikeP_TRD_5tls_pion) hLikeP_TRD_5tls_pion->Draw("colz");
    
    fCanvas->Update();
    fCanvas->Clear();
  }
  
  delete qaList;

  p.Close();
  delete fCanvas;

  if(detList.Contains("ALL")){
    if(fTree){
      TFile* fouttree=new TFile("trending.root","recreate");
      fTree->Write();
      fouttree->Close();
      delete fouttree;
    }
  }else{
    if(origf && fTree){
      origf->cd();
      fTree->Write("", TObject::kOverwrite);
      origf->Close();
      delete origf;
    }
  }
 
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
  StoreTrendingVars(hMean,hSigma);
}

void PublishCanvas(TList *qaList, const char* det, const char* name, TString nadd, TString outputformat)
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
  if(outputformat.Contains("png")) fCanvas->SaveAs(Form("%snSigmaPID.png",det));
  if(outputformat.Contains("eps"))   fCanvas->SaveAs(Form("%snSigmaPID.eps",det));
  if(outputformat.Contains("root"))   fCanvas->SaveAs(Form("%snSigmaPID.root",det));
  fCanvas->Update();
  fCanvas->Clear();

}

void CheckPIDInTracking(TList* qaListTPC,const char* listname){

  TList *histolist = (TList *)qaListTPC->FindObject(listname);
  if (!histolist) {printf(" list not found \n");  return; }
  TCanvas* ctrpid=new TCanvas("ctrpid","PID-in-track",1000,700);
  TCanvas* ctrpidall=new TCanvas("ctrpidall","PID-in-track",1500,700);
  TLegend* leg=new TLegend(0.7,0.35,0.89,0.87);
  leg->SetHeader("PID in tracking");
  ctrpidall->Divide(3,3);
  Bool_t saveCan=kFALSE;
  Int_t cols[9]={kGreen+2,kGray,1,2,4,kMagenta,kOrange+1,kYellow,kCyan};
  for(Int_t j=0; j<9; j++){
    TString histoname=Form("hSigP_TPC_TrackedAs_%s",AliPID::ParticleName(j));
    TH2* histo = (TH2*)histolist->FindObject(histoname.Data());
    if(histo){
      saveCan=kTRUE;
      TH2* histo2=(TH2*)histo->Clone(Form("%scolor",histoname.Data()));
      histo2->SetTitle(" ");
      histo2->SetStats(0);
      histo2->SetMarkerColor(cols[j]);
      leg->AddEntry(histo2,Form("%s",AliPID::ParticleName(j)),"")->SetTextColor(histo2->GetMarkerColor());
      ctrpid->cd();
      SetupPadStyle();
      if(j==0) histo2->Draw();
      else histo2->Draw("same");
      ctrpidall->cd(j+1);
      SetupPadStyle();
      histo->Draw("colz");
    }
  }
  if(saveCan){
    ctrpid->cd();
    leg->Draw();
    ctrpidall->SaveAs("TPCdEdx-PIDinTracking-9pads.png");
    ctrpid->SaveAs("TPCdEdx-PIDinTracking.png");
  }
}


void StoreTrendingVars(TH1* hMean, TH1* hSigma){

  // interface method to store trending variables
  
  if(!fTree) return;
  TString hname=hMean->GetName();
  cout<<hname.Data()<<endl;
  Int_t theDet=-1;
  for(Int_t iDet=0; iDet<nDetsForTrending; iDet++){
    if(hname.Contains(Form("hNsigmaP_%s",nameDetsForTrending[iDet].Data()))) theDet=iDet;
  }
  Int_t theSpec=-1;
  for(Int_t iSpec=0; iSpec<nSpecForTrending; iSpec++){
    if(hname.Contains(nameSpecForTrending[iSpec].Data())) theSpec=iSpec;
  }
  if(theDet>=0 && theSpec>=0){
    for(Int_t iMom=0; iMom<nMomsForTrending; iMom++){
      
    // Int_t firstIndex=theDet*(nSpecForTrending*nMomsForTrending*nVarsForTrending)+theSpec*(nMomsForTrending*nVarsForTrending);
    // Int_t binLow=hMean->FindBin(momForTrending[theDet][0]);
    // Int_t binHigh=hMean->FindBin(momForTrending[theDet][1]);
    // vecForTrend[firstIndex]=hMean->GetBinContent(binLow);
    // vecForTrend[++firstIndex]=hMean->GetBinContent(binHigh);
    // vecForTrend[++firstIndex]=hSigma->GetBinContent(binLow);
    // vecForTrend[++firstIndex]=hSigma->GetBinContent(binHigh);
      Int_t hbin=hMean->FindBin(momForTrending[theDet][iMom]);
      Float_t mu=hMean->GetBinContent(hbin);
      Float_t sig=hSigma->GetBinContent(hbin);
      TString bnam1=Form("meannSigma%s_%s_p%dMeV",nameDetsForTrending[theDet].Data(),nameSpecForTrending[theSpec].Data(),TMath::Nint(momForTrending[theDet][iMom]*1000.));
      TString bnam2=Form("signSigma%s_%s_p%dMeV",nameDetsForTrending[theDet].Data(),nameSpecForTrending[theSpec].Data(),TMath::Nint(momForTrending[theDet][iMom]*1000.));
      TBranch *br1 = fTree->Branch(bnam1.Data(),&mu,Form("%s/F",bnam1.Data()));
      TBranch *br2 = fTree->Branch(bnam2.Data(),&sig,Form("%s/F",bnam2.Data()));
      br1->Fill();
      br2->Fill();
    }
  }
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


