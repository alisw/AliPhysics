#if !defined(__CINT__) || defined(__MAKECINT__)

#include <stdio.h>
#include <stdlib.h>
#include <TROOT.h>
#include <Riostream.h>
#include <TClassTable.h>
#include <TStyle.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TLine.h>
#include <TGrid.h>
#include <TBits.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TBranch.h>
#include <TFileMerger.h>
#include <TGridResult.h>
#include <TSystem.h>
#include <TLegend.h>

#endif

void MakePlots(TString ntupleFileName){
  //***********************************************************
  // #### ROOT initialization
  gROOT->Reset();
  gStyle->SetCanvasColor(10);
  gStyle->SetFrameFillColor(10);
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetTitleTextColor(4);
  gStyle->SetStatTextColor(4);
  gStyle->SetStatX(0.92);
  gStyle->SetStatY(0.92);
  gStyle->SetLineColor(1);
  gStyle->SetPalette(1);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadRightMargin(0.04);
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadLeftMargin(0.06); 
  gStyle->SetTitleOffset(1.2,"X");
  gStyle->SetTitleOffset(0.7,"Y");
  // *************************************************************

  TFile* fil = new TFile(ntupleFileName.Data(),"read");
  if(!fil){
    printf("File with ntuple does not exist\n");
    return;
  }
  TNtuple* ntzdc = (TNtuple*)fil->Get("ntzdc");

  Float_t nrun;
  Float_t meanZNC,meanZPC,meanZNA,meanZPA,meanZEM1,meanZEM2;
  Float_t emeanZNC,emeanZPC,emeanZNA,emeanZPA,emeanZEM1,emeanZEM2;
  Float_t ZNClg,ZNAlg,pmcZNClg,pmcZNAlg;
  Float_t eZNClg,eZNAlg,epmcZNClg,epmcZNAlg;
  Float_t xZNC,yZNC,xZNA,yZNA,tdcSum,tdcDiff;
  
  ntzdc->SetBranchAddress("nrun",&nrun);
  ntzdc->SetBranchAddress("meanZNC",&meanZNC);
  ntzdc->SetBranchAddress("meanZPC",&meanZPC);
  ntzdc->SetBranchAddress("meanZNA",&meanZNA);
  ntzdc->SetBranchAddress("meanZPA",&meanZPA);
  ntzdc->SetBranchAddress("meanZEM1",&meanZEM1);
  ntzdc->SetBranchAddress("meanZEM2",&meanZEM2);
  ntzdc->SetBranchAddress("emeanZNC",&emeanZNC);
  ntzdc->SetBranchAddress("emeanZPC",&emeanZPC);
  ntzdc->SetBranchAddress("emeanZNA",&emeanZNA);
  ntzdc->SetBranchAddress("emeanZPA",&emeanZPA);
  ntzdc->SetBranchAddress("emeanZEM1",&emeanZEM1);
  ntzdc->SetBranchAddress("emeanZEM2",&emeanZEM2);
  ntzdc->SetBranchAddress("ZNClg",&ZNClg);
  ntzdc->SetBranchAddress("ZNAlg",&ZNAlg);
  ntzdc->SetBranchAddress("eZNClg",&eZNClg);
  ntzdc->SetBranchAddress("eZNAlg",&eZNAlg);
  ntzdc->SetBranchAddress("pmcZNClg",&pmcZNClg);
  ntzdc->SetBranchAddress("pmcZNAlg",&pmcZNAlg);
  ntzdc->SetBranchAddress("epmcZNClg",&epmcZNClg);
  ntzdc->SetBranchAddress("epmcZNAlg",&epmcZNAlg);
  ntzdc->SetBranchAddress("xZNC",&xZNC);
  ntzdc->SetBranchAddress("yZNC",&yZNC);
  ntzdc->SetBranchAddress("xZNA",&xZNA);
  ntzdc->SetBranchAddress("yZNA",&yZNA);
  ntzdc->SetBranchAddress("tdcSum",&tdcSum);
  ntzdc->SetBranchAddress("tdcDiff",&tdcDiff);

  TH1F *hznc      = new TH1F("hznc","",(Int_t)ntzdc->GetEntries(),0.,ntzdc->GetEntries());
  TH1F *hzna      = new TH1F("hzna","",(Int_t)ntzdc->GetEntries(),0.,ntzdc->GetEntries());
  TH1F *hzpc      = new TH1F("hzpc","",(Int_t)ntzdc->GetEntries(),0.,ntzdc->GetEntries());
  TH1F *hzpa      = new TH1F("hzpa","",(Int_t)ntzdc->GetEntries(),0.,ntzdc->GetEntries());
  TH1F *hzem1     = new TH1F("hzem1","",(Int_t)ntzdc->GetEntries(),0.,ntzdc->GetEntries());
  TH1F *hzem2     = new TH1F("hzem2","",(Int_t)ntzdc->GetEntries(),0.,ntzdc->GetEntries());
  TH1F *hznclg    = new TH1F("hznclg","",(Int_t)ntzdc->GetEntries(),0.,ntzdc->GetEntries());
  TH1F *hznalg    = new TH1F("hznalg","",(Int_t)ntzdc->GetEntries(),0.,ntzdc->GetEntries());
  TH1F *hpmcznclg = new TH1F("hpmcznclg","",(Int_t)ntzdc->GetEntries(),0.,ntzdc->GetEntries());
  TH1F *hpmcznalg = new TH1F("hpmcznalg","",(Int_t)ntzdc->GetEntries(),0.,ntzdc->GetEntries());
  TH1F *hxznc     = new TH1F("hxznc","",(Int_t)ntzdc->GetEntries(),0.,ntzdc->GetEntries());
  TH1F *hyznc     = new TH1F("hyznc","",(Int_t)ntzdc->GetEntries(),0.,ntzdc->GetEntries());
  TH1F *hxzna     = new TH1F("hxzna","",(Int_t)ntzdc->GetEntries(),0.,ntzdc->GetEntries());
  TH1F *hyzna     = new TH1F("hyzna","",(Int_t)ntzdc->GetEntries(),0.,ntzdc->GetEntries());
  TH1F *htdcsum   = new TH1F("htdcsum","",(Int_t)ntzdc->GetEntries(),0.,ntzdc->GetEntries());
  TH1F *htdcdiff  = new TH1F("htdcdiff","",(Int_t)ntzdc->GetEntries(),0.,ntzdc->GetEntries());

  for(Int_t i = 0; i<ntzdc->GetEntries();i++){
    ntzdc->GetEvent(i);
    //
    hznc->SetBinContent(i+1, meanZNC);
    hznc->SetBinError(i+1, emeanZNC);
    hznc->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hznc->GetXaxis()->SetTitle("RUN #");
    hznc->GetYaxis()->SetTitle("ZN mean signal");
    hzpc->SetBinContent(i+1, meanZPC);
    hzpc->SetBinError(i+1, emeanZPC);
    hzpc->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hzpc->GetXaxis()->SetTitle("RUN #");
    hzpc->GetYaxis()->SetTitle("ZP mean signal");
    hzna->SetBinContent(i+1, meanZNA);
    hzna->SetBinError(i+1, emeanZNA);
    hzna->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hzna->GetXaxis()->SetTitle("RUN #");
    hzna->GetYaxis()->SetTitle("ZN mean signal");
    hzpa->SetBinContent(i+1, meanZPA);
    hzpa->SetBinError(i+1, emeanZPA);
    hzpa->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hzpa->GetXaxis()->SetTitle("RUN #");
    hzpa->GetYaxis()->SetTitle("ZP mean signal");
    hzem1->SetBinContent(i+1, meanZEM1);
    hzem1->SetBinError(i+1, emeanZEM1);
    hzem1->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hzem1->GetXaxis()->SetTitle("RUN #");
    hzem1->GetYaxis()->SetTitle("ZEM mean signal");
    hzem2->SetBinContent(i+1, meanZEM2);
    hzem2->SetBinError(i+1, emeanZEM2);
    hzem2->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hzem2->GetXaxis()->SetTitle("RUN #");
    hzem2->GetYaxis()->SetTitle("ZEM mean signal");
    hznclg->SetBinContent(i+1, ZNClg);
    hznclg->SetBinError(i+1, eZNClg);
    hznclg->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hznclg->GetXaxis()->SetTitle("RUN #");
    hznclg->GetYaxis()->SetTitle("ZN LG mean signal");
    hznalg->SetBinContent(i+1, ZNAlg);
    hznalg->SetBinError(i+1, eZNAlg);
    hznalg->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hznalg->GetXaxis()->SetTitle("RUN #");
    hznalg->GetYaxis()->SetTitle("ZN LG mean signal");
    hpmcznclg->SetBinContent(i+1, pmcZNClg);
    hpmcznclg->SetBinError(i+1, epmcZNClg);
    hpmcznclg->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hpmcznclg->GetXaxis()->SetTitle("RUN #");
    hpmcznclg->GetYaxis()->SetTitle("ZN LG PMC mean signal");
    hpmcznalg->SetBinContent(i+1, pmcZNAlg);
    hpmcznalg->SetBinError(i+1, epmcZNAlg);
    hpmcznalg->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hpmcznalg->GetXaxis()->SetTitle("RUN #");
    hpmcznalg->GetYaxis()->SetTitle("ZN LG PMC mean signal");
    hxznc->SetBinContent(i+1, xZNC);
    hxznc->SetBinError(i+1, 0.1*xZNC);
    hxznc->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hxznc->GetXaxis()->SetTitle("RUN #");
    hxznc->GetYaxis()->SetTitle("X_{ZN} (cm)");
    hyznc->SetBinContent(i+1, yZNC);
    hyznc->SetBinError(i+1, 0.1*yZNC);
    hyznc->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hyznc->GetXaxis()->SetTitle("RUN #");
    hyznc->GetYaxis()->SetTitle("Y_{ZN} (cm)");
    hxzna->SetBinContent(i+1, xZNA);
    hxzna->SetBinError(i+1, 0.1*xZNA);
    hxzna->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hxzna->GetXaxis()->SetTitle("RUN #");
    hxzna->GetYaxis()->SetTitle("X_{ZN} (cm)");
    hyzna->SetBinContent(i+1, yZNA);
    hyzna->SetBinError(i+1, 0.1*yZNA);
    hyzna->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    hyzna->GetXaxis()->SetTitle("RUN #");
    hyzna->GetYaxis()->SetTitle("Y_{ZN} (cm)");
    htdcsum->SetBinContent(i+1, tdcSum);
    htdcsum->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    htdcsum->GetXaxis()->SetTitle("RUN #");
    htdcsum->GetYaxis()->SetTitle("TDC Sum (ns)");
    htdcdiff->SetBinContent(i+1, tdcDiff);
    htdcdiff->GetXaxis()->SetBinLabel(i+1,Form("%d",(Int_t)nrun));
    htdcdiff->GetXaxis()->SetTitle("RUN #");
    htdcdiff->GetYaxis()->SetTitle("TDC Diff (ns)");
  }


  TCanvas *c1 = new TCanvas("c1", "Mean value ZNs", 0, 0, 1200, 1000);
  c1->Divide(1,2);
  //
  c1->cd(1); 
  hznc->SetMarkerColor(kAzure+6); hznc->SetLineColor(kAzure+6); 
  hznc->SetMarkerStyle(21);
  hznc->Draw(""); 
  hzna->SetMarkerColor(kPink-2);  hzna->SetLineColor(kPink-2);
  hzna->SetMarkerStyle(20);
  hzna->Draw("SAME"); 
  //
  TLegend *l1 = new TLegend(0.44,0.18,0.54,0.32);
  l1->SetFillColor(kWhite);
  l1->AddEntry(hznc," ZNC " ,"P");
  l1->AddEntry(hzna," ZNA " ,"P");
  l1->Draw("");
  //
  c1->cd(2);
  hznclg->SetMarkerColor(kPink+5); hznclg->SetLineColor(kPink+5); 
  hznclg->SetMarkerStyle(21);
  hznclg->Draw(""); 
  hznalg->SetMarkerColor(kBlue+1); hznalg->SetLineColor(kBlue+1); 
  hznalg->SetMarkerStyle(20);  hznalg->SetMinimum(0);
  hznalg->Draw("SAME");
  //
  TLegend *l2 = new TLegend(0.44,0.18,0.54,0.32);
  l2->SetFillColor(kWhite);
  l2->AddEntry(hznclg," ZNC LG " ,"P");
  l2->AddEntry(hznalg," ZNA LG " ,"P");
  l2->Draw("");
  

  TCanvas *c1b = new TCanvas("c1b", "Mean value ZEMs ZPs", 200, 0, 1200, 1000);
  c1b->Divide(1,2);
  c1b->cd(1); 
  hzpc->SetMarkerColor(kAzure+6); hzpc->SetLineColor(kAzure+6); 
  hzpc->SetMarkerStyle(21);
  hzpc->Draw(""); 
  hzpa->SetMarkerColor(kPink-2);  hzpa->SetLineColor(kPink-2);
  hzpa->SetMarkerStyle(20);
  hzpa->Draw("SAME"); 
  //
  TLegend *l3 = new TLegend(0.44,0.18,0.54,0.32);
  l3->SetFillColor(kWhite);
  l3->AddEntry(hzpc," ZPC " ,"P");
  l3->AddEntry(hzpa," ZPA " ,"P");
  l3->Draw("");
  
  c1b->cd(2);
  hzem1->SetMarkerColor(kTeal-7); hzem1->SetLineColor(kTeal-7);
  hzem1->SetMarkerStyle(29);
  hzem2->SetMarkerColor(kTeal+5); hzem2->SetLineColor(kTeal+5);
  hzem2->SetMarkerStyle(30);
  hzem2->Draw(""); 
  hzem1->Draw("SAME"); 
  //
  TLegend *l4 = new TLegend(0.44,0.18,0.54,0.32);
  l4->SetFillColor(kWhite);
  l4->AddEntry(hzem1," ZEM1 " ,"P");
  l4->AddEntry(hzem2," ZEM2 " ,"P");
  l4->Draw("");
 
  
  /*TCanvas *c2 = new TCanvas("c2", "ZN centroids", 400, 400, 1200, 1000);
  c2->Divide(1,2);
  //
  c2->cd(1);
  hxznc->SetMarkerColor(kAzure+5); hxznc->SetLineColor(kAzure+5);
  hxznc->SetMarkerStyle(21); 
  hxzna->SetMarkerColor(kPink+5); hxzna->SetLineColor(kPink+5);
  hxzna->SetMarkerStyle(20); 
  hxznc->Draw(""); 
  hxzna->Draw("SAME");
  c2->cd(2);
  hyznc->SetMarkerColor(kAzure+5); hyznc->SetLineColor(kAzure+5);
  hyznc->SetMarkerStyle(21); 
  hyzna->SetMarkerColor(kPink+5); hyzna->SetLineColor(kPink+5);
  hyzna->SetMarkerStyle(20); 
  hyznc->Draw(""); 
  hyzna->Draw("SAME");*/

  
  TCanvas *c3 = new TCanvas("c3", "LG signals", 600, 0, 1200, 1000);
  c3->Divide(1,2);
  //
  c3->cd(1);
  hznclg->SetMarkerColor(kPink+5); hznclg->SetLineColor(kPink+5); 
  hznclg->SetMarkerStyle(21);
  hznclg->Draw(""); 
  hznalg->SetMarkerColor(kBlue+1); hznalg->SetLineColor(kBlue+1); 
  hznalg->SetMarkerStyle(20);
  hznalg->Draw("SAME"); 
  TLegend *l7 = new TLegend(0.44,0.18,0.54,0.32);
  l7->SetFillColor(kWhite);
  l7->AddEntry(hznclg," ZNC LG " ,"P");
  l7->AddEntry(hznalg," ZNA LG " ,"P");
  l7->Draw("");
  //
  c3->cd(2);
  hpmcznclg->SetMarkerColor(kPink+6); hpmcznclg->SetLineColor(kPink+6); 
  hpmcznclg->SetMarkerStyle(21); hpmcznclg->SetMinimum(0);
  hpmcznclg->Draw(""); 
  hpmcznalg->SetMarkerColor(kAzure+7); hpmcznalg->SetLineColor(kAzure+7); 
  hpmcznalg->SetMarkerStyle(20); hpmcznalg->SetMinimum(0);
  hpmcznalg->Draw("SAME"); 
  TLegend *l8 = new TLegend(0.44,0.15,0.58,0.32);
  l8->SetFillColor(kWhite);
  l8->AddEntry(hpmcznclg," ZNC LG PMC" ,"P");
  l8->AddEntry(hpmcznalg," ZNA LG PMC" ,"P");
  l8->Draw("");
}
