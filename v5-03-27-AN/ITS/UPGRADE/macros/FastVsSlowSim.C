#ifndef __CINT__
#include <TTree.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TStyle.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliRawReader.h"
#include <AliTPCRawStreamV3.h>
#include "AliFMDRawReader.h"
#include "AliFMDParameters.h"
#include "AliFMDDigit.h"
#include "AliTPCAltroMapping.h"
#include "AliTriggerConfiguration.h"
#include "TAlienCollection.h"
#include "AliTriggerClass.h"
#include "TGrid.h"
#include "TPRegexp.h"
#include "AliVZERORawStream.h"

#include "AliITSsegmentationUpgrade.h"
#include "AliRunLoader.h"
#include "AliStack.h"
#include "AliITSRecPointU.h"
#include "AliITSDigitUpgrade.h"
#include "TParticle.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliTrackReference.h"
#include <TROOT.h>
#include <TStyle.h>
#include "TGraphErrors.h"
#include <TF1.h>
#include <TGeoManager.h>

#endif
#include <iostream>
#include <iomanip>
#include <cstdio>

/*
 .L ~/ITSupgrade/CDRpics/FastVsSlowSim.C
 ExtractOutputHistos(0,1);


 .L ~/ITSupgrade/BuildDetector/DetectorK.cxx++
 .L FastVsSlowSim.C

 FastVsSlowSimRes();
 FastVsSlowSimPtRes();
 FastVsSlowSimEff(0,1);
 FastVsSlowSimEff(1,1);
 FastVsSlowSimEff(2,1);
 FastVsSlowSimEff(3,1);

 FastVsSlowSimEff(0);

*/


Int_t atLeastcorr = 3;


void FastVsSlowSim(Bool_t extract=1) {

  if (extract)
    ExtractOutputHistos(0,1);
  else
    plotMerged();
}

TGraph *gr[5];
Int_t colors[5]={1,2,3,4,6};
Int_t width =2;


// new ideal Pixel properties?

Double_t etaCut = 0.9;
Double_t X0     = 0.003;
Double_t resRPhi = 0.0004;
Double_t resZ   = 0.0004;

void FastVsSlowSimRes() {

  Int_t plusTPC =0;

  gROOT->LoadMacro("~/fig_template.C"); // figure style
  myOptions(0);
  gROOT->ForceStyle();

  TCanvas *myCan = new TCanvas("myCan");
  myCan->Draw();
  myCan->cd();
  
  TPad *myPad = new TPad("myPad", "The pad",0,0,1,1);
  myPadSetUp(myPad,0.15,0.04,0.04,0.15);
  myPad->Draw();   myPad->cd();
  myPad->SetGridx();   myPad->SetGridy();  myPad->SetLogx(); 

  //  TLegend *leg = new TLegend(0.7,160,20,290,"","brCDN"); 
  TLegend *leg = new TLegend(0.44,160,1.7,290,"","brCDN"); 
 
  leg->SetFillColor(0);

  
  // Current ITS +++++++++++++++++++++++++++++++++++++++++

  DetectorK its("ALICE","ITS");
  its.MakeAliceCurrent(0,plusTPC);
  its.SetMaxRadiusOfSlowDetectors(0.1);
  its.SolveViaBilloir(0);
  Int_t color=1; Int_t linewidth=2;

  TGraph *c[6];
  TGraph *d[6];

  Int_t pi =0;
  d[pi] = its.GetGraphPointingResolution(1,color,linewidth);
  d[pi]->SetLineStyle(2);
  //  d[pi]->GetYaxis()->SetTitle("Pointing resolution #sigma [#mum]");
  //  d[pi]->SetTitle("Pointing resolution .vs. Pt");
  //  d[pi]->Draw("AC");
  
  c[pi] = its.GetGraphPointingResolution(0,color,linewidth);
  c[pi]->SetMinimum(-1);
  c[pi]->Draw("AC");

  leg->AddEntry(c[pi],"FastTool:  Current ITS","l");
  //  leg->AddEntry(d[pi],"in z  - Current ITS","l");

 

 
  // Current ITS +++++++++++++++++++++++++++++++++++++++++

  Int_t color=3; Int_t linewidth=2;
  Int_t pi =2;

  DetectorK its("ALICE","ITS");
  its.MakeAliceCurrent(0,plusTPC);
  
  its.SetRadius("bpipe",2.0);
  its.AddLayer("spd0", 2.2,1,1,1);  

  its.SetRadius("spd0",2.2); its.SetRadiationLength("spd0",X0); its.SetResolution("spd0",resRPhi,resZ);
  its.SetRadius("spd1",4.8);   its.SetRadiationLength("spd1",X0); its.SetResolution("spd1",resRPhi,resZ);
  its.SetRadius("spd2",9.1);   its.SetRadiationLength("spd2",X0); its.SetResolution("spd2",resRPhi,resZ);

  its.SetMaxRadiusOfSlowDetectors(0.1);
  its.SolveViaBilloir(0);

  d[pi] = its.GetGraphPointingResolution(1,color,linewidth);
  d[pi]->SetLineStyle(2);
  //  d[pi]->Draw("C");

  c[pi] = its.GetGraphPointingResolution(0,color,linewidth);
  c[pi]->Draw("C");

  leg->AddEntry(c[pi],"FastTool: \"New SPDs\"","l");
  //  leg->AddEntry(d[pi],"in z  - \"New SPDs\"","l");



  // ALL NEW +++++++++++++++++++++++++++++++++++++++++++

  color=2; Int_t linewidth=2;
  Int_t pi =1; 


  // for a 0.8,0.2 weight configuration
  
  DetectorK *itsU = new DetectorK((char*)"ALICE",(char*)"ITS");
  
  itsU->AddLayer((char*)"bpipe", 2.0,0.0022); // beam pipe
  itsU->AddLayer((char*)"vertex",  0,     0); // dummy vertex for matrix calculation
  
  itsU->AddLayer("ddd1",  2.2 ,  X0, resRPhi, resZ); 
  itsU->AddLayer("ddd2",  3.8 ,  X0, resRPhi, resZ); 
  itsU->AddLayer("ddd3",  6.8 ,  X0, resRPhi, resZ); 
  itsU->AddLayer("ddd4", 12.4 ,  X0, resRPhi, resZ); 
  itsU->AddLayer("ddd5", 23.5 ,  X0, resRPhi, resZ); 
  itsU->AddLayer("ddd6", 39.6 ,  X0, resRPhi, resZ); 
  itsU->AddLayer("ddd7", 43.0 ,  X0, resRPhi, resZ); 
 
  if(plusTPC) itsU->AddTPC(0.1,0.1);
  itsU->SetMaxRadiusOfSlowDetectors(0.1);
  itsU->SolveViaBilloir(0);
  itsU->PrintLayout();

  
  d[pi] = itsU->GetGraphPointingResolution(1,color,linewidth);
  d[pi]->SetLineStyle(2);
  //  d[pi]->Draw("C");

  c[pi] = itsU->GetGraphPointingResolution(0,color,linewidth);
  c[pi]->SetMaximum(150);
  c[pi]->Draw("C");

  leg->AddEntry(c[pi],"FastTool: \"All New\" ","l");
  //  leg->AddEntry(d[pi],"in z  - \"All New\" ","l");


  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


 
  TFile f1("root/FastVsSlow_CurrentITS-PbPb-fran.root");
  TFile f2("root/FastVsSlow_NewSPDs-PbPb-fran.root");
  TFile f3("root/FastVsSlow_AllNew-PbPb-fran.root");
  TGraphErrors *dca1 = (TGraphErrors*)f1.Get("dca");
  TGraphErrors *dca2 = (TGraphErrors*)f2.Get("dca");
  TGraphErrors *dca3 = (TGraphErrors*)f3.Get("dca");
  
  dca1->SetMarkerStyle(21); dca1->SetMarkerColor(1);
  dca2->SetMarkerStyle(21); dca2->SetMarkerColor(3);
  dca3->SetMarkerStyle(21); dca3->SetMarkerColor(2);

  leg->AddEntry(dca1,"FullMC: Current ITS","PE");
  leg->AddEntry(dca2,"FullMC: \"New SPDs\"","PE");
  leg->AddEntry(dca3,"FullMC: \"All New\" ","PE");

  dca1->Draw("APE"); dca1->SetMinimum(-1); dca1->SetMaximum(300);
  dca2->Draw("PE");
  dca3->Draw("PE");
  c[0]->Draw("C");
  c[1]->Draw("C");
  c[2]->Draw("C");

  leg->Draw();

  myCan->SaveAs(Form("FastVsSlowSim-Res-%d.pdf",plusTPC));
  myCan->SaveAs(Form("FastVsSlowSim-Res-%d.eps",plusTPC));


}

// ==============================================================================================
// ==============================================================================================

void FastVsSlowSimPtRes() {

  Int_t plusTPC =0;

  gROOT->LoadMacro("~/fig_template.C"); // figure style
  myOptions(0);
  gROOT->ForceStyle();

  TCanvas *myCan = new TCanvas("myCan");
  myCan->Draw();
  myCan->cd();
  
  TPad *myPad = new TPad("myPad", "The pad",0,0,1,1);
  myPadSetUp(myPad,0.15,0.04,0.04,0.15);
  myPad->Draw();   myPad->cd();
  myPad->SetGridx();   myPad->SetGridy();  myPad->SetLogx(); myPad->SetLogy();


  TLegend *leg = new TLegend(0.44,0.13,0.1.7,0.9,"","brCDN"); 
  leg->SetFillColor(0);

  TGraph *c[6];
  
  // Current ITS +++++++++++++++++++++++++++++++++++++++++
  Int_t color=1; Int_t linewidth=2;

  Int_t pi =0;
 
  DetectorK its("ALICE","ITS");
  its.MakeAliceCurrent(0,plusTPC);
  its.SetMaxRadiusOfSlowDetectors(0.1);
  its.SolveViaBilloir(0);
  Int_t color=1; Int_t linewidth=2;

  c[pi] = its.GetGraphMomentumResolution(color,linewidth);
  c[pi]->Draw("AC");

  leg->AddEntry(c[pi],"FastTool: Current ITS","l");


  // Current ITS +++++++++++++++++++++++++++++++++++++++++

  Int_t color=3; Int_t linewidth=2;
  Int_t pi =2;

  DetectorK its("ALICE","ITS");
  its.MakeAliceCurrent(0,plusTPC);
  
  its.SetRadius("bpipe",2.0);
  its.AddLayer("spd0", 2.2,1,1,1);  

  its.SetRadius("spd0",2.2); its.SetRadiationLength("spd0",X0); its.SetResolution("spd0",resRPhi,resZ);
  its.SetRadius("spd1",4.8);   its.SetRadiationLength("spd1",X0); its.SetResolution("spd1",resRPhi,resZ);
  its.SetRadius("spd2",9.1);   its.SetRadiationLength("spd2",X0); its.SetResolution("spd2",resRPhi,resZ);

  its.SetMaxRadiusOfSlowDetectors(0.1);
  its.SolveViaBilloir(0);

  c[pi] = its.GetGraphMomentumResolution(color,linewidth);
  c[pi]->Draw("C");

  leg->AddEntry(c[pi],"FastTool: \"New SPDs\"","l");


  // ALL NEW +++++++++++++++++++++++++++++++++++++++++++

  color=2; Int_t linewidth=2;
  Int_t pi =1; 


  // for a 0.8,0.2 weight configuration
  
  DetectorK *itsU = new DetectorK((char*)"ALICE",(char*)"ITS");
  
  itsU->AddLayer((char*)"bpipe", 2.0,0.0022); // beam pipe
  itsU->AddLayer((char*)"vertex",  0,     0); // dummy vertex for matrix calculation
  
  itsU->AddLayer("ddd1",  2.2 ,  X0, resRPhi, resZ); 
  itsU->AddLayer("ddd2",  3.8 ,  X0, resRPhi, resZ); 
  itsU->AddLayer("ddd3",  6.8 ,  X0, resRPhi, resZ); 
  itsU->AddLayer("ddd4", 12.4 ,  X0, resRPhi, resZ); 
  itsU->AddLayer("ddd5", 23.5 ,  X0, resRPhi, resZ); 
  itsU->AddLayer("ddd6", 39.6 ,  X0, resRPhi, resZ); 
  itsU->AddLayer("ddd7", 43.0 ,  X0, resRPhi, resZ); 
 
  if(plusTPC) itsU->AddTPC(0.1,0.1);
  itsU->SetMaxRadiusOfSlowDetectors(0.1);
  itsU->SolveViaBilloir(0);
  itsU->PrintLayout();

  
  c[pi] = itsU->GetGraphMomentumResolution(color,linewidth);
  c[pi]->Draw("C");

  leg->AddEntry(c[pi],"FastTool: \"All New\" ","l");

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



  TFile f1("root/FastVsSlow_CurrentITS-PbPb-fran.root");
  TFile f2("root/FastVsSlow_NewSPDs-PbPb-fran.root");
  TFile f3("root/FastVsSlow_AllNew-PbPb-fran.root");

  TGraphErrors *gpt1 = (TGraphErrors*)f1.Get("dPt");
  TGraphErrors *gpt2 = (TGraphErrors*)f2.Get("dPt");
  TGraphErrors *gpt3 = (TGraphErrors*)f3.Get("dPt");
  
  gpt1->SetMarkerStyle(21); gpt1->SetMarkerColor(1);
  gpt2->SetMarkerStyle(21); gpt2->SetMarkerColor(3);
  gpt3->SetMarkerStyle(21); gpt3->SetMarkerColor(2);

  leg->AddEntry(gpt1,"FullMC: Current ITS","PE");
  leg->AddEntry(gpt2,"FullMC: \"New SPDs\"","PE");
  leg->AddEntry(gpt3,"FullMC: \"All New\" ","PE");

  gpt1->Draw("APE"); gpt1->SetMinimum(0.1); gpt1->SetMaximum(20);
  gpt2->Draw("PE");
  gpt3->Draw("PE");
  c[0]->Draw("C");
  c[1]->Draw("C");
  c[2]->Draw("C");

  leg->Draw();
 
  myCan->SaveAs(Form("FastVsSlowSim-PtRes-%d.pdf",plusTPC));
  myCan->SaveAs(Form("FastVsSlowSim-PtRes-%d.eps",plusTPC));



}



// ==============================================================================================
// ==============================================================================================

void FastVsSlowSimEff(Int_t id=0,Int_t PbPb=0) {

  Int_t mult = 2400; // 2800  // deducted from "Frackable"
  if (PbPb) mult=2800;


  Int_t plusTPC =0;

  gROOT->LoadMacro("~/fig_template.C"); // figure style
  myOptions(0);
  gROOT->ForceStyle();

  TCanvas *myCan = new TCanvas("myCan");
  myCan->Draw();
  myCan->cd();
  
  TPad *myPad = new TPad("myPad", "The pad",0,0,1,1);
  myPadSetUp(myPad,0.15,0.04,0.04,0.15);
  myPad->Draw();   myPad->cd();
  myPad->SetGridx();   myPad->SetGridy();//  myPad->SetLogx();


  TLegend *leg = new TLegend(0.9,30,1.7,70,"","brCDN"); 
  leg->SetFillColor(0);

  TGraph *c[6];
  if (id!=2) {
  
  
    // Current ITS +++++++++++++++++++++++++++++++++++++++++
    Int_t color=1; Int_t linewidth=2;

    Int_t pi =0;
 
    DetectorK its("ALICE","ITS");
    its.MakeAliceCurrent(0,plusTPC);
    its.SetMaxRadiusOfSlowDetectors(0.01);
    its.SetAtLeastCorr(atLeastcorr);
    if (PbPb) its.SetdNdEtaCent(mult);
    its.SolveViaBilloir(0);
    Int_t color=1; Int_t linewidth=2;

    if (id==0)
      c[pi] = its.GetGraphRecoEfficiency(0,color,linewidth);
    else if (id==1)
      c[pi] = its.GetGraphRecoPurity(0,color,linewidth);
    else 
      c[pi] = its.GetGraphRecoFakes(0,color,linewidth);

    c[pi]->Draw("AC");

    leg->AddEntry(c[pi],"FastTool: Current ITS","l");


    // NEW SPD  +++++++++++++++++++++++++++++++++++++++++

    Int_t color=3; Int_t linewidth=2;
    Int_t pi =2;

    DetectorK its("ALICE","ITS");
    its.MakeAliceCurrent(0,plusTPC);
    its.SetAtLeastCorr(atLeastcorr);
    if (PbPb) its.SetdNdEtaCent(mult);
    its.SetRadius("bpipe",2.0);
    its.AddLayer("spd0", 2.2,1,1,1);  

    its.SetRadius("spd0",2.2); its.SetRadiationLength("spd0",X0); its.SetResolution("spd0",resRPhi,resZ);
    its.SetRadius("spd1",4.8);   its.SetRadiationLength("spd1",X0); its.SetResolution("spd1",resRPhi,resZ);
    its.SetRadius("spd2",9.1);   its.SetRadiationLength("spd2",X0); its.SetResolution("spd2",resRPhi,resZ);

    its.SetMaxRadiusOfSlowDetectors(0.1);
    its.SolveViaBilloir(0);

    if (id==0)
      c[pi] = its.GetGraphRecoEfficiency(0,color,linewidth);
    else if (id==1)
      c[pi] = its.GetGraphRecoPurity(0,color,linewidth);
    else 
      c[pi] = its.GetGraphRecoFakes(0,color,linewidth);
   
    c[pi]->Draw("C");

    leg->AddEntry(c[pi],"FastTool: \"New SPDs\"","l");


    // ALL NEW +++++++++++++++++++++++++++++++++++++++++++

    color=4; Int_t linewidth=2;
    Int_t pi =1; 


    // for a 0.8,0.2 weight configuration
  
    DetectorK *itsU = new DetectorK((char*)"ALICE",(char*)"ITS");
    itsU->SetAtLeastCorr(atLeastcorr);
    itsU->AddLayer((char*)"bpipe", 2.0,0.0022); // beam pipe
    itsU->AddLayer((char*)"vertex",  0,     0); // dummy vertex for matrix calculation
  
    itsU->AddLayer("ddd1",  2.2 ,  X0, resRPhi, resZ); 
    itsU->AddLayer("ddd2",  3.8 ,  X0, resRPhi, resZ); 
    itsU->AddLayer("ddd3",  6.8 ,  X0, resRPhi, resZ); 
    itsU->AddLayer("ddd4", 12.4 ,  X0, resRPhi, resZ); 
    itsU->AddLayer("ddd5", 23.5 ,  X0, resRPhi, resZ); 
    itsU->AddLayer("ddd6", 39.6 ,  X0, resRPhi, resZ); 
    //    itsU->AddLayer("ddd6", 42.6 ,  X0, resRPhi, resZ); 
    itsU->AddLayer("ddd7", 43.0 ,  X0, resRPhi, resZ); 
    //    itsU->AddLayer("ddd8", 43.4 ,  X0, resRPhi, resZ); 


    if (PbPb) itsU->SetdNdEtaCent(mult);
    if(plusTPC) itsU->AddTPC(0.1,0.1);
    itsU->SetMaxRadiusOfSlowDetectors(0.1);
    itsU->SolveViaBilloir(0);
    itsU->PrintLayout();

    if (id==0)
      c[pi] = itsU->GetGraphRecoEfficiency(0,color,linewidth);
    else if (id==1)
      c[pi] = itsU->GetGraphRecoPurity(0,color,linewidth);
    else 
      c[pi] = itsU->GetGraphRecoFakes(0,color,linewidth);
    c[pi]->Draw("C");

    leg->AddEntry(c[pi],"FastTool: \"All New\" ","l");

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 // ALL NEW - double outer layer +++++++++++++++++++++++++++++++++++++

    color=2; Int_t linewidth=2;
    Int_t pi =3; 


    // for a 0.8,0.2 weight configuration
  
    DetectorK *itsU = new DetectorK((char*)"ALICE",(char*)"ITS");
    itsU->SetAtLeastCorr(atLeastcorr);
    itsU->AddLayer((char*)"bpipe", 2.0,0.0022); // beam pipe
    itsU->AddLayer((char*)"vertex",  0,     0); // dummy vertex for matrix calculation
  
    itsU->AddLayer("ddd1",  2.2 ,  X0, resRPhi, resZ); 
    itsU->AddLayer("ddd2",  3.8 ,  X0, resRPhi, resZ); 
    itsU->AddLayer("ddd3",  6.8 ,  X0, resRPhi, resZ); 
    itsU->AddLayer("ddd4", 12.4 ,  X0, resRPhi, resZ); 
    itsU->AddLayer("ddd5", 23.5 ,  X0, resRPhi, resZ); 
    itsU->AddLayer("ddd6", 39.6 ,  X0, resRPhi, resZ); 
    itsU->AddLayer("ddd8", 40.0 ,  X0, resRPhi, resZ); 
    itsU->AddLayer("ddd7", 43.0 ,  X0, resRPhi, resZ); 
    itsU->AddLayer("ddd9", 43.4 ,  X0, resRPhi, resZ); 


    if (PbPb) itsU->SetdNdEtaCent(mult);
    if(plusTPC) itsU->AddTPC(0.1,0.1);
    itsU->SetMaxRadiusOfSlowDetectors(0.1);
    itsU->SolveViaBilloir(0);
    itsU->PrintLayout();

    if (id==0)
      c[pi] = itsU->GetGraphRecoEfficiency(0,color,linewidth);
    else if (id==1)
      c[pi] = itsU->GetGraphRecoPurity(0,color,linewidth);
    else 
      c[pi] = itsU->GetGraphRecoFakes(0,color,linewidth);
    c[pi]->Draw("C");

    leg->AddEntry(c[pi],"FastTool: \"All New\" (2x double layer)","l");

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  }

  char h[100];
  if (PbPb==0) 
    sprintf(h,"-fran");
  else 
    sprintf(h,"-Anna");

  TFile f1(Form("root/FastVsSlow_CurrentITS-PbPb%s.root",h));
  TFile f2(Form("root/FastVsSlow_NewSPDs-PbPb%s.root",h));
  TFile f3(Form("root/FastVsSlow_AllNew-PbPb%s.root",h));
  TFile f4(Form("root/FastVsSlow_AllNew-9-PbPb%s.root",h));

  //  TFile f1(Form("root/FastVsSlow_CurrentITS%s-fran.root",h));
  //  TFile f2(Form("root/FastVsSlow_NewSPDs%s-fran.root",h));
  //  TFile f3(Form("root/FastVsSlow_AllNew%s-fran.root",h));

  TH1F *eff1 = 0;
  TH1F *eff2 = 0;
  TH1F *eff3 = 0;
  TH1F *eff4 = 0;
  if (id==0) {
    eff1 = (TH1F*)f1.Get("efficiency");
    eff2 = (TH1F*)f2.Get("efficiency");
    eff3 = (TH1F*)f3.Get("efficiency");
    eff4 = (TH1F*)f4.Get("efficiency");
    eff1->GetYaxis()->SetTitle("efficiency (%)");
  } else if (id==1) {
    eff1 = (TH1F*)f1.Get("purity");
    eff2 = (TH1F*)f2.Get("purity");
    eff3 = (TH1F*)f3.Get("purity");
    eff4 = (TH1F*)f4.Get("purity");
      eff1->GetYaxis()->SetTitle("purity (%)");
  } else if (id==2) {
    eff1 = (TH1F*)f1.Get("annaEff");
    eff2 = (TH1F*)f2.Get("annaEff");
    eff3 = (TH1F*)f3.Get("annaEff");
    eff4 = (TH1F*)f4.Get("annaEff");
    eff1->GetYaxis()->SetTitle("Overall efficiency (%)");
  } else if (id==3) {
    eff1 = (TH1F*)f1.Get("fake");
    eff2 = (TH1F*)f2.Get("fake");
    eff3 = (TH1F*)f3.Get("fake");
    eff4 = (TH1F*)f4.Get("fake");
    eff1->GetYaxis()->SetTitle("Fake ratio (%)");
  }

  eff1->SetMarkerStyle(21); eff1->SetMarkerColor(1);
  eff2->SetMarkerStyle(21); eff2->SetMarkerColor(3);
  eff3->SetMarkerStyle(21); eff3->SetMarkerColor(4);
  eff4->SetMarkerStyle(21); eff4->SetMarkerColor(2);

  leg->AddEntry(eff1,"FullMC: Current ITS","PE");
  leg->AddEntry(eff2,"FullMC: \"New SPDs\"","PE");
  leg->AddEntry(eff3,"FullMC: \"All New\" ","PE");
  leg->AddEntry(eff4,"FullMC: \"All New\" (2x double layer)","PE");

  eff1->SetMinimum(0.4); eff1->SetMaximum(100);
  eff1->DrawCopy("E");
  eff2->DrawCopy("sameE");
  eff4->DrawCopy("sameE");
  eff3->DrawCopy("sameE");
  if (id!=2) {
    c[0]->Draw("C");
    c[1]->Draw("C");
    c[2]->Draw("C");
    c[3]->Draw("C");
  }
  eff2->DrawCopy("sameE");
  eff4->DrawCopy("sameE");
  eff3->DrawCopy("sameE");

  
  leg->Draw();
 



  TPaveText *pt = 0;
  if (id!=3) 
   pt = new TPaveText(0.4,0.1,1.76,30);
  else
   pt = new TPaveText(0.4,70,1.76,100);
    
  pt->SetBorderSize(1); // no shadow
  pt->SetTextFont(12);
  TText *t1 = pt->AddText("FastTool settings: "); t1->SetTextFont(32); // bold

  pt->AddText(Form("   Tracked particles: Pions;   Average rapidity: 0.45; dN_{ch}/d#eta = %d ",mult));

  //  pt->AddText("\"New SPDs\": layer radii: r = {2.2,4.8,9.1} cm");
  //  pt->AddText("\"All New\: layer radii: r = {2.2,3.8,6.8,...} cm");
  //  pt->AddText(Form("    New layer prop.: X/X_{0}=%1.1lf%%;  #sigma_{r#phi,z}=%1.0lf#mum",X0*100,resZ*1e4));
  
  TText *t2 = pt->AddText("FullMC settings: "); t2->SetTextFont(32); // bold
  if (PbPb==0) {
    pt->AddText("   Generator: AliGenHIJINGpara (parametrized PbPb event)");
    pt->AddText("   dN_{ch.pr.}/d#eta = 2070");
    pt->AddText("   Track selection: Pions, |#eta|<0.9");
  } else {
    pt->AddText("   Generator: AliGenHijing (modified);  #sqrt{s_{NN}} = 5.5 TeV");
    pt->AddText("   dN_{ch.pr.}/d#eta = 2410; Impactparameter range: b#in(0,5)  #rightarrow central PbPb ");
    pt->AddText("   Track selection: Pions, |#eta|<0.9");
  }
  //  pt->SetLabel("Settings");
  pt->SetTextAlign(12);
  pt->SetFillColor(0);
  pt->Draw();





  if (PbPb==0) {    
    myCan->SaveAs(Form("FastVsSlowSim-Eff-%d.pdf",id));
    myCan->SaveAs(Form("FastVsSlowSim-Eff-%d.eps",id));
  }else{
    myCan->SaveAs(Form("FastVsSlowSim-Eff-PbPb-%d.pdf",id));
    myCan->SaveAs(Form("FastVsSlowSim-Eff-PbPb-%d.eps",id));
  }




}



// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Extraction


void GetDetectorRadii(TArrayD *rmin,TArrayD *rmax) {
  // Loads geometry of the ITS upgrade
  AliITSsegmentationUpgrade *seg=new AliITSsegmentationUpgrade();
  Int_t nlayers = seg->GetNLayers();
  rmin->Set(nlayers);   rmax->Set(nlayers);
  for (Int_t i=0; i<nlayers; i++) {
    rmin->AddAt(seg->GetRadius(i),i);
    rmax->AddAt(seg->GetRadius(i)+seg->GetThickness(i),i);
  }
}  

void PrintDetectorGeometry() {
  gSystem->Load("libITSUpgradeBase");
  gSystem->Load("libITSUpgradeSim");
  
  TArrayD rmin(0); 
  TArrayD rmax(0); 
  GetDetectorRadii(&rmin,&rmax);
  for (Int_t i=0; i<rmin.GetSize(); i++) {
    cout<<i<<": (rmin,rmax)=("<<rmin.At(i)<<","<<rmax.At(i)<<")cm"<<endl;
  }
}


void CountTrackableMCs(TH1F *hAllMC=0, Bool_t onlyPrims=0,Bool_t onlyPion=0);
void CountPrimaries(TH1F *hMultCount);


void ExtractOutputHistos(Bool_t onlyPrims=0,Bool_t onlyPion=0,Int_t plotFlag=0) {

  //  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);

  const Int_t nbins=20;
  Double_t ptmin=0.06;//04;
  Double_t ptmax=2.0;//GeV
  Double_t logxmin = TMath::Log10(ptmin);
  Double_t logxmax = TMath::Log10(ptmax);
  Double_t binwidth = (logxmax-logxmin)/(nbins+1);
  enum {nb=nbins+1};
  Double_t xbins[nb];
  xbins[0] = ptmin;
  for (Int_t i=1;i<=nbins;i++) {
    xbins[i] = ptmin + TMath::Power(10,logxmin+(i)*binwidth);
    //    cout<<xbins[i]<<endl;
  }
  //  TH1F *h = new TH1F("h","hist with log x axis",nbins,xbins);

  TH1F *hMultCount = new TH1F("mult","averaged multiplicity (charg. prim)",80,-4.,4.);
  hMultCount->GetXaxis()->SetTitle("eta");
  hMultCount->GetYaxis()->SetTitle("N/d#eta");

  TH1F *hAllMC = new TH1F("allMC","All Tracks MC primaries",nbins,xbins);
  TH1F *hAllFound = new TH1F("allFound","All Tracks found",nbins,xbins);
  TH1F *hImperfect = new TH1F("imperfect","Imperfect tracks",nbins,xbins);
  TH1F *hPerfect = new TH1F("perfect","Perfect tracks",nbins,xbins);
  TH1F *hEff = new TH1F("efficiency","Efficiency (Perfect tracks in \"ALL MC\")",nbins,xbins);
  TH1F *hFake = new TH1F("fake","Fake tracks (Inperfect tracks in \"ALL MC\")",nbins,xbins);
  TH1F *hPurity = new TH1F("purity","Purity (Perfect tracks in \"All Found\")",nbins,xbins);
  TH1F *hAnna = new TH1F("annaEff","AnnalisaEff ",nbins,xbins);
  TH1F *hNoMCTrack = new TH1F("noMCtrack","noMCtrack ",nbins,xbins);

  TH1F *hEta = new TH1F("","",50,-2,2);
  //  TH1F *hEtaMC = new TH1F("","",50,-2,2);

  TH2D *h2Ddca = new TH2D("dca2D","DCAvsPt2D",nbins,xbins,50,-0.05,0.05);
  TH2D *h2Dpt = new TH2D("dPt2D","dPtdvsPt2D",nbins,xbins,50,-25,25);

  // open run loader and load gAlice, kinematics and header
  AliRunLoader* runLoader = AliRunLoader::Open("galice.root");
  if (!runLoader) {
    Error("Check kine", "getting run loader from file %s failed",
          "galice.root");
    return;
  }
  runLoader->LoadgAlice();
  gAlice = runLoader->GetAliRun();
  if (!gAlice) {
    Error("Check kine", "no galice object found");
    return;
  }
  runLoader->LoadHeader();
  runLoader->LoadKinematics();

  TFile* esdFile = TFile::Open("AliESDs.root");
  if (!esdFile || !esdFile->IsOpen()) {
    Error("CheckESD", "opening ESD file %s failed", "AliESDs.root");
    return;
  }
  AliESDEvent *esd = new AliESDEvent();
  TTree* tree = (TTree*) esdFile->Get("esdTree");
  if (!tree) {
    Error("CheckESD", "no ESD tree found");
    return;
  }
  esd->ReadFromTree(tree);
  
  Int_t nTrackTotalMC = 0;
  Int_t nTrackFound = 0;
  Int_t nTrackImperfect = 0;
  Int_t nTrackPerfect = 0;
  Int_t nNoMCTrack = 0;

  
  for(Int_t iEv =0; iEv<tree->GetEntries(); iEv++){
    tree->GetEvent(iEv);
    runLoader->GetEvent(iEv);
    
    printf("+++ event %i (of %lld) +++++++++++++++++++++++  # ESDtracks: %d \n",iEv,tree->GetEntries()-1,esd->GetNumberOfTracks());
    Int_t nESDtracks = esd->GetNumberOfTracks();
    for (Int_t iTrack = 0; iTrack < nESDtracks; iTrack++) {
      AliESDtrack* track = esd->GetTrack(iTrack);
      if (!(iTrack%1000)) printf("event %i: ESD track count %d (of %d)\n",iEv,iTrack,nESDtracks);

      Int_t label = track->GetLabel();
  
      Int_t idx[12];
      //      Int_t ncl = track->GetITSclusters(idx);
   
      if(label<0) {
	//	cout<< " ESD track label " << label;
	//	cout<<"  ---> imperfect track (label "<<label<<"<0) !! -> track Pt: "<< track->Pt() << endl;
      }

      AliStack* stack = runLoader->Stack();
      //     nTrackTotalMC += stack->GetNprimary();
    

      TParticle* particle = stack->Particle(TMath::Abs(label)); 
      Double_t pt = track->Pt();
      
      if(particle) {

	if (TMath::Abs(particle->Eta())>etaCut) continue;

	Double_t ptMC = particle->Pt();

	// Efficiencies
	if (onlyPion && TMath::Abs(particle->GetPdgCode())!=211) continue;

	if ( (!onlyPrims) || stack->IsPhysicalPrimary(TMath::Abs(label))) {
	  //  cout<<" # clusters "<<ncl<<endl;

	  nTrackFound++;
	  hAllFound->Fill(ptMC);
	  hEta->Fill(track->Eta());
	  
	  if (label<0) {
	    nTrackImperfect++;
	    hImperfect->Fill(ptMC);
	  } else {
	    nTrackPerfect++;
	    hPerfect->Fill(ptMC);
	  }

	}


	// following only for "true tracks, pions

	if(particle->Pt() < 0.001)continue;
	if (TMath::Abs(particle->GetPdgCode())!=211) continue;
	if (label>0) {
	  
	  // Impact parameters for Pions only
	  Double_t dca = track->GetD(0,0,0.5);
	  h2Ddca->Fill(ptMC,dca);
	  
	  // Pt resolution for Pions only
	  Double_t dPt = (pt-ptMC)/ptMC*100;
	  h2Dpt->Fill(ptMC,dPt);
	}

      } else {
	nNoMCTrackFound++;
	hNoMCTrack->Fill(pt);
	cout<<" according MC particle not found"<<endl;
      }
      
    } //entries track esd
  

  }//entries tree
  runLoader->UnloadHeader();
  runLoader->UnloadKinematics();
  delete runLoader;

 
  // Count trackable MC tracks
  CountTrackableMCs(hAllMC, onlyPrims, onlyPion);


  // Count trackable MC tracks
  CountPrimaries(hMultCount);

 


  // Get Errors right
  hMultCount->Sumw2();
  hAllMC->Sumw2();   
  hAllFound->Sumw2();
  hPerfect->Sumw2(); 
  hImperfect->Sumw2(); 
  h2Dpt->Sumw2();
  h2Ddca->Sumw2();

  // -- Global efficienies

  nTrackTotalMC = hAllMC->GetEntries();
  Double_t eff = ((Double_t)nTrackPerfect)/nTrackTotalMC;
  printf("-> Total number of events: %lld -> MCtracks %d -> nPerfect %d  -> Eff: %3.2lf \n",
	 tree->GetEntries(),nTrackTotalMC,nTrackPerfect,eff);

  Double_t purity = ((Double_t)nTrackPerfect)/nTrackFound;
  printf("-> Total number of events: %lld -> FoundTracks %d -> nPerfect %d  -> Purity: %3.2lf \n",
	 tree->GetEntries(),nTrackFound,nTrackPerfect,purity);

  // Efficiencies - and normalize to 100%

  TF1 f1("f1","100+x*0",0.,1.e3);

  hPurity->Divide(hPerfect,hAllFound,1,1,"b"); 
  hPurity->Multiply(&f1);
  hPurity->SetMarkerColor(kGreen);
  hPurity->SetMarkerStyle(21);
  hPurity->GetXaxis()->SetTitle("transverse momentum p_{t} (GeV)");
  hPurity->SetStats(0);

  hPurity->GetYaxis()->SetRangeUser(0,100);
  hPurity->SetTitle("Efficiency & Purity");

  hEff->Divide(hPerfect,hAllMC,1,1,"b");
  hEff->Multiply(&f1);
  hEff->GetXaxis()->SetTitle("transverse momentum p_{t} (GeV)");
  hEff->SetMarkerColor(kBlue);
  hEff->SetMarkerStyle(21);
  hEff->SetStats(0);

  hFake->Divide(hImperfect,hAllMC,1,1,"b");
  hFake->Multiply(&f1);
  hFake->GetXaxis()->SetTitle("transverse momentum p_{t} (GeV)");
  hFake->SetMarkerColor(kRed);
  hFake->SetMarkerStyle(21);
  hFake->SetStats(0);


  hAnna->Divide(hAllFound,hAllMC,1,1,"b");
  hAnna->Multiply(&f1);
  hAnna->GetXaxis()->SetTitle("transverse momentum p_{t} (GeV)");
  hAnna->SetMarkerColor(kBlack);
  hAnna->SetMarkerStyle(21);
  hAnna->SetStats(0);

  TCanvas *c1 = new TCanvas("c1","NoMCTrackFound");//,200,10,900,900);
  TVirtualPad *pad =   c1->cd();
  pad->SetGridx();   pad->SetGridy();
  hNoMCTrack->Draw();

  TCanvas *c2 = new TCanvas("c2","Eff&Purity");//,200,10,900,900);
  TVirtualPad *pad =   c2->cd();
  pad->SetGridx();   pad->SetGridy();
  //  pad->SetLogx();

  hPurity->Draw("E");
  hEff->Draw("Same E");
  hFake->Draw("Same E");
  hAnna->Draw("Same E");

  TLegend *leg = new TLegend(0.1,0.8,0.6,0.9);leg->SetFillColor(0);
  leg->AddEntry(hPurity,"Purity (\"Perfect tracks\" within \"Found Tracks\")","PE");
  leg->AddEntry(hEff,"Efficiency (\"Perfect tracks\" within \"MC findable Tracks\")","PE");
  leg->AddEntry(hFake,"Fake (\"Inperfect tracks\" within \"MC findable Tracks\")","PE");
  leg->AddEntry(hAnna,"AnnaLisa - Efficiency (\"Found tracks\" within \"MC findable Tracks\")","PE");
  leg->Draw();


  if (plotFlag==1){
    hAllMC->GetXaxis()->SetTitle("transverse momentum p_{t} (GeV)");
    hAllMC->Draw();  // MC pt distribution
    hAllFound->SetLineColor(2);
    hAllFound->Draw("same");  // MC pt distribution
  }
 
 
  /*

  .L ~/ITSupgrade/BuildDetector/DetectorK.cxx+
  
  // All NEW
  DetectorK its("ALICE","ITS");
  its.MakeAliceAllNew(0);
  its.SetMaxRadiusOfSlowDetectors(0.01);
  its.SolveViaBilloir(0);
  TGraph *c = its.GetGraphRecoEfficiency(0,3,2);
  c->Draw("C");


  // Current
  DetectorK its("ALICE","ITS");
  its.MakeAliceCurrent(0,0);
  its.SetMaxRadiusOfSlowDetectors(0.01);
  its.SolveViaBilloir(0);
  TGraph *c = its.GetGraphRecoEfficiency(0,4,2);
  c->Draw("C");

  */

  TCanvas *c3 = new TCanvas("c3","impact");//,200,10,900,900);
  c3->Divide(2,1); c3->cd(1);
  // Impact parameter

  // Impact parameter resolution ---------------
  h2Ddca->Draw("colz");
  h2Ddca->FitSlicesY() ;
  TH2D *dcaM = (TH2D*)gDirectory->Get("dca2D_1"); dcaM->Draw("same");
  TH2D *dcaRMS = (TH2D*)gDirectory->Get("dca2D_2"); //dcaRMS->Draw();
  TGraphErrors *d0 = new TGraphErrors(); 
  for (Int_t ibin =1; ibin<=dcaRMS->GetXaxis()->GetNbins(); ibin++) {
    d0->SetPoint(     ibin-1,dcaRMS->GetBinCenter(ibin),dcaRMS->GetBinContent(ibin)*1e4); // microns
    d0->SetPointError(ibin-1,0,dcaRMS->GetBinError(ibin)*1e4); // microns
  }
  d0->SetMarkerStyle(21);
  d0->SetMaximum(200);  d0->SetMinimum(0);
  d0->GetXaxis()->SetTitle("transverse momentum p_{t} (GeV)");
  d0->GetYaxis()->SetTitle("R-#phi Pointing Resolution (#mum)");
  d0->SetName("dca");  d0->SetTitle("DCAvsPt");

  c3->cd(1);  h2Ddca->Draw("surf2");
  c3->cd(2);  d0->Draw("APE");

  // PT RESOLUTION ------------
  TCanvas *c4 = new TCanvas("c4","pt resolution");//,200,10,900,900);
  c4->Divide(2,1); c4->cd(1);
  // Impact parameter
  h2Dpt->Draw("colz");
  h2Dpt->FitSlicesY() ;
  TH2D *dPtM = (TH2D*)gDirectory->Get("dPt2D_1"); dPtM->Draw("same");
  TH2D *dPtRMS = (TH2D*)gDirectory->Get("dPt2D_2"); // dPtRMS->Draw("");
  TGraphErrors *gPt = new TGraphErrors(); 
  for (Int_t ibin =1; ibin<=dPtRMS->GetXaxis()->GetNbins(); ibin++) {
    gPt->SetPoint(     ibin-1,dPtRMS->GetBinCenter(ibin),dPtRMS->GetBinContent(ibin)); 
    gPt->SetPointError(ibin-1,0,dPtRMS->GetBinError(ibin)); 
  }
  gPt->SetMarkerStyle(21);
  gPt->SetMaximum(20);  gPt->SetMinimum(0);
  gPt->GetXaxis()->SetTitle("transverse momentum p_{t} (GeV)");
  gPt->GetYaxis()->SetTitle("relative momentum resolution (%)");
  gPt->SetName("dPt");  gPt->SetTitle("DPTvsPt");

  c4->cd(1);  h2Dpt->Draw("surf2");
  c4->cd(2);  gPt->Draw("APE");


  // EXPORT --------

  TFile f("histos.root","RECREATE");

  hMultCount->Write();
  hAllMC->Write();
  hAllFound->Write();
  hImperfect->Write();
  hPerfect->Write();
  hNoMCTrack->Write();

  hPurity->Write();
  hEff->Write();
  hFake->Write();
  hAnna->Write();

  h2Ddca->Write();
  d0->Write();

  h2Dpt->Write();
  gPt->Write();

  f.Close();

  return;

}

void CountTrackableMCs(TH1F *hAllMC, Bool_t onlyPrims,Bool_t onlyPion) {
  
  gSystem->Load("libITSUpgradeBase");
  gSystem->Load("libITSUpgradeSim");

 // open run loader and load gAlice, kinematics and header
  AliRunLoader* runLoader = AliRunLoader::Open("galice.root");
  if (!runLoader) {
    Error("Check kine", "getting run loader from file %s failed",
          "galice.root");
    return;
  }
  runLoader->LoadHeader();
  runLoader->LoadKinematics();
  runLoader->LoadTrackRefs();

  AliLoader *dl = runLoader->GetDetectorLoader("ITS");

  //Trackf
  TTree *trackRefTree = 0x0; 
  TClonesArray *trackRef = new TClonesArray("AliTrackReference",1000);

  //  TH1F *hRef = new TH1F("","",100,0,100);
  TH1F *hR = new TH1F("","",100,0,100);
  if (hAllMC==0) hAllMC = new TH1F("","",100,0.1,2);
  Float_t ptmin = hAllMC->GetBinCenter(1)-hAllMC->GetBinWidth(1)/2;
  Float_t ptmax = hAllMC->GetBinCenter(hAllMC->GetNbinsX())+hAllMC->GetBinWidth(hAllMC->GetNbinsX())/2;
  //  Int_t nAllMC = 0;

  // Detector geometry
  TArrayD rmin(0);   TArrayD rmax(0); 
  GetDetectorRadii(&rmin,&rmax);
  TArrayI nLaySigs(rmin.GetSize());

  printf("Counting trackable MC tracks ...\n");
  
  for(Int_t iEv =0; iEv<runLoader->GetNumberOfEvents(); iEv++){
    Int_t nTrackableTracks = 0;
    runLoader->GetEvent(iEv);
    AliStack* stack = runLoader->Stack();  
    printf("+++ event %i (of %d) +++++++++++++++++++++++  # total MCtracks: %d \n",iEv,runLoader->GetNumberOfEvents()-1,stack->GetNtrack());

    trackRefTree=runLoader->TreeTR();
    TBranch *br = trackRefTree->GetBranch("TrackReferences");
    if(!br) {
      printf("no TR branch available , exiting \n");
      return;
    }
    br->SetAddress(&trackRef);

    // init the trackRef tree 
    trackRefTree=runLoader->TreeTR();
    trackRefTree->SetBranchAddress("TrackReferences",&trackRef);
 
    // Count trackable MC tracks
    for (Int_t iMC=0; iMC<stack->GetNtrack(); iMC++) {

      TParticle* particle = stack->Particle(iMC); 
      if (TMath::Abs(particle->Eta())>etaCut) continue;
      if (onlyPrims && !stack->IsPhysicalPrimary(iMC)) continue;
      if (onlyPion && TMath::Abs(particle->GetPdgCode())!=211) continue;


      Bool_t isTrackable = 0;
      nLaySigs.Reset(0);
 
      trackRefTree->GetEntry(stack->TreeKEntry(iMC));
      Int_t nref=trackRef->GetEntriesFast();
      for(Int_t iref =0; iref<nref; iref++){
	AliTrackReference *trR = (AliTrackReference*)trackRef->At(iref);
	if(!trR) continue;
	if(trR->DetectorId()!=AliTrackReference::kITS) continue;
	Float_t radPos = trR->R();
	hR->Fill(radPos);
	for (Int_t il=0; il<rmin.GetSize();il++) {
	  if (radPos>=rmin.At(il)-0.1 && radPos<=rmax.At(il)+0.1) {
	    //	    cout<<"  in Layer "<<il<<" "<<radPos;
	    nLaySigs.AddAt(1.,il);
	    //	    cout<<" "<<nLaySigs.At(il)<<endl;
	  }
	}
      }

      if (nLaySigs.GetSum()>=3) {
	isTrackable =1;
	//	cout<<nLaySigs.GetSum()<<endl;
      }
      
      if (isTrackable) {
	Double_t ptMC = particle->Pt();
	//	Double_t etaMC = particle->Eta();
	//	if (ptMC>ptmin&&ptMC<ptmax) {nTrackableTracks++;hAllMC->Fill(ptMC);}
	if (ptMC>ptmin) {nTrackableTracks++;hAllMC->Fill(ptMC);}

      }

      
    } // entries tracks MC
    printf(" -> trackable MC tracks: %d (%d)\n",nTrackableTracks,hAllMC->GetEntries());
  }//entries Events
  

  hR->DrawCopy();
  hAllMC->DrawCopy();
  runLoader->UnloadHeader();
  runLoader->UnloadKinematics();
  delete runLoader;

}

void CountPrimaries(TH1F *hMultCount) {

  if (hMultCount==0) hMultCount = new TH1F("mult","averaged multiplicity (charg. prim)",80,-4.,4.);
  
  AliRunLoader *rl = AliRunLoader::Open("galice.root");
  rl->SetKineFileName("Kinematics.root");
  rl->LoadHeader();
  rl->LoadKinematics(); 
  Int_t nEvents = rl->GetNumberOfEvents();
  cout<< "N events "<<nEvents<<endl;
  for(Int_t iEv=0; iEv<nEvents; iEv++){
    rl->GetEvent(iEv);
    AliStack *s = rl->Stack();
    for(Int_t iP=0; iP<s->GetNtrack(); iP++ ){
      TParticle *p = s->Particle(iP);
      if (!(s->IsPhysicalPrimary(iP))) continue;
      Float_t eta = p->Eta();
      if (p->Pt()>0.06) {
	hMultCount->Fill(eta);
      }
    }
  }

  hMultCount->DrawCopy();
  rl->UnloadHeader();
  rl->UnloadKinematics();
  delete rl;



}


void plotMerged(Bool_t onlyPlot=0) {

  gStyle->SetPalette(1);
 
  TFile f("histoSum.root","UPDATE");

  TH1F* hAllMC = f.Get("allMC");
  TH1F* hAllFound= f.Get("allFound");
  TH1F* hImperfect= f.Get("imperfect");
  TH1F* hPerfect= f.Get("perfect");
  TH1F* hNoMCTrack= f.Get("noMCtrack");
  
  
  // have to be recalculated
  TH1F* hPurity = f.Get("purity");
  TH1F* hEff= f.Get("efficiency");
  TH1F* hFake= f.Get("fake");
  TH1F* hAnna= f.Get("annaEff");

  TH2D* h2Ddca= f.Get("dca2D");
  TGraphErrors *d0= f.Get("dca");

  TH2D* h2Dpt= f.Get("dPt2D");
  TGraphErrors *gPt= f.Get("dPt");


  if (!onlyPlot) {
    /*    // Get Errors right
    hAllMC->Sumw2();   
    hAllFound->Sumw2();
    hPerfect->Sumw2(); 
    hImperfect->Sumw2(); 
    h2Dpt->Sumw2();
    h2Ddca->Sumw2();
    */

    // Efficiencies - and normalize to 100%
    
    TF1 f1("f1","100+x*0",0.,1.e3);
    
    hPurity->Divide(hPerfect,hAllFound,1,1,"b"); 
    hPurity->Multiply(&f1);
    hPurity->SetMarkerColor(kGreen);
    hPurity->SetMarkerStyle(21);
    hPurity->GetXaxis()->SetTitle("transverse momentum p_{t} (GeV)");
    hPurity->SetStats(0);
    
    hPurity->GetYaxis()->SetRangeUser(0,100);
    hPurity->SetTitle("Efficiency & Purity");
    
    hEff->Divide(hPerfect,hAllMC,1,1,"b");
    hEff->Multiply(&f1);
    hEff->GetXaxis()->SetTitle("transverse momentum p_{t} (GeV)");
    hEff->SetMarkerColor(kBlue);
    hEff->SetMarkerStyle(21);
    hEff->SetStats(0);
    
    hFake->Divide(hImperfect,hAllMC,1,1,"b");
    hFake->Multiply(&f1);
    hFake->GetXaxis()->SetTitle("transverse momentum p_{t} (GeV)");
    hFake->SetMarkerColor(kRed);
    hFake->SetMarkerStyle(21);
    hFake->SetStats(0);
    
    hAnna->Divide(hAllFound,hAllMC,1,1,"b");
    hAnna->Multiply(&f1);
    hAnna->GetXaxis()->SetTitle("transverse momentum p_{t} (GeV)");
    hAnna->SetMarkerColor(kBlack);
    hAnna->SetMarkerStyle(21);
    hAnna->SetStats(0);
    
    
    // Impact parameter resolution ---------------
    TCanvas *c3 = new TCanvas("c3","impact");//,200,10,900,900);
    c3->Divide(2,1); c3->cd(1);
    h2Ddca->DrawCopy("colz");
    h2Ddca->FitSlicesY() ;
    TH2D *dcaM = (TH2D*)gDirectory->Get("dca2D_1"); dcaM->Draw("same");
    TH2D *dcaRMS = (TH2D*)gDirectory->Get("dca2D_2"); //dcaRMS->Draw();
    TGraphErrors *d0 = new TGraphErrors(); 
    for (Int_t ibin =1; ibin<=dcaRMS->GetXaxis()->GetNbins(); ibin++) {
      d0->SetPoint(     ibin-1,dcaRMS->GetBinCenter(ibin),dcaRMS->GetBinContent(ibin)*1e4); // microns
      d0->SetPointError(ibin-1,0,dcaRMS->GetBinError(ibin)*1e4); // microns
    }
    d0->SetMarkerStyle(21);
    d0->SetMaximum(200);  d0->SetMinimum(0);
    d0->GetXaxis()->SetTitle("transverse momentum p_{t} (GeV)");
    d0->GetYaxis()->SetTitle("R-#phi Pointing Resolution (#mum)");
    d0->SetName("dca");  d0->SetTitle("DCAvsPt");
    //  c3->cd(1);  h2Ddca->Draw("surf2");
    c3->cd(2);  d0->Draw("APE");
    
    // PT RESOLUTION ------------
    TCanvas *c4 = new TCanvas("c4","pt resolution");//,200,10,900,900);  
    c4->Divide(2,1); c4->cd(1);
    h2Dpt->DrawCopy("colz");
    h2Dpt->FitSlicesY() ;
    TH2D *dPtM = (TH2D*)gDirectory->Get("dPt2D_1"); dPtM->Draw("same");
    TH2D *dPtRMS = (TH2D*)gDirectory->Get("dPt2D_2"); // dPtRMS->Draw("");
    TGraphErrors *gPt = new TGraphErrors(); 
    for (Int_t ibin =1; ibin<=dPtRMS->GetXaxis()->GetNbins(); ibin++) {
      gPt->SetPoint(     ibin-1,dPtRMS->GetBinCenter(ibin),dPtRMS->GetBinContent(ibin)); 
      gPt->SetPointError(ibin-1,0,dPtRMS->GetBinError(ibin)); 
    }
    gPt->SetMarkerStyle(21);
    gPt->SetMaximum(20);  gPt->SetMinimum(0);
    gPt->GetXaxis()->SetTitle("transverse momentum p_{t} (GeV)");
    gPt->GetYaxis()->SetTitle("relative momentum resolution (%)");
    gPt->SetName("dPt");  gPt->SetTitle("DPTvsPt");
    //  c4->cd(1);  h2Dpt->Draw("surf2");
    c4->cd(2);  gPt->Draw("APE");


    // overwrite with normalized graphs
    hPurity->Write();
    hEff->Write();
    hFake->Write();
    hAnna->Write();
    h2Ddca->Write();
    d0->Write();
    h2Dpt->Write();
    gPt->Write();
   
  }
  
  // Plots

  TCanvas *c2 = new TCanvas("c2","Eff&Purity");//,200,10,900,900);
  TVirtualPad *pad =   c2->cd();
  pad->SetGridx();   pad->SetGridy();
  //  pad->SetLogx();

  TLegend *leg = new TLegend(0.1,0.8,0.6,0.9);leg->SetFillColor(0);
  leg->AddEntry(hPurity,"Purity (\"Perfect tracks\" within \"Found Tracks\")","PE");
  leg->AddEntry(hEff,"Efficiency (\"Perfect tracks\" within \"MC findable Tracks\")","PE");
  leg->AddEntry(hFake,"Fake (\"Inperfect tracks\" within \"MC findable Tracks\")","PE");
  leg->AddEntry(hAnna,"AnnaLisa - Efficiency (\"Found tracks\" within \"MC findable Tracks\")","PE");
 

  hPurity->DrawCopy("E");
  hEff->DrawCopy("Same E");
  hFake->DrawCopy("Same E");
  hAnna->DrawCopy("Same E");
  leg->Draw();

  c2->SaveAs("EffPlot.png");

  f.Close();



}

