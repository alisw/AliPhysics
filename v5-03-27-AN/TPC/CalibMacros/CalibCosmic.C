/*
Draw result of perfomance test:

aliroot -b -q  $ALICE_ROOT/TPC/scripts/loadTPCcalib.C $ALICE_ROOT/TPC/CalibMacros/CalibCosmic.C

  //gROOT->Macro("~/NimStyle.C"); 
  gSystem->AddIncludePath("-I$ALICE_ROOT/STAT")
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC")
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros")

  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
  .L $ALICE_ROOT/TPC/CalibMacros/CalibCosmic.C
  // init
  Init();
  SetDefaultCut();  // check defualt cut 
  //
  MakeDefaultPlots();

  gROOT->Macro("~/NimStyle.C");
  TFile f("cosmicPlots.root");
  TBrowser b
  b.Add(CosmicPlots,"CosmicPlot");

*/  

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "THnSparse.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TFile.h"
#include "TChain.h"
#include "TCut.h"
#include "TH3.h"
#include "TH2F.h"
#include "TProfile3D.h"
#include "TMath.h" 
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TTreeStream.h"
#include "AliExternalTrackParam.h"
#include "AliESDfriend.h"
#include "AliTPCcalibCosmic.h"
#include "TROOT.h"
#include "TPostScript.h"
#include "TStyle.h"
#include "AliTrackerBase.h"
#include "AliTPCExBEffective.h"
#include "TEntryList.h"
#include "TLegend.h"
#endif


class AliTPCcalibCosmic;
AliTPCcalibCosmic * cosmicScan =0;
TObjArray fitArr;
Int_t colors[3]={1,2,4};
TObjArray *picArray = new TObjArray();
const char * chsign[3]={"all", "Positive","Negative"};
void Init();
void SetRangeAll(Int_t axis, Float_t xmin, Float_t xmax);
void SetDefaultCut();
void  MakeDefaultPlots();



void CalibCosmic(){
  // init
  Init();
  SetDefaultCut(); 
  //
  MakeDefaultPlots();
}

void Init(){
  //
  //
  TH1::AddDirectory(0);
  TFile fcalib("TPCCosmicObjects.root");
  cosmicScan = ( AliTPCcalibCosmic *)fcalib.Get("cosmicTPC");
  TString axisName[9];
  axisName[0]  ="#Delta"; axisName[1]  ="N_{cl}";
  axisName[2]  ="DCA_{r}(cm)";
  axisName[3]  ="z (cm)"; axisName[4]  ="sin(#phi)";
  axisName[5]  ="tan(#theta)"; axisName[6]  ="1/p_{t} (1/GeV)";
  axisName[7]  ="p_{t} (GeV)"; axisName[8]  ="alpha";

  {
  for (Int_t ivar=0;ivar<6;ivar++){
    for (Int_t ivar2=0;ivar2<9;ivar2++){      
      cosmicScan->fHistoDelta[ivar]->GetAxis(ivar2)->SetName(axisName[ivar2].Data());
      cosmicScan->fHistoDelta[ivar]->GetAxis(ivar2)->SetTitle(axisName[ivar2].Data());
      cosmicScan->fHistoPull[ivar]->GetAxis(ivar2)->SetName(axisName[ivar2].Data());
      cosmicScan->fHistoPull[ivar]->GetAxis(ivar2)->SetTitle(axisName[ivar2].Data());
    }
  }
  }
}

void SetRangeAll(Int_t axis, Float_t xmin, Float_t xmax){

  for (Int_t i=0;i<6;i++){
    //
    cosmicScan->fHistoDelta[i]->GetAxis(axis)->SetRangeUser(xmin,xmax);
    cosmicScan->fHistoPull[i]->GetAxis(axis)->SetRangeUser(xmin,xmax);
  }
}


void SetDefaultCut(){
  for (Int_t i=0;i<6;i++){
    //
    //cut on number of clusters
    cosmicScan->fHistoDelta[i]->GetAxis(1)->SetRangeUser(130,200);
    cosmicScan->fHistoPull[i]->GetAxis(1)->SetRangeUser(130,200);
    //cut on DCA r
    cosmicScan->fHistoDelta[i]->GetAxis(2)->SetRangeUser(0,15);
    cosmicScan->fHistoPull[i]->GetAxis(2)->SetRangeUser(0,15);
    //cut on z at
    cosmicScan->fHistoDelta[i]->GetAxis(3)->SetRangeUser(-20,20);
    cosmicScan->fHistoPull[i]->GetAxis(3)->SetRangeUser(-20,20);
  }
  cosmicScan->fHistoDelta[0]->GetAxis(0)->SetRangeUser(-1,1);
  cosmicScan->fHistoDelta[1]->GetAxis(0)->SetRangeUser(-1,1);
  cosmicScan->fHistoDelta[4]->GetAxis(0)->SetRangeUser(-0.1,0.1);
}

TH2 * GetDelta2D(Int_t type, Int_t var){
  TH2 * his = cosmicScan->fHistoDelta[type]->Projection(0,var);
  his->SetXTitle(cosmicScan->fHistoDelta[type]->GetAxis(var)->GetName());
  his->SetYTitle(cosmicScan->fHistoDelta[type]->GetAxis(0)->GetName());
  return his;
}


TH1* GetFit2D(Int_t type, Int_t var, Bool_t resol){
  
  TH2 * his = cosmicScan->fHistoDelta[type]->Projection(0,var);
  his->SetXTitle(cosmicScan->fHistoDelta[type]->GetAxis(var)->GetName());
  his->SetYTitle(cosmicScan->fHistoDelta[type]->GetAxis(0)->GetName());
  his->FitSlicesY(0,0,-1,0,"QNR",&fitArr);
  TH1 * hres = 0;
  if (resol) hres = (TH1*)(fitArr.At(2)->Clone());
  if (!resol) hres = (TH1*)(fitArr.At(1)->Clone());
  hres->SetMarkerStyle(20);
  hres->SetMarkerColor(2);
  hres->GetYaxis()->SetTitleOffset(1.8);
  hres->GetYaxis()->SetDecimals(kTRUE);
  return hres;
}


TH1 * GetDelta(Int_t type){
  TH1 * his = cosmicScan->fHistoDelta[type]->Projection(0);
  his->SetXTitle(cosmicScan->fHistoDelta[type]->GetAxis(0)->GetName());
  return his;
}

TH2 * GetPull2D(Int_t type, Int_t var){
  TH2 * his = cosmicScan->fHistoPull[type]->Projection(0,var);
  his->SetXTitle(cosmicScan->fHistoPull[type]->GetAxis(var)->GetName());
  his->SetYTitle(cosmicScan->fHistoPull[type]->GetAxis(0)->GetName());
  return his;
}

TH1* GetPull2DSigma(Int_t type, Int_t var){
  
  TH2 * his = cosmicScan->fHistoPull[type]->Projection(0,var);
  his->SetXTitle(cosmicScan->fHistoPull[type]->GetAxis(var)->GetName());
  his->SetYTitle(cosmicScan->fHistoPull[type]->GetAxis(0)->GetName());
  his->FitSlicesY(0,0,-1,0,"QNR",&fitArr);
  TH1 * hres = (TH1*)(fitArr.At(2)->Clone());
  return hres;
}



TH1 * GetPull(Int_t type){
  TH1 * his = cosmicScan->fHistoPull[type]->Projection(0);
  his->SetXTitle(cosmicScan->fHistoPull[type]->GetAxis(0)->GetName());
  return his;
}


void DrawResoldEdx(){
  //
  //
  //
  Int_t kmicolors[10]={1,2,3,6,7,8,9,10,11,12};
  Int_t kmimarkers[10]={21,22,23,24,25,26,27,28,29,30};
  TH2 *htemp;
  TObjArray arr;
  TH1 * hResolMax[4];
  TH1 * hResolTot[4];
  //  
  for (Int_t ipad=0;ipad<4;ipad++){
    cosmicScan->fHistodEdxTot[ipad]->GetAxis(4)->SetRangeUser(-0.6,0.6);
    cosmicScan->fHistodEdxMax[ipad]->GetAxis(4)->SetRangeUser(-0.6,0.6);
  }
  cosmicScan->fHistodEdxTot[0]->GetAxis(1)->SetRangeUser(30,62);
  cosmicScan->fHistodEdxTot[1]->GetAxis(1)->SetRangeUser(30,62);
  cosmicScan->fHistodEdxTot[2]->GetAxis(1)->SetRangeUser(10,35);
  cosmicScan->fHistodEdxTot[3]->GetAxis(1)->SetRangeUser(10,150);
  cosmicScan->fHistodEdxMax[0]->GetAxis(1)->SetRangeUser(30,62);
  cosmicScan->fHistodEdxMax[1]->GetAxis(1)->SetRangeUser(30,62);
  cosmicScan->fHistodEdxMax[2]->GetAxis(1)->SetRangeUser(10,35);
  cosmicScan->fHistodEdxMax[3]->GetAxis(1)->SetRangeUser(10,150);
  //

  for (Int_t ipad=0;ipad<4;ipad++){
    htemp = cosmicScan->fHistodEdxTot[ipad]->Projection(0,1);
    htemp->FitSlicesY(0,0,-1,0,"QNR",&arr);
    hResolTot[ipad] = (TH1*)(arr.At(2)->Clone());
    delete htemp;
    arr.SetOwner(kTRUE);
    arr.Delete();
    hResolTot[ipad]->Scale(1./TMath::Sqrt(2.));
    //
    htemp = cosmicScan->fHistodEdxMax[ipad]->Projection(0,1);
    htemp->FitSlicesY(0,0,-1,0,"QNR",&arr);
    hResolMax[ipad] = (TH1*)(arr.At(2)->Clone());
    delete htemp;
    arr.SetOwner(kTRUE);
    arr.Delete();
    hResolMax[ipad]->Scale(1./TMath::Sqrt(2.));    
  }
  hResolTot[3]->GetXaxis()->SetRangeUser(0,160);
  hResolTot[3]->SetXTitle("N_{cl}");
  hResolTot[3]->SetYTitle("#sigma(dEdx/dEdx_{d})/#sqrt{2.}");
  hResolTot[3]->SetTitle("Relative dEdx resolution");
  for (Int_t ipad=3;ipad>=0;ipad--){
    hResolTot[ipad]->SetMaximum(0.1);
    hResolTot[ipad]->SetMinimum(0.);
    hResolTot[ipad]->SetMarkerColor(kmicolors[ipad]+0);
    hResolTot[ipad]->SetMarkerStyle(kmimarkers[ipad]+1);
    if (ipad==3)    hResolTot[ipad]->Draw();
    hResolTot[ipad]->Draw("same");
    //
    hResolMax[ipad]->SetMaximum(0.1);
    hResolMax[ipad]->SetMinimum(0.);
    hResolMax[ipad]->SetMarkerColor(kmicolors[ipad]+0);
    hResolMax[ipad]->SetMarkerStyle(kmimarkers[ipad]+4);
    hResolMax[ipad]->Draw("same");
  }
  
}

void DrawStat(Int_t coord, TObjArray *array=0){
  //
  //
  //
  TCanvas *cStat = new TCanvas(Form("Cosmic stat%d",coord), Form("CosmicStat%d",coord),1000,800);
  Float_t mx0=0.2, mx1=0.05, my0=0.15, my1=0.1;
  cStat->SetMargin(mx0,mx1,my0,my1);
  cStat->Divide(3,3);
  for (Int_t i=1; i<8; i++){
    cStat->cd(i+1);
    cosmicScan->fHistoDelta[0]->Projection(i)->Draw();
  }
  if (array) array->AddLast(cStat);
}

void SetStylePad(TVirtualPad *pad){
  Float_t mx0=0.2, mx1=0.05, my0=0.15, my1=0.1;
  pad->SetMargin(mx0,mx1,my0,my1);
  pad->SetTicks(1,1);
  pad->SetGrid(1,1); 
  
}

void MakePlotPt(TObjArray * array){
  //
  //
  TCanvas *cptRes = new TCanvas("TPCPtResol","TPCPtResol",900,600);
  cptRes->Divide(2,1);
  SetStylePad(cptRes->cd(1));
  SetStylePad(cptRes->cd(2)); 
  cptRes->cd(1);
  TLegend *legend = new TLegend(0.2,0.7,0.6,0.9,"");

  //
  TH1 * hisRes=0;
  TH1 * hisMean=0;
  for (Int_t i=0; i<3; i++){
    if (i==0) cosmicScan->fHistoDelta[5]->GetAxis(6)->SetRangeUser(-1,1);  
    if (i==1) cosmicScan->fHistoDelta[5]->GetAxis(6)->SetRangeUser(0.001,1);
    if (i==2) cosmicScan->fHistoDelta[5]->GetAxis(6)->SetRangeUser(-1,-0.001);
    hisRes  = (TH1*)GetFit2D(5,7,kTRUE)->Clone();
    hisMean = (TH1*)GetFit2D(5,7,kFALSE)->Clone();
    hisMean->SetName(Form("#Deltap_{t}/p_{t} %s",chsign[i]));
    hisRes->SetName(Form("#sigma_{p_{t}}/p_{t} %s",chsign[i]));
    hisMean->SetTitle(Form("#Delta_{p_{t}}/p_{t} %s",chsign[i]));
    hisRes->SetTitle(Form("#sigma_{p_{t}}/p_{t} %s",chsign[i]));

    hisRes->SetMarkerStyle(20);
    hisMean->SetMarkerStyle(20);
    hisRes->SetMarkerColor(colors[i]);
    hisMean->SetMarkerColor(colors[i]);
    hisRes->Scale(100);
    hisMean->Scale(100);
    hisRes->SetMaximum(30);
    hisRes->SetMinimum(0);
    hisMean->SetMaximum(20);
    hisMean->SetMinimum(-20);
    hisRes->SetYTitle("#sigma_{p_{t}}/p_{t} (%)");
    hisMean->SetYTitle("#Delta_{p_{t}}/p_{t} (%)");
    hisRes->GetXaxis()->SetRangeUser(0,50);
    hisMean->GetXaxis()->SetRangeUser(0,50);  
    cptRes->cd(2); 
    hisRes->Draw("same");
    if (i==0) hisRes->Draw("");
    cptRes->cd(1);
    hisMean->Draw("same");
    if (i==0) hisMean->Draw("");
    legend->AddEntry(hisMean);
  }
  legend->Draw();
  if (array) array->AddLast(cptRes);
}


void MakePlotP4(TObjArray * array){
  //
  //
  TCanvas *cptRes = new TCanvas("TPCP4Resol","TPCP4Resol",900,600);
  cptRes->Divide(2,1);
  SetStylePad(cptRes->cd(1));
  SetStylePad(cptRes->cd(2));
  cptRes->cd(1);
  TLegend *legend = new TLegend(0.2,0.7,0.6,0.9,"");

  //
  TH1 * hisRes  =0;
  TH1 * hisMean=0;

  for (Int_t i=0; i<3; i++){
    if (i==0) cosmicScan->fHistoDelta[4]->GetAxis(6)->SetRangeUser(-1,1);
    if (i==1) cosmicScan->fHistoDelta[4]->GetAxis(6)->SetRangeUser(0.001,1);
    if (i==2) cosmicScan->fHistoDelta[4]->GetAxis(6)->SetRangeUser(-1,-0.001);
    hisRes  = (TH1*)GetFit2D(4,7,kTRUE)->Clone();
    hisMean = (TH1*)GetFit2D(4,7,kFALSE)->Clone();
    hisMean->SetName(Form("#Delta_{1/p_{t}} %s",chsign[i]));
    hisRes->SetName(Form("#sigma_{1/p_{t}} %s",chsign[i]));
    hisMean->SetTitle(Form("#Delta_{1/p_{t}} %s",chsign[i]));
    hisRes->SetTitle(Form("#sigma_{1/p_{t}} %s",chsign[i]));

    hisRes->SetMarkerStyle(20+i);
    hisMean->SetMarkerStyle(20+i);
    hisMean->SetMarkerColor(colors[i]);
    hisRes->SetMarkerColor(colors[i]);
    hisRes->SetMaximum(0.02);
    hisRes->SetMinimum(-0.0);
    hisMean->SetMaximum(0.02);
    hisMean->SetMinimum(-0.02);
    hisRes->SetYTitle("#sigma_{1/pt} (1/GeV)");
    hisMean->SetYTitle("#Delta_{1/pt} (1/GeV)");
    hisMean->SetXTitle("p_{t} (GeV)");
    hisRes->SetXTitle("p_{t} (GeV)");
    hisRes->GetXaxis()->SetRangeUser(0,50);
    hisMean->GetXaxis()->SetRangeUser(0,50);      
    cptRes->cd(2); 
    hisRes->Draw("same");
    if (i==0) hisRes->Draw("");
    cptRes->cd(1);
    hisMean->Draw("same");
    if (i==0) hisMean->Draw("");
    legend->AddEntry(hisMean);
  }
  legend->Draw();
  if (array) array->AddLast(cptRes);
}






void MakePlotPosY(TObjArray * array){
  //
  //
  TCanvas *cptRes = new TCanvas("TPCPosResolY","TPCPosResolY",900,600);
  cptRes->Divide(2,1);
  SetStylePad(cptRes->cd(1));
  SetStylePad(cptRes->cd(2));
  cptRes->cd(1);
  TLegend *legend = new TLegend(0.2,0.7,0.6,0.9,"");

  //
  TH1 * hisRes=0;
  TH1 * hisMean=0;
  for (Int_t i=0; i<3; i++){
    if (i==1) cosmicScan->fHistoDelta[0]->GetAxis(6)->SetRangeUser(0.0001,1);
    if (i==2) cosmicScan->fHistoDelta[0]->GetAxis(6)->SetRangeUser(-1,-0.0001);
    if (i==0) cosmicScan->fHistoDelta[0]->GetAxis(6)->SetRangeUser(-1,1);
    hisRes  = (TH1*)GetFit2D(0,7,kTRUE)->Clone();
    hisMean = (TH1*)GetFit2D(0,7,kFALSE)->Clone();
    hisMean->SetName(Form("#Delta_{r#phi} %s",chsign[i]));
    hisRes->SetName(Form("#sigma_{r#phi} %s",chsign[i]));
    hisMean->SetTitle(Form("#Delta_{r#phi} %s",chsign[i]));
    hisRes->SetTitle(Form("#sigma_{r#phi} %s",chsign[i]));

    hisRes->SetMarkerStyle(20+i);
    hisMean->SetMarkerStyle(20+i);
    hisMean->SetMarkerColor(colors[i]);
    hisRes->SetMarkerColor(colors[i]);
    
    //
    hisRes->SetMaximum(0.4);
    hisRes->SetMinimum(0.0);
    hisMean->SetMaximum(0.4);
    hisMean->SetMinimum(-0.4);
    hisRes->SetName("Y resol");
    hisRes->SetName("Y resolution");
    hisRes->SetYTitle("#sigma_{y} (cm)");
    hisMean->SetYTitle("#Delta_{y} (cm)");
    hisRes->GetXaxis()->SetRangeUser(0,50);
    hisMean->GetXaxis()->SetRangeUser(0,50);  
    cptRes->cd(2); 
    hisRes->Draw("same");
    if (i==0) hisRes->Draw("");
    cptRes->cd(1);
    hisMean->Draw("same");
    if (i==0) hisMean->Draw("");
    legend->AddEntry(hisMean);   
  }
  legend->Draw();
  if (array) array->AddLast(cptRes);
}

void MakePlotSnp(TObjArray * array){
  //
  //
  TCanvas *cptRes = new TCanvas("TPCSnp","TPCSnp",900,600);
  cptRes->Divide(2,1);
  SetStylePad(cptRes->cd(1));
  SetStylePad(cptRes->cd(2));
  cptRes->cd(1);
  TLegend *legend = new TLegend(0.2,0.7,0.6,0.9,"");

  //
  TH1 * hisRes=0;
  TH1 * hisMean=0;
  for (Int_t i=0; i<3; i++){
    if (i==1) cosmicScan->fHistoDelta[2]->GetAxis(6)->SetRangeUser(0.0001,1);
    if (i==2) cosmicScan->fHistoDelta[2]->GetAxis(6)->SetRangeUser(-1,-0.0001);
    if (i==0) cosmicScan->fHistoDelta[2]->GetAxis(6)->SetRangeUser(-1,1);
    hisRes  = (TH1*)GetFit2D(2,7,kTRUE)->Clone();
    hisMean = (TH1*)GetFit2D(2,7,kFALSE)->Clone();
    hisMean->SetName(Form("#Delta_{#phi} %s",chsign[i]));
    hisRes->SetName(Form("#sigma_{#phi} %s",chsign[i]));
    hisMean->SetTitle(Form("#Delta_{#phi} %s",chsign[i]));
    hisRes->SetTitle(Form("#sigma_{#phi} %s",chsign[i]));

    hisRes->SetMarkerStyle(20+i);
    hisMean->SetMarkerStyle(20+i);
    hisMean->SetMarkerColor(colors[i]);
    hisRes->SetMarkerColor(colors[i]);
    
    //
    hisRes->SetMaximum(4);
    hisRes->SetMinimum(4);
    hisMean->SetMaximum(4);
    hisMean->SetMinimum(-4);
    hisRes->SetYTitle("#sigma_{#phi} (mrad)");
    hisMean->SetYTitle("#Delta_{#phi} (mrad)");
    hisRes->GetXaxis()->SetRangeUser(0,50);
    hisMean->GetXaxis()->SetRangeUser(0,50);  
    cptRes->cd(2); 
    hisRes->Draw("same");
    if (i==0) hisRes->Draw("");
    cptRes->cd(1);
    hisMean->Draw("same");
    if (i==0) hisMean->Draw("");
    legend->AddEntry(hisMean);   
  }
  legend->Draw();
  if (array) array->AddLast(cptRes);
}
void MakePlotTgl(TObjArray * array){
  //
  //
  TCanvas *cptRes = new TCanvas("TPCtgl","TPCtgl",900,600);
  cptRes->Divide(2,1);
  SetStylePad(cptRes->cd(1));
  SetStylePad(cptRes->cd(2));
  cptRes->cd(1);
  TLegend *legend = new TLegend(0.2,0.7,0.6,0.9,"");

  //
  TH1 * hisRes=0;
  TH1 * hisMean=0;
  for (Int_t i=0; i<3; i++){
    if (i==1) cosmicScan->fHistoDelta[3]->GetAxis(6)->SetRangeUser(0.0001,1);
    if (i==2) cosmicScan->fHistoDelta[3]->GetAxis(6)->SetRangeUser(-1,-0.0001);
    if (i==0) cosmicScan->fHistoDelta[3]->GetAxis(6)->SetRangeUser(-1,1);
    hisRes  = (TH1*)GetFit2D(3,7,kTRUE)->Clone();
    hisMean = (TH1*)GetFit2D(3,7,kFALSE)->Clone();
    hisMean->SetName(Form("#Delta_{#theta} %s",chsign[i]));
    hisRes->SetName(Form("#sigma_{#theta} %s",chsign[i]));
    hisMean->SetTitle(Form("#Delta_{#theta} %s",chsign[i]));
    hisRes->SetTitle(Form("#sigma_{#theta} %s",chsign[i]));

    hisRes->SetMarkerStyle(20+i);
    hisMean->SetMarkerStyle(20+i);
    hisMean->SetMarkerColor(colors[i]);
    hisRes->SetMarkerColor(colors[i]);
    
    //
    hisRes->SetMaximum(4);
    hisRes->SetMinimum(4);
    hisMean->SetMaximum(4);
    hisMean->SetMinimum(-4);
    hisRes->SetYTitle("#sigma_{#theta} (mrad)");
    hisMean->SetYTitle("#Delta_{#theta} (mrad)");
    hisRes->GetXaxis()->SetRangeUser(0,50);
    hisMean->GetXaxis()->SetRangeUser(0,50);  
    cptRes->cd(2); 
    hisRes->Draw("same");
    if (i==0) hisRes->Draw("");
    cptRes->cd(1);
    hisMean->Draw("same");
    if (i==0) hisMean->Draw("");
    legend->AddEntry(hisMean);   
  }
  legend->Draw();
  if (array) array->AddLast(cptRes);
}

void MakePlotPosZ(TObjArray * array){
  //
  //
  TCanvas *cptRes = new TCanvas("TPCPosResolZ","TPCPosResolZ",900,600);
  cptRes->Divide(2,1);
  SetStylePad(cptRes->cd(1));
  SetStylePad(cptRes->cd(2));
  cptRes->cd(1);
  TLegend *legend = new TLegend(0.2,0.7,0.6,0.9,"");

  //
  TH1 * hisRes=0;
  TH1 * hisMean=0;
  for (Int_t i=0; i<3; i++){
    if (i==1) cosmicScan->fHistoDelta[1]->GetAxis(6)->SetRangeUser(0.001,1.);
    if (i==2) cosmicScan->fHistoDelta[1]->GetAxis(6)->SetRangeUser(-1,-0.001);
    if (i==0) cosmicScan->fHistoDelta[1]->GetAxis(6)->SetRangeUser(-1,1);
    
    hisRes  = (TH1*)GetFit2D(1,7,kTRUE)->Clone();
    hisMean = (TH1*)GetFit2D(1,7,kFALSE)->Clone();   
    hisMean->SetName(Form("#Delta_{z} %s",chsign[i]));
    hisRes->SetName(Form("#sigma_{z} %s",chsign[i]));
    hisMean->SetTitle(Form("#Delta_{z} %s",chsign[i]));
    hisRes->SetTitle(Form("#sigma_{z} %s",chsign[i]));


    hisRes->SetMaximum(0.4);
    hisRes->SetMinimum(0.0);
    hisMean->SetMaximum(0.2);
    hisMean->SetMinimum(-0.2); 
    hisRes->SetMarkerStyle(20);
    hisMean->SetMarkerStyle(20);
    hisRes->SetMarkerColor(colors[i]);
    hisMean->SetMarkerColor(colors[i]);
 
    hisRes->SetName("Z resol");
    hisRes->SetName("Z resolution");
    hisRes->SetYTitle("#sigma_{z} (cm)");
    hisMean->SetYTitle("#Delta_{z} (cm)");
    hisRes->GetXaxis()->SetRangeUser(0,50);
    hisMean->GetXaxis()->SetRangeUser(0,50);  
    cptRes->cd(2); 
    hisRes->Draw("same");
    if (i==0) hisRes->Draw();
    cptRes->cd(1);
    hisMean->Draw("same");
    if (i==0) hisMean->Draw();
    legend->AddEntry(hisMean);
  }
  legend->Draw();
  if (array) array->AddLast(cptRes);
}


void  MakeDefaultPlots(){
  //
  //
  //
  gStyle->SetOptStat(1100);
  DrawStat(0,picArray);
  gStyle->SetOptStat(0);
  MakePlotPosY(picArray);
  MakePlotPosZ(picArray);
  MakePlotSnp(picArray);
  MakePlotTgl(picArray);
  MakePlotP4(picArray);
  MakePlotPt(picArray);
  //

  TFile f("cosmicPlots.root","recreate");
  picArray->Write("CosmicPlots",TObject::kSingleKey);
  f.Close();
  TPostScript *ps=new TPostScript("cosmicPerformance.ps", 112);
  ps->NewPage();
  for (Int_t ipad=0;ipad<picArray->GetEntries();ipad++){
    TCanvas *c =dynamic_cast<TCanvas*> (picArray->At(ipad));
    if (c) {
      c->Draw();
      c->Update();
      ps->NewPage();
    }
  } 
  ps->Close();
  delete ps;
}
