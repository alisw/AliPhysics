//**************************************************************************//
//WARNING: RUN AFTER COMPILING WITH .L++ TO AVOID SEGMENTATION FALUT CRASHES!!
//**************************************************************************//

#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TF1.h>
#include <Riostream.h>
#include <TBufferFile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TPaveLabel.h>
#include <TVirtualPad.h>
#include <TMath.h>
#include <TLatex.h>
#include <TColor.h>
#include <TClass.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <sstream>

TString strPtAss[3]={"0.3to1.0","1.0to2.0","2.0to3.0"};
TString strPtDCanvas[4]={"3 < #it{p}_{T}^{D} < 5 GeV/#it{c}","5 < #it{p}_{T}^{D} < 8 GeV/#it{c}","8 < #it{p}_{T}^{D} < 16 GeV/#it{c}","16 < #it{p}_{T}^{D} < 24 GeV/#it{c}"};
TString strSystem[2]={"pp, 5 TeV","p-Pb, 5 TeV"};
Color_t colSystem[2]={kBlue,kRed};
Int_t markerStyle[2]={20,21};
Bool_t useLegendForData=kTRUE;
Bool_t plotv2unc=kFALSE;
TString strFitResultPPb[2]={"./pp_5TeV/FitResults/Trends_pp","./pPb_5TeV/FitResults/Trends_pPb"}; //  "/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2015June7finalPlots/ReflectedPlots/StdRebin/AllPlots/Averages/FitResults";
Double_t canvasheight=1000;
Double_t resizeTextFactor=1.3;//1.*canvasheight/1200.; // size was tuned for canvasheight =1800. 
Double_t referencePadHeight=0.48; // Do not touch unless the canvas division is changed from 2x3, change canvasheight and resizeTextFactor instead
Double_t innerPadHeight;// not touch, set internally
Double_t innerPadWidth;// not touch, set internally
Int_t style=1;
Double_t scaleHeightPads=1;// do not touch this, it is regulated automatically in the macro
Double_t scaleWidthPads=1;// do not touch this, it is regulated automatically in the macro
TString strFitResultMC[4]={"","","",""};
TString basicdir="$PWD";
//TString strFitResultMC[2]={"/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2015June7finalPlots/MCTemplates/Templates_pp_12May15/FitResults/",
//			   "/Users/administrator/ALICE/CHARM/HFCJ/DCorrelations_Test/2015May19UseScriptPWGHF/MCTemplates/Templates_pPb_12May15/FitResults/"};

TString strquantityFile[5]={"NSYield","NSSigma","Pedestal","ASYield","ASSigma"};

Double_t maxRangePPb[6][5]={
         {1.85,0.64,13,1.95,1.09},
         {1.85,0.64,13,1.95,1.09},
         {1.85,0.64,13,1.95,1.09},
         {1.35,0.64,13,2,1.28},
         {1.35,0.64,13,2,1.28},
         {1.35,0.64,13,2,1.28}};

Double_t maxRangePPbSinglePanel[6][5]={
         {3.2,0.8,9.0,4.3,1.48},
         {1.6,0.8,6.0,2.5,1.48},
         {2.1,0.8,3.0,2.5,1.48},
         {1.5,0.8,2.0,1.5,1.48},
         {0.6,0.8,0.5,1.0,1.48},
         {0.6,0.8,0.3,1.0,1.48}};

Double_t minptMC=0.,maxptMC=99.;
Double_t minptData=0.5,maxptData=3.5; //you have to be tighter than real limits...
Int_t ncolumns=4;
Int_t nrows=2;
Double_t ytitleoffset=3,xtitleoffset=2.4;
Bool_t skip3to5=kFALSE;
Double_t markersize=1.2;
Double_t markersizeMC=1.2;
Bool_t drawSystMC=kTRUE;
Bool_t runonpPb2016=kTRUE;

void SetMinPtDisplayData(Double_t ptmin){minptMC=ptmin;}
void SetMaxPtDisplayData(Double_t ptmax){maxptMC=ptmax;}
void SetMinPtDisplayMC(Double_t ptmin){minptData=ptmin;}
void SetMaxPtDisplayMC(Double_t ptmax){maxptData=ptmax;}
void SetRunOn2016(Bool_t run){runonpPb2016=run;}

TCanvas* Compare(Int_t binass,Int_t quantity,TPad *pd=0x0,Int_t textlegendOptions=0);

void AdaptRangeHist(TH1D *h,Double_t minpt,Double_t maxpt){
  for(Int_t i=1;i<=h->GetNbinsX();i++){
    if(h->GetXaxis()->GetBinUpEdge(i)<minpt*1.001 ||  h->GetXaxis()->GetBinLowEdge(i)>maxpt*0.999 ){
      h->SetBinContent(i,-1);
      h->SetBinError(i,0);
    }
  }
}

void AdaptRangeTGraph(TGraph *gr,Double_t minpt,Double_t maxpt){
  Double_t *x=gr->GetX();
  Double_t *exl=gr->GetEXlow();
  Double_t *exh=gr->GetEXhigh();
  Int_t np=gr->GetN();
  for(Int_t i=0;i<np;i++){
    if(x[i]+exh[i]<minpt*1.0001){
      gr->RemovePoint(i);
      AdaptRangeTGraph(gr,minpt,maxpt);
      return;
    }
    if(x[i]-exl[i]>maxpt*0.999){
      gr->RemovePoint(i);
      AdaptRangeTGraph(gr,minpt,maxpt);
      return;
    }
  }
  return;
}

void SetDrawSystMC(Bool_t drawsystmc){
  drawSystMC=drawsystmc;
}
void Init3x3Settings(){
  referencePadHeight=0.35;
  resizeTextFactor=0.85;
  ytitleoffset=4.8;
  xtitleoffset=4;
  ncolumns=3;
  nrows=3;
  markersize=1.;
  markersizeMC=1.;
}

TString yaxisTitle[5]={"Associated yield","Peak width (rad)","Baseline (rad^{-1})","Associated yield","Peak width (rad)"};
Double_t leftMarginCanvas=0.17;
Double_t rightMarginCanvas=0.055;
Double_t bottomMarginCanvas=0.13;
Double_t topMarginCanvas=0.07;
const Int_t nmodels=8;
Bool_t includemodel[nmodels]={kTRUE,kTRUE,kTRUE,kTRUE,kFALSE,kFALSE,kTRUE,kFALSE};
TString strModelDir[nmodels]={"Perugia0","Perugia2010","Perugia2011","PYTHIA8","HERWIG","POWHEG","POWHEG","EPOS3"};
TString strModelDirLeg[nmodels]={"PYTHIA6, Perugia 0","PYTHIA6, Perugia 2010","PYTHIA6, Perugia 2011","PYTHIA8, Tune 4C","HERWIG","POWHEG+PYTHIA6","POWHEG+PYTHIA6 EPS09","EPOS 3.117"};
Color_t modelColors[nmodels]={kRed+2,kCyan,kGreen+2,kMagenta+1,kViolet,kBlue,kBlue,kOrange+1};
Bool_t includeinlegend[nmodels]={kTRUE,kTRUE,kTRUE,kTRUE,kFALSE,kFALSE,kTRUE,kFALSE};// this is also used to split the legend in 2!!
Int_t modelMarkerStyle[nmodels]={kOpenSquare,kOpenCircle,kOpenDiamond,28,26,3,3,33};

TH1D **hMC;
TGraphAsymmErrors **grMC;

void SetSkip3to5(Bool_t skip){
  skip3to5=skip;
}
void IncludePerugia0(Bool_t incl){
  includemodel[0]=incl;
  includeinlegend[0]=incl;
}
void IncludePerugia2010(Bool_t incl){
  includemodel[1]=incl;
  includeinlegend[1]=incl;
}
void IncludePerugia2011(Bool_t incl){
  includemodel[2]=incl;
  includeinlegend[2]=incl;
}
void IncludePythia8(Bool_t incl){
  includemodel[3]=incl;
  includeinlegend[3]=incl;
}
void IncludePowheg(Bool_t incl){
  includemodel[4]=incl;
  includeinlegend[4]=incl;
}
void IncludeEPOS(Bool_t incl){
  includemodel[5]=incl;
  includeinlegend[5]=incl;
}
void IncludeModel(Int_t imod,Bool_t incl){
  includemodel[imod]=incl;
  includeinlegend[imod]=incl;
}
Int_t GetNPythia6Tunes(){
  Int_t ncount=0;
  for(Int_t j=0;j<3;j++){
    if(includemodel[j])ncount++;
  }
  return ncount;
}

Int_t CountNmodels(){
  Int_t ncount=0;
  for(Int_t j=0;j<nmodels;j++){
    if(includemodel[j])ncount++;
  }
  return ncount;
}

Int_t CountNmodelsInLegend(){
  Int_t ncount=0;
  for(Int_t j=0;j<nmodels;j++){
    if(includeinlegend[j])ncount++;
    if(j==6&&includeinlegend[j])ncount++;
  }
  return ncount;
}

void SetBasicDir(TString dir){
  basicdir=dir;
}

void SetDirectoryFitResult(TString str1,TString str2){
  strFitResultPPb[0]=str1;
  strFitResultPPb[1]=str2;
}

void SetDirectoryFitResultsMCpp(TString str){
  strFitResultMC[0]=str;
}

void SetDirectoryFitResultsMCpPb(TString str){
  strFitResultMC[1]=str;
}

TCanvas *CreateCanvasWithDefaultStyle(TString name){
  gStyle->SetOptStat(0000);
  TString nameNoPoints=name;
  nameNoPoints.ReplaceAll(".","");
  TCanvas *cout=new TCanvas(nameNoPoints.Data(),name.Data(),800,800);
  cout->SetTicky();
  cout->SetTickx();
  cout->SetFrameBorderMode(0);
  cout->SetLeftMargin(leftMarginCanvas);
  cout->SetRightMargin(rightMarginCanvas);
  cout->SetBottomMargin(bottomMarginCanvas);
  cout->SetTopMargin(topMarginCanvas);

  return cout;
}

void ConvertTH1ToTGraphAsymmError(TH1D* h,TGraphAsymmErrors *&gr, Double_t shift) {
  const Int_t nbinsxx=4;
  Double_t x[nbinsxx], y[nbinsxx], ex1[nbinsxx], ex2[nbinsxx], ey1[nbinsxx], ey2[nbinsxx];

  for(int i=0; i<nbinsxx; i++) {
    x[i] = h->GetBinCenter(i+2)+shift;
    y[i] = h->GetBinContent(i+2);
    ex1[i] = h->GetBinCenter(i+2)-h->GetBinLowEdge(i+2)+shift;
    ex2[i] = ex1[i]-2*shift;
    ey1[i] = h->GetBinError(i+2);
    ey2[i] = h->GetBinError(i+2);
  }

  gr = new TGraphAsymmErrors(nbinsxx,x,y,ex1,ex2,ey1,ey2);
  return;

}

void ConvertTH1ToTGraphAsymmError2016(TH1D* h,TGraphAsymmErrors *&gr, Double_t shift, Int_t system) {
 
  const Int_t nbinsxx=4;
  Double_t x[nbinsxx], y[nbinsxx], ex1[nbinsxx], ex2[nbinsxx], ey1[nbinsxx], ey2[nbinsxx];
 
 if(system==0) {
  for(int i=0 ; i<nbinsxx; i++) {
    x[i] = h->GetBinCenter(i+3)+shift;
    y[i] = h->GetBinContent(i+3);
    ex1[i] = h->GetBinCenter(i+3)-h->GetBinLowEdge(i+3)+shift;
    ex2[i] = ex1[i]-2*shift;
    ey1[i] = h->GetBinError(i+3);
    ey2[i] = h->GetBinError(i+3);
  }
}
if(system==1) {
  for(int i=0; i<nbinsxx; i++) {
    x[i] = h->GetBinCenter(i+2)+shift;
    y[i] = h->GetBinContent(i+2);
    ex1[i] = h->GetBinCenter(i+2)-h->GetBinLowEdge(i+2)+shift;
    ex2[i] = ex1[i]-2*shift;
    ey1[i] = h->GetBinError(i+2);
    ey2[i] = h->GetBinError(i+2);
  }  
}  

  gr = new TGraphAsymmErrors(nbinsxx,x,y,ex1,ex2,ey1,ey2);
  return;

}

void SetPadStyle(TPad *p){
  p->SetTicky();
  p->SetTickx();
  p->SetFrameBorderMode(0);
  p->SetFillStyle(0);
  p->SetFrameFillStyle(4000);
  //  p->SetBottomMargin(0);
  //  p->SetTopMargin(0);
  //  p->Range(-2.6,-1.5,5.4,4.4);
  //  p->SetLeftMargin(leftMarginCanvas);
  //  p->SetRightMargin(rightMarginCanvas);
  //  p->SetBottomMargin(bottomMarginCanvas);
  //  p->SetTopMargin(topMarginCanvas);
}

void Set4x2PadPositions(TCanvas* c){
    
  TPad * pd1 =(TPad*) c->GetPad(1);
  TPad * pd2 =(TPad*) c->GetPad(2);
  TPad * pd3 =(TPad*) c->GetPad(3);
  TPad * pd4 =(TPad*) c->GetPad(4);
  TPad * pd5 =(TPad*) c->GetPad(5);
  TPad * pd6 =(TPad*) c->GetPad(6);
  TPad * pd7 =(TPad*) c->GetPad(7);
  TPad * pd8 =(TPad*) c->GetPad(8);

  SetPadStyle(pd1);
  SetPadStyle(pd2);
  SetPadStyle(pd3);
  SetPadStyle(pd4);
  SetPadStyle(pd5);
  SetPadStyle(pd6);
  SetPadStyle(pd7);
  SetPadStyle(pd8);

  Double_t xl,xu,yl,yu;
  Double_t marginLeft=0.1;
  Double_t marginRight=0.02;
  Double_t marginTop=0.04;
  Double_t marginBottom=0.08;
  innerPadWidth=(1-marginLeft-marginRight)/4.;// this is the width w/o margin, not the real pad width!!
  innerPadHeight=(1-marginTop-marginBottom)/2.;// this is the height w/o margin, not the real pad height, which differs between inner pads and pads at the "boarders"!!
  Printf("innerPadHeight: %f",innerPadHeight);
  Printf("innerPadWidth: %f",innerPadWidth);
  Double_t marginLeftForXAxis=0.02;
  Double_t marginBottomForYAxis=0.02;

 // Bottom row

    pd5->GetPadPar(xl,yl,xu,yu);
    Printf("PAD 7 Original values: xl %f  xu %f  yl  %f  yu %f, gPad->GetwNDC=%f, gPad->GetHNDC=%f",xl,xu,yl,yu,gPad->GetWNDC(),gPad->GetHNDC());
    pd5->SetPad(0.,0.,innerPadWidth+marginLeft,innerPadHeight+marginBottom);
    pd5->GetPadPar(xl,yl,xu,yu);
    Printf("New values: xl %f  xu %f  yl  %f  yu %f, gPad->GetwNDC=%f, gPad->GetHNDC=%f",xl,xu,yl,yu,gPad->GetWNDC(),gPad->GetHNDC());
 
    pd5->SetLeftMargin(marginLeft/(marginLeft+innerPadWidth));
    pd5->SetRightMargin(0);
    pd5->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
    pd5->SetTopMargin(0.);
    pd5->Modified();
    pd5->SetFillStyle(0);

    pd6->GetPadPar(xl,yl,xu,yu);
    Printf("PAD 10 Original values: xl %f  xu %f  yl  %f  yu %f, gPad->GetwNDC=%f, gPad->GetHNDC=%f",xl,xu,yl,yu,gPad->GetWNDC(),gPad->GetHNDC());
    pd6->SetPad(innerPadWidth+marginLeft-marginLeftForXAxis,0,2*innerPadWidth+marginLeft,innerPadHeight+marginBottom);
    pd6->GetPadPar(xl,yl,xu,yu);
    Printf("New values: xl %f  xu %f  yl  %f  yu %f, gPad->GetwNDC=%f, gPad->GetHNDC=%f",xl,xu,yl,yu,gPad->GetWNDC(),gPad->GetHNDC());
 
    pd6->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis));
    //pd8->SetLeftMargin(0.004/(1.-innerPadWidth-marginLeft+0.004));
    pd6->SetRightMargin(0.);
    pd6->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
    pd6->SetTopMargin(0);
    pd6->SetFillStyle(0);

    pd6->Modified();
    pd6->Update();

    pd7->GetPadPar(xl,yl,xu,yu);
    Printf("PAD 11 Original values: xl %f  xu %f  yl  %f  yu %f, gPad->GetwNDC=%f, gPad->GetHNDC=%f",xl,xu,yl,yu,gPad->GetWNDC(),gPad->GetHNDC());
  
    pd7->SetPad(2*innerPadWidth+marginLeft-marginLeftForXAxis,0,3*innerPadWidth+marginLeft,innerPadHeight+marginBottom);
    pd7->GetPadPar(xl,yl,xu,yu);
    Printf("New values: xl %f  xu %f  yl  %f  yu %f, gPad->GetwNDC=%f, gPad->GetHNDC=%f",xl,xu,yl,yu,gPad->GetWNDC(),gPad->GetHNDC());
    pd7->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis));
    pd7->SetRightMargin(0);
    pd7->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
    pd7->SetTopMargin(0.);
    pd7->SetFillStyle(0);

    pd7->Modified();
    pd7->Update();

    pd8->GetPadPar(xl,yl,xu,yu);
    Printf("PAD 12 Original values: xl %f  xu %f  yl  %f  yu %f, gPad->GetwNDC=%f, gPad->GetHNDC=%f",xl,xu,yl,yu,gPad->GetWNDC(),gPad->GetHNDC());
    pd8->SetPad(3*innerPadWidth+marginLeft-marginLeftForXAxis,0,1.,innerPadHeight+marginBottom);
    pd8->GetPadPar(xl,yl,xu,yu);
    Printf("New values: xl %f  xu %f  yl  %f  yu %f, gPad->GetwNDC=%f, gPad->GetHNDC=%f",xl,xu,yl,yu,gPad->GetWNDC(),gPad->GetHNDC());
 
    pd8->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis+marginRight));
    pd8->SetRightMargin(marginRight/(innerPadWidth+marginRight));
    pd8->SetBottomMargin(marginBottom/(marginBottom+innerPadHeight));
    pd8->SetTopMargin(0.);
    pd8->SetFillStyle(0);

    pd8->Modified();
    pd8->Update();

    // Top Row
    pd1->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd1->SetPad(0.,innerPadHeight+marginBottom-marginBottomForYAxis,innerPadWidth+marginLeft,2.*innerPadHeight+marginBottom);
    pd1->SetLeftMargin(marginLeft/(marginLeft+innerPadWidth));
    pd1->SetRightMargin(0.);
    pd1->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis));
    pd1->SetTopMargin(0.);
    pd1->SetFillStyle(0);

    pd1->Modified();
    pd1->Update();

    pd2->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd2->SetPad(innerPadWidth+marginLeft-marginLeftForXAxis,innerPadHeight+marginBottom-marginBottomForYAxis,2*innerPadWidth+marginLeft,2.*innerPadHeight+marginBottom);
    pd2->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis));
    pd2->SetRightMargin(0.);
    pd2->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis));
    pd2->SetTopMargin(0.);
    pd2->SetFillStyle(0);

    pd2->Modified();
    pd2->Update();

    pd3->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd3->SetPad(2*innerPadWidth+marginLeft-marginLeftForXAxis,innerPadHeight+marginBottom-marginBottomForYAxis,3*innerPadWidth+marginLeft,2.*innerPadHeight+marginBottom);
    pd3->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis));
    pd3->SetRightMargin(0.);
    pd3->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis));
    pd3->SetTopMargin(0.);

    pd3->Modified();
    pd3->Update();

    pd4->GetPadPar(xl,yl,xu,yu);
    Printf("Original values: xl %f  xu %f  yl  %f  yu %f",xl,xu,yl,yu);
    pd4->SetPad(3*innerPadWidth+marginLeft-marginLeftForXAxis,innerPadHeight+marginBottom-marginBottomForYAxis,1.,2.*innerPadHeight+marginBottom);
    pd4->SetLeftMargin(marginLeftForXAxis/(innerPadWidth+marginLeftForXAxis+marginRight));
    pd4->SetRightMargin(marginRight/(innerPadWidth+marginRight));
    pd4->SetBottomMargin(marginBottomForYAxis/(innerPadHeight+marginBottomForYAxis));
    pd4->SetTopMargin(0.);
    
    pd4->Modified();
    pd4->Update();

    scaleHeightPads=pd1->GetHNDC();
    scaleWidthPads=pd1->GetWNDC();

}

TCanvas* CompareNSyieldPPtoPPb(Int_t binass){
  return Compare(binass,0);
}

TCanvas* CompareNSsigmaPPtoPPb(Int_t binass){
  return  Compare(binass,1);

}

TCanvas* ComparePedestalPPtoPPb(Int_t binass){
  return Compare(binass,2);
}

TLatex *GetDRapForSystem(Int_t collSystem,Int_t identifier,Int_t includeDEta=0){
  
  TLatex *tlrap;
  Double_t x=0.015,y=0.28;
  if(ncolumns==3){
    x=0.015;// this and the following numbers does not make too much sense, they come just from an optimization
    if(includeDEta)x=0.060; // was 0.065
    if(collSystem==1&&includeDEta)x=0.026;// was 0.048
  }
  if(nrows==3){// these are hard coded number from an optimization
    if(includeDEta)y=0.25;
    else y=0.18;
  }
  if(nrows==2){// these are hard coded number from an optimization
    y=0.3850;// was 0.390 at round 1    
  }
  TString str;
  if(collSystem==0)str="";//"|#it{y}^{D}_{cms}| < 0.5";
  if(collSystem==1)str="";//"-0.96 < #it{y}^{D}_{cms} < 0.04";
  if(includeDEta==1){
    str.Append("|#Delta#eta| < 1");
  }
  if(style==-1){
    if(collSystem==0)tlrap=new TLatex(0.24,0.75,"|#it{y}^{D}_{cms}| < 0.5");
    else if(collSystem==1)tlrap=new TLatex(0.24,0.75,"-0.96 < #it{y}^{D}_{cms} < 0.04");
    else return 0x0;
    tlrap->SetNDC();
    tlrap->SetTextFont(42);
    tlrap->SetTextSize(0.03);
  }
  else{
    if(collSystem==0)tlrap=new TLatex(0.016/gPad->GetWNDC()+gPad->GetLeftMargin(),0.303/gPad->GetHNDC()+gPad->GetBottomMargin(),str.Data()); 
    else if(collSystem==1)tlrap=new TLatex(x/gPad->GetWNDC()+gPad->GetLeftMargin(),y/gPad->GetHNDC()+gPad->GetBottomMargin(),str.Data()); 
    tlrap->SetNDC();
    tlrap->SetTextFont(43);
    tlrap->SetTextSize(24*innerPadHeight/referencePadHeight*resizeTextFactor);//0.06/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor);
  }  
  tlrap->SetTextAlign(11);
  tlrap->SetName(Form("tlrap_%d",identifier));
  return tlrap;
}

TLatex *GetDEtaD(Int_t identifier){
  
  TLatex *tlDEta;
  Double_t x=0.015,y=0.28;
  if(ncolumns==3){
    x=0.19;
  }
  if(nrows==3){// these are hard coded number from an optimization
    y=0.25;
  }
  if(nrows==2){// these are hard coded number from an optimization
    y=0.28;
  }
  
  tlDEta=new TLatex(x/gPad->GetWNDC()+gPad->GetLeftMargin(),y/gPad->GetHNDC()+gPad->GetBottomMargin(),"|#Delta#eta| < 1"); 

  tlDEta->SetNDC();
  tlDEta->SetTextFont(43);
  tlDEta->SetTextSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);//0.06/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor);
  
  tlDEta->SetTextAlign(12);
  tlDEta->SetName(Form("tlDEta_%d",identifier));
  return tlDEta;
}


TLatex *GetCollSystem(Int_t collSystem,Int_t identifier){
  Double_t x=0.007,y=0.390;
  if(ncolumns==3){
    x=0.007;
  }
  if(nrows==3){// these are hard coded number from an optimization; has precedence on ncolumns
    x=0.15;
    y=0.25;
  }
  if(nrows==2){// these are hard coded number from an optimization
    y=0.3850; // was 0.390 at round 1   
    if(collSystem==0){
      x=0.15;
    }
    else if(collSystem==1){
      x=0.15;
    }
  }

  TLatex *tlsystem;
  if(style==-1){
    if(collSystem==0)tlsystem=new TLatex(0.24,0.8,"pp, #sqrt{#it{s}} = 5.02 TeV");
    else if(collSystem==1)tlsystem=new TLatex(0.24,0.8,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
    else return 0x0;
    tlsystem->SetNDC();
    tlsystem->SetTextFont(42);
    tlsystem->SetTextSize(0.03);
    tlsystem->SetTextAlign(12);
  }
  else{
    if(collSystem==0)tlsystem=new TLatex(x/gPad->GetWNDC()+gPad->GetLeftMargin(),y/gPad->GetHNDC()+gPad->GetBottomMargin(),"pp, #sqrt{#it{s}} = 5.02 TeV"); 
    else if(collSystem==1)tlsystem=new TLatex(x/gPad->GetWNDC()+gPad->GetLeftMargin(),y/gPad->GetHNDC()+gPad->GetBottomMargin(),"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV"); 
    tlsystem->SetNDC();
    tlsystem->SetTextFont(43);
    tlsystem->SetTextSize(24*innerPadHeight/referencePadHeight*resizeTextFactor);//0.06/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor);
    tlsystem->SetTextAlign(21);
  }
  

  tlsystem->SetName(Form("tlSystem_%d",identifier));
  return tlsystem;
}

TLatex *GetALICEtext(Int_t identifier){
  TLatex *alice;
Double_t x=0.21,y=0.390;
  if(ncolumns==3){
    x=0.21;
  }
  if(nrows==3){// these are hard coded number from an optimization
    y=0.25;
    x=0.20;// draft 2 was not present -> above value 0.21
  }
  if(nrows==2){// these are hard coded number from an optimization
    y=0.185;// draft 2 was 0.39
    x=0.20;// draft 2 was not present -> above value 0.21
  }

  if(style==-1){
    alice=new TLatex(0.75,0.85,"ALICE Preliminary");
    alice->SetNDC();
    alice->SetTextFont(42);
    alice->SetTextSize(0.03);
    alice->SetTextAlign(11);
  }
  else{
    alice= new TLatex(0.0145/gPad->GetWNDC()+gPad->GetLeftMargin(),0.395/gPad->GetHNDC()+gPad->GetBottomMargin(),"ALICE Preliminary"); 
    alice->SetNDC();
    alice->SetTextFont(43);
    alice->SetTextSize(19*innerPadHeight/referencePadHeight*resizeTextFactor);//0.06/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor);// draft 2 was: 28 *...
    alice->SetTextAlign(11);
  }


  //  TPaveText *alice = new TPaveText(0.012/gPad->GetWNDC()+gPad->GetLeftMargin(),0.26/gPad->GetHNDC()+gPad->GetBottomMargin(),0.3/gPad->GetWNDC()+gPad->GetLeftMargin(),0.28/gPad->GetHNDC()+gPad->GetBottomMargin(),"NDC");
  alice->SetName(Form("tlALICE_%d",identifier));
  return alice;
}

TLatex* GetDPtText(Int_t binD,Int_t identifier,Int_t addDEta=1){
  TLatex *tlasspt;
  Double_t x=0.025,y=0.35;
  if(ncolumns==3){
    //    x=0.035;
    x=0.15;
  }
  if(nrows==3){// these are hard coded number from an optimization
    y=0.215;
    //    if(addDEta==0)x=0.055;
    x=0.15;
  }
  if(nrows==2){// these are hard coded number from an optimization
    y=0.34;// was 0.35 in draft 2
    x=0.115;
  }

  if(style==-1){
    tlasspt=new TLatex(0.25,0.78,Form("%s, |#Delta#eta| < 1",strPtDCanvas[binD].Data()));
    tlasspt->SetNDC();
    tlasspt->SetTextFont(42);
    tlasspt->SetTextSize(0.03);
  }
  else{
    TString strTot=strPtDCanvas[binD];
    if(addDEta==1)strTot.Append(", |#Delta#eta| < 1");
    tlasspt= new TLatex(x/gPad->GetWNDC()+gPad->GetLeftMargin(),0.355/gPad->GetHNDC()+gPad->GetBottomMargin(),strTot.Data()); 
    tlasspt->SetNDC();
    tlasspt->SetTextFont(43);
    tlasspt->SetTextSize(24*innerPadHeight/referencePadHeight*resizeTextFactor);//  if font 42 is used try this: 0.06/(gPad->GetHNDC())*scaleHeightPads*resizeTextFactor) but see notes on top
    Printf("Height pad: %f, scal Height = %f",gPad->GetHNDC(),scaleHeightPads);
  }
  //  TPaveText *tlasspt = new TPaveText(0.012/gPad->GetWNDC()+gPad->GetLeftMargin(),0.26/gPad->GetHNDC()+gPad->GetBottomMargin(),0.3/gPad->GetWNDC()+gPad->GetLeftMargin(),0.28/gPad->GetHNDC()+gPad->GetBottomMargin(),"NDC");
  tlasspt->SetTextAlign(22);
  tlasspt->SetName(Form("tlAssocPt_%d",identifier));
  return tlasspt;

}

TLatex* GetTextSide(Int_t variable,Int_t identifier){
  
  TLatex *tlSide=new TLatex();    
  Double_t x=0.015,y=0.390;
  if(ncolumns==3){
    x=0.015;
  }
  if(nrows==3){// these are hard coded number from an optimization
    y=0.250;
  }
  if(nrows==2){// these are hard coded number from an optimization
    y=0.385;// draft 2 was 0.39, round 1 was 0.378
  }
  
  if(variable==2)return 0x0;
  else if(variable==0 || variable ==1)
    {
      if(style==-1){
	tlSide=new TLatex(0.25,0.85,"Near side");
	tlSide->SetNDC();
	tlSide->SetTextFont(42);
	tlSide->SetTextAlign(12);
	tlSide->SetTextSize(0.03*resizeTextFactor);
      }
      else{
	tlSide=new TLatex(0.15/gPad->GetWNDC()+gPad->GetLeftMargin(),0.395/gPad->GetHNDC()+gPad->GetBottomMargin(),"Near side");
	tlSide->SetNDC();
	tlSide->SetTextAlign(11);
	tlSide->SetTextFont(43);
	tlSide->SetTextSize(19*innerPadHeight/referencePadHeight*resizeTextFactor);// draft 2 was 28*...
      }
    }
  else {
    if(style==-1){
      tlSide=new TLatex(0.25,0.85,"Away side");
      tlSide->SetNDC();
      tlSide->SetTextFont(42);
	tlSide->SetTextAlign(12);
      tlSide->SetTextSize(0.03*resizeTextFactor);
    }
    else{
      tlSide=new TLatex(0.146/gPad->GetWNDC()+gPad->GetLeftMargin(),0.395/gPad->GetHNDC()+gPad->GetBottomMargin(),"Away side");
      tlSide->SetNDC();
      tlSide->SetTextFont(43);
      tlSide->SetTextAlign(11);
      tlSide->SetTextSize(19*innerPadHeight/referencePadHeight*resizeTextFactor);
    }
    

  }
  
  tlSide->SetName(Form("tlSideName_%d",identifier));
  return tlSide;
}

TH1D *GetAndPreparePP(Int_t binD,Int_t quantity,Int_t numsyst,TGraphAsymmErrors *&gr){


  TFile *f1=TFile::Open(Form("%s/CanvasFinalTrend%s_pthad%s.root",strFitResultPPb[numsyst].Data(),strquantityFile[quantity].Data(),strPtAss[0].Data()),"READ");
  TFile *f2=TFile::Open(Form("%s/CanvasFinalTrend%s_pthad%s.root",strFitResultPPb[numsyst].Data(),strquantityFile[quantity].Data(),strPtAss[1].Data()),"READ");
  TFile *f3=TFile::Open(Form("%s/CanvasFinalTrend%s_pthad%s.root",strFitResultPPb[numsyst].Data(),strquantityFile[quantity].Data(),strPtAss[2].Data()),"READ");  
  TCanvas *c1=(TCanvas*)f1->Get(Form("CanvasFinalTrend%s",strquantityFile[quantity].Data()));
  TCanvas *c2=(TCanvas*)f2->Get(Form("CanvasFinalTrend%s",strquantityFile[quantity].Data()));
  TCanvas *c3=(TCanvas*)f3->Get(Form("CanvasFinalTrend%s",strquantityFile[quantity].Data()));
  TGraphAsymmErrors *gr1=(TGraphAsymmErrors*)c1->FindObject(Form("fFullSystematics%s",strquantityFile[quantity].Data()));
  TGraphAsymmErrors *gr2=(TGraphAsymmErrors*)c2->FindObject(Form("fFullSystematics%s",strquantityFile[quantity].Data()));
  TGraphAsymmErrors *gr3=(TGraphAsymmErrors*)c3->FindObject(Form("fFullSystematics%s",strquantityFile[quantity].Data()));
  
  Double_t x[3] = {0.65,1.5,2.5};
  Double_t ex[3] = {0.35,0.5,0.5};
  Double_t y[3];
  Double_t eyl[3];  
  Double_t eyh[3];  
  Double_t *y1 = gr1->GetY();
  Double_t *y2 = gr2->GetY();
  Double_t *y3 = gr3->GetY();
  if(numsyst==0) {
    y[0] = y1[binD+1]; y[1] = y2[binD+1]; y[2] = y3[binD+1];
    eyl[0] = gr1->GetErrorYlow(binD+1); eyl[1] = gr2->GetErrorYlow(binD+1); eyl[2] = gr3->GetErrorYlow(binD+1);
    eyh[0] = gr1->GetErrorYhigh(binD+1); eyh[1] = gr2->GetErrorYhigh(binD+1); eyh[2] = gr3->GetErrorYhigh(binD+1);
  }    
  if(numsyst==1) {
    y[0] = y1[binD]; y[1] = y2[binD]; y[2] = y3[binD];
    eyl[0] = gr1->GetErrorYlow(binD); eyl[1] = gr2->GetErrorYlow(binD); eyl[2] = gr3->GetErrorYlow(binD);
    eyh[0] = gr1->GetErrorYhigh(binD); eyh[1] = gr2->GetErrorYhigh(binD); eyh[2] = gr3->GetErrorYhigh(binD);
  }

  gr = new TGraphAsymmErrors(3,x,y,ex,ex,eyl,eyh); //the real one!
  gr->SetFillColor(kGreen+2);
  gr->SetFillStyle(0);
  gr->SetName(Form("%sPPb",gr->GetName()));
  gr->SetMarkerColor(colSystem[numsyst]);
  gr->SetLineColor(colSystem[numsyst]);
  gr->SetLineWidth(1);
  gr->SetMarkerStyle(markerStyle[numsyst]);
  gr->SetMarkerSize(markersize);
  for(Int_t iPoint=0;iPoint<3;iPoint++) {
    gr->SetPointError(iPoint,0.3*gr->GetErrorXlow(iPoint),0.3*gr->GetErrorXhigh(iPoint),gr->GetErrorYlow(iPoint),gr->GetErrorYhigh(iPoint));
  }

  TH1D *hPPbInput1=(TH1D*)c1->FindObject(Form("FinalTrend%s",strquantityFile[quantity].Data()));
  TH1D *hPPbInput2=(TH1D*)c2->FindObject(Form("FinalTrend%s",strquantityFile[quantity].Data()));
  TH1D *hPPbInput3=(TH1D*)c3->FindObject(Form("FinalTrend%s",strquantityFile[quantity].Data()));
  Double_t xaxis[6] = {0.,0.3,1.,2.,3.,4.};
  TH1D *hPPb=new TH1D(Form("%sPPb",hPPbInput1->GetName()),"Hist_obs_vsPtAss",5,xaxis); //the real one!
  if(numsyst==1) {
    hPPb->SetBinContent(1,-1);
    hPPb->SetBinContent(2,hPPbInput1->GetBinContent(binD+2));
    hPPb->SetBinContent(3,hPPbInput2->GetBinContent(binD+2));
    hPPb->SetBinContent(4,hPPbInput3->GetBinContent(binD+2));
    hPPb->SetBinContent(5,-1);
    hPPb->SetBinError(1,0);
    hPPb->SetBinError(2,hPPbInput1->GetBinError(binD+2));
    hPPb->SetBinError(3,hPPbInput2->GetBinError(binD+2));
    hPPb->SetBinError(4,hPPbInput3->GetBinError(binD+2));
    hPPb->SetBinError(5,0);
    //printf("BinD = %d\n, pPb in 0.3-1 is %f\n",binD,hPPbInput1->GetBinContent(binD+2));
  } else {  //pp has also additional bin 2-3 to be skipped, so bin 3-5 is #3 (0-2, 2-3, 3-5)
    hPPb->SetBinContent(1,-1);
    hPPb->SetBinContent(2,hPPbInput1->GetBinContent(binD+3));
    hPPb->SetBinContent(3,hPPbInput2->GetBinContent(binD+3));
    hPPb->SetBinContent(4,hPPbInput3->GetBinContent(binD+3));
    hPPb->SetBinContent(5,-1);
    hPPb->SetBinError(1,0);
    hPPb->SetBinError(2,hPPbInput1->GetBinError(binD+3));
    hPPb->SetBinError(3,hPPbInput2->GetBinError(binD+3));
    hPPb->SetBinError(4,hPPbInput3->GetBinError(binD+3));
    hPPb->SetBinError(5,0);
    //printf("BinD = %d\n, pp in 0.3-1 is %f\n",binD,hPPbInput1->GetBinContent(binD+3));
  }
  hPPb->SetLineColor(colSystem[numsyst]);
  hPPb->SetLineWidth(1);
  hPPb->SetMarkerColor(colSystem[numsyst]);
  hPPb->SetMarkerStyle(markerStyle[numsyst]);
  hPPb->SetMarkerSize(markersize);

  hPPb->SetXTitle("Associated track #it{p}_{T} (GeV/#it{c})");
  hPPb->SetYTitle(yaxisTitle[quantity].Data());
  if(style==-1){
    hPPb->GetYaxis()->SetTitleSize(0.04);
    hPPb->GetYaxis()->SetTitleOffset(1.2);
    hPPb->GetYaxis()->SetLabelSize(0.04);
    hPPb->GetXaxis()->SetTitleSize(0.04);
    hPPb->GetXaxis()->SetLabelSize(0.04);
  }
  else {
    hPPb->GetYaxis()->SetTitle("");      

    hPPb->GetXaxis()->SetRangeUser(0,3.4);
    hPPb->GetYaxis()->SetTitleFont(43);
    hPPb->GetYaxis()->SetLabelFont(43);
    hPPb->GetXaxis()->SetTitleFont(43);
    hPPb->GetXaxis()->CenterTitle();
    hPPb->GetYaxis()->CenterTitle();
    hPPb->GetXaxis()->SetLabelFont(43);
    hPPb->GetYaxis()->SetTitleSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
    hPPb->GetYaxis()->SetTitleOffset(ytitleoffset*(gPad->GetHNDC())/scaleHeightPads/resizeTextFactor);
    hPPb->GetXaxis()->SetTitleOffset(xtitleoffset*(gPad->GetHNDC())/scaleHeightPads/resizeTextFactor);
    hPPb->GetYaxis()->SetLabelSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
    hPPb->GetXaxis()->SetTitleSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
    hPPb->GetXaxis()->SetLabelSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
  }
  if(style<=0){
    hPPb->GetYaxis()->SetRangeUser(0,maxRangePPb[binD][quantity]);
  }

//  AdaptRangeHist(hPPb,minptData,maxptData); //NOT TO BE CALLED VS PTASS!!
//  AdaptRangeTGraph(gr,minptData,maxptData); //NOT TO BE CALLED VS PTASS!!

  return hPPb;
}

TH1D *GetAndPreparePPb(Int_t binD,Int_t quantity,Int_t numsyst,TGraphAsymmErrors *&gr, TGraphAsymmErrors *&grV2tot){


  TFile *f1=TFile::Open(Form("%s/CanvasFinalTrend%s_pthad%s.root",strFitResultPPb[numsyst].Data(),strquantityFile[quantity].Data(),strPtAss[0].Data()),"READ");
  TFile *f2=TFile::Open(Form("%s/CanvasFinalTrend%s_pthad%s.root",strFitResultPPb[numsyst].Data(),strquantityFile[quantity].Data(),strPtAss[1].Data()),"READ");
  TFile *f3=TFile::Open(Form("%s/CanvasFinalTrend%s_pthad%s.root",strFitResultPPb[numsyst].Data(),strquantityFile[quantity].Data(),strPtAss[2].Data()),"READ");  
  TCanvas *c1=(TCanvas*)f1->Get(Form("CanvasFinalTrend%s",strquantityFile[quantity].Data()));
  TCanvas *c2=(TCanvas*)f2->Get(Form("CanvasFinalTrend%s",strquantityFile[quantity].Data()));
  TCanvas *c3=(TCanvas*)f3->Get(Form("CanvasFinalTrend%s",strquantityFile[quantity].Data()));
  TGraphAsymmErrors *gr1=(TGraphAsymmErrors*)c1->FindObject(Form("fFullSystematics%s",strquantityFile[quantity].Data()));
  TGraphAsymmErrors *gr2=(TGraphAsymmErrors*)c2->FindObject(Form("fFullSystematics%s",strquantityFile[quantity].Data()));
  TGraphAsymmErrors *gr3=(TGraphAsymmErrors*)c3->FindObject(Form("fFullSystematics%s",strquantityFile[quantity].Data()));
  TGraphAsymmErrors *grV1=(TGraphAsymmErrors*)c1->FindObject(Form("fv2Systematics%s",strquantityFile[quantity].Data()));
  TGraphAsymmErrors *grV2=(TGraphAsymmErrors*)c2->FindObject(Form("fv2Systematics%s",strquantityFile[quantity].Data()));
  TGraphAsymmErrors *grV3=(TGraphAsymmErrors*)c3->FindObject(Form("fv2Systematics%s",strquantityFile[quantity].Data()));

  Double_t x[3] = {0.65,1.5,2.5};
  Double_t ex[3] = {0.35,0.5,0.5};
  Double_t y[3];
  Double_t eyl[3];  
  Double_t eyh[3];  
  Double_t *y1 = gr1->GetY();
  Double_t *y2 = gr2->GetY();
  Double_t *y3 = gr3->GetY();
  Double_t yV[3];
  Double_t eyVl[3];
  Double_t eyVh[3];
  Double_t *yV1 = grV1->GetY();
  Double_t *yV2 = grV2->GetY();
  Double_t *yV3 = grV3->GetY();
  if(numsyst==0) {
    y[0] = y1[binD+1]; y[1] = y2[binD+1]; y[2] = y3[binD+1];
    eyl[0] = gr1->GetErrorYlow(binD+1); eyl[1] = gr2->GetErrorYlow(binD+1); eyl[2] = gr3->GetErrorYlow(binD+1);
    eyh[0] = gr1->GetErrorYhigh(binD+1); eyh[1] = gr2->GetErrorYhigh(binD+1); eyh[2] = gr3->GetErrorYhigh(binD+1);
  }    
  if(numsyst==1) {
    y[0] = y1[binD]; y[1] = y2[binD]; y[2] = y3[binD];
    eyl[0] = gr1->GetErrorYlow(binD); eyl[1] = gr2->GetErrorYlow(binD); eyl[2] = gr3->GetErrorYlow(binD);
    eyh[0] = gr1->GetErrorYhigh(binD); eyh[1] = gr2->GetErrorYhigh(binD); eyh[2] = gr3->GetErrorYhigh(binD);
    if(plotv2unc) {
      yV[0] = yV1[binD]; yV[1] = yV2[binD]; yV[2] = yV3[binD];
      eyVl[0] = grV1->GetErrorYlow(binD); eyVl[1] = grV2->GetErrorYlow(binD); eyVl[2] = grV3->GetErrorYlow(binD);
      eyVh[0] = grV1->GetErrorYhigh(binD); eyVh[1] = grV2->GetErrorYhigh(binD); eyVh[2] = grV3->GetErrorYhigh(binD);    
    }
  }

  gr = new TGraphAsymmErrors(3,x,y,ex,ex,eyl,eyh); //the real one!
  gr->SetFillColor(kGreen+2);
  gr->SetFillStyle(0);
  gr->SetName(Form("%sPPb",gr->GetName()));
  gr->SetMarkerColor(colSystem[numsyst]);
  gr->SetLineColor(colSystem[numsyst]);
  gr->SetLineWidth(1);
  gr->SetMarkerStyle(markerStyle[numsyst]);
  gr->SetMarkerSize(markersize);
  for(Int_t iPoint=0;iPoint<3;iPoint++) {
    gr->SetPointError(iPoint,0.3*gr->GetErrorXlow(iPoint),0.3*gr->GetErrorXhigh(iPoint),gr->GetErrorYlow(iPoint),gr->GetErrorYhigh(iPoint));
  }
  if(plotv2unc) {
    grV2tot = new TGraphAsymmErrors(3,x,yV,ex,ex,eyVl,eyVh); //the real one!
    grV2tot->SetFillColor(grV1->GetFillColor());
    grV2tot->SetFillStyle(grV1->GetFillStyle());
    grV2tot->SetName(Form("%sPPb",gr->GetName()));
    grV2tot->SetMarkerColor(colSystem[numsyst]);
    grV2tot->SetLineColor(grV1->GetLineColor());
    grV2tot->SetLineWidth(grV1->GetLineWidth());
    grV2tot->SetMarkerStyle(markerStyle[numsyst]);
    grV2tot->SetMarkerSize(markersize);
    for(Int_t iPoint=0;iPoint<3;iPoint++) {
      grV2tot->SetPointError(iPoint,0.3*grV2tot->GetErrorXlow(iPoint),0.3*grV2tot->GetErrorXhigh(iPoint),grV2tot->GetErrorYlow(iPoint),grV2tot->GetErrorYhigh(iPoint));
      printf("vals y and errs %f, %f, %f - %f %f %f\n",yV[0],yV[1],yV[2],eyVh[0],eyVh[1],eyVh[2]);
    }
  }


  TH1D *hPPbInput1=(TH1D*)c1->FindObject(Form("FinalTrend%s",strquantityFile[quantity].Data()));
  TH1D *hPPbInput2=(TH1D*)c2->FindObject(Form("FinalTrend%s",strquantityFile[quantity].Data()));
  TH1D *hPPbInput3=(TH1D*)c3->FindObject(Form("FinalTrend%s",strquantityFile[quantity].Data()));
  Double_t xaxis[6] = {0.,0.3,1.,2.,3.,4.};
  TH1D *hPPb=new TH1D(Form("%sPPb",hPPbInput1->GetName()),"Hist_obs_vsPtAss",5,xaxis); //the real one!
  if(numsyst==1) {
    hPPb->SetBinContent(1,-1);
    hPPb->SetBinContent(2,hPPbInput1->GetBinContent(binD+2));
    hPPb->SetBinContent(3,hPPbInput2->GetBinContent(binD+2));
    hPPb->SetBinContent(4,hPPbInput3->GetBinContent(binD+2));
    hPPb->SetBinContent(5,-1);
    hPPb->SetBinError(1,0);
    hPPb->SetBinError(2,hPPbInput1->GetBinError(binD+2));
    hPPb->SetBinError(3,hPPbInput2->GetBinError(binD+2));
    hPPb->SetBinError(4,hPPbInput3->GetBinError(binD+2));
    hPPb->SetBinError(5,0);
    //printf("BinD = %d\n, pPb in 0.3-1 is %f\n",binD,hPPbInput1->GetBinContent(binD+2));
  } else {  //pp has also additional bin 2-3 to be skipped, so bin 3-5 is #3 (0-2, 2-3, 3-5)
    hPPb->SetBinContent(1,-1);
    hPPb->SetBinContent(2,hPPbInput1->GetBinContent(binD+3));
    hPPb->SetBinContent(3,hPPbInput2->GetBinContent(binD+3));
    hPPb->SetBinContent(4,hPPbInput3->GetBinContent(binD+3));
    hPPb->SetBinContent(5,-1);
    hPPb->SetBinError(1,0);
    hPPb->SetBinError(2,hPPbInput1->GetBinError(binD+3));
    hPPb->SetBinError(3,hPPbInput2->GetBinError(binD+3));
    hPPb->SetBinError(4,hPPbInput3->GetBinError(binD+3));
    hPPb->SetBinError(5,0);
    //printf("BinD = %d\n, pp in 0.3-1 is %f\n",binD,hPPbInput1->GetBinContent(binD+3));
  }
  hPPb->SetLineColor(colSystem[numsyst]);
  hPPb->SetLineWidth(1);
  hPPb->SetMarkerColor(colSystem[numsyst]);
  hPPb->SetMarkerStyle(markerStyle[numsyst]);
  hPPb->SetMarkerSize(markersize);

  hPPb->SetXTitle("Associated track #it{p}_{T} (GeV/#it{c})");
  hPPb->SetYTitle(yaxisTitle[quantity].Data());
  if(style==-1){
    hPPb->GetYaxis()->SetTitleSize(0.04);
    hPPb->GetYaxis()->SetTitleOffset(1.2);
    hPPb->GetYaxis()->SetLabelSize(0.04);
    hPPb->GetXaxis()->SetTitleSize(0.04);
    hPPb->GetXaxis()->SetLabelSize(0.04);
  }
  else {
    hPPb->GetYaxis()->SetTitle("");      

    hPPb->GetXaxis()->SetRangeUser(0,3.4);
    hPPb->GetYaxis()->SetTitleFont(43);
    hPPb->GetYaxis()->SetLabelFont(43);
    hPPb->GetXaxis()->SetTitleFont(43);
    hPPb->GetXaxis()->CenterTitle();
    hPPb->GetYaxis()->CenterTitle();
    hPPb->GetXaxis()->SetLabelFont(43);
    hPPb->GetYaxis()->SetTitleSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
    hPPb->GetYaxis()->SetTitleOffset(ytitleoffset*(gPad->GetHNDC())/scaleHeightPads/resizeTextFactor);
    hPPb->GetXaxis()->SetTitleOffset(xtitleoffset*(gPad->GetHNDC())/scaleHeightPads/resizeTextFactor);
    hPPb->GetYaxis()->SetLabelSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
    hPPb->GetXaxis()->SetTitleSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
    hPPb->GetXaxis()->SetLabelSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
  }
  if(style<=0){
    hPPb->GetYaxis()->SetRangeUser(0,maxRangePPb[binD][quantity]);
  }

  return hPPb;
}

/*
TH1D *GetAndPreparePPb(Int_t binass,Int_t quantity,Int_t numsyst,TGraphAsymmErrors *&gr){


  TFile *f=TFile::Open(Form("%s/CanvasFinalTrend%s_pthad%s.root",strFitResultPPb[numsyst].Data(),strquantityFile[quantity].Data(),strPtAss[binass].Data()),"READ");
  TCanvas *c=(TCanvas*)f->Get(Form("CanvasFinalTrend%s",strquantityFile[quantity].Data()));
  gr=(TGraphAsymmErrors*)c->FindObject(Form("fFullSystematics%s",strquantityFile[quantity].Data()));
  gr->SetName(Form("%sPPb",gr->GetName()));
  for(Int_t iPoint=0;iPoint<5;iPoint++) {
    gr->SetPointError(iPoint,0.7*gr->GetErrorXlow(iPoint),0.7*gr->GetErrorXhigh(iPoint),gr->GetErrorYlow(iPoint),gr->GetErrorYhigh(iPoint));
    double xx,yy; gr->GetPoint(iPoint,xx,yy);
    printf("syst %d; point %d; x = %f, y = %f\n",numsyst,iPoint,xx,yy);
  }
  if(numsyst==1) gr->RemovePoint(4);

  TH1D *hPPb=(TH1D*)c->FindObject(Form("FinalTrend%s",strquantityFile[quantity].Data()));
  hPPb->SetName(Form("%sPPb",hPPb->GetName()));
  hPPb->SetLineColor(colSystem[numsyst]);
  hPPb->SetLineWidth(1);
  hPPb->SetMarkerColor(colSystem[numsyst]);
  hPPb->SetMarkerStyle(markerStyle[numsyst]);
  hPPb->SetMarkerSize(markersize);
  gr->SetMarkerColor(colSystem[numsyst]);
  gr->SetLineColor(colSystem[numsyst]);
  gr->SetLineWidth(1);
  gr->SetMarkerStyle(markerStyle[numsyst]);
  gr->SetMarkerSize(markersize);

  hPPb->SetXTitle("D meson #it{p}_{T} (GeV/#it{c})");
  hPPb->SetYTitle(yaxisTitle[quantity].Data());
  if(style==-1){
    hPPb->GetYaxis()->SetTitleSize(0.04);
    hPPb->GetYaxis()->SetTitleOffset(1.2);
    hPPb->GetYaxis()->SetLabelSize(0.04);
    hPPb->GetXaxis()->SetTitleSize(0.04);
    hPPb->GetXaxis()->SetLabelSize(0.04);
  }
  else {
    hPPb->GetYaxis()->SetTitle("");      

    hPPb->GetXaxis()->SetRangeUser(0,24.2);
    hPPb->GetYaxis()->SetTitleFont(43);
    hPPb->GetYaxis()->SetLabelFont(43);
    hPPb->GetXaxis()->SetTitleFont(43);
    hPPb->GetXaxis()->CenterTitle();
    hPPb->GetYaxis()->CenterTitle();
    hPPb->GetXaxis()->SetLabelFont(43);
    hPPb->GetYaxis()->SetTitleSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
    hPPb->GetYaxis()->SetTitleOffset(ytitleoffset*(gPad->GetHNDC())/scaleHeightPads/resizeTextFactor);
    hPPb->GetXaxis()->SetTitleOffset(xtitleoffset*(gPad->GetHNDC())/scaleHeightPads/resizeTextFactor);
    hPPb->GetYaxis()->SetLabelSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
    hPPb->GetXaxis()->SetTitleSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
    hPPb->GetXaxis()->SetLabelSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
  }
  if(style<=0){
    hPPb->GetYaxis()->SetRangeUser(0,maxRangePPb[binass][quantity]);
  }

  AdaptRangeHist(hPPb,minptData,maxptData);
  AdaptRangeTGraph(gr,minptData,maxptData);

  return hPPb;
}
*/
TCanvas* Compare(Int_t binD,Int_t quantity,TPad *pd,Int_t textlegendOptions){

  printf("Preparing pp... (%d)\n",quantity);
  TGraphAsymmErrors *grPPb_1;
  TH1D *hPPb_1=GetAndPreparePP(binD,quantity,0,grPPb_1);

  printf("Preparing p-Pb...(%d)\n",quantity);
  TGraphAsymmErrors *grPPb_2;
  TGraphAsymmErrors *grV2tot;
  TH1D *hPPb_2=GetAndPreparePPb(binD,quantity,1,grPPb_2,grV2tot);

  TCanvas *cout=0x0;
  if(!pd){
    cout=CreateCanvasWithDefaultStyle(Form("%sComparisonPPtoPPbBinD%d",strquantityFile[quantity].Data(),binD));
    cout->cd();
  }
  else {
    pd->cd();
  }
  //new TCanvas(Form("NSyieldComparisonBinAss%d",binass),Form("NSyieldComparisonBinAss%d",binass),800,800);

  TH2D *hDraw;
  //  hPP->SetXTitle("D meson #it{p}_{T} (GeV/#it{c})");
  //  hPP->SetYTitle(yaxisTitle[quantity].Data());
  if(style==-1){
//     hPP->GetYaxis()->SetTitleSize(0.04);
//     hPP->GetYaxis()->SetTitleOffset(1.2);
//     hPP->GetYaxis()->SetLabelSize(0.04);
//     hPP->GetXaxis()->SetTitleSize(0.04);
//     hPP->GetXaxis()->SetLabelSize(0.04);
    hPPb_1->Draw();
  }
  else {
    hDraw=new TH2D(Form("hDraw%d",10*quantity+binD),"",100,0,4.5,200,0,10);
    hDraw->GetYaxis()->SetTitle("");      
    //    hPP->GetYaxis()->SetTitle("");      

    hDraw->GetXaxis()->SetRangeUser(0,3.4);
    hDraw->GetYaxis()->SetTitleFont(43);
    hDraw->GetYaxis()->SetLabelFont(43);
    hDraw->GetXaxis()->SetTitleFont(43);
    hDraw->GetXaxis()->CenterTitle();
    hDraw->GetYaxis()->CenterTitle();
    hDraw->GetXaxis()->SetLabelFont(43);
    hDraw->GetYaxis()->SetTitleSize(24*innerPadHeight/referencePadHeight*resizeTextFactor);
    hDraw->GetYaxis()->SetTitleOffset(ytitleoffset*(gPad->GetHNDC())/scaleHeightPads/resizeTextFactor);
    hDraw->GetXaxis()->SetTitleOffset(xtitleoffset*(gPad->GetHNDC())/scaleHeightPads/resizeTextFactor);
    hDraw->GetYaxis()->SetLabelSize(22*innerPadHeight/referencePadHeight*resizeTextFactor);
    hDraw->GetXaxis()->SetTitleSize(22*innerPadHeight/referencePadHeight*resizeTextFactor);
    hDraw->GetXaxis()->SetLabelSize(22*innerPadHeight/referencePadHeight*resizeTextFactor);

    if(textlegendOptions%10==2||textlegendOptions%10==3){
      hDraw->SetYTitle(yaxisTitle[quantity].Data());
    }
    else {
      hDraw->GetYaxis()->SetLabelSize(0);
    }
    if(textlegendOptions%10==1||textlegendOptions%10==3){
      hDraw->SetXTitle("Assoc. track #it{p}_{T} (GeV/#it{c})");
    }
    else {
      hDraw->GetXaxis()->SetLabelSize(0);
    }

    // Set same drawing settings also to hPP: not needed, just to have it also in that histo
//     hPP->GetXaxis()->SetRangeUser(0,16.9);
//     hPP->GetYaxis()->SetTitleFont(43);
//     hPP->GetYaxis()->SetLabelFont(43);
//     hPP->GetXaxis()->SetTitleFont(43);
//     hPP->GetXaxis()->CenterTitle();
//     hPP->GetYaxis()->CenterTitle();
//     hPP->GetXaxis()->SetLabelFont(43);
//     hPP->GetYaxis()->SetTitleSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
//     hPP->GetYaxis()->SetTitleOffset(ytitleoffset*(gPad->GetHNDC())/scaleHeightPads/resizeTextFactor);
//     hPP->GetXaxis()->SetTitleOffset(xtitleoffset*(gPad->GetHNDC())/scaleHeightPads/resizeTextFactor);
//     hPP->GetYaxis()->SetLabelSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
//     hPP->GetXaxis()->SetTitleSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);
//     hPP->GetXaxis()->SetLabelSize(28*innerPadHeight/referencePadHeight*resizeTextFactor);

    hDraw->Draw();
  //  hPPb_1->Draw("same");
  }

  //  hPP->Draw("E0X0");// to avoid plotting the error along x
  if(style>0){
    hDraw->GetYaxis()->SetRangeUser(0,maxRangePPb[2][quantity]);
  }
  else {
    hPPb_1->GetYaxis()->SetRangeUser(0,maxRangePPb[binD][quantity]);
  }
 
//  grPPb_1->Draw("E2");
//  grPPbV2_1->Draw("E2");

//  hPPb->Draw("same");

// p-Pb: displace the syst errors along x axis (Fabio)
  Double_t shift=0.05;
  Double_t binsx=grPPb_1->GetN();
  for(int i=0;i<binsx;i++) {
    Double_t x,y,ex1,ex2,ey1,ey2;
    
    grPPb_1->GetPoint(i,x,y);
    //printf("pp ) Bin %d, point x = %f, y = %f\n",i,x,y);
    ex1=grPPb_1->GetErrorXlow(i); ex2=grPPb_1->GetErrorXhigh(i); ey1=grPPb_1->GetErrorYlow(i); ey2=grPPb_1->GetErrorYhigh(i);
    grPPb_1->SetPoint(i,x+shift*(-1),y);
    grPPb_1->SetPointError(i,ex1,ex2,ey1,ey2);

    grPPb_2->GetPoint(i,x,y);
    //printf("pPb) Bin %d, point x = %f, y = %f\n",i,x,y);
    ex1=grPPb_2->GetErrorXlow(i); ex2=grPPb_2->GetErrorXhigh(i); ey1=grPPb_2->GetErrorYlow(i); ey2=grPPb_2->GetErrorYhigh(i);
    grPPb_2->SetPoint(i,x+shift*(0),y);
    grPPb_2->SetPointError(i,ex1,ex2,ey1,ey2);
  }

/****************************************/
//REMOVE AD HOC POINTS WITH LARGE ERRORS!!
/****************************************/

  if(quantity==3 || quantity==4) { //remove NSy and NSw of 0.3-1 in 3-5 and 16-24
     if(binD==0) grPPb_1->RemovePoint(0);
     if(binD==3) grPPb_1->RemovePoint(0);
     printf("Removing NSy and NSw of 0.3-1 in 3-5 and 16-24\n");
   }

  if(plotv2unc) grV2tot->Draw("E2");
  grPPb_1->Draw("E2");
  grPPb_2->Draw("E2");

//Conversion of pPb TH1F to TGraph to displace the points along x axis (Fabio)
  TGraphAsymmErrors *gr_points_PPb1, *gr_pointCount_PPb1, *gr_points_PPb2, *gr_pointCount_PPb2, *gr_points_PPb3, *gr_pointCount_PPb3, *gr_points_PPb4, *gr_pointCount_PPb4;
    ConvertTH1ToTGraphAsymmError(hPPb_1,gr_points_PPb1,shift*(-1));
    ConvertTH1ToTGraphAsymmError(hPPb_1,gr_pointCount_PPb1,shift*(-1));    
    ConvertTH1ToTGraphAsymmError(hPPb_2,gr_points_PPb2,shift*(0));
    ConvertTH1ToTGraphAsymmError(hPPb_2,gr_pointCount_PPb2,shift*(0));  
  
  if(quantity==3 || quantity==4) { //remove NSy and NSw of 0.3-1 in 3-5 and 16-24
     if(binD==0) {gr_points_PPb1->RemovePoint(0); gr_pointCount_PPb1->RemovePoint(0);}
     if(binD==3) {gr_points_PPb1->RemovePoint(0); gr_pointCount_PPb1->RemovePoint(0);}
     printf("Removing NSy and NSw of 0.3-1 in 3-5 and 16-24\n");
   }


  gr_points_PPb1->SetLineColor(colSystem[0]);
  gr_points_PPb1->SetLineWidth(1);
  gr_points_PPb1->SetMarkerColor(colSystem[0]);
  gr_points_PPb1->SetMarkerStyle(markerStyle[0]);
  gr_points_PPb1->SetMarkerSize(markersize);
  gr_points_PPb1->Draw("samePZ");
  
  gr_pointCount_PPb1->SetLineColor(colSystem[0]);
  gr_pointCount_PPb1->SetLineWidth(1);
  gr_pointCount_PPb1->SetMarkerStyle(colSystem[0]);
  gr_pointCount_PPb1->SetMarkerColor(kBlue);
  gr_pointCount_PPb1->SetMarkerSize(markersize);
  gr_pointCount_PPb1->Draw("samePZ");

  gr_points_PPb2->SetLineColor(colSystem[1]);
  gr_points_PPb2->SetLineWidth(1);
  gr_points_PPb2->SetMarkerColor(colSystem[1]);
  gr_points_PPb2->SetMarkerStyle(markerStyle[1]);
  gr_points_PPb2->SetMarkerSize(markersize);
  gr_points_PPb2->Draw("samePZ");
  
  gr_pointCount_PPb2->SetLineColor(colSystem[1]);
  gr_pointCount_PPb2->SetLineWidth(1);
  gr_pointCount_PPb2->SetMarkerStyle(colSystem[1]);
  gr_pointCount_PPb2->SetMarkerColor(kRed);
  gr_pointCount_PPb2->SetMarkerSize(markersize);
  gr_pointCount_PPb2->Draw("samePZ");

  if(style==-1){
        
    TLatex *tlAssYieldPt=GetDPtText(binD,10*quantity+binD,0);
    tlAssYieldPt->Draw();

    TLatex *tlALICE=GetALICEtext(10*quantity+binD);
    tlALICE->Draw();
    
    if(quantity!=2){
      TLatex *tlSide=GetTextSide(quantity,10*quantity+binD);
      tlSide->Draw();
    }
  }
  else{

    if(textlegendOptions%100>=10){
      TLatex *tlAssYieldPt=GetDPtText(binD,10*quantity+binD,0);
      tlAssYieldPt->Draw();
    }
    if(textlegendOptions%1000>=100){
      TLegend * legend;
      legend = new TLegend(0.005/gPad->GetWNDC()+gPad->GetLeftMargin(),0.175/gPad->GetHNDC()+gPad->GetBottomMargin(),0.18/gPad->GetWNDC()+gPad->GetLeftMargin(),0.332/gPad->GetHNDC()+gPad->GetBottomMargin());// draft 2 (2 lines only, rapidity on the same line also for p-Pb): 0.002/gPad->GetWNDC()+gPad->GetLeftMargin(),0.23/gPad->GetHNDC()+gPad->GetBottomMargin(),0.15/gPad->GetWNDC()+gPad->GetLeftMargin(),0.30/gPad->GetHNDC()+gPad->GetBottomMargin()
      legend->SetTextFont(43);
      legend->SetTextAlign(12);
      legend->SetEntrySeparation(0.3);
      legend->SetLineColor(kWhite);
      legend->SetTextSize(17*innerPadHeight/referencePadHeight*resizeTextFactor);
      legend->AddEntry(hPPb_1,"pp, #sqrt{s} = 5.02 TeV,","lep");
      legend->AddEntry((TObject*)0,"|#it{y}^{D}_{cms}| < 0.5","");
      //legend->AddEntry((TObject*)0,"","");
      legend->AddEntry(hPPb_2,"p-Pb, #sqrt{s_{NN}} = 5.02 TeV,","lep");
      legend->AddEntry((TObject*)0,"-0.96 < #it{y}^{D}_{cms} < 0.04","");
      if(plotv2unc) legend->AddEntry(grV2tot,"Syst uncertainty from v_{2}","f");
      legend->Draw();

     /*TLegend *legendSuperimp=GetLegendDataPointsFake(hPP,hPPbSuperimp,10*quantity+binD);
      legendSuperimp->Draw("same");*/
    }
    if(textlegendOptions%10000>=1000){
      TLatex *tlALICE=GetALICEtext(10*quantity+binD);
      tlALICE->Draw();
    }
    if(textlegendOptions%100000>=10000){
      if(quantity!=2){
  TLatex *tlSide=GetTextSide(quantity,10*quantity+binD);
  tlSide->Draw();
      }
    }
  }
  
  return cout;
  
}


void InitMCobjects(){
  hMC=new TH1D*[nmodels];
  grMC=new TGraphAsymmErrors*[nmodels];
}

void CompareFitResults_vsPtAss_UniqueCanvas(){
  gStyle->SetOptStat(0000);
  TCanvas *cFinalPaperStyle;
    cFinalPaperStyle=new TCanvas("cPPvsPPbFitResultsFinalPaperStyle","pp vs.pPb fit results ",1800./1200.*canvasheight,canvasheight);
    //SetPadStyle(cFinalPaperStyle);
    //    cFinalPaperStyle->SetBottomMargin(0);
    //    cFinalPaperStyle->SetTopMargin(0);
    cFinalPaperStyle->Divide(4,2,0.0,0.0,0);
    cFinalPaperStyle->SetTicky();
    cFinalPaperStyle->SetTickx();
    cFinalPaperStyle->SetFrameBorderMode(0);

    Set4x2PadPositions(cFinalPaperStyle);
    cFinalPaperStyle->Modified();
    cFinalPaperStyle->Update();
    //    Convert3x2Matrix(cFinalPaperStyle,kTRUE);

  Int_t orderD[4]={0,1,2,3};  
  for(Int_t jp=0;jp<=3;jp++){// First loop for NS yield
    Int_t needTitle=0;
    if(jp==0)needTitle=2;
    gStyle->SetOptStat(0000);
    TPad *pd=(TPad*)cFinalPaperStyle->cd(jp+1);
    SetPadStyle(pd);
    pd->cd();    
    if(jp==1){// identifier set as: 10*quantity+binass; NS --> 0 ; binass == jp (not orderAssoc[jp])
      Compare(orderD[jp],0,pd,100+needTitle);    
    }
    else if(jp==0){// identifier set as: 10*quantity+binass; NS:
      Compare(orderD[jp],0,pd,needTitle);    
      TLatex *tlALICE=GetALICEtext(orderD[jp]);
      tlALICE->Draw();
      TLatex *tlSide=GetTextSide(0,orderD[jp]);
      tlSide->Draw();
    }
    else{
      Compare(orderD[jp],0,pd,0+needTitle);    
    }
    TLatex *tlAssYieldPt=GetDPtText(orderD[jp],orderD[jp],0);
    tlAssYieldPt->Draw();
    if(jp==0) {
      TLatex *tlDeta=GetDRapForSystem(0,orderD[jp],1);
      tlDeta->Draw();
    }

    gStyle->SetOptStat(0000);   
  }

  for(Int_t jp=0;jp<=3;jp++){// second loop --> sigma
    Int_t needTitle=0;
    if(jp==0)needTitle=3;
    else needTitle=1;
    gStyle->SetOptStat(0000);
    TPad *pd=(TPad*)cFinalPaperStyle->cd(jp+5);
    SetPadStyle(pd);
    pd->cd();
    gStyle->SetOptStat(0000);
    Compare(orderD[jp],1,pd,needTitle); 

    gStyle->SetOptStat(0000);   
  }
  for(Int_t j=8;j>=1;j--){
    TPad *pd=(TPad*)cFinalPaperStyle->cd(j);
    pd->Draw();
  }

  TString outdir = Form("%s/ScriptOutput/ReflectedPlots/StdRebin/AllPlots/Averages/ComparisonTopPb",basicdir.Data());
  
  cFinalPaperStyle->Modified();
  cFinalPaperStyle->Update();
  cFinalPaperStyle->SaveAs(Form("%s/CompareFitResults_ppVspPb_5TeV_VsAssPt.root",outdir.Data()));
  cFinalPaperStyle->SaveAs(Form("%s/CompareFitResults_ppVspPb_5TeV_VsAssPt.eps",outdir.Data()));
  cFinalPaperStyle->SaveAs(Form("%s/CompareFitResults_ppVspPb_5TeV_VsAssPt.png",outdir.Data()));
  cFinalPaperStyle->SaveAs(Form("%s/CompareFitResults_ppVspPb_5TeV_VsAssPt.pdf",outdir.Data()));

  return;

}

void CompareFitResults_vsPtAss_UniqueCanvas_AwaySide(){
  gStyle->SetOptStat(0000);
  TCanvas *cFinalPaperStyle;
    cFinalPaperStyle=new TCanvas("cPPvsPPbFitResultsFinalPaperStyle","pp vs.pPb fit results ",1800./1200.*canvasheight,canvasheight);
    //SetPadStyle(cFinalPaperStyle);
    //    cFinalPaperStyle->SetBottomMargin(0);
    //    cFinalPaperStyle->SetTopMargin(0);
    cFinalPaperStyle->Divide(4,2,0.0,0.0,0);
    cFinalPaperStyle->SetTicky();
    cFinalPaperStyle->SetTickx();
    cFinalPaperStyle->SetFrameBorderMode(0);

    Set4x2PadPositions(cFinalPaperStyle);
    cFinalPaperStyle->Modified();
    cFinalPaperStyle->Update();
    //    Convert3x2Matrix(cFinalPaperStyle,kTRUE);

  Int_t orderD[4]={0,1,2,3};  
  for(Int_t jp=0;jp<=3;jp++){// First loop for NS yield
    Int_t needTitle=0;
    if(jp==0)needTitle=2;
    gStyle->SetOptStat(0000);
    TPad *pd=(TPad*)cFinalPaperStyle->cd(jp+1);
    SetPadStyle(pd);
    pd->cd();    
    if(jp==1){// identifier set as: 10*quantity+binass; NS --> 0 ; binass == jp (not orderAssoc[jp])
      Compare(orderD[jp],3,pd,100+needTitle);    
    }
    else if(jp==0){// identifier set as: 10*quantity+binass; NS:
      Compare(orderD[jp],3,pd,needTitle);    
      TLatex *tlALICE=GetALICEtext(orderD[jp]);
      tlALICE->Draw();
      TLatex *tlSide=GetTextSide(3,orderD[jp]);
      tlSide->Draw();
    }
    else{
      Compare(orderD[jp],3,pd,0+needTitle);    
    }
    TLatex *tlAssYieldPt=GetDPtText(orderD[jp],orderD[jp],0);
    tlAssYieldPt->Draw();
    if(jp==0) {
      TLatex *tlDeta=GetDRapForSystem(0,orderD[jp],1);
      tlDeta->Draw();
    }

    gStyle->SetOptStat(0000);   
  }

  for(Int_t jp=0;jp<=3;jp++){// second loop --> sigma
    Int_t needTitle=0;
    if(jp==0)needTitle=3;
    else needTitle=1;
    gStyle->SetOptStat(0000);
    TPad *pd=(TPad*)cFinalPaperStyle->cd(jp+5);
    SetPadStyle(pd);
    pd->cd();
    gStyle->SetOptStat(0000);
    Compare(orderD[jp],4,pd,needTitle); 

    gStyle->SetOptStat(0000);   
  }
  for(Int_t j=8;j>=1;j--){
    TPad *pd=(TPad*)cFinalPaperStyle->cd(j);
    pd->Draw();
  }

  TString outdir = Form("%s/ScriptOutput/ReflectedPlots/StdRebin/AllPlots/Averages/ComparisonTopPb",basicdir.Data());
  
  cFinalPaperStyle->Modified();
  cFinalPaperStyle->Update();
  cFinalPaperStyle->SaveAs(Form("%s/CompareFitResults_ppVspPb_5TeV_VsAssPt_AwaySide.root",outdir.Data()));
  cFinalPaperStyle->SaveAs(Form("%s/CompareFitResults_ppVspPb_5TeV_VsAssPt_AwaySide.eps",outdir.Data()));
  cFinalPaperStyle->SaveAs(Form("%s/CompareFitResults_ppVspPb_5TeV_VsAssPt_AwaySide.png",outdir.Data()));
  cFinalPaperStyle->SaveAs(Form("%s/CompareFitResults_ppVspPb_5TeV_VsAssPt_AwaySide.pdf",outdir.Data()));

  return;

}

