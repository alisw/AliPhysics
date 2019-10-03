//gROOT->Reset();
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <Riostream.h>
#include "RooFit.h"
#include "RooVoigtian.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooComplex.h"
#include "RooMath.h"

#include "TVector2.h"
#include "TFile.h"
#include "TString.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "TText.h"
#include "TRandom3.h"
#include "TArray.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMinuit.h"

#include "TCanvas.h"
#include "TLatex.h"
#include "TImage.h"
#include "TGaxis.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"

using namespace std;

bool SaveToFile=kFALSE;
int VariationType=0;// 0 = CutVariations, 1 = FitVariations
//
int CutVariation=0;// 0=Standard Cuts, 1-12 are Systematic deviations
int FitVariation=1;// 1 for VoigtOnly with BkgSub, 2 for VoigtOnly without BkgSub, 3 for RightNorm, 4 Alternate BinCounting Range
TString *outfilename = new TString();
//
// below are default (standard) settings
int POLdegree=1;// 1 with BkgSub, 5 without 
bool XiStarCase=kTRUE;// Standard: kTRUE
bool VoigtFitOnly=kFALSE;// Standard: kFALSE
bool bkgSubtract=kTRUE;// Standard: kTRUE
bool LeftNorm=kTRUE;// Standard: kTRUE
bool MCAssocCase=kTRUE;// Standard: kTRUE
short ParticleCase=0;// 0 for particle+anti-particle. 1 for particle. 2 for anti-particle


const int max_ptbins = 8; // 18 for Xi, 8 for XiStar
double raprange[2]={-0.5,0.5};// for both Xi and XiStar (bin size = 0.1)

double XiyieldRange[2]={1.313, 1.332}; 
double FitRangeXi[2]={1.27,1.36};

double XiStaryieldRange[2]={1.51, 1.56};
double FitRangeXiStar[2]={1.48,1.59};

double XiCountBound[4]={1.304,1.313,1.331,1.340};// -9sigma, -4.5sigma, +4.5sigma, +9sigma 
double XiStarCountBound[6]={1.481, 1.490, 1.523, 1.541, 1.570, 1.579};// 1 gamma range (Standard)
double XiStarCountBound_SysRange[6]={1.481, 1.490, 1.528, 1.537, 1.570, 1.579};// 1/2*gamma range (Systematic Variation)

double yieldRange[2]={0};
bool reject;


double PolFunction(double *x, double *par);
double PolFunctionSpecialXi(double *x, double *par);
double BWFunction(double *x, double *par);
double BWplusPol(double *x, double *par);
double GausFunction(double *x, double *par);
double GausplusPol(double *x, double *par);
double GausplusPolSpecialXi(double *x, double *par);
double myLevyPt(Double_t *x, Double_t *par);
double myExpFit(Double_t *x, Double_t *par);
double Voigt(Double_t *x, Double_t *par);
double VoigtplusPol(double *x, double *par);



void Plot_macroXiStar(){

  if(VariationType==0) {
    outfilename->Append("CutVariation_");
    *outfilename += CutVariation;
  }
  if(VariationType==1) {
    outfilename->Append("FitVariation_");
    *outfilename += FitVariation;
    // 1 for VoigtOnly with BkgSub, 2 for VoigtOnly without BkgSub, 3 for RightNorm, 4 Alternate BinCounting Range
    if(FitVariation==1) {VoigtFitOnly=kTRUE;}
    if(FitVariation==2) {VoigtFitOnly=kTRUE; bkgSubtract=kFALSE; POLdegree=5;}
    if(FitVariation==3) {LeftNorm=kFALSE;}
    if(FitVariation==4) {for(int ii=0; ii<6; ii++) XiStarCountBound[ii]=XiStarCountBound_SysRange[ii];}
  }
  outfilename->Append(".root");

  double FitRange[2]={0};
  if(XiStarCase) {
    FitRange[0] = FitRangeXiStar[0]; FitRange[1] = FitRangeXiStar[1];
    yieldRange[0] = XiStaryieldRange[0]; yieldRange[1] = XiStaryieldRange[1];
  }
  else {
    FitRange[0] = FitRangeXi[0]; FitRange[1] = FitRangeXi[1];
    yieldRange[0] = XiyieldRange[0]; yieldRange[1] = XiyieldRange[1];
  }
  


  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111);
  //gStyle->SetOptFit(1111);
  
  
  TCanvas *can = new TCanvas("can","can",13,34,700,500);
  gStyle->SetOptFit(1111);
  can->Range(-1.25,-0.2625,11.25,2.3625);
  can->SetFillColor(10);
  can->SetBorderMode(0);
  can->SetBorderSize(2);
  //can->SetGridx();
  //can->SetGridy();
  can->SetFrameFillColor(0);
  can->SetFrameBorderMode(0);
  can->SetFrameBorderMode(0);
  can->Divide(4,2);

  TLegend *legend = new TLegend(.1,.7,.35,.9,NULL,"brNDC");
  legend->SetBorderSize(1);
  legend->SetTextSize(.04);// small .03; large .036 
  //legend->SetLineColor(0);
  legend->SetFillColor(0);
  TLegend *legend2 = new TLegend(.35,.62,.9,.72,NULL,"brNDC");
  legend2->SetBorderSize(1);
  legend2->SetTextSize(.04);// small .03; large .036 
  //legend2->SetLineColor(0);
  legend2->SetFillColor(0);



  int PARBINS_temp;
  PARBINS_temp = POLdegree+1+4;
  const int PARBINS=PARBINS_temp;
  int offset=0;
  offset=4;
  
  double bkg_params[POLdegree+1];

  TString *polname=new TString();
  if(POLdegree==1) polname->Append("pol1(0)");
  if(POLdegree==2) polname->Append("pol2(0)");
  if(POLdegree==3) polname->Append("pol3(0)");
  if(POLdegree==4) polname->Append("pol4(0)");
  if(POLdegree==5) polname->Append("pol5(0)");

 
  TF1 *myPol[max_ptbins];
  for(int i=0; i<max_ptbins; i++){
    TString *name = new TString("myPol");
    *name +=i+1;
    myPol[i] = new TF1(name->Data(),polname->Data(),0,2);
    myPol[i]->SetLineColor(3);
    myPol[i]->SetRange(FitRange[0],FitRange[1]);
  }

  
  //
  
  
  // use 10b + 10c + 10d for final result
  //TFile *myfile = new TFile("Results/Real_LHC10b_lego.root","READ");
  //TFile *myfile = new TFile("Results/Real_LHC10c_lego.root","READ");
  //TFile *myfile = new TFile("Results/Real_LHC10d_lego.root","READ");
  TFile *myfile = new TFile("Results/RealMerge_10b_10c_10d_lego.root","READ");// accepted version
  

  // use 10d1 + 10d4 + 10f6a for final result
  //TFile *myfile2 = new TFile("Results/MC_LHC10d1_lego.root","READ");
  //TFile *myfile2 = new TFile("Results/MC_LHC10d4_lego.root","READ");
  //TFile *myfile2 = new TFile("Results/MC_LHC10f6a_lego.root","READ");// Pythia
  TFile *myfile2 = new TFile("Results/MCMerge_10d1_10d4_10f6a_lego.root","READ");// accepted version
  //TFile *myfile2 = new TFile("Results/MC_LHC10f6_lego.root","READ");// Phojet
  
  double David_XiMinus_eff[18]={0.012631, 0.0321763, 0.0473262, 0.059645, 0.0758325, 0.0926316, 0.109662, 0.120194, 0.145389, 0.17972, 0.218773, 0.248567, 0.287159, 0.298501, 0.287736, 0.278715, 0.241559, 0.168857};
  double David_XiMinus_eff_e[18]={0.00038781, 0.000951899, 0.00123674, 0.00147704, 0.0017635, 0.00209177, 0.00243421, 0.00272509, 0.00236159, 0.0030364, 0.00329984, 0.00389785, 0.00508514, 0.00597671, 0.00825827, 0.0124455, 0.0167743, 0.0190201};
  double David_XiPlus_eff[18]={0.010763, 0.0294655, 0.0429606, 0.0558855, 0.0722606, 0.0901825, 0.0982448, 0.114756, 0.140581, 0.172746, 0.208525, 0.247291, 0.28522, 0.29594, 0.28497, 0.276188, 0.249522, 0.210375};
  double David_XiPlus_eff_e[18]={0.000382058, 0.000962961, 0.00122875, 0.00147933, 0.00178481, 0.00212508, 0.00236694, 0.00273767, 0.00238592, 0.00306165, 0.00329165, 0.0039834, 0.00516792, 0.00603728, 0.00836173, 0.0124634, 0.0168862, 0.0220455};
  
  TDirectoryFile *tdir=(TDirectoryFile*)myfile->Get("PWGLF.outputXiStarAnalysis.root");
  TList *List1=(TList*)tdir->Get("XiStarOutput");
  //TList *List1=(TList*)myfile->Get("MyList");
  myfile->Close();
  
  TDirectoryFile *tdir2=(TDirectoryFile*)myfile2->Get("PWGLF.outputXiStarAnalysis.root");
  TList *List2=(TList*)tdir2->Get("XiStarOutput");
  //TList *List2=(TList*)myfile2->Get("MyList");
  myfile2->Close();
  
  TH1F *Events_PhysSel=(TH1F*)(List1->FindObject("fMultDist1"));
  TH1F *Events_postPV=(TH1F*)(List1->FindObject("fMultDist3"));
  TH1F *EventsMC_PhysSel=(TH1F*)(List2->FindObject("fMultDist1"));
  TH1F *EventsMC_postPV=(TH1F*)(List2->FindObject("fMultDist3"));
  
  cout<<"# Events passing PhysSelection Cuts = "<<Events_PhysSel->GetEntries()<<endl;
  cout<<"# Events passing PhysSel + Vz Cut + PileUp Cut = "<<Events_postPV->GetEntries()<<endl;
  cout<<"# MC Events passing Zvertex and PhysSelection Cuts = "<<EventsMC_postPV->GetEntries()<<endl;
  
  TString *name = new TString("fXi_0");
  TH3F *Xi3 = (TH3F*)(List1->FindObject(name->Data()));
  name = new TString("fXibar_0");
  TH3F *Xi_bar3 = (TH3F*)(List1->FindObject(name->Data()));
  
  name = new TString("fXiMinusPiPlus_"); *name += CutVariation;
  TH3F *XiMinus_piPlus3 = (TH3F*)(List1->FindObject(name->Data()));
  name = new TString("fXiMinusPiMinus_"); *name += CutVariation;
  TH3F *XiMinus_piMinus3 = (TH3F*)(List1->FindObject(name->Data()));
  name = new TString("fXiPlusPiPlus_"); *name += CutVariation;
  TH3F *XiPlus_piPlus3 = (TH3F*)(List1->FindObject(name->Data()));
  name = new TString("fXiPlusPiMinus_"); *name += CutVariation;
  TH3F *XiPlus_piMinus3 = (TH3F*)(List1->FindObject(name->Data()));
  
  name = new TString("fXiMinusPiPlusbkg_"); *name += CutVariation;
  TH3F *XiMinus_piPlus3_bkg = (TH3F*)(List1->FindObject(name->Data()));
  name = new TString("fXiMinusPiMinusbkg_"); *name += CutVariation;
  TH3F *XiMinus_piMinus3_bkg = (TH3F*)(List1->FindObject(name->Data()));
  name = new TString("fXiPlusPiPlusbkg_"); *name += CutVariation;
  TH3F *XiPlus_piPlus3_bkg = (TH3F*)(List1->FindObject(name->Data()));
  name = new TString("fXiPlusPiMinusbkg_"); *name += CutVariation;
  TH3F *XiPlus_piMinus3_bkg = (TH3F*)(List1->FindObject(name->Data()));
  
  
  TH3F *MCinput_Xi_prePV = (TH3F*)(List2->FindObject("fMCinputTotalXi1"));
  TH3F *MCinput_Xi_bar_prePV = (TH3F*)(List2->FindObject("fMCinputTotalXibar1"));
  TH3F *MCinput_XiStar_prePV = (TH3F*)(List2->FindObject("fMCinputTotalXiStar1"));
  TH3F *MCinput_XiStar_bar_prePV = (TH3F*)(List2->FindObject("fMCinputTotalXiStarbar1"));
  TH3F *MCinput_XiStar = (TH3F*)(List2->FindObject("fMCinputTotalXiStar3"));
  TH3F *MCinput_XiStar_bar = (TH3F*)(List2->FindObject("fMCinputTotalXiStarbar3"));
  TH3F *MCinput_Xi = (TH3F*)(List2->FindObject("fMCinputTotalXi3"));
  TH3F *MCinput_Xi_bar = (TH3F*)(List2->FindObject("fMCinputTotalXibar3"));

  int binYL = MCinput_XiStar->GetYaxis()->FindBin(raprange[0]+0.01);
  int binYH = MCinput_XiStar->GetYaxis()->FindBin(raprange[1]-0.01);
  int binML = MCinput_XiStar->GetZaxis()->FindBin(1.51);
  int binMH = MCinput_XiStar->GetZaxis()->FindBin(1.55);
  
  TH1D *MCinput_Spectrum;
  TH1D *MCinput_Spectrum_bar;
  if(XiStarCase){
    MCinput_Spectrum= (TH1D*)MCinput_XiStar->ProjectionX("MCinput_Spectrum",binYL,binYH, binML,binMH);
    MCinput_Spectrum_bar = (TH1D*)MCinput_XiStar_bar->ProjectionX("MCinput_Spectrum_bar",binYL,binYH, binML,binMH);
  }else{
    MCinput_Spectrum= (TH1D*)MCinput_Xi->ProjectionX("MCinput_Spectrum",binYL,binYH, binML,binMH);
    MCinput_Spectrum_bar = (TH1D*)MCinput_Xi_bar->ProjectionX("MCinput_Spectrum_bar",binYL,binYH, binML,binMH);
  }
  MCinput_Spectrum->Add(MCinput_Spectrum_bar);

  name = new TString("fMCrecXiMinusPiPlus_"); *name += CutVariation;
  TH3F *MCrec_XiStar = (TH3F*)(List2->FindObject(name->Data()));
  name = new TString("fMCrecXiPlusPiMinus_"); *name += CutVariation;
  TH3F *MCrec_XiStar_bar = (TH3F*)(List2->FindObject(name->Data()));
  
  
  TH3F *MCrec_Xi;
  TH3F *MCrec_Xi_bar;
  if(MCAssocCase){
    MCrec_Xi = (TH3F*)(List2->FindObject("fMCrecXi"));
    MCrec_Xi_bar = (TH3F*)(List2->FindObject("fMCrecXibar"));
  }else{
    name = new TString("fXi_0");
    MCrec_Xi = (TH3F*)(List2->FindObject(name->Data()));
    name = new TString("fXibar_0");
    MCrec_Xi_bar = (TH3F*)(List2->FindObject(name->Data()));
  }


  
  
  double ptedges_Dhevan[9]={0.8, 1.2, 1.6, 2.0, 2.4, 3.2, 4.0, 4.8, 5.6};// standard
  double pt_points_Dhevan[8]={1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.4, 5.2};// standard
  double pt_points_e_Dhevan[8]={.2, .2, .2, .2, .4, .4, .4, .4};// standard
  //double ptedges_Dhevan[17]={0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4.0, 4.4, 4.8, 5.2, 5.6};
  //double pt_points_Dhevan[16]={0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.6, 3.0, 3.4, 3.8, 4.2, 4.6, 5.0, 5.4};
  //double pt_points_e_Dhevan[16]={.1, .1, .1, .1, .1, .1, .1, .1, .2, .2, .2, .2, .2, .2, .2, .2};
  
  double ptedges_David[19]={0.6, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.7, 1.9, 2.2, 2.6, 3.1, 3.9, 4.9, 6, 7.2, 8.5};
  double pt_points_David[18]={0.7, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.6, 1.8, 2.05, 2.4, 2.85, 3.5, 4.4, 5.45, 6.6, 7.85};
  double pt_points_e_David[18]={0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.3, 0.4, 0.5, 0.8, 1, 1.1, 1.2, 1.3};
  
  
  double ptedges[max_ptbins+1+1]={0};
  double pt_points[max_ptbins]={0};
  double pt_points_e[max_ptbins]={0};
  TString *names[max_ptbins];
  
  for(int i=0; i<max_ptbins; i++){
    if(XiStarCase){
      //ptedges[i] = ptedges_Dhevan[i];
      //ptedges[max_ptbins] = ptedges_Dhevan[max_ptbins];
      pt_points[i] = pt_points_Dhevan[i];
      pt_points_e[i] = pt_points_e_Dhevan[i];
      names[i] = new TString("pt bin ");
      *names[i] += i+1;
    }
    else{
      //ptedges[i] = ptedges_David[i];
      //ptedges[max_ptbins] = ptedges_David[max_ptbins];
      pt_points[i] = pt_points_David[i];
      pt_points_e[i] = pt_points_e_David[i]/2.;
      names[i] = new TString();
      *names[i] += float(ptedges_David[i]);
      names[i]->Append(" < pt < ");
      *names[i] += float(ptedges_David[i+1]);
    }
  }
  for(int i=0; i<max_ptbins+1; i++){
    if(XiStarCase){
      if(i==0) ptedges[i]=0;
      else {
	ptedges[i] = ptedges_Dhevan[i-1];
	ptedges[max_ptbins+1] = ptedges_Dhevan[max_ptbins];
      }
    }else {
      if(i==0) ptedges[i]=0;
      else {
	ptedges[i] = ptedges_David[i-1];
	ptedges[max_ptbins+1] = ptedges_David[max_ptbins];
      }
    }
  }

  double Eff_corr_Xiplus[max_ptbins]={0};
  double Eff_corr_Ximinus[max_ptbins]={0};
  
  if(!XiStarCase){
    // David binning
    double geant3fluka_corr_Xiplus[18]={0.8646, 0.8958, 0.9063, 0.9135, 0.9190, 0.9241, 0.9285, 0.9321, 0.9370, 0.9428, 0.9485, 0.9548, 0.9612, 0.9677, 0.9745, 0.9802, 0.9847, 0.9884};
    double geant3fluka_corr_Ximinus[18]={0.9843, 0.9903, 0.9906, 0.9906, 0.9906, 0.9906, 0.9906, 0.9906, 0.9906, 0.9906, 0.9906, 0.9906, 0.9906, 0.9906, 0.9906, 0.9906, 0.9906, 0.9906};
    for(int i=0; i<max_ptbins; i++) {
      Eff_corr_Xiplus[i] = geant3fluka_corr_Xiplus[i];
      Eff_corr_Ximinus[i] = geant3fluka_corr_Ximinus[i];
    }
  }
  if(XiStarCase){
    double geant3fluka_corr_Xiplus[8]={0.9118, 0.9305, 0.9427, 0.9514, 0.9601, 0.9686, 0.9746, 0.9792};// standard
    double geant3fluka_corr_Ximinus[8]={0.9906, 0.9906, 0.9906, 0.9906, 0.9906, 0.9906, 0.9906, 0.9906};// standard
    //double geant3fluka_corr_Xiplus[16]={0.9025, 0.9166, 0.9264, 0.9339, 0.9400, 0.9452, 0.9496, 0.9534, 0.9579, 0.9631, 0.9672, 0.9706, 0.9737, 0.9760, 0.9784, 0.9803};
    //double geant3fluka_corr_Ximinus[16]={0.9905, 0.9906, 0.9906, 0.9906, 0.9906, 0.9906, 0.9906, 0.9906, 0.9906, 0.9906, 0.9906, 0.9906, 0.9906, 0.9906, 0.9906, 0.9906};
    for(int i=0; i<max_ptbins; i++) {
      Eff_corr_Xiplus[i] = geant3fluka_corr_Xiplus[i];
      Eff_corr_Ximinus[i] = geant3fluka_corr_Ximinus[i];
    }
  }
  
  double Eff[max_ptbins]={0}; 
  double Eff_e[max_ptbins]={0};
  
  
  double spectrum[max_ptbins]={0};
  double spectrum_e[max_ptbins]={0};
  double chi2perNDF[max_ptbins]={0};
  double yield[max_ptbins]={0};
  double yield_e[max_ptbins]={0};
  double width[max_ptbins]={0};
  double width_e[max_ptbins]={0};
  double mass[max_ptbins]={0};
  double mass_e[max_ptbins]={0};
  double MCwidth[max_ptbins]={0};
  double MCwidth_e[max_ptbins]={0};
  double MCmass[max_ptbins]={0};
  double MCmass_e[max_ptbins]={0};
  double MissedYield[max_ptbins]={0};
  double MissedYield_e[max_ptbins]={0};

  TH1F *MCrec_1d[max_ptbins];
  TH1F *MCrec_1d_bar[max_ptbins];
  TF1 *MCfit = new TF1("MCfit",Voigt,FitRange[0],FitRange[1],4);
  double LostRatio=0, LostRatio_e=0;
  
  for(int ptbin=0; ptbin<max_ptbins; ptbin++){
    can->cd(ptbin+1);
    cout<<"+++++++++++++++ ptbin = "<<ptbin<<"  +++++++++++++++++"<<endl;
    int rapbinFirst = Xi3->GetYaxis()->FindBin(raprange[0]+0.01);// bin 16
    int rapbinLast = Xi3->GetYaxis()->FindBin(raprange[1]-0.01);// bin 25
    int ptbinFirst = Xi3->GetXaxis()->FindBin(ptedges[ptbin+1]+0.01);
    int ptbinLast = Xi3->GetXaxis()->FindBin(ptedges[ptbin+1+1]-0.01);
    
    ////////////////////////////
    // Efficiency calcualtion
    double rec_num=0, input_num=0;
    double rec_num1=0, input_num1=0;
    double rec_num2=0, input_num2=0;
    if(XiStarCase){
      int massbinFirst = 1;
      int massbinLast = MCrec_XiStar->GetNbinsZ();
      rec_num1 = MCrec_XiStar->Integral(ptbinFirst,ptbinLast, rapbinFirst,rapbinLast, massbinFirst,massbinLast);
      input_num1 = MCinput_XiStar->Integral(ptbinFirst,ptbinLast, rapbinFirst,rapbinLast, massbinFirst,massbinLast);
      rec_num2 = MCrec_XiStar_bar->Integral(ptbinFirst,ptbinLast, rapbinFirst,rapbinLast, massbinFirst,massbinLast);
      input_num2 = MCinput_XiStar_bar->Integral(ptbinFirst,ptbinLast, rapbinFirst,rapbinLast, massbinFirst,massbinLast);
    
      // set ranges for projections later on
      MCrec_XiStar->GetXaxis()->SetRange(ptbinFirst,ptbinLast);
      MCrec_XiStar->GetYaxis()->SetRange(rapbinFirst,rapbinLast);
      MCrec_XiStar->GetZaxis()->SetRange(massbinFirst,massbinLast);
      MCrec_XiStar_bar->GetXaxis()->SetRange(ptbinFirst,ptbinLast);
      MCrec_XiStar_bar->GetYaxis()->SetRange(rapbinFirst,rapbinLast);
      MCrec_XiStar_bar->GetZaxis()->SetRange(massbinFirst,massbinLast);
      // Estimate Lost number of Xi* from PV cut 
      double num=0, den=0;
      if(ParticleCase==0){
	num = MCinput_XiStar->Integral(ptbinFirst,ptbinLast, rapbinFirst,rapbinLast, massbinFirst,massbinLast) + MCinput_XiStar_bar->Integral(ptbinFirst,ptbinLast, rapbinFirst,rapbinLast, massbinFirst,massbinLast);
	den = MCinput_XiStar_prePV->Integral(ptbinFirst,ptbinLast, rapbinFirst,rapbinLast, massbinFirst,massbinLast) + MCinput_XiStar_bar_prePV->Integral(ptbinFirst,ptbinLast, rapbinFirst,rapbinLast, massbinFirst,massbinLast);
      }else if(ParticleCase==1){
	num = MCinput_XiStar->Integral(ptbinFirst,ptbinLast, rapbinFirst,rapbinLast, massbinFirst,massbinLast);
	den = MCinput_XiStar_prePV->Integral(ptbinFirst,ptbinLast, rapbinFirst,rapbinLast, massbinFirst,massbinLast);
      }else {
	num = MCinput_XiStar_bar->Integral(ptbinFirst,ptbinLast, rapbinFirst,rapbinLast, massbinFirst,massbinLast);
	den = MCinput_XiStar_bar_prePV->Integral(ptbinFirst,ptbinLast, rapbinFirst,rapbinLast, massbinFirst,massbinLast);
      }
      LostRatio = num/den;
      LostRatio_e = sqrt(pow(sqrt(num)/den,2)+pow(sqrt(den)*num/pow(den,2),2));
      
      ///////////////////////////////////////
      // MC fitting for mass and sigma values
      TString *nameMCrec=new TString("MCrec_1d_");
      TString *nameMCrecbar=new TString("MCrecbar_1d_");
      *nameMCrec += ptbin;
      *nameMCrecbar += ptbin;
      MCrec_1d[ptbin] = (TH1F*)MCrec_XiStar->ProjectionZ(nameMCrec->Data(),ptbinFirst,ptbinLast, rapbinFirst,rapbinLast);
      MCrec_1d_bar[ptbin] = (TH1F*)MCrec_XiStar->ProjectionZ(nameMCrecbar->Data(),ptbinFirst,ptbinLast, rapbinFirst,rapbinLast);
      MCrec_1d[ptbin]->Add(MCrec_1d_bar[ptbin]);
      MCrec_1d[ptbin]->GetXaxis()->SetRangeUser(1.46, 1.62);
      MCrec_1d[ptbin]->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
      MCrec_1d[ptbin]->SetTitle("p_{T} bin 8");
      MCfit->SetParameter(0,.5);// mean
      MCfit->SetParameter(1,1.532);// mean
      MCfit->SetParameter(2,.002);// Gaussian sigma
      MCfit->SetParLimits(0,0.01,10.);// 0.0001 to 0.01 for BW only; 0.01 to 10 for true Voigtian
      MCfit->SetParLimits(1,1.52,1.54);
      MCfit->SetParLimits(2,0.00001,.01);
      MCfit->FixParameter(3,.0091);// BW width
      MCrec_1d[ptbin]->Fit(MCfit,"IMQ","",FitRange[0],FitRange[1]);
      MCwidth[ptbin] = MCfit->GetParameter(2);
      MCwidth_e[ptbin] = MCfit->GetParError(2);
      MCmass[ptbin] = MCfit->GetParameter(1);
      MCmass_e[ptbin] = MCfit->GetParError(1);
      //continue;// use this only to plot MC yields one by one
      ///////////////////////////////////////////

    }else {// Xi
      rec_num1 = MCrec_Xi->Integral(ptbinFirst,ptbinLast, rapbinFirst,rapbinLast, 1,MCrec_Xi->GetZaxis()->GetNbins());
      input_num1 = MCinput_Xi->Integral(ptbinFirst,ptbinLast, rapbinFirst,rapbinLast, 1,MCinput_Xi->GetZaxis()->GetNbins());
      rec_num2 = MCrec_Xi_bar->Integral(ptbinFirst,ptbinLast, rapbinFirst,rapbinLast, 1,MCrec_Xi_bar->GetZaxis()->GetNbins());
      input_num2 = MCinput_Xi_bar->Integral(ptbinFirst,ptbinLast, rapbinFirst,rapbinLast, 1,MCinput_Xi_bar->GetZaxis()->GetNbins());
      //
      MCrec_Xi->GetXaxis()->SetRange(ptbinFirst,ptbinLast);
      MCrec_Xi->GetYaxis()->SetRange(rapbinFirst,rapbinLast);
      MCrec_Xi_bar->GetXaxis()->SetRange(ptbinFirst,ptbinLast);
      MCrec_Xi_bar->GetYaxis()->SetRange(rapbinFirst,rapbinLast);
      MCrec_1d[ptbin]=(TH1F*)MCrec_Xi->Project3D("z");
      MCrec_1d_bar[ptbin]=(TH1F*)MCrec_Xi_bar->Project3D("z");
      
      rec_num1 = MCrec_1d[ptbin]->Integral(MCrec_1d[ptbin]->GetXaxis()->FindBin(XiCountBound[1]),MCrec_1d[ptbin]->GetXaxis()->FindBin(XiCountBound[2]-.0001));
      rec_num1 -= MCrec_1d[ptbin]->Integral(MCrec_1d[ptbin]->GetXaxis()->FindBin(XiCountBound[0]),MCrec_1d[ptbin]->GetXaxis()->FindBin(XiCountBound[1]-.0001));
      rec_num1 -= MCrec_1d[ptbin]->Integral(MCrec_1d[ptbin]->GetXaxis()->FindBin(XiCountBound[2]),MCrec_1d[ptbin]->GetXaxis()->FindBin(XiCountBound[3]-.0001));
      
      rec_num2 = MCrec_1d_bar[ptbin]->Integral(MCrec_1d_bar[ptbin]->GetXaxis()->FindBin(XiCountBound[1]),MCrec_1d_bar[ptbin]->GetXaxis()->FindBin(XiCountBound[2]-.0001));
      rec_num2 -= MCrec_1d_bar[ptbin]->Integral(MCrec_1d_bar[ptbin]->GetXaxis()->FindBin(XiCountBound[0]),MCrec_1d_bar[ptbin]->GetXaxis()->FindBin(XiCountBound[1]-.0001));
      rec_num2 -= MCrec_1d_bar[ptbin]->Integral(MCrec_1d_bar[ptbin]->GetXaxis()->FindBin(XiCountBound[2]),MCrec_1d_bar[ptbin]->GetXaxis()->FindBin(XiCountBound[3]-.0001));

      // Estimate Lost number of Xi from PV cut 
      if(ParticleCase==0){
	LostRatio = MCinput_Xi->Integral(ptbinFirst,ptbinLast, rapbinFirst,rapbinLast, 1,MCrec_Xi->GetZaxis()->GetNbins()) + MCinput_Xi_bar->Integral(ptbinFirst,ptbinLast, rapbinFirst,rapbinLast, 1,MCrec_Xi->GetZaxis()->GetNbins());
	LostRatio /= MCinput_Xi_prePV->Integral(ptbinFirst,ptbinLast, rapbinFirst,rapbinLast, 1,MCrec_Xi->GetZaxis()->GetNbins()) + MCinput_Xi_bar_prePV->Integral(ptbinFirst,ptbinLast, rapbinFirst,rapbinLast, 1,MCrec_Xi->GetZaxis()->GetNbins());
      }else if(ParticleCase==1){
	LostRatio = MCinput_Xi->Integral(ptbinFirst,ptbinLast, rapbinFirst,rapbinLast, 1,MCrec_Xi->GetZaxis()->GetNbins());
	LostRatio /= MCinput_Xi_prePV->Integral(ptbinFirst,ptbinLast, rapbinFirst,rapbinLast, 1,MCrec_Xi->GetZaxis()->GetNbins());
      }else {
	LostRatio = MCinput_Xi_bar->Integral(ptbinFirst,ptbinLast, rapbinFirst,rapbinLast, 1,MCrec_Xi->GetZaxis()->GetNbins());
	LostRatio /= MCinput_Xi_bar_prePV->Integral(ptbinFirst,ptbinLast, rapbinFirst,rapbinLast, 1,MCrec_Xi->GetZaxis()->GetNbins());
      }
      

    }// else Xi case
    cout<<"LostRatio = "<<LostRatio<<" +- "<<LostRatio_e<<endl;
    
    
    // Apply Geant3/fluka correction
    rec_num1 /= Eff_corr_Ximinus[ptbin];
    rec_num2 /= Eff_corr_Xiplus[ptbin];
    

    if(ParticleCase == 0){
      rec_num = rec_num1 + rec_num2;
      input_num = input_num1 + input_num2;
    }
    if(ParticleCase == 1){
      rec_num = rec_num1;
      input_num = input_num1;
    }
    if(ParticleCase == 2){
      rec_num = rec_num2;
      input_num = input_num2;
    }
    

    Eff[ptbin] = rec_num/input_num;
    Eff_e[ptbin] = sqrt( pow(sqrt(rec_num)/input_num,2) + pow(sqrt(input_num)*rec_num/pow(input_num,2),2));
    
    cout<<"Efficiency = "<<Eff[ptbin]<<"  +- "<<Eff_e[ptbin]<<"   MC_rec_count = "<<rec_num<<endl;
    //cout<<"Reassigning Eff in last pt bin to second to last!!!!!!!!!!!!!"<<endl; Eff[7] = Eff[6];

   

    //cout<<"Re-assigning to David's Eff"<<endl; Eff[ptbin] = David_XiMinus_eff[ptbin];
    
        
    reject=kTRUE;
    TF1 *myBkg=new TF1("myBkg",PolFunction,FitRange[0],FitRange[1],POLdegree+1);
    //TF1 *myBkg=new TF1("myBkg",PolFunctionSpecialXi,1.27,1.36,POLdegree+1);// special bkg for Xi
    myBkg->SetParameters(6,0);
    myBkg->SetLineColor(3);
    
   
    Xi3->GetXaxis()->SetRange(ptbinFirst,ptbinLast);// pt bins
    Xi_bar3->GetXaxis()->SetRange(ptbinFirst,ptbinLast);
    Xi3->GetYaxis()->SetRange(rapbinFirst,rapbinLast);
    Xi_bar3->GetYaxis()->SetRange(rapbinFirst,rapbinLast);
    
    XiMinus_piPlus3->GetXaxis()->SetRange(ptbinFirst,ptbinLast);
    XiPlus_piMinus3->GetXaxis()->SetRange(ptbinFirst,ptbinLast);
    XiMinus_piMinus3->GetXaxis()->SetRange(ptbinFirst,ptbinLast);
    XiPlus_piPlus3->GetXaxis()->SetRange(ptbinFirst,ptbinLast);
    
    XiMinus_piPlus3->GetYaxis()->SetRange(rapbinFirst,rapbinLast);
    XiPlus_piMinus3->GetYaxis()->SetRange(rapbinFirst,rapbinLast);
    XiMinus_piMinus3->GetYaxis()->SetRange(rapbinFirst,rapbinLast);
    XiPlus_piPlus3->GetYaxis()->SetRange(rapbinFirst,rapbinLast);
    
    
    XiMinus_piPlus3_bkg->GetXaxis()->SetRange(ptbinFirst,ptbinLast);
    XiPlus_piMinus3_bkg->GetXaxis()->SetRange(ptbinFirst,ptbinLast);
    XiMinus_piMinus3_bkg->GetXaxis()->SetRange(ptbinFirst,ptbinLast);
    XiPlus_piPlus3_bkg->GetXaxis()->SetRange(ptbinFirst,ptbinLast);

    XiMinus_piPlus3_bkg->GetYaxis()->SetRange(rapbinFirst,rapbinLast);
    XiPlus_piMinus3_bkg->GetYaxis()->SetRange(rapbinFirst,rapbinLast);
    XiMinus_piMinus3_bkg->GetYaxis()->SetRange(rapbinFirst,rapbinLast);
    XiPlus_piPlus3_bkg->GetYaxis()->SetRange(rapbinFirst,rapbinLast);

  
    TH1F *Xi=(TH1F*)Xi3->Project3D("z");
    TH1F *Xibar=(TH1F*)Xi_bar3->Project3D("z");
    if(ParticleCase==0) Xi->Add(Xibar);// to add particle and anti-particle
 
    TH1F *XiMinus_piPlus = (TH1F*)XiMinus_piPlus3->Project3D("z");
    TH1F *XiPlus_piMinus = (TH1F*)XiPlus_piMinus3->Project3D("z");
    TH1F *XiMinus_piMinus = (TH1F*)XiMinus_piMinus3->Project3D("z");
    TH1F *XiPlus_piPlus = (TH1F*)XiPlus_piPlus3->Project3D("z");
    if(ParticleCase==0) {
      XiMinus_piPlus->Add(XiPlus_piMinus);
      XiMinus_piMinus->Add(XiPlus_piPlus);
    }
    

    TH1F *XiMinus_piPlus_bkg = (TH1F*)XiMinus_piPlus3_bkg->Project3D("z");
    TH1F *XiPlus_piMinus_bkg = (TH1F*)XiPlus_piMinus3_bkg->Project3D("z");
    TH1F *XiMinus_piMinus_bkg = (TH1F*)XiMinus_piMinus3_bkg->Project3D("z");
    TH1F *XiPlus_piPlus_bkg = (TH1F*)XiPlus_piPlus3_bkg->Project3D("z");
    if(ParticleCase==0) {
      XiMinus_piPlus_bkg->Add(XiPlus_piMinus_bkg); 
      XiMinus_piMinus_bkg->Add(XiPlus_piPlus_bkg);
    }

    TH1F *os=0x0;
    TH1F *os_bkg=0x0;
    if(XiStarCase) {
      if(ParticleCase==0 || ParticleCase==1) {
	os=(TH1F*)XiMinus_piPlus->Clone();
	os_bkg=(TH1F*)XiMinus_piPlus_bkg->Clone();
      }
      if(ParticleCase==2) os=(TH1F*)XiPlus_piMinus->Clone();
    }else {
      if(ParticleCase==0 || ParticleCase==1) os=(TH1F*)Xi->Clone();
      if(ParticleCase==2) os=(TH1F*)Xibar->Clone();
    }
    
  
    os->Sumw2();
    TH1F *ss=(TH1F*)XiMinus_piMinus->Clone();
    TH1F *ss_bkg=(TH1F*)XiMinus_piMinus_bkg->Clone();
    ss_bkg->SetLineColor(2);
    
    ss->SetLineColor(4);
    ss->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
    ss->GetYaxis()->SetTitle("Counts");
    //ss->SetTitle(name[ptbin]);
    ss->SetTitle(names[ptbin]->Data());
    
    os->SetLineColor(4);
    os->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
    os->GetYaxis()->SetTitle("Counts");
    //os->SetTitle(name[ptbin]);
    os->SetTitle(names[ptbin]->Data());
  
    
    double scalefactor=1;
    if(XiStarCase) {
      os_bkg->SetLineColor(2);
      os_bkg->Sumw2();
      double leftS=0, rightS=0;
      if(LeftNorm) {leftS=1.49; rightS=1.51;}
      else {leftS=1.56; rightS=1.58;}
     
      scalefactor = os->Integral(os->GetXaxis()->FindBin(leftS),os->GetXaxis()->FindBin(rightS));
      scalefactor /= os_bkg->Integral(os_bkg->GetXaxis()->FindBin(leftS),os_bkg->GetXaxis()->FindBin(rightS));
      os_bkg->Scale(scalefactor);
      if(bkgSubtract) os->Add(os_bkg,-1);
    }
    
   
    /////////////////////////////////////
    // 1st Iteration
    reject=kTRUE;
    if(bkgSubtract) {
      for(int polbin=0; polbin<POLdegree+1-2; polbin++) myBkg->FixParameter(polbin+2,0);
      myBkg->SetParameter(0,0);
      myBkg->SetParameter(1,0);
    }
    /////////////
    os->Fit(myBkg,"NQ","",FitRange[0],FitRange[1]);// 1st fit
    /////////////

   
    //ss->Draw();
    //ss_bkg->Scale(ss->Integral()/ss_bkg->Integral());
    //ss_bkg->Draw("same");
  
    TF1 *FullFit_t2;
    TF1 *FullFit_t6;
    
  
    
    // Voigt
    FullFit_t2=new TF1("FullFit_t2",VoigtplusPol,FitRange[0],FitRange[1],PARBINS+1);
    FullFit_t2->SetParameter(0,1.0);// mean
    FullFit_t2->SetParLimits(0,0.01,10.);// 0.0001 to 0.01 for BW only; 0.01 to 10 for true Voigtian
    FullFit_t2->SetParameter(1,1.532);// mean
    FullFit_t2->SetParLimits(1,1.5,1.56);
    //FullFit_t2->FixParameter(2,0);// Gaussian sigma
    FullFit_t2->SetParameter(2,.003);// Gaussian sigma
    FullFit_t2->SetParLimits(2,0.0001,.015);
    FullFit_t2->FixParameter(3,.0091);// BW width
    //FullFit_t2->SetParameter(3,.0091);// BW width
    //FullFit_t2->SetParLimits(3,.001,.015);
    FullFit_t2->SetParName(0,"Norm");
    FullFit_t2->SetParName(1,"Mass");
    FullFit_t2->SetParName(2,"#sigma");
    FullFit_t2->SetParName(3,"#Gamma");
    
    FullFit_t2->SetParameter(4,0);
    FullFit_t2->SetParameter(5,0);
    FullFit_t2->SetParameter(6,0);
    FullFit_t2->SetParameter(7,0);
    FullFit_t2->SetParameter(8,0);
    FullFit_t2->SetParameter(9,0);
    
  
  
       
    if(XiStarCase) os->GetXaxis()->SetRange(1000*(FitRangeXiStar[0]-1.4),1000*(FitRangeXiStar[1]-1.4));
    else os->GetXaxis()->SetRange(1000*(FitRangeXi[0]-1.2),1000*(FitRangeXi[1]+-1.2));
   
    /////////////////////////////////////////////
    // 2nd interation 
    reject=kFALSE;
    double pars[PARBINS];
    double pars_e[PARBINS];

    for(int polbin=0; polbin<PARBINS-offset; polbin++) bkg_params[polbin] = myBkg->GetParameter(polbin);
    
    
    for(int polbin=0; polbin<PARBINS-offset; polbin++) FullFit_t2->FixParameter(offset+polbin, myBkg->GetParameter(polbin));
    //////////////////
    os->Fit(FullFit_t2,"IMQ","",FitRange[0],FitRange[1]);// 2nd fit
    //////////////////
    for(int polbin=0; polbin<PARBINS; polbin++) pars[polbin] = FullFit_t2->GetParameter(polbin);
    ////////////////////////////////////////

   
  
    ///////////////////////////////////////
    // 3rd iteration
    FullFit_t6=new TF1("FullFit_t6",VoigtplusPol,FitRange[0],FitRange[1],PARBINS);
    FullFit_t6->SetParLimits(0,0.01,10.);// 0.0001 to 0.01 for BW only; 0.01 to 10 for true Voigtian
    FullFit_t6->SetParLimits(1,1.5,1.56);
    FullFit_t6->SetParLimits(2,0.0001,.015);
    FullFit_t6->FixParameter(3,.0091);// BW width
    FullFit_t6->SetParName(0,"Norm");
    FullFit_t6->SetParName(1,"Mass");
    FullFit_t6->SetParName(2,"#sigma");
    FullFit_t6->SetParName(3,"#Gamma");
        
    
    for(int polbin=0; polbin<PARBINS; polbin++) FullFit_t6->SetParameter(polbin, pars[polbin]);
    if(bkgSubtract) {
      for(int polbin=0; polbin<PARBINS-offset-2; polbin++) FullFit_t6->FixParameter(polbin+offset+2,0);
      FullFit_t6->SetParameter(offset,0);
      FullFit_t6->SetParameter(offset+1,0);
    }
    //////////////////
    os->Fit(FullFit_t6,"IMEQ","",FitRange[0],FitRange[1]);// 3rd fit
    //////////////////
    for(int polbin=0; polbin<PARBINS-offset; polbin++) bkg_params[polbin] = FullFit_t6->GetParameter(polbin+offset);
    for(int polbin=0; polbin<PARBINS; polbin++) pars[polbin] = FullFit_t6->GetParameter(polbin);
    for(int polbin=0; polbin<PARBINS; polbin++) pars_e[polbin] = FullFit_t6->GetParError(polbin);
  
    

    /////////////////////////////////////////
    if(!bkgSubtract) os_bkg->Draw("same");// XiStar
    for(int polbin=0; polbin<POLdegree+1; polbin++) myBkg->SetParameter(polbin, bkg_params[polbin]);
    
    for(int polbin=0; polbin<POLdegree+1; polbin++) myPol[ptbin]->FixParameter(polbin, bkg_params[polbin]);
    myPol[ptbin]->DrawCopy("same");
    if((ptbin+1)==max_ptbins){
      legend->AddEntry(os,"opp-charge","l");
      legend->AddEntry(os_bkg,"opp-charge event mixing","l");
      //legend->Draw("same");
    }
    
   
    double temp_yield=0, temp_yield_e=0;
    MissedYield[ptbin] = 1000.*(FullFit_t6->Integral(0.,XiStarCountBound[2]) + FullFit_t6->Integral(XiStarCountBound[3],100.) - myBkg->Integral(0.,XiStarCountBound[2]) - myBkg->Integral(XiStarCountBound[3],100.));
    MissedYield_e[ptbin] = MissedYield[ptbin]*FullFit_t6->GetParError(0)/FullFit_t6->GetParameter(0);
    //cout<<"Signal/Bkg = "<<temp_yield/(1000.*myBkg->Integral(XiStarYieldCountRange[0],XiStarYieldCountRange[1]))<<endl;
    //cout<<"Significance = "<<temp_yield/sqrt(temp_yield)<<endl;

    if(VoigtFitOnly) {
      if(!bkgSubtract) temp_yield = 1000.*(FullFit_t6->Integral(1.48,1.59) - myBkg->Integral(1.48,1.59));// must have fit range limits for high degree polynomial bkg
      else temp_yield = 1000.*(FullFit_t6->Integral(0.,100.) - myBkg->Integral(0.,100.));// yield from Voigt fit (1.48-1.59 is function fit range)
      temp_yield_e = temp_yield*FullFit_t6->GetParError(0)/FullFit_t6->GetParameter(0);
    }else {// yield from bin counting
      if(XiStarCase){
	temp_yield = os->Integral(os->GetXaxis()->FindBin(XiStarCountBound[2]),os->GetXaxis()->FindBin(XiStarCountBound[3]-.0001));
	temp_yield -= myBkg->Integral(XiStarCountBound[2],XiStarCountBound[3]);
	cout<<"Missed yield fraction = "<<MissedYield[ptbin]/(temp_yield+MissedYield[ptbin])<<endl;
	temp_yield += MissedYield[ptbin];
	for(int massBin=os->GetXaxis()->FindBin(XiStarCountBound[2]); massBin<=os->GetXaxis()->FindBin(XiStarCountBound[3]); massBin++){
	  temp_yield_e += pow(os->GetBinError(massBin),2);
	}
	temp_yield_e += pow(MissedYield_e[ptbin],2);
	temp_yield_e = sqrt(temp_yield_e);
	cout<<"Included/Total = "<<(temp_yield-MissedYield[ptbin])/temp_yield<<endl;
      }else {// Xi case
	temp_yield = os->Integral(os->GetXaxis()->FindBin(XiCountBound[1]),os->GetXaxis()->FindBin(XiCountBound[2]-.0001));
	temp_yield -= os->Integral(os->GetXaxis()->FindBin(XiCountBound[0]),os->GetXaxis()->FindBin(XiCountBound[1]-.0001)); 
	temp_yield -= os->Integral(os->GetXaxis()->FindBin(XiCountBound[2]),os->GetXaxis()->FindBin(XiCountBound[3]-.0001));
      }
    }
    cout<<"yield = "<<temp_yield<<" +- "<<temp_yield_e<<endl;
       
    
    
    double base = Eff[ptbin]*LostRatio*(raprange[1]-raprange[0])*(2*pt_points_e[ptbin])*Events_PhysSel->GetEntries()/0.852;
    
    spectrum[ptbin] = temp_yield;
    spectrum[ptbin] /= base;
    //
    spectrum_e[ptbin] = pow(temp_yield_e/base,2);
    spectrum_e[ptbin] += pow(temp_yield/base*Eff_e[ptbin]/Eff[ptbin],2);
    spectrum_e[ptbin] += pow(temp_yield/base*LostRatio_e/LostRatio,2);
    spectrum_e[ptbin] += pow(temp_yield/base*sqrt(Events_PhysSel->GetEntries())/Events_PhysSel->GetEntries(),2);
    spectrum_e[ptbin] = sqrt(spectrum_e[ptbin]);
    //spectrum_e[ptbin] = temp_yield/base*Eff_e[ptbin]/Eff[ptbin];// for eff errors only
    if(ParticleCase==0) spectrum[ptbin] /= 2.;// (particle + cc)/2
    if(ParticleCase==0) spectrum_e[ptbin] /= 2.;// (particle + cc)/2
    //
   
    
    cout<<"spectrum point = "<<spectrum[ptbin]<<"  +- "<<spectrum_e[ptbin]<<endl;
    
    yield[ptbin] = temp_yield;
    yield_e[ptbin] = temp_yield_e;
    mass[ptbin] = pars[1]; 
    mass_e[ptbin] = pars_e[1];
    width[ptbin] = pars[2]; 
    width_e[ptbin] = pars_e[2];
    
    chi2perNDF[ptbin] = FullFit_t6->GetChisquare()/double(FullFit_t6->GetNDF());
    
  }
  
  
  
  TH1D *RawYields = new TH1D("RawYields","",8,0.5,8.5);
  for(int ii=0; ii<8; ii++){ RawYields->SetBinContent(ii+1, yield[ii]); RawYields->SetBinError(ii+1, yield_e[ii]);}
  TH1D *Efficiency = new TH1D("Efficiency","",8,0.5,8.5);
  for(int ii=0; ii<8; ii++){ Efficiency->SetBinContent(ii+1, Eff[ii]); Efficiency->SetBinError(ii+1, Eff_e[ii]);}
  
  /*
  TLatex *tex2 = new TLatex(1.6,1.5,"#Xi(1530)^{0} #rightarrow #Xi^{+-} + #pi^{-+}");
  tex2->SetTextSize(.06);
  tex2->Draw();

  TLatex *tex3 = new TLatex(1.6,1.5,"ALICE Performance");
  tex3->SetTextSize(.06);
  tex3->SetTextColor(2);
  tex3->Draw();

  TLatex *tex4 = new TLatex(1.6,1.5,"12/09/2011");
  tex4->SetTextSize(.04);
  tex4->Draw();

  TLatex *tex5 = new TLatex(1.6,1.5,"pp, #sqrt{s} = 7 TeV");
  tex5->SetTextSize(.06);
  tex5->Draw();

  TPad *c1_1 = new TPad("c1_1", "c1_1",0.25,0.3,0.4,0.5);
  c1_1->SetBorderMode(0);
  c1_1->SetLineColor(0);
  c1_1->SetFillColor(0);
  c1_1->Draw();
  
  c1_1->cd();
  TImage *alicelogo = TImage::Open("../../Pictures/alice_logo_transparent2.png");
  //alicelogo->SetConstRatio(kTRUE);
  //alicelogo->SetImageQuality(TAttImage::kImgBest);
  alicelogo->Draw("same");
  */

  
  TCanvas *can2 = new TCanvas("can2","can2",13,34,700,500);
  gStyle->SetOptFit(111);
  can2->Range(-1.25,-0.2625,11.25,2.3625);
  can2->SetFillColor(10);
  can2->SetBorderMode(0);
  can2->SetBorderSize(2);
  can2->SetGridx();
  can2->SetGridy();
  //can2->SetLogy();
  can2->SetFrameFillColor(0);
  can2->SetFrameBorderMode(0);
  can2->SetFrameBorderMode(0);
  
 
  
  TH1D *h_Xispectrum = new TH1D("h_Xispectrum","Xi spectrum",max_ptbins+1,ptedges);
  h_Xispectrum->SetMarkerStyle(20);
  h_Xispectrum->SetMarkerSize(1.0);
  if(XiStarCase) h_Xispectrum->SetTitle("#Xi(1530) spectrum");
  else h_Xispectrum->SetTitle("#Xi spectrum");
  h_Xispectrum->SetTitle("pp #sqrt{s}= 7 TeV");
  h_Xispectrum->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_Xispectrum->GetYaxis()->SetTitle("1/N_{E}d^{2}N/dydp_{T} (GeV/c)^{-1}");
  //h_Xispectrum->GetYaxis()->SetTitle("(Real N_Xi-)/(MC N_Xi-)");
  h_Xispectrum->SetMinimum(3.0e-6);
  h_Xispectrum->SetMaximum(0.02);
  h_Xispectrum->GetXaxis()->SetLimits(0,5.6);

  h_Xispectrum->SetBinContent(1,+100);// first bin is just there to visualize pt=0
  for(int i=0; i<max_ptbins; i++) {    
    h_Xispectrum->SetBinContent(i+1+1, spectrum[i]);
    h_Xispectrum->SetBinError(i+1+1, spectrum_e[i]);
  }
  
  //h_Xispectrum->Draw();
    

  double minFitpoint=0, maxFitpoint=0;
  if(XiStarCase) {minFitpoint = 0.8; maxFitpoint = 5.6;}// 0.8 to 5.6
  else {minFitpoint = 0.6; maxFitpoint = 8.5;}

  
  TF1 *myLevy=new TF1("myLevy",myLevyPt,0,8.5,3);
  myLevy->SetParName(0,"dN/dy");
  myLevy->SetParName(1,"C");
  myLevy->SetParName(2,"n");
  myLevy->SetParameter(0,.008);
  myLevy->SetParameter(1,.3);
  myLevy->SetParameter(2,15);
  myLevy->SetParLimits(0,.001,.2);
  myLevy->SetParLimits(1,.1,1);
  myLevy->SetParLimits(2,1,500);
  
  
  h_Xispectrum->Draw();
  h_Xispectrum->Fit(myLevy,"IME","",minFitpoint,maxFitpoint);
  //h_Xispectrum->Fit(myExp,"IME","",minFitpoint,maxFitpoint);
  //myLevy->Draw("same");
  //TH1D *temp = (TH1D*)myLevy->GetHistogram();

  //legend2->AddEntry(h_Xispectrum,"#Xi^{+}, This analysis","p");
  //legend2->AddEntry(h_Xispectrum,"10b+10c+10d","p");
  legend2->AddEntry(h_Xispectrum,"(#Xi(1530) + cc)/2","p");
  
  cout<<"dN/dy = "<<myLevy->Integral(0,10)<<"  , Covered range = "<<myLevy->Integral(0.8,5.6)<<endl;
 

  double Data_Levy_Ratio[8]={0};
  double Data_Levy_Ratio_e[8]={0};
  for(int i=0; i<8; i++){
    Data_Levy_Ratio[i] = (ptedges[i+2]-ptedges[i+1])*spectrum[i]/myLevy->Integral(ptedges[i+1],ptedges[i+2]);
    Data_Levy_Ratio_e[i] = (ptedges[i+2]-ptedges[i+1])*spectrum_e[i]/myLevy->Integral(ptedges[i+1],ptedges[i+2]);
  }
  TGraphErrors *gr_Data_Levy_Ratio = new TGraphErrors(8, pt_points, Data_Levy_Ratio, pt_points_e, Data_Levy_Ratio_e);
  gr_Data_Levy_Ratio->SetMarkerStyle(20);
  gr_Data_Levy_Ratio->SetTitle("");
  gr_Data_Levy_Ratio->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  gr_Data_Levy_Ratio->GetYaxis()->SetTitle("(Measured Yield)/(Levy fit integral)");
  //gr_Data_Levy_Ratio->Draw("AP");

  // Reference Spectra
  double LHC10b_XiStar[8]={0.00175909, 0.000893716, 0.000957646, 0.000399652, 0.000217706, 6.86327e-05, 3.24599e-05, 1.28201e-05};
  double LHC10b_XiStar_e[8]={0.000460783, 0.000152526, 0.000136327, 6.06295e-05, 2.79773e-05, 1.2783e-05, 9.14534e-06, 4.461e-06};
  double LHC10c_XiStar[8]={0.00153422, 0.00106629, 0.000752749, 0.000495777, 0.000255897, 6.88883e-05, 2.41039e-05, 1.09414e-05};
  double LHC10c_XiStar_e[8]={0.000262538, 0.000106795, 6.6372e-05, 4.51931e-05, 2.07905e-05, 8.01244e-06, 4.03196e-06, 2.75332e-06};
  double LHC10d_XiStar[8]={0.00137741, 0.00118096, 0.000813771, 0.000465871, 0.000223755, 6.5902e-05, 2.2662e-05, 7.11972e-06};
  double LHC10d_XiStar_e[8]={7.81918e-05, 4.79133e-05, 2.92694e-05, 1.16203e-05, 5.00591e-06, 2.60709e-06, 1.18255e-06};
  TGraphErrors *gr_10bXiStar = new TGraphErrors(8, pt_points, LHC10b_XiStar, pt_points_e, LHC10b_XiStar_e);
  gr_10bXiStar->SetMarkerStyle(21);
  gr_10bXiStar->SetMarkerColor(2);
  gr_10bXiStar->SetLineColor(2);
  TGraphErrors *gr_10cXiStar = new TGraphErrors(8, pt_points, LHC10c_XiStar, pt_points_e, LHC10c_XiStar_e);
  gr_10cXiStar->SetMarkerStyle(22);
  gr_10cXiStar->SetMarkerColor(4);
  gr_10cXiStar->SetLineColor(4);
  TGraphErrors *gr_10dXiStar = new TGraphErrors(8, pt_points, LHC10d_XiStar, pt_points_e, LHC10d_XiStar_e);
  gr_10dXiStar->SetMarkerStyle(23);
  gr_10dXiStar->SetMarkerColor(6);
  gr_10dXiStar->SetLineColor(6);
  //gr_10bXiStar->Draw("P");
  //gr_10cXiStar->Draw("P");
  //gr_10dXiStar->Draw("P");
  //legend2->AddEntry(gr_10bXiStar,"LHC10b","p");
  //legend2->AddEntry(gr_10cXiStar,"LHC10c","p");
  //legend2->AddEntry(gr_10dXiStar,"LHC10d","p");

  // Old QM2011 values for (Xi* + cc)/2
  double Old_XiStarPoints[8]={0.001259, 0.0008586, 0.0006316, 0.0003604, 1.856e-4, 5.4601e-5, 2.0187e-5, 7.6559e-6};
  double Old_XiStarPoints_e[8]={2.44e-4, 1.224e-4, 7.136e-5, 3.85e-5, 1.59e-5, 7.49e-6, 3.51e-6, 1.822e-6};
  TGraphErrors *gr_OldXiStar = new TGraphErrors(8, pt_points, Old_XiStarPoints, pt_points_e, Old_XiStarPoints_e);
  gr_OldXiStar->SetMarkerStyle(21);
  gr_OldXiStar->SetMarkerColor(2);
  gr_OldXiStar->SetLineColor(2);
  //gr_OldXiStar->Draw("P");
  //legend2->AddEntry(gr_OldXiStar,"(#Xi(1530) + cc)/2: QM11 Results","p");
  /*
  // To fit the Old QM results
  TH1D *h_OldXiStar = (TH1D*)h_Xispectrum->Clone();
  h_OldXiStar->Add(h_OldXiStar,-1);
  for(int i=0; i<8; i++){ 
    h_OldXiStar->SetBinContent(i+2, Old_XiStarPoints[i]);
    h_OldXiStar->SetBinError(i+2, Old_XiStarPoints_e[i]);
  }
  h_OldXiStar->SetBinContent(1,+100);// first bin has no data
  h_OldXiStar->Fit(myLevy,"IME","",minFitpoint,maxFitpoint);
  cout<<"old mean pt = "<<myLevy->Moment(1,0.,10)<<endl;
  */
  //double pt_points_shifted[max_ptbins]={0};
  //TF1 *AvgX = new TF1("AvgX","x*myExp",0.6,5.6);
  
  //AvgX->SetParName(0,"A");
  //AvgX->SetParName(1,"T");
  //for(int i=0; i<8; i++){
  //cout<<AvgX->Integral(ptedges[i],ptedges[i+1])/myExp->Integral(ptedges[i],ptedges[i+1])<<endl;
  //}
   
  
 

                      
  double David_points_pt[18]={0.7, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.6, 1.8, 2.05, 2.4, 2.85, 3.5, 4.4, 5.45, 6.6, 7.85};
  
  
  // David's new points
  // David's Xi-
  //double David_points_y[18]={0.00470846, 0.00494042, 0.00474707, 0.0046308, 0.00426987, 0.00390345, 0.00345295, 0.0032779, 0.00269152, 0.00202204, 0.00145725, 0.00088381, 0.000460706, 0.000196437, 6.31383e-05, 1.69709e-05, 5.72763e-06, 2.22189e-06};
  //double David_points_y_e[18]={0.000313742, 0.000316252, 0.000294082, 0.000292834, 0.000255195, 0.000233124, 0.00020646, 0.000207611, 0.000160407, 0.000116738, 8.30453e-05, 5.03115e-05, 2.65867e-05, 1.14573e-05, 4.13348e-06, 1.4395e-06, 6.63351e-07, 4.01286e-07};
  // David's Xi+
  //double David_points_y[18]={0.00476135, 0.00480974, 0.00469503, 0.0046272, 0.0040248, 0.00364938, 0.00356068, 0.00306262, 0.0026221, 0.00196555, 0.00139214, 0.000856611, 0.000443817, 0.000192649, 6.06012e-05, 1.77073e-05, 5.04921e-06, 2.01156e-06};
  //double David_points_y_e[18]={0.000324969, 0.000323376, 0.000320961, 0.00029042, 0.000248456, 0.000219039, 0.000231559, 0.000187812, 0.000153252, 0.000114861, 7.95886e-05, 4.92256e-05, 2.57394e-05, 1.15276e-05, 4.04731e-06, 1.43029e-06, 6.10334e-07, 3.50592e-07};
  // Published ALICE points
  // Xi-
  double David_points_y[18]={4.734067e-03, 4.967289e-03, 4.772892e-03, 4.655989e-03, 4.293091e-03, 3.924680e-03, 3.471732e-03, 3.295728e-03, 2.706158e-03, 2.033039e-03, 1.465172e-03, 8.886171e-04, 4.632116e-04, 1.975056e-04, 6.348164e-05, 1.706320e-05, 5.758779e-06, 2.233978e-06};
  double David_points_y_e[18]={3.203305e-04, 3.230481e-04, 3.006813e-04, 2.992160e-04, 2.610787e-04, 2.384269e-04, 2.111132e-04, 2.118811e-04, 1.638975e-04, 1.193448e-04, 8.491450e-05, 5.149418e-05, 2.718932e-05, 1.174048e-05, 4.233300e-06, 1.475196e-06, 6.953602e-07, 4.110991e-07};
  // Xi+
  //double David_points_y[18]={4.787242e-03, 4.835899e-03, 4.720567e-03, 4.652362e-03, 4.046691e-03, 3.669227e-03, 3.580047e-03, 3.079273e-03, 2.636364e-03, 1.976242e-03, 1.399711e-03, 8.612696e-04, 4.462308e-04, 1.936968e-04, 6.093081e-05, 1.780359e-05, 5.076670e-06, 2.022504e-06};
  //double David_points_y_e[18]={3.315195e-04, 3.301755e-04, 3.272160e-04, 2.971086e-04, 2.538620e-04, 2.236666e-04, 2.363366e-04, 1.920731e-04, 1.566324e-04, 1.174641e-04, 8.139283e-05, 5.033151e-05, 2.632723e-05, 1.179350e-05, 4.145917e-06, 1.471487e-06, 6.243488e-07, 3.638050e-07};

  TGraphErrors *gr_DavidXi = new TGraphErrors(18, David_points_pt, David_points_y, 0, David_points_y_e); 
  gr_DavidXi->SetMarkerStyle(24);
  gr_DavidXi->SetMarkerSize(2);
  gr_DavidXi->SetMarkerColor(4);
  //gr_DavidXi->Draw("P");
  //legend2->AddEntry(gr_DavidXi,"ALICE Published #Xi^{+}","p");
  //legend2->AddEntry(gr_DavidXi,"#Xi^{+-}","p");

 

  //cout<<"Efficiencies:"<<endl;
  //for(int i=0; i<max_ptbins; i++) cout<<Eff[i]<<", ";
  //cout<<endl;
  //cout<<"Corresponding errors:"<<endl;
  //for(int i=0; i<max_ptbins; i++) cout<<Eff_e[i]<<", ";
  //cout<<endl;
  //cout<<"Raw Yields:"<<endl;
  //for(int i=0; i<max_ptbins; i++) cout<<yield[i]<<", ";
  //cout<<endl;
  //cout<<"Ratio of My spectrum to Davids:"<<endl;
  //for(int i=0; i<max_ptbins; i++) cout<<spectrum[i]/David_points_y[i]<<", ";
  //cout<<endl;
  //cout<<"Corresponding ratio errors:"<<endl;
  //for(int i=0; i<max_ptbins; i++) cout<<sqrt( pow(spectrum_e[i]/David_points_y[i],2) + pow(spectrum[i]*David_points_y_e[i]/(David_points_y[i]*David_points_y[i]),2))<<", ";
  //cout<<endl;
  //cout<<"Ratio of My Eff to Davids:"<<endl;
  //for(int i=0; i<max_ptbins; i++) cout<<Eff[i]/David_XiMinus_eff[i]<<", ";
  //cout<<endl;
  //cout<<"Corresponding ratio errors:"<<endl;
  //for(int i=0; i<max_ptbins; i++) cout<<sqrt( pow(Eff_e[i]/David_XiMinus_eff[i],2) + pow(Eff_e[i]*David_XiMinus_eff_e[i]/(David_XiMinus_eff[i]*David_XiMinus_eff[i]),2))<<", ";
  //cout<<endl;


  //Xi 10b+10c+10d with first eff bug fix, second bug fix is a very tiny change
  //double Eff_ref[18]={0.0155246, 0.0366893, 0.0522083, 0.0653027, 0.0795129, 0.0942177, 0.109902, 0.125146, 0.150523, 0.184772, 0.225384, 0.27233, 0.324995, 0.371011, 0.382669, 0.354245, 0.344473, 0.299919};
  // My eff with David's cuts XiMinus
  //double Eff_refXiMinus[18]={0.0105995, 0.0312385, 0.0477496, 0.0645244, 0.0812671, 0.0996255, 0.117602, 0.13177, 0.158334, 0.19623, 0.230584, 0.280132, 0.324593, 0.351765, 0.353174, 0.321661, 0.285389, 0.20558};
  //double Eff_refXiMinus_e[18]={0.000430506, 0.00115503, 0.00152918, 0.00189702, 0.00227815, 0.00268649, 0.0031259, 0.00357281, 0.00307359, 0.00400903, 0.00422553, 0.0052229, 0.00685152, 0.00835609, 0.0118039, 0.0168519, 0.0234869, 0.0274467};
  //double Eff_refXiPlus[18]={0.0103308, 0.0291685, 0.0473052, 0.0586698, 0.0741367, 0.0936288, 0.110488, 0.129917, 0.15185, 0.192956, 0.226617, 0.286435, 0.318134, 0.355376, 0.358586, 0.325878, 0.249803, 0.256949};
  //double Eff_refXiPlus_e[18]={0.000429664, 0.00112569, 0.00153199, 0.00180478, 0.00218061, 0.00261963, 0.00307867, 0.00358387, 0.00304106, 0.00398421, 0.00422447, 0.00531824, 0.00683356, 0.00842619, 0.012285, 0.0171912, 0.0224039, 0.0320204};
  
  // 10d1 eff
  double Eff_ref1[8]={0.00585138, 0.0271111, 0.0370528, 0.0767083, 0.0970953, 0.129362, 0.144039, 0.130414};
  double Eff_ref1_e[8]={0.00110308, 0.00301515, 0.00433263, 0.00845994, 0.00972208, 0.0175748, 0.0318935, 0.0367763};
  // 10d4 eff
  double Eff_ref2[8]={0.00646719, 0.0236728, 0.0454735, 0.0692855, 0.0922256, 0.13542, 0.161154, 0.126502};
  double Eff_ref2_e[8]={0.000740606, 0.00176465, 0.00315074, 0.00498191, 0.00600009, 0.0117801, 0.0203019, 0.0247307};
  // 10f6a eff
  double Eff_ref3[8]={0.00589799, 0.0228898, 0.043247, 0.0692314, 0.104212, 0.14611, 0.1587, 0.176539};
  double Eff_ref3_e[8]={0.000465501, 0.00114491, 0.00199253, 0.00333088, 0.00421532, 0.00822058, 0.0133099, 0.0202604};
  // 10f6a eff doubled binning
  //double Eff_ref3[16]={0.00288254, 0.00957083, 0.0177512, 0.0293962, 0.0382567, 0.0496768, 0.0650758, 0.074635, 0.096423, 0.116841, 0.142147, 0.15214, 0.169098, 0.143671, 0.185977, 0.162155};
  //double Eff_ref3_e[16]={0.000438375, 0.000885026, 0.00134525, 0.00195975, 0.00249222, 0.00323944, 0.00428832, 0.00525628, 0.00513653, 0.00726998, 0.0104245, 0.0133468, 0.0179437, 0.0196858, 0.0268477, 0.030696};

  // Pythia to Phojet eff ratio
  for(int ii=0; ii<max_ptbins; ii++){
    Eff_e[ii] = sqrt(pow(Eff_e[ii]/Eff_ref3[ii],2) + pow(Eff[ii]*Eff_ref3_e[ii]/pow(Eff_ref3[ii],2),2));
    Eff[ii] /= Eff_ref3[ii];
  }
  
  TGraphErrors *gr_Eff = new TGraphErrors(max_ptbins,pt_points, Eff, pt_points_e, Eff_e);
  gr_Eff->SetMarkerStyle(20);
  gr_Eff->SetMinimum(0);
  gr_Eff->SetMaximum(.4);
  //gr_Eff->GetXaxis()->SetLimits(0,6.1);
  gr_Eff->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  gr_Eff->GetYaxis()->SetTitle("Efficiency");
  //gr_Eff->SetTitle("#Xi(1530)^{0}+cc Efficiency");
  //gr_Eff->Draw("AP");
  //legend2->AddEntry(gr_Eff,"10d1+10d4+10f6a","p");
  //legend2->AddEntry(gr_Eff,"10f6 (Phojet)","p");
  //gr_Eff->Fit("pol0","IME","",1.6,4.8);

  TGraphErrors *gr_Eff_ref1 = new TGraphErrors(max_ptbins,pt_points, Eff_ref1, pt_points_e, Eff_ref1_e);
  gr_Eff_ref1->SetMarkerStyle(20);
  gr_Eff_ref1->SetMarkerColor(2);
  gr_Eff_ref1->SetLineColor(2);
  //gr_Eff_ref1->Draw("P");
  //legend2->AddEntry(gr_Eff_ref1,"10d1 MC","p");
  TGraphErrors *gr_Eff_ref2 = new TGraphErrors(max_ptbins,pt_points, Eff_ref2, pt_points_e, Eff_ref2_e);
  gr_Eff_ref2->SetMarkerStyle(20);
  gr_Eff_ref2->SetMarkerColor(4);
  gr_Eff_ref2->SetLineColor(4);
  //gr_Eff_ref2->Draw("P");
  //legend2->AddEntry(gr_Eff_ref2,"10d4 MC","p");
  TGraphErrors *gr_Eff_ref3 = new TGraphErrors(max_ptbins,pt_points, Eff_ref3, pt_points_e, Eff_ref3_e);
  gr_Eff_ref3->SetMarkerStyle(20);
  gr_Eff_ref3->SetMarkerColor(3);
  gr_Eff_ref3->SetLineColor(3);
  //gr_Eff_ref3->Draw("P");
  //legend2->AddEntry(gr_Eff_ref3,"10f6a (Pythia Perugia)","p");


  //TGraphErrors *gr_Eff_xiplus = new TGraphErrors(max_ptbins,pt_points, Eff_refXiPlus, 0, Eff_refXiPlus_e);
  //gr_Eff_xiplus->SetMarkerStyle(21);
  //gr_Eff_xiplus->SetMarkerColor(4);
  //gr_Eff_xiplus->GetXaxis()->SetTitle("pt (GeV/c)");
  //gr_Eff_xiplus->GetYaxis()->SetTitle("Efficiency");
  //gr_Eff_xiplus->SetTitle("#Xi(1530)^0 Efficiency");
  //gr_Eff_xiplus->Draw("P");


  
  TGraphErrors *gr_David_Eff_XiMinus = new TGraphErrors(18,pt_points_David,David_XiMinus_eff,0,David_XiMinus_eff_e);
  gr_David_Eff_XiMinus->SetMarkerStyle(24);
  gr_David_Eff_XiMinus->SetMarkerColor(2);
  //gr_David_Eff_XiMinus->Draw("P");
  //legend2->AddEntry(gr_David_Eff_XiMinus,"David/Antonin Xi- Eff","p");

  TGraphErrors *gr_David_Eff_XiPlus = new TGraphErrors(18,pt_points_David,David_XiPlus_eff,0,David_XiPlus_eff_e);
  gr_David_Eff_XiPlus->SetMarkerStyle(25);
  gr_David_Eff_XiPlus->SetMarkerColor(4);
  //gr_David_Eff_XiPlus->Draw("P");


  legend2->Draw("same");

  
  TLegend *legend3 = new TLegend(.35,.62,.9,.72,NULL,"brNDC");
  legend3->SetBorderSize(1);
  legend3->SetTextSize(.04);// small .03; large .036 
  //legend3->SetLineColor(0);
  legend3->SetFillColor(0);
  

  TGraphErrors *gr_RealWidth = new TGraphErrors(8, pt_points, width, pt_points_e, width_e);
  gr_RealWidth->SetMarkerStyle(20);
  gr_RealWidth->SetMinimum(0);
  gr_RealWidth->SetMaximum(0.0035);
  gr_RealWidth->GetYaxis()->SetTitle("#sigma (GeV/c)");
  gr_RealWidth->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  gr_RealWidth->SetTitle("");
  //gr_RealWidth->Draw("AP");
  //legend3->AddEntry(gr_RealWidth,"Real data (10b+10c+10d)","p");
  

  TGraphErrors *gr_MCWidth = new TGraphErrors(8, pt_points, MCwidth, pt_points_e, MCwidth_e);
  gr_MCWidth->SetMarkerStyle(20);
  gr_MCWidth->SetMarkerColor(2);
  gr_MCWidth->SetLineColor(2);
  gr_MCWidth->GetYaxis()->SetTitle("#sigma (GeV/c)");
  gr_MCWidth->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  gr_MCWidth->SetTitle("");
  //gr_MCWidth->Draw("P");
  //legend3->AddEntry(gr_MCWidth,"Pythia+ALICE (10d1+10d4+10f6a)","p");
  



  TGraphErrors *gr_RealMass = new TGraphErrors(8, pt_points, mass, pt_points_e, mass_e);
  gr_RealMass->SetMarkerStyle(20);
  gr_RealMass->SetMinimum(1.529);
  gr_RealMass->SetMaximum(1.5345);
  gr_RealMass->SetTitle("");
  gr_RealMass->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  gr_RealMass->GetYaxis()->SetTitleOffset(1.3);
  gr_RealMass->GetYaxis()->SetTitle("Mass (GeV/c^{2})");
  //gr_RealMass->Draw("AP");
  //legend3->AddEntry(gr_RealMass,"Real data (10b+10c+10d)","p");

  TGraphErrors *gr_MCMass = new TGraphErrors(8, pt_points, MCmass, pt_points_e, MCmass_e);
  gr_MCMass->SetMarkerStyle(20);
  gr_MCMass->SetMarkerColor(2);
  gr_MCMass->SetTitle("");
  gr_MCMass->GetYaxis()->SetTitle("Mass (GeV/c^{2})");
  gr_MCMass->GetYaxis()->SetTitleOffset(1.3);
  gr_MCMass->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  //gr_MCMass->Draw("P");
  //legend3->AddEntry(gr_MCMass,"Pythia+ALICE (10d1+10d4+10f6a)","p");
  //TF1 *pdg_massline=new TF1("pdg_massline","1.53178",0,10);
  //pdg_massline->Draw("same");
  
  //legend3->Draw("same");
  TF1 *myLevy2=new TF1("myLevy2",myLevyPt,0,8.5,3);
  myLevy2->SetLineColor(4);
  myLevy2->SetParName(0,"dN/dy");
  myLevy2->SetParName(1,"C");
  myLevy2->SetParName(2,"n");
  myLevy2->SetParameter(0,.003);
  myLevy2->SetParameter(1,.3);
  myLevy2->SetParameter(2,15);
  myLevy2->SetParLimits(0,0.0001,0.01);
  myLevy2->SetParLimits(1,.1,1);
  myLevy2->SetParLimits(2,1,500);
  TF1 *myExp2 = new TF1("myExp2",myExpFit,0,8.5,2);
  myExp2->SetLineColor(4);
  myExp2->SetParName(0,"dN/dy");
  myExp2->SetParName(1,"C");
  myExp2->SetParameter(0,.003);
  myExp2->SetParameter(1,.3);
  myExp2->SetParLimits(0,0.0001,0.1);
  myExp2->SetParLimits(1,.1,1);
  

  MCinput_Spectrum->SetMarkerStyle(20);
  MCinput_Spectrum->SetMarkerColor(2);
  MCinput_Spectrum->SetLineColor(2);
  MCinput_Spectrum->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  MCinput_Spectrum->GetYaxis()->SetTitle("1/N_{E}d^{2}N/dydp_{T} (GeV/c)^{-1}");
  //MCinput_Spectrum->SetTitle("Pythia Perugia");
  //MCinput_Spectrum->SetTitle("Phojet");

  MCinput_Spectrum->Sumw2();
  MCinput_Spectrum->Scale(1/MCinput_Spectrum->GetXaxis()->GetBinWidth(1));
  MCinput_Spectrum->Scale(1/EventsMC_postPV->GetEntries());
 
  // Scale Pythia to Data fit Tsallis
  MCinput_Spectrum->Scale(myLevy->Integral(0.8,1.2)/(MCinput_Spectrum->Integral(9,12)*MCinput_Spectrum->GetXaxis()->GetBinWidth(1)));
  //
  cout<<"low-pt extrapolation from Levy = "<<myLevy->Integral(0,0.8)<<"    Pythia Histogram Integral = "<<MCinput_Spectrum->Integral(1,8)*MCinput_Spectrum->GetXaxis()->GetBinWidth(1)<<endl;
  cout<<"high pt extrapolation over total Levy = "<<myLevy->Integral(5.6,10.)/myLevy->Integral(0,10.)<<endl;
  cout<<"high pt Pythia over total Pythia = "<<MCinput_Spectrum->Integral(56,100)/MCinput_Spectrum->Integral(1,100)<<endl;
  cout<<"Total yield Pythia over total Levy = "<<MCinput_Spectrum->Integral(1,100)*MCinput_Spectrum->GetXaxis()->GetBinWidth(1)/myLevy->Integral(0,10.)<<endl;
  
  //MCinput_Spectrum->Draw();
  //minFitpoint=0.; maxFitpoint=10.;
  //MCinput_Spectrum->Fit(myLevy2,"IME","",minFitpoint,maxFitpoint);
  //MCinput_Spectrum->Fit(myExp,"IME","",minFitpoint,maxFitpoint);
  //myExp->Draw("same");
  //myLevy2->Draw("same");
  //myLevy->Draw("same");
  //legend3->AddEntry(MCinput_Spectrum,"Pythia Perugia (Scaled to Real Data)","p");
  //legend3->AddEntry(myLevy2,"Phojet Levy fit (fit range from pt=0.8 to 5.6 GeV/c","l");
  //legend3->AddEntry(myLevy,"Real data Levy fit","l");
  
  double xaxis[60]={0};
  double yaxis[60]={0};
  for(int ii=0; ii<60; ii++){
    xaxis[ii] = MCinput_Spectrum->GetXaxis()->GetBinCenter(ii+1);
    yaxis[ii] = myLevy->Eval(xaxis[ii])/myLevy2->Eval(xaxis[ii]);
  }
  TGraph *LevyRatio = new TGraph(60, xaxis, yaxis);
  LevyRatio->SetTitle("Levy Ratio");
  LevyRatio->GetXaxis()->SetTitle("p_{T}");
  LevyRatio->GetYaxis()->SetTitle("(Real Data Levy)/(Pythia fit Levy)");
  LevyRatio->SetMarkerStyle(20);
  //LevyRatio->Draw("AP");
  

  //legend3->Draw("same");

  //cout<<"low-pt extrapolation from Levy = "<<myLevy2->Integral(0,0.8)<<"   Histogram Integral = "<<MCinput_Spectrum->Integral(1,8)*MCinput_Spectrum->GetXaxis()->GetBinWidth(1)<<endl;
  //cout<<"Total Bin Count yield = "<<MCinput_Spectrum->Integral(1,MCinput_Spectrum->GetNbinsX())*MCinput_Spectrum->GetXaxis()->GetBinWidth(1)<<endl;
  //cout<<"low-pt extrapolation from Exp = "<<myExp->Integral(0,0.8)<<"   Histogram Integral = "<<MCinput_Spectrum->Integral(1,8)*MCinput_Spectrum->GetXaxis()->GetBinWidth(1)<<endl;


  if(SaveToFile){  
    TFile *outputfile = new TFile(outfilename->Data(),"RECREATE");
    can->Write("can");
    can2->Write("can2");
    myLevy->Write("myLevy");
    h_Xispectrum->Write("h_Spectrum");
    RawYields->Write("RawYields");
    Efficiency->Write("Efficiency");

    outputfile->Close();
  }
  
  cout<<endl;
  
}

//________________________________________________________________________
double PolFunction(double *x, double *par){
 
  if (reject && x[0] > yieldRange[0] && x[0] < yieldRange[1]) {
    TF1::RejectPoint();
    return 0;
  }
  
  if(POLdegree==1) return par[0] + par[1]*x[0];
  else if(POLdegree==2) return par[0] + par[1]*x[0] + par[2]*pow(x[0],2);
  else if(POLdegree==3) return par[0] + par[1]*x[0] + par[2]*pow(x[0],2) + par[3]*pow(x[0],3);
  else if(POLdegree==4) return par[0] + par[1]*x[0] + par[2]*pow(x[0],2) + par[3]*pow(x[0],3) + par[4]*pow(x[0],4);
  else if(POLdegree==5) return par[0] + par[1]*x[0] + par[2]*pow(x[0],2) + par[3]*pow(x[0],3) + par[4]*pow(x[0],4) + par[5]*pow(x[0],5);
  else return 0;
  
}
//________________________________________________________________________
double PolFunctionSpecialXi(double *x, double *par){
 
  if (x[0] > 1.28 && x[0] < 1.3) {
    TF1::RejectPoint();
    return 0;
  }
  if (reject && x[0] > 1.31 && x[0] < 1.334) {
    TF1::RejectPoint();
    return 0;
  }
  
  if(POLdegree==1) return par[0] + par[1]*x[0];
  else if(POLdegree==2) return par[0] + par[1]*x[0] + par[2]*pow(x[0],2);
  else if(POLdegree==3) return par[0] + par[1]*x[0] + par[2]*pow(x[0],2) + par[3]*pow(x[0],3);
  else if(POLdegree==4) return par[0] + par[1]*x[0] + par[2]*pow(x[0],2) + par[3]*pow(x[0],3) + par[4]*pow(x[0],4);
  else if(POLdegree==5) return par[0] + par[1]*x[0] + par[2]*pow(x[0],2) + par[3]*pow(x[0],3) + par[4]*pow(x[0],4) + par[5]*pow(x[0],5);
  else return 0;
  
}
//________________________________________________________________________
double BWFunction(double *x, double *par){
  return (par[0]*par[2]/(1000.*2.*3.1415926))/( pow(x[0]-par[1],2) + par[2]*par[2]/4.);
  //return (par[0]*par[1]/(1000.*2.*3.1415926))/( pow(x[0]-par[2],2) + par[1]*par[1]/4.);
}
//________________________________________________________________________
double BWplusPol(double *x, double *par){
  return BWFunction(x,par) + PolFunction(x,&par[3]);
}
//________________________________________________________________________
double GausFunction(double *x, double *par){
  return (par[0]*.001)/(sqrt(2*3.14159*par[2]*par[2]))*exp(-pow((x[0]-par[1])/(sqrt(2)*par[2]),2));
}
//________________________________________________________________________
double GausplusPol(double *x, double *par){
  return GausFunction(x,par) + PolFunction(x,&par[3]);
}
//________________________________________________________________________
double GausplusPolSpecialXi(double *x, double *par){
  return GausFunction(x,par) + PolFunctionSpecialXi(x,&par[3]);
}
//________________________________________________________________________
double myLevyPt(Double_t *x, Double_t *par)
{
  Double_t lMass=0;
  if(XiStarCase) lMass = 1.5318; //Xi mass
  else lMass = 1.32171; //Xi mass

  Double_t ldNdy  = par[0]; // dN/dy
  Double_t l2pi   = 2*TMath::Pi(); // 2pi
  Double_t lTemp = par[1]; // Temperature
  Double_t lPower = par[2]; // power=n
  
  Double_t lBigCoef = ((lPower-1)*(lPower-2)) / (l2pi*lPower*lTemp*(lPower*lTemp+lMass*(lPower-2)));
  Double_t lInPower = 1 + (TMath::Sqrt(x[0]*x[0]+lMass*lMass)-lMass) / (lPower*lTemp);
  
  return l2pi * ldNdy * x[0] * lBigCoef * TMath::Power(lInPower,(-1)*lPower);
}
//________________________________________________________________________
double myExpFit(Double_t *x, Double_t *par)
{
  Double_t lMass=0;
  if(XiStarCase) lMass = 1.5318; //Xi mass
  else lMass = 1.32171; //Xi mass

  return 2*TMath::Pi()*x[0]*par[0]*exp(-TMath::Sqrt(x[0]*x[0]+lMass*lMass)/par[1]);
}
//________________________________________________________________________
double Voigt(Double_t *x, Double_t *par)
{// code taken directly from RooVoigtian in ROOT
  
  Double_t Norm = par[0];
  Double_t mean = par[1];
  Double_t s = fabs(par[2]);// sigma; Gaussian sigma
  Double_t w = fabs(par[3]);// Width; BW width

  Double_t arg = x[0] - mean;
  Double_t _invRootPi = 1./sqrt(atan2(0.,-1.));

  // return constant for zero width and sigma
  if (s==0. && w==0.) return 1.;

  // Breit-Wigner for zero sigma
  if (s==0.) return (Norm*1./(arg*arg+0.25*w*w));
  Double_t coef= -0.5/(s*s);
  
  // Gauss for zero width
  if (w==0.) return Norm*exp(coef*arg*arg);

  // actual Voigtian for non-trivial width and sigma
  Double_t c = 1./(sqrt(2.)*s);
  Double_t a = 0.5*c*w;
  Double_t u = c*arg;
  RooComplex z(u,a);
  RooComplex v(0.);


  v = RooMath::ComplexErrFunc(z);
  
  return Norm*c*_invRootPi*v.re();

}
//________________________________________________________________________
 double VoigtplusPol(double *x, double *par){
   return Voigt(x,par) + PolFunction(x,&par[4]);
}
