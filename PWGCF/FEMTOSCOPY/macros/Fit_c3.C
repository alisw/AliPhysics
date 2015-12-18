#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <Riostream.h>
#include <assert.h>

#include "TVector2.h"
#include "TFile.h"
#include "TString.h"
#include "TF1.h"
#include "TF3.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TMath.h"
#include "TText.h"
#include "TRandom3.h"
#include "TArray.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMinuit.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TPaveStats.h"
#include "TASImage.h"
#include "TGraph.h"
#include "TSpline.h"
#include "TVirtualFitter.h"

#define BohrR 1963.6885 // Bohr Radius for pions
#define FmToGeV 0.19733 // conversion to fm
#define PI 3.1415926
#define masspiC 0.1395702 // pi+ mass (GeV/c^2)

using namespace std;

// Fit Normalization is fixed to 1.0, explains why Chi2 higher here than in Three-pion radii paper
bool PrintData=0;
//
int CollisionType_def=0;// PbPb, pPb, pp
int FitType=0;// (0)Edgeworth, (1)Laguerre, (2)Levy
bool FixNorm=1;// Fix Normalization (PbPb only)
//
int Mbin=0;// 0-9: centrality bin in widths of 5%
int CHARGE=-1;// -1 or +1: + or - pions for same-charge case, --+ or -++,  ---+ or -+++
//
int EDbin=0;// 0 or 1: Kt3 bin
double G_def = 0;// coherent fraction (percentage)
double Rcoh_def = 100.;// Radius of Gaussian coherent component in fm, 100 for full sized source
//
bool MRC=1;// Momentum Resolution Corrections?
bool MuonCorrection=1;// correct for Muons?
//
float TwoFrac=0.7;// Lambda parameter
int f_choice=0;// 0(Core/Halo), 1(60fm), 2(80fm), 3(100fm)
//
//
//
//
bool CumulantFit=kFALSE;// (0) C3, (1) c3
bool SaveToFile_def=0;
int fFSIindex=0;
float OneFrac;// Lambda^{1/2}
float ThreeFrac;// Lambda^{3/2}
double Qcut_pp = 0.6;// 0.6
double Qcut_PbPb = 0.1;
double NormQcutLow_pp = 1.0;
double NormQcutHigh_pp = 1.5;
double NormQcutLow_PbPb = .15;
double NormQcutHigh_PbPb = .2;// was .175

int TextFont=42;// 63, or 42
float SizeLabel=0.06;// 20(63 font), 0.08(42 font)
float SizeLegend=0.06;// .08
float SizeTitle=0.06;// 
float SizeSpecif=0.045;// 
float SF1=2/3.*0.95;
float SF2=1/2.*0.95;
float MarkerSize=1.25;

double RightMargin=0.004;// 0.002

const int BINRANGE_3=60;// q12,q13,q23
int BINLIMIT_3;
double A_3[BINRANGE_3][BINRANGE_3][BINRANGE_3];// q12,q13,q23
double A_3_e[BINRANGE_3][BINRANGE_3][BINRANGE_3];// q12,q13,q23
double B_3[BINRANGE_3][BINRANGE_3][BINRANGE_3];// q12,q13,q23
double a_3[BINRANGE_3][BINRANGE_3][BINRANGE_3];// q12,q13,q23
double a_3_e[BINRANGE_3][BINRANGE_3][BINRANGE_3];// q12,q13,q23
double b_3[BINRANGE_3][BINRANGE_3][BINRANGE_3];// q12,q13,q23
double BinCenters[400];
double Chi2_c3;
double NFitPoints_c3;
double Q3LimitLow;
double Q3LimitHigh;
void fcn_c3(int&, double*, double&, double[], int);

int RbinMRC;

TH1D *fFSIss[12];
TH1D *fFSIos[12];

void SetFSICorrelations();
void SetFSIindex(Float_t);
Float_t FSICorrelation(Int_t, Int_t, Float_t);
void SetMuonCorrections();
void SetMomResCorrections();
double Gamov(int, double);

//
float fact(float);



TH1D *MRC_SC_3[3];
TH1D *C3muonCorrectionSC[2];

double AvgQinvSS[30];
double AvgQinvOS[30];
double BinCentersQ4[20];



void Fit_c3(bool SaveToFile=SaveToFile_def, int CollisionType=CollisionType_def, double G=G_def, double Rcoh=Rcoh_def){

  gStyle->SetErrorX(0.0001);
  
  SaveToFile_def=SaveToFile;
  CollisionType_def=CollisionType;
  G /= 100.;
  G_def=G;
  Rcoh_def=Rcoh;

  if(CollisionType!=0) FixNorm=1;

  
  TFile *_file0;
  if(CollisionType==0){// PbPb
    //_file0 = new TFile("Results/PDC_11h_3Dhistos.root","READ");
    //_file0 = new TFile("Results/PDC_11h_RcohFinite.root","READ");
    _file0 = new TFile("Results/PDC_11h_LowNorm_HighNorm.root","READ");
  }else if(CollisionType==1){// pPb
    //_file0 = new TFile("Results/PDC_13bc_kINT7_LowNorm.root","READ");
    _file0 = new TFile("Results/PDC_13bc_kINT7_LowNorm_HighNorm.root","READ");
  }else{// pp
    //_file0 = new TFile("Results/PDC_10bcde_kMB_3Dhisto_LowNorm_HighNorm.root","READ");// used for paper proposal
    _file0 = new TFile("Results/PDC_10bcde_kMB_4kappas_Norm0p9to1p2.root","READ");
  }

  TList *MyList;
  TDirectoryFile *tdir = (TDirectoryFile*)_file0->Get("PWGCF.outputFourPionAnalysis.root");
  if(CollisionType==0 || CollisionType==2) MyList=(TList*)tdir->Get("FourPionOutput_1");
  else MyList=(TList*)tdir->Get("FourPionOutput_2");
  //MyList=(TList*)_file0->Get("MyList");

  _file0->Close();


  if(CollisionType==0) {Q3LimitLow = 0.01; Q3LimitHigh = 0.08;}// 0.01 and 0.08 
  else {Q3LimitLow = 0.01; Q3LimitHigh = 0.25;}// 0.01 and 0.25
  
  // Pb-Pb points from main analysis
  double C3_pointsPbPb_accepted[10]={0, 3.80885, 2.82493, 2.10035, 1.66763, 1.42446, 1.28301, 1.19701, 1.14319, 1.10818};
  double c3_pointsPbPb_accepted[10]={0, 1.8686, 1.45941, 1.21784, 1.08782, 1.03138, 1.00964, 1.00001, 0.996496, 0.99515};
  

  //
  OneFrac = sqrt(TwoFrac);
  ThreeFrac = pow(TwoFrac, 1.5);
  
  //// Core/Halo, 60fm, 80fm, 100fm
  float ThermShift_f33[4]={0., 0.05884, 0.04487, 0.04132};// % shift from core/halo
  float ThermShift_f32[4]={0., -0.01622, -0.01259, -0.01169};
  float ThermShift_f31[4]={0., -0.01019, -0.007113, -0.006259};
  
  float f33prime = ThreeFrac;
  float f32prime = TwoFrac*(1-OneFrac);
  float f31prime = pow(1-OneFrac,3) + 3*OneFrac*pow(1-OneFrac,2);
  f33prime += ThermShift_f33[f_choice];
  f32prime += ThermShift_f32[f_choice];
  f31prime += ThermShift_f31[f_choice];
  float f33 = f33prime;
  float f32 = f32prime/TwoFrac;
  float f31 = f31prime - 3*((1-TwoFrac)*(1-OneFrac) + ThermShift_f32[f_choice]*(1-TwoFrac)/TwoFrac);
  //cout<<f33 + 3*f32 + f31<<endl;

 
  cout<<"Mbin = "<<Mbin<<"   KT3 = "<<EDbin<<endl;
  
  if(CollisionType==0){
    if(Mbin==0) {RbinMRC=10;}
    else if(Mbin==1) {RbinMRC=9;}
    else if(Mbin<=3) {RbinMRC=8;}
    else if(Mbin<=5) {RbinMRC=7;}
    else {RbinMRC=6;}
  }else{
    RbinMRC=2;
  }


  if(CollisionType==0) BINLIMIT_3=20;
  else BINLIMIT_3=30;
 

  TString *System=new TString("");
  if(CollisionType==0) System->Append("ALICE 0-5% Pb-Pb #\sqrt{#font[12]{s}_{NN}}=2.76 TeV");
  else if(CollisionType==1) System->Append("ALICE p-Pb #\sqrt{#font[12]{s}_{NN}}=5.02 TeV");
  else System->Append("ALICE pp #\sqrt{#font[12]{s}}=7 TeV");
  double XstartALICE=0.43;
  if(CollisionType==1) XstartALICE=0.54;
  if(CollisionType==2) XstartALICE=0.65;
  TLatex *ALICEspecif = new TLatex(XstartALICE,.93,System->Data());// ALICE specifications
  ALICEspecif->SetNDC(1);
  ALICEspecif->SetTextFont(TextFont);
  ALICEspecif->SetTextSize(SizeSpecif);
  //
  TString *KT=new TString("");
  if(EDbin==0) KT->Append("0.16<#font[12]{K}_{T3}<0.3 GeV/#font[12]{c}");
  else KT->Append("0.3<#font[12]{K}_{T3}<1.0 GeV/#font[12]{c}");
  TLatex *KTspecif = new TLatex(0.64,0.5,KT->Data());
  KTspecif->SetNDC(1);
  KTspecif->SetTextFont(TextFont);
  KTspecif->SetTextSize(SizeSpecif);
    
  TString *Centname=new TString("0-5%");
  TLatex *Centspecif = new TLatex(0.52,0.8,Centname->Data());
  Centspecif->SetNDC(1);
  Centspecif->SetTextFont(TextFont);
  Centspecif->SetTextSize(SizeSpecif);

  TString *ChargeCombSt3 = new TString("");
  ChargeCombSt3->Append("#pi^{-}#pi^{-}#pi^{-}");
  float xCCT=0.77, yCCT=0.4;
  TLatex *ChargeCombTitle3 = new TLatex(xCCT,yCCT,ChargeCombSt3->Data());// 0.7,yCCT
  ChargeCombTitle3->SetNDC(1);
  ChargeCombTitle3->SetTextFont(TextFont);
  ChargeCombTitle3->SetTextSize(2*SizeSpecif);

  // bin centers from QS+FSI
  double BinCenterPbPbCentral[40]={0.00206385, 0.00818515, 0.0129022, 0.0177584, 0.0226881, 0.027647, 0.032622, 0.0376015, 0.042588, 0.0475767, 0.0525692, 0.0575625, 0.0625569, 0.0675511, 0.0725471, 0.0775436, 0.0825399, 0.0875364, 0.0925339, 0.0975321, 0.102529, 0.107527, 0.112525, 0.117523, 0.122522, 0.12752, 0.132519, 0.137518, 0.142516, 0.147515, 0.152514, 0.157513, 0.162513, 0.167512, 0.172511, 0.177511, 0.18251, 0.187509, 0.192509, 0.197509};
  double BinCenterpPbAndpp[40]={0.00359275, 0.016357, 0.0257109, 0.035445, 0.045297, 0.0552251, 0.0651888, 0.0751397, 0.0851088, 0.0951108, 0.105084, 0.115079, 0.12507, 0.135059, 0.145053, 0.155049, 0.16505, 0.175038, 0.185039, 0.195034, 0.205023, 0.215027, 0.225024, 0.235023, 0.245011, 0.255017, 0.265017, 0.275021, 0.285021, 0.295017, 0.305018, 0.315018, 0.325013, 0.335011, 0.345016, 0.355019, 0.365012, 0.375016, 0.385017, 0.395016};
  if(CollisionType==0){
    for(int i=0; i<40; i++) BinCenters[i] = BinCenterPbPbCentral[i];
  }else{
    for(int i=0; i<40; i++) BinCenters[i] = BinCenterpPbAndpp[i];
  }
  
  // extend BinCenters for high q
  for(int index=40; index<400; index++){
    if(CollisionType==0) BinCenters[index] = (index+0.5)*(0.005);
    else BinCenters[index] = (index+0.5)*(0.010);
  }
  // Set 0's to 3-particle fit arrays
  for(int i=1; i<=BINLIMIT_3; i++){// bin number
    for(int j=1; j<=BINLIMIT_3; j++){// bin number
      for(int k=1; k<=BINLIMIT_3; k++){// bin number
	A_3[i-1][j-1][k-1]=0;
	A_3_e[i-1][j-1][k-1]=0;
	B_3[i-1][j-1][k-1]=0;
      }
    }
  }


  //
  TH3D *ThreeParticle[2][2][2][5];// ch1,ch2,ch3,term
  TProfile3D *K3avg[2][2][2][4];
  double norm_3[5]={0};
  //

  gStyle->SetOptStat(0);
  gStyle->SetOptDate(0);
  gStyle->SetOptFit(0111);


  
  SetFSIindex(10.);
  SetFSICorrelations();
  SetMomResCorrections();
  SetMuonCorrections();
  //
  /////////////////////////////////////////////////////////
  

 
  TH1D *Events = (TH1D*)MyList->FindObject("fEvents2");
  //

  cout<<"#Events = "<<Events->Integral(Mbin+1,Mbin+1)<<endl;

  
  
  ///////////////////////////////////
  // Get Histograms
  
  for(int term=0; term<5; term++){
    
    TString *name3 = new TString("ThreeParticle_Charge1_");
    *name3 += 0;
    name3->Append("_Charge2_");
    *name3 += 0;
    name3->Append("_Charge3_");
    *name3 += 0;
    name3->Append("_M_");
    *name3 += Mbin;
    name3->Append("_ED_");
    *name3 += EDbin;
    name3->Append("_Term_");
    *name3 += term+1;
    
    TString *nameNorm3=new TString(name3->Data());
    nameNorm3->Append("_Norm");
    //
    TString *nameK3=new TString(name3->Data());
    nameK3->Append("_Kfactor3D");
    //
    name3->Append("_3D");
    
    
    
    norm_3[term] = ((TH1D*)MyList->FindObject(nameNorm3->Data()))->Integral();
    ThreeParticle[0][0][0][term] = (TH3D*)MyList->FindObject(name3->Data());
    ThreeParticle[0][0][0][term]->Sumw2();
    //if(0==0 && 0==0) cout<<"3-pion norms  "<<norm_3[term]<<endl;
    ThreeParticle[0][0][0][term]->Scale(norm_3[0]/norm_3[term]);
    ThreeParticle[0][0][0][term]->SetMarkerStyle(20);
    ThreeParticle[0][0][0][term]->SetTitle("");
    //
    
    //
    if(term<4){
      K3avg[0][0][0][term] = (TProfile3D*)MyList->FindObject(nameK3->Data());
    }
    
    //
  }// term 
  
  
  

  cout<<"Done getting Histograms"<<endl;
  
  TF1 *Unity = new TF1("Unity","1",0,100);
  Unity->SetLineStyle(2);
  Unity->SetLineColor(1);

  int ch1=0,ch2=0,ch3=0;
  
  
  
  ///////////////////////////////////////////////////////////
  // 3-pion 
  cout<<"3-pion section"<<endl;
 
  TCanvas *can1 = new TCanvas("can1", "can1",10,0,700,800);// 11,53,700,500
  can1->SetHighLightColor(2);
  gStyle->SetOptFit(0111);
  can1->SetFillColor(0);//10
  can1->SetBorderMode(0);
  can1->SetBorderSize(2);
  can1->SetFrameFillColor(0);
  can1->SetFrameBorderMode(0);
  can1->SetFrameBorderMode(0);
    
  TPad *pad1 = new TPad("pad1","pad1",0.,0.3,1.,1.);
  TPad *pad1_2 = new TPad("pad1_2","pad1_2",0.,0.,1.,0.3);
  pad1->Draw();
  pad1_2->Draw();
  pad1->cd();
  gPad->SetLeftMargin(0.11); gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.03); 
  gPad->SetBottomMargin(0.01);
    

  TLegend *legend1 = new TLegend(.65,.55, .85,.9,NULL,"brNDC");//.45 or .4 for x1
  legend1->SetBorderSize(0);
  legend1->SetFillColor(0);
  legend1->SetTextFont(42);
  legend1->SetTextSize(SizeLegend);
  TLegend *legend2 = (TLegend*)legend1->Clone();
  legend2->SetY1(0.6); legend2->SetY2(0.95);
  legend2->SetTextSize(2.5*SizeLegend);
  
  int Q3BINS=12;
  float Q3HistoLimit=0.12;
  if(CollisionType!=0){ Q3BINS=60; Q3HistoLimit=0.6;}

  TH1D *C3_hist = new TH1D("C3_hist","",Q3BINS,0,Q3HistoLimit);
  TH1D *c3_hist = new TH1D("c3_hist","",Q3BINS,0,Q3HistoLimit);
  TH1D *Combinatorics_1d = new TH1D("Combinatorics_1d","",Q3BINS,0,Q3HistoLimit);
  TH1D *GenSignalExpected_num[2];
  GenSignalExpected_num[0]=new TH1D("GenSignalExpected_num1","",Q3BINS,0,Q3HistoLimit);
  GenSignalExpected_num[1]=new TH1D("GenSignalExpected_num2","",Q3BINS,0,Q3HistoLimit);
  TH1D *GenSignalExpected_den[2];
  GenSignalExpected_den[0]=new TH1D("GenSignalExpected_den1","",Q3BINS,0,Q3HistoLimit);
  GenSignalExpected_den[1]=new TH1D("GenSignalExpected_den2","",Q3BINS,0,Q3HistoLimit);
  double c3_e[100]={0};
  double C3_e[100]={0};
  
  double value_num; 
  double value_num_e;
  double N3_QS;
  double N3_QS_e;
  
  for(int i=2; i<=ThreeParticle[0][0][0][0]->GetNbinsX(); i++){// bin number
    double Q12 = BinCenters[i-1];// true center
    //int Q12bin = int(Q12/0.01)+1;
    //
    for(int j=2; j<=ThreeParticle[0][0][0][0]->GetNbinsY(); j++){// bin number
      double Q13 = BinCenters[j-1];// true center
      //int Q13bin = int(Q13/0.01)+1;
      //
      for(int k=2; k<=ThreeParticle[0][0][0][0]->GetNbinsZ(); k++){// bin number
	double Q23 = BinCenters[k-1];// true center
	//int Q23bin = int(Q23/0.01)+1;
	//
		
	if(Q12 < sqrt(pow(Q13,2)+pow(Q23,2) - 2*Q13*Q23)) continue;// not all configurations are possible
	if(Q12 > sqrt(pow(Q13,2)+pow(Q23,2) + 2*Q13*Q23)) continue;// not all configurations are possible
	
	double Q3 = sqrt(pow(Q12,2) + pow(Q13,2) + pow(Q23,2));
       	int Q3bin = c3_hist->GetXaxis()->FindBin(Q3);
	
       	//
	double K3 = K3avg[0][0][0][0]->GetBinContent(i,j,k);
	double K2_12 = K3avg[0][0][0][1]->GetBinContent(i,j,k);
	double K2_13 = K3avg[0][0][0][2]->GetBinContent(i,j,k);
	double K2_23 = K3avg[0][0][0][3]->GetBinContent(i,j,k);

	
	if(K3==0) continue;

	double TERM1=ThreeParticle[ch1][ch2][ch3][0]->GetBinContent(i,j,k);// N3 (3-pion yield per q12,q13,q23 cell). 3-pions from same-event
	double TERM2=ThreeParticle[ch1][ch2][ch3][1]->GetBinContent(i,j,k);// N2*N1. 1 and 2 from same-event
	double TERM3=ThreeParticle[ch1][ch2][ch3][2]->GetBinContent(i,j,k);// N2*N1. 1 and 3 from same-event
	double TERM4=ThreeParticle[ch1][ch2][ch3][3]->GetBinContent(i,j,k);// N2*N1. 2 and 3 from same-event
	double TERM5=ThreeParticle[ch1][ch2][ch3][4]->GetBinContent(i,j,k);// N1*N1*N1. All from different events (pure combinatorics)
	
	
	if(TERM1==0 && TERM2==0 && TERM3==0 && TERM4==0 && TERM5==0) continue;
	if(TERM1==0 || TERM2==0 || TERM3==0 || TERM4==0 || TERM5==0) continue;
	
     	//
	if(MRC){
	  TERM1 *= MRC_SC_3[0]->GetBinContent(MRC_SC_3[0]->GetXaxis()->FindBin(Q3));
	  TERM2 *= MRC_SC_3[1]->GetBinContent(MRC_SC_3[1]->GetXaxis()->FindBin(Q3));
	  TERM3 *= MRC_SC_3[1]->GetBinContent(MRC_SC_3[1]->GetXaxis()->FindBin(Q3));
	  TERM4 *= MRC_SC_3[1]->GetBinContent(MRC_SC_3[1]->GetXaxis()->FindBin(Q3));
	  TERM5 *= MRC_SC_3[2]->GetBinContent(MRC_SC_3[2]->GetXaxis()->FindBin(Q3));
	}
	double MuonCorr1=1, MuonCorr2=1, MuonCorr3=1, MuonCorr4=1;
	if(MuonCorrection){
	  MuonCorr1 = C3muonCorrectionSC[0]->GetBinContent(C3muonCorrectionSC[0]->GetXaxis()->FindBin(Q3));
	  MuonCorr2 = C3muonCorrectionSC[1]->GetBinContent(C3muonCorrectionSC[0]->GetXaxis()->FindBin(Q3));
	  MuonCorr3 = MuonCorr2;
	  MuonCorr4 = MuonCorr2;
	}

		
	
	// Purify.  Isolate pure 3-pion QS correlations using Lambda and K3 (removes lower order correlations)
	N3_QS = TERM1;
	N3_QS -= f31 * TERM5;
	N3_QS -= f32 * (TERM2 + TERM3 + TERM4);
	N3_QS /= f33;
	N3_QS *= K3;
	N3_QS *=  MuonCorr1;

	
	// Isolate 3-pion cumulant
	value_num = N3_QS;
	value_num -= (TERM2 - (1-TwoFrac)*TERM5)*K2_12/TwoFrac * MuonCorr2;
	value_num -= (TERM3 - (1-TwoFrac)*TERM5)*K2_13/TwoFrac * MuonCorr3;
	value_num -= (TERM4 - (1-TwoFrac)*TERM5)*K2_23/TwoFrac * MuonCorr4;
	value_num += 2*TERM5;
	
	
	
	
	// errors
	N3_QS_e = TERM1;
	N3_QS_e += pow(( pow(1-OneFrac,3) +3*OneFrac*pow(1-OneFrac,2) )*sqrt(TERM5),2);
	N3_QS_e += (pow((1-OneFrac),2)*(TERM2 + TERM3 + TERM4) + pow((1-OneFrac)*3*(1-TwoFrac)*sqrt(TERM5),2));
	N3_QS_e /= pow(ThreeFrac,2);
	N3_QS_e *= pow(K3,2);
	//
	value_num_e = N3_QS_e;
	value_num_e += (pow(K2_12/TwoFrac*sqrt(TERM2),2) + pow((1-TwoFrac)*K2_12/TwoFrac*sqrt(TERM5),2));
	value_num_e += (pow(K2_13/TwoFrac*sqrt(TERM3),2) + pow((1-TwoFrac)*K2_13/TwoFrac*sqrt(TERM5),2));
	value_num_e += (pow(K2_23/TwoFrac*sqrt(TERM4),2) + pow((1-TwoFrac)*K2_23/TwoFrac*sqrt(TERM5),2));
	value_num_e += pow(2*sqrt(TERM5),2);
	
	C3_e[Q3bin-1] += N3_QS_e;
	c3_e[Q3bin-1] += value_num_e + TERM5;// add baseline stat error


	// Fill histograms
	C3_hist->Fill(Q3, N3_QS);
	c3_hist->Fill(Q3, value_num + TERM5);// for cumulant correlation function
	Combinatorics_1d->Fill(Q3, TERM5);
	
	//

	a_3[i-1][j-1][k-1] = value_num + TERM5;
	b_3[i-1][j-1][k-1] = TERM5;
	a_3_e[i-1][j-1][k-1] = sqrt(value_num_e + TERM5);
	
	A_3[i-1][j-1][k-1] = N3_QS;
	B_3[i-1][j-1][k-1] = TERM5;
	A_3_e[i-1][j-1][k-1] = sqrt(N3_QS_e);

	/*int fst=4, scd=5, trd=6;
	if(i==fst && j==scd && k==trd) cout<<A_3[i-1][j-1][k-1]<<"  "<<A_3_e[i-1][j-1][k-1]<<endl;
	if(i==scd && j==fst && k==trd) cout<<A_3[i-1][j-1][k-1]<<"  "<<A_3_e[i-1][j-1][k-1]<<endl;
	if(i==trd && j==scd && k==fst) cout<<A_3[i-1][j-1][k-1]<<"  "<<A_3_e[i-1][j-1][k-1]<<endl;
	if(i==fst && j==trd && k==scd) cout<<A_3[i-1][j-1][k-1]<<"  "<<A_3_e[i-1][j-1][k-1]<<endl;
	if(i==trd && j==fst && k==scd) cout<<A_3[i-1][j-1][k-1]<<"  "<<A_3_e[i-1][j-1][k-1]<<endl;
	if(i==scd && j==trd && k==fst) cout<<A_3[i-1][j-1][k-1]<<"  "<<A_3_e[i-1][j-1][k-1]<<endl;*/
	///////////////////////////////////////////////////////////
	
      }
    }
  }

  // symmetrize i-j-k
  /*for(int i=2; i<=ThreeParticle[0][0][0][0]->GetNbinsX(); i++){// bin number
    for(int j=i; j<=ThreeParticle[0][0][0][0]->GetNbinsY(); j++){// bin number
      for(int k=j; k<=ThreeParticle[0][0][0][0]->GetNbinsZ(); k++){// bin number
	a_3[i-1][j-1][k-1] = (a_3[i-1][j-1][k-1] + a_3[i-1][k-1][j-1] + a_3[j-1][i-1][k-1] + a_3[k-1][j-1][i-1]);
	a_3[i-1][j-1][k-1] += (a_3[k-1][i-1][j-1] + a_3[j-1][k-1][i-1]);
	a_3[i-1][j-1][k-1] /= 6.;
	b_3[i-1][j-1][k-1] = (b_3[i-1][j-1][k-1] + b_3[i-1][k-1][j-1] + b_3[j-1][i-1][k-1] + b_3[k-1][j-1][i-1]);
	b_3[i-1][j-1][k-1] += (b_3[k-1][i-1][j-1] + b_3[j-1][k-1][i-1]);
	b_3[i-1][j-1][k-1] /= 6.;
	//
	int degCount=1;
	if(i!=j) {a_3_e[i-1][j-1][k-1] += pow(a_3_e[j-1][i-1][k-1],2); degCount++;}
	if(i!=k) {a_3_e[i-1][j-1][k-1] += pow(a_3_e[k-1][j-1][i-1],2); degCount++;}
	if(j!=k) {a_3_e[i-1][j-1][k-1] += pow(a_3_e[i-1][k-1][j-1],2); degCount++;}
	if(i!=j && i!=k && j!=k) {
	  a_3_e[i-1][j-1][k-1] += pow(a_3_e[k-1][i-1][j-1],2);
	  a_3_e[i-1][j-1][k-1] += pow(a_3_e[j-1][k-1][i-1],2);
	  degCount += 2;
	}
	a_3_e[i-1][j-1][k-1] = sqrt(a_3_e[i-1][j-1][k-1] / degCount);
	///////
	//
	A_3[i-1][j-1][k-1] = (A_3[i-1][j-1][k-1] + A_3[i-1][k-1][j-1] + A_3[j-1][i-1][k-1] + A_3[k-1][j-1][i-1]);
	A_3[i-1][j-1][k-1] += (A_3[k-1][i-1][j-1] + A_3[j-1][k-1][i-1]);
	A_3[i-1][j-1][k-1] /= 6.;
	B_3[i-1][j-1][k-1] = (B_3[i-1][j-1][k-1] + B_3[i-1][k-1][j-1] + B_3[j-1][i-1][k-1] + B_3[k-1][j-1][i-1]);
	B_3[i-1][j-1][k-1] += (B_3[k-1][i-1][j-1] + B_3[j-1][k-1][i-1]);
	B_3[i-1][j-1][k-1] /= 6.;
	//
	degCount=1;
	if(i!=j) {A_3_e[i-1][j-1][k-1] += pow(A_3_e[j-1][i-1][k-1],2); degCount++;}
	if(i!=k) {A_3_e[i-1][j-1][k-1] += pow(A_3_e[k-1][j-1][i-1],2); degCount++;}
	if(j!=k) {A_3_e[i-1][j-1][k-1] += pow(A_3_e[i-1][k-1][j-1],2); degCount++;}
	if(i!=j && i!=k && j!=k) {
	  A_3_e[i-1][j-1][k-1] += pow(A_3_e[k-1][i-1][j-1],2);
	  A_3_e[i-1][j-1][k-1] += pow(A_3_e[j-1][k-1][i-1],2);
	  degCount += 2;
	}
	A_3_e[i-1][j-1][k-1] = sqrt(A_3_e[i-1][j-1][k-1] / degCount);
      }
    }
    }*/


  ////////////////////////////

  // Intermediate steps
  for(int i=0; i<Q3BINS; i++) {
    C3_hist->SetBinError(i+1, sqrt(C3_e[i]));
    c3_hist->SetBinError(i+1, sqrt(c3_e[i]));
  }

  C3_hist->Divide(Combinatorics_1d);
  c3_hist->Divide(Combinatorics_1d);

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  
  
 
  c3_hist->GetXaxis()->SetTitle("#font[12]{Q}_{3} (GeV/#font[12]{c})");
  c3_hist->GetYaxis()->SetTitle("Three pion correlation");
  c3_hist->SetMarkerSize(MarkerSize);
  c3_hist->SetMarkerStyle(24);
  c3_hist->GetYaxis()->SetTitleSize(0.045); c3_hist->GetXaxis()->SetTitleSize(0.045);
  c3_hist->GetYaxis()->SetTitleOffset(0.8);
  c3_hist->GetXaxis()->SetRangeUser(0,Q3LimitHigh);
  c3_hist->GetYaxis()->SetRangeUser(0.9,4);
  c3_hist->SetMarkerColor(1);
  c3_hist->SetLineColor(1);
  c3_hist->SetMaximum(3);
  c3_hist->SetMinimum(.7);
  //

  C3_hist->GetXaxis()->SetTitle("#font[12]{Q}_{3} (GeV/#font[12]{c})");
  C3_hist->GetYaxis()->SetTitle("Three pion correlation");
  C3_hist->SetMarkerSize(MarkerSize);
  C3_hist->GetXaxis()->SetTitleSize(SizeTitle);
  C3_hist->GetXaxis()->SetLabelSize(SizeLabel);
  C3_hist->GetYaxis()->SetTitleSize(SizeTitle);
  C3_hist->GetYaxis()->SetLabelSize(SizeLabel);
  C3_hist->GetXaxis()->SetTitleOffset(1.05);
  C3_hist->GetYaxis()->SetTitleOffset(0.8);
  C3_hist->GetXaxis()->SetNdivisions(606);
  C3_hist->GetYaxis()->SetNdivisions(505);
  C3_hist->GetXaxis()->SetRangeUser(0,Q3LimitHigh);
  C3_hist->GetYaxis()->SetRangeUser(0.9,4);
  C3_hist->SetMarkerStyle(20);
  C3_hist->SetMarkerColor(4);
  C3_hist->SetLineColor(4);
  if(CollisionType==0) C3_hist->SetMaximum(4);
  else C3_hist->SetMaximum(6);
  C3_hist->SetMinimum(.8);
  //
 
  TGraph *gr_c3Spline[2];// c3/C3
  TMinuit *MyMinuit_fit[2];// c3/C3
  TF1 *c3Fit1D_Expan[2];
  if(FitType==0) {
    c3Fit1D_Expan[0]=new TF1("c3Fit1D_Expan1","[0]*(1+[1]*exp(-pow([2]*x/0.19733,2)/2.) * pow(1 + ([3]/(6.*pow(2.,1.5))*(8.*pow([2]*x/sqrt(3.)/0.19733,3) - 12.*pow([2]*x/sqrt(3.)/0.19733,1))) + ([4]/(24.*pow(2.,2))*(16.*pow([2]*x/sqrt(3.)/0.19733,4) -48.*pow([2]*x/sqrt(3.)/0.19733,2) + 12)) + [5]/(120.*pow(2.,2.5))*(32.*pow(x/sqrt(3.)*[2]/0.19733,5) - 160.*pow(x/sqrt(3.)*[2]/0.19733,3) + 120*x/sqrt(3.)*[2]/0.19733)   +  [6]/pow(2.,6/2.)/(720.)*(64*pow(x/sqrt(3.)*[2]/0.19733,6) - 480*pow(x/sqrt(3.)*[2]/0.19733,4) + 720*pow(x/sqrt(3.)*[2]/0.19733,2) - 120),3))",0,1);
    c3Fit1D_Expan[1]=new TF1("c3Fit1D_Expan2","[0]*(1+[1]*exp(-pow([2]*x/0.19733,2)/2.) * pow(1 + ([3]/(6.*pow(2.,1.5))*(8.*pow([2]*x/sqrt(3.)/0.19733,3) - 12.*pow([2]*x/sqrt(3.)/0.19733,1))) + ([4]/(24.*pow(2.,2))*(16.*pow([2]*x/sqrt(3.)/0.19733,4) -48.*pow([2]*x/sqrt(3.)/0.19733,2) + 12)) + [5]/(120.*pow(2.,2.5))*(32.*pow(x/sqrt(3.)*[2]/0.19733,5) - 160.*pow(x/sqrt(3.)*[2]/0.19733,3) + 120*x/sqrt(3.)*[2]/0.19733)  +  [6]/pow(2.,6/2.)/(720.)*(64*pow(x/sqrt(3.)*[2]/0.19733,6) - 480*pow(x/sqrt(3.)*[2]/0.19733,4) + 720*pow(x/sqrt(3.)*[2]/0.19733,2) - 120),3))",0,1);
  }else if(FitType==1){
    c3Fit1D_Expan[0]=new TF1("c3Fit1D_Expan1","[0]*(1+[1]*exp(-[2]*x/0.19733 * sqrt(3.)/2.) * pow(1 + [3]*(-[2]*x/sqrt(3.)/0.19733 + 1) + [4]/2*(pow([2]*x/sqrt(3.)/0.19733,2) - 4*[2]*x/sqrt(3.)/0.19733 + 2) + [5]/6.*(-pow(x/sqrt(3.)*[1]/0.19733,3) + 9*pow(x/sqrt(3.)*[1]/0.19733,2) - 18*x/sqrt(3.)*[1]/0.19733 + 6)  +  [6]/24.*(pow(x/sqrt(3.)*[1]/0.19733,4) - 16*pow(x/sqrt(3.)*[1]/0.19733,3) + 72*pow(x/sqrt(3.)*[1]/0.19733,2) - 96*x/sqrt(3.)*[1]/0.19733 + 24),3))",0,1);
    c3Fit1D_Expan[1]=new TF1("c3Fit1D_Expan2","[0]*(1+[1]*exp(-[2]*x/0.19733 * sqrt(3.)/2.) * pow(1 + [3]*(-[2]*x/sqrt(3.)/0.19733 + 1) + [4]/2*(pow([2]*x/sqrt(3.)/0.19733,2) - 4*[2]*x/sqrt(3.)/0.19733 + 2) + [5]/6.*(-pow(x/sqrt(3.)*[1]/0.19733,3) + 9*pow(x/sqrt(3.)*[1]/0.19733,2) - 18*x/sqrt(3.)*[1]/0.19733 + 6)  +  [6]/24.*(pow(x/sqrt(3.)*[1]/0.19733,4) - 16*pow(x/sqrt(3.)*[1]/0.19733,3) + 72*pow(x/sqrt(3.)*[1]/0.19733,2) - 96*x/sqrt(3.)*[1]/0.19733 + 24),3))",0,1);
  }else{
    c3Fit1D_Expan[0]=new TF1("c3Fit1D_Expan1","[0]*(1+[1]*exp(-pow([2]*x/0.19733, [3])))",0,1);
    c3Fit1D_Expan[1]=new TF1("c3Fit1D_Expan2","[0]*(1+[1]*exp(-pow([2]*x/0.19733, [3])))",0,1);
  }
  
  TF1 *testEq[2];
  if(FitType==0) {
    testEq[0] = new TF1("testEq_c3","[0]*exp(-pow(x*[1]/0.19733,2)/2.) * ( 1 + [2]/(6.*pow(2.,1.5))*(8*pow(x*[1]/0.19733,3) - 12*pow(x*[1]/0.19733,1)) + [3]/(24.*pow(2.,2))*(16*pow(x*[1]/0.19733,4) -48*pow(x*[1]/0.19733,2) + 12) + [4]/(120.*pow(2.,2.5))*(32.*pow(x*[1]/0.19733,5) - 160.*pow(x*[1]/0.19733,3) + 120*x*[1]/0.19733)  +  [5]/pow(2.,6/2.)/(720.)*(64*pow(x*[1]/0.19733,6) - 480*pow(x*[1]/0.19733,4) + 720*pow(x*[1]/0.19733,2) - 120) )",0,1);
    testEq[1] = new TF1("testEq_C3","[0]*exp(-pow(x*[1]/0.19733,2)/2.) * ( 1 + [2]/(6.*pow(2.,1.5))*(8*pow(x*[1]/0.19733,3) - 12*pow(x*[1]/0.19733,1)) + [3]/(24.*pow(2.,2))*(16*pow(x*[1]/0.19733,4) -48*pow(x*[1]/0.19733,2) + 12) + [4]/(120.*pow(2.,2.5))*(32.*pow(x*[1]/0.19733,5) - 160.*pow(x*[1]/0.19733,3) + 120*x*[1]/0.19733) +  [5]/pow(2.,6/2.)/(720.)*(64*pow(x*[1]/0.19733,6) - 480*pow(x*[1]/0.19733,4) + 720*pow(x*[1]/0.19733,2) - 120) )",0,1);
  }else {
    testEq[0] = new TF1("testEq_c3","[0]*exp(-x*[1]/0.19733/2.) * ( 1 + [2]*(-x*[1]/0.19733 + 1) + [3]/2.*(pow(x*[1]/0.19733,2) - 4*x*[1]/0.19733 + 2) + [4]/6.*(-pow(x*[1]/0.19733,3) + 9*pow(x*[1]/0.19733,2) - 18*x*[1]/0.19733 + 6)  + [5]/24.*(pow(x*[1]/0.19733,4) - 16*pow(x*[1]/0.19733,3) + 72*pow(x*[1]/0.19733,2) - 96*x*[1]/0.19733 + 24))",0,1);
    testEq[1] = new TF1("testEq_C3","[0]*exp(-x*[1]/0.19733/2.) * ( 1 + [2]*(-x*[1]/0.19733 + 1) + [3]/2.*(pow(x*[1]/0.19733,2) - 4*x*[1]/0.19733 + 2) + [4]/6.*(-pow(x*[1]/0.19733,3) + 9*pow(x*[1]/0.19733,2) - 18*x*[1]/0.19733 + 6)  + [5]/24.*(pow(x*[1]/0.19733,4) - 16*pow(x*[1]/0.19733,3) + 72*pow(x*[1]/0.19733,2) - 96*x*[1]/0.19733 + 24))",0,1);
  }
  
  
  for(int fitIt=0; fitIt<2; fitIt++){
    if(fitIt==0) CumulantFit=kTRUE;
    else CumulantFit=kFALSE;
    
    const int npar_c3=8;
    TMinuit MyMinuit_c3(npar_c3);
    MyMinuit_c3.SetFCN(fcn_c3);
    int ierflg_c3=0;
    double arglist_c3 = 0;
    MyMinuit_c3.mnexcm("SET NOWarnings",&arglist_c3,1, ierflg_c3);  
    arglist_c3 = -1;
    MyMinuit_c3.mnexcm("SET PRint",&arglist_c3,1, ierflg_c3);
    //arglist_c3=2;// improve Minimization Strategy (1 is default)
    //MyMinuit_c3.mnexcm("SET STR",&arglist_c3,1,ierflg_c3);
    //arglist_c3 = 0;
    //MyMinuit_c3.mnexcm("SCAN", &arglist_c3,1,ierflg_c3);
    arglist_c3 = 5000;
    MyMinuit_c3.mnexcm("MIGRAD", &arglist_c3 ,1,ierflg_c3);
    
   

  double OutputPar_c3[npar_c3]={0};
  double OutputPar_c3_e[npar_c3]={0};
  
  double par_c3[npar_c3];               // the start values
  double stepSize_c3[npar_c3];          // step sizes 
  double minVal_c3[npar_c3];            // minimum bound on parameter 
  double maxVal_c3[npar_c3];            // maximum bound on parameter
  string parName_c3[npar_c3];
  //          1.0              1.0              10.              0.              0.              0.             1.5  
  par_c3[0] = 1.0; par_c3[1] = 1.0; par_c3[2] = 10.; par_c3[3] = 0.; par_c3[4] = 0.; par_c3[5] = 0; par_c3[6] = 0; par_c3[7] = 1.5;
  stepSize_c3[0] = 0.01; stepSize_c3[1] = 0.1; stepSize_c3[2] = 0.1; stepSize_c3[3] = 0.01; stepSize_c3[4] = 0.01; stepSize_c3[5] = 0.01; stepSize_c3[6] = 0.01; stepSize_c3[7] = 0.1;
  minVal_c3[0] = 0.8; minVal_c3[1] = 0.2; minVal_c3[2] = 4.; minVal_c3[3] = -2; minVal_c3[4] = -2; minVal_c3[5] = -5; minVal_c3[6] = -5; minVal_c3[6] = 0.5;
  maxVal_c3[0] = 1.1; maxVal_c3[1] = 1000.; maxVal_c3[2] = 100.; maxVal_c3[3] = +2; maxVal_c3[4] = +2; maxVal_c3[5] = +5; maxVal_c3[6] = +5; maxVal_c3[7] = 2.5;
  parName_c3[0] = "N"; parName_c3[1] = "s"; parName_c3[2] = "R_{inv}"; parName_c3[3] = "coeff_{3}"; parName_c3[4] = "coeff_{4}"; parName_c3[5] = "coeff_{5}"; parName_c3[6] = "coeff_{6}"; parName_c3[7] = "#alpha";

  if(CollisionType==0){ 
    if(FitType!=0) {
      par_c3[2]=15.; minVal_c3[2] = 8.;
    }
  }else{
    if(FitType==0) {par_c3[2] = 2.0; minVal_c3[2] = 1.0;}
    else {
      par_c3[1] = 4.0; minVal_c3[1] = .2; maxVal_c3[1] = 3.0;
      //
      par_c3[2] = 1.3; minVal_c3[2] = .9; maxVal_c3[2] = 4.;
      minVal_c3[3] = -1; minVal_c3[4] = -1;
      maxVal_c3[3] = 1; maxVal_c3[4] = 1;
    }
  }
  
  if(FitType==0) {par_c3[7] = 2.;}
  if(FitType==1) {par_c3[7] = 1.;}
  if(FitType==2) {par_c3[3] = 0; par_c3[4] = 0; par_c3[5] = 0; par_c3[6] = 0;}

  if(FitType==2) {maxVal_c3[1] = 10.;}

  for (int i=0; i<npar_c3; i++){
    MyMinuit_c3.DefineParameter(i, parName_c3[i].c_str(), par_c3[i], stepSize_c3[i], minVal_c3[i], maxVal_c3[i]);
  }
  if(FitType==0 || FitType==1) { MyMinuit_c3.FixParameter(7);}
  if(FitType==2){
    MyMinuit_c3.FixParameter(3);
    MyMinuit_c3.FixParameter(4);
    MyMinuit_c3.FixParameter(5);
    MyMinuit_c3.FixParameter(6);
  }
  if(FixNorm) MyMinuit_c3.FixParameter(0);
  //MyMinuit_c3.FixParameter(1);
  //MyMinuit_c3.FixParameter(2);
  //MyMinuit_c3.FixParameter(3);
  //MyMinuit_c3.FixParameter(4);
  //MyMinuit_c3.FixParameter(5);
  //MyMinuit_c3.FixParameter(6);
  /////////////////////////////////////////////////////////////
  // Do the minimization!
  cout<<"Start Three-d fit"<<endl;
  MyMinuit_c3.Migrad();// Minuit's best minimization algorithm
  cout<<"End Three-d fit"<<endl;
  /////////////////////////////////////////////////////////////
  MyMinuit_c3.mnexcm("SHOw PARameters", &arglist_c3, 1, ierflg_c3);
  if(fitIt==0) cout<<"c3 Fit: Chi2 = ";
  else cout<<"C3 Fit: Chi2 = ";
  cout<<Chi2_c3<<"   NDF = "<<(NFitPoints_c3-MyMinuit_c3.GetNumFreePars())<<endl;
  cout<<" Chi2/NDF = "<<Chi2_c3 / (NFitPoints_c3-MyMinuit_c3.GetNumFreePars())<<endl;

  for(int i=0; i<npar_c3; i++){
    MyMinuit_c3.GetParameter(i,OutputPar_c3[i],OutputPar_c3_e[i]);
  }
  
  
  if(FitType!=2){
    c3Fit1D_Expan[fitIt]->SetParameter(0,OutputPar_c3[0]);
    c3Fit1D_Expan[fitIt]->SetParameter(1,OutputPar_c3[1]);
    c3Fit1D_Expan[fitIt]->SetParameter(2,OutputPar_c3[2]);
    c3Fit1D_Expan[fitIt]->SetParameter(3,OutputPar_c3[3]);
    c3Fit1D_Expan[fitIt]->SetParameter(4,OutputPar_c3[4]);
    c3Fit1D_Expan[fitIt]->SetParameter(5,OutputPar_c3[5]);
    c3Fit1D_Expan[fitIt]->SetParameter(6,OutputPar_c3[6]);
  }else{// Levy
    c3Fit1D_Expan[fitIt]->SetParameter(0,OutputPar_c3[0]);
    c3Fit1D_Expan[fitIt]->SetParameter(1,OutputPar_c3[1]);
    c3Fit1D_Expan[fitIt]->SetParameter(2,OutputPar_c3[2]);
    c3Fit1D_Expan[fitIt]->SetParameter(3,OutputPar_c3[6]);
  }
  c3Fit1D_Expan[fitIt]->SetLineStyle(1);
  //c3Fit1D_Expan->Draw("same");
  
  for(int i=0; i<6; i++) testEq[fitIt]->SetParameter(i,OutputPar_c3[i+1]);
  cout<<"EA = "<<testEq[fitIt]->Eval(0.025)<<endl;

  // project 3D EW fit onto 1D histogram
  for(int i=2; i<=BINLIMIT_3; i++){// bin number
    double Q12 = BinCenters[i-1];// true center
    for(int j=2; j<=BINLIMIT_3; j++){// bin number
      double Q13 = BinCenters[j-1];// true center
      for(int k=2; k<=BINLIMIT_3; k++){// bin number
	double Q23 = BinCenters[k-1];// true center
	//	
	double Q3 = sqrt(pow(Q12,2) + pow(Q13,2) + pow(Q23,2));
	
	if(Q12 < sqrt(pow(Q13,2)+pow(Q23,2) - 2*Q13*Q23)) continue;// not all configurations are possible
	if(Q12 > sqrt(pow(Q13,2)+pow(Q23,2) + 2*Q13*Q23)) continue;// not all configurations are possible
	
	double TERM5=ThreeParticle[ch1][ch2][ch3][4]->GetBinContent(i,j,k);// N1*N1*N1. All from different events (pure combinatorics)
	if(TERM5==0) continue;
	
	
	if(MRC) TERM5 *= MRC_SC_3[2]->GetBinContent(MRC_SC_3[2]->GetXaxis()->FindBin(Q3));
	//
	double radius = OutputPar_c3[2]/FmToGeV;
	double arg12 = Q12*radius;
	double arg13 = Q13*radius;
	double arg23 = Q23*radius;
	double Expan12=1, Expan13=1, Expan23=1;
	if(FitType==0){
	  Expan12 += OutputPar_c3[3]/pow(2.,3/2.)/(6.)*(8*pow(arg12,3) - 12*pow(arg12,1));
	  Expan12 += OutputPar_c3[4]/pow(2.,4/2.)/(24.)*(16*pow(arg12,4) -48*pow(arg12,2) + 12);
	  Expan12 += OutputPar_c3[5]/pow(2.,5/2.)/(120.)*(32*pow(arg12,5) - 160*pow(arg12,3) + 120*arg12);
	  Expan12 += OutputPar_c3[6]/pow(2.,6/2.)/(720.)*(64*pow(arg12,6) - 480*pow(arg12,4) + 720*pow(arg12,2) - 120);
	  //
	  Expan13 += OutputPar_c3[3]/pow(2.,3/2.)/(6.)*(8*pow(arg13,3) - 12*pow(arg13,1));
	  Expan13 += OutputPar_c3[4]/pow(2.,4/2.)/(24.)*(16*pow(arg13,4) -48*pow(arg13,2) + 12);
	  Expan13 += OutputPar_c3[5]/pow(2.,5/2.)/(120.)*(32*pow(arg13,5) - 160*pow(arg13,3) + 120*arg13);
	  Expan13 += OutputPar_c3[6]/pow(2.,6/2.)/(720.)*(64*pow(arg13,6) - 480*pow(arg13,4) + 720*pow(arg13,2) - 120);
	  //
	  Expan23 += OutputPar_c3[3]/pow(2.,3/2.)/(6.)*(8*pow(arg23,3) - 12*pow(arg23,1));
	  Expan23 += OutputPar_c3[4]/pow(2.,4/2.)/(24.)*(16*pow(arg23,4) -48*pow(arg23,2) + 12);
	  Expan23 += OutputPar_c3[5]/pow(2.,5/2.)/(120.)*(32*pow(arg23,5) - 160*pow(arg23,3) + 120*arg23);
	  Expan23 += OutputPar_c3[6]/pow(2.,6/2.)/(720.)*(64*pow(arg23,6) - 480*pow(arg23,4) + 720*pow(arg23,2) - 120);
	}else if(FitType==1){
	  Expan12 += OutputPar_c3[3]*(-arg12 + 1);
	  Expan12 += OutputPar_c3[4]/2.*(pow(arg12,2) - 4*arg12 + 2);
	  Expan12 += OutputPar_c3[5]/6.*(-pow(arg12,3) + 9*pow(arg12,2) - 18*arg12 + 6);
	  Expan12 += OutputPar_c3[6]/24.*(pow(arg12,4) - 16*pow(arg12,3) + 72*pow(arg12,2) - 96*arg12 + 24);
	  //
	  Expan13 += OutputPar_c3[3]*(-arg13 + 1);
	  Expan13 += OutputPar_c3[4]/2.*(pow(arg13,2) - 4*arg13 + 2);
	  Expan13 += OutputPar_c3[5]/6.*(-pow(arg13,3) + 9*pow(arg13,2) - 18*arg13 + 6);
	  Expan13 += OutputPar_c3[6]/24.*(pow(arg13,4) - 16*pow(arg13,3) + 72*pow(arg13,2) - 96*arg13 + 24);
	  //
	  Expan23 += OutputPar_c3[3]*(-arg23 + 1);
	  Expan23 += OutputPar_c3[4]/2.*(pow(arg23,2) - 4*arg23 + 2);
	  Expan23 += OutputPar_c3[5]/6.*(-pow(arg23,3) + 9*pow(arg23,2) - 18*arg23 + 6);
	  Expan23 += OutputPar_c3[6]/24.*(pow(arg23,4) - 16*pow(arg23,3) + 72*pow(arg23,2) - 96*arg23 + 24);
	}else{}
	
	//
	double t12_coh = exp(-pow(Rcoh/FmToGeV * Q12,2)/2.);
	double t23_coh = exp(-pow(Rcoh/FmToGeV * Q23,2)/2.);
	double t13_coh = exp(-pow(Rcoh/FmToGeV * Q13,2)/2.);
	if(Rcoh>=100){
	  t12_coh = exp(-pow(arg12,OutputPar_c3[7])/2.) * Expan12;
	  t23_coh = exp(-pow(arg23,OutputPar_c3[7])/2.) * Expan23;
	  t13_coh = exp(-pow(arg13,OutputPar_c3[7])/2.) * Expan13;
	}
	double C = 2*pow(OutputPar_c3[1],3) * pow(1-G,3)*exp(-(pow(arg12,OutputPar_c3[7])+pow(arg13,OutputPar_c3[7])+pow(arg23,OutputPar_c3[7]))/2.)*Expan12*Expan13*Expan23;
	C += 2*pow(OutputPar_c3[1], 2) * G*pow(1-G,2)*exp(-(pow(arg12,OutputPar_c3[7])+pow(arg13,OutputPar_c3[7]))/2.)*Expan12*Expan13*t23_coh;
	C += 2*pow(OutputPar_c3[1], 2) * G*pow(1-G,2)*exp(-(pow(arg12,OutputPar_c3[7])+pow(arg23,OutputPar_c3[7]))/2.)*Expan12*Expan23*t13_coh;
	C += 2*pow(OutputPar_c3[1], 2) * G*pow(1-G,2)*exp(-(pow(arg13,OutputPar_c3[7])+pow(arg23,OutputPar_c3[7]))/2.)*Expan13*Expan23*t12_coh;
	if(!CumulantFit){
	  C += pow(OutputPar_c3[1],2) * pow(1-G,2)*exp(-pow(arg12,OutputPar_c3[7]))*pow(Expan12,2);
	  C += pow(OutputPar_c3[1],2) * pow(1-G,2)*exp(-pow(arg13,OutputPar_c3[7]))*pow(Expan13,2);
	  C += pow(OutputPar_c3[1],2) * pow(1-G,2)*exp(-pow(arg23,OutputPar_c3[7]))*pow(Expan23,2);
	  C += 2*OutputPar_c3[1] * G*(1-G)*exp(-pow(arg12,OutputPar_c3[7])/2.)*Expan12 * t12_coh;
	  C += 2*OutputPar_c3[1] * G*(1-G)*exp(-pow(arg13,OutputPar_c3[7])/2.)*Expan13 * t13_coh;
	  C += 2*OutputPar_c3[1] * G*(1-G)*exp(-pow(arg23,OutputPar_c3[7])/2.)*Expan23 * t23_coh;
	}
	C += 1.0;
	C *= TERM5;
	C *= OutputPar_c3[0];
	//if(Q3<0.018) continue;
	GenSignalExpected_num[fitIt]->Fill(Q3, C);
	GenSignalExpected_den[fitIt]->Fill(Q3, TERM5);
	//if(Q3<0.02) cout<<Q3<<"  "<<TERM5<<endl;
	///////////////////////////////////////////////////////////
	
	
	
      }
    }
  }
  

  GenSignalExpected_num[fitIt]->Sumw2();
  GenSignalExpected_num[fitIt]->Divide(GenSignalExpected_den[fitIt]);
  

  TSpline3 *c3Fit1D_ExpanSpline = new TSpline3(GenSignalExpected_num[fitIt]);
  c3Fit1D_ExpanSpline->SetLineWidth(2);
  const int Npoints=500;
  double xpoints[Npoints], ypoints[Npoints];
  bool splineOnly=kFALSE;
  double Q3Start=0.01;
  if(CollisionType!=0) Q3Start = 0.02;
  for(int iii=0; iii<Npoints; iii++){
    if(CollisionType==0) xpoints[iii] = Q3Start + (iii+0.5)*0.01;
    else xpoints[iii] = Q3Start + (iii+0.5)*0.03;
    //
    //ypoints[iii] = c3Fit1D_ExpanSpline->Eval(xpoints[iii]);// to skip spline
    splineOnly=kTRUE;// to skip 1D approximation
    //if(CollisionType==0) splineOnly=kTRUE;
    //if(CollisionType!=0 && xpoints[iii] > 0.06) splineOnly=kTRUE;
    if(!splineOnly){
      ypoints[iii] = c3Fit1D_Expan[fitIt]->Eval(xpoints[iii]);
      if(c3Fit1D_Expan[fitIt]->Eval(xpoints[iii])<2. && fabs(c3Fit1D_Expan[fitIt]->Eval(xpoints[iii])-c3Fit1D_ExpanSpline->Eval(xpoints[iii])) < 0.01) splineOnly=kTRUE;
    }
    else {ypoints[iii] = c3Fit1D_ExpanSpline->Eval(xpoints[iii]); splineOnly=kTRUE;}
    //
    if(CollisionType==0){
      int Q3Bin = int(xpoints[iii]/0.01) + 1;
      if(Q3Bin<=10){// to fix the Q3 problem with offline 1D projections of 3D histo
	if(fitIt==0) ypoints[iii] *= c3_pointsPbPb_accepted[Q3Bin-1] / c3_hist->GetBinContent(Q3Bin);
	else ypoints[iii] *= C3_pointsPbPb_accepted[Q3Bin-1] / C3_hist->GetBinContent(Q3Bin);
      }
    } 
  }
  gr_c3Spline[fitIt] = new TGraph(Npoints,xpoints,ypoints);
  gr_c3Spline[fitIt]->SetLineWidth(2); 
  if(fitIt==0) {
    gr_c3Spline[fitIt]->SetLineColor(1);
    gr_c3Spline[fitIt]->SetLineStyle(4);
  }else{
    gr_c3Spline[fitIt]->SetLineColor(4);
  }
  MyMinuit_fit[fitIt] = (TMinuit*)MyMinuit_c3.Clone();
  
  }// fitIt

  double CorrectionFactor[2][10]={{0}};
  if(CollisionType==0){
    for(int i=0; i<10; i++){
      if(c3_hist->GetBinContent(i+1) > 0) CorrectionFactor[0][i] = c3_pointsPbPb_accepted[i] / c3_hist->GetBinContent(i+1);
      if(C3_hist->GetBinContent(i+1) > 0) CorrectionFactor[1][i] = C3_pointsPbPb_accepted[i] / C3_hist->GetBinContent(i+1);
      C3_hist->SetBinContent(i+1, C3_pointsPbPb_accepted[i]);
      c3_hist->SetBinContent(i+1, c3_pointsPbPb_accepted[i]);
    }
  }

  if(CollisionType!=0){
    int binH=4;
    if(CollisionType==2) binH=5;
    for(int bin=1; bin<=binH; bin++){
      C3_hist->SetBinContent(bin,100);
      c3_hist->SetBinContent(bin,100);
    }
  }
  C3_hist->Draw();
  c3_hist->Draw("same");
  legend1->AddEntry(C3_hist,"#font[12]{C}_{3}^{QS}","p");
  legend1->AddEntry(c3_hist,"#font[12]{c}_{3}^{QS}","p");
 
  Unity->Draw("same");
  
  ALICEspecif->Draw("same");
  KTspecif->Draw("same");
  

  gr_c3Spline[0]->Draw("same");
  gr_c3Spline[1]->Draw("same");
  // Systematics
  TH1D *C3QS_Syst = (TH1D*)C3_hist->Clone();
  C3QS_Syst->SetBinContent(1,100); 
  C3QS_Syst->SetMarkerSize(0); C3QS_Syst->SetFillColor(kBlue-10);
  C3QS_Syst->SetLineColor(kBlue-10); C3QS_Syst->SetLineWidth(5);
  C3QS_Syst->SetMarkerColor(kBlue-10);
  for(int bin=1; bin<=C3QS_Syst->GetNbinsX(); bin++){
    double q3 = C3QS_Syst->GetXaxis()->GetBinCenter(bin);
    if(CollisionType!=0) q3 /= 2.5;// Scale factor for pp and p-Pb
    double syst1 = pow(0.001,2);// cc
    syst1 += pow(0.002 - 0.002*q3/0.1,2);// 11h to 10h
    syst1 += pow(0.01*(1 - q3/0.1),2);// f coefficients, r*<70. was pow(0.9913 - 0.2231*q3 - 1,2)
    syst1 += pow(0.9847 + 0.358*q3 - 2.133*q3*q3 - 1,2);// MRC
    syst1 += pow(0.975 + 0.4189*q3 - 2.055*q3*q3 - 1,2);// Muon, 92%
    syst1 += pow(0.936 + 1.194*q3 - 5.912*q3*q3 - 1,2);// fc2 scale
    syst1 = sqrt(syst1);
    syst1 *= C3_hist->GetBinContent(bin);
    C3QS_Syst->SetBinError(bin, syst1);
  }  
  
  C3QS_Syst->Draw("same");
  C3_hist->Draw("same");
  //
  legend1->AddEntry(gr_c3Spline[0],"Fit to #font[12]{c}_{3}^{QS}","l");
  legend1->AddEntry(gr_c3Spline[1],"Fit to #font[12]{C}_{3}^{QS}","l");
  legend1->Draw("same");

  ChargeCombTitle3->Draw("same");

  //cout<<c3_hist->GetBinContent(3)<<"  "<<C3_hist->GetBinContent(3)<<endl;
  // Ratio plot
  pad1_2->cd();
  gPad->SetLeftMargin(0.11); gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.03); gPad->SetBottomMargin(0.32);
  TH1D *Ratio_C3 = (TH1D*)C3_hist->Clone();
  TH1D *Ratio_c3 = (TH1D*)c3_hist->Clone();
  Ratio_C3->GetYaxis()->SetTitle("Data / fit   ");
  Ratio_C3->GetXaxis()->SetTitleSize(2.2*SizeTitle);
  Ratio_C3->GetXaxis()->SetLabelSize(2.2*SizeLabel);
  Ratio_C3->GetYaxis()->SetTitleSize(2.2*SizeTitle);
  Ratio_C3->GetYaxis()->SetLabelSize(2.2*SizeLabel);
  Ratio_C3->GetXaxis()->SetTitleOffset(1.05);
  Ratio_C3->GetYaxis()->SetTitleOffset(0.4);
  Ratio_C3->GetXaxis()->SetNdivisions(606);
  Ratio_C3->GetYaxis()->SetNdivisions(202);
  Ratio_C3->SetMinimum(0.96); Ratio_C3->SetMaximum(1.054);
  if(CollisionType!=0) {Ratio_C3->SetMinimum(0.86); Ratio_C3->SetMaximum(1.056);}
  int BinKill=1;
  if(CollisionType!=0) BinKill=10;
  for(int bin=1; bin<=C3_hist->GetNbinsX()-BinKill; bin++){
    //double q3 = Ratio_c3->GetXaxis()->GetBinCenter(bin);
    double value = Ratio_c3->GetBinContent(bin);
    double value_e = Ratio_c3->GetBinError(bin);
    if(value < 0.9) continue;
    if(GenSignalExpected_num[0]->GetBinContent(bin) < 0.9) continue;
    if(GenSignalExpected_num[1]->GetBinContent(bin) < 0.9) continue;
    value /= GenSignalExpected_num[0]->GetBinContent(bin);
    value_e /= GenSignalExpected_num[0]->GetBinContent(bin);
    if(CorrectionFactor[0][bin-1]>0 && CollisionType==0) {
      value /= CorrectionFactor[0][bin-1];
      value_e /=  CorrectionFactor[0][bin-1];
    }
    Ratio_c3->SetBinContent(bin, value);
    Ratio_c3->SetBinError(bin, value_e);
    //
    value = Ratio_C3->GetBinContent(bin);
    value_e = Ratio_C3->GetBinError(bin);
    value /= GenSignalExpected_num[1]->GetBinContent(bin);
    value_e /= GenSignalExpected_num[1]->GetBinContent(bin);
    if(CorrectionFactor[1][bin-1]>0 && CollisionType==0) {
      value /= CorrectionFactor[1][bin-1];
      value_e /=  CorrectionFactor[1][bin-1];
    }
    Ratio_C3->SetBinContent(bin, value);
    Ratio_C3->SetBinError(bin, value_e);
  }
  Ratio_C3->Draw();
  Ratio_c3->Draw("same");
  legend2->AddEntry(Ratio_C3,"#font[12]{C}_{3}^{QS}","p");
  legend2->AddEntry(Ratio_c3,"#font[12]{c}_{3}^{QS}","p");
  //legend2->Draw("same");
  Unity->Draw("same");


  if(PrintData){
    int Q3bins=9;
    if(CollisionType!=0) Q3bins=25;
    cout.precision(4);
    cout<<"C3QS  C3QS/Fit  c3QS  c3QS/Fit"<<endl;
    for(int qbin=1; qbin<=Q3bins; qbin++){
      if(C3_hist->GetBinContent(qbin)>99) continue;
      if(C3_hist->GetBinContent(qbin)==0) continue;
      cout<<C3_hist->GetXaxis()->GetBinLowEdge(qbin)<<" TO "<<C3_hist->GetXaxis()->GetBinUpEdge(qbin)<<"; "<<C3_hist->GetBinContent(qbin)<<" +- "<<C3_hist->GetBinError(qbin)<<" (DSYS=+-"<<C3QS_Syst->GetBinError(qbin)<<"); ";
      cout<<Ratio_C3->GetBinContent(qbin)<<" +- "<<Ratio_C3->GetBinError(qbin)<<"; ";
      cout<<c3_hist->GetBinContent(qbin)<<" +- "<<c3_hist->GetBinError(qbin)<<" (DSYS=+-"<<C3QS_Syst->GetBinError(qbin) * c3_hist->GetBinContent(qbin)/C3_hist->GetBinContent(qbin)<<"); ";
      cout<<Ratio_c3->GetBinContent(qbin)<<" +- "<<Ratio_c3->GetBinError(qbin)<<";"<<endl;
    }
    cout<<endl;
  }


  /*
  double ChiSqSum_1D=0, SumNDF_1D=0;
  for(int bin=1; bin<=300; bin++){
    double binCenter = c3_hist->GetXaxis()->GetBinCenter(bin);
    if(binCenter > Q3Limit) continue;
    if(c3_hist->GetBinError(bin)==0) continue;
    if(binCenter < 0.01) continue;
    int grIndex=1;
    for(int gr=0; gr<999; gr++){
      if(binCenter > xpoints[gr] && (binCenter < xpoints[gr+1])) {grIndex=gr; break;}
    }

    ChiSqSum_1D += pow((c3_hist->GetBinContent(bin)-ypoints[grIndex]) / c3_hist->GetBinError(bin),2);
    //cout<<c3_hist->GetBinContent(bin)<<"  "<<ypoints[grIndex]<<"  "<<c3_hist->GetBinError(bin)<<endl;
    cout<<pow((c3_hist->GetBinContent(bin)-ypoints[grIndex]) / c3_hist->GetBinError(bin),2)<<endl;
    SumNDF_1D++;
  }
  cout<<"1D Chi2/NDF = "<<ChiSqSum_1D / (SumNDF_1D-5.)<<endl;
  */
  
  TCanvas *can2 = new TCanvas("can2", "can2",800,0,700,700);// 11,53,700,500
  can2->SetHighLightColor(2);
  gStyle->SetOptFit(0111);
  
  if(CollisionType==0) testEq[0]->GetXaxis()->SetRangeUser(0,0.1);
  else testEq[0]->GetXaxis()->SetRangeUser(0,0.5);
  testEq[0]->Draw();
  testEq[1]->SetLineColor(4);
  testEq[1]->Draw("same");

  if(SaveToFile){
    TString *savefilename = new TString("FitFiles/");
    savefilename->Append("FitFile_CT");
    *savefilename += CollisionType;
    savefilename->Append("_FT");
    *savefilename += FitType;
    savefilename->Append("_R");
    *savefilename += Rcoh;
    savefilename->Append("_G");
    *savefilename += int((G+0.001)/0.02);
    savefilename->Append(".root");
    TFile *savefile = new TFile(savefilename->Data(),"RECREATE");
    MyMinuit_fit[0]->Write("MyMinuit_c3");
    MyMinuit_fit[1]->Write("MyMinuit_C3");
    //
    savefile->Close();
  }
  
  
}

//________________________________________________________________________
void SetFSICorrelations(){
  // read in 2-particle and 3-particle FSI correlations = K2 & K3
  // 2-particle input histo from file is binned in qinv.  3-particle in qinv of each pair
  TFile *fsifile = new TFile("KFile.root","READ");
  if(!fsifile->IsOpen()) {
    cout<<"No FSI file found!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
  }
  
  TH1D *temphistoSS[12];
  TH1D *temphistoOS[12];
  for(Int_t MB=0; MB<12; MB++) {
    TString *nameK2SS = new TString("K2ss_");
    *nameK2SS += MB;
    temphistoSS[MB] = (TH1D*)fsifile->Get(nameK2SS->Data());
    //
    TString *nameK2OS = new TString("K2os_");
    *nameK2OS += MB;
    temphistoOS[MB] = (TH1D*)fsifile->Get(nameK2OS->Data());
    //
    fFSIss[MB] = (TH1D*)temphistoSS[MB]->Clone();
    fFSIos[MB] = (TH1D*)temphistoOS[MB]->Clone();
    fFSIss[MB]->SetDirectory(0);
    fFSIos[MB]->SetDirectory(0);
  }
  //
  
  fsifile->Close();
  
}

double Gamov(int chargeProduct, double qinv){
  
  double arg = chargeProduct*2.*PI/(BohrR*qinv/2.);
  
  return arg/(exp(arg)-1);
}

void SetMomResCorrections(){
 
  TString *momresfilename = new TString("MomResFile");
  if(CollisionType_def!=0) momresfilename->Append("_ppAndpPb");
  momresfilename->Append(".root");
  
  TFile *MomResFile = new TFile(momresfilename->Data(),"READ");
  
  TString *proName[28];
  for(int ii=0; ii<28; ii++){
    proName[ii] = new TString("MRC_pro_");
    *proName[ii] += ii;
  }

  
  //
  MRC_SC_3[0] = (TH1D*)(((TH2D*)MomResFile->Get("MRC_3_SC_term1"))->ProjectionY(proName[4]->Data(), RbinMRC, RbinMRC));
  MRC_SC_3[1] = (TH1D*)(((TH2D*)MomResFile->Get("MRC_3_SC_term2"))->ProjectionY(proName[5]->Data(), RbinMRC, RbinMRC));
  MRC_SC_3[2] = (TH1D*)(((TH2D*)MomResFile->Get("MRC_3_SC_term3"))->ProjectionY(proName[6]->Data(), RbinMRC, RbinMRC));
  MRC_SC_3[0]->SetDirectory(0); MRC_SC_3[1]->SetDirectory(0); MRC_SC_3[2]->SetDirectory(0);
  //
  
  if(!MRC){
    for(int bin=1; bin<=MRC_SC_3[0]->GetNbinsX(); bin++){
      MRC_SC_3[0]->SetBinContent(bin, 1.0); MRC_SC_3[1]->SetBinContent(bin, 1.0); MRC_SC_3[2]->SetBinContent(bin, 1.0);
    }
    
  }
  MomResFile->Close();
  
}


float fact(float n){
  return (n < 1.00001) ? 1 : fact(n - 1) * n;
}
//________________________________________________________________________________________
void SetMuonCorrections(){
  TString *name = new TString("MuonCorrection");
  if(CollisionType_def!=0) name->Append("_ppAndpPb");
  
  name->Append(".root");
  TFile *MuonFile=new TFile(name->Data(),"READ");
  TString *proName[22];
  for(int ii=0; ii<22; ii++){
    proName[ii] = new TString("MuonCorr_pro_");
    *proName[ii] += ii;
  }
 
  //
  C3muonCorrectionSC[0] = (TH1D*)(((TH2D*)MuonFile->Get("C3muonCorrection_SC_term1"))->ProjectionY(proName[3]->Data(), RbinMRC, RbinMRC));
  C3muonCorrectionSC[1] = (TH1D*)(((TH2D*)MuonFile->Get("C3muonCorrection_SC_term2"))->ProjectionY(proName[4]->Data(), RbinMRC, RbinMRC));
  
  C3muonCorrectionSC[0]->SetDirectory(0); C3muonCorrectionSC[1]->SetDirectory(0);
 
  //
  if(!MuonCorrection){
    for(int bin=1; bin<=C3muonCorrectionSC[0]->GetNbinsX(); bin++){
      C3muonCorrectionSC[0]->SetBinContent(bin, 1.0); C3muonCorrectionSC[1]->SetBinContent(bin, 1.0);
    }
  }
  
  MuonFile->Close();
}
//________________________________________________________________________
void SetFSIindex(Float_t R){
  if(CollisionType_def==0){
    if(Mbin==0) fFSIindex = 0;//0-5%
    else if(Mbin==1) fFSIindex = 1;//5-10%
    else if(Mbin<=3) fFSIindex = 2;//10-20%
    else if(Mbin<=5) fFSIindex = 3;//20-30%
    else if(Mbin<=7) fFSIindex = 4;//30-40%
    else if(Mbin<=9) fFSIindex = 5;//40-50%
    else if(Mbin<=12) fFSIindex = 6;//40-50%
    else if(Mbin<=15) fFSIindex = 7;//40-50%
    else if(Mbin<=18) fFSIindex = 8;//40-50%
    else fFSIindex = 8;//90-100%
  }else fFSIindex = 9;// pp and pPb
  
}
//________________________________________________________________________
Float_t FSICorrelation(Int_t charge1, Int_t charge2, Float_t qinv){
  // returns 2-particle Coulomb correlations = K2
  Int_t qbinL = fFSIss[fFSIindex]->GetXaxis()->FindBin(qinv-fFSIss[fFSIindex]->GetXaxis()->GetBinWidth(1)/2.);
  Int_t qbinH = qbinL+1;
  if(qbinL <= 0) return 1.0;
  if(qbinH > fFSIss[fFSIindex]->GetNbinsX()) {// Scaled Gamov approximation 
    int chargeproduct = 1;
    if(charge1!=charge2) {
      chargeproduct = -1;
      Float_t ScaleFac = (fFSIos[fFSIindex]->GetBinContent(fFSIos[fFSIindex]->GetNbinsX()-1) - 1);
      ScaleFac /= (Gamov(chargeproduct, fFSIos[fFSIindex]->GetXaxis()->GetBinCenter(fFSIos[fFSIindex]->GetNbinsX()-1)) - 1);
      return ( (Gamov(chargeproduct, qinv)-1)*ScaleFac + 1); 
    }else{
      Float_t ScaleFac = (fFSIss[fFSIindex]->GetBinContent(fFSIss[fFSIindex]->GetNbinsX()-1) - 1);
      ScaleFac /= (Gamov(chargeproduct, fFSIss[fFSIindex]->GetXaxis()->GetBinCenter(fFSIss[fFSIindex]->GetNbinsX()-1)) - 1);
      return ( (Gamov(chargeproduct, qinv)-1)*ScaleFac + 1);
    }
  }
  
  Float_t slope=0;
  if(charge1==charge2){
    slope = fFSIss[fFSIindex]->GetBinContent(qbinL) - fFSIss[fFSIindex]->GetBinContent(qbinH);
    slope /= fFSIss[fFSIindex]->GetXaxis()->GetBinCenter(qbinL) - fFSIss[fFSIindex]->GetXaxis()->GetBinCenter(qbinH);
    return (slope*(qinv - fFSIss[fFSIindex]->GetXaxis()->GetBinCenter(qbinL)) + fFSIss[fFSIindex]->GetBinContent(qbinL));
  }else {
    slope = fFSIos[fFSIindex]->GetBinContent(qbinL) - fFSIos[fFSIindex]->GetBinContent(qbinH);
    slope /= fFSIos[fFSIindex]->GetXaxis()->GetBinCenter(qbinL) - fFSIos[fFSIindex]->GetXaxis()->GetBinCenter(qbinH);
    return (slope*(qinv - fFSIos[fFSIindex]->GetXaxis()->GetBinCenter(qbinL)) + fFSIos[fFSIindex]->GetBinContent(qbinL));
  }
}
//__________________________________________________________________________
void fcn_c3(int& npar, double* deriv, double& f, double par[], int flag){

  double q12=0, q13=0, q23=0;
  double Expan12=0, Expan13=0, Expan23=0;
  double C=0;
  double Rch=par[2]/FmToGeV;
  double SumChi2=0;
  //double lnL=0, term1=0, term2=0;
  NFitPoints_c3=0;
  //double SumChi2_test=0;

  for(int i=1; i<=BINLIMIT_3; i++){// q12
    for(int j=1; j<=BINLIMIT_3; j++){// q13
      for(int k=1; k<=BINLIMIT_3; k++){// q23
	
	if(B_3[i][j][k] == 0) continue;
	if(A_3[i][j][k] == 0) continue;
	if(A_3_e[i][j][k] == 0) continue;

	q12 = BinCenters[i];
	q13 = BinCenters[j];
	q23 = BinCenters[k];
	double q3 = sqrt(pow(q12,2)+pow(q13,2)+pow(q23,2));
	if(q3 > Q3LimitHigh) continue;
	if(q3 < Q3LimitLow) continue;
	//
	double arg12 = q12*Rch;
	double arg13 = q13*Rch;
	double arg23 = q23*Rch;
	if(FitType==0){// Edgeworth expansion
	  Expan12 = 1;
	  Expan12 += par[3]/pow(2.,3/2.)/(6.)*(8*pow(arg12,3) - 12*pow(arg12,1));
	  Expan12 += par[4]/pow(2.,4/2.)/(24.)*(16*pow(arg12,4) -48*pow(arg12,2) + 12);
	  Expan12 += par[5]/pow(2.,5/2.)/(120.)*(32.*pow(arg12,5) - 160.*pow(arg12,3) + 120*arg12);
	  Expan12 += par[6]/pow(2.,6/2.)/(720.)*(64*pow(arg12,6) - 480*pow(arg12,4) + 720*pow(arg12,2) - 120);
	  //
	  Expan13 = 1;
	  Expan13 += par[3]/pow(2.,3/2.)/(6.)*(8*pow(arg13,3) - 12*pow(arg13,1));
	  Expan13 += par[4]/pow(2.,4/2.)/(24.)*(16*pow(arg13,4) -48*pow(arg13,2) + 12);
	  Expan13 += par[5]/pow(2.,5/2.)/(120.)*(32.*pow(arg13,5) - 160.*pow(arg13,3) + 120*arg13);
	  Expan13 += par[6]/pow(2.,6/2.)/(720.)*(64*pow(arg13,6) - 480*pow(arg13,4) + 720*pow(arg13,2) - 120);
	  //
	  Expan23 = 1;
	  Expan23 += par[3]/pow(2.,3/2.)/(6.)*(8*pow(arg23,3) - 12*pow(arg23,1));
	  Expan23 += par[4]/pow(2.,4/2.)/(24.)*(16*pow(arg23,4) -48*pow(arg23,2) + 12);
	  Expan23 += par[5]/pow(2.,5/2.)/(120.)*(32.*pow(arg23,5) - 160.*pow(arg23,3) + 120*arg23);
	  Expan23 += par[6]/pow(2.,6/2.)/(720.)*(64*pow(arg23,6) - 480*pow(arg23,4) + 720*pow(arg23,2) - 120);
	}else if(FitType==1){// Laguerre expansion
	  Expan12 = 1;
	  Expan12 += par[3]*(-arg12 + 1);
	  Expan12 += par[4]/2.*(pow(arg12,2) - 4*arg12 + 2);
	  Expan12 += par[5]/6.*(-pow(arg12,3) + 9*pow(arg12,2) - 18*arg12 + 6);
	  Expan12 += par[6]/24.*(pow(arg12,4) - 16*pow(arg12,3) + 72*pow(arg12,2) - 96*arg12 + 24);
	  //
	  Expan13 = 1;
	  Expan13 += par[3]*(-arg13 + 1);
	  Expan13 += par[4]/2.*(pow(arg13,2) - 4*arg13 + 2);
	  Expan13 += par[5]/6.*(-pow(arg13,3) + 9*pow(arg13,2) - 18*arg13 + 6);
	  Expan13 += par[6]/24.*(pow(arg13,4) - 16*pow(arg13,3) + 72*pow(arg13,2) - 96*arg13 + 24);
	  //
	  Expan23 = 1;
	  Expan23 += par[3]*(-arg23 + 1);
	  Expan23 += par[4]/2.*(pow(arg23,2) - 4*arg23 + 2);
	  Expan23 += par[5]/6.*(-pow(arg23,3) + 9*pow(arg23,2) - 18*arg23 + 6);
	  Expan23 += par[6]/24.*(pow(arg23,4) - 16*pow(arg23,3) + 72*pow(arg23,2) - 96*arg23 + 24);
	}else{Expan12=1.0; Expan13=1.0; Expan23=1.0;}
	//
	double t12_coh = exp(-pow(Rcoh_def/FmToGeV * q12,2)/2.);
	double t23_coh = exp(-pow(Rcoh_def/FmToGeV * q23,2)/2.);
	double t13_coh = exp(-pow(Rcoh_def/FmToGeV * q13,2)/2.);
	if(Rcoh_def>=100){
	  t12_coh = exp(-pow(arg12,par[7])/2.) * Expan12;
	  t23_coh = exp(-pow(arg23,par[7])/2.) * Expan23;
	  t13_coh = exp(-pow(arg13,par[7])/2.) * Expan13;
	}
	C = 2*pow(par[1],3) * pow(1-G_def,3)*exp(-(pow(arg12,par[7])+pow(arg13,par[7])+pow(arg23,par[7]))/2.)*Expan12*Expan13*Expan23;
	C += 2*pow(par[1],2) * G_def*pow(1-G_def,2)*exp(-(pow(arg12,par[7])+pow(arg13,par[7]))/2.)*Expan12*Expan13*t23_coh;
	C += 2*pow(par[1],2) * G_def*pow(1-G_def,2)*exp(-(pow(arg12,par[7])+pow(arg23,par[7]))/2.)*Expan12*Expan23*t13_coh;
	C += 2*pow(par[1],2) * G_def*pow(1-G_def,2)*exp(-(pow(arg13,par[7])+pow(arg23,par[7]))/2.)*Expan13*Expan23*t12_coh;
	if(!CumulantFit){
	  C += pow(par[1],2) * pow(1-G_def,2)*exp(-pow(arg12,par[7]))*pow(Expan12,2);
	  C += pow(par[1],2) * pow(1-G_def,2)*exp(-pow(arg13,par[7]))*pow(Expan13,2);
	  C += pow(par[1],2) * pow(1-G_def,2)*exp(-pow(arg23,par[7]))*pow(Expan23,2);
	  C += 2*par[1] * G_def*(1-G_def)*exp(-pow(arg12,par[7])/2.)*Expan12 * t12_coh;
	  C += 2*par[1] * G_def*(1-G_def)*exp(-pow(arg13,par[7])/2.)*Expan13 * t13_coh;
	  C += 2*par[1] * G_def*(1-G_def)*exp(-pow(arg23,par[7])/2.)*Expan23 * t23_coh;
	}
	C += 1.0;
	C *= par[0];// Norm
	

	double error = 0;
	if(!CumulantFit){
	  error = pow(A_3_e[i][j][k]/B_3[i][j][k],2);
	  error += pow(sqrt(B_3[i][j][k])*A_3[i][j][k]/pow(B_3[i][j][k],2),2);
	  error = sqrt(error);
	  SumChi2 += pow( (C - A_3[i][j][k]/B_3[i][j][k])/error,2);
	}else{
	  error = pow(a_3_e[i][j][k]/b_3[i][j][k],2);
	  error += pow(sqrt(b_3[i][j][k])*a_3[i][j][k]/pow(b_3[i][j][k],2),2);
	  error = sqrt(error);
	  SumChi2 += pow( (C - a_3[i][j][k]/b_3[i][j][k])/error,2);
	}
	//if(q3<0.05) SumChi2_test += pow( (C - A_3[i][j][k]/B_3[i][j][k])/error,2);
	//
	/*if(A_3[i][j][k] >= 1) term1 = C*(A_3[i][j][k]+B_3[i][j][k])/(A_3[i][j][k]*(C+1));
	else term1 = 0;
	term2 = (A_3[i][j][k]+B_3[i][j][k])/(B_3[i][j][k]*(C+1));
	
	if(term1 > 0.0 && term2 > 0.0){
	  lnL += A_3[i][j][k]*log(term1) + B_3[i][j][k]*log(term2);
	}else if(term1==0 && term2 > 0.0){
	  lnL += B_3[i][j][k]*log(term2);
	}else {cout<<"WARNING -- term1 and term2 are negative"<<endl;}
	*/

	NFitPoints_c3++;
	
      }
    }
  }
  //f = -2.0*lnL;// log-liklihood minimization
  f = SumChi2;// Chi2 minimization
  Chi2_c3 = f;
  //Chi2_c3 = SumChi2_test;
  
}
