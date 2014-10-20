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

//
int CollisionType_def=1;// PbPb, pPb, pp
int FitType=0;// (0)Edgeworth, (1)Laguerre, (2)Levy
//
int Mbin=0;// 0-9: centrality bin in widths of 5%
int CHARGE=-1;// -1 or +1: + or - pions for same-charge case, --+ or -++,  ---+ or -+++
//
int EDbin=0;// 0 or 1: Kt3 bin
double G_def = 90.;// coherent %
//
bool MRC=1;// Momentum Resolution Corrections?
bool MuonCorrection=1;// correct for Muons?
//
int f_choice=0;// 0(Core/Halo), 1(40fm), 2(70fm), 3(100fm)
//
//
//
//
bool SaveToFile_def=0;
int fFSIindex=0;
float TwoFrac;// Lambda parameter
float OneFrac;// Lambda^{1/2}
float ThreeFrac;// Lambda^{3/2}
double Qcut_pp = 0.6;// 0.6
double Qcut_PbPb = 0.1;
double NormQcutLow_pp = 1.0;
double NormQcutHigh_pp = 1.5;
double NormQcutLow_PbPb = .15;
double NormQcutHigh_PbPb = .2;// was .175


const int BINRANGE_3=60;// q12,q13,q23
int BINLIMIT_3;
double A_3[BINRANGE_3][BINRANGE_3][BINRANGE_3];// q12,q13,q23
double A_3_e[BINRANGE_3][BINRANGE_3][BINRANGE_3];// q12,q13,q23
double B_3[BINRANGE_3][BINRANGE_3][BINRANGE_3];// q12,q13,q23
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



void Fit_c3(bool SaveToFile=SaveToFile_def, int CollisionType=CollisionType_def, double G=G_def){
  
  SaveToFile_def=SaveToFile;
  CollisionType_def=CollisionType;
  G /= 100.;
  G_def=G;

  TFile *_file0;
  if(CollisionType==0){// PbPb
    _file0 = new TFile("Results/PDC_11h_3Dhistos.root","READ");
  }else if(CollisionType==1){// pPb
    _file0 = new TFile("Results/PDC_13bc_kINT7.root","READ");
  }else{// pp
    _file0 = new TFile("Results/PDC_10bcde_kMB.root","READ");
  }

  TList *MyList;
  TDirectoryFile *tdir = (TDirectoryFile*)_file0->Get("PWGCF.outputFourPionAnalysis.root");
  MyList=(TList*)tdir->Get("FourPionOutput_1");
  //MyList=(TList*)_file0->Get("MyList");

  _file0->Close();


  if(CollisionType==0) {Q3LimitLow = 0.02; Q3LimitHigh = 0.06;}// 0.01 and 0.08 
  else {Q3LimitLow = 0.08; Q3LimitHigh = 0.2;}// 0.01 and 0.25
  
  
  //
  TwoFrac=0.7;
  OneFrac = sqrt(TwoFrac);
  ThreeFrac = pow(TwoFrac, 1.5);
  
  //// Core/Halo, 40fm, 70fm, 100fm
  float ThermShift_f33[4]={0., 0.06933, 0.01637, 0.006326};
  float ThermShift_f32[4]={0., -0.0185, -0.005301, -0.002286};
  float ThermShift_f31[4]={0., -0.01382, -0.0004682, 0.0005337};
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
  
  TF1 *Unity=new TF1("Unity","1",0,100);
  Unity->SetLineStyle(2);


  int ch1=0,ch2=0,ch3=0;
  
  
  
  ///////////////////////////////////////////////////////////
  // 3-pion 
  cout<<"3-pion section"<<endl;
 
  TCanvas *can2 = new TCanvas("can2", "can2",600,53,700,500);
  can2->SetHighLightColor(2);
  can2->Range(-0.7838115,-1.033258,0.7838115,1.033258);
  gStyle->SetOptFit(0111);
  can2->SetFillColor(10);
  can2->SetBorderMode(0);
  can2->SetBorderSize(2);
  can2->SetGridx();
  can2->SetGridy();
  can2->SetFrameFillColor(0);
  can2->SetFrameBorderMode(0);
  can2->SetFrameBorderMode(0);
  gPad->SetRightMargin(0.02); gPad->SetTopMargin(0.02);
 
  int Q3BINS=12;
  float Q3HistoLimit=0.12;
  if(CollisionType!=0){ Q3BINS=60; Q3HistoLimit=0.6;}

  TH1D *c3_hist = new TH1D("c3_hist","",Q3BINS,0,Q3HistoLimit);
  TH1D *Combinatorics_1d = new TH1D("Combinatorics_1d","",Q3BINS,0,Q3HistoLimit);
  TH1D *GenSignalExpected_num=new TH1D("GenSignalExpected_num","",Q3BINS,0,Q3HistoLimit);
  TH1D *GenSignalExpected_den=new TH1D("GenSignalExpected_den","",Q3BINS,0,Q3HistoLimit);
  double c3_e[100]={0};
  
  
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
	N3_QS -= ( pow(1-OneFrac,3) + 3*OneFrac*pow(1-OneFrac,2) )*TERM5;
	N3_QS -= (1-OneFrac)*(TERM2 + TERM3 + TERM4 - 3*(1-TwoFrac)*TERM5);
	N3_QS /= ThreeFrac;
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
	
	c3_e[Q3bin-1] += value_num_e + TERM5;// add baseline stat error


	// Fill histograms
	c3_hist->Fill(Q3, value_num + TERM5);// for cumulant correlation function
	Combinatorics_1d->Fill(Q3, TERM5);
	
	//
	A_3[i-1][j-1][k-1] = value_num + TERM5;
	B_3[i-1][j-1][k-1] = TERM5;
	A_3_e[i-1][j-1][k-1] = sqrt(value_num_e + TERM5);
	//if(i==j && i==k && i==4) cout<<A_3[i-1][j-1][k-1]<<"  "<<B_3[i-1][j-1][k-1]<<"  "<<A_3_e[i-1][j-1][k-1]<<endl;
	///////////////////////////////////////////////////////////
	
      }
    }
  }
  
  
 
  ////////////////////////////

  // Intermediate steps
  for(int i=0; i<Q3BINS; i++) {c3_hist->SetBinError(i+1, sqrt(c3_e[i]));}

  c3_hist->Divide(Combinatorics_1d);

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  
  
 
  c3_hist->GetXaxis()->SetTitle("#font[12]{Q}_{3} (GeV/#font[12]{c})");
  c3_hist->GetYaxis()->SetTitle("#font[12]{#bf{c}}_{3}");
  c3_hist->GetYaxis()->SetTitleSize(0.045); c3_hist->GetXaxis()->SetTitleSize(0.045);
  c3_hist->GetYaxis()->SetTitleOffset(1.1);
  c3_hist->GetXaxis()->SetRangeUser(0,Q3LimitHigh);
  c3_hist->GetYaxis()->SetRangeUser(0.9,4);
  c3_hist->SetMarkerStyle(25);
  c3_hist->SetMarkerColor(2);
  c3_hist->SetLineColor(2);
  c3_hist->SetMaximum(3);
  c3_hist->SetMinimum(.7);
  c3_hist->Draw();
  //legend2->AddEntry(c3_hist,"#font[12]{#bf{c}}_{3} (cumulant correlation)","p");
  
  

  const int npar_c3=7;
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
  

  TF1 *c3Fit1D_Expan;
  if(FitType==0) {
    c3Fit1D_Expan=new TF1("c3Fit1D_Expan","[0]*(1+[1]*exp(-pow([2]*x/0.19733,2)/2.) * pow(1 + ([3]/(6.*pow(2.,1.5))*(8.*pow([2]*x/sqrt(3.)/0.19733,3) - 12.*pow([2]*x/sqrt(3.)/0.19733,1))) + ([4]/(24.*pow(2.,2))*(16.*pow([2]*x/sqrt(3.)/0.19733,4) -48.*pow([2]*x/sqrt(3.)/0.19733,2) + 12)) + [5]/(120.*pow(2.,2.5))*(32.*pow(x/sqrt(3.)*[2]/0.19733,5) - 160.*pow(x/sqrt(3.)*[2]/0.19733,3) + 120*x/sqrt(3.)*[2]/0.19733) ,3))",0,1);
  }else if(FitType==1){
    c3Fit1D_Expan=new TF1("c3Fit1D_Expan","[0]*(1+[1]*exp(-[2]*x/0.19733 * sqrt(3.)/2.) * pow(1 + [3]*([2]*x/sqrt(3.)/0.19733 - 1) + [4]/2*(pow([2]*x/sqrt(3.)/0.19733,2) - 4*[2]*x/sqrt(3.)/0.19733 + 2) + [5]/6.*(-pow(x/sqrt(3.)*[1]/0.19733,3) + 9*pow(x/sqrt(3.)*[1]/0.19733,2) - 18*x/sqrt(3.)*[1]/0.19733 + 6),3))",0,1);
  }else{
    c3Fit1D_Expan=new TF1("c3Fit1D_Expan","[0]*(1+[1]*exp(-pow([2]*x/0.19733, [3])))",0,1);
  }
  

  double OutputPar_c3[npar_c3]={0};
  double OutputPar_c3_e[npar_c3]={0};
  
  double par_c3[npar_c3];               // the start values
  double stepSize_c3[npar_c3];          // step sizes 
  double minVal_c3[npar_c3];            // minimum bound on parameter 
  double maxVal_c3[npar_c3];            // maximum bound on parameter
  string parName_c3[npar_c3];
  //          1.0              1.5              10.              0.              0.              0.             1.5  
  par_c3[0] = 1.0; par_c3[1] = 1.5; par_c3[2] = 10.; par_c3[3] = 0.; par_c3[4] = 0.; par_c3[5] = 0; par_c3[6] = 1.5;
  stepSize_c3[0] = 0.01; stepSize_c3[1] = 0.1; stepSize_c3[2] = 0.1; stepSize_c3[3] = 0.01; stepSize_c3[4] = 0.01; stepSize_c3[5] = 0.01; stepSize_c3[6] = 0.1;
  minVal_c3[0] = 0.8; minVal_c3[1] = 0.4; minVal_c3[2] = 4.; minVal_c3[3] = -2; minVal_c3[4] = -1; minVal_c3[5] = -1; minVal_c3[6] = 0.5;
  maxVal_c3[0] = 1.1; maxVal_c3[1] = 2000.; maxVal_c3[2] = 100.; maxVal_c3[3] = +2; maxVal_c3[4] = +1; maxVal_c3[5] = +1; maxVal_c3[6] = 2.5;
  parName_c3[0] = "N"; parName_c3[1] = "#lambda_{3}"; parName_c3[2] = "R_{inv}"; parName_c3[3] = "coeff_{3}"; parName_c3[4] = "coeff_{4}"; parName_c3[5] = "coeff_{5}"; parName_c3[6] = "#alpha";

  if(CollisionType==0){ 
    if(FitType!=0) {
      par_c3[2]=15.; minVal_c3[2] = 8.;
    }
  }else{
    if(FitType==0) {par_c3[2] = 2.0; minVal_c3[2] = 1.0;}
    else {
      par_c3[1] = 4.0; minVal_c3[1] = 1.0;
      //
      par_c3[2] = 1.3; minVal_c3[2] = .9; maxVal_c3[2] = 10.;
    }
  }
  
  if(FitType==0) {par_c3[6] = 2.;}
  if(FitType==1) {par_c3[6] = 1.;}
  if(FitType==2) {par_c3[3] = 0; par_c3[4] = 0; par_c3[5] = 0;}

  if(FitType==2) {maxVal_c3[1] = 10.;}

  for (int i=0; i<npar_c3; i++){
    MyMinuit_c3.DefineParameter(i, parName_c3[i].c_str(), par_c3[i], stepSize_c3[i], minVal_c3[i], maxVal_c3[i]);
  }
  if(FitType==0 || FitType==1) { MyMinuit_c3.FixParameter(6);}
  if(FitType==2){
    MyMinuit_c3.FixParameter(3);
    MyMinuit_c3.FixParameter(4);
    MyMinuit_c3.FixParameter(5);
  }
  MyMinuit_c3.FixParameter(0);
  //MyMinuit_c3.FixParameter(1);
  //MyMinuit_c3.FixParameter(2);
  //MyMinuit_c3.FixParameter(3);
  //MyMinuit_c3.FixParameter(4);
  MyMinuit_c3.FixParameter(5);
  
  /////////////////////////////////////////////////////////////
  // Do the minimization!
  cout<<"Start Three-d fit"<<endl;
  MyMinuit_c3.Migrad();// Minuit's best minimization algorithm
  cout<<"End Three-d fit"<<endl;
  /////////////////////////////////////////////////////////////
  MyMinuit_c3.mnexcm("SHOw PARameters", &arglist_c3, 1, ierflg_c3);
  cout<<"c3 Fit: Chi2 = "<<Chi2_c3<<"   NDF = "<<(NFitPoints_c3-MyMinuit_c3.GetNumFreePars())<<endl;
  cout<<" Chi2/NDF = "<<Chi2_c3 / (NFitPoints_c3-MyMinuit_c3.GetNumFreePars())<<endl;

  for(int i=0; i<npar_c3; i++){
    MyMinuit_c3.GetParameter(i,OutputPar_c3[i],OutputPar_c3_e[i]);
  }
  
  cout<<"Tij Norm = "<<pow(OutputPar_c3[1]/2., 1/3.)<<endl;
  
  if(FitType!=2){
    c3Fit1D_Expan->FixParameter(0,OutputPar_c3[0]);
    c3Fit1D_Expan->FixParameter(1,OutputPar_c3[1]);
    c3Fit1D_Expan->FixParameter(2,OutputPar_c3[2]);
    c3Fit1D_Expan->FixParameter(3,OutputPar_c3[3]);
    c3Fit1D_Expan->FixParameter(4,OutputPar_c3[4]);
    c3Fit1D_Expan->FixParameter(5,OutputPar_c3[5]);
  }else{// Levy
    c3Fit1D_Expan->FixParameter(0,OutputPar_c3[0]);
    c3Fit1D_Expan->FixParameter(1,OutputPar_c3[1]);
    c3Fit1D_Expan->FixParameter(2,OutputPar_c3[2]);
    c3Fit1D_Expan->FixParameter(3,OutputPar_c3[6]);
  }
  c3Fit1D_Expan->SetLineStyle(1);
  //c3Fit1D_Expan->Draw("same");
  
 
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
	  Expan12 += OutputPar_c3[5]/pow(2.,5/2.)/(120.)*(32.*pow(arg12,5) - 160.*pow(arg12,3) + 120*arg12);
	  //
	  Expan13 += OutputPar_c3[3]/pow(2.,3/2.)/(6.)*(8*pow(arg13,3) - 12*pow(arg13,1));
	  Expan13 += OutputPar_c3[4]/pow(2.,4/2.)/(24.)*(16*pow(arg13,4) -48*pow(arg13,2) + 12);
	  Expan13 += OutputPar_c3[5]/pow(2.,5/2.)/(120.)*(32.*pow(arg13,5) - 160.*pow(arg13,3) + 120*arg13);
	  //
	  Expan23 += OutputPar_c3[3]/pow(2.,3/2.)/(6.)*(8*pow(arg23,3) - 12*pow(arg23,1));
	  Expan23 += OutputPar_c3[4]/pow(2.,4/2.)/(24.)*(16*pow(arg23,4) -48*pow(arg23,2) + 12);
	  Expan23 += OutputPar_c3[5]/pow(2.,5/2.)/(120.)*(32.*pow(arg23,5) - 160.*pow(arg23,3) + 120*arg23);
	}else if(FitType==1){
	  Expan12 += OutputPar_c3[3]*(arg12 - 1);
	  Expan12 += OutputPar_c3[4]/2.*(pow(arg12,2) - 4*arg12 + 2);
	  Expan12 += OutputPar_c3[5]/6.*(-pow(arg12,3) + 9*pow(arg12,2) - 18*arg12 + 6);
	  //
	  Expan13 += OutputPar_c3[3]*(arg13 - 1);
	  Expan13 += OutputPar_c3[4]/2.*(pow(arg13,2) - 4*arg13 + 2);
	  Expan13 += OutputPar_c3[5]/6.*(-pow(arg13,3) + 9*pow(arg13,2) - 18*arg13 + 6);
	  //
	  Expan23 += OutputPar_c3[3]*(arg23 - 1);
	  Expan23 += OutputPar_c3[4]/2.*(pow(arg23,2) - 4*arg23 + 2);
	  Expan23 += OutputPar_c3[5]/6.*(-pow(arg23,3) + 9*pow(arg23,2) - 18*arg23 + 6);
	}else{}
	
	//
	double C = OutputPar_c3[1] * pow(1-G,3)*exp(-(pow(arg12,OutputPar_c3[6])+pow(arg13,OutputPar_c3[6])+pow(arg23,OutputPar_c3[6]))/2.)*Expan12*Expan13*Expan23;
	C += pow(OutputPar_c3[1], 2/3.) * G*pow(1-G,2)*exp(-(pow(arg12,OutputPar_c3[6])+pow(arg13,OutputPar_c3[6]))/2.)*Expan12*Expan13;
	C += pow(OutputPar_c3[1], 2/3.) * G*pow(1-G,2)*exp(-(pow(arg12,OutputPar_c3[6])+pow(arg23,OutputPar_c3[6]))/2.)*Expan12*Expan23;
	C += pow(OutputPar_c3[1], 2/3.) * G*pow(1-G,2)*exp(-(pow(arg13,OutputPar_c3[6])+pow(arg23,OutputPar_c3[6]))/2.)*Expan13*Expan23;
	C += 1.0;
	//if(Q3<0.026) cout<<Q3<<"  "<<C<<endl;
	C *= TERM5;
	C *= OutputPar_c3[0];
	//if(Q3<0.018) continue;
	GenSignalExpected_num->Fill(Q3, C);
	GenSignalExpected_den->Fill(Q3, TERM5);
	//if(Q3<0.02) cout<<Q3<<"  "<<TERM5<<endl;
	///////////////////////////////////////////////////////////
	
	
	
      }
    }
  }
  

  GenSignalExpected_num->Sumw2();
  GenSignalExpected_num->Divide(GenSignalExpected_den);
 
  TSpline3 *c3Fit1D_ExpanSpline = new TSpline3(GenSignalExpected_num);
  c3Fit1D_ExpanSpline->SetLineWidth(2);
  double xpoints[1000], ypoints[1000];
  bool splineOnly=kFALSE;
  for(int iii=0; iii<1000; iii++){
    xpoints[iii] = 0 + (iii+0.5)*0.001;
    //ypoints[iii] = c3Fit1D_ExpanSpline->Eval(xpoints[iii]);// to skip spline
    splineOnly=kTRUE;// to skip 1D approximation
    if(CollisionType==0) splineOnly=kTRUE;
    if(CollisionType!=0 && xpoints[iii] > 0.06) splineOnly=kTRUE;
    if(!splineOnly){
      ypoints[iii] = c3Fit1D_Expan->Eval(xpoints[iii]);
      if(c3Fit1D_Expan->Eval(xpoints[iii])<2. && fabs(c3Fit1D_Expan->Eval(xpoints[iii])-c3Fit1D_ExpanSpline->Eval(xpoints[iii])) < 0.01) splineOnly=kTRUE;
    }
    else {ypoints[iii] = c3Fit1D_ExpanSpline->Eval(xpoints[iii]); splineOnly=kTRUE;}
  }
  TGraph *gr_c3Spline = new TGraph(1000,xpoints,ypoints);
  gr_c3Spline->SetLineWidth(2);
  if(CollisionType==0) gr_c3Spline->SetLineColor(1);
  if(CollisionType==1) gr_c3Spline->SetLineColor(2);
  if(CollisionType==2) gr_c3Spline->SetLineColor(4);
  
  gr_c3Spline->Draw("c same");
  
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
  
  
 

  if(SaveToFile){
    TString *savefilename = new TString("FitFiles/FitFile_CT");
    *savefilename += CollisionType;
    savefilename->Append("_FT");
    *savefilename += FitType;
    savefilename->Append("_G");
    *savefilename += int((G+0.001)/0.02);
    savefilename->Append(".root");
    TFile *savefile = new TFile(savefilename->Data(),"RECREATE");
    MyMinuit_c3.Write("MyMinuit_c3");
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

  for(int i=0; i<=BINLIMIT_3; i++){// q12
    for(int j=0; j<=BINLIMIT_3; j++){// q13
      for(int k=0; k<=BINLIMIT_3; k++){// q23
	
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
	  //
	  Expan13 = 1;
	  Expan13 += par[3]/pow(2.,3/2.)/(6.)*(8*pow(arg13,3) - 12*pow(arg13,1));
	  Expan13 += par[4]/pow(2.,4/2.)/(24.)*(16*pow(arg13,4) -48*pow(arg13,2) + 12);
	  Expan13 += par[5]/pow(2.,5/2.)/(120.)*(32.*pow(arg13,5) - 160.*pow(arg13,3) + 120*arg13);
	  //
	  Expan23 = 1;
	  Expan23 += par[3]/pow(2.,3/2.)/(6.)*(8*pow(arg23,3) - 12*pow(arg23,1));
	  Expan23 += par[4]/pow(2.,4/2.)/(24.)*(16*pow(arg23,4) -48*pow(arg23,2) + 12);
	  Expan23 += par[5]/pow(2.,5/2.)/(120.)*(32.*pow(arg23,5) - 160.*pow(arg23,3) + 120*arg23);
	}else if(FitType==1){// Laguerre expansion
	  Expan12 = 1;
	  Expan12 += par[3]*(arg12 - 1);
	  Expan12 += par[4]/2.*(pow(arg12,2) - 4*arg12 + 2);
	  Expan12 += par[5]/6.*(-pow(arg12,3) + 9*pow(arg12,2) - 18*arg12 + 6);
	  //
	  Expan13 = 1;
	  Expan13 += par[3]*(arg13 - 1);
	  Expan13 += par[4]/2.*(pow(arg13,2) - 4*arg13 + 2);
	  Expan13 += par[5]/6.*(-pow(arg13,3) + 9*pow(arg13,2) - 18*arg13 + 6);
	  //
	  Expan23 = 1;
	  Expan23 += par[3]*(arg23 - 1);
	  Expan23 += par[4]/2.*(pow(arg23,2) - 4*arg23 + 2);
	  Expan23 += par[5]/6.*(-pow(arg23,3) + 9*pow(arg23,2) - 18*arg23 + 6);
	}else{Expan12=1.0; Expan13=1.0; Expan23=1.0;}
	//

	C = par[1] * pow(1-G_def,3)*exp(-(pow(arg12,par[6])+pow(arg13,par[6])+pow(arg23,par[6]))/2.)*Expan12*Expan13*Expan23;
	C += pow(par[1],2/3.) * G_def*pow(1-G_def,2)*exp(-(pow(arg12,par[6])+pow(arg13,par[6]))/2.)*Expan12*Expan13;
	C += pow(par[1],2/3.) * G_def*pow(1-G_def,2)*exp(-(pow(arg12,par[6])+pow(arg23,par[6]))/2.)*Expan12*Expan23;
	C += pow(par[1],2/3.) * G_def*pow(1-G_def,2)*exp(-(pow(arg13,par[6])+pow(arg23,par[6]))/2.)*Expan13*Expan23;
	C += 1.0;
	C *= par[0];// Norm
	

	double error = pow(A_3_e[i][j][k]/B_3[i][j][k],2);
	error += pow(sqrt(B_3[i][j][k])*A_3[i][j][k]/pow(B_3[i][j][k],2),2);
	error = sqrt(error);
	SumChi2 += pow( (C - A_3[i][j][k]/B_3[i][j][k])/error,2);
	
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
