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

#define BohrR 1963.6885 // Bohr Radius for pions
#define FmToGeV 0.19733 // conversion to fm
#define PI 3.1415926
#define masspiC 0.1395702 // pi+ mass (GeV/c^2)

using namespace std;

int FileSetting=0;// 0(standard), 1(r3 Lambda), 2(TTC), 3(PID), 4(Pos B), 5(Neg B), 6(GaussR8to5), 7(GaussR11to6), 8(4 kt bins old TTC cuts), 9 (FB5and7overlap), 10 (muon Variation: -0.5<NsigmaPion<2.0)
bool MCcase_def=kFALSE;// MC data?
int CHARGE_def=-1;// -1 or +1: + or - pions for same-charge case, --+ or ++- for mixed-charge case
bool SameCharge_def=1;// 3-pion same-charge?
int EDbin_def=0;// 0 or 1: Kt3 bin
//
bool MuonCorrection=1;// correct for Muons?
bool IncludeG_def=kFALSE;// Include coherence?
bool IncludeEW_def=kTRUE;// Include EdgeWorth coefficients?
bool IncludeEWfromTherm_def=kFALSE;// Include EdgeWorth coefficients from Therminator?
bool GRS_def=kTRUE;// Generalized RiverSide 3-body FSI correction?
int Mbin_def=0;// 0-9: centrality bin in widths of 5%
int Ktbin_def=10;// 1(0.2-0.3),..., 6(0.7-0.8), 10 = Full Range
double MRCShift=1.0;// 1.0=full Momentum Resolution Correction. 1.1 for 10% systematic increase
double KShift=1.0;// K factor shift for testing (1.0 for default).  TURN OFF STRONGSC FOR THIS TESTING
//
//
bool IncludeMJcorrection=kTRUE;// linear Mini-Jet correction for denominator of r3?
bool IncludeStrongSC=kTRUE;// same-charge strong FSI
bool SaveToFile_def=kFALSE;// Save outputs to file?
int SourceType=0;// 0=Therminator, 1=Therminator with 50fm cut (keep at 0 for default), 2=Gaussian
bool ConstantFSI=kFALSE;// Constant FSI's for each kt bin?
bool GofP=kFALSE;// Include momentum dependence of coherent fraction?
bool ChargeConstraint=kFALSE;// Include Charge Constraint for coherent states?
bool LinkN=kFALSE;// Make N(++)=N(+-)?  
bool GeneratedSignal=kFALSE;// Was the QS+FSI signal generated in MC? 
bool Tricubic=kFALSE;// Tricubic or Trilinear interpolation?  Tricubic was found to be unstable.
bool IncludeSS=kTRUE;// Include same-charge two-pion correlations in the fit?
bool IncludeOS=kTRUE;// Include mixed-charge two-pion correlations in the fit?
//
//
//
//
const int Sbins_2=1;// 2-particle Species bins. 1=Pion-Pion only. max=6 (pi-pi, pi-k, pi-p, k-p, k-k, p-p)
const int Sbins_3=1;// 3-particle Species bins. 1=Pion-Pion-Pion only. max=10
const int fitlimitLowBin = 3;
const int fitlimitHighBin = 14;// max = 20 (100MeV)
bool LHC11h=kTRUE;// always kTRUE
bool ExplicitLoop=kFALSE;// always kFALSE
bool PairCut=kTRUE;// always kTRUE
bool PbPbcase=kTRUE;// always kTRUE
int KtbinFSI;
int Ktlowbin;
int Kthighbin;
int bBin;// impact parameter bin
int MomResCentBin;// momentum resolution centrality bin
int RValue;// Radius setting for Gaussian source profile
int bValue;// impact parameter for Therminator source profile (default)
float TwoFrac;// Lambda parameter
float OneFrac;// Lambda^{1/2}
float ThreeFrac;// Lambda^{3/2}
float TwoFracMomRes;// Lambda parameter for momentum resolution corrections
float RValueMomRes;// Gaussian radius value for momentum resolution corrections
float TwoFracr3;// Lambda parameter for r3
int TPN_bin;// Two Particle Normalization bin (always 0)
double Qcut_pp = 0.6;// 0.6
double Qcut_PbPb = 0.1;
double NormQcutLow_pp = 1.0;
double NormQcutHigh_pp = 1.5;
double NormQcutLow_PbPb = .15;//1.06
double NormQcutHigh_PbPb = .175;//1.1

// single-pion emitter size (for fits with momentum dependence of coherence)
double RT = 0.72/FmToGeV;
double Temp = 0.135;// 0.135 GeV


const int Nlines = 50;
TH3D *CoulCorrOmega0SS;
TH3D *CoulCorrOmega0OS;
TH1D *CoulCorr2SS;
TH1D *CoulCorr2OS;
TF1 *StrongSC;// same-charge pion strong FSI


void ReadCoulCorrections(int, int, int);
void ReadCoulCorrections_Omega0();
void ReadMuonCorrections(int);
void ReadMomResFile(int, double);
double CoulCorr2(int, double);
double CoulCorrOmega0(bool, double, double, double);
double CoulCorrGRS(bool, double, double, double);
double OSLfitfunction(Double_t *x, Double_t *par);
double D3fitfunction_3(Double_t *x, Double_t *par);
double r3_fitfunction(Double_t *x, Double_t *par);
double Gamov(int, double);
//double LednickyCorr(double);
double C2ssFitFunction(Double_t *x, Double_t *par);
double C2osFitFunction(Double_t *x, Double_t *par);
double cubicInterpolate(double[4], double);
double bicubicInterpolate(double[4][4], double, double);
double tricubicInterpolate(double[4][4][4], double, double, double);
void DrawALICELogo(Bool_t, Float_t, Float_t, Float_t, Float_t);
//
void fcnC2_Global(int&, double*, double&, double[], int);
void fcnOSL(int&, double*, double&, double[], int);
float fact(float);

const int BINRANGE_C2global=20;
double C2ssFitting[BINRANGE_C2global];
double C2osFitting[BINRANGE_C2global];
double C2ssFitting_e[BINRANGE_C2global];
double C2osFitting_e[BINRANGE_C2global];
double K2SS[BINRANGE_C2global];
double K2OS[BINRANGE_C2global];
double K2SS_e[BINRANGE_C2global];
double K2OS_e[BINRANGE_C2global];
double K3SS_e[BINRANGE_C2global][BINRANGE_C2global][BINRANGE_C2global];
double K3OS_e[BINRANGE_C2global][BINRANGE_C2global][BINRANGE_C2global];
double Chi2_C2global;
double NFitPoints_C2global;
double Chi2_C3global;
double NFitPoints_C3global;

const int BINS_OSL=40;// out,side,long bins
double A[BINS_OSL][BINS_OSL][BINS_OSL];// out,side,long
double B[BINS_OSL][BINS_OSL][BINS_OSL];// out,side,long
double A_e[BINS_OSL][BINS_OSL][BINS_OSL];// out,side,long
double B_e[BINS_OSL][BINS_OSL][BINS_OSL];// out,side,long
double avg_q[BINS_OSL][BINS_OSL][BINS_OSL];// out,side,long
double K_OSL[BINS_OSL][BINS_OSL][BINS_OSL];// out,side,long
double K_OSL_e[BINS_OSL][BINS_OSL][BINS_OSL];// out,side,long
double Chi2_OSL;
double NFitPoints_OSL;

const int BINRANGE_3=50;// q12,q13,q23
double A_3[BINRANGE_3][BINRANGE_3][BINRANGE_3];// q12,q13,q23
double B_3[BINRANGE_3][BINRANGE_3][BINRANGE_3];// q12,q13,q23

double C3_base[BINRANGE_3][BINRANGE_3][BINRANGE_3];
double N_term1[BINRANGE_3][BINRANGE_3][BINRANGE_3];
double N_term2[BINRANGE_3][BINRANGE_3][BINRANGE_3];
double N_term3[BINRANGE_3][BINRANGE_3][BINRANGE_3];
double N_term4[BINRANGE_3][BINRANGE_3][BINRANGE_3];
double N_term5[BINRANGE_3][BINRANGE_3][BINRANGE_3];

double MomRes_C2_pp[40];// Qinv bins
double MomRes_C2_mp[40];//
double MomRes_term1_pp[40];
double MomRes_term2_pp[40];
double MomRes_term1_mp[40];
double MomRes_term2_mp[40];

TH3D *MomRes3d[2][5];// SC/MC, term#
TH1D *MomRes1d[2][5];// SC/MC, term#
TH3D *MomRes3d_FVP[2][5];// sum#, term#
TH3D *AvgK3[2];// SC/MC
double AvgP[6][20];// kt bin, qinv bin
int Ktbin_GofP;
//
double MomResDev_C2_pp[40];// Qinv bins
double MomResDev_C2_mp[40];
double MomResDev_term1_pp[40];
double MomResDev_term2_pp[40];


double r3_3d_arr[20][20][20];
double r3_3d_arr_e[20][20][20];
double AvgQinvSS[30];
double AvgQinvOS[30];
//
TH1D *C2muonCorrection_SC;
TH1D *C2muonCorrection_MC;
TH1D *WeightmuonCorrection;
TH1D *C3muonCorrection;




void Plot_PDCumulants(bool SaveToFile=SaveToFile_def, bool MCcase=MCcase_def, bool IncludeEWfromTherm=IncludeEWfromTherm_def, bool SameCharge=SameCharge_def, bool IncludeG=IncludeG_def, bool IncludeEW=IncludeEW_def, bool GRS=GRS_def, int EDbin=EDbin_def, int CHARGE=CHARGE_def, int Mbin=Mbin_def, int Ktbin=Ktbin_def){
  
  EDbin_def=EDbin;
  SaveToFile_def=SaveToFile;
  MCcase_def=MCcase;
  IncludeEWfromTherm_def=IncludeEWfromTherm;
  CHARGE_def=CHARGE;
  IncludeG_def=IncludeG;
  IncludeEW_def=IncludeEW;
  SameCharge_def=SameCharge;// 3-pion same-charge
  Mbin_def=Mbin;
  Ktbin_def=Ktbin;
  if(Ktbin_def==10) Ktbin_GofP=2;
  else Ktbin_GofP=Ktbin_def;
  //
  Ktlowbin=(Ktbin)*2+3;// kt bins are 0.5 GeV/c wide (0-0.5, 0.5-1.0 ...)
  Kthighbin=(Ktbin)*2+4;
  
  //
  if(ConstantFSI) KtbinFSI=10;
  else KtbinFSI=Ktbin;
  

  TwoFracr3 = 0.7;// standard lambda is 0.7 for r3
  if(FileSetting==1) TwoFracr3 = 0.6;

  if(Mbin == 0) {// 0-5%
    RValue=8;
    TwoFrac=0.68;
    TwoFracMomRes=0.68;
    RValueMomRes=10;// for C2 
    bValue=2; 
    bBin=1;
    MomResCentBin=1;// for C3
  }else if(Mbin == 1) {// 5-10%
    RValue=8;
    TwoFrac=0.68;
    TwoFracMomRes=0.68;
    RValueMomRes=9;
    bValue=3;
    bBin=2;
    MomResCentBin=2;
  }else if(Mbin <= 3) {// 10-20%
    RValue=7;
    TwoFrac=0.68;
    TwoFracMomRes=0.68;
    RValueMomRes=8;
    bValue=5;
    bBin=3;
    MomResCentBin=3;
  }else if(Mbin <= 5) {// 20-30%
    RValue=6;
    TwoFrac=0.68;
    TwoFracMomRes=0.68;
    RValueMomRes=7;
    bValue=7;
    bBin=4;
    MomResCentBin=4;
  }else if(Mbin <= 7){// 30-40%
    RValue=6;
    TwoFrac=0.68;
    TwoFracMomRes=0.68;
    RValueMomRes=6;
    bValue=8;
    bBin=5;
    MomResCentBin=5;
  }else {// 40-50%
    RValue=5;
    TwoFrac=0.68;
    TwoFracMomRes=0.68;
    RValueMomRes=5;
    bValue=9;
    bBin=6;
    MomResCentBin=6;
  }

  OneFrac = sqrt(TwoFrac);
  ThreeFrac = pow(TwoFrac,3/2.);

  // same-charge pion strong FSI (obtained from ratio of CoulStrong to Coul at R=7 Gaussian from Lednicky's code)
  StrongSC=new TF1("StrongSC","[0]+[1]*exp(-pow([2]*x,2))",0,50);
  StrongSC->FixParameter(0,9.99192e-01);
  StrongSC->FixParameter(1,-8.01113e-03);
  StrongSC->FixParameter(2,5.35912e-02);

  
  if(SourceType==2 && RValue > 8) {cout<<"Radius value too large!!!"<<endl; return;}

  cout<<"Mbin = "<<Mbin<<"   Kt = "<<Ktbin<<"   SameCharge = "<<SameCharge<<endl;
  
  // Core-Halo modeling of single-particle and triple-particle core fraction 
  //float AvgN[10]={472.56, 390.824, 319.597, 265.933, 218.308, 178.865, 141.2, 115.88, 87.5556, 69.3235};// 10h (avg total pion mult)/2. 2.0sigma
  //float AvgN[10]={571.2, 472.7, 400.2, 325.2, 268.721, 220.3, 178.9, 143.4, 113.412, 88.0};// 11h (avg total pion mult)/2.  2.0sigma
 
  //
  // avg Qinv within the 5 MeV bins (0-5, 5-10,...) for true bin mean values.  From unsmeared HIJING 0-5% with input signal
  double AvgQinvSS_temp[30]={0.00421164, 0.00796583, 0.0128473, 0.0178262, 0.0228416, 0.0276507, 0.0325368, 0.0376114, 0.0424707, 0.0476274, 0.0526015, 0.0575645, 0.0625478, 0.0675416, 0.0725359, 0.0775219, 0.0825521, 0.0875211, 0.0925303, 0.0975319, 0.102544, 0.10753, 0.112499, 0.117545, 0.122522, 0.127522, 0.132499, 0.137514, 0.142495, 0.147521};
  double AvgQinvOS_temp[30]={0.00352789, 0.00797596, 0.0128895, 0.0177198, 0.0226397, 0.0276331, 0.0326309, 0.0376471, 0.0426217, 0.047567, 0.0525659, 0.0575472, 0.0625886, 0.0675261, 0.0725543, 0.077564, 0.0825617, 0.0875465, 0.092539, 0.0975348, 0.102529, 0.107527, 0.112506, 0.117531, 0.122536, 0.1275, 0.132508, 0.137508, 0.14251, 0.147511};
  for(int i=0; i<30; i++){
    AvgQinvSS[i] = AvgQinvSS_temp[i];
    AvgQinvOS[i] = AvgQinvOS_temp[i];
  }
 

  // average total momentum in each kt and qinv bin.  Used for a test of momentum dependence of fits including coherence
  double AvgP_kt1[20]={0.25, 0.254421, 0.270286, 0.27346, 0.274017, 0.274252, 0.274438, 0.274675, 0.27497, 0.275354, 0.275844, 0.276433, 0.27708, 0.277768, 0.278469, 0.279179, 0.279871, 0.280536, 0.281167, 0.281749};
  double AvgP_kt2[20]={0.34, 0.34282, 0.370731, 0.38009, 0.381977, 0.382246, 0.382231, 0.38212, 0.381979, 0.381813, 0.381653, 0.381502, 0.381358, 0.381227, 0.381111, 0.381008, 0.380923, 0.380864, 0.380835, 0.380824};
  double AvgP_kt3[20]={0.47, 0.47, 0.469358, 0.486368, 0.492106, 0.49335, 0.493497, 0.49341, 0.493263, 0.493048, 0.492821, 0.492587, 0.49234, 0.492093, 0.491842, 0.491607, 0.491363, 0.491126, 0.490887, 0.490645};
  double AvgP_kt4[20]={0.59, 0.59, 0.59, 0.588732, 0.598683, 0.601864, 0.602459, 0.602302, 0.602005, 0.601617, 0.601192, 0.600794, 0.600365, 0.599961, 0.599571, 0.599208, 0.598841, 0.598511, 0.59816, 0.597836};
  double AvgP_kt5[20]={0.67, 0.67, 0.67, 0.676098, 0.694986, 0.701248, 0.703524, 0.704307, 0.704402, 0.704432, 0.704255, 0.703942, 0.703521, 0.703031, 0.702484, 0.701869, 0.701206, 0.700521, 0.699822, 0.699084};
  double AvgP_kt6[20]={0.79, 0.79, 0.79, 0.79, 0.790909, 0.799252, 0.801984, 0.801945, 0.800688, 0.799177, 0.797658, 0.796156, 0.794652, 0.793256, 0.791942, 0.790697, 0.789588, 0.788608, 0.787698, 0.786941};
  for(int ktbin=1; ktbin<=6; ktbin++){
    for(int qbin=1; qbin<=20; qbin++){
      if(ktbin==1) AvgP[ktbin-1][qbin-1] = AvgP_kt1[qbin-1];
      else if(ktbin==2) AvgP[ktbin-1][qbin-1] = AvgP_kt2[qbin-1];
      else if(ktbin==3) AvgP[ktbin-1][qbin-1] = AvgP_kt3[qbin-1];
      else if(ktbin==4) AvgP[ktbin-1][qbin-1] = AvgP_kt4[qbin-1];
      else if(ktbin==5) AvgP[ktbin-1][qbin-1] = AvgP_kt5[qbin-1];
      else AvgP[ktbin-1][qbin-1] = AvgP_kt6[qbin-1];
    }
  }

  
  
  gStyle->SetOptStat(0);
  gStyle->SetOptDate(0);
  gStyle->SetOptFit(0111);

  TFile *myfile;
  if(PbPbcase){// PbPb
    if(MCcase) {
      if(Mbin<=1){
	//myfile = new TFile("Results/PDC_HIJING_12a17ad_fix_genSignal_Rinv11.root","READ");
	//myfile = new TFile("Results/PDC_HIJING_12a17ad_fix_genSignal_Rinv11_Smeared.root","READ");
	//myfile = new TFile("Results/PDC_HIJING_12a17ad_fix_KtgenSignal_Rinv11.root","READ");
	myfile = new TFile("Results/PDC_pureHIJING_12a17a_fix_KtgenSignal_Rinv11.root","READ");
	//myfile = new TFile("Results/PDC_HIJING_10h8.root","READ");
      }else{
	myfile = new TFile("Results/PDC_HIJING_12a17b_myRun_L0p68R11_C2plots.root","READ");
      }
    }else{
      //myfile = new TFile("Results/PDC_11h_v1paper.root","READ");// v1 paper
      //if(FileSetting==0) myfile = new TFile("Results/PDC_11h_standard_and_Raw0p04TTC.root","READ");// Old standard (bad r3 Interp, 0.035TTC)
      //if(FileSetting==0) myfile = new TFile("Results/PDC_11h_2sigmaTTC_3sigmaTTC.root","READ");// v3 Standard 
      //
      //if(FileSetting==0) myfile = new TFile("Results/PDC_11h_2Kt3bins_FullTPC.root","READ");// FullTPC runlist
      //if(FileSetting==0) myfile = new TFile("Results/PDC_11h_HighNorm.root","READ");// High Norm variation 1.06<qinv<1.1
      if(FileSetting==0) myfile = new TFile("Results/PDC_11h_2Kt3bins.root","READ");// Standard 
      if(FileSetting==1) myfile = new TFile("Results/PDC_11h_Lam0p7_and_Lam0p6.root","READ");// Lam=0.7 and Lam=0.6 (3ktbins, 0.035TTC)
      if(FileSetting==2) myfile = new TFile("Results/PDC_11h_3sigmaTTC.root","READ");// TTC variation
      if(FileSetting==3) myfile = new TFile("Results/PDC_11h_PID1p5.root","READ");// PID variation (0.035TTC)
      if(FileSetting==4) myfile = new TFile("Results/PDC_11h_PosB.root","READ");// Positive B field for r3num (0.035TTC)
      if(FileSetting==5) myfile = new TFile("Results/PDC_11h_NegB.root","READ");// Negative B field for r3num (0.035TTC)
      if(FileSetting==6) myfile = new TFile("Results/PDC_11h_GaussR8to5.root","READ");// Gaussian
      if(FileSetting==7) myfile = new TFile("Results/PDC_11h_GaussR11to6.root","READ");// Gaussian
      if(FileSetting==8) myfile = new TFile("Results/PDC_11h_4ktbins.root","READ");// 4 kt bins (0.035TTC)
      if(FileSetting==9) myfile = new TFile("Results/PDC_11h_FB5and7overlap_l0p8.root","READ");// FB5and7overlap
      if(FileSetting==10) myfile = new TFile("Results/PDC_11h_muonVar.root","READ");
    }
  }else{// pp
    cout<<"pp not currently supported"<<endl;
    return;
  }

  ReadCoulCorrections(RValue, bValue, KtbinFSI);
  //ReadLednickyFile(RValue);
  ReadMomResFile(RValueMomRes, TwoFracMomRes);
  ReadCoulCorrections_Omega0();
  ReadMuonCorrections(Mbin);
  //
  /////////////////////////////////////////////////////////
  

  double NormQcutLow;
  double NormQcutHigh;
  if(PbPbcase) {
    NormQcutLow = NormQcutLow_PbPb;
    NormQcutHigh = NormQcutHigh_PbPb;
  }else {
    NormQcutLow = NormQcutLow_pp;
    NormQcutHigh = NormQcutHigh_pp;
  }

 
  TList *MyList;
  if(!MCcase){
    TDirectoryFile *tdir = (TDirectoryFile*)myfile->Get("PWGCF.outputChaoticityAnalysis.root");
    if(FileSetting!=1) MyList=(TList*)tdir->Get("ChaoticityOutput_1");
    else MyList=(TList*)tdir->Get("ChaoticityOutput_2");
  }else{
    MyList=(TList*)myfile->Get("MyList");
  }
  myfile->Close();

 
  TH1D *Events = (TH1D*)MyList->FindObject("fEvents2");
  //

  cout<<"#Events = "<<Events->Integral(Mbin+1,Mbin+1)<<endl;

  
  TH1D *ChiSquaredNDF = new TH1D("ChiSquaredNDF","",2,0.5,2.5);// Chi^2/NDF records

  // Explicit Loop Histos
  TH2D *Two_ex_2d[2][2][Sbins_2][2];
  TH1D *Two_ex[2][2][Sbins_2][2];
      
  // Normalizations
  double NormH_pc[5]={0};
  TH3D *PionNorm[2];
  //TH3D *PionNorm_e[2];

  double norm_ex_2[6][2]={{0}};
  
  // 3d method histograms
  TH3D *Three_3d[2][2][2][Sbins_3][5];
  // 4-vect Product Sum
  TH3D *TPNRejects;
  ///////////////////////////////////
  // Create Histograms
  
  // 2-particle
  for(int c1=0; c1<2; c1++){
    for(int c2=0; c2<2; c2++){
      for(int sc=0; sc<Sbins_2; sc++){
	for(int term=0; term<2; term++){
	  
	  TString *name2_ex = new TString("Explicit2_Charge1_");
	  *name2_ex += c1;
	  name2_ex->Append("_Charge2_");
	  *name2_ex += c2;
	  name2_ex->Append("_SC_");
	  *name2_ex += sc;
	  name2_ex->Append("_M_");
	  *name2_ex += Mbin;
	  name2_ex->Append("_ED_");
	  *name2_ex += 0;// EDbin
	  name2_ex->Append("_Term_");
	  *name2_ex += term+1;
	  
	  if(sc==0 || sc==3 || sc==5){
	    if( (c1+c2)==1 ) {if(c1!=0) continue;}// skip degenerate histogram
	  }
	  
	      
	  
	  Two_ex_2d[c1][c2][sc][term] = (TH2D*)MyList->FindObject(name2_ex->Data());
	  Two_ex_2d[c1][c2][sc][term]->Sumw2();
	  TString *proname = new TString(name2_ex->Data());
	  proname->Append("_pro");
	  
	  if(Ktbin==10) {Ktlowbin=1; Kthighbin=Two_ex_2d[c1][c2][sc][term]->GetNbinsX();}
	  Two_ex[c1][c2][sc][term] = (TH1D*)Two_ex_2d[c1][c2][sc][term]->ProjectionY(proname->Data(),Ktlowbin,Kthighbin);// bins:5-6 (kt:0.2-0.3)
	  
	  norm_ex_2[sc][term] = Two_ex[c1][c2][sc][term]->Integral(Two_ex[c1][c2][sc][term]->GetXaxis()->FindBin(NormQcutLow),Two_ex[c1][c2][sc][term]->GetXaxis()->FindBin(NormQcutHigh));
	  Two_ex[c1][c2][sc][term]->Scale(norm_ex_2[sc][0]/norm_ex_2[sc][term]);
	  
	  Two_ex[c1][c2][sc][term]->SetMarkerStyle(20);
	  Two_ex[c1][c2][sc][term]->GetXaxis()->SetTitle("q_{inv}");
	  Two_ex[c1][c2][sc][term]->GetYaxis()->SetTitle("C_{2}");
	  Two_ex[c1][c2][sc][term]->SetTitle("");
	  
	}// term
      }// sc
    
      // 3-particle
      for(int c3=0; c3<2; c3++){
	for(int sc=0; sc<Sbins_3; sc++){
	  for(int term=0; term<5; term++){
	   
	    TString *name3_ex = new TString("Explicit3_Charge1_");
	    *name3_ex += c1;
	    name3_ex->Append("_Charge2_");
	    *name3_ex += c2;
	    name3_ex->Append("_Charge3_");
	    *name3_ex += c3;
	    name3_ex->Append("_SC_");
	    *name3_ex += sc;
	    name3_ex->Append("_M_");
	    *name3_ex += Mbin;
	    name3_ex->Append("_ED_");
	    *name3_ex += EDbin;
	    name3_ex->Append("_Term_");
	    *name3_ex += term+1;
	    
	    
	    TString *name3_pc = new TString("PairCut3_Charge1_");
	    *name3_pc += c1;
	    name3_pc->Append("_Charge2_");
	    *name3_pc += c2;
	    name3_pc->Append("_Charge3_");
	    *name3_pc += c3;
	    name3_pc->Append("_SC_");
	    *name3_pc += sc;
	    name3_pc->Append("_M_");
	    *name3_pc += Mbin;
	    name3_pc->Append("_ED_");
	    *name3_pc += EDbin;
	    name3_pc->Append("_Term_");
	    *name3_pc += term+1;
	    
	    ///////////////////////////////////////
	    // skip degenerate histograms
	    if(sc==0 || sc==6 || sc==9){// Identical species
	      if( (c1+c2+c3)==1) {if(c1!=0 || c2!=0 || c3!=1) continue;}
	      if( (c1+c2+c3)==2) {if(c1!=0) continue;}
	    }else if(sc!=5){
	      if( (c1+c2)==1) {if(c1!=0) continue;}
	    }else {}// do nothing for pi-k-p case
	    
	    /////////////////////////////////////////
	    
	    
	    if(ExplicitLoop) {
	      /*
		Three_ex[c1][c2][c3][sc][term] = (TH1D*)MyList->FindObject(name3_ex->Data());
		
		Three_ex[c1][c2][c3][sc][term]->SetMarkerStyle(20);
		Three_ex[c1][c2][c3][sc][term]->SetLineColor(2);
		Three_ex[c1][c2][c3][sc][term]->SetMarkerColor(2);
		Three_ex[c1][c2][c3][sc][term]->GetXaxis()->SetTitle("Q_{3}");
		Three_ex[c1][c2][c3][sc][term]->GetYaxis()->SetTitle("C_{3}");
		Three_ex[c1][c2][c3][sc][term]->SetTitle("");
	      
		TString *name = new TString(name3_ex->Data());
		name->Append("_Norm");
		NormH_ex[term] = ((TH1D*)MyList->FindObject(name->Data()))->Integral();
		  		  
		if(NormH_ex[term] > 0){
		  Three_ex[c1][c2][c3][sc][term]->Scale(NormH_ex[0]/NormH_ex[term]);
		}else{ cout<<"normalization = 0.  Skipping this SC."<<endl;}
	      */		
	    }
	    if(PairCut){
	      
	      TString *nameNorm=new TString("PairCut3_Charge1_");
	      *nameNorm += c1;
	      nameNorm->Append("_Charge2_");
	      *nameNorm += c2;
	      nameNorm->Append("_Charge3_");
	      *nameNorm += c3;
	      nameNorm->Append("_SC_");
	      *nameNorm += sc;
	      nameNorm->Append("_M_");
	      *nameNorm += Mbin;
	      nameNorm->Append("_ED_");
	      *nameNorm += 0;// EDbin
	      nameNorm->Append("_Term_");
	      *nameNorm += term+1;
	      nameNorm->Append("_Norm");
	      NormH_pc[term] = ((TH1D*)MyList->FindObject(nameNorm->Data()))->Integral();
	      	      
	      if(NormH_pc[term] > 0){
		
		if(sc<=2) {
		 
		  TString *name_3d = new TString(name3_pc->Data());
		  name_3d->Append("_3D");
		  Three_3d[c1][c2][c3][sc][term] = (TH3D*)MyList->FindObject(name_3d->Data());
		  Three_3d[c1][c2][c3][sc][term]->Sumw2();
		  Three_3d[c1][c2][c3][sc][term]->Scale(NormH_pc[0]/NormH_pc[term]);
		  //cout<<"Scale factor = "<<NormH_pc[0]/NormH_pc[term]<<endl;

		  // 4-vect Product Sum
		  if(c1==c2 && c1==c3 && sc==0 && term==4){
		    
		    TString *nameDenType=new TString("PairCut3_Charge1_");
		    *nameDenType += c1;
		    nameDenType->Append("_Charge2_");
		    *nameDenType += c2;
		    nameDenType->Append("_Charge3_");
		    *nameDenType += c3;
		    nameDenType->Append("_SC_");
		    *nameDenType += sc;
		    nameDenType->Append("_M_");
		    *nameDenType += Mbin;
		    nameDenType->Append("_ED_");
		    *nameDenType += EDbin;
		    nameDenType->Append("_TPN_");
		    *nameDenType += 0;
		    
		    PionNorm[c1] = (TH3D*)MyList->FindObject(nameDenType->Data());
		    PionNorm[c1]->Sumw2();
		    PionNorm[c1]->Scale(NormH_pc[0]/NormH_pc[term]);
		    if(c1==0) {
		      TPNRejects = (TH3D*)MyList->FindObject("fTPNRejects4");
		      TPNRejects->Scale(NormH_pc[0]/NormH_pc[term]);
		    }
		    		    
		  } 
		}
		
	      }else{
		cout<<"normalization = 0.  Skipping this SC.  SC="<<sc<<"  c1="<<c1<<"  c2="<<c2<<"  c3="<<c3<<endl;
	      }	      
	      
	      
	     
	    }// pair cut
	  }// term
	  
	}// SC
	
	
      }// c3
    }// c2
  }// c1
    

  cout<<"Done getting Histograms"<<endl;
  
  
  TCanvas *can1 = new TCanvas("can1", "can1",11,53,700,500);
  can1->SetHighLightColor(2);
  can1->Range(-0.7838115,-1.033258,0.7838115,1.033258);
  gStyle->SetOptFit(0111);
  can1->SetFillColor(10);
  can1->SetBorderMode(0);
  can1->SetBorderSize(2);
  can1->SetGridx();
  can1->SetGridy();
  can1->SetFrameFillColor(0);
  can1->SetFrameBorderMode(0);
  can1->SetFrameBorderMode(0);
  
  TLegend *legend1 = new TLegend(.6,.67,.87,.87,NULL,"brNDC");
  legend1->SetBorderSize(1);
  legend1->SetTextSize(.04);// small .03; large .036 
  legend1->SetFillColor(0);
  
  /////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////
  // This part for plotting track splitting/merging effects in MC data only
  /*
  TH3F *Merge3d_num=(TH3F*)MyList->FindObject("fPairsDetaDPhiNum");
  TH3F *Merge3d_den=(TH3F*)MyList->FindObject("fPairsDetaDPhiDen");
  //TH3F *Merge3d_num=(TH3F*)MyList->FindObject("Pairs_dEtadPhi_UNI_num");
  //TH3F *Merge3d_den=(TH3F*)MyList->FindObject("Pairs_dEtadPhi_UNI_den");
 
  TH1D *Merge1d_num[10];
  TH1D *Merge1d_den[10];
  TString *newnamenum[10];
  TString *newnameden[10];
  TF1 *MergedGaus=new TF1("MergedGaus","1-[0]*exp(-pow(x/[1],2))",-0.1, 0.1);
  MergedGaus->SetParName(0,"Amplitude");
  MergedGaus->SetParName(1,"width");
  MergedGaus->SetParameter(0,0.06);
  MergedGaus->SetParameter(1,0.01);
  MergedGaus->SetParLimits(0,0.001,0.5);
  MergedGaus->SetParLimits(1,0.001,0.1);
  
  for(int i=2; i<10; i++){
    if(i!=5 && i!=8) continue;// 5 and 8
    newnamenum[i]=new TString("namenum_");
    *newnamenum[i] += i;
    newnameden[i]=new TString("nameden_");
    *newnameden[i] += i;
  
    Merge1d_num[i]=(TH1D*)Merge3d_num->ProjectionZ(newnamenum[i]->Data(),i+1,i+1,90,110);//90,110 (phi)
    Merge1d_den[i]=(TH1D*)Merge3d_den->ProjectionZ(newnameden[i]->Data(),i+1,i+1,90,110);// (phi)
    //Merge1d_num[i]=(TH1D*)Merge3d_num->ProjectionY(newnamenum[i]->Data(),i+1,i+1,190,410);// 190,410 (eta)
    //Merge1d_den[i]=(TH1D*)Merge3d_den->ProjectionY(newnameden[i]->Data(),i+1,i+1,190,410);// (eta)
    //Merge1d_num[i]->Rebin(2);
    //Merge1d_den[i]->Rebin(2);
    Merge1d_num[i]->Sumw2();
    Merge1d_den[i]->Sumw2();
    Merge1d_num[i]->SetMarkerStyle(20);

    if(Merge1d_den[i]->Integral(1,100)<=0) continue;
    double SF_merge = Merge1d_num[i]->Integral(1,100)/Merge1d_den[i]->Integral(1,100);// Z projection (phi)
    //double SF_merge = Merge1d_num[i]->Integral(1,50)/Merge1d_den[i]->Integral(1,50);// Y projection (eta)
    Merge1d_den[i]->Scale(SF_merge);
    Merge1d_num[i]->Divide(Merge1d_den[i]);
   
    
    if(i<9){
      Merge1d_num[i]->SetLineColor(i+1);
      Merge1d_num[i]->SetMarkerColor(i+1);
    }else{
      Merge1d_num[i]->SetLineColor(11);
      Merge1d_num[i]->SetMarkerColor(11);
    }
    if(i==4) {
      Merge1d_num[i]->SetLineColor(2);
      Merge1d_num[i]->SetMarkerColor(2);
    }
    if(i==5) {
      Merge1d_num[i]->GetXaxis()->SetTitle("#Delta#phi*");
      //Merge1d_num[i]->GetXaxis()->SetTitle("#Delta#eta");
      Merge1d_num[i]->GetYaxis()->SetTitle("C_{2}^{HIJING}");
      Merge1d_num[i]->GetXaxis()->SetRangeUser(-.1,.1);
      Merge1d_num[i]->SetMinimum(.91);
      Merge1d_num[i]->SetMaximum(1.1);
      Merge1d_num[i]->DrawCopy();
           
      //Merge1d_num[i]->Fit(MergedGaus,"IME","",-0.1,0.1);
    }else{
      Merge1d_num[i]->DrawCopy("same");
    }
    
    TString *Dname=new TString("D=0.2*");
    *Dname +=i;
    Dname->Append(" m");
    legend1->AddEntry(newnamenum[i]->Data(),Dname->Data(),"lpf");
  }
  legend1->Draw("same");
  gStyle->SetOptFit(111);
  Merge1d_num[8]->Fit(MergedGaus,"IME","",-0.1,0.1);
  MergedGaus->Draw("same");
  */
  /*TH3D *PadRowNum3= (TH3D*)MyList->FindObject("fPairsPadRowNum");// kt, shfrac, qinv
  TH3D *PadRowDen3= (TH3D*)MyList->FindObject("fPairsPadRowDen");// kt, shfrac, qinv
  PadRowDen3->Scale(PadRowNum3->Integral(1,20,1,159, 35,40)/PadRowDen3->Integral(1,20,1,159, 35,40));
  PadRowNum3->GetYaxis()->SetRangeUser(0,0.01);
  PadRowDen3->GetYaxis()->SetRangeUser(0,0.01);
  TH1D *PadRowNum=(TH1D*)PadRowNum3->Project3D("z");
  TH1D *PadRowDen=(TH1D*)PadRowDen3->Project3D("z");
  PadRowNum->Divide(PadRowDen);
  PadRowNum->Draw();*/
  /*
  TH3D *PadRowNum3= (TH3D*)MyList->FindObject("fPairsShareFracDPhiNum");// r, shfrac, deltaphi
  TH3D *PadRowDen3= (TH3D*)MyList->FindObject("fPairsShareFracDPhiDen");// r, shfrac, deltaphi
  PadRowDen3->Scale(PadRowNum3->Integral(1,10,1,159, 90,100)/PadRowDen3->Integral(1,10,1,159, 90,100));
  PadRowNum3->GetXaxis()->SetRange(5,5);
  PadRowDen3->GetXaxis()->SetRange(5,5);
  TH2D *PadRowNum=(TH2D*)PadRowNum3->Project3D("zy");
  TH2D *PadRowDen=(TH2D*)PadRowDen3->Project3D("zy");
  PadRowNum->Divide(PadRowDen);
  PadRowNum->Draw("lego");
  */
  /////////////////////////
  // 2-particle legend
  // 0 = pi-pi
  // 1 = pi-k
  // 2 = pi-p
  // 3 = k-k
  // 4 = k-p
  // 5 = p-p
  /////////////////////////
  TH1D *Two_pipi_mp = (TH1D*)Two_ex[0][1][0][0]->Clone();
  Two_pipi_mp->Divide(Two_ex[0][1][0][1]);

  // Normalization range counting.  Just to evaluate required normalization width in qinv.
  //cout<<Two_ex[0][0][0][0]->Integral(Two_ex[0][0][0][0]->GetXaxis()->FindBin(1.06), Two_ex[0][0][0][0]->GetXaxis()->FindBin(1.1))<<endl;
  //cout<<Two_ex[0][0][0][0]->Integral(Two_ex[0][0][0][0]->GetXaxis()->FindBin(0.15), Two_ex[0][0][0][0]->GetXaxis()->FindBin(0.175))<<endl;
  //cout<<Two_ex[0][0][0][0]->Integral(Two_ex[0][0][0][0]->GetXaxis()->FindBin(0.15), Two_ex[0][0][0][0]->GetXaxis()->FindBin(0.175))<<endl;


  const int SCOI_2=0;
 
  TH1D *Two_ex_clone_mm=(TH1D*)Two_ex[0][0][SCOI_2][0]->Clone();
  Two_ex_clone_mm->Divide(Two_ex[0][0][SCOI_2][1]);
  TH1D *Two_ex_clone_mp=(TH1D*)Two_ex[0][1][SCOI_2][0]->Clone();
  Two_ex_clone_mp->Divide(Two_ex[0][1][SCOI_2][1]);
  TH1D *Two_ex_clone_pp=(TH1D*)Two_ex[1][1][SCOI_2][0]->Clone();
  Two_ex_clone_pp->Divide(Two_ex[1][1][SCOI_2][1]);
  
  // Mini-jet ++ background linear estimation.
  TF1 *MJ_bkg_ss=new TF1("MJ_bkg_ss","pol1",0,1);
  Two_ex_clone_mm->Fit(MJ_bkg_ss,"IMENQ","",0.2,0.4);
  cout<<"Non-femto bkg C2(q=0.01) = "<<MJ_bkg_ss->Eval(0.01)<<endl;

  Two_ex_clone_mm->GetYaxis()->SetTitle("C_{2}");
  Two_ex_clone_mm->SetTitle("");
  Two_ex_clone_mp->GetYaxis()->SetTitle("C_{2}");
  Two_ex_clone_mm->SetMarkerColor(2);
  Two_ex_clone_mm->SetLineColor(2);
  Two_ex_clone_mp->SetMarkerColor(1);
  Two_ex_clone_pp->SetMarkerColor(4);
  Two_ex_clone_pp->SetLineColor(4);
  Two_ex_clone_mm->GetXaxis()->SetRangeUser(0,0.1);
  Two_ex_clone_mm->SetMinimum(0.95);
  Two_ex_clone_mm->SetMaximum(1.4);
  
  if(MCcase){
    Two_ex_clone_mp->SetMarkerColor(4);
    Two_ex_clone_mp->SetLineColor(4);
    Two_ex_clone_mm->Add(Two_ex_clone_pp);
    Two_ex_clone_mm->Scale(0.5);
    Two_ex_clone_mm->GetYaxis()->SetTitle("C_{2}^{HIJING}");
    Two_ex_clone_mm->GetYaxis()->SetTitleOffset(1.3);
    Two_ex_clone_mm->SetMinimum(0.97);
    Two_ex_clone_mm->SetMaximum(1.02);
    Two_ex_clone_mm->DrawCopy();
    //Two_ex_clone_pp->DrawCopy("same");
    Two_ex_clone_mp->DrawCopy("same");
    legend1->AddEntry(Two_ex_clone_mm,"same-charge","p");
    //legend1->AddEntry(Two_ex_clone_pp,"++","p");
    legend1->AddEntry(Two_ex_clone_mp,"mixed-charge","p");
    legend1->Draw("same");
  }

 
    
  ////////////////////////////////////////////
  // Get Therminator EdgeWorth coefficients
  TString *ThermName = new TString("../ThermData/Therm_FSI_b");
  *ThermName += bValue;
  ThermName->Append(".root");
  TFile *Thermfile = new TFile(ThermName->Data(),"READ");
  TH2D *thermNum2d_ss=(TH2D*)Thermfile->Get("Num_CosFSI_ss");
  TH2D *thermNumSq2d_ss=(TH2D*)Thermfile->Get("NumSq_CosFSI_ss");
  TH2D *thermDen2d_ss=(TH2D*)Thermfile->Get("Den_ss");
  TH2D *thermLargeR2D_ss=(TH2D*)Thermfile->Get("LargeRpairs_ss");
  TH1D *thermNum_ss=(TH1D*)thermNum2d_ss->ProjectionY("thermNum_ss",Ktlowbin,Kthighbin);
  TH1D *thermNumSq_ss=(TH1D*)thermNumSq2d_ss->ProjectionY("thermNumSq_ss",Ktlowbin,Kthighbin);
  TH1D *thermDen_ss=(TH1D*)thermDen2d_ss->ProjectionY("thermDen_ss",Ktlowbin,Kthighbin);
  TH1D *thermLargeR_ss=(TH1D*)thermLargeR2D_ss->ProjectionY("thermLargeR_ss",Ktlowbin,Kthighbin);
  //
  TH2D *thermNum2d_os=(TH2D*)Thermfile->Get("Num_CosFSI_os");
  TH2D *thermNumSq2d_os=(TH2D*)Thermfile->Get("NumSq_CosFSI_os");
  TH2D *thermDen2d_os=(TH2D*)Thermfile->Get("Den_os");
  TH2D *thermLargeR2D_os=(TH2D*)Thermfile->Get("LargeRpairs_os");
  TH1D *thermNum_os=(TH1D*)thermNum2d_os->ProjectionY("thermNum_os",Ktlowbin,Kthighbin);
  TH1D *thermNumSq_os=(TH1D*)thermNumSq2d_os->ProjectionY("thermNumSq_os",Ktlowbin,Kthighbin);
  TH1D *thermDen_os=(TH1D*)thermDen2d_os->ProjectionY("thermDen_os",Ktlowbin,Kthighbin);
  TH1D *thermLargeR_os=(TH1D*)thermLargeR2D_os->ProjectionY("thermLargeR_os",Ktlowbin,Kthighbin);
  TH1D *C2Therm_ss = (TH1D*)thermNum_ss->Clone();
  TH1D *C2Therm_os = (TH1D*)thermNum_os->Clone();
  TH1D *C2Den_ss = (TH1D*)thermDen_ss->Clone();
  TH1D *C2Den_os = (TH1D*)thermDen_os->Clone();
  C2Therm_ss->Add(thermLargeR_ss);
  C2Den_ss->Add(thermLargeR_ss);
  C2Therm_os->Add(thermLargeR_os);
  C2Den_os->Add(thermLargeR_os);
  C2Therm_ss->Sumw2();
  C2Therm_os->Sumw2();

  C2Therm_ss->Divide(C2Den_ss);
  C2Therm_os->Divide(C2Den_os);
 

  for(int i=1; i<=thermDen_ss->GetNbinsX(); i++){
    if(thermDen_ss->GetBinContent(i) > 0){
      double err = thermNumSq_ss->GetBinContent(i)/thermDen_ss->GetBinContent(i);
      err -= pow(thermNum_ss->GetBinContent(i)/thermDen_ss->GetBinContent(i),2);
      err /= thermDen_ss->GetBinContent(i);
      err = sqrt(fabs(err));
      C2Therm_ss->SetBinError(i, err);
    }
    if(thermDen_os->GetBinContent(i) > 0){
      double err = thermNumSq_os->GetBinContent(i)/thermDen_os->GetBinContent(i);
      err -= pow(thermNum_os->GetBinContent(i)/thermDen_os->GetBinContent(i),2);
      err /= thermDen_os->GetBinContent(i);
      err = sqrt(fabs(err));
      C2Therm_os->SetBinError(i, err);
    }
  }
  
  C2Therm_ss->SetMarkerStyle(20);
  C2Therm_os->SetMarkerStyle(25);
  C2Therm_ss->SetMarkerColor(2);
  C2Therm_os->SetMarkerColor(4);
  C2Therm_ss->SetLineColor(2);
  C2Therm_os->SetLineColor(4);
  C2Therm_ss->SetMarkerSize(1.5);
  C2Therm_os->SetMarkerSize(1.5);
  C2Therm_ss->SetDirectory(0);
  C2Therm_os->SetDirectory(0);
  C2Therm_ss->GetXaxis()->SetTitle("q_{inv} (GeV/c)");
  C2Therm_ss->GetYaxis()->SetTitle("1+<cos(qx)>");
  C2Therm_os->GetXaxis()->SetTitle("q_{inv} (GeV/c)");
  C2Therm_os->GetYaxis()->SetTitle("C_{2}^{+-}");
  C2Therm_os->SetMinimum(0.98);
  C2Therm_os->SetMaximum(1.25);
  //C2Therm_ss->DrawCopy();
  //C2Therm_os->DrawCopy("same");
  C2Therm_ss->SetDirectory(0);
  C2Therm_os->SetDirectory(0);
  //C2Therm->Fit(GaussEW,"IME","",0.005,0.1);
  Thermfile->Close();
  


  /////////////////////////////////////////////////////
  // Global fitting C2os and C2ss
  double ThermEW[4]={0};
  double C2ssSys_e[BINRANGE_C2global]={0};
  double C2osSys_e[BINRANGE_C2global]={0};
  TF1 *fitC2ssTherm;
  TF1 *fitC2osTherm;
  TF1 *fitC2ssThermGaus;
  TF1 *fitC2osThermGaus;
  //
  const int npar=11;// 10 
  TMinuit MyMinuit(npar);
  MyMinuit.SetFCN(fcnC2_Global);
  double OutputPar[npar]={0};
  double OutputPar_e[npar]={0};

  double par[npar];               // the start values
  double stepSize[npar];          // step sizes 
  double minVal[npar];            // minimum bound on parameter 
  double maxVal[npar];            // maximum bound on parameter
  string parName[npar];

  TH1D *temp_mm=(TH1D*)Two_ex[0][0][SCOI_2][0]->Clone();
  temp_mm->Divide(Two_ex[0][0][SCOI_2][1]);
  TH1D *temp_mp=(TH1D*)Two_ex[0][1][SCOI_2][0]->Clone();
  temp_mp->Divide(Two_ex[0][1][SCOI_2][1]);
  TH1D *temp_pp=(TH1D*)Two_ex[1][1][SCOI_2][0]->Clone();
  temp_pp->Divide(Two_ex[1][1][SCOI_2][1]);
  TH1D *C2ssRaw=(TH1D*)temp_mm->Clone();// MRC uncorrected
  TH1D *C2osRaw=(TH1D*)temp_mp->Clone();// MRC uncorrected
  C2ssRaw->SetMarkerStyle(24);
  C2osRaw->SetMarkerStyle(21);//21
  C2ssRaw->SetMarkerColor(2);
  C2osRaw->SetMarkerColor(4);
 
  
  int iteration=0;
  while(iteration<3){// 0=pure Gaussian Therminator fits, 1=EW Therminator fits, 2=real data fits with Therm EW
    if(!IncludeEWfromTherm) iteration=2;// skip Therminator 
   
    for(int i=0; i<BINRANGE_C2global; i++){
      
      C2ssFitting[i] = (temp_mm->GetBinContent(i+1) + temp_pp->GetBinContent(i+1))/2.;
      if(!GeneratedSignal && !MCcase) C2ssFitting[i] *= (MomRes_C2_pp[i]-1)*MRCShift+1;
      C2ssFitting_e[i] = pow(MomRes_C2_pp[i]*sqrt(pow(temp_mm->GetBinError(i+1),2) + pow(temp_pp->GetBinError(i+1),2))/sqrt(2.),2);
      C2ssRaw->SetBinContent(i+1, (temp_mm->GetBinContent(i+1) + temp_pp->GetBinContent(i+1))/2.);
      C2ssRaw->SetBinError(i+1, pow(sqrt(pow(temp_mm->GetBinError(i+1),2) + pow(temp_pp->GetBinError(i+1),2))/sqrt(2.),2));
      C2ssFitting_e[i] += pow((MomRes_C2_pp[i]-1)*0.1 * (temp_mm->GetBinContent(i+1) + temp_pp->GetBinContent(i+1))/2.,2);
      C2ssFitting_e[i] = sqrt(C2ssFitting_e[i]);
      C2osFitting[i] = temp_mp->GetBinContent(i+1);
      if(!GeneratedSignal && !MCcase) C2osFitting[i] *= (MomRes_C2_mp[i]-1)*MRCShift+1;
      C2osFitting_e[i] = pow(MomRes_C2_mp[i]*temp_mp->GetBinError(i+1),2);
      C2osRaw->SetBinContent(i+1, temp_mp->GetBinContent(i+1));
      C2osRaw->SetBinError(i+1, temp_mp->GetBinError(i+1));
      C2osFitting_e[i] += pow((MomRes_C2_mp[i]-1)*0.1 * temp_mp->GetBinContent(i+1),2);
      C2osFitting_e[i] = sqrt(C2osFitting_e[i]);
      //
      C2ssSys_e[i] = pow((MomRes_C2_pp[i]-1)*0.1 * (temp_mm->GetBinContent(i+1) + temp_pp->GetBinContent(i+1))/2.,2);
      C2ssSys_e[i] = sqrt(C2ssSys_e[i]);
      C2osSys_e[i] = pow((MomRes_C2_mp[i]-1)*0.1 * temp_mp->GetBinContent(i+1),2);
      C2osSys_e[i] = sqrt(C2osSys_e[i]);
      //
      
      K2SS[i] = CoulCorr2(+1, AvgQinvSS[i]);
      K2OS[i] = CoulCorr2(-1, AvgQinvOS[i]);
      //K2SS[i] = 1;
      //K2OS[i] = 1;
      //
      K2SS_e[i] = sqrt(pow((K2SS[i]-1)*0.02,2) + pow((K2SS[i]-Gamov(+1, AvgQinvSS[i]))*0.02,2));
      K2OS_e[i] = sqrt(pow((K2OS[i]-1)*0.02,2) + pow((K2OS[i]-Gamov(-1, AvgQinvOS[i]))*0.02,2));
      //K2SS_e[i] = 0.0001;
      //K2OS_e[i] = 0.0001;
      
      //
      if(iteration<2){
	C2ssFitting[i] = C2Therm_ss->GetBinContent(i+1);// Therminator fit
	C2ssFitting_e[i] = C2Therm_ss->GetBinError(i+1);// Therminator fit
	C2osFitting[i] = C2Therm_os->GetBinContent(i+1);// Therminator fit
	C2osFitting_e[i] = C2Therm_os->GetBinError(i+1);// Therminator fit
	K2SS_e[i] = 0.;
	K2OS_e[i] = 0.;
      }
      
    }
    
    
    
    C2ssFitting[0]=-1000; C2osFitting[0]=-1000;
    C2ssFitting_e[0]=1000; C2osFitting_e[0]=1000;
    K2SS_e[0]=1000; K2OS_e[0]=1000;
    
    
    
    par[0] = 1.0; par[1]=0.5; par[2]=0.5; par[3]=9.2; par[4] = .1; par[5] = .2; par[6] = .0; par[7] = 0.; par[8] = 0.; par[9] = 1.; par[10] = 0.01;
    stepSize[0] = 0.01; stepSize[1] = 0.01;  stepSize[2] = 0.02; stepSize[3] = 0.2; stepSize[4] = 0.01; stepSize[5] = 0.001; stepSize[6] = 0.001; stepSize[7] = 0.001; stepSize[8]=0.001; stepSize[9]=0.01; stepSize[10]=0.01; 
    minVal[0] = 0.995; minVal[1] = 0.2; minVal[2] = 0.; minVal[3] = 0.1; minVal[4] = 0.001; minVal[5] = -10.; minVal[6] = -10.; minVal[7] = -10.; minVal[8]=-10; minVal[9] = 0.995; minVal[10] = 0; 
    maxVal[0] = 1.03; maxVal[1] = 1.0; maxVal[2] = 0.99; maxVal[3] = 15.; maxVal[4] = 2.; maxVal[5] = 10.; maxVal[6] = 10.; maxVal[7] = 10.; maxVal[8]=10.; maxVal[9]=1.03; maxVal[10]=1.;
    parName[0] = "Norm"; parName[1] = "#Lambda"; parName[2] = "G"; parName[3] = "Rch"; parName[4] = "Rcoh"; 
    parName[5] = "kappa_3"; parName[6] = "kappa_4"; parName[7] = "kappa_5"; parName[8]="kappa_6"; parName[9]="Norm_2"; parName[10]="P_{coh}";
    
    MyMinuit.DefineParameter(10, parName[10].c_str(), 0., stepSize[10], minVal[10], maxVal[10]); MyMinuit.FixParameter(10);// extra parameter
    
    for (int i=0; i<npar; i++){
      MyMinuit.DefineParameter(i, parName[i].c_str(), par[i], stepSize[i], minVal[i], maxVal[i]);
    }
    //MyMinuit.DefineParameter(0, parName[0].c_str(), .998, stepSize[0], minVal[0], maxVal[0]); MyMinuit.FixParameter(0);// N
    //MyMinuit.DefineParameter(1, parName[1].c_str(), .68, stepSize[1], minVal[1], maxVal[1]); MyMinuit.FixParameter(1);// lambda
    //MyMinuit.DefineParameter(2, parName[2].c_str(), 0.9999, stepSize[2], minVal[2], maxVal[2]); MyMinuit.FixParameter(2);// G
    
    if(IncludeG==kFALSE || iteration<2) {
      MyMinuit.DefineParameter(2, parName[2].c_str(), 0., stepSize[2], minVal[2], maxVal[2]); MyMinuit.FixParameter(2);// G
      MyMinuit.DefineParameter(4, parName[4].c_str(), 0., stepSize[4], minVal[4], maxVal[4]); MyMinuit.FixParameter(4);// Rcoh
    }else{
      Int_t err=0;
      Double_t tmp[1];
      tmp[0] = 2+1;
      MyMinuit.mnexcm( "RELease", tmp,  1,  err );// G
      tmp[0] = 4+1;
      MyMinuit.mnexcm( "RELease", tmp,  1,  err );// Rcoh
    }
    
    if(!IncludeEW){
      if(!IncludeEWfromTherm){
	// kappa3=0.24, kappa4=0.16 for testing
	MyMinuit.DefineParameter(5, parName[5].c_str(), 0, stepSize[5], minVal[5], maxVal[5]); MyMinuit.FixParameter(5);// k3 
	MyMinuit.DefineParameter(6, parName[6].c_str(), 0, stepSize[6], minVal[6], maxVal[6]); MyMinuit.FixParameter(6);// k4
	MyMinuit.DefineParameter(7, parName[7].c_str(), 0, stepSize[7], minVal[7], maxVal[7]); MyMinuit.FixParameter(7);// k5
	MyMinuit.DefineParameter(8, parName[8].c_str(), 0, stepSize[8], minVal[8], maxVal[8]); MyMinuit.FixParameter(8);// k6
      }else{
	if(iteration==0){
	  MyMinuit.DefineParameter(5, parName[5].c_str(), 0, stepSize[5], minVal[5], maxVal[5]); MyMinuit.FixParameter(5);// k3 
	  MyMinuit.DefineParameter(6, parName[6].c_str(), 0, stepSize[6], minVal[6], maxVal[6]); MyMinuit.FixParameter(6);// k4
	}else if(iteration==1){
	  Int_t err=0;
	  Double_t tmp[1];
	  tmp[0] = 5+1;
	  MyMinuit.mnexcm( "RELease", tmp,  1,  err );// k3
	  tmp[0] = 6+1;
	  MyMinuit.mnexcm( "RELease", tmp,  1,  err );// k4
	}else{// iteration=2
	  MyMinuit.DefineParameter(5, parName[5].c_str(), ThermEW[0], stepSize[5], minVal[5], maxVal[5]); MyMinuit.FixParameter(5);// k3 
	  MyMinuit.DefineParameter(6, parName[6].c_str(), ThermEW[1], stepSize[6], minVal[6], maxVal[6]); MyMinuit.FixParameter(6);// k4
	  MyMinuit.DefineParameter(7, parName[7].c_str(), ThermEW[2], stepSize[7], minVal[7], maxVal[7]); MyMinuit.FixParameter(7);// k5
	  MyMinuit.DefineParameter(8, parName[8].c_str(), ThermEW[3], stepSize[8], minVal[8], maxVal[8]); MyMinuit.FixParameter(8);// k6
	}
      }// IncludeEWfromTherm
    }// IncludeEW
    
    if(IncludeSS==kFALSE){
      MyMinuit.DefineParameter(3, parName[3].c_str(), 9.1, stepSize[3], minVal[3], maxVal[3]); MyMinuit.FixParameter(3);// Rch
      MyMinuit.DefineParameter(0, parName[0].c_str(), .998, stepSize[0], minVal[0], maxVal[0]); MyMinuit.FixParameter(0);// N
      MyMinuit.DefineParameter(5, parName[5].c_str(), 0, stepSize[5], minVal[5], maxVal[5]); MyMinuit.FixParameter(5);// k3 
      MyMinuit.DefineParameter(6, parName[6].c_str(), 0, stepSize[6], minVal[6], maxVal[6]); MyMinuit.FixParameter(6);// k4
    }
    if(IncludeOS==kFALSE){
      MyMinuit.DefineParameter(9, parName[9].c_str(), 1.002, stepSize[9], minVal[9], maxVal[9]); MyMinuit.FixParameter(9);// N_2
    }
    
    //MyMinuit.DefineParameter(3, parName[3].c_str(), 10.5, stepSize[3], minVal[3], maxVal[3]); MyMinuit.FixParameter(3);// Rch
    //MyMinuit.DefineParameter(4, parName[4].c_str(), 5., stepSize[4], minVal[4], maxVal[4]); MyMinuit.FixParameter(4);// Rcoh
    //MyMinuit.DefineParameter(5, parName[5].c_str(), 0, stepSize[5], minVal[5], maxVal[5]); MyMinuit.FixParameter(5);// k3 
    //MyMinuit.DefineParameter(6, parName[6].c_str(), 0, stepSize[6], minVal[6], maxVal[6]); MyMinuit.FixParameter(6);// k4
    MyMinuit.DefineParameter(7, parName[7].c_str(), 0, stepSize[7], minVal[7], maxVal[7]); MyMinuit.FixParameter(7);// k5
    MyMinuit.DefineParameter(8, parName[8].c_str(), 0, stepSize[8], minVal[8], maxVal[8]); MyMinuit.FixParameter(8);// k6
    //
    int ierflg=0;
    double arglist[10];
    //arglist[0]=2;// improve Minimization Strategy (1 is default)
    //MyMinuit.mnexcm("SET STR",arglist,1,ierflg);
    //arglist[0] = 0;
    //MyMinuit.mnexcm("SCAN", arglist,1,ierflg);
    arglist[0] = 5000;
    MyMinuit.mnexcm("MIGRAD", arglist ,1,ierflg);
    // Do the minimization!
    cout<<"Start C2 Global fit"<<endl;
    MyMinuit.Migrad();// generally the best minimization algorithm
    cout<<"End C2 Global fit"<<endl;
    
    for (int i=0; i<npar; i++){
      MyMinuit.GetParameter(i,OutputPar[i],OutputPar_e[i]);
    }
    
    cout<<"Chi2/NDF = "<<Chi2_C2global/(NFitPoints_C2global - MyMinuit.GetNumFreePars())<<"   Chi^2="<<Chi2_C2global<<endl;
    if(iteration==2) {
      ChiSquaredNDF->Fill(1, Chi2_C2global/(NFitPoints_C2global - MyMinuit.GetNumFreePars()));
    }if(iteration==1) {
      ThermEW[0]=OutputPar[5]; ThermEW[1]=OutputPar[6]; ThermEW[2]=OutputPar[7]; ThermEW[3]=OutputPar[8];
      fitC2ssTherm = new TF1("fitC2ssTherm",C2ssFitFunction, 0.005,0.2, npar);
      for(int i=0; i<npar; i++) {
	fitC2ssTherm->FixParameter(i,OutputPar[i]);
	fitC2ssTherm->SetParError(i,OutputPar_e[i]);
      }
      fitC2osTherm = new TF1("fitC2osTherm",C2osFitFunction, 0.005,0.2, npar);
      for(int i=0; i<npar; i++) {
	fitC2osTherm->FixParameter(i,OutputPar[i]);
	fitC2osTherm->SetParError(i,OutputPar_e[i]);
      }
      fitC2osTherm->SetLineColor(4);
    }
    if(iteration==0){
      fitC2ssThermGaus = new TF1("fitC2ssThermGaus",C2ssFitFunction, 0.005,0.2, npar);
      for(int i=0; i<npar; i++) {
	fitC2ssThermGaus->FixParameter(i,OutputPar[i]);
	fitC2ssThermGaus->SetParError(i,OutputPar_e[i]);
      }
      fitC2osThermGaus = new TF1("fitC2osThermGaus",C2osFitFunction, 0.005,0.2, npar);
      for(int i=0; i<npar; i++) {
	fitC2osThermGaus->FixParameter(i,OutputPar[i]);
	fitC2osThermGaus->SetParError(i,OutputPar_e[i]);
      }
      fitC2osThermGaus->SetLineColor(4);
    }
    
    iteration++;
    
  }// iteration loop
  
 
  TH1D *C2_ss=(TH1D*)Two_ex_clone_mm->Clone();
  TH1D *C2_os=(TH1D*)Two_ex_clone_mp->Clone();
  TH1D *C2_ss_Momsys=(TH1D*)Two_ex_clone_mm->Clone();
  TH1D *C2_os_Momsys=(TH1D*)Two_ex_clone_mp->Clone();
  TH1D *C2_ss_Ksys=(TH1D*)Two_ex_clone_mm->Clone();
  TH1D *C2_os_Ksys=(TH1D*)Two_ex_clone_mp->Clone();
  TH1D *K2_ss = (TH1D*)Two_ex_clone_mm->Clone();
  TH1D *K2_os = (TH1D*)Two_ex_clone_mp->Clone();

  for(int i=0; i<BINRANGE_C2global; i++){ 
    C2_ss->SetBinContent(i+1, C2ssFitting[i]);
    C2_os->SetBinContent(i+1, C2osFitting[i]);
    C2_ss_Momsys->SetBinContent(i+1, C2ssFitting[i]);
    C2_os_Momsys->SetBinContent(i+1, C2osFitting[i]);
    C2_ss_Ksys->SetBinContent(i+1, C2ssFitting[i]);
    C2_os_Ksys->SetBinContent(i+1, C2osFitting[i]);
    double Toterror_ss = sqrt(fabs(pow(C2ssFitting_e[i],2)-pow(C2ssSys_e[i],2)));
    double Toterror_os = sqrt(fabs(pow(C2osFitting_e[i],2)-pow(C2osSys_e[i],2)));
    C2_ss->SetBinError(i+1, Toterror_ss);
    C2_os->SetBinError(i+1, Toterror_os);
    C2_ss_Momsys->SetBinError(i+1, C2ssSys_e[i]);
    C2_os_Momsys->SetBinError(i+1, C2osSys_e[i]);
    C2_ss_Ksys->SetBinError(i+1, OutputPar[1]*K2SS_e[i]);
    C2_os_Ksys->SetBinError(i+1, OutputPar[1]*K2OS_e[i]);
    //
    K2_ss->SetBinContent(i+1, K2SS[i]); K2_ss->SetBinError(i+1, 0);
    K2_os->SetBinContent(i+1, K2OS[i]); K2_os->SetBinError(i+1, 0);
  }

  C2_ss_Momsys->SetMarkerSize(0);
  C2_ss_Momsys->SetFillStyle(1000);
  C2_ss_Momsys->SetFillColor(kRed-10);
  C2_os_Momsys->SetMarkerSize(0);
  C2_os_Momsys->SetFillStyle(1000);
  C2_os_Momsys->SetFillColor(kBlue-10);
  C2_ss_Ksys->SetMarkerSize(0);
  C2_ss_Ksys->SetFillStyle(3004);
  C2_ss_Ksys->SetFillColor(kRed);
  C2_os_Ksys->SetMarkerSize(0);
  C2_os_Ksys->SetFillStyle(3004);
  C2_os_Ksys->SetFillColor(kBlue);
  
  
  
  C2_ss->GetXaxis()->SetRangeUser(0,0.098);
  C2_ss->GetYaxis()->SetRangeUser(0.98,1.3);
  C2_ss->GetXaxis()->SetTitleOffset(.8);
  C2_ss->GetYaxis()->SetTitleOffset(.8);
  C2_ss->GetXaxis()->SetTitle("q_{inv} (GeV/c)");
  C2_ss->GetXaxis()->SetTitleSize(0.05);
  C2_ss->GetYaxis()->SetTitleSize(0.05);
  C2_os->GetXaxis()->SetRangeUser(0,0.098);
  C2_os->GetYaxis()->SetRangeUser(0.98,1.3);
  C2_os->GetXaxis()->SetTitleOffset(.8);
  C2_os->GetYaxis()->SetTitleOffset(.8);
  C2_os->GetXaxis()->SetTitle("q_{inv} (GeV/c)");
  C2_os->GetXaxis()->SetTitleSize(0.05);
  C2_os->GetYaxis()->SetTitleSize(0.05);

  C2_ss->SetMarkerSize(1.5);
  C2_os->SetMarkerSize(1.5);
  C2_os->SetMarkerStyle(25);
  C2_os->SetMarkerColor(4);
 
  TF1 *fitC2ss = new TF1("fitC2ss",C2ssFitFunction, 0.001,0.2, npar);
  TF1 *fitC2os = new TF1("fitC2os",C2osFitFunction, 0.001,0.2, npar);
  for(int i=0; i<npar; i++) {
    fitC2ss->FixParameter(i,OutputPar[i]);
    fitC2ss->SetParError(i,OutputPar_e[i]);
  }
  for(int i=0; i<npar; i++) {
    fitC2os->FixParameter(i,OutputPar[i]);
    fitC2os->SetParError(i,OutputPar_e[i]);
  }
  TH1D *fitC2ss_h=new TH1D("fitC2ss_h","",C2_ss->GetNbinsX(),0,C2_ss->GetXaxis()->GetBinUpEdge(C2_ss->GetNbinsX()));
  TH1D *fitC2os_h=new TH1D("fitC2os_h","",C2_os->GetNbinsX(),0,C2_os->GetXaxis()->GetBinUpEdge(C2_os->GetNbinsX()));
  fitC2ss_h->SetLineWidth(2); fitC2os_h->SetLineWidth(2);
  fitC2ss_h->SetLineColor(2); fitC2os_h->SetLineColor(4);

  for(int bin=1; bin<=fitC2ss_h->GetNbinsX(); bin++){
    double qinv=(bin-0.5)*0.005;
    fitC2ss_h->SetBinContent(bin, fitC2ss->Eval(qinv));
    fitC2os_h->SetBinContent(bin, fitC2os->Eval(qinv));
  }

  if(!MCcase){
    C2_ss->DrawCopy();
    legend1->AddEntry(C2_ss,"same-charge","p");
    C2_os->DrawCopy("same");
    legend1->AddEntry(C2_os,"opp-charge","p");
    C2_ss_Momsys->DrawCopy("E2 same");
    C2_os_Momsys->DrawCopy("E2 same");
    C2_ss_Ksys->DrawCopy("E2 same");
    C2_os_Ksys->DrawCopy("E2 same");
    C2_ss->DrawCopy("same");
    C2_os->DrawCopy("same");
    fitC2os->SetLineColor(4);
    fitC2ss_h->DrawCopy("C same");
    fitC2os_h->DrawCopy("C same");
    MJ_bkg_ss->SetLineColor(1);
    MJ_bkg_ss->Draw("same");
    legend1->AddEntry(MJ_bkg_ss,"Non-femto bkg","l");
    legend1->Draw("same");
  }
  //C2Therm_ss->DrawCopy();
  //C2Therm_os->DrawCopy("same");

  /*
  ////////// C2 systematics
  // C2 -- base, M0, Pos B (0.035 TTC for all)
  //double C2ss_base[20]={0, 1.17121, 1.17272, 1.13947, 1.09819, 1.06512, 1.04115, 1.02563, 1.01558, 1.00937, 1.00551, 1.00309, 1.00157, 1.00065, 1.00023, 0.999911, 0.99976, 0.999733, 0.999593, 0.999641};
  //double C2ss_base_e[20]={0, 0.00135764, 0.000430769, 0.000242626, 0.000167895, 0.00012913, 0.000105724, 9.01743e-05, 7.90587e-05, 7.07036e-05, 6.41813e-05, 5.89418e-05, 5.46406e-05, 5.10458e-05, 4.80009e-05, 4.53857e-05, 4.31177e-05, 4.11375e-05, 3.93869e-05, 3.78413e-05};
  // --, M0, 
  double C2ss_base[20]={0, 1.16926, 1.17306, 1.13932, 1.09834, 1.06475, 1.04096, 1.02546, 1.01547, 1.00931, 1.00542, 1.00311, 1.00155, 1.00071, 1.00018, 0.999786, 0.999742, 0.999612, 0.999548, 0.999631};
  double C2ss_base_e[20]={0, 0.00128558, 0.000410643, 0.000231512, 0.000160337, 0.000123314, 0.000100987, 8.61441e-05, 7.55245e-05, 6.75439e-05, 6.13117e-05, 5.63119e-05, 5.21976e-05, 4.87656e-05, 4.58515e-05, 4.33492e-05, 4.11861e-05, 3.92891e-05, 3.76189e-05, 3.61428e-05};
  // -+, M0, 
  //double C2ss_base[20]={1.97333, 1.17447, 1.09412, 1.05561, 1.03613, 1.02507, 1.01809, 1.01363, 1.01048, 1.00816, 1.00645, 1.00526, 1.00427, 1.0035, 1.00289, 1.00241, 1.00209, 1.0018, 1.00148, 1.00127};
  //double C2ss_base_e[20]={1.61122, 0.00030054, 0.000174753, 0.000123051, 9.53599e-05, 7.81537e-05, 6.64323e-05, 5.79583e-05, 5.15436e-05, 4.65277e-05, 4.25025e-05, 3.92115e-05, 3.64672e-05, 3.4149e-05, 3.21673e-05, 3.04575e-05, 2.89698e-05, 2.76642e-05, 2.65087e-05, 2.54847e-05};
  // --, M0, kt4
  //double C2ss_base[20]={0, 0, 0, 1.15387, 1.14361, 1.10213, 1.06552, 1.04357, 1.02804, 1.01826, 1.01182, 1.00767, 1.00509, 1.00294, 1.00169, 1.0008, 1.00063, 1.00022, 1.00001, 0.999986};
  //double C2ss_base_e[20]={0, 0, 0, 0.00229526, 0.000813871, 0.000518539, 0.00039803, 0.000328676, 0.000281681, 0.00024771, 0.000221807, 0.000201394, 0.000184824, 0.000171048, 0.000159468, 0.000149554, 0.000141058, 0.000133583, 0.00012702, 0.000121218};
  // --, M0, kt5
  //double C2ss_base[20]={0, 0, 0, 1.11133, 1.13824, 1.11806, 1.08288, 1.05648, 1.03774, 1.02516, 1.01728, 1.01166, 1.00841, 1.00553, 1.00357, 1.00194, 1.00132, 1.00104, 1.00058, 1.00039};
  //double C2ss_base_e[20]={0, 0, 0, 0.0147639, 0.00188068, 0.000910766, 0.000637117, 0.000510753, 0.000432397, 0.000377855, 0.000337624, 0.000306333, 0.000281495, 0.00026099, 0.000243881, 0.000229331, 0.000216961, 0.000206243, 0.000196785, 0.000188487};
  // --, M0, kt6
  //double C2ss_base[20]={0, 0, 0, 0, 1.08106, 1.12109, 1.09838, 1.07193, 1.05087, 1.03534, 1.02433, 1.01785, 1.01339, 1.00928, 1.00662, 1.00441, 1.00326, 1.00192, 1.00175, 1.00143};
  //double C2ss_base_e[20]={0, 0, 0, 0, 0.00722887, 0.00211091, 0.00123591, 0.000924715, 0.000765679, 0.000662277, 0.000588573, 0.000533669, 0.000490831, 0.000456049, 0.000427965, 0.000404612, 0.000385285, 0.000368671, 0.00035474, 0.000342708};

  cout<<"C2 values:"<<endl;
  for(int ii=0; ii<20; ii++){
    cout<<Two_ex_clone_mm->GetBinContent(ii+1)<<", ";
  }
  cout<<endl;
  cout<<"C2 errors:"<<endl;
  for(int ii=0; ii<20; ii++){
    cout<<Two_ex_clone_mm->GetBinError(ii+1)<<", ";
  }
  cout<<endl;
  TH1D *C2ss_Sys=(TH1D*)Two_ex_clone_mm->Clone();
  for(int ii=1; ii<20; ii++){
    if(C2ss_base[ii] ==0) {
      C2ss_Sys->SetBinContent(ii+1, 100.);
      C2ss_Sys->SetBinError(ii+1, 100.);
      continue;
    }
    C2ss_Sys->SetBinContent(ii+1, (C2ss_base[ii]-C2ss_Sys->GetBinContent(ii+1))/C2ss_base[ii]);
    C2ss_Sys->SetBinError(ii+1, sqrt(fabs(pow(C2ss_Sys->GetBinError(ii+1),2) - pow(C2ss_base_e[ii],2))));
  }
  C2ss_Sys->GetXaxis()->SetRangeUser(0,0.099);
  C2ss_Sys->GetYaxis()->SetTitle("(C_{2}^{s}-C_{2}^{v})/C_{2}^{s}");
  C2ss_Sys->GetYaxis()->SetTitleOffset(1.3);
  C2ss_Sys->SetMinimum(-0.01);
  C2ss_Sys->SetMaximum(0.01);
  C2ss_Sys->Draw();
  TF1 *C2lineSys=new TF1("C2lineSys","pol0",0,0.099);
  //C2ss_Sys->Fit(C2lineSys,"MEQ","",0,0.099);
  */

  // Momentum resolution uncorrected C2
  /*C2ssRaw->Draw("same");
  C2osRaw->Draw("same");
  legend1->AddEntry(C2ssRaw,"same-charge, Raw","p");
  legend1->AddEntry(C2osRaw,"opp-charge, Raw","p");
  legend1->Draw("same");
  */

  // FSI corrected C2+-
  /*TH1D *C2_os_FSIcorrected=(TH1D*)C2_os->Clone();
  for(int ii=2; ii<20; ii++){
    double value = (C2_os_FSIcorrected->GetBinContent(ii) - (1-OutputPar[1]))/(OutputPar[1]*K2OS[ii-1]);
    C2_os_FSIcorrected->SetBinContent(ii,value);
  }
  C2_os_FSIcorrected->SetMinimum(0.95);
  C2_os_FSIcorrected->SetMarkerStyle(24);
  C2_os_FSIcorrected->Draw();
  */

  
  /*Two_ex_clone_mm->SetMarkerStyle(25);
  Two_ex_clone_mp->SetMarkerStyle(25);
  Two_ex_clone_mm->SetMarkerColor(2);
  Two_ex_clone_mp->SetMarkerColor(4);
  Two_ex_clone_mm->Draw("same");
  Two_ex_clone_mp->Draw("same");
  legend1->AddEntry(C2_ss,"same-charge; Momentum Resolution Corrected","p");
  legend1->AddEntry(C2_os,"opp-charge; Momentum Resolution Corrected","p");
  legend1->AddEntry(Two_ex_clone_mm,"same-charge; Raw","p");
  legend1->AddEntry(Two_ex_clone_mp,"opp-charge; Raw","p");
  legend1->Draw("same");
  */
  //MyMinuit.SetErrorDef(4); //note 4 and not 2!
  //TGraph *gr2 = (TGraph*)MyMinuit.Contour(15,1,2);
  //gr2->GetXaxis()->SetTitle("lambda");
  //gr2->GetYaxis()->SetTitle("G");
  //gr2->SetTitle("");
  //gr2->SetFillColor(kYellow);
  //gr2->Draw("alf");
  //Get contour for parameter 0 versus parameter 2 for ERRDEF=1  
  //MyMinuit.SetErrorDef(1);
  //TGraph *gr1 = (TGraph*)MyMinuit.Contour(15,1,2);
  //gr1->SetFillColor(kGreen);//38
  //gr1->Draw("lf");
  

 
  ///////////////////////////
  // STAR comparison
  /* double C2_star_mm[50]={0, 0, 0.566654, 1.11195, 1.09964, 1.12566, 1.1479, 1.15596, 1.15658, 1.14583, 1.13182, 1.11544, 1.09851, 1.08326, 1.0701, 1.05686, 1.048, 1.03924, 1.03122, 1.02544, 1.02019, 1.01617, 1.01239, 1.00942, 1.00736, 1.00485, 1.00296, 1.00152, 1.00039, 0.999461, 0.998703, 0.997563, 0.997294, 0.996726, 0.996739, 0.996355, 0.996426, 0.996581, 0.996336, 0.996157, 0.996312, 0.996413, 0.996501, 0.996204, 0.996816, 0.996654, 0.99695, 0.997494, 0.997135, 0.997188};
  double C2_star_mm_e[50]={0, 0, 0.454334, 0.0216601, 0.00590084, 0.00297232, 0.001939, 0.00142608, 0.00112851, 0.000929702, 0.000792097, 0.000690265, 0.000613054, 0.000552755, 0.000504494, 0.000464031, 0.000431371, 0.000403169, 0.000378672, 0.000357685, 0.000339133, 0.00032286, 0.000308251, 0.000295149, 0.000283433, 0.000272615, 0.000262784, 0.000253889, 0.000245693, 0.00023813, 0.000231135, 0.00022455, 0.000218593, 0.000212961, 0.000207822, 0.000202906, 0.000198339, 0.000194087, 0.000189971, 0.000186181, 0.000182575, 0.000179202, 0.000175989, 0.000172894, 0.000170072, 0.000167293, 0.000164706, 0.000162286, 0.000159867, 0.000157606};
  TH1D *Two_star_mm=(TH1D*)Two_ex_clone_mm->Clone();
  Two_star_mm->Add(Two_star_mm,-1);
  for(int i=0; i<50; i++) {Two_star_mm->SetBinContent(i+1, C2_star_mm[i]); Two_star_mm->SetBinError(i+1, C2_star_mm_e[i]);}
  Two_star_mm->SetMarkerColor(4);
  Two_star_mm->SetLineColor(4);
  //Two_star_mm->Draw("same");
  //legend1->AddEntry(Two_star_mm,"--, 200 GeV","p");

  //Two_ex_clone_mp->Multiply(Two_ex_clone_pp);
  //Two_ex_clone_mp->Draw();
  //Two_ex_clone_mm->Multiply(Two_ex_clone_mp);
  //Two_ex_clone_mm->Draw();
  //Two_ex_clone_pp->Draw("same");
  
  //legend1->Draw("same");
  */


  
 

  /////////////////////////////////////////////////////////////////////////
  // 3-d fitting (out,side,long).  Not used for paper.
  /*
  TString *name1 = new TString("Explicit2_Charge1_0_Charge2_0_SC_0_M_");
  *name1 += Mbin;
  TString *name2 = new TString(name1->Data());
  TString *name3 = new TString(name1->Data());
  name1->Append("_ED_0_Term_1_osl_b2");// b1 (0.2<kt<0.3). b2 (0.6<kt<0.7)
  name2->Append("_ED_0_Term_2_osl_b2");
  name3->Append("_ED_0_Term_1_osl_b2_QW");
  

  TH3D *num_osl = (TH3D*)MyList->FindObject(name1->Data());
  TH3D *den_osl = (TH3D*)MyList->FindObject(name2->Data());
  den_osl->Scale(num_osl->Integral(28,40, 28,40, 28,40)/den_osl->Integral(28,40, 28,40, 28,40));
  TH3D *num_osl_QW = (TH3D*)MyList->FindObject(name3->Data());
  
  for(int i=0; i<BINS_OSL; i++){
    for(int j=0; j<BINS_OSL; j++){
      for(int k=0; k<BINS_OSL; k++){
	if(num_osl->GetBinContent(i+1,j+1,k+1) < 1 || den_osl->GetBinContent(i+1,j+1,k+1) < 1) continue;
	
	avg_q[i][j][k] = num_osl_QW->GetBinContent(i+1,j+1,k+1)/num_osl->GetBinContent(i+1,j+1,k+1);
	int QinvBin = int((avg_q[i][j][k])/0.005);
	if(QinvBin >=20) QinvBin=19;
	if(MomRes_term1_pp[QinvBin] ==0) continue;
	if(MomRes_term2_pp[QinvBin] ==0) continue;
	//
	num_osl->SetBinContent(i+1,j+1,k+1, num_osl->GetBinContent(i+1,j+1,k+1)*MomRes_term1_pp[QinvBin]);
	den_osl->SetBinContent(i+1,j+1,k+1, den_osl->GetBinContent(i+1,j+1,k+1)*MomRes_term2_pp[QinvBin]);
	//
	A[i][j][k] = num_osl->GetBinContent(i+1,j+1,k+1);
	B[i][j][k] = den_osl->GetBinContent(i+1,j+1,k+1);
	//
	A_e[i][j][k] = num_osl->GetBinContent(i+1,j+1,k+1);
	A_e[i][j][k] += pow(num_osl->GetBinContent(i+1,j+1,k+1)/MomRes_term1_pp[QinvBin] *fabs(MomRes_term1_pp[QinvBin]-1)*0.1,2);
	A_e[i][j][k] = sqrt(A_e[i][j][k]);
	B_e[i][j][k] = den_osl->GetBinContent(i+1,j+1,k+1);
	B_e[i][j][k] += pow(den_osl->GetBinContent(i+1,j+1,k+1)/MomRes_term2_pp[QinvBin] *fabs(MomRes_term2_pp[QinvBin]-1)*0.1,2);
	B_e[i][j][k] = sqrt(B_e[i][j][k]);
	//
	K_OSL[i][j][k] = CoulCorr2(+1, avg_q[i][j][k]);
	K_OSL_e[i][j][k] = sqrt(pow((K_OSL[i][j][k]-1)*0.02,2) + pow((K_OSL[i][j][k]-Gamov(+1,avg_q[i][j][k]))*0.01,2));
	//K_OSL_e[i][j][k] = 0;
      }
    }
  }
  
  
  const int npar_osl=6;//5
  TMinuit MyMinuit_osl(npar_osl);
  MyMinuit_osl.SetFCN(fcnOSL);
  double OutputPar_osl[npar_osl]={0};
  double OutputPar_osl_e[npar_osl]={0};

  double par_osl[npar_osl];               // the start values
  double stepSize_osl[npar_osl];          // step sizes 
  double minVal_osl[npar_osl];            // minimum bound on parameter 
  double maxVal_osl[npar_osl];            // maximum bound on parameter
  string parName_osl[npar_osl];
  
  par_osl[0] = 1.0; par_osl[1]=0.5; par_osl[2]=0.; par_osl[3]=5.; par_osl[4] = 5.; par_osl[5]=5.;
  stepSize_osl[0] = 0.01; stepSize_osl[1] = 0.01; stepSize_osl[2] = 0.02; stepSize_osl[3] = 0.2; stepSize_osl[4] = 0.2; stepSize_osl[5] = 0.2;
  minVal_osl[0] = 0.8; minVal_osl[1] = 0.01; minVal_osl[2] = 0.0; minVal_osl[3] = 1.; minVal_osl[4] = 1.; minVal_osl[5] = 1.;
  maxVal_osl[0] = 1.2; maxVal_osl[1] = 1.0; maxVal_osl[2] = 0.99; maxVal_osl[3] = 20.; maxVal_osl[4] = 20.; maxVal_osl[5] = 20.;
  parName_osl[0] = "Norm"; parName_osl[1] = "Lambda"; parName_osl[2] = "G"; parName_osl[3] = "R_out"; parName_osl[4] = "R_side"; parName_osl[5] = "R_long";
  
  for (int i=0; i<npar_osl; i++){
    MyMinuit_osl.DefineParameter(i, parName_osl[i].c_str(), par_osl[i], stepSize_osl[i], minVal_osl[i], maxVal_osl[i]);
  }
  MyMinuit_osl.FixParameter(2);// G
  //MyMinuit.FixParameter(1);// lambda
  cout<<"here!!"<<endl;

  // Do the minimization!
  cout<<"Start Three-d fit"<<endl;
  MyMinuit_osl.Migrad();// Minuit's best minimization algorithm
  cout<<"End Three-d fit"<<endl;
  cout<<"Chi2/NDF = "<<Chi2_OSL/(NFitPoints_OSL-MyMinuit_osl.GetNumFreePars())<<endl;

  for (int i=0; i<npar_osl; i++){
    MyMinuit_osl.GetParameter(i,OutputPar_osl[i],OutputPar_osl_e[i]);
  }
  
  
  TF3 *fit3d = new TF3("fit3d",OSLfitfunction, 0,0.2, 0,0.2, 0,0.2, npar_osl);
  for(int i=0; i<npar_osl; i++) {fit3d->FixParameter(i,OutputPar_osl[i]);}
  Int_t BinsOSL = num_osl->GetNbinsX();
  double LimitOSL = num_osl->GetXaxis()->GetBinUpEdge(BinsOSL);
  TH3D *den_timesFit = new TH3D("den_timesFit","",BinsOSL,0,LimitOSL, BinsOSL,0,LimitOSL, BinsOSL,0,LimitOSL);
  for(int i=0; i<BinsOSL; i++){
    for(int j=0; j<BinsOSL; j++){
      for(int k=0; k<BinsOSL; k++){
	
	double qout=(i+.5)*0.005;
	double qside=(j+.5)*0.005;
	double qlong=(k+.5)*0.005;
	den_timesFit->SetBinContent(i+1,j+1,k+1, fit3d->Eval(qout,qside,qlong)*den_osl->GetBinContent(i+1,j+1,k+1));
	den_timesFit->SetBinError(i+1,j+1,k+1, 0);
      }
    }
  }
  
  int binL=1, binH=4;
  TH1D *num_pro=(TH1D*)num_osl->ProjectionX("num_pro",binL,binH, binL,binH);
  TH1D *den_pro=(TH1D*)den_osl->ProjectionX("den_pro",binL,binH, binL,binH);
  TH1D *dentimesFit_pro=(TH1D*)den_timesFit->ProjectionX("dentimesFit_pro",binL,binH, binL,binH);
  num_pro->Sumw2();
  den_pro->Sumw2();
  num_pro->Divide(den_pro);
  num_pro->SetMarkerStyle(20);
  num_pro->SetTitle("");
  num_pro->GetXaxis()->SetTitle("q_{out}");
  num_pro->GetYaxis()->SetTitle("C_{2}");
  num_pro->SetMinimum(0.97);
  num_pro->SetMaximum(1.48);
  num_pro->DrawCopy();
  //
  dentimesFit_pro->Divide(den_pro);
  dentimesFit_pro->SetLineColor(2);
  dentimesFit_pro->DrawCopy("same");
  //
  */
  //MyMinuit.SetErrorDef(4); //note 4 and not 2!
  //TGraph *gr2 = (TGraph*)MyMinuit.Contour(10,1,2);
  //gr2->SetFillColor(kYellow);
  //gr2->Draw("alf");
  //Get contour for parameter 0 versus parameter 2 for ERRDEF=1  
  //MyMinuit.SetErrorDef(1);
  //TGraph *gr1 = (TGraph*)MyMinuit.Contour(10,1,2);
  //gr1->SetFillColor(kGreen);//38
  //gr1->Draw("lf");
  

 
  //////////////////////////////////////////////////////////////////////////////////////
  
  
  // To visualize the Qcut and Norm Regions
  //TH1D *QcutRegion = new TH1D("QcutRegion","",400,0,2);
  //TH1D *NormRegion1 = new TH1D("NormRegion1","",400,0,2);  
  //TH1D *NormRegion2 = new TH1D("NormRegion2","",400,0,2);  
  //for(int bin=1; bin<=20; bin++) QcutRegion->SetBinContent(bin,Two_ex[0][0][0][0]->GetBinContent(bin));
  //for(int bin=213; bin<=220; bin++) NormRegion1->SetBinContent(bin,Two_ex[0][0][0][0]->GetBinContent(bin));
  //for(int bin=31; bin<=35; bin++) NormRegion2->SetBinContent(bin,Two_ex[0][0][0][0]->GetBinContent(bin));
  //QcutRegion->SetFillColor(4);
  //NormRegion1->SetFillColor(2);
  //NormRegion2->SetFillColor(3);
  //Two_ex[0][0][0][0]->Draw();
  //QcutRegion->Draw("same");
  //NormRegion1->Draw("same");
  //NormRegion2->Draw("same");

 
 

  TCanvas *can2 = new TCanvas("can2", "can2",800,0,900,900);//800,0,600,900
  can2->SetHighLightColor(2);
  gStyle->SetOptFit(0111);
  can2->SetFillColor(10);
  can2->SetBorderMode(0);
  can2->SetBorderSize(2);
  can2->SetGridx();
  can2->SetGridy();
  can2->SetFrameFillColor(0);
  can2->SetFrameBorderMode(0);
  can2->SetFrameBorderMode(0);
  can2->cd();
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.02);
  
  TLegend *legend2 = new TLegend(.4,.67,1.,.87,NULL,"brNDC");
  legend2->SetBorderSize(1);
  legend2->SetTextSize(.06);// small .03; large .06
  legend2->SetFillColor(0);

  /////////////////////////////////////////////
  TH3D *C3QS_3d = new TH3D("C3QS_3d","",20,0,.98, 20,0,.1, 20,0,.1);
  TH3D *Combinatorics_3d = new TH3D("Combinatorics_3d","",20,0,.1, 20,0,.1, 20,0,.1);
  //
  const int Q3BINS = 20;
  TH1D *num_withRS = new TH1D("num_withRS","",Q3BINS,0,0.2);
  TH1D *num_withGRS = new TH1D("num_withGRS","",Q3BINS,0,0.2);
  TH1D *num_withTM = new TH1D("num_withTM","",Q3BINS,0,0.2);
  TH1D *Cterm1 = new TH1D("Cterm1","",Q3BINS,0,0.2);
  TH1D *Cterm1_noMRC = new TH1D("Cterm1_noMRC","",Q3BINS,0,0.2);
  TH1D *Cterm2 = new TH1D("Cterm2","",Q3BINS,0,0.2);
  TH1D *Cterm3 = new TH1D("Cterm3","",Q3BINS,0,0.2);
  TH1D *Cterm4 = new TH1D("Cterm4","",Q3BINS,0,0.2);
  TH1D *num_QS = new TH1D("num_QS","",Q3BINS,0,0.2);
  TH1D *Combinatorics_1d = new TH1D("Combinatorics_1d","",Q3BINS,0,0.2);
  TH1D *Combinatorics_1d_noMRC = new TH1D("Combinatorics_1d_noMRC","",Q3BINS,0,0.2);
  TH1D *Coul_Riverside = new TH1D("Coul_Riverside","",Q3BINS,0,0.2);
  TH1D *Coul_GRiverside = new TH1D("Coul_GRiverside","",Q3BINS,0,0.2);
  TH1D *Coul_Omega0 = new TH1D("Coul_Omega0","",Q3BINS,0,0.2);
  TH1D *c3_hist = new TH1D("c3_hist","",Q3BINS,0,0.2);
  TH1D *c3_hist_STAR = new TH1D("c3_hist_STAR","",Q3BINS,0,0.2);
  //TProfile *MomRes_2d = new TProfile("MomRes_2d","",Q3BINS,0,0.2, 0,20.,"");
  TProfile *MomRes_3d_term1 = new TProfile("MomRes_3d_term1","",Q3BINS,0,0.2, 0,20.,"");
  TProfile *MomRes_3d_term2 = new TProfile("MomRes_3d_term2","",Q3BINS,0,0.2, 0,20.,"");
  TProfile *MomRes_3d_term5 = new TProfile("MomRes_3d_term5","",Q3BINS,0,0.2, 0,20.,"");
  TH1D *r3_num_Q3 = new TH1D("r3_num_Q3","",Q3BINS,0,0.2);
  TH1D *r3_den_Q3 = new TH1D("r3_den_Q3","",Q3BINS,0,0.2);
  r3_num_Q3->GetXaxis()->SetTitle("Q_{3} (GeV/c)");
  r3_num_Q3->GetYaxis()->SetTitle("r_{3}");
  r3_num_Q3->GetYaxis()->SetTitleOffset(1.3);
  r3_num_Q3->GetXaxis()->SetRangeUser(0,0.1);
  r3_num_Q3->SetMinimum(0);
  r3_num_Q3->SetMaximum(2.4);
  r3_den_Q3->SetLineColor(2);
  r3_num_Q3->SetMarkerStyle(20);
  Coul_Omega0->GetXaxis()->SetRangeUser(0,0.099);
  Coul_Omega0->GetXaxis()->SetLabelSize(0.04);
  Coul_Omega0->GetYaxis()->SetLabelSize(0.04);
  c3_hist_STAR->GetXaxis()->SetRangeUser(0,0.099);
  c3_hist_STAR->SetMinimum(0.8); c3_hist_STAR->SetMaximum(1.02);
  c3_hist_STAR->GetXaxis()->SetTitle("Q_{3} (GeV/c)");
  c3_hist_STAR->GetYaxis()->SetTitle("c_{3}*^{++-}");
  c3_hist_STAR->GetYaxis()->SetTitleOffset(1.6);
  if(SameCharge) {Cterm1_noMRC->Sumw2(); Combinatorics_1d_noMRC->Sumw2();}
  //
  double num_QS_e[Q3BINS]={0};
  double c3_e[Q3BINS]={0};
  double r3_e[Q3BINS]={0};

  // CB=Charge Bin. 0 means pi-, 1 means pi+
  int CB1=0, CB2=0, CB3=0;
  int CP12=1, CP13=1, CP23=1;
  if(CHARGE==-1) {
    if(SameCharge) {CB1=0; CB2=0; CB3=0;}
    else {CB1=0; CB2=0; CB3=1; CP12=+1, CP13=-1, CP23=-1;}
  }else {
    if(SameCharge) {CB1=1; CB2=1; CB3=1;}
    else {CB1=0; CB2=1; CB3=1; CP12=-1, CP13=-1, CP23=+1;}
  }
  
  cout<<"here"<<endl;
  // SC = species combination.  Always 0 meaning pi-pi-pi.
  int SCBin=0;
  //
  ReadCoulCorrections(RValue, bValue, 10);// switch to full kt range, 10.
  //ReadCoulCorrections(0, 5, 2, 10);// Change to Gaussian R=5 fm calculation (STAR method testing)
  TH1D *GenSignalExpected_num=new TH1D("GenSignalExpected_num","",20,0,0.2);
  TH1D *GenSignalExpected_den=new TH1D("GenSignalExpected_den","",20,0,0.2);
  //
  double value_num[2]={0}; 
  double value_num_e[2]={0};
  double value_den[2]={0};
  double N3_QS[2]={0};
  double N3_QS_e[2]={0};
  double OutTriplets=0, InTriplets=0;
  for(int i=2; i<=20; i++){// bin number
    //double Q12_m = (i-0.5)*(0.005);// geometric center
    double Q12 = AvgQinvSS[i-1]; if(!SameCharge && CHARGE==+1) {Q12 = AvgQinvOS[i-1];}// true center
    int Q12bin = int(Q12/0.005)+1;
    //
    for(int j=2; j<=20; j++){// bin number
      //double Q13_m = (j-0.5)*(0.005);// geometric center
      double Q13 = AvgQinvSS[j-1]; if(!SameCharge) {Q13 = AvgQinvOS[j-1];}// true center
      int Q13bin = int(Q13/0.005)+1;
      //
      for(int k=2; k<=20; k++){// bin number
	//double Q23_m = (k-0.5)*(0.005);// geometric center
	double Q23 = AvgQinvSS[k-1]; if(!SameCharge && CHARGE==-1) {Q23 = AvgQinvOS[k-1];}// true center
	int Q23bin = int(Q23/0.005)+1;
	//
	//if(fabs(i-j)>3 || fabs(i-k)>3 || fabs(j-k)>3) continue;// testing
	
	double Q3 = sqrt(pow(Q12,2) + pow(Q13,2) + pow(Q23,2));
       	int Q3bin = Cterm1->GetXaxis()->FindBin(Q3);
	
	if(Q12 < sqrt(pow(Q13,2)+pow(Q23,2) - 2*Q13*Q23)) continue;// not all configurations are possible
	if(Q12 > sqrt(pow(Q13,2)+pow(Q23,2) + 2*Q13*Q23)) continue;// not all configurations are possible
	
	
       	// point-source Coulomb correlation
	double G_12 = Gamov(CP12, Q12);
	double G_13 = Gamov(CP13, Q13);
	double G_23 = Gamov(CP23, Q23);
	// finite-source Coulomb+Strong correlation from Therminator.
	double K2_12 = CoulCorr2(CP12, Q12);
	double K2_13 = CoulCorr2(CP13, Q13);
	double K2_23 = CoulCorr2(CP23, Q23);
	double K3 = 1.;// 3-body Coulomb+Strong correlation
	if(GRS) {// GRS approach
	  if(SameCharge || CHARGE==-1) K3 = CoulCorrGRS(SameCharge, Q12, Q13, Q23);
	  else K3 = CoulCorrGRS(SameCharge, Q23, Q12, Q13);
	}else {
	  if(SameCharge || CHARGE==-1) K3 = CoulCorrOmega0(SameCharge, Q12, Q13, Q23);
	  else K3 = CoulCorrOmega0(SameCharge, Q23, Q12, Q13);
	}
	if(MCcase && !GeneratedSignal) { K2_12=1.; K2_13=1.; K2_23=1.; K3=1.;}
	if(K3==0) continue;
	if(GeneratedSignal) K3 = K2_12*K2_13*K2_23;// no interpolation for Generated signal

	double TERM1=Three_3d[CB1][CB2][CB3][SCBin][0]->GetBinContent(i,j,k);// N3 (3-pion yield per q12,q13,q23 cell). 3-pions from same-event
	double TERM2=Three_3d[CB1][CB2][CB3][SCBin][1]->GetBinContent(i,j,k);// N2*N1. 1 and 2 from same-event
	double TERM3=Three_3d[CB1][CB2][CB3][SCBin][2]->GetBinContent(i,j,k);// N2*N1. 1 and 3 from same-event
	double TERM4=Three_3d[CB1][CB2][CB3][SCBin][3]->GetBinContent(i,j,k);// N2*N1. 2 and 3 from same-event
	double TERM5=Three_3d[CB1][CB2][CB3][SCBin][4]->GetBinContent(i,j,k);// N1*N1*N1. All from different events (pure combinatorics)
	if(TERM1==0 && TERM2==0 && TERM3==0 && TERM4==0 && TERM5==0) continue;
	if(TERM1==0 || TERM2==0 || TERM3==0 || TERM4==0 || TERM5==0) continue;
	
	double muonCorr_C3=1.0, muonCorr_C2_12=1.0, muonCorr_C2_13=1.0, muonCorr_C2_23=1.0;
	double muonCorr_W12=1.0, muonCorr_W13=1.0, muonCorr_W23=1.0;
	if(MuonCorrection){
	  if(SameCharge){
	    muonCorr_C3 = C3muonCorrection->GetBinContent(Q3bin);
	    muonCorr_C2_12 = C2muonCorrection_SC->GetBinContent(Q12bin);
	    muonCorr_C2_13 = C2muonCorrection_SC->GetBinContent(Q13bin);
	    muonCorr_C2_23 = C2muonCorrection_SC->GetBinContent(Q23bin);
	    muonCorr_W12 = WeightmuonCorrection->GetBinContent(Q12bin);
	    muonCorr_W13 = WeightmuonCorrection->GetBinContent(Q13bin);
	    muonCorr_W23 = WeightmuonCorrection->GetBinContent(Q23bin);
	  }else{
	    muonCorr_C3 = C3muonCorrection->GetBinContent(Q3bin);
	    if(CHARGE==-1){// pi-pi-pi+
	      muonCorr_C2_12 = C2muonCorrection_SC->GetBinContent(Q12bin);
	      muonCorr_C2_13 = C2muonCorrection_MC->GetBinContent(Q13bin);
	      muonCorr_C2_23 = C2muonCorrection_MC->GetBinContent(Q23bin);
	    }else{// pi-pi+pi+
	      muonCorr_C2_12 = C2muonCorrection_MC->GetBinContent(Q12bin);
	      muonCorr_C2_13 = C2muonCorrection_MC->GetBinContent(Q13bin);
	      muonCorr_C2_23 = C2muonCorrection_SC->GetBinContent(Q23bin);
	    }
	  }
	}	

	if(Q3>0.08 && Q3<0.09){// just for testing
	  if(Q12>0.08 || Q13>0.08 || Q23>0.08) OutTriplets++;
	  else InTriplets++;
	}
	
	// apply momentum resolution correction
	if(!MCcase && !GeneratedSignal){
	  if(SameCharge){
	    // 3d momentum resolution corrections
	    TERM1 *= (MomRes3d[0][0]->GetBinContent(Q12bin, Q13bin, Q23bin)-1)*MRCShift+1;
	    TERM2 *= (MomRes3d[0][1]->GetBinContent(Q12bin, Q13bin, Q23bin)-1)*MRCShift+1;
	    TERM3 *= (MomRes3d[0][2]->GetBinContent(Q12bin, Q13bin, Q23bin)-1)*MRCShift+1;
	    TERM4 *= (MomRes3d[0][3]->GetBinContent(Q12bin, Q13bin, Q23bin)-1)*MRCShift+1;
	    TERM5 *= (MomRes3d[0][4]->GetBinContent(Q12bin, Q13bin, Q23bin)-1)*MRCShift+1;
	    MomRes_3d_term1->Fill(Q3, MomRes3d[0][0]->GetBinContent(Q12bin, Q13bin, Q23bin),TERM1);
	    MomRes_3d_term2->Fill(Q3, MomRes3d[0][1]->GetBinContent(Q12bin, Q13bin, Q23bin),TERM2);
	    MomRes_3d_term5->Fill(Q3, MomRes3d[0][4]->GetBinContent(Q12bin, Q13bin, Q23bin),TERM5);
	  }else{
	    if(CHARGE==-1){// pi-pi-pi+
	      TERM1 *= (MomRes3d[1][0]->GetBinContent(Q12bin, Q13bin, Q23bin)-1)*MRCShift+1;
	      TERM2 *= (MomRes3d[1][1]->GetBinContent(Q12bin, Q13bin, Q23bin)-1)*MRCShift+1;
	      TERM3 *= (MomRes3d[1][2]->GetBinContent(Q12bin, Q13bin, Q23bin)-1)*MRCShift+1;
	      TERM4 *= (MomRes3d[1][3]->GetBinContent(Q12bin, Q13bin, Q23bin)-1)*MRCShift+1;
	      TERM5 *= (MomRes3d[1][4]->GetBinContent(Q12bin, Q13bin, Q23bin)-1)*MRCShift+1;
	    }else {// pi-pi+pi+
	      TERM1 *= (MomRes3d[1][0]->GetBinContent(Q23bin, Q13bin, Q12bin)-1)*MRCShift+1;
	      TERM2 *= (MomRes3d[1][3]->GetBinContent(Q23bin, Q13bin, Q12bin)-1)*MRCShift+1;
	      TERM3 *= (MomRes3d[1][2]->GetBinContent(Q23bin, Q13bin, Q12bin)-1)*MRCShift+1;
	      TERM4 *= (MomRes3d[1][1]->GetBinContent(Q23bin, Q13bin, Q12bin)-1)*MRCShift+1;
	      TERM5 *= (MomRes3d[1][4]->GetBinContent(Q23bin, Q13bin, Q12bin)-1)*MRCShift+1;
	    }
	    MomRes_3d_term1->Fill(Q3, MomRes3d[1][0]->GetBinContent(Q12bin, Q13bin, Q23bin), TERM1);
	    MomRes_3d_term2->Fill(Q3, MomRes3d[1][1]->GetBinContent(Q12bin, Q13bin, Q23bin),TERM2);
	    MomRes_3d_term5->Fill(Q3, MomRes3d[1][4]->GetBinContent(Q12bin, Q13bin, Q23bin),TERM5);
	  }
	}
	//
	
	for(int LamType=0; LamType<2; LamType++){
	  
	  if(LamType==1) TwoFrac=TwoFracr3;// r3
	  else TwoFrac = 0.7;//OutputPar[1];// c3 (newly fixed to 0.7)
	  
	  if(FileSetting==9) TwoFrac=0.8;// for FB5and7overlap choose a fixed higher lambda

	  OneFrac=pow(TwoFrac,0.5); // 0.495 to bring baseline up
	  ThreeFrac=pow(TwoFrac,1.5);
	  
	  // finite-multiplicity method
	  //OneFrac = (1+sqrt(1 + 4*AvgN[Mbin]*TwoFrac*(AvgN[Mbin]-1)))/(2*AvgN[Mbin]);
	  //ThreeFrac = (OneFrac*AvgN[Mbin]*(OneFrac*AvgN[Mbin]-1)*(OneFrac*AvgN[Mbin]-2))/(AvgN[Mbin]*(AvgN[Mbin]-1)*(AvgN[Mbin]-2));
	  

	  // Purify.  Isolate pure 3-pion QS correlations using Lambda and K3 (removes lower order correlations)
	  N3_QS[LamType] = TERM1;
	  N3_QS[LamType] -= ( pow(1-OneFrac,3) + 3*OneFrac*pow(1-OneFrac,2) )*TERM5;
	  N3_QS[LamType] -= (1-OneFrac)*(TERM2 + TERM3 + TERM4 - 3*(1-TwoFrac)*TERM5);
	  N3_QS[LamType] /= ThreeFrac;
	  N3_QS[LamType] /= K3;
	  N3_QS[LamType] *= muonCorr_C3;
	  
	  if(LamType==0) num_QS->Fill(Q3, N3_QS[LamType]);
	  
	  // Isolate 3-pion cumulant
	  value_num[LamType] = N3_QS[LamType];
	  value_num[LamType] -= (TERM2 - (1-TwoFrac)*TERM5)/K2_12/TwoFrac * muonCorr_C2_12;
	  value_num[LamType] -= (TERM3 - (1-TwoFrac)*TERM5)/K2_13/TwoFrac * muonCorr_C2_13;
	  value_num[LamType] -= (TERM4 - (1-TwoFrac)*TERM5)/K2_23/TwoFrac * muonCorr_C2_23;
	  value_num[LamType] += 2*TERM5;
	  
	  //value_num[LamType] += 0.004*TERM5;
	  	  
	  // r3 denominator
	  if(LamType==1 && SameCharge) {
	    value_den[LamType] = PionNorm[CB1]->GetBinContent(i,j,k);// Raw r3 denominator
	    if(!MCcase){
	      // Momentum Resolution Correction Systematic variation. Only important when MRCShift != 1.0.
	      double denMRC1 = (C2ssRaw->GetBinContent(i)*MomRes_C2_pp[i-1] - TwoFrac*K2_12 - (1-TwoFrac))/(TwoFrac*K2_12);
	      denMRC1 *= (C2ssRaw->GetBinContent(j)*MomRes_C2_pp[j-1] - TwoFrac*K2_13 - (1-TwoFrac))/(TwoFrac*K2_13);
	      denMRC1 *= (C2ssRaw->GetBinContent(k)*MomRes_C2_pp[k-1] - TwoFrac*K2_23 - (1-TwoFrac))/(TwoFrac*K2_23);
	      double denMRC2 = (C2ssRaw->GetBinContent(i)*((MomRes_C2_pp[i-1]-1)*MRCShift+1) - TwoFrac*K2_12 - (1-TwoFrac))/(TwoFrac*K2_12);
	      denMRC2 *= (C2ssRaw->GetBinContent(j)*((MomRes_C2_pp[j-1]-1)*MRCShift+1) - TwoFrac*K2_13 - (1-TwoFrac))/(TwoFrac*K2_13);
	      denMRC2 *= (C2ssRaw->GetBinContent(k)*((MomRes_C2_pp[k-1]-1)*MRCShift+1) - TwoFrac*K2_23 - (1-TwoFrac))/(TwoFrac*K2_23);
	      // Non-femto background correction (Mini-Jet).
	      double denMJ = (C2ssRaw->GetBinContent(i)*MomRes_C2_pp[i-1] - (MJ_bkg_ss->Eval(Q12)-1) - TwoFrac*K2_12 - (1-TwoFrac))/(TwoFrac*K2_12);
	      denMJ *= (C2ssRaw->GetBinContent(j)*MomRes_C2_pp[j-1] - (MJ_bkg_ss->Eval(Q13)-1) - TwoFrac*K2_13 - (1-TwoFrac))/(TwoFrac*K2_13);
	      denMJ *= (C2ssRaw->GetBinContent(k)*MomRes_C2_pp[k-1] - (MJ_bkg_ss->Eval(Q23)-1) - TwoFrac*K2_23 - (1-TwoFrac))/(TwoFrac*K2_23);
	      // Strong same-charge correction
	      double Coul_12 = K2_12/StrongSC->Eval(Q12*1000/2.);// Coulomb used online
	      double Coul_13 = K2_13/StrongSC->Eval(Q13*1000/2.);// Coulomb used online
	      double Coul_23 = K2_23/StrongSC->Eval(Q23*1000/2.);// Coulomb used online
	      double denCoulOnly = (C2ssRaw->GetBinContent(i)*MomRes_C2_pp[i-1] - TwoFrac*Coul_12 - (1-TwoFrac))/(TwoFrac*Coul_12);
	      denCoulOnly *= (C2ssRaw->GetBinContent(j)*MomRes_C2_pp[j-1] - TwoFrac*Coul_13 - (1-TwoFrac))/(TwoFrac*Coul_13);
	      denCoulOnly *= (C2ssRaw->GetBinContent(k)*MomRes_C2_pp[k-1] - TwoFrac*Coul_23 - (1-TwoFrac))/(TwoFrac*Coul_23);
	      // K shift testing
	      double K2back_12 = (K2_12-1)/KShift+1;
	      double K2back_13 = (K2_13-1)/KShift+1;
	      double K2back_23 = (K2_23-1)/KShift+1;
	      double denKshiftBack = (C2ssRaw->GetBinContent(i)*MomRes_C2_pp[i-1] - TwoFrac*K2back_12 - (1-TwoFrac))/(TwoFrac*K2back_12);
	      denKshiftBack *= (C2ssRaw->GetBinContent(j)*MomRes_C2_pp[j-1] - TwoFrac*K2back_13 - (1-TwoFrac))/(TwoFrac*K2back_13);
	      denKshiftBack *= (C2ssRaw->GetBinContent(k)*MomRes_C2_pp[k-1] - TwoFrac*K2back_23 - (1-TwoFrac))/(TwoFrac*K2back_23);
	      
	      double den_ratio = sqrt(fabs(denMRC2))/sqrt(fabs(denMRC1));
	      value_den[LamType] *= den_ratio;// apply Momentum Resolution correction systematic variation if any.
	      double den_ratioMJ = sqrt(fabs(denMJ))/sqrt(fabs(denMRC1));
	      if(IncludeMJcorrection) value_den[LamType] *= den_ratioMJ;// Non-femto bkg correction
	      if(IncludeStrongSC){
		double den_ratioStrong = sqrt(fabs(denMRC1))/sqrt(fabs(denCoulOnly));
		value_den[LamType] *= den_ratioStrong;
	      }
	      double den_ratioKShift = sqrt(fabs(denMRC1))/sqrt(fabs(denKshiftBack));
	      value_den[LamType] *= den_ratioKShift;
	      //
	      value_den[LamType] *= sqrt(muonCorr_W12*muonCorr_W13*muonCorr_W23);// muon correction
	      //
	      //value_den[LamType] -= TPNRejects->GetBinContent(i,j,k);// Testing for 0-5% only to estimate the effect when C2^QS < 1.0
	      value_den[LamType] *= (MomRes3d[0][4]->GetBinContent(Q12bin, Q13bin, Q23bin)-1)*MRCShift+1;// momentum resolution correction to combinatorics
	    }// !MCcase
	  }
	 

	  // errors
	  N3_QS_e[LamType] = fabs(TERM1);
	  N3_QS_e[LamType] += pow(( pow(1-OneFrac,3) +3*OneFrac*pow(1-OneFrac,2) )*sqrt(fabs(TERM5)),2);
	  N3_QS_e[LamType] += (pow((1-OneFrac),2)*fabs(TERM2 + TERM3 + TERM4) + pow((1-OneFrac)*3*(1-TwoFrac)*sqrt(fabs(TERM5)),2));
	  N3_QS_e[LamType] /= pow(ThreeFrac,2);
	  N3_QS_e[LamType] /= pow(K3,2);
	  if(LamType==0) num_QS_e[Q3bin-1] += N3_QS_e[LamType];
	  value_num_e[LamType] = N3_QS_e[LamType];
	  value_num_e[LamType] += (pow(1/K2_12/TwoFrac*sqrt(fabs(TERM2)),2) + pow((1-TwoFrac)/K2_12/TwoFrac*sqrt(fabs(TERM5)),2));
	  value_num_e[LamType] += (pow(1/K2_13/TwoFrac*sqrt(fabs(TERM3)),2) + pow((1-TwoFrac)/K2_13/TwoFrac*sqrt(fabs(TERM5)),2));
	  value_num_e[LamType] += (pow(1/K2_23/TwoFrac*sqrt(fabs(TERM4)),2) + pow((1-TwoFrac)/K2_23/TwoFrac*sqrt(fabs(TERM5)),2));
	  value_num_e[LamType] += pow(2*sqrt(fabs(TERM5)),2);
	  if(LamType==0) c3_e[Q3bin-1] += value_num_e[LamType] + TERM5;// add baseline stat error
	  else r3_e[Q3bin-1] += value_num_e[LamType];
	  
	}

	// Fill histograms
	c3_hist->Fill(Q3, value_num[0] + TERM5);// for cumulant correlation function
	C3QS_3d->SetBinContent(i,j,k, N3_QS[0]);
	C3QS_3d->SetBinError(i,j,k, N3_QS_e[0]);
	//
	Coul_Riverside->Fill(Q3, G_12*G_13*G_23*TERM5);
	Coul_GRiverside->Fill(Q3, TERM5*CoulCorrGRS(SameCharge, Q12, Q13, Q23));
	Coul_Omega0->Fill(Q3, K3*TERM5);
	num_withGRS->Fill(Q3, N3_QS[0]);
	Cterm1->Fill(Q3, TERM1);
	Cterm2->Fill(Q3, TERM2);
	Cterm3->Fill(Q3, TERM3);
	Cterm4->Fill(Q3, TERM4);
	Combinatorics_1d->Fill(Q3, TERM5);
	Combinatorics_3d->Fill(Q12,Q13,Q23, TERM5);
	r3_num_Q3->Fill(Q3, value_num[1]);
	r3_den_Q3->Fill(Q3, value_den[1]);
	double Q3_m = sqrt(pow((i-0.5)*(0.005),2) + pow((j-0.5)*(0.005),2) + pow((k-0.5)*(0.005),2));
	Cterm1_noMRC->Fill(Q3_m, Three_3d[CB1][CB2][CB3][SCBin][0]->GetBinContent(i,j,k));
	Combinatorics_1d_noMRC->Fill(Q3_m, Three_3d[CB1][CB2][CB3][SCBin][4]->GetBinContent(i,j,k));
	double cumulant_STAR = Three_3d[CB1][CB2][CB3][SCBin][0]->GetBinContent(i,j,k)/(K2_12*K2_13*K2_23);
	cumulant_STAR -= Three_3d[CB1][CB2][CB3][SCBin][1]->GetBinContent(i,j,k)/(K2_12);
	cumulant_STAR -= Three_3d[CB1][CB2][CB3][SCBin][2]->GetBinContent(i,j,k)/(K2_13);
	cumulant_STAR -= Three_3d[CB1][CB2][CB3][SCBin][3]->GetBinContent(i,j,k)/(K2_23);
	c3_hist_STAR->Fill(Q3_m, cumulant_STAR + 3*Three_3d[CB1][CB2][CB3][SCBin][4]->GetBinContent(i,j,k));
	
	///////////////////////////////////////////////////////////
	// Edgeworth 3-pion expection.  Not important for r3.
	//double radius_exp = (11-(MomResCentBin-1))/FmToGeV;
	//TwoFrac=0.68; OneFrac=pow(TwoFrac,.5), ThreeFrac=pow(TwoFrac,1.5);
	double radius_exp = (OutputPar[3])/FmToGeV;
	TwoFrac=OutputPar[1]; OneFrac=pow(TwoFrac,.5), ThreeFrac=pow(TwoFrac,1.5);
	double arg12 = Q12*radius_exp;
	double arg13 = Q13*radius_exp;
	double arg23 = Q23*radius_exp;
	Float_t EW12 = 1 + 0.2/(6.*pow(2.,1.5))*(8.*pow(arg12,3) - 12.*arg12);
	EW12 += 0.45/(24.*pow(2.,2))*(16.*pow(arg12,4) -48.*pow(arg12,2) + 12);
	Float_t EW13 = 1 + 0.2/(6.*pow(2.,1.5))*(8.*pow(arg13,3) - 12.*arg13);
	EW13 += 0.45/(24.*pow(2.,2))*(16.*pow(arg13,4) -48.*pow(arg13,2) + 12);
	Float_t EW23 = 1 + 0.2/(6.*pow(2.,1.5))*(8.*pow(arg23,3) - 12.*arg23);
	EW23 += 0.45/(24.*pow(2.,2))*(16.*pow(arg23,4) -48.*pow(arg23,2) + 12);
	
	//Float_t EW12 = 1 + OutputPar[5]/(6.*pow(2.,1.5))*(8.*pow(arg12,3) - 12.*arg12);
	//EW12 += OutputPar[6]/(24.*pow(2.,2))*(16.*pow(arg12,4) -48.*pow(arg12,2) + 12);
	//Float_t EW13 = 1 + OutputPar[5]/(6.*pow(2.,1.5))*(8.*pow(arg13,3) - 12.*arg13);
	//EW13 += OutputPar[6]/(24.*pow(2.,2))*(16.*pow(arg13,4) -48.*pow(arg13,2) + 12);
	//Float_t EW23 = 1 + OutputPar[5]/(6.*pow(2.,1.5))*(8.*pow(arg23,3) - 12.*arg23);
	//EW23 += OutputPar[6]/(24.*pow(2.,2))*(16.*pow(arg23,4) -48.*pow(arg23,2) + 12);
	//
	double cumulant_exp=0, term1_exp=0;
	if(SameCharge) {
	  cumulant_exp = (exp(-pow(radius_exp*Q12,2))*pow(EW12,2)+exp(-pow(radius_exp*Q13,2))*pow(EW13,2)+exp(-pow(radius_exp*Q23,2))*pow(EW23,2))*TERM5;
	  cumulant_exp += 2*EW12*EW13*EW23*exp(-pow(radius_exp,2)/2. * (pow(Q12,2)+pow(Q13,2)+pow(Q23,2)))*TERM5 + TERM5;
	  term1_exp = ( pow(1-OneFrac,3) + 3*OneFrac*pow(1-OneFrac,2) )*TERM5;
	  term1_exp += TwoFrac*(1-OneFrac)*(K2_12*(1+exp(-pow(radius_exp*Q12,2))*pow(EW12,2))+K2_13*(1+exp(-pow(radius_exp*Q13,2))*pow(EW13,2))+K2_23*(1+exp(-pow(radius_exp*Q23,2))*pow(EW23,2)))*TERM5;
	  term1_exp += ThreeFrac*K3*cumulant_exp;
	  //term1_exp = ((1-TwoFrac) + TwoFrac*K2_12*(1+exp(-pow(radius_exp*Q12,2))*pow(EW12,2)))*TERM5;
	}else {
	  cumulant_exp = (exp(-pow(radius_exp*Q12,2))*pow(EW12,2))*TERM5 + TERM5;
	  term1_exp = ( pow(1-OneFrac,3) + 3*OneFrac*pow(1-OneFrac,2) )*TERM5;
	  term1_exp += TwoFrac*(1-OneFrac)*(K2_12*(1+exp(-pow(radius_exp*Q12,2))*pow(EW12,2)) + K2_13 + K2_23)*TERM5;
	  term1_exp += ThreeFrac*K3*cumulant_exp;
	  term1_exp = ((1-TwoFrac) + TwoFrac*K2_12*(1+exp(-pow(radius_exp*Q12,2))*pow(EW12,2)))*TERM5;
	  //term1_exp = ((1-TwoFrac) + TwoFrac*K2_13)*TERM5;
	}
	
	GenSignalExpected_num->Fill(Q3, term1_exp);
	GenSignalExpected_den->Fill(Q3, TERM5);
	///////////////////////////////////////////////////////////
	
      }
    }
  }
  
  cout<<"OutTriplets: "<<OutTriplets<<"   InTriplets: "<<InTriplets<<endl;
  ////////////////////////////

  // Intermediate steps
  num_withRS->Sumw2();
  num_withGRS->Sumw2();
  num_withTM->Sumw2();
  Cterm1->Sumw2();
  Cterm2->Sumw2();
  Cterm3->Sumw2();
  Cterm4->Sumw2();
  Combinatorics_1d->Sumw2();
  Combinatorics_3d->Sumw2();
  r3_num_Q3->Sumw2();
  r3_den_Q3->Sumw2();
  for(int i=0; i<Q3BINS; i++) {c3_hist->SetBinError(i+1, sqrt(c3_e[i])); num_QS->SetBinError(i+1, sqrt(num_QS_e[i]));}
  for(int i=0; i<Q3BINS; i++) {c3_hist_STAR->SetBinError(i+1, sqrt(c3_e[i]));}// approximate error
  num_withRS->Divide(Combinatorics_1d);
  num_withGRS->Divide(Combinatorics_1d);
  num_withTM->Divide(Combinatorics_1d);
  for(int q3bin=1; q3bin<=r3_num_Q3->GetNbinsX(); q3bin++){
    r3_num_Q3->SetBinError(q3bin, sqrt(r3_e[q3bin-1]));
  }
  Cterm1->Divide(Combinatorics_1d);
  Cterm1_noMRC->Divide(Combinatorics_1d_noMRC);
  Cterm2->Divide(Combinatorics_1d);
  Cterm3->Divide(Combinatorics_1d);
  Cterm4->Divide(Combinatorics_1d);
  c3_hist->Divide(Combinatorics_1d);
  c3_hist_STAR->Divide(Combinatorics_1d_noMRC);
  num_QS->Divide(Combinatorics_1d);
  r3_num_Q3->Divide(r3_den_Q3);
  GenSignalExpected_num->Sumw2();
  GenSignalExpected_num->Divide(GenSignalExpected_den);
  
  //
  //
  Coul_Riverside->Divide(Combinatorics_1d);
  Coul_GRiverside->Divide(Combinatorics_1d);
  Coul_Omega0->Divide(Combinatorics_1d);
  for(int ii=1; ii<=Coul_Omega0->GetNbinsX(); ii++){
    Coul_Riverside->SetBinError(ii,0.000001);
    Coul_GRiverside->SetBinError(ii,0.000001);
    Coul_Omega0->SetBinError(ii,0.000001);
  }
  ////////////////////////////
 
  
 
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  
  Coul_Riverside->SetLineColor(1);
  Coul_GRiverside->SetLineColor(2);
  Coul_Omega0->SetLineColor(4);

  //Coul_Riverside->Draw();
  //legend1->AddEntry(Coul_Riverside,"Riverside","l");
  //Coul_GRiverside->Draw("same");
  //legend1->AddEntry(Coul_GRiverside,"Generalized Riverside","l");
  //Coul_Omega0->Draw("same");
  //legend1->AddEntry(Coul_Omega0,"Omega0","l");
 
  
  Cterm1->SetMarkerStyle(20);
  Cterm1->SetMarkerColor(4);
  Cterm1->SetLineColor(4);
  Cterm1->GetXaxis()->SetTitle("Q_{3} (GeV/c)");
  Cterm1->GetYaxis()->SetTitle("C_{3}");
  //Cterm1->SetTitle("#pi^{-}#pi^{-}#pi^{-}");
  Cterm1->SetMinimum(0.99);
  Cterm1->SetMaximum(1.08);// 6.1 for same-charge
  Cterm1->GetXaxis()->SetRangeUser(0,.095);
  Cterm1->GetYaxis()->SetTitleOffset(1.4);
  Cterm1->Draw();
  legend2->AddEntry(Cterm1,"C_{3} raw","p");
  //
  Cterm2->SetMarkerStyle(20);
  Cterm2->SetMarkerColor(7);
  
  Cterm3->SetMarkerStyle(25);
  Cterm3->SetMarkerColor(8);
  Cterm4->SetMarkerStyle(26);
  Cterm4->SetMarkerColor(9);
  //Cterm2->Draw("same");
  //Cterm3->Draw("same");
  //Cterm4->Draw("same");
  //legend2->AddEntry(Cterm1,"++-","p");
  
  
  if(MCcase){
    double C3points_HIJING_mmm[10]={0, 0.961834, 1.01827, 0.999387, 1.00202, 1.00081, 1.00082, 1.00037, 0.999074, 0.999099};
    double C3points_HIJING_mmm_e[10]={0, 0.0833895, 0.015158, 0.0047978, 0.00235067, 0.00138155, 0.00087485, 0.000612203, 0.000450162, 0.000346943};
    double C3points_HIJING_ppp[10]={0, 1.13015, 1.00623, 1.00536, 1.00293, 0.999964, 1.00116, 1.0007, 1.00037, 1.00105};
    double C3points_HIJING_ppp_e[10]={0, 0.0977058, 0.0150694, 0.0048196, 0.00235518, 0.00138376, 0.000877562, 0.000614069, 0.000452051, 0.00034856};
    TH1D *C3_HIJING_mmm=(TH1D*)Cterm1->Clone();
    TH1D *C3_HIJING_ppp=(TH1D*)Cterm1->Clone();
    for(int i=0; i<10; i++){
      C3_HIJING_mmm->SetBinContent(i+1, C3points_HIJING_mmm[i]);
      C3_HIJING_mmm->SetBinError(i+1, C3points_HIJING_mmm_e[i]);
      C3_HIJING_ppp->SetBinContent(i+1, C3points_HIJING_ppp[i]);
      C3_HIJING_ppp->SetBinError(i+1, C3points_HIJING_ppp_e[i]);
    }
    C3_HIJING_mmm->SetMarkerColor(2);
    C3_HIJING_mmm->SetLineColor(2);
    C3_HIJING_ppp->SetMarkerColor(4);
    C3_HIJING_ppp->SetLineColor(4);
    //C3_HIJING_mmm->Draw("same");
    //C3_HIJING_ppp->Draw("same");
    //legend2->AddEntry(C3_HIJING_mmm,"---","p");
    //legend2->AddEntry(C3_HIJING_ppp,"+++","p");
  }

  num_QS->SetMarkerStyle(24);
  num_QS->SetMarkerColor(4);
  num_QS->SetLineColor(4);
  num_QS->GetXaxis()->SetTitle("Q_{3}");
  num_QS->GetYaxis()->SetTitle("C_{3}^{QS}");
  num_QS->GetXaxis()->SetRangeUser(0,.12);
  num_QS->SetMaximum(6);
  //num_QS->Draw("same");
  //legend2->AddEntry(num_QS,"C_{3}^{QS}","p");
 
  c3_hist->GetXaxis()->SetTitle("Q_{3}");
  c3_hist->GetYaxis()->SetTitle("c_{3}^{++-}");
  c3_hist->GetXaxis()->SetRangeUser(0,.15);
  c3_hist->SetMarkerStyle(25);
  c3_hist->SetMarkerColor(2);
  c3_hist->SetLineColor(2);
  c3_hist->SetMaximum(3);
  c3_hist->SetMinimum(.9);
  if(!MCcase && !GeneratedSignal) c3_hist->Draw("same");
  //legend2->AddEntry(c3_hist,"#font[12]{c}_{3} (cumulant correlation)","p");
  
  if(SameCharge){
    //for(int ii=0; ii<10; ii++) cout<<c3_hist->GetBinContent(ii+1)<<", ";
    TF1 *c3Fit=new TF1("c3Fit","[0]*(1+[1]*exp(-pow([2]*x/0.19733,2)/2.))",0,0.2);
    //TF1 *c3Fit=new TF1("c3Fit","[0]*(1+[1]*exp(-pow([2]*x/0.19733,2)/2.) * (1 + (-0.12/(6.*pow(2,1.5))*(8.*pow([2]*x/0.19733,3) - 12.*pow([2]*x/0.19733,1))) + (0.43/(24.*pow(2,2))*(16.*pow([2]*x/0.19733,4) -48.*pow([2]*x/0.19733,2) + 12))))",0,1);
    c3Fit->SetParameter(0,1);
    c3Fit->SetParameter(1,2);
    c3Fit->SetParameter(2,6);
    //c3Fit->FixParameter(1,0.8545);
    //c3Fit->FixParameter(2,8.982);
    //c3_hist->Fit(c3Fit,"IME","",0,0.11);
    //cout<<"Chi2/NDF = "<<c3Fit->GetChisquare()/c3Fit->GetNDF()<<endl;
  }

  GenSignalExpected_num->SetMarkerStyle(20);
  //GenSignalExpected_num->Draw("same");
  //legend2->AddEntry(GenSignalExpected_num,"#kappa_{3}=0.24, #kappa_{4}=0.16, #lambda=0.68, R=6 fm","p");
  //legend2->AddEntry(GenSignalExpected_num,"Edeworth expectation (fully chaotic)","p");

  //MomRes_2d->SetMarkerStyle(20);
  //MomRes_3d->SetMarkerStyle(20);
  //MomRes_3d->SetMarkerColor(4);
  //MomRes_2d->GetXaxis()->SetTitle("Q_{3} (GeV/c)");
  //MomRes_2d->GetYaxis()->SetTitle("< MRC >");
  //MomRes_2d->SetTitle("");
  //MomRes_2d->Draw();
  //legend2->AddEntry(MomRes_2d,"2D: Triple pair product","p");
  MomRes_3d_term2->SetLineColor(2);
  MomRes_3d_term5->SetLineColor(4);
  //MomRes_3d_term1->Draw();
  //MomRes_3d_term2->Draw("same");
  //MomRes_3d_term5->Draw("same");
  //legend2->AddEntry(MomRes_3d,"3D","p");
  
  //legend2->Draw("same");
  
  //cout<<c3_hist->Integral(8,10)/3.<<endl;
  
  /*
  ////////// C3 systematics
  // C3 --+ base, M0, Kt3 0
  double C3_base[10]={0, 1.64715, 1.40709, 1.24344, 1.15809, 1.1071, 1.07544, 1.0534, 1.03881, 1.02974};
  double C3_base_e[10]={0, 0.00348782, 0.000841611, 0.00031699, 0.000163779, 9.95652e-05, 6.43811e-05, 4.60347e-05, 3.45978e-05, 2.68396e-05};
  //
  // HIJING C3 --- base, M0
  //double C3_base[10]={0, 0.970005, 1.00655, 1.00352, 1.00291, 1.00071, 1.0002, 0.999524, 0.999404, 0.999397};
  //double C3_base_e[10]={0, 0.050845, 0.0099602, 0.00334862, 0.00138008, 0.000841743, 0.000531776, 0.000371712, 0.000274716, 0.00021};
  // HIJING C3 --- base, M1
  //double C3_base[10]={0, 1.03657, 1.00199, 0.997984, 1.00067, 1.0006, 0.999901, 0.999967, 0.999792, 0.999777};
  //double C3_base_e[10]={0, 0.0634232, 0.0117204, 0.0039446, 0.00163131, 0.000996638, 0.000629662, 0.000440266, 0.00032534, 0.000249};
  // HIJING C3 --- base, M2
  //double C3_base[10]={0, 1.34345, 1.04226, 1.0278, 0.99582, 1.00554, 1.00296, 1.00057, 1.00271, 1.00152};
  //double C3_base_e[10]={0, 0.363559, 0.0503531, 0.0170117, 0.00679185, 0.00419035, 0.00264603, 0.00184587, 0.00136663, 0.00104772};
  // HIJING C3 --- base, M3
  //double C3_base[10]={0, 0.890897, 1.05222, 1.00461, 1.01364, 0.998981, 1.00225, 1.00305, 1.00235, 1.00043};
  //double C3_base_e[10]={0, 0.297124, 0.0604397, 0.0195066, 0.00812906, 0.00490835, 0.00310751, 0.00217408, 0.00160575, 0.00123065};
  TH1D *C3baseHisto=(TH1D*)Cterm1->Clone();
  for(int ii=0; ii<10; ii++){
    C3baseHisto->SetBinContent(ii+1, C3_base[ii]);
    C3baseHisto->SetBinError(ii+1, C3_base_e[ii]);
  }

  cout<<"C3 values:"<<endl;
  for(int ii=0; ii<10; ii++){
    cout<<Cterm1->GetBinContent(ii+1)<<", ";
  }
  cout<<endl;
  cout<<"C3 errors:"<<endl;
  for(int ii=0; ii<10; ii++){
    cout<<Cterm1->GetBinError(ii+1)<<", ";
  }
  cout<<endl;
  TH1D *C3_Sys=(TH1D*)Cterm1->Clone();
  for(int ii=1; ii<10; ii++){
    if(C3_base[ii] ==0) {
      C3_Sys->SetBinContent(ii+1, 100.);
      C3_Sys->SetBinError(ii+1, 100.);
      continue;
    }
    C3_Sys->SetBinContent(ii+1, (C3_base[ii]-C3_Sys->GetBinContent(ii+1))/C3_base[ii]);
    C3_Sys->SetBinError(ii+1, sqrt(fabs(pow(C3_Sys->GetBinError(ii+1),2) - pow(C3_base_e[ii],2))));
  }
  gStyle->SetOptFit(111);
  C3_Sys->GetXaxis()->SetRangeUser(0,0.099);
  C3_Sys->GetYaxis()->SetTitle("(C_{3}^{t1}-C_{3}^{t2})/C_{3}^{t1}");
  C3_Sys->GetYaxis()->SetTitleOffset(2);
  C3_Sys->SetMinimum(-0.01);
  C3_Sys->SetMaximum(0.01);
  //C3_Sys->Draw();
  TF1 *C3lineSys=new TF1("C3lineSys","pol3",0,0.099);
  C3_Sys->Fit(C3lineSys,"MEQ","",0,0.099);
  
  if(MCcase){
    // Plotting +++ added with --- for HIJING
    TLegend *legendC3Hijing = new TLegend(.5,.15,.9,.3,NULL,"brNDC");
    legendC3Hijing->SetBorderSize(1);
    legendC3Hijing->SetTextSize(.03);// small .03; large .06
    legendC3Hijing->SetFillColor(0);
    
    Cterm1->Add(C3baseHisto);
    Cterm1->Scale(0.5);
    Cterm1->SetMinimum(0.9);
    Cterm1->SetMaximum(1.1);
    Cterm1->Draw();
    legendC3Hijing->AddEntry(Cterm1,"same-charge, HIJING","p");
    legendC3Hijing->Draw("same");
  }
  */
  /*
  ////////// c3 systematics
  // c3 --- base, M0, (0.04 TTC )
  //double c3_base[10]={0, 1.86014, 1.47533, 1.23733, 1.09944, 1.04145, 1.01693, 1.00715, 1.00253, 1.00111};
  //double c3_base_e[10]={0, 0.104645, 0.0120917, 0.00333303, 0.00118126, 0.0006692, 0.000405246, 0.000274163, 0.000198507, 0.000150258};
  // c3 --- base, M0, Kt3 0 (New Standard)
  double c3_base[10]={0, 1.00906, 1.00013, 1.00203, 1.00001, 1.00017, 1.00001, 0.999721, 0.999778, 0.999844};
  double c3_base_e[10]={0, 0.00726286, 0.00196376, 0.00081047, 0.00044086, 0.000276571, 0.000182574, 0.000132467, 0.000100557, 7.85009e-05};

  cout<<"c3 values:"<<endl;
  for(int ii=0; ii<10; ii++){
    cout<<c3_hist->GetBinContent(ii+1)<<", ";
  }
  cout<<endl;
  cout<<"c3 errors:"<<endl;
  for(int ii=0; ii<10; ii++){
    cout<<c3_hist->GetBinError(ii+1)<<", ";
  }
  cout<<endl;
  TH1D *c3_Sys=(TH1D*)c3_hist->Clone();
  for(int ii=1; ii<10; ii++){
    if(c3_base[ii] ==0) {
      c3_Sys->SetBinContent(ii+1, 100.);
      c3_Sys->SetBinError(ii+1, 100.);
      continue;
    }
    c3_Sys->SetBinContent(ii+1, (c3_base[ii]-c3_Sys->GetBinContent(ii+1))/c3_base[ii]);
    c3_Sys->SetBinError(ii+1, sqrt(fabs(pow(c3_Sys->GetBinError(ii+1),2) - pow(c3_base_e[ii],2))));
  }
  gStyle->SetOptFit(111);
  c3_Sys->GetXaxis()->SetRangeUser(0,0.099);
  c3_Sys->GetYaxis()->SetTitle("(C_{3}^{t1}-C_{3}^{t2})/C_{3}^{t1}");
  c3_Sys->GetYaxis()->SetTitleOffset(2);
  c3_Sys->SetMinimum(-0.01);
  c3_Sys->SetMaximum(0.01);
  c3_Sys->Draw();
  TF1 *c3lineSys=new TF1("c3lineSys","pol3",0,0.099);
  c3_Sys->Fit(c3lineSys,"MEQ","",0,0.099);
  */

  
  /*TPad *pad1 = new TPad("pad1","pad1",0.0,0.0,1.,1.);
  gPad->SetGridx(0);
  gPad->SetGridy(0);
  gPad->SetTickx();
  gPad->SetTicky();
  pad1->SetTopMargin(0.02);//0.05
  //pad1->SetLeftMargin(0.13);//.14 for wide title, .10 for narrow title, 0.08 for DeltaQ
  pad1->SetRightMargin(0.01);//1e-2
  pad1->SetBottomMargin(0.06);//0.12
  pad1->Divide(1,2,0,0);
  pad1->Draw();
  pad1->cd(1);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.02);
  if(SameCharge){
    Coul_Omega0->SetMinimum(0.48);
    Coul_Omega0->SetMaximum(1.01);
  }else{
    Coul_Omega0->SetMinimum(0.99);
    Coul_Omega0->SetMaximum(1.32);
  }
  Coul_Omega0->SetMarkerStyle(20);
  Coul_Omega0->SetMarkerColor(4);
  Coul_GRiverside->SetMarkerStyle(25);
  Coul_GRiverside->SetMarkerColor(2);
  Coul_Riverside->SetMarkerStyle(23);
  Coul_Omega0->Draw();
  //Coul_Riverside->Draw("same");
  Coul_GRiverside->Draw("same");
  legend2->AddEntry(Coul_Omega0,"#Omega_{0}","p");
  legend2->AddEntry(Coul_GRiverside,"Generalized Riverside","p");
  //legend2->AddEntry(Coul_Riverside,"Riverside","p");
  legend2->Draw("same");
  TLatex *K3Label = new TLatex(-0.011,0.92,"K_{3}");// -0.011,0.92 (ss), 1.26 (os)
  K3Label->SetTextSize(0.08);
  K3Label->SetTextAngle(90);
  K3Label->Draw();
  //
  pad1->cd(2);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.02);
  gPad->SetBottomMargin(0.13);
  TH1D *K3_compOmega0 = (TH1D*)Coul_Omega0->Clone();
  TH1D *K3_compGRS = (TH1D*)Coul_GRiverside->Clone();
  for(int ii=0; ii<K3_compOmega0->GetNbinsX(); ii++){ 
    K3_compOmega0->SetBinContent(ii+1, K3_compOmega0->GetBinContent(ii+1)-1.0);
    K3_compGRS->SetBinContent(ii+1, K3_compGRS->GetBinContent(ii+1)-1.0);
  }
  K3_compOmega0->Add(K3_compGRS,-1);
  K3_compOmega0->Divide(K3_compGRS);
  K3_compOmega0->SetMarkerStyle(20);
  K3_compOmega0->SetMarkerColor(1);
  K3_compOmega0->SetLineColor(1);
  K3_compOmega0->SetMinimum(-0.12);// -0.021
  K3_compOmega0->SetMaximum(0.12);// 0.021
  K3_compOmega0->SetBinContent(1,-100);
  K3_compOmega0->Draw();
  TLatex *RatioLabel = new TLatex(-.011,-0.05,"#Delta (K_{3}-1)/(K_{3}(#Omega_{0})-1)");// -0.011,0.04
  RatioLabel->SetTextSize(0.08);
  RatioLabel->SetTextAngle(90);
  RatioLabel->Draw();
  TLatex *Q3Label = new TLatex(.065,-0.147,"Q_{3} (GeV/c)");//.065,-0.147
  Q3Label->SetTextSize(0.08);
  Q3Label->Draw();
  */


  /////////////////////////////
  // 3d C3QS fitting
  /*
  for(int i=0; i<BINRANGE_C2global; i++){
    for(int j=0; j<BINRANGE_C2global; j++){
      for(int k=0; k<BINRANGE_C2global; k++){
	C3ssFitting[i][j][k] = C3QS_3d->GetBinContent(i+1,j+1,k+1);
	C3ssFitting_e[i][j][k] = C3QS_3d->GetBinError(i+1,j+1,k+1);
      }
    }
  }
  
  TF3 *C3ssFitResult=new TF3("C3ssFitResult",C3ssFitFunction,0,0.2,0,0.2,0,0.2,npar);
  for(int i=0; i<npar; i++) {
    if(i<=5) C3ssFitResult->FixParameter(i, OutputPar[i]);
    else C3ssFitResult->FixParameter(i, OutputPar[i]-OutputPar_e[i]);
  }
  
  
  const int NbinsC3_1d=Q3BINS;
  TH1D *C3fit_1d=new TH1D("C3fit_1d","",num_QS->GetNbinsX(),num_QS->GetXaxis()->GetBinLowEdge(1),num_QS->GetXaxis()->GetBinUpEdge(num_QS->GetNbinsX()));
  TH1D *C3fit_1d_den=new TH1D("C3fit_1d_den","",num_QS->GetNbinsX(),num_QS->GetXaxis()->GetBinLowEdge(1),num_QS->GetXaxis()->GetBinUpEdge(num_QS->GetNbinsX()));
  double C3_err[NbinsC3_1d]={0};
  
     
  for(int i=0; i<C3QS_3d->GetNbinsX(); i++){// bin number
    for(int j=0; j<C3QS_3d->GetNbinsY(); j++){// bin number
      for(int k=0; k<C3QS_3d->GetNbinsZ(); k++){// bin number
	
	if(i==0 || j==0 || k==0) continue;
	//if(i==1 || j==1 || k==1) continue;
	double Q3 = sqrt(pow(AvgQinvSS[i],2) + pow(AvgQinvSS[j],2) + pow(AvgQinvSS[k],2));
	C3fit_1d->Fill(Q3, C3ssFitResult->Eval(AvgQinvSS[i],AvgQinvSS[j],AvgQinvSS[k])*Combinatorics_3d->GetBinContent(i+1,j+1,k+1));
	C3fit_1d_den->Fill(Q3, Combinatorics_3d->GetBinContent(i+1,j+1,k+1));
      }
    }
  }
  
  C3fit_1d->Divide(C3fit_1d_den);
  for(int q3bin=1; q3bin<=num_QS->GetNbinsX(); q3bin++) C3fit_1d->SetBinError(q3bin, 0.001);// non-zero to display marker
  C3fit_1d->SetMarkerStyle(22);
  C3fit_1d->SetLineColor(4);
  //C3fit_1d->Draw("same");
 
  TF1 *lineC3 = new TF1("lineC3","6",0,10);
  lineC3->SetLineColor(4);
  lineC3->SetLineStyle(2);
  //lineC3->Draw("same");
  //legend2->AddEntry(lineC3,"Chaotic Limit C_{3}^{QS}","l");
  
  TF1 *linec3 = new TF1("linec3","3",0,10);
  linec3->SetLineColor(2);
  linec3->SetLineStyle(2);
  //linec3->Draw("same");
  //legend2->AddEntry(linec3,"Chaotic Limit #font[12]{c}_{3}","l"); 
  */
  
  //legend2->Draw("same");
 

  /////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////
  
  TCanvas *can3 = new TCanvas("can3", "can3",1600,0,700,700);
  can3->SetHighLightColor(2);
  can3->Range(-0.7838115,-1.033258,0.7838115,1.033258);
  //gStyle->SetOptFit(0111);
  can3->SetFillColor(10);
  can3->SetBorderMode(0);
  can3->SetBorderSize(2);
  can3->SetGridx();
  can3->SetGridy();
  can3->SetFrameFillColor(0);
  can3->SetFrameBorderMode(0);
  can3->SetFrameBorderMode(0);
  can3->cd();
  TPad *pad3 = new TPad("pad3","pad3",0.0,0.0,1.,1.);
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  gPad->SetTickx();
  gPad->SetTicky();
  pad3->SetTopMargin(0.02);//0.05
  //pad3->SetLeftMargin(0.13);//.14 for wide title, .10 for narrow title, 0.08 for DeltaQ
  pad3->SetRightMargin(0.01);//1e-2
  pad3->SetBottomMargin(0.06);//0.12
  pad3->Draw();
  pad3->cd(1);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.02);

  TF1 *QuarticFit = new TF1("QuarticFit","[0]*(1-[1]*pow(x,4))",0,.1);
  QuarticFit->SetParameter(0,2); QuarticFit->SetParameter(1,0);
  QuarticFit->SetLineColor(4);
  QuarticFit->SetParName(0,"I^{Quartic}"); QuarticFit->SetParName(1,"a^{Quartic}"); 
  TF1 *QuadraticFit = new TF1("QuadraticFit","[0]*(1-[1]*pow(x,2))",0,.1);
  QuadraticFit->SetParameter(0,2); QuadraticFit->SetParameter(1,0);
  QuadraticFit->SetParName(0,"I^{Quadratic}"); QuadraticFit->SetParName(1,"a^{Quadratic}"); 
  TLegend *legend3 = new TLegend(.2,.85,.5,.95,NULL,"brNDC");
  legend3->SetBorderSize(1);
  legend3->SetTextSize(.04);// small .03; large .06
  legend3->SetFillColor(0);


  if(SameCharge){
    
    r3_num_Q3->SetMinimum(0); r3_num_Q3->SetMaximum(2.5);// 0 to 2.5
    r3_num_Q3->GetXaxis()->SetRangeUser(0.0,0.1);
    r3_num_Q3->Draw();


    // HIJING standalone
    // Kt3 1
    /*double Gen_r3_m_M0[10]={0, 1.97963, 2.10248, 2.04465, 1.96697, 2.02295, 1.92281, 2.10031, 2.22338, 2.49729};
    double Gen_r3_m_M0_e[10]={0, 0.132726, 0.0668652, 0.0451347, 0.042838, 0.0546967, 0.0754362, 0.133624, 0.286307, 0.628381};
    double Gen_r3_p_M0[10]={0, 1.91771, 1.9653, 2.00742, 2.02393, 1.90624, 1.93554, 1.66, 1.79584, 0.301761};
    double Gen_r3_p_M0_e[10]={0, 0.133654, 0.0667213, 0.0451512, 0.0428925, 0.0547591, 0.0754764, 0.13365, 0.28638, 0.628441};*/
    // Kt3 2
    /*double Gen_r3_m_M0[10]={0, 1.12993, 2.09715, 1.91886, 2.08493, 2.10931, 2.00286, 1.99898, 1.78549, 1.91861};
    double Gen_r3_m_M0_e[10]={0, 0.237903, 0.115454, 0.0725178, 0.0630867, 0.0730326, 0.0886163, 0.12885, 0.213654, 0.379273};
    double Gen_r3_p_M0[10]={0, 2.05766, 1.97408, 2.05182, 2.02431, 2.11783, 1.93294, 1.97525, 2.21833, 2.31318};
    double Gen_r3_p_M0_e[10]={0, 0.24238, 0.11515, 0.0725077, 0.0629938, 0.0729842, 0.0886386, 0.12884, 0.213737, 0.379196};*/
    
    
    // HIJING + ALICE
    // Kt3 1
    /*double Gen_r3_m_M0[10]={0, 2.03468, 2.05783, 1.97757, 2.03809, 1.95703, 2.02915, 1.80055, 1.97664, 1.49573};
    double Gen_r3_m_M0_e[10]={0, 0.284164, 0.0923288, 0.0543837, 0.045952, 0.0565222, 0.0748221, 0.128994, 0.271954, 0.599077};
    double Gen_r3_p_M0[10]={0, 1.74611, 1.92369, 2.11024, 2.02823, 2.06235, 2.00127, 1.91551, 1.90576, 2.521};
    double Gen_r3_p_M0_e[10]={0, 0.279962, 0.0928866, 0.0548168, 0.0462253, 0.0569129, 0.0752683, 0.129803, 0.2736, 0.602535};
    double Gen_r3_m_M1[10]={0, 1.71341, 1.8735, 1.9352, 1.96466, 1.93852, 1.69509, 1.51402, 1.6014, -0.802941};
    double Gen_r3_m_M1_e[10]={0, 0.334562, 0.108368, 0.0640739, 0.0540579, 0.0663174, 0.0876532, 0.149302, 0.307982, 0.640548};
    double Gen_r3_p_M1[10]={0, 1.50861, 2.09508, 2.05338, 1.96178, 1.97999, 2.07162, 1.8607, 1.91805, 1.10468};
    double Gen_r3_p_M1_e[10]={0, 0.3399, 0.109449, 0.0646241, 0.0544419, 0.0667711, 0.0882695, 0.150195, 0.309844, 0.643923};*/
    // Kt3 2
    /*double Gen_r3_m_M0[10]={0, 1.04272, 2.05961, 1.93025, 2.08203, 2.07565, 2.29192, 1.93009, 2.68715, 1.71175};
    double Gen_r3_m_M0_e[10]={0, 3.09343, 0.30016, 0.117634, 0.0757208, 0.0794904, 0.090823, 0.128235, 0.21308, 0.40023};
    double Gen_r3_p_M0[10]={0, 1.83276, 1.88969, 2.03778, 2.0907, 2.11919, 2.04992, 2.02593, 1.95209, 1.68264};
    double Gen_r3_p_M0_e[10]={0, 3.59, 0.298056, 0.118354, 0.0762849, 0.079909, 0.0912179, 0.128866, 0.213999, 0.402173};
    double Gen_r3_m_M1[10]={0, 5.40628, 1.52822, 1.93258, 2.13338, 2.05811, 2.02963, 2.1204, 2.04906, 1.9021};
    double Gen_r3_m_M1_e[10]={0, 4.3025, 0.350163, 0.13871, 0.0897272, 0.0938973, 0.107045, 0.151557, 0.25029, 0.446325};
    double Gen_r3_p_M1[10]={0, 3.15883, 1.86772, 2.24914, 2.12118, 2.09175, 2.01447, 2.36802, 2.57239, 2.68729};
    double Gen_r3_p_M1_e[10]={0, 4.6622, 0.348237, 0.140294, 0.0903003, 0.09446, 0.107692, 0.152575, 0.251939, 0.448919};*/
    
    //double r3Sys_M1[10]={0, 0.097, 0.056, 0.063, 0.097, 0.17, 0.32, 0.66, 1.4, 3.4};
    // muon variation
    /*double r3_muonVar_M1[10]={0, 1.73244, 1.80921, 1.77852, 1.7192, 1.62059, 1.50122, 1.37656, 1.01344, 0.755781};
    double r3_muonVar_M1_e[10]={0, 0.160786, 0.051075, 0.031982, 0.0293467, 0.0390241, 0.0514494, 0.0760959, 0.112944, 0.154592};
    TH1D *r3_muonVar=(TH1D*)r3_num_Q3->Clone();
    for(int i=0; i<10; i++){ 
      r3_muonVar->SetBinContent(i+1,r3_muonVar_M1[i]); 
      r3_muonVar->SetBinError(i+1,sqrt(pow(r3_muonVar_M1_e[i],2)+pow(r3Sys_M1[i],2)));
    }
    r3_muonVar->SetMarkerStyle(21);
    r3_muonVar->SetMarkerColor(2); r3_muonVar->SetLineColor(2);
    r3_muonVar->SetFillStyle(1000); r3_muonVar->SetFillColor(kRed-10);
    r3_muonVar->Draw("E2 same");
    r3_num_Q3->Draw("same");
    legend3->AddEntry(r3_num_Q3,"-2.0<N#sigma_{Pion}<2.0 (default)","p");
    legend3->AddEntry(r3_muonVar,"-0.5<N#sigma_{Pion}<2.0","p");
    legend3->Draw("same");*/
    
    // muon correction 
    /*double r3_muonCorr_M1[10]={0, 1.8462, 1.84371, 1.79934, 1.77217, 1.76725, 1.79545, 1.7986, 2.11717, 2.86177};
    double r3_muonCorr_M1_e[10]={0, 0.0838277, 0.0266269, 0.0168108, 0.0156166, 0.0211333, 0.0286836, 0.0450595, 0.075618, 0.141419};
    TH1D *r3_muonCorr=(TH1D*)r3_num_Q3->Clone();
    for(int i=0; i<10; i++){ 
      r3_muonCorr->SetBinContent(i+1,r3_muonCorr_M1[i]); 
      r3_muonCorr->SetBinError(i+1,sqrt(pow(r3_muonCorr_M1_e[i],2) + pow(r3Sys_M1[i],2)));
    }
    r3_muonCorr->SetMarkerStyle(21);
    r3_muonCorr->SetMarkerColor(2); r3_muonCorr->SetLineColor(2);
    r3_muonCorr->SetFillStyle(1000); r3_muonCorr->SetFillColor(kRed-10);
    r3_muonCorr->Draw("E2 same");
    r3_num_Q3->Draw("same");
    legend3->AddEntry(r3_num_Q3,"No Muon Correction","p");
    legend3->AddEntry(r3_muonCorr,"Correction Applied","p");
    legend3->Draw("same");*/
    
    // muon correction for muon variation data
    /*double r3_muonCorr_M1[10]={0, 1.71115, 1.8128, 1.79761, 1.76112, 1.70699, 1.63908, 1.63417, 1.49861, 1.75565};
    double r3_muonCorr_M1_e[10]={0, 0.155228, 0.0491558, 0.0308125, 0.0283357, 0.0378035, 0.0501422, 0.0756057, 0.119169, 0.202093};
    TH1D *r3_muonCorr=(TH1D*)r3_num_Q3->Clone();
    for(int i=0; i<10; i++){ 
      r3_muonCorr->SetBinContent(i+1,r3_muonCorr_M1[i]); 
      r3_muonCorr->SetBinError(i+1,sqrt(pow(r3_muonCorr_M1_e[i],2) + pow(r3Sys_M1[i],2)));
    }
    r3_muonCorr->SetMarkerStyle(21);
    r3_muonCorr->SetMarkerColor(2); r3_muonCorr->SetLineColor(2);
    r3_muonCorr->SetFillStyle(1000); r3_muonCorr->SetFillColor(kRed-10);
    r3_muonCorr->Draw("E2 same");
    r3_num_Q3->Draw("same");
    legend3->AddEntry(r3_num_Q3,"Muon Corrected.  Default PID","p");
    legend3->AddEntry(r3_muonCorr,"Muon Corrected.  Varied PID","p");
    legend3->Draw("same");*/

    /*for(int i=1; i<=10; i++){
      cout<<r3_num_Q3->GetBinContent(i)<<", ";
    }
    cout<<endl;
    for(int i=1; i<=10; i++){
      cout<<r3_num_Q3->GetBinError(i)<<", ";
      }*/

    /*
    TH1D *Merged_SanityCheck=(TH1D*)r3_num_Q3->Clone();
    for(int i=1; i<=10; i++){
      double value = (Gen_r3_m_M0[i-1]+Gen_r3_p_M0[i-1])/2.;
      double value_e = sqrt(pow(Gen_r3_m_M0_e[i-1],2)+pow(Gen_r3_p_M0_e[i-1],2))/2.;
      //double value = (Gen_r3_m_M0[i-1]+Gen_r3_p_M0[i-1]+Gen_r3_m_M1[i-1]+Gen_r3_p_M1[i-1])/4.;
      //double value_e = sqrt(pow(Gen_r3_m_M0_e[i-1],2)+pow(Gen_r3_p_M0_e[i-1],2)+pow(Gen_r3_m_M1_e[i-1],2)+pow(Gen_r3_p_M1_e[i-1],2))/4.;
      Merged_SanityCheck->SetBinContent(i,value);
      Merged_SanityCheck->SetBinError(i,value_e);
    }
    gPad->SetTopMargin(0.02); gPad->SetLeftMargin(0.1);
    gPad->SetGridx(0); gPad->SetGridy(0);
    Merged_SanityCheck->GetXaxis()->SetTitleOffset(1.2);
    Merged_SanityCheck->GetYaxis()->SetTitleOffset(1.3);
    Merged_SanityCheck->SetMinimum(0.3); Merged_SanityCheck->SetMaximum(2.68);
    Merged_SanityCheck->Draw();
    
    
    //r3_num_Q3->Fit(QuarticFit,"IME","",0,0.1);
    Merged_SanityCheck->Fit(QuarticFit,"IME","",0,0.1);
    gPad->Update();
    //TPaveStats *Quartic_stats=(TPaveStats*)r3_num_Q3->FindObject("stats");
    TPaveStats *Quartic_stats=(TPaveStats*)Merged_SanityCheck->FindObject("stats");
    Quartic_stats->SetX1NDC(0.15);
    Quartic_stats->SetX2NDC(0.52);
    Quartic_stats->SetY1NDC(.2);
    Quartic_stats->SetY2NDC(.3);
    
    //TH1D *r3_clone=(TH1D*)r3_num_Q3->Clone();
    TH1D *r3_clone=(TH1D*)Merged_SanityCheck->Clone();
    r3_clone->Fit(QuadraticFit,"IME","",0,0.1);
    gPad->Update();
    TPaveStats *Quadratic_stats=(TPaveStats*)r3_clone->FindObject("stats");
    Quadratic_stats->SetX1NDC(0.54);
    Quadratic_stats->SetX2NDC(0.91);
    Quadratic_stats->SetY1NDC(.2);
    Quadratic_stats->SetY2NDC(.3);
    
    QuarticFit->Draw("same");
    QuadraticFit->Draw("same");
    Quartic_stats->Draw("same");
    Quadratic_stats->Draw("same");
    
    TF1 *ChaoticLimit = new TF1("ChaoticLimit","2.0",0,1);
    ChaoticLimit->SetLineStyle(3);
    ChaoticLimit->Draw("same");
    
    TLatex *Specif_SanityCheck=new TLatex(0.2,0.35,"HIJING With Simulated QS+FSI Weights");
    Specif_SanityCheck->SetNDC();
    Specif_SanityCheck->SetTextFont(42);
    Specif_SanityCheck->SetTextSize(0.04);
    Specif_SanityCheck->Draw();
    //TLatex *Specif_Kt3=new TLatex(0.2,0.45,"0.16<K_{t,3}<0.3 GeV/c");
    TLatex *Specif_Kt3=new TLatex(0.2,0.45,"0.3<K_{t,3}<1.0 GeV/c");
    Specif_Kt3->SetNDC();
    Specif_Kt3->SetTextFont(42);
    Specif_Kt3->SetTextSize(0.04);
    Specif_Kt3->Draw();
    
    legend3->AddEntry(QuarticFit,"Quartic fit","l");
    legend3->AddEntry(QuadraticFit,"Quadratic fit","l");
    legend3->Draw("same");
    //DrawALICELogo(kFALSE, .72, .4, .87, .55);
    */
    /*
    ////////// 
    // r3 Sys
    cout<<"r3 values:"<<endl;
    for(int ii=0; ii<10; ii++){
      cout<<r3_num_Q3->GetBinContent(ii+1)<<", ";
    }
    cout<<endl;
    cout<<"r3 errors:"<<endl;
    for(int ii=0; ii<10; ii++){
      cout<<r3_num_Q3->GetBinError(ii+1)<<", ";
    }
    cout<<endl;
    
    // M0,pi- base values with TTC=0.035 (old standard)
    //double r3base[10]={0, 2.02584, 1.82659, 1.85384, 1.84294, 1.82431, 1.72205, 1.5982, 1.17379, 0.881535};
    //double r3base_e[10]={0, 0.15115, 0.038311, 0.0231936, 0.0209217, 0.0291552, 0.0406093, 0.0624645, 0.0933104, 0.122683};
    // M0,pi+ base values with TTC=0.035 (old standard)
    //double r3base[10]={0, 1.84738, 1.81805, 1.80057, 1.82563, 1.76585, 1.6749, 1.37454, 1.37558, 0.74562};
    //double r3base_e[10]={0, 0.147285, 0.0381688, 0.0231643, 0.0209571, 0.0292255, 0.0407515, 0.0627228, 0.0937758, 0.123392};
    // M1,pi- base values with TTC=0.035 (old standard)
    //double r3base[10]={0, 1.72641, 1.82547, 1.79747, 1.84141, 1.83557, 1.7994, 1.67511, 1.44238, 1.03004};
    //double r3base_e[10]={0, 0.192385, 0.0471262, 0.0270385, 0.0229849, 0.0302041, 0.0397806, 0.0595701, 0.0902909, 0.123976};
    // M2,pi- base values with TTC=0.035 (old standard)
    //double r3base[10]={0, 1.61726, 1.58973, 1.7834, 1.8914, 1.75445, 1.76511, 1.81367, 1.55617, 1.17797};
    //double r3base_e[10]={0, 0.417924, 0.0973277, 0.0527296, 0.0421998, 0.0522792, 0.065144, 0.0923891, 0.135867, 0.190967};
    // M9,pi- base values with TTC=0.035 (old standard)
    //double r3base[10]={0, 2.51235, 2.68042, 1.85333, 1.66794, 1.39832, 1.717, 1.5337, 1.80686, 1.42286};
    //double r3base_e[10]={0, 3.09403, 0.650924, 0.272827, 0.160562, 0.145494, 0.137767, 0.146212, 0.164592, 0.189126};

    // M0 PosB, pi- base values with TTC=0.035 (old standard)
    //double r3base[10]={0, 1.81361, 1.84779, 1.81693, 1.8252, 1.78366, 1.71822, 1.41243, 1.31656, 0.719842};
    //double r3base_e[10]={0, 0.158656, 0.0405, 0.0244075, 0.0219925, 0.0306193, 0.0426494, 0.0655837, 0.0980259, 0.129011};
    // M1 PosB, pi- base values with TTC=0.035 (old standard)
    //double r3base[10]={0, 1.7208, 1.79016, 1.76858, 1.85008, 1.80416, 1.67274, 1.6501, 1.45098, 1.1052};
    //double r3base_e[10]={0, 0.202391, 0.0496228, 0.0282103, 0.0239023, 0.0313499, 0.041234, 0.061679, 0.0934039, 0.128249};

    // M1 PID1p5 special case (M0 MomRes used),pi- base values with TTC=0.035 (old standard)
    //double r3base[10]={0, 1.77035, 1.82471, 1.77934, 1.80615, 1.77819, 1.71105, 1.60505, 1.44563, 1.09473};
    //double r3base_e[10]={0, 0.193326, 0.0471058, 0.0269928, 0.0229581, 0.0301777, 0.0397553, 0.0595472, 0.0902768, 0.123968};

    // M0,pi- base values with TTC=0.04 (New standard, r3 Interp fix, GRS)
    //double r3base[10]={0, 1.93236, 1.77265, 1.71324, 1.63426, 1.56167, 1.38348, 1.25064, 0.897615, 0.691928};
    //double r3base_e[10]={0, 0.221991, 0.0433103, 0.023234, 0.0187678, 0.0243522, 0.0318668, 0.0459759, 0.066614, 0.0873847};
    // M1,pi- base values with TTC=0.04 (New standard, r3 Interp fix, GRS)
    //double r3base[10]={0, 1.46815, 1.76361, 1.68282, 1.6334, 1.56039, 1.44628, 1.27343, 1.01443, 0.691937};
    //double r3base_e[10]={0, 0.277164, 0.053134, 0.0271048, 0.0207043, 0.0253691, 0.0315642, 0.0443283, 0.0641466, 0.0858159};
    // M2,pi- base values with TTC=0.04 (New standard, r3 Interp fix, GRS)
    //double r3base[10]={0, 1.43331, 1.72759, 1.71155, 1.68957, 1.53846, 1.45209, 1.40698, 1.16202, 1.04672};
    //double r3base_e[10]={0, 0.594776, 0.109832, 0.0531013, 0.0383514, 0.0444564, 0.0525963, 0.0705969, 0.0998679, 0.135742};
    // M9,pi- base values with TTC=0.04 (New standard, r3 Interp fix, GRS)
    //double r3base[10]={0, 6.98181, 2.42845, 1.71962, 1.57641, 1.23102, 1.49543, 1.48063, 1.50938, 1.27714};
    //double r3base_e[10]={0, 7.09254, 0.731841, 0.277336, 0.152199, 0.130775, 0.118981, 0.121864, 0.132875, 0.148409};
    //
    // M0,pi- base values with 3sigma
    //double r3base[10]={0, 1.71114, 1.72142, 1.74749, 1.69154, 1.60411, 1.47805, 1.27561, 1.02355, 0.769853};
    //double r3base_e[10]={0, 0.251731, 0.0398122, 0.0196354, 0.0153904, 0.0197609, 0.0258499, 0.0377931, 0.0560878, 0.0774703};
    // M1,pi- base values with 3sigma
    //double r3base[10]={0, 2.02471, 1.78737, 1.68681, 1.69639, 1.62634, 1.57337, 1.42229, 1.28352, 0.958121};
    //double r3base_e[10]={0, 0.314475, 0.048507, 0.022855, 0.0169939, 0.0207322, 0.0260424, 0.0375715, 0.0569927, 0.0830874};
    // M2,pi- base values with 3sigma
    //double r3base[10]={0, 1.50807, 1.75138, 1.71876, 1.76974, 1.6322, 1.5883, 1.5079, 1.4551, 1.49369};
    //double r3base_e[10]={0, 0.658863, 0.0979719, 0.0436736, 0.0308217, 0.035654, 0.0425953, 0.0585648, 0.0865817, 0.128056};
    // M9,pi- base values with 3sigma
    //double r3base[10]={0, 1.17509, 1.5117, 1.98348, 1.63032, 1.57195, 1.48828, 1.48218, 1.53689, 1.39533};
    //double r3base_e[10]={0, 4.89181, 0.60212, 0.230712, 0.11983, 0.102686, 0.0937373, 0.0971281, 0.107949, 0.123923};
    //
    // M0,pi-,EDbin=0, Femto Plus+Minus
    //double r3base[10]={0, 1.77055, 1.75926, 1.7711, 1.67159, 1.57851, 1.41187, 1.11004, 0.870851, 0.442019};
    //double r3base_e[10]={0, 0.0688956, 0.0232182, 0.0156323, 0.0155991, 0.0225749, 0.0324584, 0.0516114, 0.0797203, 0.112197};
    // M1,pi-,EDbin=0, Femto Plus+Minus
    //double r3base[10]={0, 1.83169, 1.79042, 1.7275, 1.67116, 1.6032, 1.5423, 1.30426, 1.13238, 0.656902};
    //double r3base_e[10]={0, 0.0891247, 0.0286553, 0.0182713, 0.0171956, 0.0235918, 0.0324304, 0.0510863, 0.0823637, 0.125496};
    // M2,pi-,EDbin=0, Femto Plus+Minus
    //double r3base[10]={0, 1.60769, 1.81565, 1.75909, 1.74026, 1.62865, 1.58101, 1.50583, 1.3503, 1.27367};
    //double r3base_e[10]={0, 0.188422, 0.0581482, 0.0348698, 0.0309363, 0.0401393, 0.0519907, 0.0771148, 0.120583, 0.187679};
    // M3,pi-,EDbin=0, Femto Plus+Minus
    double r3base[10]={0, 1.52778, 1.81442, 1.78982, 1.69978, 1.70934, 1.59516, 1.56605, 1.45676, 1.17468};
    double r3base_e[10]={0, 0.288329, 0.0882351, 0.0510576, 0.0432015, 0.0536294, 0.0667149, 0.0950229, 0.143021, 0.217322};
    // M8,pi-,EDbin=0, Femto Plus+Minus
    //double r3base[10]={0, 2.27921, 2.40794, 1.72409, 1.66853, 1.70966, 1.52681, 1.64425, 1.18747, 1.35988};
    //double r3base_e[10]={0, 1.08315, 0.289559, 0.136951, 0.0919134, 0.0900869, 0.0901754, 0.101887, 0.121289, 0.148063};
    // M9,pi-,EDbin=0, Femto Plus+Minus
    //double r3base[10]={0, 2.1271, 2.12837, 1.95071, 1.4832, 1.50308, 1.41068, 1.29826, 1.10739, 0.755899};
    //double r3base_e[10]={0, 1.52325, 0.385705, 0.181328, 0.116165, 0.10896, 0.105705, 0.115163, 0.131994, 0.156299};

    //
    // M9,pi+ (No MRC)
    double r3base_noMRC[10]={0, 4.16425, 0.429681, 1.56506, 1.64596, 1.73785, 1.57181, 1.51971, 1.59096, 1.54403};
    double r3base_noMRC_e[10]={0, 2.00855, 0.437928, 0.187872, 0.1087, 0.097367, 0.0893871, 0.0933358, 0.103603, 0.11839};
    // M9,pi+ (MRC normal)
    double r3base_MRC[10]={0, 3.84986, 0.558354, 1.60733, 1.67855, 1.75362, 1.59637, 1.54972, 1.62045, 1.57307};
    double r3base_MRC_e[10]={0, 1.93428, 0.402333, 0.171664, 0.0996222, 0.0903509, 0.0840615, 0.0889731, 0.100023, 0.115704};
    // M9,pi+ (MRC 50% larger)
    double r3base_MRC_over[10]={0, 3.70355, 0.614999, 1.62641, 1.69371, 1.76209, 1.60931, 1.56565, 1.63663, 1.58935};
    double r3base_MRC_over_e[10]={0, 1.90276, 0.387084, 0.164703, 0.0956857, 0.0872585, 0.0816823, 0.0870033, 0.0984005, 0.114503};

    

    TH1D *r3_noMRC=(TH1D*)r3_num_Q3->Clone();
    TH1D *r3_MRC=(TH1D*)r3_num_Q3->Clone();
    TH1D *r3_MRC_over=(TH1D*)r3_num_Q3->Clone();

    TH1D *r3Sys=(TH1D*)r3_num_Q3->Clone();
    for(int ii=1; ii<10; ii++){
      double Stat_RelativeError = sqrt(fabs(pow(r3_num_Q3->GetBinError(ii+1),2) - pow(r3base_e[ii],2)));// + for independent stat's
      r3Sys->SetBinContent(ii+1, (r3base[ii]-r3_num_Q3->GetBinContent(ii+1))/r3base[ii]);
      r3Sys->SetBinError(ii+1, Stat_RelativeError);
      //
      r3_noMRC->SetBinContent(ii+1, r3base_noMRC[ii]);
      r3_MRC->SetBinContent(ii+1, r3base_MRC[ii]);
      r3_MRC_over->SetBinContent(ii+1, r3base_MRC_over[ii]);
      r3_noMRC->SetBinError(ii+1, r3base_noMRC_e[ii]);
      r3_MRC->SetBinError(ii+1, r3base_MRC_e[ii]);
      r3_MRC_over->SetBinError(ii+1, r3base_MRC_over_e[ii]);
    }
    r3_noMRC->SetMarkerStyle(20); r3_noMRC->SetMarkerColor(1); r3_noMRC->SetLineColor(1);
    r3_MRC->SetMarkerStyle(23); r3_MRC->SetMarkerColor(4); r3_MRC->SetLineColor(4);
    r3_MRC_over->SetMarkerStyle(24); r3_MRC_over->SetMarkerColor(2); r3_MRC_over->SetLineColor(2);
    
    //r3_MRC->Draw();
    //r3_noMRC->Draw("same");
    //r3_MRC_over->Draw("same");
    TLegend *legendMRC = new TLegend(.45,.15,.95,.3,NULL,"brNDC");
    legendMRC->SetBorderSize(1);
    legendMRC->SetTextSize(.03);// small .03; large .06
    legendMRC->SetFillColor(0);
    //legendMRC->AddEntry(r3_MRC,"Momentum Resolution Corrected","p");
    //legendMRC->AddEntry(r3_noMRC,"No Correction","p");
    //legendMRC->AddEntry(r3_MRC_over,"50% Increased Correction","p");
    //legendMRC->Draw("same");

    r3Sys->GetYaxis()->SetTitle("(r_{3}^{t1}-r_{3}^{t2})/r_{3}^{t1}");
    r3Sys->GetYaxis()->SetTitleOffset(1.8);
    r3Sys->SetMinimum(-0.1);
    r3Sys->SetMaximum(0.1);
    r3Sys->GetXaxis()->SetRangeUser(0,0.099);
    r3Sys->Draw();
    TF1 *lineFit=new TF1("lineFit","pol0",0,0.02);
    gStyle->SetOptFit(111);
    r3Sys->Fit(lineFit,"MEQ","",0,0.1);
    */

  
    /*
    ////////////////////////////////////////////////////////////////////////
    // The STAR method
    ReadCoulCorrections(0, 5, 2, 10);// Change to Gaussian R=5 fm calculation
    
    TLegend *legendSTAR = new TLegend(.45,.15,.95,.3,NULL,"brNDC");
    legendSTAR->SetBorderSize(1);
    legendSTAR->SetTextSize(.03);// small .03; large .06
    legendSTAR->SetFillColor(0);
    
    TH1D *r3_STAR = new TH1D("r3_STAR","",Q3BINS,0,0.2);
    TH1D *r3_den_STAR = new TH1D("r3_den_STAR","",Q3BINS,0,0.2);
    double r3_STAR_num_e[20]={0};
    double r3_STAR_den_e[20]={0};

    for(int i=2; i<=20; i++){// bin number
      double Q12 = (i-0.5)*(0.005);
      int Q12bin = int(Q12/0.005)+1;
      //
      for(int j=2; j<=20; j++){// bin number
	double Q13 = (j-0.5)*(0.005);
	int Q13bin = int(Q13/0.005)+1;
	//
	for(int k=2; k<=20; k++){// bin number
	  double Q23 = (k-0.5)*(0.005);
	  int Q23bin = int(Q23/0.005)+1;
	  //
	  double Q3 = sqrt(pow(Q12,2) + pow(Q13,2) + pow(Q23,2));
	  int Q3bin = Cterm1_noMRC->GetXaxis()->FindBin(Q3);
	  if(Q3bin>19) continue;

	  if(Q12 < sqrt(pow(Q13,2)+pow(Q23,2) - 2*Q13*Q23)) continue;
	  if(Q12 > sqrt(pow(Q13,2)+pow(Q23,2) + 2*Q13*Q23)) continue;
	  
	  double K2_12 = CoulCorr2(CP12, Q12);
	  double K2_13 = CoulCorr2(CP13, Q13);
	  double K2_23 = CoulCorr2(CP23, Q23);
	  
	  double r3num = Cterm1_noMRC->GetBinContent(Q3bin)/(K2_12*K2_13*K2_23) - 1;
	  r3num -= C2ssRaw->GetBinContent(Q12bin)/(K2_12) - 1;
	  r3num -= C2ssRaw->GetBinContent(Q13bin)/(K2_13) - 1;
	  r3num -= C2ssRaw->GetBinContent(Q23bin)/(K2_23) - 1;
	  double r3den = (C2ssRaw->GetBinContent(Q12bin)/(K2_12) - 1);
	  r3den *= (C2ssRaw->GetBinContent(Q13bin)/(K2_13) - 1);
	  r3den *= (C2ssRaw->GetBinContent(Q23bin)/(K2_23) - 1);
	  if(r3den<0) continue;
	  r3den = sqrt(r3den);
	  //
	  r3_STAR_num_e[Q3bin-1] += pow(Cterm1_noMRC->GetBinError(Q3bin)/(K2_12*K2_13*K2_23),2);
	  r3_STAR_num_e[Q3bin-1] += pow(C2ssRaw->GetBinError(Q12bin)/(K2_12),2);
	  r3_STAR_num_e[Q3bin-1] += pow(C2ssRaw->GetBinError(Q13bin)/(K2_13),2);
	  r3_STAR_num_e[Q3bin-1] += pow(C2ssRaw->GetBinError(Q23bin)/(K2_23),2);
	  r3_STAR_den_e[Q3bin-1] += pow(0.5*C2ssRaw->GetBinError(Q12bin)/(K2_12) * (C2ssRaw->GetBinContent(Q13bin)/(K2_13) - 1)*(C2ssRaw->GetBinContent(Q23bin)/(K2_23) - 1) /r3den,2);
	  r3_STAR_den_e[Q3bin-1] += pow(0.5*C2ssRaw->GetBinError(Q13bin)/(K2_13) * (C2ssRaw->GetBinContent(Q12bin)/(K2_12) - 1)*(C2ssRaw->GetBinContent(Q23bin)/(K2_23) - 1) /r3den,2);
	  r3_STAR_den_e[Q3bin-1] += pow(0.5*C2ssRaw->GetBinError(Q23bin)/(K2_23) * (C2ssRaw->GetBinContent(Q12bin)/(K2_12) - 1)*(C2ssRaw->GetBinContent(Q13bin)/(K2_13) - 1) /r3den,2);
	  //
	  r3_STAR->Fill(Q3, r3num);
	  r3_den_STAR->Fill(Q3, r3den);
	}
      }
    }
    
    
    for(int i=0; i<20; i++){
      r3_STAR->SetBinError(i+1, sqrt( r3_STAR_num_e[i]));
      r3_den_STAR->SetBinError(i+1, sqrt( r3_STAR_den_e[i]));
    }
    r3_STAR->Divide(r3_den_STAR);
    r3_STAR->SetMaximum(5);
    r3_STAR->SetMinimum(-10);
    r3_STAR->GetXaxis()->SetRangeUser(0,0.099);
    r3_STAR->GetXaxis()->SetTitle("Q_{3} (GeV/c)");
    r3_STAR->GetYaxis()->SetTitle("r_{3}*");
    r3_STAR->SetLineColor(2); r3_STAR->SetMarkerColor(2);
    r3_STAR->SetMarkerStyle(kFullStar);
    r3_STAR->SetMarkerSize(1.2);
    r3_STAR->Draw();
    r3_STAR->Fit(QuarticFit,"IME","",0.01,0.1);
    gPad->Update();
    TPaveStats *Quartic_stats_STAR=(TPaveStats*)r3_STAR->FindObject("stats");
    Quartic_stats_STAR->SetX1NDC(0.2);
    Quartic_stats_STAR->SetX2NDC(0.5);
    Quartic_stats_STAR->SetY1NDC(.8);
    Quartic_stats_STAR->SetY2NDC(.9);
    
    TH1D *r3_clone_STAR=(TH1D*)r3_STAR->Clone();
    r3_STAR->Fit(QuadraticFit,"IME","",0,0.1);
    gPad->Update();
    TPaveStats *Quadratic_stats_STAR=(TPaveStats*)r3_clone_STAR->FindObject("stats");
    Quadratic_stats_STAR->SetX1NDC(0.55);
    Quadratic_stats_STAR->SetX2NDC(0.85);
    Quadratic_stats_STAR->SetY1NDC(.8);
    Quadratic_stats_STAR->SetY2NDC(.9);
    
    QuarticFit->Draw("same");
    QuadraticFit->Draw("same");
    Quartic_stats_STAR->Draw("same");
    Quadratic_stats_STAR->Draw("same");
    //Cterm1_noMRC->Draw();
    
    /////////////////////////////////////////////////////////////////////////
    if(SameCharge) r3_num_Q3->Draw("same");
    legendSTAR->AddEntry(r3_STAR,"STAR method","p");
    legendSTAR->AddEntry(r3_num_Q3,"ALICE method","p");
    legendSTAR->Draw("same");
    */
  }// SameCharge

  //c3_hist_STAR->Draw();
  //Cterm1_noMRC->Draw();

  //r3_num_Q3->SetMarkerStyle(20);
  //r3_num_Q3->GetXaxis()->SetRangeUser(0,0.1);
  //r3_num_Q3->SetMinimum(-7);
  //r3_num_Q3->SetMaximum(2.2);
 

  /*
    double C3_star_noCoul[20]={0, 1.59478, 1.49099, 1.4239, 1.31166, 1.21337, 1.14344, 1.09671, 1.06525, 1.04391, 1.0306, 1.02045, 1.01359, 1.0035, 0.992264, 0.989587, 0.989968, 0.996795, 0, 0};
  double C3_star_noCoul_e[20]={0, 0.260834, 0.0191839, 0.00524928, 0.0021307, 0.00110448, 0.000661211, 0.000438894, 0.000311247, 0.00023452, 0.000182819, 0.000147496, 0.000126547, 0.000130428, 0.000159494, 0.000240658, 0.000489113, 0.00387954, 0, 0};
  TH1D *C3_star=(TH1D*)Cterm1->Clone();
  C3_star->Add(C3_star,-1);
  for(int i=0; i<20; i++) {C3_star->SetBinContent(i+1, C3_star_noCoul[i]); C3_star->SetBinError(i+1, C3_star_noCoul_e[i]);}
  C3_star->SetMarkerColor(4);
  C3_star->SetLineColor(4);
  //C3_star->Draw("same");
  //legend1->AddEntry(C3_star,"---; 200 GeV","p");

  //legend1->Draw("same");

  //TGraph *gr_MomRes_term1=new TGraph(10,
  */
  


  
  
  //////////////////////////////////////////////////////////////////////////////////////

  


 
  
  if(SaveToFile){
    TString *savefilename = new TString("OutFiles/OutFile");
    if(SameCharge) savefilename->Append("SC");
    else savefilename->Append("MC");
    //
    if(CHARGE==1) savefilename->Append("Pos");
    else savefilename->Append("Neg");
    //
    if(IncludeG) savefilename->Append("YesG");
    else savefilename->Append("NoG");
    //
    if(IncludeEWfromTherm) savefilename->Append("EWfromTherm");
    if(IncludeEW) savefilename->Append("EW");
    if(GRS) savefilename->Append("GRS");
    else savefilename->Append("Omega0");

    if(EDbin==0) savefilename->Append("Kt3_1");
    else savefilename->Append("Kt3_2");
    
    savefilename->Append("_Kt");
    *savefilename += Ktbin;
    savefilename->Append("_M");
    *savefilename += Mbin;
    savefilename->Append(".root");
    TFile *savefile = new TFile(savefilename->Data(),"RECREATE");
    can1->Write("can1");
    C2_ss->Write("C2_ss");
    C2_os->Write("C2_os");
    C2_ss_Momsys->Write("C2_ss_Momsys");
    C2_os_Momsys->Write("C2_os_Momsys");
    C2_ss_Ksys->Write("C2_ss_Ksys");
    C2_os_Ksys->Write("C2_os_Ksys");
    fitC2ss->Write("fitC2ss");
    fitC2os->Write("fitC2os");
    fitC2ss_h->Write("fitC2ss_h");
    fitC2os_h->Write("fitC2os_h");
    if(IncludeEWfromTherm) {
      fitC2ssTherm->Write("fitC2ssTherm");
      fitC2osTherm->Write("fitC2osTherm");
      fitC2ssThermGaus->Write("fitC2ssThermGaus");
      fitC2osThermGaus->Write("fitC2osThermGaus");
    }
    C2Therm_ss->Write("C2Therm_ss");
    C2Therm_os->Write("C2Therm_os");
    K2_ss->Write("K2_ss");
    K2_os->Write("K2_os");
    Cterm1->Write("C3");
    Coul_Omega0->Write("Coul_Omega0");
    Coul_GRiverside->Write("Coul_GRiverside");
    GenSignalExpected_num->Write("C3_EWexpectation");
    num_QS->Write("num_QS");
    c3_hist->Write("c3");
    Combinatorics_1d->Write("Combinatorics_1d");
    ChiSquaredNDF->Write("ChiSquaredNDF");
    r3_num_Q3->Write("r3_Q3");
    
    //
    savefile->Close();
  }
  

}

void ReadCoulCorrections(int RVal, int bVal, int kt){
  ///////////////////////
  
  TString *fname2;
  if(FileSetting==6) fname2 = new TString("KFile_GaussR8to5.root");
  else if(FileSetting==7) fname2 = new TString("KFile_GaussR11to6.root");
  else{
    if(SourceType==0) fname2 = new TString("KFile.root");
    if(SourceType==1) fname2 = new TString("KFile_50fmCut.root");
    if(SourceType==2) fname2 = new TString("KFile_GaussR11to6.root");
  }
  
  TFile *File=new TFile(fname2->Data(),"READ");
  if(bVal!=2 && bVal!=3 && bVal!=5 && bVal!=7 && bVal!=8 && bVal!=9) cout<<"Therminator bVal not acceptable in 2-particle Coulomb read"<<endl;
  
  if(kt==10){// kt integrated
    TH2D *tempT_ss = (TH2D*)File->Get("K2ssT");
    TH2D *tempT_os = (TH2D*)File->Get("K2osT");
    CoulCorr2SS = (TH1D*)tempT_ss->ProjectionY("CoulCorr2SS",bBin, bBin);
    CoulCorr2OS = (TH1D*)tempT_os->ProjectionY("CoulCorr2OS",bBin, bBin);
  }else{
    if(kt < 1 || kt > 6) cout<<"kt bin out of range in 2-particle Coulomb read"<<endl;
    TH3D *tempT3_ss = (TH3D*)File->Get("K2ssT_kt");
    TH3D *tempT3_os = (TH3D*)File->Get("K2osT_kt");
    CoulCorr2SS = (TH1D*)tempT3_ss->ProjectionZ("CoulCorr2SS",bBin, bBin, kt,kt);
    CoulCorr2OS = (TH1D*)tempT3_os->ProjectionZ("CoulCorr2OS",bBin, bBin, kt,kt);
  }
  
  if(IncludeStrongSC){// include strong FSI for pi+pi+ as per Lednicky code with R=7 (ratio of Coul+Strong to Coul) 
    for(int ii=1; ii<=CoulCorr2SS->GetNbinsX(); ii++){
      double newValue = CoulCorr2SS->GetBinContent(ii) * StrongSC->Eval(CoulCorr2SS->GetBinCenter(ii)*1000/2.);// k* not qinv
      CoulCorr2SS->SetBinContent(ii, newValue);
    }
  }
  // K factor shift for testing
  for(int ii=1; ii<=CoulCorr2SS->GetNbinsX(); ii++){
    double newValueSS = (CoulCorr2SS->GetBinContent(ii)-1)*KShift+1;// k* not qinv
    CoulCorr2SS->SetBinContent(ii, newValueSS);
    double newValueOS = (CoulCorr2OS->GetBinContent(ii)-1)*KShift+1;// k* not qinv
    CoulCorr2OS->SetBinContent(ii, newValueOS);
  }

  CoulCorr2SS->SetDirectory(0);
  CoulCorr2OS->SetDirectory(0);
  File->Close(); 
}
double CoulCorr2(int chargeproduct, double Q2){// returns K2
  int indexL=0;
  int indexH=0;
  double slope=0;
 
  
  indexL = int(fabs(Q2 - CoulCorr2SS->GetXaxis()->GetBinWidth(1)/2.)/(CoulCorr2SS->GetXaxis()->GetBinWidth(1)));
  indexH = indexL+1;
 
  if(indexH >= CoulCorr2SS->GetNbinsX()) return 1.0;
  if(chargeproduct==1){
    slope = (CoulCorr2SS->GetBinContent(indexL+1) - CoulCorr2SS->GetBinContent(indexH+1));
    slope /= (CoulCorr2SS->GetXaxis()->GetBinCenter(indexL+1) - CoulCorr2SS->GetXaxis()->GetBinCenter(indexH+1));
    return slope*(Q2 - CoulCorr2SS->GetXaxis()->GetBinCenter(indexL+1)) + CoulCorr2SS->GetBinContent(indexL+1);
  }else {
    slope = (CoulCorr2OS->GetBinContent(indexL+1) - CoulCorr2OS->GetBinContent(indexH+1));
    slope /= (CoulCorr2OS->GetXaxis()->GetBinCenter(indexL+1) - CoulCorr2OS->GetXaxis()->GetBinCenter(indexH+1));
    return slope*(Q2 - CoulCorr2OS->GetXaxis()->GetBinCenter(indexL+1)) + CoulCorr2OS->GetBinContent(indexL+1);
  }
  
}

//----------------------------------------------------------------------


//________________________________________________________________________
void fcnC2_Global(int& npar, double* deriv, double& f, double par[], int flag){
  
  double qinvSS=0;
  
  double Rch=par[3]/FmToGeV;
  double Rcoh=par[4]/FmToGeV;
  double Dp=0;
  double CSS=0, COS=0;
  double SumChi2=0;
  double EW=0;
  double MT = sqrt(pow(0.25 + 0.1*Ktbin_GofP,2) + pow(masspiC,2));
  NFitPoints_C2global=0;
  if(LinkN) par[9]=par[0];// Link N factors

  for(int i=1; i<BINRANGE_C2global; i++){
    
    qinvSS = AvgQinvSS[i];
    //if(qinvSS>0.08) continue;

    if(!GofP) Dp=fabs(par[2])/(1-fabs(par[2]));// p independent
    else Dp = fabs(par[2])/(1-fabs(par[2]));
    
    //double Dp1 = fabs(par[2])/(1-fabs(par[2])) * exp(-2*(pow(Rcoh,2) - pow(RT,2))*pow(AvgP[Ktbin_GofP-1][i]-qinvSS/2.,2));
    //double Dp2 = fabs(par[2])/(1-fabs(par[2])) * exp(-2*(pow(Rcoh,2) - pow(RT,2))*pow(AvgP[Ktbin_GofP-1][i]+qinvSS/2.,2));
    //double Dp1 = fabs(par[2])/(1-fabs(par[2])) * exp(-2*pow(Rcoh,2)*pow(AvgP[Ktbin_GofP-1][i]-qinvSS/2.,2) - 2*MT/Temp);
    //double Dp2 = fabs(par[2])/(1-fabs(par[2])) * exp(-2*pow(Rcoh,2)*pow(AvgP[Ktbin_GofP-1][i]+qinvSS/2.,2) - 2*MT/Temp);
    //double Dp1 = fabs(par[2])/(1-fabs(par[2])) * exp(-2*pow(Rcoh,2)*pow(par[10]-qinvSS/2.,2) - 2*MT/Temp);
    //double Dp2 = fabs(par[2])/(1-fabs(par[2])) * exp(-2*pow(Rcoh,2)*pow(par[10]+qinvSS/2.,2) - 2*MT/Temp);
    //double Dp1 = fabs(par[2])/(1-fabs(par[2])) * exp(-2*pow(Rcoh,2)*pow(qinvSS/2.,2) - 2*MT/Temp);
    //double Dp2 = fabs(par[2])/(1-fabs(par[2])) * exp(-2*pow(Rcoh,2)*pow(qinvSS/2.,2) - 2*MT/Temp);
    double Dp1 = fabs(par[2])/(1-fabs(par[2])) * exp(2*(pow(Rcoh,2)-pow(RT,2))*pow(qinvSS/2.,2));
    double Dp2 = fabs(par[2])/(1-fabs(par[2])) * exp(2*(pow(Rcoh,2)-pow(RT,2))*pow(qinvSS/2.,2));
    
    if(!GofP) {Dp1=Dp; Dp2=Dp;}
    
    //
    double arg=qinvSS*Rch;
    EW = 1 + par[5]/(6.*pow(2,1.5))*(8.*pow(arg,3) - 12.*pow(arg,1));
    EW += par[6]/(24.*pow(2,2))*(16.*pow(arg,4) -48.*pow(arg,2) + 12);
    EW += par[7]/(120.*pow(2,2.5))*(32.*pow(arg,5) - 160.*pow(arg,3) + 120*arg);
    EW += par[8]/(720.*pow(2,3))*(64.*pow(arg,6) - 480.*pow(arg,4) + 720.*pow(arg,2) - 120);
    //
    double Gaus_coh = exp(-pow(Rcoh*qinvSS,2)/2.);
    double Gaus_ch = exp(-pow(Rch*qinvSS,2)/2.) * EW;
    CSS = 1 + pow(1 + Dp*Gaus_coh/Gaus_ch,2)/((1+Dp1)*(1+Dp2)) * exp(-pow(Rch*qinvSS,2))*pow(EW,2);
    if(ChargeConstraint) CSS -= 4/5.* Dp1/(1+Dp1) * Dp2/(1+Dp2);
    else CSS -= pow(Gaus_coh*Dp,2)/((1+Dp1)*(1+Dp2));
    //else CSS -= Dp1/(1+Dp1) * Dp2/(1+Dp2);

    CSS *= par[1]*K2SS[i];
    if(MuonCorrection) CSS /= C2muonCorrection_SC->GetBinContent(i+1);
    CSS += 1-par[1];
    CSS *= par[0];
    //
    COS = 1;
    if(ChargeConstraint && GofP) COS += 1/5.* Dp1/(1+Dp1) * Dp2/(1+Dp2);
    COS *= par[1]*K2OS[i];
    if(MuonCorrection) COS /= C2muonCorrection_MC->GetBinContent(i+1);
    COS += 1-par[1];
    COS *= par[9];// different Norm
    //
    if(C2ssFitting[i] > 0){
      if(IncludeSS) {
	double errorSS = sqrt(pow( (CSS/par[0] - (1-par[1]))/K2SS[i] * K2SS_e[i] * par[0],2) + pow(C2ssFitting_e[i],2));
	//double errorSS = C2ssFitting_e[i];
	SumChi2 += pow((CSS - C2ssFitting[i])/errorSS,2);
	NFitPoints_C2global++;
      }
    }
    if(IncludeOS) {
      double errorOS = sqrt(pow( (COS/par[9] - (1-par[1]))/K2OS[i] * K2OS_e[i] * par[9],2) + pow(C2osFitting_e[i],2));
      //double errorOS = C2osFitting_e[i];
      SumChi2 += pow((COS - C2osFitting[i])/errorOS,2);
      NFitPoints_C2global++;
    }
  }
    
  
  f = SumChi2;
  Chi2_C2global = f;
  
}
//________________________________________________________________________
double C2ssFitFunction(Double_t *x, Double_t *par)
{
  double Rch=par[3]/FmToGeV;
  double Rcoh=par[4]/FmToGeV;
  double Dp=0;
  int qbin = int(fabs(x[0]/0.005));
  if(qbin >= BINRANGE_C2global) return 1.0;
  double qinvSS = AvgQinvSS[qbin];
  double MT = sqrt(pow(0.25 + 0.1*Ktbin_GofP,2) + pow(masspiC,2));

  if(!GofP) Dp = fabs(par[2])/(1-fabs(par[2]));// p independent
  else Dp = fabs(par[2])/(1-fabs(par[2]));
 
  //double Dp1 = fabs(par[2])/(1-fabs(par[2])) * exp(-2*(pow(Rcoh,2)-pow(RT,2))*pow(AvgP[Ktbin_GofP-1][qbin]-qinvSS/2.,2));
  //double Dp2 = fabs(par[2])/(1-fabs(par[2])) * exp(-2*(pow(Rcoh,2)-pow(RT,2))*pow(AvgP[Ktbin_GofP-1][qbin]+qinvSS/2.,2));
  //double Dp1 = fabs(par[2])/(1-fabs(par[2])) * exp(-2*pow(Rcoh,2)*pow(AvgP[Ktbin_GofP-1][qbin]-qinvSS/2.,2) - 2*MT/Temp);
  //double Dp2 = fabs(par[2])/(1-fabs(par[2])) * exp(-2*pow(Rcoh,2)*pow(AvgP[Ktbin_GofP-1][qbin]+qinvSS/2.,2) - 2*MT/Temp);
  //double Dp1 = fabs(par[2])/(1-fabs(par[2])) * exp(-2*pow(Rcoh,2)*pow(qinvSS/2.,2) - 2*MT/Temp);
  //double Dp2 = fabs(par[2])/(1-fabs(par[2])) * exp(-2*pow(Rcoh,2)*pow(qinvSS/2.,2) - 2*MT/Temp);
  double Dp1 = fabs(par[2])/(1-fabs(par[2])) * exp(2*(pow(Rcoh,2)-pow(RT,2))*pow(qinvSS/2.,2));
  double Dp2 = fabs(par[2])/(1-fabs(par[2])) * exp(2*(pow(Rcoh,2)-pow(RT,2))*pow(qinvSS/2.,2));
  
  if(!GofP) {Dp1=Dp; Dp2=Dp;}
  double arg=qinvSS*Rch;
  double EW = 1 + par[5]/(6.*pow(2,1.5))*(8.*pow(arg,3) - 12.*pow(arg,1));
  EW += par[6]/(24.*pow(2,2))*(16.*pow(arg,4) -48.*pow(arg,2) + 12);
  EW += par[7]/(120.*pow(2,2.5))*(32.*pow(arg,5) - 160.*pow(arg,3) + 120*arg);
  EW += par[8]/(720.*pow(2,3))*(64.*pow(arg,6) - 480.*pow(arg,4) + 720.*pow(arg,2) - 120);
  double Gaus_coh = exp(-pow(Rcoh*qinvSS,2)/2.);
  double Gaus_ch = exp(-pow(Rch*qinvSS,2)/2.) * EW;
  double CSS = 1 + pow(1 + Dp*Gaus_coh/Gaus_ch,2)/((1+Dp1)*(1+Dp2)) * exp(-pow(Rch*qinvSS,2))*pow(EW,2);
  if(ChargeConstraint) CSS -= 4/5.* Dp1/(1+Dp1) * Dp2/(1+Dp2);
  else CSS -= pow(Gaus_coh*Dp,2)/((1+Dp1)*(1+Dp2));
  //else CSS -= Dp1/(1+Dp1) * Dp2/(1+Dp2);
  CSS *= par[1]*K2SS[qbin];
  if(MuonCorrection) CSS /= C2muonCorrection_SC->GetBinContent(qbin+1);
  CSS += 1-par[1];
  CSS *= par[0];

  //cout<<qinvSS<<"  "<<Dp1/(1+Dp1) * Dp2/(1+Dp2)<<endl;

  return CSS;
}
//________________________________________________________________________
double C2osFitFunction(Double_t *x, Double_t *par)
{
  if(LinkN) par[9]=par[0];// Link N factors
  double Rcoh=par[4]/FmToGeV;
  double Dp=0;
  int qbin = int(fabs(x[0]/0.005));
  if(qbin >= BINRANGE_C2global) return 1.0;
  double qinvOS = AvgQinvOS[qbin];
  
  if(!GofP) Dp = fabs(par[2])/(1-fabs(par[2]));// p independent
  else Dp = fabs(par[2])/(1-fabs(par[2])) * exp(-2*(pow(Rcoh,2)-pow(RT,2))*pow(AvgP[Ktbin_GofP-1][qbin],2));
  //Dp = fabs(par[2])/(1-fabs(par[2]));
  double Dp1 = fabs(par[2])/(1-fabs(par[2])) * exp(-2*(pow(Rcoh,2)-pow(RT,2))*pow(AvgP[Ktbin_GofP-1][qbin]-qinvOS/2.,2));
  double Dp2 = fabs(par[2])/(1-fabs(par[2])) * exp(-2*(pow(Rcoh,2)-pow(RT,2))*pow(AvgP[Ktbin_GofP-1][qbin]+qinvOS/2.,2));
  //double Dp1 = fabs(par[2])/(1-fabs(par[2])) * exp(-2*(pow(Rcoh,2)-pow(RT,2))*pow(qinvOS/2.,2));
  //double Dp2 = fabs(par[2])/(1-fabs(par[2])) * exp(-2*(pow(Rcoh,2)-pow(RT,2))*pow(qinvOS/2.,2));

  if(!GofP) {Dp1=Dp; Dp2=Dp;}
  double COS = 1;
  if(ChargeConstraint && GofP) COS += 1/5.* Dp1/(1+Dp1) * Dp2/(1+Dp2);
  COS *= par[1]*K2OS[qbin];
  if(MuonCorrection) COS /= C2muonCorrection_MC->GetBinContent(qbin+1);
  COS += 1-par[1];
  COS *= par[9];
  return COS;
}
//________________________________________________________________________
void fcnC3_Global(int& npar, double* deriv, double& f, double par[], int flag){

  double qinvSS_12=0;
  double qinvSS_13=0;
  double qinvSS_23=0;

  //double D=fabs(par[2])/(1-fabs(par[2]));
  double Rch=par[3]/FmToGeV;
  //double Rcoh=par[4]/FmToGeV;
  double CSS=0;
  double SumChi2=0;
  double EW_12=0, EW_13=0, EW_23=0;
  NFitPoints_C3global=0;

  for(int i=0; i<BINRANGE_C2global; i++){
    qinvSS_12 = AvgQinvSS[i];
    double arg_12=qinvSS_12*Rch;
    EW_12 = 1 + par[5]/(6.*pow(2,1.5))*(8.*pow(arg_12,3) - 12.*pow(arg_12,1));
    EW_12 += par[6]/(24.*pow(2,2))*(16.*pow(arg_12,4) -48.*pow(arg_12,2) + 12);
    EW_12 += par[7]/(120.*pow(2,2.5))*(32.*pow(arg_12,5) - 160.*pow(arg_12,3) + 120*arg_12);
    EW_12 += par[8]/(720.*pow(2,3))*(64.*pow(arg_12,6) - 480.*pow(arg_12,4) + 720.*pow(arg_12,2) - 120);

    for(int j=0; j<BINRANGE_C2global; j++){
      qinvSS_13 = AvgQinvSS[j];
      double arg_13=qinvSS_13*Rch;
      EW_13 = 1 + par[5]/(6.*pow(2,1.5))*(8.*pow(arg_13,3) - 12.*pow(arg_13,1));
      EW_13 += par[6]/(24.*pow(2,2))*(16.*pow(arg_13,4) -48.*pow(arg_13,2) + 12);
      EW_13 += par[7]/(120.*pow(2,2.5))*(32.*pow(arg_13,5) - 160.*pow(arg_13,3) + 120*arg_13);
      EW_13 += par[8]/(720.*pow(2,3))*(64.*pow(arg_13,6) - 480.*pow(arg_13,4) + 720.*pow(arg_13,2) - 120);

      for(int k=0; k<BINRANGE_C2global; k++){
	qinvSS_23 = AvgQinvSS[k];
	double arg_23=qinvSS_23*Rch;
	EW_23 = 1 + par[5]/(6.*pow(2,1.5))*(8.*pow(arg_23,3) - 12.*pow(arg_23,1));
	EW_23 += par[6]/(24.*pow(2,2))*(16.*pow(arg_23,4) -48.*pow(arg_23,2) + 12);
	EW_23 += par[7]/(120.*pow(2,2.5))*(32.*pow(arg_23,5) - 160.*pow(arg_23,3) + 120*arg_23);
	EW_23 += par[8]/(720.*pow(2,3))*(64.*pow(arg_23,6) - 480.*pow(arg_23,4) + 720.*pow(arg_23,2) - 120);
	//
	if(i==0 || j==0 || k==0) continue;
	CSS = 1 + pow(EW_12,2)*exp(-pow(Rch*qinvSS_12,2));
	CSS += pow(EW_13,2)*exp(-pow(Rch*qinvSS_13,2));
	CSS += pow(EW_23,2)*exp(-pow(Rch*qinvSS_23,2));
	CSS += 2*EW_12*EW_13*EW_23 * exp(-1/2.*pow(Rch,2)*(pow(qinvSS_12,2)+pow(qinvSS_13,2)+pow(qinvSS_23,2)));
	//
	if(IncludeSS) {
	  /*double errorSS = sqrt(pow(C3ssFitting_e[i],2));
	  if(errorSS > 0){
	    SumChi2 += pow((CSS - C3ssFitting[i])/errorSS,2);
	    NFitPoints_C3global++;
	    }*/
	}
	if(IncludeOS) {
	  //double errorOS = sqrt(pow(C3osFitting_e[i],2));
	  //SumChi2 += pow((COS - C3osFitting[i])/errorOS,2);
	  //NFitPoints_C3global++;
	}
      
	
      }// k
    }// j
  }// i 

  f = SumChi2;
  Chi2_C3global = f;

}
//______________________________________________________________________
/*void C3ssFitFunction(Double_t *x, Double_t *par){

  double qinvSS_12=0, qinvOS_12=0;
  double qinvSS_13=0, qinvOS_13=0;
  double qinvSS_23=0, qinvOS_23=0;

  double D=fabs(par[2])/(1-fabs(par[2]));
  double Rch=par[3]/FmToGeV;
  double Rcoh=par[4]/FmToGeV;
  double CSS=0, COS=0;
  double term1=0, term2=0;
  double SumChi2=0;
  double EW_12=0, EW_13=0, EW_23=0;
  NFitPoints_C3global=0;

  
  int qbin_12 = x[0]*1000/5;
  qinvSS_12 = AvgQinvSS[qbin_12];
  qinvOS_12 = AvgQinvOS[qbin_12];
  double arg_12=qinvSS_12*Rch;
  EW_12 = 1 + par[5]/(6.*pow(2,1.5))*(8.*pow(arg_12,3) - 12.*pow(arg_12,1));
  EW_12 += par[6]/(24.*pow(2,2))*(16.*pow(arg_12,4) -48.*pow(arg_12,2) + 12);
  EW_12 += par[7]/(120.*pow(2,2.5))*(32.*pow(arg_12,5) - 160.*pow(arg_12,3) + 120*arg_12);
  EW_12 += par[8]/(720.*pow(2,3))*(64.*pow(arg_12,6) - 480.*pow(arg_12,4) + 720.*pow(arg_12,2) - 120);
  //
  int qbin_13 = x[1]*1000/5;
  qinvSS_13 = AvgQinvSS[qbin_13];
  qinvOS_13 = AvgQinvOS[qbin_13];
  double arg_13=qinvSS_13*Rch;
  EW_13 = 1 + par[5]/(6.*pow(2,1.5))*(8.*pow(arg_13,3) - 12.*pow(arg_13,1));
  EW_13 += par[6]/(24.*pow(2,2))*(16.*pow(arg_13,4) -48.*pow(arg_13,2) + 12);
  EW_13 += par[7]/(120.*pow(2,2.5))*(32.*pow(arg_13,5) - 160.*pow(arg_13,3) + 120*arg_13);
  EW_13 += par[8]/(720.*pow(2,3))*(64.*pow(arg_13,6) - 480.*pow(arg_13,4) + 720.*pow(arg_13,2) - 120);
  //
  int qbin_23 = x[2]*1000/5;
  qinvSS_23 = AvgQinvSS[qbin_23];
  qinvOS_23 = AvgQinvOS[qbin_23];
  double arg_23=qinvSS_23*Rch;
  EW_23 = 1 + par[5]/(6.*pow(2,1.5))*(8.*pow(arg_23,3) - 12.*pow(arg_23,1));
  EW_23 += par[6]/(24.*pow(2,2))*(16.*pow(arg_23,4) -48.*pow(arg_23,2) + 12);
  EW_23 += par[7]/(120.*pow(2,2.5))*(32.*pow(arg_23,5) - 160.*pow(arg_23,3) + 120*arg_23);
  EW_23 += par[8]/(720.*pow(2,3))*(64.*pow(arg_23,6) - 480.*pow(arg_23,4) + 720.*pow(arg_23,2) - 120);
  //
  if(qbin_12==0 || qbin_13==0 || qbin_23==0) return 0;
  CSS = 1 + pow(EW_12,2)*exp(-pow(Rch*qinvSS_12,2));
  CSS += pow(EW_13,2)*exp(-pow(Rch*qinvSS_13,2));
  CSS += pow(EW_23,2)*exp(-pow(Rch*qinvSS_23,2));
  CSS += 2*EW_12*EW_13*EW_23 * exp(-1/2.*pow(Rch,2)*(pow(qinvSS_12,2)+pow(qinvSS_13,2)+pow(qinvSS_23,2)));
  //
  
  return CSS;
  }*/
//----------------------------------------------------------------------
void fcnOSL(int& npar, double* deriv, double& f, double par[], int flag){

  double qout=0, qside=0, qlong=0;
  double Gaus_ch=0;
  double C=0;
  //double lnL=0;
  double SumChi2=0;
  NFitPoints_OSL=0;

  for(int i=0; i<BINS_OSL-10; i++){// out, 0 to BINS_OSL-10
    for(int j=0; j<BINS_OSL-10; j++){// side, 0 to BINS_OSL-10
      for(int k=0; k<BINS_OSL-10; k++){// long, 0 to BINS_OSL-10
	
	if(B[i][j][k] < 1. ) continue;
	
	qout = (i+0.5)*0.005;
	qside = (j+0.5)*0.005;
	qlong = (k+0.5)*0.005;
	
	Gaus_ch = exp(-pow(par[3]*qout/FmToGeV,2) - pow(par[4]*qside/FmToGeV,2) - pow(par[5]*qlong/FmToGeV,2));
	
	//C = par[0]*(1 + par[1]*(K_OSL[i][j][k]-1) + K_OSL[i][j][k]*par[1]*Gaus_ch);
	C = par[0]*( (1-par[1]) + par[1]*K_OSL[i][j][k]*(1+Gaus_ch));
		
	if(C < 0) continue;
	
	double error = pow(A_e[i][j][k]/B[i][j][k],2) + pow(B_e[i][j][k]*A[i][j][k]/pow(B[i][j][k],2),2);
	error += pow((C/par[0] - (1-par[1]))/K_OSL[i][j][k] * K_OSL_e[i][j][k]*par[0],2);
	error = sqrt(error);
	SumChi2 += pow( (C - A[i][j][k]/B[i][j][k])/error,2);

	/*
	if(A[i][j][k] >= 1) term1 = C*(A[i][j][k]+B[i][j][k])/(A[i][j][k]*(C+1));
	else term1 = 0;
	term2 = (A[i][j][k]+B[i][j][k])/(B[i][j][k]*(C+1));
	
	if(term1 > 0.0 && term2 > 0.0){
	  lnL += A[i][j][k]*log(term1) + B[i][j][k]*log(term2);
	}else if(term1==0 && term2 > 0.0){
	  lnL += B[i][j][k]*log(term2);
	}else {cout<<"WARNING -- term1 and term2 are negative"<<endl;}
	*/
	NFitPoints_OSL++;
      }
    }
  }
  f = SumChi2;
  //f = -2.0*lnL;
  Chi2_OSL = f;
      
}
//________________________________________________________________________
double OSLfitfunction(Double_t *x, Double_t *par)
{
  int bin_out = int(1000*(x[0])/5.);
  int bin_side = int(1000*(x[1])/5.);
  int bin_long = int(1000*(x[2])/5.);
  
  double K = CoulCorr2(+1, avg_q[bin_out][bin_side][bin_long]);
  double Gaus_ch = exp(-pow(par[3]*x[0]/FmToGeV,2) - pow(par[4]*x[1]/FmToGeV,2) - pow(par[5]*x[2]/FmToGeV,2));
  return par[0]*(1 + par[1]*(K-1) + K*par[1]*Gaus_ch);
  
}
//__________________________________________________________________________

void fcn_3(int& npar, double* deriv, double& f, double par[], int flag){

  double q12=0, q13=0, q23=0;
  double K12=0, K13=0, K23=0, K123=0;
  double G1=0, G2=0, G3=0, G4=0;
  double EW12=0, EW13=0, EW23=0;
  double C=0;
  double term1=0, term2=0;
  double lnL=0;
  
  // start from 2-4 MeV bins
  for(int i=1; i<BINRANGE_3-10; i++){// q12
    for(int j=1; j<BINRANGE_3-10; j++){// q13
      for(int k=1; k<BINRANGE_3-10; k++){// q23
	
	if(B_3[i][j][k] < 1. ) continue;
 
	q12 = (i+0.5)*0.002;
	q13 = (j+0.5)*0.002;
	q23 = (k+0.5)*0.002;
		
	//K12 = 1/CoulCorr_CsorgoMate(2, q12);
	//K13 = 1/CoulCorr_CsorgoMate(2, q13);
	//K23 = 1/CoulCorr_CsorgoMate(2, q23);
	//K12 = 1/Coul_pipi_data(q12);
	//K13 = 1/Coul_pipi_data(q13);
	//K23 = 1/Coul_pipi_data(q23);
	K12=1.0; K13=1.0; K23=1.0;
	K123 = K12*K13*K23;
	//
	//K123 = 1 - pow(par[1],3)*(1-K12*K13*K23);
	//K12 = 1 - pow(par[1],3)*(1-K12);
	//K13 = 1 - pow(par[1],3)*(1-K13);
	//K23 = 1 - pow(par[1],3)*(1-K23);
	//
	EW12 = 1 + par[4]/(6.*pow(2,1.5))*(8*pow(q12*par[2]/FmToGeV,3) - 12*pow(q12*par[2]/FmToGeV,1)) + par[5]/(24.*pow(2,2))*(16*pow(q12*par[2]/FmToGeV,4) -48*pow(q12*par[2]/FmToGeV,2) + 12);
	EW13 = 1 + par[4]/(6.*pow(2,1.5))*(8*pow(q13*par[2]/FmToGeV,3) - 12*pow(q13*par[2]/FmToGeV,1)) + par[5]/(24.*pow(2,2))*(16*pow(q13*par[2]/FmToGeV,4) -48*pow(q13*par[2]/FmToGeV,2) + 12);
	EW23 = 1 + par[4]/(6.*pow(2,1.5))*(8*pow(q23*par[2]/FmToGeV,3) - 12*pow(q23*par[2]/FmToGeV,1)) + par[5]/(24.*pow(2,2))*(16*pow(q23*par[2]/FmToGeV,4) -48*pow(q23*par[2]/FmToGeV,2) + 12);
	//EW12=1; EW13=1; EW23=1;
	//
	G1 = K12*exp(-pow(par[2]*q12/FmToGeV,2)/2.)*EW12;
	G1 += K13*exp(-pow(par[2]*q13/FmToGeV,2)/2.)*EW13;
	G1 += K23*exp(-pow(par[2]*q23/FmToGeV,2)/2.)*EW23;
	//
	G2 = K12*exp(-pow(par[2]*q12/FmToGeV,2))*pow(EW12,2);
	G2 += K13*exp(-pow(par[2]*q13/FmToGeV,2))*pow(EW13,2);
	G2 += K23*exp(-pow(par[2]*q23/FmToGeV,2))*pow(EW23,2);
	//
	G3 = exp(-pow(par[2]/FmToGeV,2)*(q12*q12 + q23*q23)/2.)*EW12*EW23;
	G3 += exp(-pow(par[2]/FmToGeV,2)*(q13*q13 + q23*q23)/2.)*EW13*EW23; 
	G3 += exp(-pow(par[2]/FmToGeV,2)*(q12*q12 + q13*q13)/2.)*EW12*EW13;
	G3 *= K123;
	//
	G4 = exp(-pow(par[2]/FmToGeV,2)*(q12*q12 + q13*q13 + q23*q23)/2.);
	G4 *= EW12*EW13*EW23;
	G4 *= K123;
	//
	C = 1 - pow(par[1],3)*(1-K123) - pow(par[1],2)*((1-K12) + (1-K13) + (1-K23));// pure Coulomb
	//C = 1 - (1-K123) - ((1-K12) + (1-K13) + (1-K23));// pure Coulomb
	C += pow(par[1],2)*2*par[3]*(1-par[3])*G1;// 2-mixed chaos/coherence
	C += pow(par[1],2)*pow(par[3],2)*G2;// 2-pure chaos
	C += pow(par[1],3)*2*pow(par[3],2)*(1-par[3])*G3;// 3-mixed chaos/coherence
	C += pow(par[1],3)*2*pow(par[3],3)*G4;// 3-pure chaos 
	C *= par[0];// Norm
	


	if(C < 0) continue;
	
	if(A_3[i][j][k] >= 1) term1 = C*(A_3[i][j][k]+B_3[i][j][k])/(A_3[i][j][k]*(C+1));
	else term1 = 0;
	term2 = (A_3[i][j][k]+B_3[i][j][k])/(B_3[i][j][k]*(C+1));
	
	if(term1 > 0.0 && term2 > 0.0){
	  lnL += A_3[i][j][k]*log(term1) + B_3[i][j][k]*log(term2);
	}else if(term1==0 && term2 > 0.0){
	  lnL += B_3[i][j][k]*log(term2);
	}else {cout<<"WARNING -- term1 and term2 are negative"<<endl;}
	
      }
    }
  }
  
  f = -2.0*lnL;
      
}
//________________________________________________________________________
double D3fitfunction_3(Double_t *x, Double_t *par)
{
  double K12=CoulCorr2(+1, x[0]);
  double K13=CoulCorr2(+1, x[1]);
  double K23=CoulCorr2(+1, x[2]);
  K12=1.0; K13=1.0; K23=1.0;
  double K123=K12*K13*K23;
  //
  //K123 = 1 - pow(par[1],3)*(1 - K123);
  //K12 = 1 - pow(par[1],3)*(1-K12);
  //K13 = 1 - pow(par[1],3)*(1-K13);
  //K23 = 1 - pow(par[1],3)*(1-K23);
  //
  double EW12 = 1 + par[4]/(6.*pow(2,1.5))*(8*pow(x[0]*par[2]/FmToGeV,3) - 12*pow(x[0]*par[2]/FmToGeV,1)) + par[5]/(24.*pow(2,2))*(16*pow(x[0]*par[2]/FmToGeV,4) -48*pow(x[0]*par[2]/FmToGeV,2) + 12);
  double EW13 = 1 + par[4]/(6.*pow(2,1.5))*(8*pow(x[1]*par[2]/FmToGeV,3) - 12*pow(x[1]*par[2]/FmToGeV,1)) + par[5]/(24.*pow(2,2))*(16*pow(x[1]*par[2]/FmToGeV,4) -48*pow(x[1]*par[2]/FmToGeV,2) + 12);
  double EW23 = 1 + par[4]/(6.*pow(2,1.5))*(8*pow(x[2]*par[2]/FmToGeV,3) - 12*pow(x[2]*par[2]/FmToGeV,1)) + par[5]/(24.*pow(2,2))*(16*pow(x[2]*par[2]/FmToGeV,4) -48*pow(x[2]*par[2]/FmToGeV,2) + 12);
  //EW12=1; EW13=1; EW23=1;
  //
  double G1=0, G2=0, G3=0, G4=0;
  
  G1 = K12*exp(-pow(par[2]*x[0]/FmToGeV,2)/2.)*EW12;
  G1 += K13*exp(-pow(par[2]*x[1]/FmToGeV,2)/2.)*EW13;
  G1 += K23*exp(-pow(par[2]*x[2]/FmToGeV,2)/2.)*EW23;
  //
  G2 = K12*exp(-pow(par[2]*x[0]/FmToGeV,2))*pow(EW12,2);
  G2 += K13*exp(-pow(par[2]*x[1]/FmToGeV,2))*pow(EW13,2);
  G2 += K23*exp(-pow(par[2]*x[2]/FmToGeV,2))*pow(EW23,2);
  //
  G3 = exp(-pow(par[2]/FmToGeV,2)*(x[0]*x[0] + x[2]*x[2])/2.)*EW12*EW23;
  G3 += exp(-pow(par[2]/FmToGeV,2)*(x[1]*x[1] + x[2]*x[2])/2.)*EW13*EW23;
  G3 += exp(-pow(par[2]/FmToGeV,2)*(x[0]*x[0] + x[1]*x[1])/2.)*EW12*EW13;
  G3 *= K123;
  //
  G4 = exp(-pow(par[2]/FmToGeV,2)*(x[0]*x[0] + x[1]*x[1] + x[2]*x[2])/2.);
  G4 *= EW12*EW13*EW23;
  G4 *= K123;
  //
  double C = 1 - pow(par[1],3)*(1-K123) - pow(par[1],2)*((1-K12) + (1-K13) + (1-K23));// pure Coulomb
  //double C = 1 - (1-K123) - ((1-K12) + (1-K13) + (1-K23));// pure Coulomb
  C += pow(par[1],2)*2*par[3]*(1-par[3])*G1;// 2-mixed chaos/coherence
  C += pow(par[1],2)*pow(par[3],2)*G2;// 2-pure chaos
  C += pow(par[1],3)*2*pow(par[3],2)*(1-par[3])*G3;// 3-mixed chaos/coherence
  C += pow(par[1],3)*2*pow(par[3],3)*G4;// 3-pure chaos 
  C *= par[0];// Norm
  
  return C;
}

void ReadCoulCorrections_Omega0(){
  // read in 3d 3-particle coulomb correlations = K3
  
  TString *fname;
  if(FileSetting==6) fname = new TString("KFile_GaussR8to5.root");
  else if(FileSetting==7) fname = new TString("KFile_GaussR11to6.root");
  else{
    if(SourceType==0) fname = new TString("KFile.root");
    if(SourceType==1) fname = new TString("KFile_50fmCut.root");
    if(SourceType==2) fname = new TString("KFile_GaussR11to6.root");
  }
  TFile *coulfile=new TFile(fname->Data(),"READ");

  TString *name=new TString("K3ss_");
  *name += bBin-1;
  CoulCorrOmega0SS = (TH3D*)coulfile->Get(name->Data());
  CoulCorrOmega0SS->SetDirectory(0);
  name=new TString("K3os_");
  *name += bBin-1;
  CoulCorrOmega0OS = (TH3D*)coulfile->Get(name->Data());
  CoulCorrOmega0OS->SetDirectory(0);
  coulfile->Close();

}
double CoulCorrGRS(bool SC, double Q_12, double Q_13, double Q_23){
  /*int Q12bin = CoulCorr2SS->GetXaxis()->FindBin(Q_12);
  int Q13bin = CoulCorr2SS->GetXaxis()->FindBin(Q_13);
  int Q23bin = CoulCorr2SS->GetXaxis()->FindBin(Q_23);*/
  int index12L = int(fabs(Q_12 - CoulCorr2SS->GetXaxis()->GetBinWidth(1)/2.)/(CoulCorr2SS->GetXaxis()->GetBinWidth(1)));
  int index12H = index12L+1;
  int index13L = int(fabs(Q_13 - CoulCorr2SS->GetXaxis()->GetBinWidth(1)/2.)/(CoulCorr2SS->GetXaxis()->GetBinWidth(1)));
  int index13H = index13L+1;
  int index23L = int(fabs(Q_23 - CoulCorr2SS->GetXaxis()->GetBinWidth(1)/2.)/(CoulCorr2SS->GetXaxis()->GetBinWidth(1)));
  int index23H = index23L+1;

  if(Tricubic){
    // Tricubic Interpolation
    double arr[4][4][4]={{{0}}};
    for(int x=0; x<4; x++){
      for(int y=0; y<4; y++){
	for(int z=0; z<4; z++){
	  if(SC){
	    arr[x][y][z] = CoulCorr2SS->GetBinContent((index12L)+x)*CoulCorr2SS->GetBinContent((index23L)+y)*CoulCorr2SS->GetBinContent((index13L)+z);
	  }else{
	    arr[x][y][z] = CoulCorr2SS->GetBinContent((index12L)+x)*CoulCorr2OS->GetBinContent((index23L)+y)*CoulCorr2OS->GetBinContent((index13L)+z);
	  }
	  
	}
      }
    }
    return tricubicInterpolate(arr, Q_12, Q_23, Q_13);
  }else{
    // Trilinear Interpolation.  See for instance: https://en.wikipedia.org/wiki/Trilinear_interpolation
    //
    double xd = (Q_12-CoulCorr2SS->GetXaxis()->GetBinCenter(index12L+1));
    xd /= (CoulCorr2SS->GetXaxis()->GetBinCenter(index12H+1)-CoulCorr2SS->GetXaxis()->GetBinCenter(index12L+1));
    double yd = (Q_13-CoulCorr2SS->GetXaxis()->GetBinCenter(index13L+1));
    yd /= (CoulCorr2SS->GetXaxis()->GetBinCenter(index13H+1)-CoulCorr2SS->GetXaxis()->GetBinCenter(index13L+1));
    double zd = (Q_23-CoulCorr2SS->GetXaxis()->GetBinCenter(index23L+1));
    zd /= (CoulCorr2SS->GetXaxis()->GetBinCenter(index23H+1)-CoulCorr2SS->GetXaxis()->GetBinCenter(index23L+1));
    //
    if(SC){
      double c00 = CoulCorr2SS->GetBinContent(index12L+1)*CoulCorr2SS->GetBinContent(index13L+1)*CoulCorr2SS->GetBinContent(index23L+1)*(1-xd);
      c00 += CoulCorr2SS->GetBinContent(index12H+1)*CoulCorr2SS->GetBinContent(index13L+1)*CoulCorr2SS->GetBinContent(index23L+1)*xd;
      double c10 = CoulCorr2SS->GetBinContent(index12L+1)*CoulCorr2SS->GetBinContent(index13H+1)*CoulCorr2SS->GetBinContent(index23L+1)*(1-xd);
      c10 += CoulCorr2SS->GetBinContent(index12H+1)*CoulCorr2SS->GetBinContent(index13H+1)*CoulCorr2SS->GetBinContent(index23L+1)*xd;
      double c01 = CoulCorr2SS->GetBinContent(index12L+1)*CoulCorr2SS->GetBinContent(index13L+1)*CoulCorr2SS->GetBinContent(index23H+1)*(1-xd);
      c01 += CoulCorr2SS->GetBinContent(index12H+1)*CoulCorr2SS->GetBinContent(index13L+1)*CoulCorr2SS->GetBinContent(index23H+1)*xd;
      double c11 = CoulCorr2SS->GetBinContent(index12L+1)*CoulCorr2SS->GetBinContent(index13H+1)*CoulCorr2SS->GetBinContent(index23H+1)*(1-xd);
      c11 += CoulCorr2SS->GetBinContent(index12H+1)*CoulCorr2SS->GetBinContent(index13H+1)*CoulCorr2SS->GetBinContent(index23H+1)*xd;
      //
      double c0 = c00*(1-yd) + c10*yd;
      double c1 = c01*(1-yd) + c11*yd;
      return (c0*(1-zd) + c1*zd);
    }else{
      double c00 = CoulCorr2SS->GetBinContent(index12L+1)*CoulCorr2OS->GetBinContent(index13L+1)*CoulCorr2OS->GetBinContent(index23L+1)*(1-xd);
      c00 += CoulCorr2SS->GetBinContent(index12H+1)*CoulCorr2OS->GetBinContent(index13L+1)*CoulCorr2OS->GetBinContent(index23L+1)*xd;
      double c10 = CoulCorr2SS->GetBinContent(index12L+1)*CoulCorr2OS->GetBinContent(index13H+1)*CoulCorr2OS->GetBinContent(index23L+1)*(1-xd);
      c10 += CoulCorr2SS->GetBinContent(index12H+1)*CoulCorr2OS->GetBinContent(index13H+1)*CoulCorr2OS->GetBinContent(index23L+1)*xd;
      double c01 = CoulCorr2SS->GetBinContent(index12L+1)*CoulCorr2OS->GetBinContent(index13L+1)*CoulCorr2OS->GetBinContent(index23H+1)*(1-xd);
      c01 += CoulCorr2SS->GetBinContent(index12H+1)*CoulCorr2OS->GetBinContent(index13L+1)*CoulCorr2OS->GetBinContent(index23H+1)*xd;
      double c11 = CoulCorr2SS->GetBinContent(index12L+1)*CoulCorr2OS->GetBinContent(index13H+1)*CoulCorr2OS->GetBinContent(index23H+1)*(1-xd);
      c11 += CoulCorr2SS->GetBinContent(index12H+1)*CoulCorr2OS->GetBinContent(index13H+1)*CoulCorr2OS->GetBinContent(index23H+1)*xd;
      //
      double c0 = c00*(1-yd) + c10*yd;
      double c1 = c01*(1-yd) + c11*yd;
      return (c0*(1-zd) + c1*zd);
    }
  }

}
double CoulCorrOmega0(bool SC, double Q_12, double Q_13, double Q_23){// 12 is same-charge pair
  int Q12bin = CoulCorrOmega0SS->GetXaxis()->FindBin(Q_12);
  int Q13bin = CoulCorrOmega0SS->GetZaxis()->FindBin(Q_13);
  int Q23bin = CoulCorrOmega0SS->GetYaxis()->FindBin(Q_23);
  int index12L = int(fabs(Q_12 - CoulCorrOmega0SS->GetXaxis()->GetBinWidth(1)/2.)/(CoulCorrOmega0SS->GetXaxis()->GetBinWidth(1)));
  int index12H = index12L+1;
  int index13L = int(fabs(Q_13 - CoulCorrOmega0SS->GetXaxis()->GetBinWidth(1)/2.)/(CoulCorrOmega0SS->GetXaxis()->GetBinWidth(1)));
  int index13H = index13L+1;
  int index23L = int(fabs(Q_23 - CoulCorrOmega0SS->GetXaxis()->GetBinWidth(1)/2.)/(CoulCorrOmega0SS->GetXaxis()->GetBinWidth(1)));
  int index23H = index23L+1;
  

  if(SC) {if(CoulCorrOmega0SS->GetBinContent(Q12bin, Q23bin, Q13bin) <=0) {cout<<"Problematic same-charge Omega0 bin"<<endl;}}
  else {if(CoulCorrOmega0OS->GetBinContent(Q12bin, Q23bin, Q13bin) <=0) {cout<<"Problematic mixed-charge Omega0 bin"<<endl;}}
  
  if(Tricubic){
    //////////////////////////
    // Tricubic Interpolation
    double arr[4][4][4]={{{0}}};
    for(int x=0; x<4; x++){
      for(int y=0; y<4; y++){
	for(int z=0; z<4; z++){
	  if(SC) arr[x][y][z] = CoulCorrOmega0SS->GetBinContent((index12L)+x, (index23L)+y, (index13L)+z);
	  else arr[x][y][z] = CoulCorrOmega0OS->GetBinContent((Q12bin)+x, (Q23bin)+y, (Q13bin)+z);
	}
      }
    }
    return tricubicInterpolate(arr, Q_12, Q_23, Q_13);
  }else{
    ///////////////////////////
    // Trilinear Interpolation.  See for instance: https://en.wikipedia.org/wiki/Trilinear_interpolation
    //
    double xd = (Q_12-CoulCorrOmega0SS->GetXaxis()->GetBinCenter(index12L+1));
    xd /= (CoulCorrOmega0SS->GetXaxis()->GetBinCenter(index12H+1)-CoulCorrOmega0SS->GetXaxis()->GetBinCenter(index12L+1));
    double yd = (Q_23-CoulCorrOmega0SS->GetYaxis()->GetBinCenter(index23L+1));
    yd /= (CoulCorrOmega0SS->GetYaxis()->GetBinCenter(index23H+1)-CoulCorrOmega0SS->GetYaxis()->GetBinCenter(index23L+1));
    double zd = (Q_13-CoulCorrOmega0SS->GetZaxis()->GetBinCenter(index13L+1));
    zd /= (CoulCorrOmega0SS->GetZaxis()->GetBinCenter(index13H+1)-CoulCorrOmega0SS->GetZaxis()->GetBinCenter(index13L+1));
    
    
    double c00=0, c10=0, c01=0, c11=0;
    if(SC){
      // interpolate along x
      if(CoulCorrOmega0SS->GetBinContent(index12L+1, index23L+1, index13L+1)>0 && CoulCorrOmega0SS->GetBinContent(index12H+1, index23L+1, index13L+1)>0) c00 = CoulCorrOmega0SS->GetBinContent(index12L+1, index23L+1, index13L+1)*(1-xd) + CoulCorrOmega0SS->GetBinContent(index12H+1, index23L+1, index13L+1)*xd;
      if(CoulCorrOmega0SS->GetBinContent(index12L+1, index23H+1, index13L+1)>0 && CoulCorrOmega0SS->GetBinContent(index12H+1, index23H+1, index13L+1)>0) c10 = CoulCorrOmega0SS->GetBinContent(index12L+1, index23H+1, index13L+1)*(1-xd) + CoulCorrOmega0SS->GetBinContent(index12H+1, index23H+1, index13L+1)*xd;
      if(CoulCorrOmega0SS->GetBinContent(index12L+1, index23L+1, index13H+1)>0 && CoulCorrOmega0SS->GetBinContent(index12H+1, index23L+1, index13H+1)>0) c01 = CoulCorrOmega0SS->GetBinContent(index12L+1, index23L+1, index13H+1)*(1-xd) + CoulCorrOmega0SS->GetBinContent(index12H+1, index23L+1, index13H+1)*xd;
      if(CoulCorrOmega0SS->GetBinContent(index12L+1, index23H+1, index13H+1)>0 && CoulCorrOmega0SS->GetBinContent(index12H+1, index23H+1, index13H+1)>0) c11 = CoulCorrOmega0SS->GetBinContent(index12L+1, index23H+1, index13H+1)*(1-xd) + CoulCorrOmega0SS->GetBinContent(index12H+1, index23H+1, index13H+1)*xd;
    }else {
      if(CoulCorrOmega0OS->GetBinContent(index12L+1, index23L+1, index13L+1)>0 && CoulCorrOmega0OS->GetBinContent(index12H+1, index23L+1, index13L+1)>0) c00 = CoulCorrOmega0OS->GetBinContent(index12L+1, index23L+1, index13L+1)*(1-xd) + CoulCorrOmega0OS->GetBinContent(index12H+1, index23L+1, index13L+1)*xd;
      if(CoulCorrOmega0OS->GetBinContent(index12L+1, index23H+1, index13L+1)>0 && CoulCorrOmega0OS->GetBinContent(index12H+1, index23H+1, index13L+1)>0) c10 = CoulCorrOmega0OS->GetBinContent(index12L+1, index23H+1, index13L+1)*(1-xd) + CoulCorrOmega0OS->GetBinContent(index12H+1, index23H+1, index13L+1)*xd;
      if(CoulCorrOmega0OS->GetBinContent(index12L+1, index23L+1, index13H+1)>0 && CoulCorrOmega0OS->GetBinContent(index12H+1, index23L+1, index13H+1)>0) c01 = CoulCorrOmega0OS->GetBinContent(index12L+1, index23L+1, index13H+1)*(1-xd) + CoulCorrOmega0OS->GetBinContent(index12H+1, index23L+1, index13H+1)*xd;
      if(CoulCorrOmega0OS->GetBinContent(index12L+1, index23H+1, index13H+1)>0 && CoulCorrOmega0OS->GetBinContent(index12H+1, index23H+1, index13H+1)>0) c11 = CoulCorrOmega0OS->GetBinContent(index12L+1, index23H+1, index13H+1)*(1-xd) + CoulCorrOmega0OS->GetBinContent(index12H+1, index23H+1, index13H+1)*xd;
    }
    //
    double c0=0, c1=0;
    if(c00>0 && c10>0) c0 = c00*(1-yd) + c10*yd;
    if(c01>0 && c11>0) c1 = c01*(1-yd) + c11*yd;
    
    if(c0>0 && c1>0) {
      double valueOm = (c0*(1-zd) + c1*zd);
      if(SC) valueOm *= StrongSC->Eval(Q_12*1000/2.) * StrongSC->Eval(Q_23*1000/2.) * StrongSC->Eval(Q_13*1000/2.);
      else valueOm *= StrongSC->Eval(Q_12*1000/2.);
      return valueOm;
    }else {
      cout<<"Not all Omega0 bins well defined.  Use GRS instead"<<endl;
      return CoulCorrGRS(kTRUE, Q_12, Q_13, Q_23);// if cell not well defined use GRS
    }
  }
}

double Gamov(int chargeProduct, double qinv){
  
  double arg = chargeProduct*2.*PI/(BohrR*qinv/2.);
  
  return arg/(exp(arg)-1);
}
/*double LednickyCorr(double qinv){
  
  int indexL = int(fabs(1000*qinv-Lednicky_qinv[0])/(2*2.027027));// starts from 2.027 MeV
  int indexH = indexL+1;
  if( (indexL >= 73)) return 1.0;
  double slope = (Lednicky_CoulStrong[indexL] - Lednicky_CoulStrong[indexH])/(Lednicky_qinv[indexL]-Lednicky_qinv[indexH]);
  return (slope*(1000*qinv - Lednicky_qinv[indexL]) + Lednicky_CoulStrong[indexL]);
  
}
void ReadLednickyFile(int Rvalue){
  
  TString *filename=new TString("../Strong_FSI/converted_CoulStrong_");
  *filename += Rvalue;
  filename->Append("fm.root");
  TFile *Ledfile = new TFile(filename->Data(),"READ");
  TH1F *LedHisto = (TH1F*)Ledfile->Get("h27");
  for(int i=0; i<74; i++){
    Lednicky_qinv[i] = 2.*(LedHisto->GetXaxis()->GetBinLowEdge(i+1) + LedHisto->GetXaxis()->GetBinUpEdge(i+1))/2.;
    Lednicky_CoulStrong[i] = LedHisto->GetBinContent(i+1);
  }
  Ledfile->Close();
  }*/
void ReadMomResFile(int rvalue, double lambda){
 
  for(int i=0; i<40; i++){
    MomRes_C2_pp[i] = 1.;
    MomRes_C2_mp[i] = 1.;
    MomRes_term1_pp[i] = 1.;
    MomRes_term2_pp[i] = 1.;
    MomRes_term1_mp[i] = 1.;
    MomRes_term2_mp[i] = 1.;
    //
    MomResDev_C2_pp[i] = 1.;
    MomResDev_C2_mp[i] = 1.;
    MomResDev_term1_pp[i] = 1.;
    MomResDev_term2_pp[i] = 1.;
  }
 
  TString *momresfilename = new TString("MomResFile");
  momresfilename->Append("_Offline_");
  if(FileSetting==3) momresfilename->Append("PID1p5");
  else if(FileSetting==9) momresfilename->Append("FB5and7overlap");
  else momresfilename->Append("11h");
  momresfilename->Append(".root");// no corresponding file for 10h
  
  TFile *MomResFile = new TFile(momresfilename->Data(),"READ");
  
  TH2D *MomResWeights_pp = (TH2D*)MomResFile->Get("MomResHisto_pp");
  TH2D *MomResWeights_mp = (TH2D*)MomResFile->Get("MomResHisto_mp");
  TH2D *MomResWeights_pp_term1 = (TH2D*)MomResFile->Get("MomResHisto_pp_term1");
  TH2D *MomResWeights_pp_term2 = (TH2D*)MomResFile->Get("MomResHisto_pp_term2");
  TH2D *MomResWeights_mp_term1 = (TH2D*)MomResFile->Get("MomResHisto_mp_term1");
  TH2D *MomResWeights_mp_term2 = (TH2D*)MomResFile->Get("MomResHisto_mp_term2");
  //
  //
  TString *names3[2][5];// SC/MC, term#
  TString *names1D[2][5];// SC/MC, term#
  TString *names3_AvgK3[2];// SC/MC
  for(int ChProd=0; ChProd<2; ChProd++){
    if(ChProd==0) names3_AvgK3[ChProd] = new TString("AvgK3_SC_M");
    else names3_AvgK3[ChProd] = new TString("AvgK3_MC_M");
    *names3_AvgK3[ChProd] += MomResCentBin-1;
    //AvgK3[ChProd] = (TH3D*)MomResFile->Get(names3_AvgK3[ChProd]->Data());
    //AvgK3[ChProd]->SetDirectory(0);
    //
    for(int term=0; term<5; term++){
      //
      if(ChProd==0) {names3[ChProd][term] = new TString("MomResHisto3_SC_term"); names1D[ChProd][term] = new TString("MomResHisto1D_SC_term");}
      else {names3[ChProd][term] = new TString("MomResHisto3_MC_term"); names1D[ChProd][term] = new TString("MomResHisto1D_MC_term");}
      *names3[ChProd][term] += term+1;
      *names1D[ChProd][term] += term+1;
      names3[ChProd][term]->Append("_M");
      names1D[ChProd][term]->Append("_M");
      *names3[ChProd][term] += MomResCentBin-1;
      *names1D[ChProd][term] += MomResCentBin-1;
      MomRes3d[ChProd][term] = (TH3D*)MomResFile->Get(names3[ChProd][term]->Data());
      MomRes1d[ChProd][term] = (TH1D*)MomResFile->Get(names1D[ChProd][term]->Data());
      MomRes3d[ChProd][term]->SetDirectory(0);
      MomRes1d[ChProd][term]->SetDirectory(0);

    }
  }
  
  
  Int_t LambdaRbin;
  if(rvalue<=5) LambdaRbin = 1 + int(float(lambda-0.5+0.0001)/0.02);
  else LambdaRbin = 1 + 16*(rvalue-5) + int(float(lambda-0.5+0.0001)/0.02);
  Int_t LambdaRbinDev;
  if(rvalue<=5) LambdaRbinDev = 1 + 16*(1) + int(float((lambda+0.04)-0.5+0.0001)/0.02);
  else LambdaRbinDev = 1 + 16*(rvalue-5 - 1) + int(float((lambda+0.04)-0.5+0.0001)/0.02);

  for(int i=0; i<40; i++){
    MomRes_C2_pp[i] = MomResWeights_pp->GetBinContent(LambdaRbin, i+1);
    MomRes_C2_mp[i] = MomResWeights_mp->GetBinContent(LambdaRbin, i+1);
    MomRes_term1_pp[i] = MomResWeights_pp_term1->GetBinContent(LambdaRbin, i+1);
    MomRes_term2_pp[i] = MomResWeights_pp_term2->GetBinContent(LambdaRbin, i+1);
    MomRes_term1_mp[i] = MomResWeights_mp_term1->GetBinContent(LambdaRbin, i+1);
    MomRes_term2_mp[i] = MomResWeights_mp_term2->GetBinContent(LambdaRbin, i+1);
    //
    MomResDev_C2_pp[i] = MomResWeights_pp->GetBinContent(LambdaRbinDev, i+1);
    MomResDev_C2_mp[i] = MomResWeights_mp->GetBinContent(LambdaRbinDev, i+1);
    MomResDev_term1_pp[i] = MomResWeights_pp_term1->GetBinContent(LambdaRbinDev, i+1);
    MomResDev_term2_pp[i] = MomResWeights_pp_term2->GetBinContent(LambdaRbinDev, i+1);
  }

  MomResFile->Close();
  
}

void fcn_r3_3d(int& npar, double* deriv, double& f, double par[], int flag){
  
  double r3=0;
  double Chi2sum=0;
  
  for(int i=fitlimitLowBin-1; i<fitlimitHighBin; i++){
    double q12=(i+0.5)*0.005;
    for(int j=fitlimitLowBin-1; j<fitlimitHighBin; j++){
      double q13=(j+0.5)*0.005;
      for(int k=fitlimitLowBin-1; k<fitlimitHighBin; k++){
	double q23=(k+0.5)*0.005;
	
	//r3 = 1 - par[2]*pow(q12*q13 + q12*q23 + q13*q23,2) + 1/5.*pow(par[0],2)*(3-par[1]*(q12*q12+q13*q13+q23*q23)) - 112/70.*pow(par[0],3);
	//r3 /= pow( (1-par[1]*q12*q12 - 4/5.*pow(par[0],2))*(1-par[1]*q13*q13 - 4/5.*pow(par[0],2))*(1-par[1]*q23*q23 - 4/5.*pow(par[0],2)),0.5);
	r3 = par[0] - par[1]*(q12*q13 + q12*q23 + q13*q23) - par[2]*pow(q12*q13 + q12*q23 + q13*q23,2);
	if(r3_3d_arr_e[i][j][k] > 0) Chi2sum += pow((r3_3d_arr[i][j][k] - r3)/r3_3d_arr_e[i][j][k],2);
	
      }
    }
  }
  
  f = Chi2sum; 
}

double r3_fitfunction(Double_t *x, Double_t *par){
  
  //double r3 = 1 - par[2]*pow(x[0]*x[1] + x[0]*x[2] + x[1]*x[2],2) + 1/5.*pow(par[0],2)*(3-par[1]*(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])) - 112/70.*pow(par[0],3);
  //r3 /= pow( (1-par[1]*x[0]*x[0] - 4/5.*pow(par[0],2))*(1-par[1]*x[1]*x[1] - 4/5.*pow(par[0],2))*(1-par[1]*x[2]*x[2] - 4/5.*pow(par[0],2)),0.5);
  double r3 = par[0] - par[1]*(x[0]*x[1] + x[0]*x[2] + x[1]*x[2]) - par[2]*pow(x[0]*x[1] + x[0]*x[2] + x[1]*x[2],2);
  
  return r3;
  
}

double cubicInterpolate (double p[4], double x) {
  return p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));// Paulinternet
}
double bicubicInterpolate (double p[4][4], double x, double y) {
	double arr[4];
	arr[0] = cubicInterpolate(p[0], y);
	arr[1] = cubicInterpolate(p[1], y);
	arr[2] = cubicInterpolate(p[2], y);
	arr[3] = cubicInterpolate(p[3], y);
	return cubicInterpolate(arr, x);
}
double tricubicInterpolate (double p[4][4][4], double x, double y, double z) {
	double arr[4];
	arr[0] = bicubicInterpolate(p[0], y, z);
	arr[1] = bicubicInterpolate(p[1], y, z);
	arr[2] = bicubicInterpolate(p[2], y, z);
	arr[3] = bicubicInterpolate(p[3], y, z);
	return cubicInterpolate(arr, x);
}
float fact(float n){
  return (n < 1.00001) ? 1 : fact(n - 1) * n;
}
void DrawALICELogo(Bool_t prel, Float_t x1, Float_t y1, Float_t x2, Float_t y2)
{
  // revision on July 24th or 25th 2012
  // correct for aspect ratio of figure plus aspect ratio of pad (coordinates are NDC!)
  x2 = x1 + (y2 - y1) * (466. / 523) * gPad->GetWh() * gPad->GetHNDC() / (gPad->GetWNDC() * gPad->GetWw());
  // Printf("%f %f %f %f", x1, x2, y1, y2);
  
  TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo", x1, y1, x2, y2);
  myPadLogo->SetLeftMargin(0);
  myPadLogo->SetTopMargin(0);
  myPadLogo->SetRightMargin(0);
  myPadLogo->SetBottomMargin(0);
  myPadLogo->Draw();
  myPadLogo->cd();
  TASImage *myAliceLogo = new TASImage((prel) ? "~/Pictures/2011-Nov-24-ALICE_PRELIMARY_logo_BLACK_small_usage_design.eps" : "~/Pictures/2011-Nov-24-ALICE_PERFORMANCE_logo_BLACK_small_usage_design.eps");
  myAliceLogo->Draw();
}
//________________________________________________________________________________________
void ReadMuonCorrections(int mbin){
  TString *name = new TString("MuonCorrection_");
  if(mbin<=1) *name += 0;
  else if(mbin<=3) *name += 1;
  else if(mbin<=5) *name += 2;
  else if(mbin<=7) *name += 3;
  else if(mbin<=9) *name += 3;
  else *name += 4;
  name->Append(".root");
  TFile *muonfile=new TFile(name->Data(),"READ");
  C2muonCorrection_SC = (TH1D*)muonfile->Get("C2muonCorrection_SC");
  C2muonCorrection_MC = (TH1D*)muonfile->Get("C2muonCorrection_MC");
  WeightmuonCorrection = (TH1D*)muonfile->Get("WeightmuonCorrection");
  if(SameCharge_def) C3muonCorrection = (TH1D*)muonfile->Get("C3muonCorrection_SC");
  else C3muonCorrection = (TH1D*)muonfile->Get("C3muonCorrection_MC");
  //
  C2muonCorrection_SC->SetDirectory(0);
  C2muonCorrection_MC->SetDirectory(0);
  WeightmuonCorrection->SetDirectory(0);
  C3muonCorrection->SetDirectory(0);
  //
  muonfile->Close();
}
