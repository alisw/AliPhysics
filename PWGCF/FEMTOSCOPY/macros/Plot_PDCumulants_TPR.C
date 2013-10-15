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

#define BohrR 1963.6885 // Bohr Radius for pions
#define FmToGeV 0.19733 // conversion to fm
#define PI 3.1415926
#define masspiC 0.1395702 // pi+ mass (GeV/c^2)
#define kappa3 0.12
#define kappa4 0.43

using namespace std;

int CollisionType_def=2;// 0=PbPb, 1=pPb, 2=pp
bool MCcase_def=kTRUE;// MC data?
int CHARGE_def=-1;// -1 or +1: + or - pions for same-charge case, --+ or ++- for mixed-charge case
bool SameCharge_def=kFALSE;// 3-pion same-charge?
bool AddCC=kTRUE;
//
bool IncludeEW_def=kFALSE;// Include EdgeWorth coefficients?
bool FixEWavg=kTRUE;
int Mbin_def=18;// 0-19 (0=1050-2000 pions, 19=0-5 pions)
int EDbin_def=0;// 0-2: Kt3 bin
int Ktbin_def=3;// 1-6, 10=full range
double MRCShift=1.0;// 1.0=full Momentum Resolution Correction. 1.1 for 10% systematic increase
//
//
bool SaveToFile_def=kFALSE;// Save outputs to file?
int SourceType=0;// 0=Therminator, 1=Therminator with 50fm cut (keep at 0 for default), 2=Gaussian
bool ConstantFSI=kFALSE;// Constant FSI's for each kt bin?
bool GofP=kFALSE;// Include momentum dependence of coherent fraction?
bool ChargeConstraint=kFALSE;// Include Charge Constraint for coherent states?
bool LinkN=kFALSE;// Make N(++)=N(+-)?
bool GeneratedSignal=kFALSE;// Was the QS+FSI signal generated in MC? 
bool Tricubic=kFALSE;// Tricubic or Trilinear interpolation?  Tricubic was found to be unstable.
bool IncludeSS=kTRUE;// Include same-charge two-pion correlations in the fit?
bool IncludeOS=kFALSE;// Include mixed-charge two-pion correlations in the fit?
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
int FSIindex;
int Ktlowbin;
int Kthighbin;
int MomResCentBin;// momentum resolution centrality bin
float TwoFrac;// Lambda parameter
float OneFrac;// Lambda^{1/2}
float ThreeFrac;// Lambda^{3/2}
float Q2Limit;
float Q3Limit;

const int Nlines = 50;
TH1D *CoulCorr2SS;
TH1D *CoulCorr2OS;
TH1D *CoulCorr2SS_2;// for interpolation in transition from Therminator to Gauss K factors
TH1D *CoulCorr2OS_2;// for interpolation in transition from Therminator to Gauss K factors
//

void ReadCoulCorrections(int);
void ReadMomResFile(int);
double CoulCorr2(int, double);
double CoulCorrGRS(bool, double, double, double);
double Dfitfunction_c3(Double_t *x, Double_t *par);
double Gamov(int, double);
double C2ssFitFunction(Double_t *x, Double_t *par);
double C2osFitFunction(Double_t *x, Double_t *par);
double cubicInterpolate(double[4], double);
double bicubicInterpolate(double[4][4], double, double);
double tricubicInterpolate(double[4][4][4], double, double, double);
//
void fcnC2_Global(int&, double*, double&, double[], int);

const int BINRANGE_C2global=40;
double C2ssFitting[BINRANGE_C2global];
double C2osFitting[BINRANGE_C2global];
double C2ssFitting_e[BINRANGE_C2global];
double C2osFitting_e[BINRANGE_C2global];
double K2SS[BINRANGE_C2global];
double K2OS[BINRANGE_C2global];
double K2SS_e[BINRANGE_C2global];
double K2OS_e[BINRANGE_C2global];
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

const int BINRANGE_3=40;// q12,q13,q23
int BINLIMIT_3=20;
double A_3[BINRANGE_3][BINRANGE_3][BINRANGE_3];// q12,q13,q23
double A_3_e[BINRANGE_3][BINRANGE_3][BINRANGE_3];// q12,q13,q23
double B_3[BINRANGE_3][BINRANGE_3][BINRANGE_3];// q12,q13,q23
double BinCenters[400];
double BinWidthQ2;
double Chi2_c3;
double NFitPoints_c3;
void fcn_c3(int&, double*, double&, double[], int);

double C3_base[BINRANGE_3][BINRANGE_3][BINRANGE_3];
double N_term1[BINRANGE_3][BINRANGE_3][BINRANGE_3];
double N_term2[BINRANGE_3][BINRANGE_3][BINRANGE_3];
double N_term3[BINRANGE_3][BINRANGE_3][BINRANGE_3];
double N_term4[BINRANGE_3][BINRANGE_3][BINRANGE_3];
double N_term5[BINRANGE_3][BINRANGE_3][BINRANGE_3];


TH3D *MomRes3d[2][5];// SC/MC, term#
TH1D *MomRes1d[2][5];// SC/MC, term#
TH1D *MomResC2[2];// SC/MC
int MRC2index;// index for C2 MRC 

double AvgP[6][20];// kt bin, qinv bin



TF1 *StrongSC;// same-charge pion strong FSI
TF1 *MixedChargeSysFit;// mixed-charge 3-pion cumulant residue obtained from Plot_plotsTPR.C
//


void Plot_PDCumulants_TPR(bool SaveToFile=SaveToFile_def, int CollisionType=CollisionType_def, bool MCcase=MCcase_def, bool SameCharge=SameCharge_def, bool IncludeEW=IncludeEW_def, int EDbin=EDbin_def, int CHARGE=CHARGE_def, int Mbin=Mbin_def, int Ktbin=Ktbin_def){
  
  EDbin_def=EDbin;
  Ktbin_def=Ktbin;
  CollisionType_def=CollisionType;
  SaveToFile_def=SaveToFile;
  MCcase_def=MCcase;
  CHARGE_def=CHARGE;
  IncludeEW_def=IncludeEW;
  SameCharge_def=SameCharge;// 3-pion same-charge
  Mbin_def=Mbin;
  //
 

  if(CollisionType==0){
    if(Mbin==0) FSIindex = 0;//0-5%
    else if(Mbin==1) FSIindex = 1;//5-10%
    else if(Mbin<=3) FSIindex = 2;//10-20%
    else if(Mbin<=5) FSIindex = 3;//20-30%
    else if(Mbin<=7) FSIindex = 4;//30-40%
    else if(Mbin<=9) FSIindex = 5;//40-50%
    else if(Mbin<=12) FSIindex = 6;//40-50%
    else if(Mbin<=15) FSIindex = 7;//6, 40-50%
    else if(Mbin<=18) FSIindex = 8;//7, 40-50%
    else FSIindex = 8;// 8, 90-100%
  }else FSIindex = 9;// 9, pp and pPb
  
  if(Mbin<=2) MRC2index=138;
  else if(Mbin<=5) MRC2index=106;
  else if(Mbin<=9) MRC2index=74;
  else if(Mbin<=9) MRC2index=42;
  else MRC2index=10;
  

  
  TwoFrac = 0.7;
  OneFrac = sqrt(TwoFrac);
  ThreeFrac = pow(TwoFrac,3/2.);
  
  if(CollisionType==0 && Mbin<10) BINLIMIT_3=20;
  else BINLIMIT_3=40;

  Ktlowbin=(Ktbin)*2+3;// kt bins are 0.5 GeV/c wide (0-0.5, 0.5-1.0 ...)
  Kthighbin=(Ktbin)*2+4;
  
  

  // Core-Halo modeling of single-particle and triple-particle core fraction 
  //float AvgN[10]={1.99266, 3.97789, 6.4624, 8.94042, 11.4194, 13.8987, 16.385, 18.8756, 21.3691, 26.0742};// 13c (avg total pion mult)/2.  2.0sigma
 
  //
  
  // bin centers from QS+FSI
  double BinCenterPbPbCentral[40]={0.00206385, 0.00818515, 0.0129022, 0.0177584, 0.0226881, 0.027647, 0.032622, 0.0376015, 0.042588, 0.0475767, 0.0525692, 0.0575625, 0.0625569, 0.0675511, 0.0725471, 0.0775436, 0.0825399, 0.0875364, 0.0925339, 0.0975321, 0.102529, 0.107527, 0.112525, 0.117523, 0.122522, 0.12752, 0.132519, 0.137518, 0.142516, 0.147515, 0.152514, 0.157513, 0.162513, 0.167512, 0.172511, 0.177511, 0.18251, 0.187509, 0.192509, 0.197509};
  double BinCenterPbPbPeripheral[40]={0.00206385, 0.00818515, 0.0129022, 0.0177584, 0.0226881, 0.027647, 0.032622, 0.0376015, 0.042588, 0.0475767, 0.0525692, 0.0575625, 0.0625569, 0.0675511, 0.0725471, 0.0775436, 0.0825399, 0.0875364, 0.0925339, 0.0975321, 0.102529, 0.107527, 0.112525, 0.117523, 0.122522, 0.12752, 0.132519, 0.137518, 0.142516, 0.147515, 0.152514, 0.157513, 0.162513, 0.167512, 0.172511, 0.177511, 0.18251, 0.187509, 0.192509, 0.197509};
  double BinCenterpPbAndpp[40]={0.00359275, 0.016357, 0.0257109, 0.035445, 0.045297, 0.0552251, 0.0651888, 0.0751397, 0.0851088, 0.0951108, 0.105084, 0.115079, 0.12507, 0.135059, 0.145053, 0.155049, 0.16505, 0.175038, 0.185039, 0.195034, 0.205023, 0.215027, 0.225024, 0.235023, 0.245011, 0.255017, 0.265017, 0.275021, 0.285021, 0.295017, 0.305018, 0.315018, 0.325013, 0.335011, 0.345016, 0.355019, 0.365012, 0.375016, 0.385017, 0.395016};
  if(CollisionType==0){
    Q3Limit = 0.1 + 0.1*Mbin/16.;
    for(int i=0; i<40; i++){
      if(Mbin<10) BinCenters[i] = BinCenterPbPbCentral[i];
      else BinCenters[i] = BinCenterPbPbPeripheral[i];
    }
    BinWidthQ2 = 0.005;
  }else{
    Q3Limit = 0.3 + 0.2*fabs(Mbin-10)/9.;
    for(int i=0; i<40; i++){
      BinCenters[i] = BinCenterpPbAndpp[i];
    }
    BinWidthQ2 = 0.01;
  }
  Q2Limit = Q3Limit/sqrt(2.);
  //Q2Limit -= 0.02;// Systematic variations
  
  // extend BinCenters for high q
  for(int index=40; index<400; index++){
    if(CollisionType==0) BinCenters[index] = (index+0.5)*(0.005);
    else BinCenters[index] = (index+0.5)*(0.010);
  }
  

  // same-charge pion strong FSI (obtained from ratio of CoulStrong to Coul at R=7 Gaussian from Lednicky's code)
  StrongSC=new TF1("StrongSC","[0]+[1]*exp(-pow([2]*x,2))",0,1000);
  StrongSC->FixParameter(0,9.99192e-01);
  StrongSC->FixParameter(1,-8.01113e-03);
  StrongSC->FixParameter(2,5.35912e-02);
  // mixed-charge 3-pion cumulant residue obtained from Plot_plotsTPR.C
  MixedChargeSysFit=new TF1("MixedChargeSysFit","[0]+[1]*exp(-pow([2]*x/0.19733,2))",0,.5);
  MixedChargeSysFit->FixParameter(0,1); 
  MixedChargeSysFit->FixParameter(1,.06);
  MixedChargeSysFit->FixParameter(2,5.5 - 6*(Mbin/19.));
  
  
  gStyle->SetOptStat(0);
  gStyle->SetOptDate(0);
  gStyle->SetOptFit(0111);

  TFile *myfile;
  if(CollisionType==0){// PbPb
    if(MCcase) {
      if(Mbin<=1){
	myfile = new TFile("Results/PDC_HIJING_12a17ad_fix_genSignal_Rinv11.root","READ");
      }else{
	myfile = new TFile("Results/PDC_HIJING_12a17b_myRun_L0p68R11_C2plots.root","READ");
      }
    }else{
      //myfile = new TFile("Results/PDC_10h_11h_0to50_50to100.root","READ"); 
      myfile = new TFile("Results/PDC_10h_11h_0to50_50to100_3Ktbins.root","READ");
    }
  }else if(CollisionType==1){// pPb
    if(!MCcase){
      //myfile = new TFile("Results/PDC_13bc.root","READ");
      myfile = new TFile("Results/PDC_13bc_3Ktbins.root","READ");
    }else{
      myfile = new TFile("Results/PDC_13b2_efix_p1_R2.root","READ");
    }
  }else{// pp
    if(!MCcase){
      //myfile = new TFile("Results/PDC_10de.root","READ");
      myfile = new TFile("Results/PDC_10cde_3Ktbins.root","READ");
    }else{
      //myfile = new TFile("Results/PDC_10f6a_R2.root","READ");
      myfile = new TFile("Results/PDC_10f6a_R2_3Ktbins.root","READ");
    }
  }
  
  ReadCoulCorrections(FSIindex);
  ReadMomResFile(Mbin);
  //
  /////////////////////////////////////////////////////////


  double NormQcutLow;
  double NormQcutHigh;
  if(CollisionType==0) {// PbPb
    if(Mbin<10){
      NormQcutLow = 0.15;
      NormQcutHigh = 0.175;
    }else{
      NormQcutLow = 0.3;
      NormQcutHigh = 0.35;
    }
  }else{// pPb and pp
    NormQcutLow = 1.0;
    NormQcutHigh = 1.2;
  }
  
  
  TList *MyList;
  TDirectoryFile *tdir;
  if(!MCcase) tdir = (TDirectoryFile*)myfile->Get("PWGCF.outputThreePionRadiiAnalysis.root");
  
  if(CollisionType==0){
    if(!MCcase){
      if(Mbin<6) MyList=(TList*)tdir->Get("ThreePionRadiiOutput_1");
      else MyList=(TList*)tdir->Get("ThreePionRadiiOutput_2");
    }else{
      MyList=(TList*)myfile->Get("MyList");
    }
  }else {
    if(!MCcase) MyList=(TList*)tdir->Get("ThreePionRadiiOutput_1");
    else MyList=(TList*)myfile->Get("MyList");
  }
  myfile->Close();
  
  TH1D *Events = (TH1D*)MyList->FindObject("fEvents2");
  //

  cout<<"#Events = "<<Events->Integral(Mbin+1,Mbin+1)<<endl;

  
  //TH1D *ChiSquaredNDF = new TH1D("ChiSquaredNDF","",2,0.5,2.5);// Chi^2/NDF records

  // Explicit Loop Histos
  TH2D *Two_ex_2d[2][2][Sbins_2][2];
  TH1D *Two_ex[2][2][Sbins_2][2];
      
  // Normalizations
  double NormH_pc[5]={0};
  
  double norm_ex_2[6][2]={{0}};
  
  // 3d method histograms
  TH3D *Three_3d[2][2][2][Sbins_3][5];
  TH1D *Three_1d[2][2][2][Sbins_3][5];
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
	  
	  if(Ktbin==10) {Ktlowbin=1; Kthighbin=Two_ex_2d[c1][c2][sc][term]->GetNbinsX();}//full kt interval
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
		  TString *name_Q3 = new TString(name3_pc->Data());
		  name_Q3->Append("_Q3");
		  Three_1d[c1][c2][c3][sc][term] = (TH1D*)MyList->FindObject(name_Q3->Data());
		  Three_1d[c1][c2][c3][sc][term]->Sumw2();
		  Three_1d[c1][c2][c3][sc][term]->Scale(NormH_pc[0]/NormH_pc[term]);
		  //cout<<"Scale factor = "<<NormH_pc[0]/NormH_pc[term]<<endl;

		  
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
  //cout<<Two_ex[0][0][0][0]->Integral(Two_ex[0][0][0][0]->GetXaxis()->FindBin(0), Two_ex[0][0][0][0]->GetXaxis()->FindBin(0.1))<<endl;
  //cout<<Two_ex[0][0][0][0]->Integral(Two_ex[0][0][0][0]->GetXaxis()->FindBin(1.06), Two_ex[0][0][0][0]->GetXaxis()->FindBin(1.1))<<endl;
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
  Two_ex_clone_mm->GetXaxis()->SetRangeUser(0,0.6);
  Two_ex_clone_mm->SetMinimum(0.95);
  Two_ex_clone_mm->SetMaximum(1.4);
  
  if(MCcase){
    Two_ex_clone_mp->SetMarkerColor(4);
    Two_ex_clone_mp->SetLineColor(4);
    Two_ex_clone_mm->Add(Two_ex_clone_pp);
    Two_ex_clone_mm->Scale(0.5);
    Two_ex_clone_mm->GetYaxis()->SetTitle("C_{2}^{MC}");
    Two_ex_clone_mm->GetYaxis()->SetTitleOffset(1.3);
    Two_ex_clone_mm->SetMinimum(0.95);
    Two_ex_clone_mm->SetMaximum(1.3);
    Two_ex_clone_mm->DrawCopy();
    //Two_ex_clone_pp->DrawCopy("same");
    Two_ex_clone_mp->DrawCopy("same");
    legend1->AddEntry(Two_ex_clone_mm,"same-charge","p");
    //legend1->AddEntry(Two_ex_clone_pp,"++","p");
    legend1->AddEntry(Two_ex_clone_mp,"mixed-charge","p");
    legend1->Draw("same");
  }


 

  
  /////////////////////////////////////////////////////
  // Global fitting C2os and C2ss
  double C2ssSys_e[BINRANGE_C2global]={0};
  double C2osSys_e[BINRANGE_C2global]={0};
  //
  const int npar=10;// 10 
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
  
  for(int i=0; i<BINRANGE_C2global; i++){
    if(i >= temp_mm->GetNbinsX()) continue;// bin limit
    
    C2ssFitting[i] = (temp_mm->GetBinContent(i+1) + temp_pp->GetBinContent(i+1))/2.;
    double MRC_SC = 1.0;
    double MRC_MC = 1.0;
    if(MomResC2[0]->GetBinContent(i+1) > 0) MRC_SC = MomResC2[0]->GetBinContent(i+1);
    if(MomResC2[1]->GetBinContent(i+1) > 0) MRC_MC = MomResC2[1]->GetBinContent(i+1);
    if(!GeneratedSignal && !MCcase) C2ssFitting[i] *= MRC_SC;
    C2ssFitting_e[i] = pow(MRC_SC*sqrt(pow(temp_mm->GetBinError(i+1),2) + pow(temp_pp->GetBinError(i+1),2))/sqrt(2.),2);
    C2ssRaw->SetBinContent(i+1, (temp_mm->GetBinContent(i+1) + temp_pp->GetBinContent(i+1))/2.);
    C2ssRaw->SetBinError(i+1, pow(sqrt(pow(temp_mm->GetBinError(i+1),2) + pow(temp_pp->GetBinError(i+1),2))/sqrt(2.),2));
    //C2ssFitting_e[i] += pow((MRC_SC-1)*0.1 * (temp_mm->GetBinContent(i+1) + temp_pp->GetBinContent(i+1))/2.,2);
    C2ssFitting_e[i] = sqrt(C2ssFitting_e[i]);
    C2osFitting[i] = temp_mp->GetBinContent(i+1);
    if(!GeneratedSignal && !MCcase) C2osFitting[i] *= MRC_MC;
    C2osFitting_e[i] = pow(MRC_MC*temp_mp->GetBinError(i+1),2);
    C2osRaw->SetBinContent(i+1, temp_mp->GetBinContent(i+1));
    C2osRaw->SetBinError(i+1, temp_mp->GetBinError(i+1));
    C2osFitting_e[i] += pow((MRC_MC-1)*0.1 * temp_mp->GetBinContent(i+1),2);
    C2osFitting_e[i] = sqrt(C2osFitting_e[i]);
    //
    C2ssSys_e[i] = pow((MRC_SC-1)*0.1 * (temp_mm->GetBinContent(i+1) + temp_pp->GetBinContent(i+1))/2.,2);
    C2ssSys_e[i] = sqrt(C2ssSys_e[i]);
    C2osSys_e[i] = pow((MRC_MC-1)*0.1 * temp_mp->GetBinContent(i+1),2);
    C2osSys_e[i] = sqrt(C2osSys_e[i]);
    //
    K2SS[i] = CoulCorr2(+1, BinCenters[i]);
    K2OS[i] = CoulCorr2(-1, BinCenters[i]);
    //K2SS[i] = 1;
    //K2OS[i] = 1;
    //
    //K2SS_e[i] = sqrt(pow((K2SS[i]-1)*0.02,2) + pow((K2SS[i]-Gamov(+1, BinCenters[i]))*0.02,2));
    //K2OS_e[i] = sqrt(pow((K2OS[i]-1)*0.02,2) + pow((K2OS[i]-Gamov(-1, BinCenters[i]))*0.02,2));
    K2SS_e[i] = 0.0;
    K2OS_e[i] = 0.0;
  }
    
    
    
    
    
  C2ssFitting[0]=-1000; C2osFitting[0]=-1000;
  C2ssFitting_e[0]=1000; C2osFitting_e[0]=1000;
  K2SS_e[0]=1000; K2OS_e[0]=1000;
  
  
  
  par[0] = 1.0; par[1]=0.5; par[2]=0.5; par[3]=9.2; par[4] = .1; par[5] = .2; par[6] = .0; par[7] = 0.; par[8] = 0.; par[9] = 0.;
  stepSize[0] = 0.01; stepSize[1] = 0.01;  stepSize[2] = 0.02; stepSize[3] = 0.2; stepSize[4] = 0.01; stepSize[5] = 0.001; stepSize[6] = 0.001; stepSize[7] = 0.001; stepSize[8]=0.001; stepSize[9]=0.01;
  minVal[0] = 0.995; minVal[1] = 0.2; minVal[2] = 0.; minVal[3] = 0.1; minVal[4] = 0.001; minVal[5] = -10.; minVal[6] = -10.; minVal[7] = -10.; minVal[8]=-10; minVal[9] = 0.995;
  maxVal[0] = 1.1; maxVal[1] = 1.0; maxVal[2] = 0.99; maxVal[3] = 15.; maxVal[4] = 2.; maxVal[5] = 10.; maxVal[6] = 10.; maxVal[7] = 10.; maxVal[8]=10.; maxVal[9]=1.1;
  parName[0] = "Norm"; parName[1] = "#Lambda"; parName[2] = "G"; parName[3] = "Rch"; parName[4] = "Rcoh"; 
  parName[5] = "kappa_3"; parName[6] = "kappa_4"; parName[7] = "kappa_5"; parName[8]="kappa_6"; parName[9]="Norm_2";
      
  for (int i=0; i<npar; i++){
    MyMinuit.DefineParameter(i, parName[i].c_str(), par[i], stepSize[i], minVal[i], maxVal[i]);
  }
  
  //MyMinuit.DefineParameter(1, parName[1].c_str(), 0.7, stepSize[1], minVal[1], maxVal[1]); MyMinuit.FixParameter(1);// lambda
  //MyMinuit.DefineParameter(0, parName[0].c_str(), .998, stepSize[0], minVal[0], maxVal[0]); MyMinuit.FixParameter(0);// N
  MyMinuit.DefineParameter(2, parName[2].c_str(), 0., stepSize[2], minVal[2], maxVal[2]); MyMinuit.FixParameter(2);// G
  MyMinuit.DefineParameter(4, parName[4].c_str(), 0., stepSize[4], minVal[4], maxVal[4]); MyMinuit.FixParameter(4);// Rcoh
  MyMinuit.DefineParameter(7, parName[7].c_str(), 0, stepSize[7], minVal[7], maxVal[7]); MyMinuit.FixParameter(7);// k5
  MyMinuit.DefineParameter(8, parName[8].c_str(), 0, stepSize[8], minVal[8], maxVal[8]); MyMinuit.FixParameter(8);// k6
  
  //
  if(!IncludeEW){
    MyMinuit.DefineParameter(5, parName[5].c_str(), 0, stepSize[5], minVal[5], maxVal[5]); MyMinuit.FixParameter(5);// k3 
    MyMinuit.DefineParameter(6, parName[6].c_str(), 0, stepSize[6], minVal[6], maxVal[6]); MyMinuit.FixParameter(6);// k4
    MyMinuit.DefineParameter(7, parName[7].c_str(), 0, stepSize[7], minVal[7], maxVal[7]); MyMinuit.FixParameter(7);// k5
    MyMinuit.DefineParameter(8, parName[8].c_str(), 0, stepSize[8], minVal[8], maxVal[8]); MyMinuit.FixParameter(8);// k6
    
  }else{// IncludeEW
    if(FixEWavg){
      MyMinuit.DefineParameter(5, parName[5].c_str(), kappa3, stepSize[5], minVal[5], maxVal[5]); MyMinuit.FixParameter(5);// k3 
      MyMinuit.DefineParameter(6, parName[6].c_str(), kappa4, stepSize[6], minVal[6], maxVal[6]); MyMinuit.FixParameter(6);// k4
    }
  }

  if(IncludeSS==kFALSE){
    MyMinuit.DefineParameter(3, parName[3].c_str(), 9.1, stepSize[3], minVal[3], maxVal[3]); MyMinuit.FixParameter(3);// Rch
    MyMinuit.DefineParameter(0, parName[0].c_str(), .998, stepSize[0], minVal[0], maxVal[0]); MyMinuit.FixParameter(0);// N
    MyMinuit.DefineParameter(5, parName[5].c_str(), 0, stepSize[5], minVal[5], maxVal[5]); MyMinuit.FixParameter(5);// k3 
    MyMinuit.DefineParameter(6, parName[6].c_str(), 0, stepSize[6], minVal[6], maxVal[6]); MyMinuit.FixParameter(6);// k4
  }
  if(IncludeOS==kFALSE){
    MyMinuit.DefineParameter(9, parName[9].c_str(), 1.002, stepSize[9], minVal[9], maxVal[9]); MyMinuit.FixParameter(9);// N_2
  }
  
  
  int ierflg=0;
  double arglist[10];
  //arglist[0]=2;// improve Minimization Strategy (1 is default)
  //MyMinuit.mnexcm("SET STR",arglist,1,ierflg);
  //arglist[0] = 0;
  //MyMinuit.mnexcm("SCAN", arglist,1,ierflg);
  arglist[0] = 5000;
  MyMinuit.mnexcm("MIGRAD", arglist ,1,ierflg);
  // Do the minimization!
  if(!MCcase){
    cout<<"Start C2 Global fit"<<endl;
    MyMinuit.Migrad();// generally the best minimization algorithm
    cout<<"End C2 Global fit"<<endl;
  }
  for (int i=0; i<npar; i++){
    MyMinuit.GetParameter(i,OutputPar[i],OutputPar_e[i]);
  }
  
  cout<<"C2 fit: Chi2/NDF = "<<Chi2_C2global/(NFitPoints_C2global - MyMinuit.GetNumFreePars())<<"   Chi^2="<<Chi2_C2global<<endl;
  
  
  
  
 
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


  C2_ss->GetXaxis()->SetRangeUser(0,Q2Limit);// 0,0.098
  C2_ss->GetYaxis()->SetRangeUser(0.95,1.6);//0.98,1.3
  C2_ss->GetXaxis()->SetTitleOffset(.8);
  C2_ss->GetYaxis()->SetTitleOffset(.8);
  C2_ss->GetXaxis()->SetTitle("q_{inv} (GeV/c)");
  C2_ss->GetXaxis()->SetTitleSize(0.05);
  C2_ss->GetYaxis()->SetTitleSize(0.05);
  C2_os->GetXaxis()->SetRangeUser(0,0.6);
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
  
  
  TF1 *fitC2ss = new TF1("fitC2ss",C2ssFitFunction, 0.005,2, npar);//0.2
  TF1 *fitC2os = new TF1("fitC2os",C2osFitFunction, 0.005,2, npar);//Q2Limit
  for(int i=0; i<npar; i++) {
    fitC2ss->FixParameter(i,OutputPar[i]);
    fitC2ss->SetParError(i,OutputPar_e[i]);
  }
  
  TH1D *fitC2ss_h = new TH1D("fitC2ss_h","",C2_ss->GetNbinsX(),C2_ss->GetXaxis()->GetBinLowEdge(1), C2_ss->GetXaxis()->GetBinUpEdge(C2_ss->GetNbinsX()));
  for(int bin=1; bin<=C2_ss->GetNbinsX(); bin++){
    double qinv = C2_ss->GetXaxis()->GetBinCenter(bin);
    if(!MCcase) fitC2ss_h->SetBinContent(bin, fitC2ss->Eval(qinv));
  }
  fitC2ss_h->SetLineWidth(2);
  fitC2ss_h->SetLineColor(2);
  
  if(!MCcase){

    C2_ss->DrawCopy();
    legend1->AddEntry(C2_ss,"same-charge","p");
    C2_os->DrawCopy("same");
    legend1->AddEntry(C2_os,"mixed-charge","p");

    /*TF1 *BkgFitC2 = new TF1("BkgFitC2","1+[0]*exp(-pow([1]*x/0.19733,2))",0,1);
    BkgFitC2->SetParameter(0,0.08);
    BkgFitC2->SetParameter(1,0.5);
    BkgFitC2->SetLineColor(1);
    C2_ss->Fit(BkgFitC2,"IME","",0.2,0.8);
    BkgFitC2->Draw("same");*/
    //C2_os->DrawCopy("same");
    //C2_ss_Momsys->DrawCopy("E2 same");
    ///C2_os_Momsys->DrawCopy("E2 same");
    //C2_ss_Ksys->DrawCopy("E2 same");
    //C2_os_Ksys->DrawCopy("E2 same");
    //C2_ss->DrawCopy("same");
    //C2_os->DrawCopy("same");
    fitC2ss_h->Draw("C same");
    //fitC2os->SetLineColor(4);
    //fitC2os->DrawCopy("same");
    //MJ_bkg_ss->SetLineColor(1);
    //MJ_bkg_ss->Draw("same");
    //legend1->AddEntry(MJ_bkg_ss,"Non-femto bkg","l");
    legend1->Draw("same");
  }
  
  

  
 
  cout<<"============================================"<<endl;
  cout<<"Start 3-pion section"<<endl;
  
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
  //gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.02);
  
  TLegend *legend2 = new TLegend(.58,.55, .97,.65,NULL,"brNDC");
  legend2->SetBorderSize(1);
  legend2->SetTextSize(.03);// small .03; large .06
  legend2->SetFillColor(0);

  /////////////////////////////////////////////
  TH3D *C3QS_3d = new TH3D("C3QS_3d","",BINRANGE_3,0,.4, BINRANGE_3,0,.4, BINRANGE_3,0,.4);
  TH3D *Combinatorics_3d = new TH3D("Combinatorics_3d","",BINRANGE_3,0,.4, BINRANGE_3,0,.4, BINRANGE_3,0,.4);
  //TH3D *c3_num_3d = new TH3D("c3_num_3d","",BINRANGE_3,0,0.4, BINRANGE_3,0,0.4, BINRANGE_3,0,0.4);
  
  //
  const float Q3HistoLimit=0.5;
  const int Q3BINS = int((Q3HistoLimit+0.0001)/0.01);
  TH1D *num_withRS = new TH1D("num_withRS","",Q3BINS,0,Q3HistoLimit);
  TH1D *num_withGRS = new TH1D("num_withGRS","",Q3BINS,0,Q3HistoLimit);
  TH1D *num_withTM = new TH1D("num_withTM","",Q3BINS,0,Q3HistoLimit);
  TH1D *Cterm1 = new TH1D("Cterm1","",Q3BINS,0,Q3HistoLimit);
  TH1D *Cterm1_noMRC = new TH1D("Cterm1_noMRC","",Q3BINS,0,Q3HistoLimit);
  TH1D *Cterm2 = new TH1D("Cterm2","",Q3BINS,0,Q3HistoLimit);
  TH1D *Cterm3 = new TH1D("Cterm3","",Q3BINS,0,Q3HistoLimit);
  TH1D *Cterm4 = new TH1D("Cterm4","",Q3BINS,0,Q3HistoLimit);
  TH1D *num_QS = new TH1D("num_QS","",Q3BINS,0,Q3HistoLimit);
  TH1D *Combinatorics_1d = new TH1D("Combinatorics_1d","",Q3BINS,0,Q3HistoLimit);
  TH1D *Combinatorics_1d_noMRC = new TH1D("Combinatorics_1d_noMRC","",Q3BINS,0,Q3HistoLimit);
  TH1D *Coul_Riverside = new TH1D("Coul_Riverside","",Q3BINS,0,Q3HistoLimit);
  TH1D *Coul_GRiverside = new TH1D("Coul_GRiverside","",Q3BINS,0,Q3HistoLimit);
  TH1D *c3_hist = new TH1D("c3_hist","",Q3BINS,0,Q3HistoLimit);
  TH1D *dentimesFit_c3 = new TH1D("den_timesFit_c3","",Q3BINS,0,Q3HistoLimit);
  if(SameCharge) {Cterm1_noMRC->Sumw2(); Combinatorics_1d_noMRC->Sumw2();}
  //
  double num_QS_e[Q3BINS];
  double c3_e[Q3BINS];
  //double c3_3d_e[BINRANGE_3][BINRANGE_3][BINRANGE_3]={{{0}}};
  for(int ii=0; ii<Q3BINS; ii++){
    num_QS_e[ii]=0.;
    c3_e[ii]=0.;
    //for(int jj=0; jj<Q3BINS; jj++){
    //for(int kk=0; kk<Q3BINS; kk++){
    //	c3_3d_e[ii][jj][kk]=0;
    //}
    //}
    
  }
  
  // CB=Charge Bin. 0 means pi-, 1 means pi+
  int CB1=0, CB2=0, CB3=0;
  int CP12=1, CP13=1, CP23=1;
  int CB1_2=0, CB2_2=0, CB3_2=0;
  if(CHARGE==-1) {
    if(SameCharge) {CB1=0; CB2=0; CB3=0;}
    else {CB1=0; CB2=0; CB3=1; CP12=+1, CP13=-1, CP23=-1;}
  }else {
    if(SameCharge) {CB1=1; CB2=1; CB3=1;}
    else {CB1=0; CB2=1; CB3=1; CP12=-1, CP13=-1, CP23=+1;}
  }
  if(AddCC){
    if(CHARGE==-1) {
      if(SameCharge) {CB1_2=1; CB2_2=1; CB3_2=1;}
      else {CB1_2=0; CB2_2=1; CB3_2=1;};
    }else {
      if(SameCharge) {CB1_2=0; CB2_2=0; CB3_2=0;}
      else {CB1_2=0; CB2_2=0; CB3_2=1;}
    }
  }

  
  // SC = species combination.  Always 0 meaning pi-pi-pi.
  int SCBin=0;
  //
  ReadCoulCorrections(FSIindex);// switch to full kt range, 10.
  
  TH1D *GenSignalExpected_num=new TH1D("GenSignalExpected_num","",Q3BINS,0,Q3HistoLimit);
  TH1D *GenSignalExpected_den=new TH1D("GenSignalExpected_den","",Q3BINS,0,Q3HistoLimit);
  //
  double value_num; 
  double value_num_e;
  double N3_QS;
  double N3_QS_e;
  double OutTriplets=0, InTriplets=0;
  
  for(int i=2; i<=BINLIMIT_3; i++){// bin number
    //double Q12 = (i-0.5)*(0.01);// geometric center
    double Q12 = BinCenters[i-1];// true center
    //int Q12bin = int(Q12/0.01)+1;
    //
    for(int j=2; j<=BINLIMIT_3; j++){// bin number
      //double Q13 = (j-0.5)*(0.01);// geometric center
      double Q13 = BinCenters[j-1];// true center
      //int Q13bin = int(Q13/0.01)+1;
      //
      for(int k=2; k<=BINLIMIT_3; k++){// bin number
	//double Q23 = (k-0.5)*(0.01);// geometric center
	double Q23 = BinCenters[k-1];// true center
	//int Q23bin = int(Q23/0.01)+1;
	//
	//if(fabs(i-j)>3 || fabs(i-k)>3 || fabs(j-k)>3) continue;// testing
	A_3[i-1][j-1][k-1] = 0;
	B_3[i-1][j-1][k-1] = 0;
	A_3_e[i-1][j-1][k-1] = 0;
	double Q3 = sqrt(pow(Q12,2) + pow(Q13,2) + pow(Q23,2));
       	int Q3bin = Cterm1->GetXaxis()->FindBin(Q3);
	//if(Q3>=Q3Limit) continue;

	if(Q12 < sqrt(pow(Q13,2)+pow(Q23,2) - 2*Q13*Q23)) continue;// not all configurations are possible
	if(Q12 > sqrt(pow(Q13,2)+pow(Q23,2) + 2*Q13*Q23)) continue;// not all configurations are possible
	
	
       	//
	double G_12 = Gamov(CP12, Q12);// point-source Coulomb correlation
	double G_13 = Gamov(CP13, Q13);
	double G_23 = Gamov(CP23, Q23);
	double K2_12 = CoulCorr2(CP12, Q12);// finite-source Coulomb+Strong correlation from Therminator.
	double K2_13 = CoulCorr2(CP13, Q13);
	double K2_23 = CoulCorr2(CP23, Q23);
	double K3 = 1.;// 3-body Coulomb+Strong correlation
	if(SameCharge || CHARGE==-1) K3 = CoulCorrGRS(SameCharge, Q12, Q13, Q23);
	else K3 = CoulCorrGRS(SameCharge, Q23, Q12, Q13);
	
	if(MCcase && !GeneratedSignal) { K2_12=1.; K2_13=1.; K2_23=1.; K3=1.;}
	if(K3==0) continue;

	double TERM1=Three_3d[CB1][CB2][CB3][SCBin][0]->GetBinContent(i,j,k);// N3 (3-pion yield per q12,q13,q23 cell). 3-pions from same-event
	double TERM2=Three_3d[CB1][CB2][CB3][SCBin][1]->GetBinContent(i,j,k);// N2*N1. 1 and 2 from same-event
	double TERM3=Three_3d[CB1][CB2][CB3][SCBin][2]->GetBinContent(i,j,k);// N2*N1. 1 and 3 from same-event
	double TERM4=Three_3d[CB1][CB2][CB3][SCBin][3]->GetBinContent(i,j,k);// N2*N1. 2 and 3 from same-event
	double TERM5=Three_3d[CB1][CB2][CB3][SCBin][4]->GetBinContent(i,j,k);// N1*N1*N1. All from different events (pure combinatorics)
	if(AddCC){
	  TERM1 += Three_3d[CB1_2][CB2_2][CB3_2][SCBin][0]->GetBinContent(i,j,k);
	  TERM2 += Three_3d[CB1_2][CB2_2][CB3_2][SCBin][1]->GetBinContent(i,j,k);
	  TERM3 += Three_3d[CB1_2][CB2_2][CB3_2][SCBin][2]->GetBinContent(i,j,k);
	  TERM4 += Three_3d[CB1_2][CB2_2][CB3_2][SCBin][3]->GetBinContent(i,j,k);
	  TERM5 += Three_3d[CB1_2][CB2_2][CB3_2][SCBin][4]->GetBinContent(i,j,k);
	}
	
	if(TERM1==0 && TERM2==0 && TERM3==0 && TERM4==0 && TERM5==0) continue;
	if(TERM1==0 || TERM2==0 || TERM3==0 || TERM4==0 || TERM5==0) continue;
	

	if(Q3>0.08 && Q3<0.09){// just for testing
	  if(Q12>0.08 || Q13>0.08 || Q23>0.08) OutTriplets++;
	  else InTriplets++;
	}
	
	// apply momentum resolution correction
	if(!MCcase && !GeneratedSignal){
	  if(SameCharge){
	    TERM1 *= (MomRes1d[0][0]->GetBinContent(Q3bin)-1)*MRCShift+1;
	    TERM2 *= (MomRes1d[0][1]->GetBinContent(Q3bin)-1)*MRCShift+1;
	    TERM3 *= (MomRes1d[0][2]->GetBinContent(Q3bin)-1)*MRCShift+1;
	    TERM4 *= (MomRes1d[0][3]->GetBinContent(Q3bin)-1)*MRCShift+1;
	    TERM5 *= (MomRes1d[0][4]->GetBinContent(Q3bin)-1)*MRCShift+1;
	  }else{
	    TERM1 *= (MomRes1d[1][0]->GetBinContent(Q3bin)-1)*MRCShift+1;
	    TERM2 *= (MomRes1d[1][1]->GetBinContent(Q3bin)-1)*MRCShift+1;
	    TERM3 *= (MomRes1d[1][2]->GetBinContent(Q3bin)-1)*MRCShift+1;
	    TERM4 *= (MomRes1d[1][3]->GetBinContent(Q3bin)-1)*MRCShift+1;
	    TERM5 *= (MomRes1d[1][4]->GetBinContent(Q3bin)-1)*MRCShift+1;
	  }
	}
	//
	
	TwoFrac=0.7;
	OneFrac=pow(TwoFrac,.5); ThreeFrac=pow(TwoFrac,1.5);
	// finite-multiplicity method
	//OneFrac = (1+sqrt(1 + 4*AvgN[Mbin]*TwoFrac*(AvgN[Mbin]-1)))/(2*AvgN[Mbin]);
	//ThreeFrac = (OneFrac*AvgN[Mbin]*(OneFrac*AvgN[Mbin]-1)*(OneFrac*AvgN[Mbin]-2))/(AvgN[Mbin]*(AvgN[Mbin]-1)*(AvgN[Mbin]-2));
	
	
	// Purify.  Isolate pure 3-pion QS correlations using Lambda and K3 (removes lower order correlations)
	N3_QS = TERM1;
	N3_QS -= ( pow(1-OneFrac,3) + 3*OneFrac*pow(1-OneFrac,2) )*TERM5;
	N3_QS -= (1-OneFrac)*(TERM2 + TERM3 + TERM4 - 3*(1-TwoFrac)*TERM5);
	N3_QS /= ThreeFrac;
	N3_QS /= K3;
	
	num_QS->Fill(Q3, N3_QS);
	
	// Isolate 3-pion cumulant
	value_num = N3_QS;
	value_num -= (TERM2 - (1-TwoFrac)*TERM5)/K2_12/TwoFrac;
	value_num -= (TERM3 - (1-TwoFrac)*TERM5)/K2_13/TwoFrac;
	value_num -= (TERM4 - (1-TwoFrac)*TERM5)/K2_23/TwoFrac;
	value_num += 2*TERM5;
	
	
	
	
	// errors
	N3_QS_e = TERM1;
	N3_QS_e += pow(( pow(1-OneFrac,3) +3*OneFrac*pow(1-OneFrac,2) )*sqrt(TERM5),2);
	N3_QS_e += (pow((1-OneFrac),2)*(TERM2 + TERM3 + TERM4) + pow((1-OneFrac)*3*(1-TwoFrac)*sqrt(TERM5),2));
	N3_QS_e /= pow(ThreeFrac,2);
	N3_QS_e /= pow(K3,2);
	num_QS_e[Q3bin-1] += N3_QS_e;
	// errors 
	value_num_e = N3_QS_e;
	value_num_e += (pow(1/K2_12/TwoFrac*sqrt(TERM2),2) + pow((1-TwoFrac)/K2_12/TwoFrac*sqrt(TERM5),2));
	value_num_e += (pow(1/K2_13/TwoFrac*sqrt(TERM3),2) + pow((1-TwoFrac)/K2_13/TwoFrac*sqrt(TERM5),2));
	value_num_e += (pow(1/K2_23/TwoFrac*sqrt(TERM4),2) + pow((1-TwoFrac)/K2_23/TwoFrac*sqrt(TERM5),2));
	value_num_e += pow(2*sqrt(TERM5),2);
	c3_e[Q3bin-1] += value_num_e;// add baseline stat error
	//c3_3d_e[i-1][j-1][k-1] += value_num_e + TERM5;// add baseline stat error
		

	// Fill histograms
	c3_hist->Fill(Q3, value_num + TERM5);// for cumulant correlation function
	//c3_hist_3d->Fill(Q12,Q13,Q23, value_num[0] + TERM5);// for cumulant correlation function
	C3QS_3d->SetBinContent(i,j,k, N3_QS);
	C3QS_3d->SetBinError(i,j,k, N3_QS_e);
	//
	Coul_Riverside->Fill(Q3, G_12*G_13*G_23*TERM5);
	Coul_GRiverside->Fill(Q3, TERM5*CoulCorrGRS(SameCharge, Q12, Q13, Q23));
	num_withGRS->Fill(Q3, N3_QS);
	Cterm1->Fill(Q3, TERM1);
	Cterm2->Fill(Q3, TERM2);
	Cterm3->Fill(Q3, TERM3);
	Cterm4->Fill(Q3, TERM4);
	Combinatorics_1d->Fill(Q3, TERM5);
	Combinatorics_3d->Fill(Q12,Q13,Q23, TERM5);
	
	double Q3_m = sqrt(pow((i-0.5)*(0.005),2) + pow((j-0.5)*(0.005),2) + pow((k-0.5)*(0.005),2));
	Cterm1_noMRC->Fill(Q3_m, Three_3d[CB1][CB2][CB3][SCBin][0]->GetBinContent(i,j,k));
	Combinatorics_1d_noMRC->Fill(Q3_m, Three_3d[CB1][CB2][CB3][SCBin][4]->GetBinContent(i,j,k));
	
	//
	A_3[i-1][j-1][k-1] = value_num + TERM5;
	B_3[i-1][j-1][k-1] = TERM5;
	A_3_e[i-1][j-1][k-1] = sqrt(value_num_e);
	
	///////////////////////////////////////////////////////////
	// Edgeworth 3-pion expection.
	double radius_exp = (2.0)/FmToGeV;
	TwoFrac=0.7; OneFrac=pow(TwoFrac,.5), ThreeFrac=pow(TwoFrac,1.5);
	double arg12 = Q12*radius_exp;
	double arg13 = Q13*radius_exp;
	double arg23 = Q23*radius_exp;
	//Float_t EW12 = 1 + 0.24/(6.*pow(2.,1.5))*(8.*pow(arg12,3) - 12.*arg12);
	//EW12 += 0.16/(24.*pow(2.,2))*(16.*pow(arg12,4) -48.*pow(arg12,2) + 12);
	//Float_t EW13 = 1 + 0.24/(6.*pow(2.,1.5))*(8.*pow(arg13,3) - 12.*arg13);
	//EW13 += 0.16/(24.*pow(2.,2))*(16.*pow(arg13,4) -48.*pow(arg13,2) + 12);
	//Float_t EW23 = 1 + 0.24/(6.*pow(2.,1.5))*(8.*pow(arg23,3) - 12.*arg23);
	//EW23 += 0.16/(24.*pow(2.,2))*(16.*pow(arg23,4) -48.*pow(arg23,2) + 12);
	
	Float_t EW12 = 1 + 0.2/(6.*pow(2.,1.5))*(8.*pow(arg12,3) - 12.*arg12);
	EW12 += 0.45/(24.*pow(2.,2))*(16.*pow(arg12,4) -48.*pow(arg12,2) + 12);
	Float_t EW13 = 1 + 0.2/(6.*pow(2.,1.5))*(8.*pow(arg13,3) - 12.*arg13);
	EW13 += 0.45/(24.*pow(2.,2))*(16.*pow(arg13,4) -48.*pow(arg13,2) + 12);
	Float_t EW23 = 1 + 0.2/(6.*pow(2.,1.5))*(8.*pow(arg23,3) - 12.*arg23);
	EW23 += 0.45/(24.*pow(2.,2))*(16.*pow(arg23,4) -48.*pow(arg23,2) + 12);
	//
	//EW12=1; EW13=1; EW23=1;
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
	  double MJWeight12 = 0.1*exp(-pow(0.8*Q12/FmToGeV,2));
	  double MJWeight13 = 0.1*exp(-pow(0.4*Q13/FmToGeV,2));
	  double MJWeight23 = 0.1*exp(-pow(0.4*Q23/FmToGeV,2));
	  cumulant_exp = ThreeFrac*K3*(1 + exp(-pow(radius_exp*Q12,2))*pow(EW12,2) + (1+MJWeight12)*(1+MJWeight13)*(1+MJWeight23) - 1)*TERM5;
	  term1_exp = ( pow(1-OneFrac,3) + 3*OneFrac*pow(1-OneFrac,2) )*TERM5;
	  term1_exp += TwoFrac*(1-OneFrac)*(K2_12*(1+exp(-pow(radius_exp*Q12,2))*pow(EW12,2)+MJWeight12) + K2_13*(1+MJWeight13) + K2_23*(1+MJWeight23))*TERM5;
	  term1_exp += cumulant_exp;
	  //term1_exp = ((1-TwoFrac) + TwoFrac*K2_12*(1+exp(-pow(radius_exp*Q12,2))*pow(EW12,2)))*TERM5 + MJWeight12*TERM5;
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
  for(int i=0; i<Q3BINS; i++) {c3_hist->SetBinError(i+1, sqrt(c3_e[i])); num_QS->SetBinError(i+1, sqrt(num_QS_e[i]));}
  /*for(int i=0; i<BINRANGE_3; i++){
    for(int j=0; j<BINRANGE_3; j++){
      for(int k=0; k<BINRANGE_3; k++){
	c3_hist->SetBinError(i+1,j+1,k+1, sqrt(c3_3d_e[i][j][k]));
      }
    }
    }*/

  num_withRS->Divide(Combinatorics_1d);
  num_withGRS->Divide(Combinatorics_1d);
  num_withTM->Divide(Combinatorics_1d);
  Cterm1->Divide(Combinatorics_1d);
  Cterm1_noMRC->Divide(Combinatorics_1d_noMRC);
  Cterm2->Divide(Combinatorics_1d);
  Cterm3->Divide(Combinatorics_1d);
  Cterm4->Divide(Combinatorics_1d);
  c3_hist->Divide(Combinatorics_1d);
  num_QS->Divide(Combinatorics_1d);
  GenSignalExpected_num->Sumw2();
  GenSignalExpected_num->Divide(GenSignalExpected_den);
  
 
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  Cterm1->SetMarkerStyle(20);
  Cterm1->SetMarkerColor(4);
  Cterm1->SetLineColor(4);
  Cterm1->GetXaxis()->SetTitle("Q_{3} (GeV/c)");
  Cterm1->GetYaxis()->SetTitle("C_{3}");
  //Cterm1->SetTitle("#pi^{-}#pi^{-}#pi^{-}");
  Cterm1->SetMinimum(0.97);
  Cterm1->SetMaximum(5.7);// 6.1 for same-charge
  Cterm1->GetXaxis()->SetRangeUser(0,Q3Limit);
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
 
  c3_hist->GetXaxis()->SetTitle("Q_{3} (GeV/c)");
  c3_hist->GetYaxis()->SetTitle("C_{3} or c_{3}");
  c3_hist->GetYaxis()->SetTitleOffset(1.4);
  c3_hist->GetXaxis()->SetRangeUser(0,Q3Limit);
  c3_hist->GetYaxis()->SetRangeUser(0.9,4);
  c3_hist->SetMarkerStyle(25);
  c3_hist->SetMarkerColor(2);
  c3_hist->SetLineColor(2);
  c3_hist->SetMaximum(3);
  c3_hist->SetMinimum(.9);
  if(!MCcase) c3_hist->Draw("same");
  legend2->AddEntry(c3_hist,"#font[12]{c}_{3} (cumulant correlation)","p");
  c3_hist->Draw();
  
  
  const int npar_c3=5;// 10 
  TMinuit MyMinuit_c3(npar_c3);
  MyMinuit_c3.SetFCN(fcn_c3);
  double OutputPar_c3[npar_c3]={0};
  double OutputPar_c3_e[npar_c3]={0};
  
  double par_c3[npar_c3];               // the start values
  double stepSize_c3[npar_c3];          // step sizes 
  double minVal_c3[npar_c3];            // minimum bound on parameter 
  double maxVal_c3[npar_c3];            // maximum bound on parameter
  string parName_c3[npar_c3];
  par_c3[0] = 1.0; par_c3[1] = 0.5; par_c3[2] = 5; par_c3[3] = 0; par_c3[4] = 0;
  stepSize_c3[0] = 0.01; stepSize_c3[1] = 0.1; stepSize_c3[2] = 0.2; stepSize_c3[3] = 0.01; stepSize_c3[4] = 0.01;
  minVal_c3[0] = 0.95; minVal_c3[1] = 0.2; minVal_c3[2] = 0.5; minVal_c3[3] = -1; minVal_c3[4] = -1;
  maxVal_c3[0] = 1.1; maxVal_c3[1] = 2.0; maxVal_c3[2] = 15.; maxVal_c3[3] = +1; maxVal_c3[4] = +1;
  parName_c3[0] = "N"; parName_c3[1] = "#lambda_{3}"; parName_c3[2] = "R_{inv}"; parName_c3[3] = "#kappa_{3}"; parName_c3[4] = "#kappa_{4}";
  //TF1 *c3Fit1D;
 
  if(SameCharge && !MCcase){
    for (int i=0; i<npar_c3; i++){
      MyMinuit_c3.DefineParameter(i, parName_c3[i].c_str(), par_c3[i], stepSize_c3[i], minVal_c3[i], maxVal_c3[i]);
    }
    // kappa_3=0.2 and kappa_4=0.45 
    if(!IncludeEW){
      MyMinuit_c3.DefineParameter(3, parName_c3[3].c_str(), 0., stepSize_c3[3], minVal_c3[3], maxVal_c3[3]); MyMinuit_c3.FixParameter(3);
      MyMinuit_c3.DefineParameter(4, parName_c3[4].c_str(), 0., stepSize_c3[4], minVal_c3[4], maxVal_c3[4]); MyMinuit_c3.FixParameter(4);
    }else{
      if(FixEWavg){// <kappa3>=0.12, <kappa4>=0.43 for all 3 systems (error weighted)
	MyMinuit_c3.DefineParameter(3, parName_c3[3].c_str(), 0.12, stepSize_c3[3], minVal_c3[3], maxVal_c3[3]); MyMinuit_c3.FixParameter(3);
	MyMinuit_c3.DefineParameter(4, parName_c3[4].c_str(), 0.43, stepSize_c3[4], minVal_c3[4], maxVal_c3[4]); MyMinuit_c3.FixParameter(4);
      }
    }

    int ierflg_c3=0;
    double arglist_c3[10];
    //arglist_c3[0]=2;// improve Minimization Strategy (1 is default)
    //MyMinuit_c3.mnexcm("SET STR",arglist_c3,1,ierflg_c3);
    //arglist_c3[0] = 0;
    //MyMinuit_c3.mnexcm("SCAN", arglist_c3,1,ierflg_c3);
    arglist_c3[0] = 5000;
    MyMinuit_c3.mnexcm("MIGRAD", arglist_c3 ,1,ierflg_c3);
    // Do the minimization!
    cout<<"Start Three-d fit"<<endl;
    MyMinuit_c3.Migrad();// Minuit's best minimization algorithm
    cout<<"End Three-d fit"<<endl;
    cout<<"Chi2/NDF = "<<Chi2_c3/(NFitPoints_c3-MyMinuit_c3.GetNumFreePars())<<endl;
    
    for (int i=0; i<npar_c3; i++){
      MyMinuit_c3.GetParameter(i,OutputPar_c3[i],OutputPar_c3_e[i]);
    }
    
    TF3 *fitcopy_c3 = new TF3("fitcopy_c3",Dfitfunction_c3, 0,0.4, 0,0.4, 0,0.4, npar_c3);
    for(int i=0; i<npar_c3; i++) {fitcopy_c3->FixParameter(i,OutputPar_c3[i]);}
    for(int i=0; i<BINLIMIT_3; i++){
      for(int j=0; j<BINLIMIT_3; j++){
	for(int k=0; k<BINLIMIT_3; k++){
	  
	  double q12=BinCenters[i];
	  double q13=BinCenters[j];
	  double q23=BinCenters[k];
	  double q3=sqrt(pow(q12,2)+pow(q13,2)+pow(q23,2));
	  if(B_3[i][j][k]==0) {
	    dentimesFit_c3->Fill(q3, fitcopy_c3->Eval(q12,q13,q23));
	    Combinatorics_1d->Fill(q3, 1);
	    continue;
	  }
	 
	  dentimesFit_c3->Fill(q3, fitcopy_c3->Eval(q12,q13,q23)*B_3[i][j][k]);
	  dentimesFit_c3->SetBinError(dentimesFit_c3->GetXaxis()->FindBin(q3), 0);
	}
      }
    }

    //
    dentimesFit_c3->Divide(Combinatorics_1d);
    for(int i=0; i<Q3BINS; i++) dentimesFit_c3->SetBinError(i+1,0);
    
    dentimesFit_c3->SetLineColor(2);
    dentimesFit_c3->SetLineWidth(2);
    dentimesFit_c3->DrawCopy("C same");
    
   


    TF1 *c3Fit1D=new TF1("c3Fit1D","[0]*(1+[1]*exp(-pow([2]*x/0.19733,2)/2.))",0,1);
    //TF1 *c3Fit1D=new TF1("c3Fit1D","[0]*(1+[1]*exp(-pow([2]*x/0.19733,2)/2.) * (1 + (-0.12/(6.*pow(2,1.5))*(8.*pow([2]*x/0.19733,3) - 12.*pow([2]*x/0.19733,1))) + (0.43/(24.*pow(2,2))*(16.*pow([2]*x/0.19733,4) -48.*pow([2]*x/0.19733,2) + 12))))",0,1);
    //TF1 *c3Fit1D=new TF1("c3Fit1D","[0]*(1+[1]*exp(-pow([2]*x/0.19733,2)/2.) * (1 + ([3]/(6.*pow(2,1.5))*(8.*pow([2]*x/0.19733,3) - 12.*pow([2]*x/0.19733,1))) + ([4]/(24.*pow(2,2))*(16.*pow([2]*x/0.19733,4) -48.*pow([2]*x/0.19733,2) + 12))))",0,1);
    //c3Fit1D->SetLineColor(4);
    c3Fit1D->SetParameter(0,1); c3Fit1D->SetParName(0,"N");
    c3Fit1D->SetParameter(1,2); c3Fit1D->SetParName(1,"#lambda_{3}");
    c3Fit1D->SetParameter(2,5); c3Fit1D->SetParName(2,"R_{inv}");
    //c3Fit1D->FixParameter(1,.74);
    //c3Fit1D->SetParName(3,"#kappa_{3}"); c3Fit1D->SetParName(4,"#kappa_{4}");
    //c3Fit1D->SetParameter(3,.24); c3Fit1D->SetParameter(3,.16);
    //c3Fit1D->SetParLimits(2,0.5,10);
    //c3_hist->Fit(c3Fit1D,"IME","",0,Q3Limit);
    //cout<<"1d fit: Chi2/NDF = "<<c3Fit1D->GetChisquare()/c3Fit1D->GetNDF()<<endl;
    //dentimesFit_c3->DrawCopy("l same");
    //
    //c3Fit1D->FixParameter(0, fitcopy_c3->GetParameter(0));
    //c3Fit1D->FixParameter(1, fitcopy_c3->GetParameter(1));
    //c3Fit1D->FixParameter(2, fitcopy_c3->GetParameter(2));

    
  }
  //for(int i=1; i<50; i++) cout<<c3_hist->GetBinContent(i)<<", ";
  //cout<<endl;
  //for(int i=1; i<50; i++) cout<<c3_hist->GetBinError(i)<<", ";
  // pp M2, pi-
  double c3RefPoints[50]={0, 0, -4.18002, 1.04882, 3.09586, 2.47307, 2.97563, 2.23708, 2.59743, 2.03775, 1.85523, 1.78085, 1.72349, 1.72578, 1.55568, 1.51951, 1.57711, 1.48621, 1.39782, 1.40689, 1.36859, 1.25888, 1.2722, 1.24943, 1.21129, 1.19069, 1.16784, 1.15226, 1.13766, 1.13617, 1.1063, 1.09591, 1.09502, 1.08399, 1.07097, 1.06779, 1.06477, 1.06908, 1.03977, 1.04485, 1.02871, 1.03256, 1.03799, 1.02454, 1.01394, 1.01878, 1.00983, 1.01341, 1.01001};
  double c3RefPoints_e[50]={0, 0, 10.0715, 1.84526, 0.966399, 0.603807, 0.360672, 0.255103, 0.187488, 0.136952, 0.104286, 0.0887452, 0.0703067, 0.0580599, 0.0475478, 0.0442814, 0.0375324, 0.0333063, 0.0293538, 0.026505, 0.0236435, 0.0221139, 0.0197568, 0.0189292, 0.0168686, 0.0159787, 0.0148127, 0.0142807, 0.0133157, 0.0125173, 0.0117914, 0.0113195, 0.010635, 0.0103195, 0.0100161, 0.00945027, 0.00895993, 0.00877647, 0.00845139, 0.00819599, 0.00778514, 0.00766359, 0.00736487, 0.00704136, 0.00689525, 0.00683821, 0.00650538, 0.00633411, 0.00618844};
  //
  TH1D *c3Ref=(TH1D*)c3_hist->Clone();
  for(int i=1; i<50; i++) {c3Ref->SetBinContent(i,c3RefPoints[i-1]); c3Ref->SetBinError(i,c3RefPoints_e[i-1]);}
  c3Ref->SetMarkerStyle(20);
  //c3Ref->Draw("same");
  TLegend *legend3 = new TLegend(.58,.55, .97,.65,NULL,"brNDC");
  legend3->SetBorderSize(1);
  legend3->SetTextSize(.03);// small .03; large .06
  legend3->SetFillColor(0);
  legend3->AddEntry(c3_hist,"p-Pb","p");
  legend3->AddEntry(c3Ref,"p-p","p");
  //legend3->Draw("same");

  Cterm1->Draw("same");
  legend2->Draw("same");
  

  //for(int ii=0; ii<10; ii++) cout<<c3_hist->GetBinContent(ii+1)<<", ";
  //Coul_GRiverside->Draw();
  //Coul_Riverside->Draw("same");



  TLegend *legend4 = new TLegend(.3,.8, .9,.9,NULL,"brNDC");
  legend4->SetBorderSize(1);
  legend4->SetTextSize(.03);// small .03; large .06
  legend4->SetFillColor(0);
  TwoFrac=0.7; OneFrac=pow(TwoFrac,.5), ThreeFrac=pow(TwoFrac,1.5);
  
  TH1D *c3_Extended = (TH1D*)Three_1d[CB1][CB2][CB3][SCBin][0]->Clone();
  c3_Extended->Add(Three_1d[CB1][CB2][CB3][SCBin][4], -( pow(1-OneFrac,3) + 3*OneFrac*pow(1-OneFrac,2) ));
  c3_Extended->Add(Three_1d[CB1][CB2][CB3][SCBin][1], -(1-OneFrac));
  c3_Extended->Add(Three_1d[CB1][CB2][CB3][SCBin][2], -(1-OneFrac));
  c3_Extended->Add(Three_1d[CB1][CB2][CB3][SCBin][3], -(1-OneFrac));
  c3_Extended->Add(Three_1d[CB1][CB2][CB3][SCBin][4], (1-OneFrac)*3*(1-TwoFrac));
  c3_Extended->Scale(1/ThreeFrac);
  //
  c3_Extended->Add(Three_1d[CB1][CB2][CB3][SCBin][1], -1/TwoFrac);
  c3_Extended->Add(Three_1d[CB1][CB2][CB3][SCBin][2], -1/TwoFrac);
  c3_Extended->Add(Three_1d[CB1][CB2][CB3][SCBin][3], -1/TwoFrac);
  c3_Extended->Add(Three_1d[CB1][CB2][CB3][SCBin][4], 3*(1-TwoFrac)/TwoFrac);
  c3_Extended->Add(Three_1d[CB1][CB2][CB3][SCBin][4], 3);
  //
  c3_Extended->Divide(Three_1d[CB1][CB2][CB3][SCBin][4]);
  c3_Extended->GetXaxis()->SetTitle("Q_{3} (GeV/c)");
  c3_Extended->GetYaxis()->SetTitle("c_{3}");
  c3_Extended->SetMinimum(0.9); c3_Extended->SetMaximum(3.0);
  c3_Extended->SetMarkerStyle(24);
  c3_Extended->SetMarkerColor(2);
  //
  TH1D *C3_Extended = (TH1D*)Three_1d[CB1][CB2][CB3][SCBin][0]->Clone();
  C3_Extended->Divide(Three_1d[CB1][CB2][CB3][SCBin][4]);
  C3_Extended->GetXaxis()->SetTitle("Q_{3} (GeV/c)");
  C3_Extended->GetYaxis()->SetTitle("C_{3}");
  C3_Extended->SetMarkerStyle(20);
  C3_Extended->SetMarkerColor(4);
  //c3_Extended->Draw("same");
  //C3_Extended->Draw("same");
  //legend4->AddEntry(c3_Extended,"c_{3}, #pi^{#pm}#pi^{#pm}#pi^{#pm}, Coulomb Uncorrected","p");
  //legend4->AddEntry(C3_Extended,"C_{3}, #pi^{#pm}#pi^{#pm}#pi^{#pm}, Coulomb Uncorrected","p");
  //legend4->Draw("same");
  /*
  int binLow=10;
  int binHigh=50;
  //TwoFrac=0.5; OneFrac=pow(TwoFrac,0.5); ThreeFrac=pow(TwoFrac,1.5);
  TH1D *K0sProjection = (TH1D*)(Three_3d[0][0][1][SCBin][0]->ProjectionY("K0sProjection",binLow,binHigh,binLow,binHigh));
  TH1D *K0sProjection_term1 = (TH1D*)(Three_3d[0][0][1][SCBin][0]->ProjectionY("K0sProjection_term1",binLow,binHigh,binLow,binHigh));
  TH1D *K0sProjection_term2 = (TH1D*)(Three_3d[0][0][1][SCBin][1]->ProjectionY("K0sProjection_term2",binLow,binHigh,binLow,binHigh));
  TH1D *K0sProjection_term3 = (TH1D*)(Three_3d[0][0][1][SCBin][2]->ProjectionY("K0sProjection_term3",binLow,binHigh,binLow,binHigh));
  TH1D *K0sProjection_term4 = (TH1D*)(Three_3d[0][0][1][SCBin][3]->ProjectionY("K0sProjection_term4",binLow,binHigh,binLow,binHigh));
  TH1D *K0sProjection_term5 = (TH1D*)(Three_3d[0][0][1][SCBin][4]->ProjectionY("K0sProjection_term5",binLow,binHigh,binLow,binHigh));
  K0sProjection_term1->Add(K0sProjection_term5, -( pow(1-OneFrac,3) + 3*OneFrac*pow(1-OneFrac,2) ));
  K0sProjection_term1->Add(K0sProjection_term2, -(1-OneFrac));
  K0sProjection_term1->Add(K0sProjection_term3, -(1-OneFrac));
  K0sProjection_term1->Add(K0sProjection_term4, -(1-OneFrac));
  K0sProjection_term1->Add(K0sProjection_term5, (1-OneFrac)*3*(1-TwoFrac));
  K0sProjection_term1->Scale(1/ThreeFrac);
  //
  K0sProjection_term1->Add(K0sProjection_term2, -1/TwoFrac);
  K0sProjection_term1->Add(K0sProjection_term3, -1/TwoFrac);
  K0sProjection_term1->Add(K0sProjection_term4, -1/TwoFrac);
  K0sProjection_term1->Add(K0sProjection_term5, 3*(1-TwoFrac)/TwoFrac);
  K0sProjection_term1->Add(K0sProjection_term5, 3);
  //
  K0sProjection->Divide(K0sProjection_term5);
  K0sProjection_term1->Divide(K0sProjection_term5);
  K0sProjection_term1->GetXaxis()->SetRangeUser(0,.5);
  K0sProjection_term1->SetMarkerStyle(24);
  K0sProjection_term1->SetMarkerColor(4);
  K0sProjection->SetMarkerStyle(20);
  K0sProjection->SetMarkerColor(4);
  K0sProjection->SetMinimum(0.95);
  K0sProjection->GetXaxis()->SetTitle("q_{13}^{#pm#mp} (GeV/c)");
  K0sProjection->GetYaxis()->SetTitle("C_{3} or c_{3}"); 
  K0sProjection->GetYaxis()->SetTitleOffset(1.3);
  K0sProjection->Draw("");
  K0sProjection_term1->Draw("same");
  legend4->AddEntry(K0sProjection,"C_{3}^{#pm#pm#mp}, 0.1<q_{12}^{#pm#pm},q_{23}^{#pm#mp}<0.5 GeV/c","p");
  legend4->AddEntry(K0sProjection,"c_{3}^{#pm#pm#mp}, 0.1<q_{12,23}<0.5 GeV/c","p");
  legend4->Draw("same");
  */

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
  //MomRes_3d->Draw("same");
  //legend2->AddEntry(MomRes_3d,"3D","p");
   
  //legend2->Draw("same");
  
  
  /*
  ////////// C3 systematics
  // C3 --- base, M0, (0.035 TTC for all)
  //double C3_base[10]={0, 1.63072, 1.6096, 1.43731, 1.28205, 1.17777, 1.11973, 1.07932, 1.05459, 1.04029};
  //double C3_base_e[10]={0, 0.022528, 0.00387504, 0.00115667, 0.000423799, 0.000238973, 0.000143993, 9.71502e-05, 7.02007e-05, 5.31267e-05};
  // C3 --- base, M0, Pos B (0.035 TTC for all)
  //double C3_base[10]={0, 1.62564, 1.60461, 1.438, 1.28234, 1.17827, 1.12009, 1.07953, 1.05474, 1.04037};
  //double C3_base_e[10]={0, 0.0239233, 0.00409002, 0.0012215, 0.000446701, 0.000251769, 0.000151651, 0.000102284, 7.38993e-05, 5.59212e-05};
  // C3 --+ base, M0, (0.035 TTC for all)
  //double C3_base[10]={0, 1.66016, 1.41961, 1.25229, 1.16459, 1.11254, 1.08012, 1.05768, 1.04265, 1.0332};
  //double C3_base_e[10]={0, 0.00539779, 0.00111398, 0.000387926, 0.000192906, 0.00011428, 7.24765e-05, 5.09126e-05, 3.76421e-05, 2.87533e-05};
  
  // C3 --- base, M0, (New Standard: 0.04 TTC )
  //double C3_base[10]={0, 1.63903, 1.60254, 1.43381, 1.2794, 1.17603, 1.11806, 1.07806, 1.05345, 1.03936};
  //double C3_base_e[10]={0, 0.0322796, 0.00438361, 0.00122249, 0.000424557, 0.000233965, 0.000139058, 9.28136e-05, 6.66159e-05, 5.01816e-05};
  // C3 --- base, M1, (New Standard: 0.04 TTC )
  //double C3_base[10]={0, 1.6127, 1.65561, 1.49508, 1.33093, 1.21458, 1.14708, 1.099, 1.06864, 1.05064};
  //double C3_base_e[10]={0, 0.0425447, 0.0061176, 0.00172948, 0.000600699, 0.000329342, 0.000194832, 0.000129427, 9.25599e-05, 6.95395e-05};
  // C3 --- base, M2, (New Standard: 0.04 TTC )
  //double C3_base[10]={0, 1.6509, 1.71863, 1.54652, 1.38092, 1.25226, 1.17549, 1.12068, 1.08408, 1.06236};
  //double C3_base_e[10]={0, 0.0981473, 0.0143699, 0.00404612, 0.00141189, 0.000770764, 0.000453944, 0.000300452, 0.000214068, 0.000160448};
  // C3 --- base, M9, (New Standard: 0.04 TTC )
  //double C3_base[10]={0, 2.41982, 2.18303, 1.93205, 1.80399, 1.60955, 1.49305, 1.38565, 1.29873, 1.23626};
  //double C3_base_e[10]={0, 1.60227, 0.177274, 0.0501455, 0.018456, 0.00998147, 0.00583719, 0.00379465, 0.00264116, 0.0019391};
  //
  // C3 --+ base, M0, (New Standard: 0.04 TTC )
  //double C3_base[10]={0, 1.66087, 1.41943, 1.25081, 1.16313, 1.11143, 1.07917, 1.05701, 1.04215, 1.0328};
  //double C3_base_e[10]={0, 0.00584743, 0.00111278, 0.000374009, 0.00018315, 0.000107523, 6.78669e-05, 4.75116e-05, 3.50489e-05, 2.67328e-05};
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
  double C3_base[10]={0, 0.890897, 1.05222, 1.00461, 1.01364, 0.998981, 1.00225, 1.00305, 1.00235, 1.00043};
  double C3_base_e[10]={0, 0.297124, 0.0604397, 0.0195066, 0.00812906, 0.00490835, 0.00310751, 0.00217408, 0.00160575, 0.00123065};
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
  //C3_Sys->Fit(C3lineSys,"MEQ","",0,0.099);
  
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
  // c3 --- base, M0, (New Standard: 0.04 TTC )
  double c3_base[10]={0, 1.86014, 1.47533, 1.23733, 1.09944, 1.04145, 1.01693, 1.00715, 1.00253, 1.00111};
  double c3_base_e[10]={0, 0.104645, 0.0120917, 0.00333303, 0.00118126, 0.0006692, 0.000405246, 0.000274163, 0.000198507, 0.000150258};

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


  

  
  
  //////////////////////////////////////////////////////////////////////////////////////

  


 
  
  if(SaveToFile){
    TString *savefilename = new TString("OutFiles/OutFile");
    if(CollisionType==0) savefilename->Append("PbPb");
    else if(CollisionType==1) savefilename->Append("pPb");
    else savefilename->Append("pp");
    //
    if(MCcase) savefilename->Append("MonteCarlo");
    //
    if(SameCharge) savefilename->Append("SC");
    else savefilename->Append("MC");
    //
    if(CHARGE==1) savefilename->Append("Pos");
    else savefilename->Append("Neg");
    //
    if(IncludeEW) savefilename->Append("EW");
    
    savefilename->Append("Kt3_");
    *savefilename += EDbin+1;
        
    savefilename->Append("_M");
    *savefilename += Mbin;
    savefilename->Append(".root");
    TFile *savefile = new TFile(savefilename->Data(),"RECREATE");
    can1->Write("can1");
    C2_ss->Write("C2_ss");
    C2_os->Write("C2_os");
    fitC2ss_h->Write("fitC2ss_h");
    fitC2ss->Write("fitC2ss");
    Cterm1->Write("C3");
    c3_hist->Write("c3");
    //c3Fit1D->Write("c3Fit1D");
    Combinatorics_1d->Write("Combinatorics_1d");
    dentimesFit_c3->Write("dentimesFit_c3");
    MyMinuit_c3.Write("MyMinuit_c3");
    //ChiSquaredNDF->Write("ChiSquaredNDF");
    //
    savefile->Close();
  }
  
  
}
void ReadCoulCorrections(int index){
  ///////////////////////
  TString *fname2 = new TString("KFile.root");
 
  TFile *File=new TFile(fname2->Data(),"READ");
  if(index>=10) cout<<"FSI index not acceptable in 2-particle Coulomb read"<<endl;
  TString *nameSS = new TString("K2ss_");
  TString *nameOS = new TString("K2os_");
  *nameSS += index;
  *nameOS += index;
  TString *nameSS_2 = new TString("K2ss_");
  TString *nameOS_2 = new TString("K2os_");
  if(index<9) {*nameSS_2 += index+1; *nameOS_2 += index+1;}
  else {*nameSS_2 += index; *nameOS_2 += index;}
  TH1D *temp_ss = (TH1D*)File->Get(nameSS->Data());
  TH1D *temp_os = (TH1D*)File->Get(nameOS->Data());
  TH1D *temp_ss_2 = (TH1D*)File->Get(nameSS_2->Data());
  TH1D *temp_os_2 = (TH1D*)File->Get(nameOS_2->Data());

  CoulCorr2SS = (TH1D*)temp_ss->Clone();
  CoulCorr2OS = (TH1D*)temp_os->Clone();
  CoulCorr2SS_2 = (TH1D*)temp_ss_2->Clone();
  CoulCorr2OS_2 = (TH1D*)temp_os_2->Clone();
  // include strong FSI for pi+pi+ as per Lednicky code with R=7 (ratio of Coul+Strong to Coul) 
  for(int ii=1; ii<=CoulCorr2SS->GetNbinsX(); ii++){
    double newValue = CoulCorr2SS->GetBinContent(ii) * StrongSC->Eval(CoulCorr2SS->GetBinCenter(ii)*1000/2.);// k* not qinv
    CoulCorr2SS->SetBinContent(ii, newValue);
  }
  for(int ii=1; ii<=CoulCorr2SS_2->GetNbinsX(); ii++){
    double newValue = CoulCorr2SS_2->GetBinContent(ii) * StrongSC->Eval(CoulCorr2SS->GetBinCenter(ii)*1000/2.);// k* not qinv
    CoulCorr2SS_2->SetBinContent(ii, newValue);
  }
  
  //
  CoulCorr2SS->SetDirectory(0);
  CoulCorr2OS->SetDirectory(0);
  CoulCorr2SS_2->SetDirectory(0);
  CoulCorr2OS_2->SetDirectory(0);
  
  File->Close(); 
}
double CoulCorr2(int chargeproduct, double Q2){// returns K2
  int indexL=0;
  int indexH=0;
  double slope=0;
  double value1=1.0, value2=1.0;

  indexL = int(fabs(Q2 - CoulCorr2SS->GetXaxis()->GetBinWidth(1)/2.)/(CoulCorr2SS->GetXaxis()->GetBinWidth(1)));
  indexH = indexL+1;
 
  if(indexH >= CoulCorr2SS->GetNbinsX()) return 1.0;
  if(chargeproduct==1){
    slope = (CoulCorr2SS->GetBinContent(indexL+1) - CoulCorr2SS->GetBinContent(indexH+1));
    slope /= (CoulCorr2SS->GetXaxis()->GetBinCenter(indexL+1) - CoulCorr2SS->GetXaxis()->GetBinCenter(indexH+1));
    value1 = slope*(Q2 - CoulCorr2SS->GetXaxis()->GetBinCenter(indexL+1)) + CoulCorr2SS->GetBinContent(indexL+1);
  }else {
    slope = (CoulCorr2OS->GetBinContent(indexL+1) - CoulCorr2OS->GetBinContent(indexH+1));
    slope /= (CoulCorr2OS->GetXaxis()->GetBinCenter(indexL+1) - CoulCorr2OS->GetXaxis()->GetBinCenter(indexH+1));
    value1 = slope*(Q2 - CoulCorr2OS->GetXaxis()->GetBinCenter(indexL+1)) + CoulCorr2OS->GetBinContent(indexL+1);
  }
  if(Mbin_def<=9 || Mbin_def>12) return value1;
  else{// interpolate in K factor transition region
    if(chargeproduct==1){
      slope = (CoulCorr2SS_2->GetBinContent(indexL+1) - CoulCorr2SS_2->GetBinContent(indexH+1));
      slope /= (CoulCorr2SS_2->GetXaxis()->GetBinCenter(indexL+1) - CoulCorr2SS_2->GetXaxis()->GetBinCenter(indexH+1));
      value2 = slope*(Q2 - CoulCorr2SS_2->GetXaxis()->GetBinCenter(indexL+1)) + CoulCorr2SS_2->GetBinContent(indexL+1);
    }else {
      slope = (CoulCorr2OS_2->GetBinContent(indexL+1) - CoulCorr2OS_2->GetBinContent(indexH+1));
      slope /= (CoulCorr2OS_2->GetXaxis()->GetBinCenter(indexL+1) - CoulCorr2OS_2->GetXaxis()->GetBinCenter(indexH+1));
      value2 = slope*(Q2 - CoulCorr2OS_2->GetXaxis()->GetBinCenter(indexL+1)) + CoulCorr2OS_2->GetBinContent(indexL+1);
    }
    
    if(value1>0 && value2>0){
      return (value1 + (Mbin_def-9)*(value2-value1)/(3.));
    }else if(value1>0){
      return value1;
    }else if(value2>0){
      return value2;
    }else return 1.0;
    
  }
  //
  
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
  //double MT = sqrt(pow(0.25 + 0.1*Ktbin_GofP,2) + pow(masspiC,2));
  NFitPoints_C2global=0;
  if(LinkN) par[9]=par[0];// Link N factors

  for(int i=1; i<BINRANGE_C2global; i++){
    
    qinvSS = BinCenters[i];
    if(qinvSS > Q2Limit) continue;
    
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
    double Dp1 = fabs(par[2])/(1-fabs(par[2]));
    double Dp2 = fabs(par[2])/(1-fabs(par[2]));
    
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
    CSS += 1-par[1];
    CSS *= par[0];
    //
    COS = 1;
    if(ChargeConstraint && GofP) COS += 1/5.* Dp1/(1+Dp1) * Dp2/(1+Dp2);
    COS *= par[1]*K2OS[i];
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
  int qbin = int(fabs(x[0]/BinWidthQ2));
  //if(qbin >= BINRANGE_C2global) return 1.0;
  
  double qinvSS = BinCenters[qbin];
  //double MT = sqrt(pow(0.25 + 0.1*Ktbin_GofP,2) + pow(masspiC,2));
  
  if(!GofP) Dp = fabs(par[2])/(1-fabs(par[2]));// p independent
  else Dp = fabs(par[2])/(1-fabs(par[2]));
 
  //double Dp1 = fabs(par[2])/(1-fabs(par[2])) * exp(-2*(pow(Rcoh,2)-pow(RT,2))*pow(AvgP[Ktbin_GofP-1][qbin]-qinvSS/2.,2));
  //double Dp2 = fabs(par[2])/(1-fabs(par[2])) * exp(-2*(pow(Rcoh,2)-pow(RT,2))*pow(AvgP[Ktbin_GofP-1][qbin]+qinvSS/2.,2));
  //double Dp1 = fabs(par[2])/(1-fabs(par[2])) * exp(-2*pow(Rcoh,2)*pow(AvgP[Ktbin_GofP-1][qbin]-qinvSS/2.,2) - 2*MT/Temp);
  //double Dp2 = fabs(par[2])/(1-fabs(par[2])) * exp(-2*pow(Rcoh,2)*pow(AvgP[Ktbin_GofP-1][qbin]+qinvSS/2.,2) - 2*MT/Temp);
  //double Dp1 = fabs(par[2])/(1-fabs(par[2])) * exp(-2*pow(Rcoh,2)*pow(qinvSS/2.,2) - 2*MT/Temp);
  //double Dp2 = fabs(par[2])/(1-fabs(par[2])) * exp(-2*pow(Rcoh,2)*pow(qinvSS/2.,2) - 2*MT/Temp);
  double Dp1 = fabs(par[2])/(1-fabs(par[2]));
  double Dp2 = fabs(par[2])/(1-fabs(par[2]));
  
  if(!GofP) {Dp1=Dp; Dp2=Dp;}
  double arg=qinvSS*Rch;
  double EW = 1 + par[5]/(6.*pow(2,1.5))*(8.*pow(arg,3) - 12.*pow(arg,1));
  EW += par[6]/(24.*pow(2,2))*(16.*pow(arg,4) -48.*pow(arg,2) + 12);
  EW += par[7]/(120.*pow(2,2.5))*(32.*pow(arg,5) - 160.*pow(arg,3) + 120*arg);
  EW += par[8]/(720.*pow(2,3))*(64.*pow(arg,6) - 480.*pow(arg,4) + 720.*pow(arg,2) - 120);
  double Gaus_coh = exp(-pow(Rcoh*qinvSS,2)/2.);
  double Gaus_ch = exp(-pow(Rch*qinvSS,2)/2.) * EW + 0.00001;// Add on a tiny amount to avoid float exception at high q
  double CSS = 1 + pow(1 + Dp*Gaus_coh/Gaus_ch,2)/((1+Dp1)*(1+Dp2)) * exp(-pow(Rch*qinvSS,2))*pow(EW,2);
  double K2=1.0;
  if(qbin < BINRANGE_C2global) K2=K2SS[qbin];
  if(ChargeConstraint) CSS -= 4/5.* Dp1/(1+Dp1) * Dp2/(1+Dp2);
  else CSS -= pow(Gaus_coh*Dp,2)/((1+Dp1)*(1+Dp2));
  //else CSS -= Dp1/(1+Dp1) * Dp2/(1+Dp2);
  CSS *= par[1]*K2;
  CSS += 1-par[1];
  CSS *= par[0];
  
  //if(qinvSS > .85 && qinvSS < .9) cout<<"In function: "<<qinvSS<<"  "<<CSS<<endl;
  
  return CSS;
}
//________________________________________________________________________
double C2osFitFunction(Double_t *x, Double_t *par)
{
  if(LinkN) par[9]=par[0];// Link N factors
  double Rcoh=par[4]/FmToGeV;
  double Dp=0;
  int qbin = int(fabs(x[0]/BinWidthQ2));
  //if(qbin >= BINRANGE_C2global) return 1.0;
  
  double qinvOS = BinCenters[qbin];
  
  if(!GofP) Dp = fabs(par[2])/(1-fabs(par[2]));// p independent
  else Dp = fabs(par[2])/(1-fabs(par[2]));
  //Dp = fabs(par[2])/(1-fabs(par[2]));
  double Dp1 = fabs(par[2])/(1-fabs(par[2]));
  double Dp2 = fabs(par[2])/(1-fabs(par[2]));
  //double Dp1 = fabs(par[2])/(1-fabs(par[2])) * exp(-2*(pow(Rcoh,2)-pow(RT,2))*pow(qinvOS/2.,2));
  //double Dp2 = fabs(par[2])/(1-fabs(par[2])) * exp(-2*(pow(Rcoh,2)-pow(RT,2))*pow(qinvOS/2.,2));

  if(!GofP) {Dp1=Dp; Dp2=Dp;}
  double COS = 1;
  double K2=1.0;
  if(qbin < BINRANGE_C2global) K2=K2OS[qbin];
  if(ChargeConstraint && GofP) COS += 1/5.* Dp1/(1+Dp1) * Dp2/(1+Dp2);
  COS *= par[1]*K2;
  COS += 1-par[1];
  COS *= par[9];
  return COS;
}
//__________________________________________________________________________

void fcn_c3(int& npar, double* deriv, double& f, double par[], int flag){

  double q12=0, q13=0, q23=0;
  double EW12=0, EW13=0, EW23=0;
  double C=0;
  double Rch=par[2]/FmToGeV;
  double SumChi2=0;
  //double lnL=0, term1=0, term2=0;
  NFitPoints_c3=0;

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
	if(q3 > Q3Limit) continue;
	//
	double arg12 = q12*Rch;
	double arg13 = q13*Rch;
	double arg23 = q23*Rch;
	EW12 = 1 + par[3]/(6.*pow(2,1.5))*(8*pow(arg12,3) - 12*pow(arg12,1)) + par[4]/(24.*pow(2,2))*(16*pow(arg12,4) -48*pow(arg12,2) + 12);
	EW13 = 1 + par[3]/(6.*pow(2,1.5))*(8*pow(arg13,3) - 12*pow(arg13,1)) + par[4]/(24.*pow(2,2))*(16*pow(arg13,4) -48*pow(arg13,2) + 12);
	EW23 = 1 + par[3]/(6.*pow(2,1.5))*(8*pow(arg23,3) - 12*pow(arg23,1)) + par[4]/(24.*pow(2,2))*(16*pow(arg23,4) -48*pow(arg23,2) + 12);
	//EW12=1; EW13=1; EW23=1;
	//
	C = 1 + par[1]*exp(-(pow(arg12,2)+pow(arg13,2)+pow(arg23,2))/2.)*EW12*EW13*EW23;
	C *= par[0];// Norm
	
	double error = pow(A_3_e[i][j][k]/B_3[i][j][k],2);
	//error += pow( fabs(MixedChargeSysFit->Eval(q3)-1) * A_3[i][j][k]/B_3[i][j][k],2);
	//error += pow( 0.1 * (MomResC2[0]->GetBinContent(MomResC2[0]->GetXaxis()->FindBin(q3))-1) * A_3[i][j][k]/B_3[i][j][k],2);
	error = sqrt(error);
	SumChi2 += pow( (C - A_3[i][j][k]/B_3[i][j][k])/error,2);
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
      
}
//________________________________________________________________________
double Dfitfunction_c3(Double_t *x, Double_t *par)
{
  double Rch = par[2]/FmToGeV;
  double arg12 = x[0]*Rch;
  double arg13 = x[1]*Rch;
  double arg23 = x[2]*Rch;
  double EW12 = 1 + par[3]/(6.*pow(2,1.5))*(8*pow(arg12,3) - 12*pow(arg12,1)) + par[4]/(24.*pow(2,2))*(16*pow(arg12,4) -48*pow(arg12,2) + 12);
  double EW13 = 1 + par[3]/(6.*pow(2,1.5))*(8*pow(arg13,3) - 12*pow(arg13,1)) + par[4]/(24.*pow(2,2))*(16*pow(arg13,4) -48*pow(arg13,2) + 12);
  double EW23 = 1 + par[3]/(6.*pow(2,1.5))*(8*pow(arg23,3) - 12*pow(arg23,1)) + par[4]/(24.*pow(2,2))*(16*pow(arg23,4) -48*pow(arg23,2) + 12);
  //
  double C = 1 + par[1]*exp(-(pow(arg12,2)+pow(arg13,2)+pow(arg23,2))/2.)*EW12*EW13*EW23;
  C *= par[0];// Norm
  
  return C;
}
//________________________________________________________________________
double CoulCorrGRS(bool SC, double Q_12, double Q_13, double Q_23){
 
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
    double value1=1.0, value2=1.0;
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
      value1 = (c0*(1-zd) + c1*zd);
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
      value1 = (c0*(1-zd) + c1*zd);
    }
    
    if(Mbin_def<=9 || Mbin_def>12) return value1;
    //
    // interpolation for K factor tansition
    if(SC){
      double c00 = CoulCorr2SS_2->GetBinContent(index12L+1)*CoulCorr2SS_2->GetBinContent(index13L+1)*CoulCorr2SS_2->GetBinContent(index23L+1)*(1-xd);
      c00 += CoulCorr2SS_2->GetBinContent(index12H+1)*CoulCorr2SS_2->GetBinContent(index13L+1)*CoulCorr2SS_2->GetBinContent(index23L+1)*xd;
      double c10 = CoulCorr2SS_2->GetBinContent(index12L+1)*CoulCorr2SS_2->GetBinContent(index13H+1)*CoulCorr2SS_2->GetBinContent(index23L+1)*(1-xd);
      c10 += CoulCorr2SS_2->GetBinContent(index12H+1)*CoulCorr2SS_2->GetBinContent(index13H+1)*CoulCorr2SS_2->GetBinContent(index23L+1)*xd;
      double c01 = CoulCorr2SS_2->GetBinContent(index12L+1)*CoulCorr2SS_2->GetBinContent(index13L+1)*CoulCorr2SS_2->GetBinContent(index23H+1)*(1-xd);
      c01 += CoulCorr2SS_2->GetBinContent(index12H+1)*CoulCorr2SS_2->GetBinContent(index13L+1)*CoulCorr2SS_2->GetBinContent(index23H+1)*xd;
      double c11 = CoulCorr2SS_2->GetBinContent(index12L+1)*CoulCorr2SS_2->GetBinContent(index13H+1)*CoulCorr2SS_2->GetBinContent(index23H+1)*(1-xd);
      c11 += CoulCorr2SS_2->GetBinContent(index12H+1)*CoulCorr2SS_2->GetBinContent(index13H+1)*CoulCorr2SS_2->GetBinContent(index23H+1)*xd;
      //
      double c0 = c00*(1-yd) + c10*yd;
      double c1 = c01*(1-yd) + c11*yd;
      value2 = (c0*(1-zd) + c1*zd);
    }else{
      double c00 = CoulCorr2SS_2->GetBinContent(index12L+1)*CoulCorr2OS_2->GetBinContent(index13L+1)*CoulCorr2OS_2->GetBinContent(index23L+1)*(1-xd);
      c00 += CoulCorr2SS_2->GetBinContent(index12H+1)*CoulCorr2OS_2->GetBinContent(index13L+1)*CoulCorr2OS_2->GetBinContent(index23L+1)*xd;
      double c10 = CoulCorr2SS_2->GetBinContent(index12L+1)*CoulCorr2OS_2->GetBinContent(index13H+1)*CoulCorr2OS_2->GetBinContent(index23L+1)*(1-xd);
      c10 += CoulCorr2SS_2->GetBinContent(index12H+1)*CoulCorr2OS_2->GetBinContent(index13H+1)*CoulCorr2OS_2->GetBinContent(index23L+1)*xd;
      double c01 = CoulCorr2SS_2->GetBinContent(index12L+1)*CoulCorr2OS_2->GetBinContent(index13L+1)*CoulCorr2OS_2->GetBinContent(index23H+1)*(1-xd);
      c01 += CoulCorr2SS_2->GetBinContent(index12H+1)*CoulCorr2OS_2->GetBinContent(index13L+1)*CoulCorr2OS_2->GetBinContent(index23H+1)*xd;
      double c11 = CoulCorr2SS_2->GetBinContent(index12L+1)*CoulCorr2OS_2->GetBinContent(index13H+1)*CoulCorr2OS_2->GetBinContent(index23H+1)*(1-xd);
      c11 += CoulCorr2SS_2->GetBinContent(index12H+1)*CoulCorr2OS_2->GetBinContent(index13H+1)*CoulCorr2OS_2->GetBinContent(index23H+1)*xd;
      //
      double c0 = c00*(1-yd) + c10*yd;
      double c1 = c01*(1-yd) + c11*yd;
      value2 = (c0*(1-zd) + c1*zd);
    }
    
    if(value1>0 && value2>0){
      return (value1 + (Mbin_def-9)*(value2-value1)/(3.));
    }else if(value1>0){
      return value1;
    }else if(value2>0){
      return value2;
    }else return 1.0;
    
  }
  
}

double Gamov(int chargeProduct, double qinv){
  
  double arg = chargeProduct*2.*PI/(BohrR*qinv/2.);
  
  return arg/(exp(arg)-1);
}
void ReadMomResFile(int Mbin){
 
  int MRCindex_L=0, MRCindex_H=0;
  float MRCindex_weight=0;
  if(Mbin<=2) {MRCindex_L=0; MRCindex_H=1; MRCindex_weight = Mbin/2.;}
  else if(Mbin<=6) {MRCindex_L=1; MRCindex_H=2; MRCindex_weight = fabs(Mbin-3)/3.;}
  else if(Mbin<=11) {MRCindex_L=2; MRCindex_H=3; MRCindex_weight = fabs(Mbin-7)/4.;}
  else {MRCindex_L=3; MRCindex_H=4; MRCindex_weight = fabs(Mbin-12)/7.;}
 
  
  TH1D *temp_L[2][5];
  TH1D *temp_H[2][5];
  TH1D *tempC2_L[2];
  TH1D *tempC2_H[2];
  TString *momresfilenameL = new TString("MomResFile_M");
  *momresfilenameL += MRCindex_L;
  momresfilenameL->Append(".root");
  TFile *MomResFileL = new TFile(momresfilenameL->Data(),"READ");
  TString *names1D_L[2][5];// SC/MC, term#
  for(int ChProd=0; ChProd<2; ChProd++){
    for(int term=0; term<5; term++){
      //
      if(ChProd==0) {names1D_L[ChProd][term] = new TString("MRC3_SC_term");}
      else {names1D_L[ChProd][term] = new TString("MRC3_MC_term");}
      *names1D_L[ChProd][term] += term+1;
      temp_L[ChProd][term] = (TH1D*)MomResFileL->Get(names1D_L[ChProd][term]->Data());
      temp_L[ChProd][term]->SetDirectory(0);
      
    }
    TString *C2MRCname=new TString("MomResHisto_");
    if(ChProd==0) C2MRCname->Append("pp");
    else C2MRCname->Append("mp");
    tempC2_L[ChProd]=(TH1D*)(((TH2D*)(MomResFileL->Get(C2MRCname->Data())))->ProjectionY("C2MRCproL",MRC2index,MRC2index));
    tempC2_L[ChProd]->SetDirectory(0);
  }
  
  //
  TString *momresfilenameH = new TString("MomResFile_M");
  *momresfilenameH += MRCindex_H;
  momresfilenameH->Append(".root");
  TFile *MomResFileH = new TFile(momresfilenameH->Data(),"READ");
  TString *names1D_H[2][5];// SC/MC, term#
  for(int ChProd=0; ChProd<2; ChProd++){
    for(int term=0; term<5; term++){
      //
      if(ChProd==0) {names1D_H[ChProd][term] = new TString("MRC3_SC_term");}
      else {names1D_H[ChProd][term] = new TString("MRC3_MC_term");}
      *names1D_H[ChProd][term] += term+1;
      temp_H[ChProd][term] = (TH1D*)MomResFileH->Get(names1D_H[ChProd][term]->Data());
      temp_H[ChProd][term]->SetDirectory(0);
      MomRes1d[ChProd][term] = (TH1D*)MomResFileH->Get(names1D_H[ChProd][term]->Data());
      MomRes1d[ChProd][term]->SetDirectory(0);
    }
    TString *C2MRCname=new TString("MomResHisto_");
    if(ChProd==0) C2MRCname->Append("pp");
    else C2MRCname->Append("mp");
    tempC2_H[ChProd]=(TH1D*)(((TH2D*)(MomResFileH->Get(C2MRCname->Data())))->ProjectionY("C2MRCproH",MRC2index,MRC2index));
    tempC2_H[ChProd]->SetDirectory(0);
    TString *name=new TString("MomResC2_");
    *name += ChProd;
    MomResC2[ChProd] = (TH1D*)(((TH2D*)(MomResFileH->Get(C2MRCname->Data())))->ProjectionY(name->Data(),MRC2index,MRC2index));
    MomResC2[ChProd]->SetDirectory(0);
  }
  
  for(int ChProd=0; ChProd<2; ChProd++){
    // C3 MRC
    for(int term=0; term<5; term++){
      for(int bin=1; bin<=temp_H[ChProd][term]->GetNbinsX(); bin++){
	double value=1;
	if(temp_L[ChProd][term]->GetBinContent(bin)>0 && temp_H[ChProd][term]->GetBinContent(bin)>0){// both have entries
	  value = temp_L[ChProd][term]->GetBinContent(bin) + MRCindex_weight * temp_H[ChProd][term]->GetBinContent(bin);
	}else if(temp_L[ChProd][term]->GetBinContent(bin)>0){
	  value = temp_L[ChProd][term]->GetBinContent(bin);
	}else if(temp_H[ChProd][term]->GetBinContent(bin)>0){
	  value = temp_H[ChProd][term]->GetBinContent(bin);
	}else value=1.0;
	
	MomRes1d[ChProd][term]->SetBinContent(bin,value);
      }
    }
    // C2 MRC
    for(int bin=1; bin<=tempC2_H[ChProd]->GetNbinsX(); bin++){
      double value=1;
      if(tempC2_L[ChProd]->GetBinContent(bin)>0 && tempC2_H[ChProd]->GetBinContent(bin)>0){// both have entries
	value = tempC2_L[ChProd]->GetBinContent(bin)*(1-MRCindex_weight) + MRCindex_weight * tempC2_H[ChProd]->GetBinContent(bin);
      }else if(tempC2_L[ChProd]->GetBinContent(bin)>0){
	value = tempC2_L[ChProd]->GetBinContent(bin);
      }else if(tempC2_H[ChProd]->GetBinContent(bin)>0){
	value = tempC2_H[ChProd]->GetBinContent(bin);
      }else value=1.0;
      MomResC2[ChProd]->SetBinContent(bin,value);
    } 
  }

  for(int ChProd=0; ChProd<2; ChProd++){
    for(int term=0; term<5; term++){
      for(int i=0; i<MomRes1d[0][0]->GetNbinsX(); i++){
	if(SourceType==0 && Mbin<10){
	  if(MomRes1d[ChProd][term]->GetXaxis()->GetBinCenter(i+1) > 0.1) MomRes1d[ChProd][term]->SetBinContent(i+1, 1.0);
	}else if(SourceType==0 && Mbin>=10){
	  if(MomRes1d[ChProd][term]->GetXaxis()->GetBinCenter(i+1) > 0.2) MomRes1d[ChProd][term]->SetBinContent(i+1, 1.0);
	}else{
	  if(MomRes1d[ChProd][term]->GetXaxis()->GetBinCenter(i+1) > 0.5) MomRes1d[ChProd][term]->SetBinContent(i+1, 1.0);
	}

      }
    }
    
  }
  
  MomResFileL->Close();
  MomResFileH->Close();

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
