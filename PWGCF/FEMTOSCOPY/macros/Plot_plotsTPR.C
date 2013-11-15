#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <Riostream.h>

#include "TVector2.h"
#include "TFile.h"
#include "TString.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TPad.h"
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
#include "TASImage.h"
#include "TGraphErrors.h"

#define BohrR 1963.6885
#define FmToGeV 0.19733 // conversion to fm
#define PI 3.1415926
#define masspiC 0.1395702 // pi+ mass (GeV/c^2)

using namespace std;

bool SaveFiles=kFALSE;
const int ChProdBOI=0;// 0=SameCharge, 1=MixedCharge
const int KT3Bin=0;// Kt3 bin. 0-2
bool IncludeEW=kFALSE;
bool AddedCC=kTRUE;// Charge Conjugate already added?
//
int MbinMaxPbPb=15;// 15
int MbinMinpPb=13;// 13
int MbinMinpp=14;// 14 
int MbinMinPbPb=0;// 0
int MbinMaxpPb=18;// 18
//
//
//
const int MaxKT3Bins=2;
int TextFont=42;// 63, or 42
float SizeLabel=0.08;// 20(63 font), 0.08(42 font)
float SizeLegend=0.07;// 
float SizeTitle=0.08;// 
float SizeSpecif=0.06;// 
float SF1=2/3.*0.95;
float SF2=1/2.*0.95;

double RightMargin=0.002;// 0.002
//
double Chi2_C2global;
double NFitPoints_C2global;
TH1D *C2_ss[20];
TH1D *C2_os[20];



void DrawALICELogo(Bool_t, Float_t, Float_t, Float_t, Float_t);

void Plot_plotsTPR(){

  gStyle->SetOptStat(0);
  gStyle->SetOptDate(0);
  
  ////////////////////////////////////
  // Get Nrec to Nch mapping
  double meanNchPbPb[20]={0};
  double meanNchPbPb_e[20]={0};
  double meanNchpPb[20]={0};
  double meanNchpPb_e[20]={0};
  double meanNchpp[20]={0};
  double meanNchpp_e[20]={0};
  TFile *NrecMapFile;
  TH2D *NrecMap;
  TList *MyList;
  //
  NrecMapFile = new TFile("Results/NrecMapping_12a17a.root","READ");
  //NrecMapFile = new TFile("Results/NrecMapping_12a11a.root","READ");// MC variation
  //NrecMapFile = new TFile("Results/NrecMapping_12a17a_FB5and7overlap.root","READ");
  MyList=(TList*)NrecMapFile->Get("MyList");
  NrecMap = (TH2D*)MyList->FindObject("fNchTrueDist");
  for(int bin=1; bin<=2; bin++){// 1 to 2 (FB7),  1 to 1 (AMPT),  1 to 4 (FB5and7overlap)
    NrecMap->GetXaxis()->SetRangeUser(bin,bin);
    if(NrecMap->GetMean(2)>0) meanNchPbPb[bin-1] = log10(NrecMap->GetMean(2));
  }
  NrecMapFile->Close();
  //
  NrecMapFile = new TFile("Results/NrecMapping_12a17e.root","READ");
  //NrecMapFile = new TFile("Results/NrecMapping_12a11b.root","READ");// MC variation
  //NrecMapFile = new TFile("Results/NrecMapping_12a17e_FB5and7overlap.root","READ");
  MyList=(TList*)NrecMapFile->Get("MyList");
  NrecMap = (TH2D*)MyList->FindObject("fNchTrueDist");
  for(int bin=3; bin<=10; bin++){// 3 to 10 (FB7),  2 to 3 (AMPT), 5 to 12 (FB5and7overlap)
    NrecMap->GetXaxis()->SetRangeUser(bin,bin);
    if(NrecMap->GetMean(2)>0) meanNchPbPb[bin-1] = log10(NrecMap->GetMean(2));
  }
  NrecMapFile->Close();
  //
  ////////////////////// extra for AMPT 
  /*NrecMapFile = new TFile("Results/NrecMapping_12a11d.root","READ");
  MyList=(TList*)NrecMapFile->Get("MyList");
  NrecMap = (TH2D*)MyList->FindObject("fNchTrueDist");
  for(int bin=4; bin<=7; bin++){// 4 to 7 (AMPT)
    NrecMap->GetXaxis()->SetRangeUser(bin,bin);
    if(NrecMap->GetMean(2)>0) meanNchPbPb[bin-1] = log10(NrecMap->GetMean(2));
  }
  NrecMapFile->Close();*/
  //////////////////////////////////////////
  //
  NrecMapFile = new TFile("Results/NrecMapping_12a17c.root","READ");
  //NrecMapFile = new TFile("Results/NrecMapping_12a11g.root","READ");// MC variation
  //NrecMapFile = new TFile("Results/NrecMapping_12a17c_FB5and7overlap.root","READ");
  MyList=(TList*)NrecMapFile->Get("MyList");
  NrecMap = (TH2D*)MyList->FindObject("fNchTrueDist");
  for(int bin=11; bin<=19; bin++){// 11 to 19 (FB7),  1 to 1 (AMPT), 13 to 19 (FB5and7overlap)
    NrecMap->GetXaxis()->SetRangeUser(bin,bin);
    if(NrecMap->GetMean(2)>0) meanNchPbPb[bin-1] = log10(NrecMap->GetMean(2));
  }
  NrecMapFile->Close();
  //
  NrecMapFile = new TFile("Results/PDC_13b2_efix_p1_R2.root","READ");
  //NrecMapFile = new TFile("Results/NrecMapping_13b3.root","READ");// MC variation
  //NrecMapFile = new TFile("Results/NrecMapping_13b2_efix_p1_FB5and7overlap.root","READ");
  MyList=(TList*)NrecMapFile->Get("MyList");
  NrecMap = (TH2D*)MyList->FindObject("fNchTrueDist");
  for(int bin=10; bin<=20; bin++){
    NrecMap->GetXaxis()->SetRangeUser(bin,bin);
    if(NrecMap->GetMean(2)>0) meanNchpPb[bin-1] = log10(NrecMap->GetMean(2));
  }
  NrecMapFile->Close();
  //
  NrecMapFile = new TFile("Results/PDC_10f6a_R2.root","READ");
  //NrecMapFile = new TFile("Results/NrecMapping_10f6.root","READ");// MC variation
  //NrecMapFile = new TFile("Results/NrecMapping_10f6a_FB5and7overlap.root","READ");
  MyList=(TList*)NrecMapFile->Get("MyList");
  NrecMap = (TH2D*)MyList->FindObject("fNchTrueDist");
  for(int bin=10; bin<=20; bin++){
    NrecMap->GetXaxis()->SetRangeUser(bin,bin);
    if(NrecMap->GetMean(2)>0) meanNchpp[bin-1] = log10(NrecMap->GetMean(2));
  }
  NrecMapFile->Close();
  
  //for(int i=0; i<20; i++) cout<<pow(10,meanNchPbPb[i])<<endl;
  //cout<<"+++++++++++++++"<<endl;
  //for(int i=0; i<20; i++) cout<<pow(10,meanNchpPb[i])<<endl;
  //cout<<"+++++++++++++++"<<endl;
  //for(int i=0; i<20; i++) cout<<pow(10,meanNchpp[i])<<endl;

  TFile *ExRangeFile=new TFile("Results/ExtendedQ3rangeM0.root","READ");
  TList *ExList=(TList*)ExRangeFile->Get("MyList");
  TH1D *ExRangeTerm1=(TH1D*)ExList->FindObject("fExtendedQ3Histo_term1");
  TH1D *ExRangeTerm2=(TH1D*)ExList->FindObject("fExtendedQ3Histo_term2");
  TH1D *ExRangeTerm5=(TH1D*)ExList->FindObject("fExtendedQ3Histo_term5");
  ExRangeTerm1->SetDirectory(0); ExRangeTerm2->SetDirectory(0); ExRangeTerm5->SetDirectory(0);
  ExRangeFile->Close();
  ExRangeTerm1->Sumw2(); ExRangeTerm2->Sumw2(); ExRangeTerm5->Sumw2();
  ExRangeTerm2->Scale(ExRangeTerm1->Integral(25,29)/ExRangeTerm2->Integral(25,29));
  ExRangeTerm5->Scale(ExRangeTerm1->Integral(25,29)/ExRangeTerm5->Integral(25,29));
  float TwoFrac=0.7;
  float OneFrac=pow(TwoFrac,.5); 
  float ThreeFrac=pow(TwoFrac,1.5);
  // Purify.  Isolate pure 3-pion QS correlations using Lambda and K3 (removes lower order correlations)
  ExRangeTerm1->Add(ExRangeTerm5, -(pow(1-OneFrac,3) + 3*OneFrac*pow(1-OneFrac,2))); 
  ExRangeTerm1->Add(ExRangeTerm2, -3*(1-OneFrac));
  ExRangeTerm1->Add(ExRangeTerm5, (1-OneFrac)*3*(1-TwoFrac));
  ExRangeTerm1->Scale(1/ThreeFrac);
  //N3_QS /= K3;
  // Isolate 3-pion cumulant
  ExRangeTerm1->Add(ExRangeTerm2, -3/TwoFrac);
  ExRangeTerm1->Add(ExRangeTerm5, 3*(1-TwoFrac)/TwoFrac);
  ExRangeTerm1->Add(ExRangeTerm5, 3);
  ExRangeTerm1->Divide(ExRangeTerm5);
  ExRangeTerm1->SetMarkerStyle(20);
  ExRangeTerm1->SetMarkerColor(1);
  
  //
  TF1 *MixedChargeSysFit=new TF1("MixedChargeSysFit","[0]+[1]*exp(-pow([2]*x/0.19733,2))",0,.5);
  //TFile *files_3[3][2][2][2][3][20];// CollType, Real/MonteCarlo, SC/MC, +/-, EDbin, MBINS
  //
  TF1 *Fit_C2[3][3][20];// CollType, EDbin, cb
  TH1D *Fit_h_C2[3][3][20];// CollType, EDbin, cb
  TH1D *Parameters_C2[3][3][5];// CollType, EDbin, Parameter#
  TH1D *C2[3][3][20];// CollType, EDbin, cb
  TH1D *C2_Sys[3][3][20];// CollType, EDbin, cb
  TH1D *C3[3][2][2][3][20];// CollType, Real/MonteCarlo, SC/MC, EDbin, cb
  TH1D *c3[3][2][2][3][20];// CollType, Real/MonteCarlo, SC/MC, EDbin, cb
  TH1D *C3_Sys[3][2][2][3][20];// CollType, Real/MonteCarlo, SC/MC, EDbin, cb
  TH1D *c3_Sys[3][2][2][3][20];// CollType, Real/MonteCarlo, SC/MC, EDbin, cb
  TH1D *c3_fit[3][2][3][20];// CollType, SC/MC, EDbin, cb
  TMinuit *Min[3][2][3][20];// CollType, +/-, EDbin, cb
  TH1D *Parameters_c3[3][3][5];// CollType, EDbin, Parameter#
  //char *labels[20]={"1050-2000","850-1050","700-850","600-700","500-600","400-500","320-400","260-320","200-260","150-200","100-150","70-100","50-70","40-50","30-40","20-30","15-20","10-15","5-10","0-5"};
  for(int ct=0; ct<3; ct++){
    for(int kt3=0; kt3<3; kt3++){
      for(int par=0; par<5; par++){
	TString *name_C2=new TString("Parameters_C2_");
	*name_C2 += ct;
	*name_C2 += kt3;
	*name_C2 += par;
	TString *name_c3=new TString("Parameters_c3_");
	*name_c3 += ct;
	*name_c3 += kt3;
	*name_c3 += par;
	//Parameters_c3[ct][par] = new TH1D(name->Data(),"",20,0.5,20.5);
	//for(int cb=0; cb<20; cb++) Parameters_c3[ct][par]->GetXaxis()->SetBinLabel(cb+1,labels[cb]);
	Parameters_C2[ct][kt3][par] = new TH1D(name_C2->Data(),"",2900,0.6,3.4);// 2900,0.6,3.4
	Parameters_c3[ct][kt3][par] = new TH1D(name_c3->Data(),"",2900,0.6,3.4);
      }
    }
  }
  //
  double N_1 = 0, N_1_e=0;
  double lambda_1 = 0, lambda_1_e=0;
  double radius_1 = 0, radius_1_e=0;
  double EW1_1 = 0, EW1_1_e=0;
  double EW2_1 = 0, EW2_1_e=0;
  double N_2 = 0, N_2_e=0;
  double lambda_2 = 0, lambda_2_e=0;
  double radius_2 = 0, radius_2_e=0;
  double EW1_2 = 0, EW1_2_e=0;
  double EW2_2 = 0, EW2_2_e=0;
 
  // Start File access
  for(int ct=0; ct<3; ct++){
    for(int dt=0; dt<1; dt++){// data type (Real or Monte Carlo)
      if(ct==0 && dt==1) continue; // no MC for PbPb
      for(int cb=0; cb<20; cb++){
	if(ct==0 && cb<MbinMinPbPb) continue;
	if(ct==0 && cb>MbinMaxPbPb) continue;
	if(ct==1 && cb>MbinMaxpPb) continue;
	if(ct==1 && cb<MbinMinpPb) continue;
	if(ct==2 && cb<MbinMinpp) continue;
	if(ct==2 && dt==1 && cb<=13) continue;// no Pythia data for cb=13
	for(int ChComb=0; ChComb<2; ChComb++) {// SC or MC
	  for(int ch=0; ch<1; ch++) {// - or +
	    for(int KT3=0; KT3<MaxKT3Bins; KT3++) {// Kt3 bin
	      if(dt==1 && KT3>0) continue;// no MC data yet for higher kt3
	      TString *name3 = new TString("OutFiles/OutFile");
	      if(ct==0) name3->Append("PbPb");
	      if(ct==1) name3->Append("pPb");
	      if(ct==2) name3->Append("pp");
	      if(dt==1) name3->Append("MonteCarlo");
	      if(ChComb==0) name3->Append("SC");
	      else name3->Append("MC");
	      if(ch==0) name3->Append("Neg");
	      else name3->Append("Pos");
	      if(IncludeEW && dt==0 && ChComb==0) name3->Append("EW");
	      name3->Append("Kt3_");
	      *name3 += KT3+1;
	      name3->Append("_M");
	      *name3 += cb;
	      name3->Append(".root");
	      //files_3[ct][dt][ChComb][ch][KT3][cb] = new TFile(name3->Data(),"READ");
	      TFile *file=new TFile(name3->Data(),"READ");
	      //
	      if(ChComb==0 && dt==0){
		Min[ct][ch][KT3][cb]=(TMinuit*)file->Get("MyMinuit_c3");
		//Min[ct][ch][KT3][cb]->SetDirectory(0);
	      }
	      
	      if(ch==0) {
		C3[ct][dt][ChComb][KT3][cb]=(TH1D*)file->Get("C3");
		C3[ct][dt][ChComb][KT3][cb]->SetDirectory(0);
		C3[ct][dt][ChComb][KT3][cb]->SetMarkerStyle(24+ct); 
		c3[ct][dt][ChComb][KT3][cb]=(TH1D*)file->Get("c3");
		c3[ct][dt][ChComb][KT3][cb]->SetDirectory(0);
		c3[ct][dt][ChComb][KT3][cb]->SetMarkerStyle(20+ct);
		if(dt==1) {c3[ct][dt][ChComb][KT3][cb]->SetMarkerStyle(28);}
		C3[ct][dt][ChComb][KT3][cb]->GetXaxis()->SetRangeUser(0,0.5);
		
		C3[ct][dt][ChComb][KT3][cb]->GetXaxis()->SetRangeUser(0,0.5);
		C3[ct][dt][ChComb][KT3][cb]->GetXaxis()->SetLabelFont(TextFont); C3[ct][dt][ChComb][KT3][cb]->GetYaxis()->SetLabelFont(TextFont); 
		C3[ct][dt][ChComb][KT3][cb]->GetXaxis()->SetTitleOffset(1.2);
		C3[ct][dt][ChComb][KT3][cb]->GetYaxis()->SetTitleOffset(1.2);
		C3[ct][dt][ChComb][KT3][cb]->GetYaxis()->SetTitleFont(TextFont);
		c3[ct][dt][ChComb][KT3][cb]->GetXaxis()->SetRangeUser(0,0.5);
		c3[ct][dt][ChComb][KT3][cb]->GetXaxis()->SetLabelFont(TextFont); c3[ct][dt][ChComb][KT3][cb]->GetYaxis()->SetLabelFont(TextFont); 
		c3[ct][dt][ChComb][KT3][cb]->GetXaxis()->SetTitleOffset(1.2);
		c3[ct][dt][ChComb][KT3][cb]->GetYaxis()->SetTitleOffset(1.2);
		c3[ct][dt][ChComb][KT3][cb]->GetYaxis()->SetTitleFont(TextFont);
		if(ct==0) {
		  C3[ct][dt][ChComb][KT3][cb]->SetMarkerColor(1); C3[ct][dt][ChComb][KT3][cb]->SetLineColor(1);
		  c3[ct][dt][ChComb][KT3][cb]->SetMarkerColor(1); c3[ct][dt][ChComb][KT3][cb]->SetLineColor(1);
		}else if(ct==1){
		  C3[ct][dt][ChComb][KT3][cb]->SetMarkerColor(2); C3[ct][dt][ChComb][KT3][cb]->SetLineColor(2);
		  c3[ct][dt][ChComb][KT3][cb]->SetMarkerColor(2); c3[ct][dt][ChComb][KT3][cb]->SetLineColor(2);
		  if(dt==1) {c3[ct][dt][ChComb][KT3][cb]->SetMarkerColor(1); c3[ct][dt][ChComb][KT3][cb]->SetLineColor(1);}
		}else {
		  C3[ct][dt][ChComb][KT3][cb]->SetMarkerColor(4); C3[ct][dt][ChComb][KT3][cb]->SetLineColor(4);
		  c3[ct][dt][ChComb][KT3][cb]->SetMarkerColor(4); c3[ct][dt][ChComb][KT3][cb]->SetLineColor(4);
		  if(dt==1) {c3[ct][dt][ChComb][KT3][cb]->SetMarkerColor(1); c3[ct][dt][ChComb][KT3][cb]->SetLineColor(1);}
		}
		if(dt==0) {
		  if(ChComb==0) {
		    C2[ct][KT3][cb]=(TH1D*)file->Get("C2_ss");
		    Fit_h_C2[ct][KT3][cb]=(TH1D*)file->Get("fitC2ss_h");
		    Fit_C2[ct][KT3][cb]=(TF1*)file->Get("fitC2ss");
		    C2_Sys[ct][KT3][cb]=(TH1D*)C2[ct][KT3][cb]->Clone();
		    if(ct==0){
		      C2[ct][KT3][cb]->SetMarkerColor(1); C2[ct][KT3][cb]->SetLineColor(1);
		      Fit_h_C2[ct][KT3][cb]->SetLineColor(1);
		      C2_Sys[ct][KT3][cb]->SetFillColor(kGray);
		    }else if(ct==1){
		      C2[ct][KT3][cb]->SetMarkerColor(2); C2[ct][KT3][cb]->SetLineColor(2);
		      Fit_h_C2[ct][KT3][cb]->SetLineColor(2);
		      C2_Sys[ct][KT3][cb]->SetFillColor(kRed-10);
		    }else{
		      C2[ct][KT3][cb]->SetMarkerColor(4); C2[ct][KT3][cb]->SetLineColor(4);
		      Fit_h_C2[ct][KT3][cb]->SetLineColor(4);
		      C2_Sys[ct][KT3][cb]->SetFillColor(kBlue-10);
		    }
		    C2_Sys[ct][KT3][cb]->SetMarkerSize(0); C2_Sys[ct][KT3][cb]->SetFillStyle(1000);
		    // C2 Systematics
		    for(int bin=1; bin<=C2_Sys[ct][KT3][cb]->GetNbinsX(); bin++){
		      C2_Sys[ct][KT3][cb]->SetBinError(bin, 0.01*C2_Sys[ct][KT3][cb]->GetBinContent(bin));
		    }
		    C2[ct][KT3][cb]->SetDirectory(0); Fit_h_C2[ct][KT3][cb]->SetDirectory(0);
		    C2_Sys[ct][KT3][cb]->SetDirectory(0);
		    C2[ct][KT3][cb]->SetMarkerStyle(24+ct);
		    C2[ct][KT3][cb]->GetXaxis()->SetRangeUser(0,0.33);
		    C2[ct][KT3][cb]->GetXaxis()->SetLabelFont(TextFont); C2[ct][KT3][cb]->GetYaxis()->SetLabelFont(TextFont); 
		    C2[ct][KT3][cb]->GetXaxis()->SetTitleOffset(1.2); C2[ct][KT3][cb]->GetYaxis()->SetTitleOffset(1.2);
		    C2[ct][KT3][cb]->GetXaxis()->SetTitleFont(TextFont); C2[ct][KT3][cb]->GetYaxis()->SetTitleFont(TextFont);
		    C2[ct][KT3][cb]->GetXaxis()->SetTitleSize(SizeTitle*SF2); //C2[ct][KT3][cb]->GetYaxis()->SetTitleFont(SizeTitle*SF2);
		    //C2[ct][KT3][cb]->GetYaxis()->SetTitle("C_{2}^{#pm#pm}");

		  }
		  c3_fit[ct][ChComb][KT3][cb]=(TH1D*)file->Get("dentimesFit_c3");
		  c3_fit[ct][ChComb][KT3][cb]->SetDirectory(0);
		  if(ct==0) c3_fit[ct][ChComb][KT3][cb]->SetLineColor(1);
		  if(ct==1) c3_fit[ct][ChComb][KT3][cb]->SetLineColor(2);
		  if(ct==2) c3_fit[ct][ChComb][KT3][cb]->SetLineColor(4);
		}
	      }else{
		TH1D *tempC3=(TH1D*)file->Get("C3");
		TH1D *tempc3=(TH1D*)file->Get("c3");
		TH1D *tempc3_fit=(TH1D*)file->Get("dentimesFit_c3");
		if(dt==0) c3_fit[ct][ChComb][KT3][cb]->Add(tempc3_fit);
		if(ChComb==0){
		  N_1 = 0; N_1_e=0;
		  lambda_1 = 0; lambda_1_e=0;
		  radius_1 = 0; radius_1_e=0;
		  EW1_1 = 0; EW1_1_e=0;
		  EW2_1 = 0; EW2_1_e=0;
		  N_2 = 0; N_2_e=0;
		  lambda_2 = 0; lambda_2_e=0;
		  radius_2 = 0; radius_2_e=0;
		  EW1_2 = 0; EW1_2_e=0;
		  EW2_2 = 0; EW2_2_e=0;
		  //
		  if(!AddedCC){
		    Min[ct][0][KT3][cb]->GetParameter(0,N_1,N_1_e); Min[ct][1][KT3][cb]->GetParameter(0,N_2,N_2_e);
		    Min[ct][0][KT3][cb]->GetParameter(1,lambda_1,lambda_1_e); Min[ct][1][KT3][cb]->GetParameter(1,lambda_2,lambda_2_e);
		    Min[ct][0][KT3][cb]->GetParameter(2,radius_1,radius_1_e); Min[ct][1][KT3][cb]->GetParameter(2,radius_2,radius_2_e);
		    Min[ct][0][KT3][cb]->GetParameter(3,EW1_1,EW1_1_e); Min[ct][1][KT3][cb]->GetParameter(3,EW1_2,EW1_2_e);
		    Min[ct][0][KT3][cb]->GetParameter(4,EW2_1,EW2_1_e); Min[ct][1][KT3][cb]->GetParameter(4,EW2_2,EW2_2_e);
		    Parameters_c3[ct][KT3][0]->SetBinContent(cb+1, (N_1+N_2)/2.); Parameters_c3[ct][KT3][0]->SetBinError(cb+1, sqrt(pow(N_1_e,2)+pow(N_2_e,2))/2.);
		    Parameters_c3[ct][KT3][1]->SetBinContent(cb+1, (lambda_1+lambda_2)/2.); Parameters_c3[ct][KT3][1]->SetBinError(cb+1, sqrt(pow(lambda_1_e,2)+pow(lambda_2_e,2))/2.);
		    Parameters_c3[ct][KT3][2]->SetBinContent(cb+1, (radius_1+radius_2)/2.); Parameters_c3[ct][KT3][2]->SetBinError(cb+1, sqrt(pow(radius_1_e,2)+pow(radius_2_e,2))/2.);
		    Parameters_c3[ct][KT3][3]->SetBinContent(cb+1, (EW1_1+EW1_2)/2.); Parameters_c3[ct][KT3][3]->SetBinError(cb+1, sqrt(pow(EW1_1_e,2)+pow(EW1_2_e,2))/2.);
		    Parameters_c3[ct][KT3][4]->SetBinContent(cb+1, (EW2_1+EW2_2)/2.); Parameters_c3[ct][KT3][4]->SetBinError(cb+1, sqrt(pow(EW2_1_e,2)+pow(EW2_2_e,2))/2.);
		  }
		}
		//
		if(!AddedCC){
		  for(int bin=1; bin<=tempC3->GetNbinsX(); bin++){
		    double valueC3 = (C3[ct][dt][ChComb][KT3][cb]->GetBinContent(bin) + tempC3->GetBinContent(bin))/2.;
		    double valueC3_e = sqrt(pow(C3[ct][dt][ChComb][KT3][cb]->GetBinError(bin),2)+pow(tempC3->GetBinError(bin),2))/sqrt(2.);
		    double valuec3 = (c3[ct][dt][ChComb][KT3][cb]->GetBinContent(bin) + tempc3->GetBinContent(bin))/2.;
		    double valuec3_e = sqrt(pow(c3[ct][dt][ChComb][KT3][cb]->GetBinError(bin),2)+pow(tempc3->GetBinError(bin),2))/sqrt(2.);
		    C3[ct][dt][ChComb][KT3][cb]->SetBinContent(bin, valueC3);
		    C3[ct][dt][ChComb][KT3][cb]->SetBinError(bin, valueC3_e);
		    c3[ct][dt][ChComb][KT3][cb]->SetBinContent(bin, valuec3);
		    c3[ct][dt][ChComb][KT3][cb]->SetBinError(bin, valuec3_e);
		    c3[ct][dt][ChComb][KT3][cb]->SetMarkerStyle(20); c3[ct][dt][ChComb][KT3][cb]->SetMarkerColor(4); c3[ct][dt][ChComb][KT3][cb]->SetLineColor(4);
		    c3[ct][dt][ChComb][KT3][cb]->GetYaxis()->SetTitleOffset(1.4);
		    if(ChComb==1){
		      c3[ct][dt][ChComb][KT3][cb]->SetMinimum(0.96);// 0.99
		      c3[ct][dt][ChComb][KT3][cb]->SetMaximum(1.02);// 1.02
		      c3[ct][dt][ChComb][KT3][cb]->GetXaxis()->SetRangeUser(0,0.13);
		      c3[ct][dt][ChComb][KT3][cb]->GetYaxis()->SetTitleOffset(1.1);
		      c3[ct][dt][ChComb][KT3][cb]->GetYaxis()->SetTitleSize(.05);
		    }
		  }
		}
	      }// ch==1
	      for(int bin=1; bin<10; bin++){// Remove large error bins
		if(C3[ct][dt][ChComb][KT3][cb]->GetBinError(bin) > 0.5*C3[ct][dt][ChComb][KT3][cb]->GetBinContent(bin)){
		  C3[ct][dt][ChComb][KT3][cb]->SetBinContent(bin,10); C3[ct][dt][ChComb][KT3][cb]->SetBinError(bin,10);
		}
		if(c3[ct][dt][ChComb][KT3][cb]->GetBinError(bin) > 0.5*c3[ct][dt][ChComb][KT3][cb]->GetBinContent(bin)){
		  c3[ct][dt][ChComb][KT3][cb]->SetBinContent(bin,10); c3[ct][dt][ChComb][KT3][cb]->SetBinError(bin,10);
		}
	      }
	      if(AddedCC && dt==0){
		Min[ct][0][KT3][cb]->GetParameter(0,N_1,N_1_e);
		Min[ct][0][KT3][cb]->GetParameter(1,lambda_1,lambda_1_e);
		Min[ct][0][KT3][cb]->GetParameter(2,radius_1,radius_1_e);
		Min[ct][0][KT3][cb]->GetParameter(3,EW1_1,EW1_1_e);
		Min[ct][0][KT3][cb]->GetParameter(4,EW2_1,EW2_1_e);
		//Parameters_c3[ct][0]->SetBinContent(cb+1, N_1); Parameters_c3[ct][0]->SetBinError(cb+1, N_1_e);
		//Parameters_c3[ct][1]->SetBinContent(cb+1, lambda_1); Parameters_c3[ct][1]->SetBinError(cb+1, lambda_1_e);
		//Parameters_c3[ct][2]->SetBinContent(cb+1, radius_1); Parameters_c3[ct][2]->SetBinError(cb+1, radius_1_e);
		//Parameters_c3[ct][3]->SetBinContent(cb+1, EW1_1); Parameters_c3[ct][3]->SetBinError(cb+1, EW1_1_e);
		//Parameters_c3[ct][4]->SetBinContent(cb+1, EW2_1); Parameters_c3[ct][4]->SetBinError(cb+1, EW2_1_e);
		
		double logNch=0;
		if(ct==0) logNch=meanNchPbPb[cb];
		else if(ct==1) logNch=meanNchpPb[cb];
		else logNch=meanNchpp[cb];
		int logNchBin = Parameters_c3[ct][KT3][0]->GetXaxis()->FindBin(logNch);
		Parameters_c3[ct][KT3][0]->SetBinContent(logNchBin, N_1); Parameters_c3[ct][KT3][0]->SetBinError(logNchBin, N_1_e);
		Parameters_c3[ct][KT3][1]->SetBinContent(logNchBin, lambda_1); Parameters_c3[ct][KT3][1]->SetBinError(logNchBin, lambda_1_e);
		Parameters_c3[ct][KT3][2]->SetBinContent(logNchBin, radius_1); Parameters_c3[ct][KT3][2]->SetBinError(logNchBin, radius_1_e);
		Parameters_c3[ct][KT3][3]->SetBinContent(logNchBin, EW1_1); Parameters_c3[ct][KT3][3]->SetBinError(logNchBin, EW1_1_e);
		Parameters_c3[ct][KT3][4]->SetBinContent(logNchBin, EW2_1); Parameters_c3[ct][KT3][4]->SetBinError(logNchBin, EW2_1_e);
		//
		Parameters_C2[ct][KT3][0]->SetBinContent(logNchBin, Fit_C2[ct][KT3][cb]->GetParameter(0));// N
		Parameters_C2[ct][KT3][0]->SetBinError(logNchBin, Fit_C2[ct][KT3][cb]->GetParError(0));
		Parameters_C2[ct][KT3][1]->SetBinContent(logNchBin, Fit_C2[ct][KT3][cb]->GetParameter(1));// lambda
		Parameters_C2[ct][KT3][1]->SetBinError(logNchBin, Fit_C2[ct][KT3][cb]->GetParError(1));
		Parameters_C2[ct][KT3][2]->SetBinContent(logNchBin, Fit_C2[ct][KT3][cb]->GetParameter(3));// R
		Parameters_C2[ct][KT3][2]->SetBinError(logNchBin, Fit_C2[ct][KT3][cb]->GetParError(3));
		Parameters_C2[ct][KT3][3]->SetBinContent(logNchBin, Fit_C2[ct][KT3][cb]->GetParameter(5));// kappa3
		Parameters_C2[ct][KT3][3]->SetBinError(logNchBin, Fit_C2[ct][KT3][cb]->GetParError(5));
		Parameters_C2[ct][KT3][4]->SetBinContent(logNchBin, Fit_C2[ct][KT3][cb]->GetParameter(6));// kappa4
		Parameters_C2[ct][KT3][4]->SetBinError(logNchBin, Fit_C2[ct][KT3][cb]->GetParError(6));
		//
		// Sys errors
		C3_Sys[ct][0][ChComb][KT3][cb] = (TH1D*)C3[ct][0][ChComb][KT3][cb]->Clone();
		c3_Sys[ct][0][ChComb][KT3][cb] = (TH1D*)c3[ct][0][ChComb][KT3][cb]->Clone();
		if(ct==0){
		  C3_Sys[ct][0][ChComb][KT3][cb]->SetMarkerSize(0); C3_Sys[ct][0][ChComb][KT3][cb]->SetFillStyle(1000); C3_Sys[ct][0][ChComb][KT3][cb]->SetFillColor(kGray);
		  c3_Sys[ct][0][ChComb][KT3][cb]->SetMarkerSize(0); c3_Sys[ct][0][ChComb][KT3][cb]->SetFillStyle(1000); c3_Sys[ct][0][ChComb][KT3][cb]->SetFillColor(kGray);
		}else if(ct==1){
		  C3_Sys[ct][0][ChComb][KT3][cb]->SetMarkerSize(0); C3_Sys[ct][0][ChComb][KT3][cb]->SetFillStyle(1000); C3_Sys[ct][0][ChComb][KT3][cb]->SetFillColor(kRed-10);
		  c3_Sys[ct][0][ChComb][KT3][cb]->SetMarkerSize(0); c3_Sys[ct][0][ChComb][KT3][cb]->SetFillStyle(1000); c3_Sys[ct][0][ChComb][KT3][cb]->SetFillColor(kRed-10);
		}else {
		  C3_Sys[ct][0][ChComb][KT3][cb]->SetMarkerSize(0); C3_Sys[ct][0][ChComb][KT3][cb]->SetFillStyle(1000); C3_Sys[ct][0][ChComb][KT3][cb]->SetFillColor(kBlue-10);
		  c3_Sys[ct][0][ChComb][KT3][cb]->SetMarkerSize(0); c3_Sys[ct][0][ChComb][KT3][cb]->SetFillStyle(1000); c3_Sys[ct][0][ChComb][KT3][cb]->SetFillColor(kBlue-10);
		}
		
		if(ChComb==1){
		  MixedChargeSysFit->SetParameter(0,1);
		  MixedChargeSysFit->SetParameter(1,.1);
		  MixedChargeSysFit->SetParameter(2,1);
		  c3[ct][0][ChComb][KT3][cb]->Fit(MixedChargeSysFit,"IMNQ","",0.01,0.5);
		  for(int i=1; i<=c3[ct][0][ChComb][KT3][cb]->GetNbinsX(); i++) {
		    float Q3=(i-0.5)*0.01;
		    // SameCharge
		    C3_Sys[ct][0][0][KT3][cb]->SetBinError(i, 0.02 * C3_Sys[ct][0][0][KT3][cb]->GetBinContent(i));
		    c3_Sys[ct][0][0][KT3][cb]->SetBinError(i, MixedChargeSysFit->Eval(Q3)-1.0);
		    // MixedCharge
		    C3_Sys[ct][0][1][KT3][cb]->SetBinError(i, 0.02 * C3_Sys[ct][0][1][KT3][cb]->GetBinContent(i));
		    c3_Sys[ct][0][1][KT3][cb]->SetBinError(i, 0.02 * c3_Sys[ct][0][1][KT3][cb]->GetBinContent(i));
		  }
		}
		C3_Sys[ct][0][ChComb][KT3][cb]->SetDirectory(0); c3_Sys[ct][0][ChComb][KT3][cb]->SetDirectory(0);
	      }// AddedCC and dt==0
	      file->Close();
	      
	    }// Kt3
	  }// ch
	}// ChComb
	
      }// cb
      
      if(dt==0){
	for(int KT3=0; KT3<3; KT3++){
	  for(int par=0; par<5; par++){
	    if(ct<2){
	      Parameters_C2[ct][KT3][par]->SetMarkerStyle(24+ct);
	      Parameters_c3[ct][KT3][par]->SetMarkerStyle(20+ct);
	    }else{
	      Parameters_C2[ct][KT3][par]->SetMarkerStyle(28);
	      Parameters_c3[ct][KT3][par]->SetMarkerStyle(34);
	    }

	    if(ct==0){
	      Parameters_C2[ct][KT3][par]->SetMarkerColor(1); 
	      Parameters_C2[ct][KT3][par]->SetLineColor(1);
	      Parameters_c3[ct][KT3][par]->SetMarkerColor(1); 
	      Parameters_c3[ct][KT3][par]->SetLineColor(1);
	    }else if(ct==1){
	      Parameters_C2[ct][KT3][par]->SetMarkerColor(2); 
	      Parameters_C2[ct][KT3][par]->SetLineColor(2);
	      Parameters_c3[ct][KT3][par]->SetMarkerColor(2); 
	      Parameters_c3[ct][KT3][par]->SetLineColor(2);
	    }else {
	      Parameters_C2[ct][KT3][par]->SetMarkerColor(4); 
	      Parameters_C2[ct][KT3][par]->SetLineColor(4);
	      Parameters_c3[ct][KT3][par]->SetMarkerColor(4); 
	      Parameters_c3[ct][KT3][par]->SetLineColor(4);
	    }
	    
	    if(par==0) Parameters_c3[ct][KT3][par]->GetYaxis()->SetTitle("N");
	    //if(par==1) Parameters_c3[ct][KT3][par]->GetYaxis()->SetTitle("#lambda_{3}");
	    if(par==1) Parameters_c3[ct][KT3][par]->GetYaxis()->SetTitle("#lambda or #lambda_{3}");
	    //if(par==2) Parameters_c3[ct][KT3][par]->GetYaxis()->SetTitle("R_{inv,3} (fm)");
	    if(par==2) Parameters_c3[ct][KT3][par]->GetYaxis()->SetTitle("R_{inv,2} or R_{inv,3} (fm)");
	    if(par==3) Parameters_c3[ct][KT3][par]->GetYaxis()->SetTitle("#kappa_{3}");
	    if(par==4) Parameters_c3[ct][KT3][par]->GetYaxis()->SetTitle("#kappa_{4}");
	    Parameters_c3[ct][KT3][par]->GetXaxis()->SetTitle("log_{10}(N_{ch})");


	  }// par
	}// KT3
      }
    }// dt
  }// ct
  
 
  TF1 *Unity = new TF1("Unity","1",0,100);
  Unity->SetLineStyle(2);
  Unity->SetLineColor(1);
  
 

  
  ////////////////////////////////////////////////////
  ////////////////////////////////////////////////////
  // Progaganda Plot
  /*
  TCanvas *can2 = new TCanvas("can2", "can2",10,0,600,600);// 11,53,700,500
  can2->SetHighLightColor(2);
  gStyle->SetOptFit(0111);
  can2->SetFillColor(10);//10
  can2->SetBorderMode(0);
  can2->SetBorderSize(2);
  can2->SetFrameFillColor(0);
  can2->SetFrameBorderMode(0);
  can2->SetFrameBorderMode(0);
  can2->cd();
  TPad *pad2 = new TPad("pad2","pad2",0.0,0.0,1.,1.);
  gPad->SetGridx(0);
  gPad->SetGridy(0);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  pad2->SetTopMargin(0.02);//0.05
  pad2->SetRightMargin(0.01);//3e-2
  pad2->SetBottomMargin(0.07);//0.12
  pad2->Draw();
  pad2->cd();
  TLegend *legend3 = new TLegend(.35,.75, .94,.95,NULL,"brNDC");//.45 or .4 for x1
  legend3->SetBorderSize(0);
  legend3->SetFillColor(0);
  legend3->SetTextFont(TextFont);
  legend3->SetTextSize(SizeLegend*SF2);

  gPad->SetRightMargin(0.03); gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.1); gPad->SetTopMargin(0.02);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  int Mb_pp=18;
  int Mb_pPb=13;
  int Mb_PbPb=0;
  c3[2][0][0][0][Mb_pp]->GetXaxis()->SetLabelSize(SizeLabel*SF2); c3[2][0][0][0][Mb_pp]->GetYaxis()->SetLabelSize(SizeLabel*SF2);
  c3[2][0][0][0][Mb_pp]->GetXaxis()->SetNdivisions(808);
  //
  c3[2][0][0][0][Mb_pp]->GetYaxis()->SetTitle("#font[12]{#bf{c}}_{3}");
  c3[2][0][0][0][Mb_pp]->GetXaxis()->SetTitle("#font[12]{Q}_{3} (GeV/#font[12]{c})");
  c3[2][0][0][0][Mb_pp]->GetYaxis()->SetTitleOffset(1.1);
  c3[2][0][0][0][Mb_pp]->GetYaxis()->SetTitleSize(0.05);
  c3[2][0][0][0][Mb_pp]->SetMinimum(0.9);
  c3[2][0][0][0][Mb_pp]->Draw();
  //legend3->AddEntry(c3[2][0][0][Mb_pp],"N_{ch} = 9#pm0.2, pp #sqrt{s}=7 TeV","p");
  //legend3->AddEntry(c3[1][0][0][Mb_pPb],"N_{ch} = 59#pm2, p-Pb #sqrt{s_{NN}}=5.02 TeV","p");
  //legend3->AddEntry(c3[0][0][0][Mb_PbPb],"N_{ch} = 1922#pm135, Pb-Pb #sqrt{s_{NN}}=2.76 TeV","p");
  legend3->AddEntry(c3[2][0][0][0][Mb_pp],"Low N_{ch}, pp #sqrt{#font[12]{s}}=7 TeV","p");
  legend3->AddEntry(c3[1][0][0][0][Mb_pPb],"Mid N_{ch}, p-Pb #sqrt{#font[12]{s_{NN}}}=5.02 TeV","p");
  legend3->AddEntry(c3[0][0][0][0][Mb_PbPb],"High N_{ch}, Pb-Pb #sqrt{#font[12]{s_{NN}}}=2.76 TeV","p");
  //
  //TH1D *FillerHist = new TH1D("FillerHist","",50,0,0.5);
  for(int i=1; i<=50; i++) {
    if(i<13) {
      ExRangeTerm1->SetBinContent(i,10.0);
      //FillerHist->SetBinContent(i,10.0); FillerHist->SetBinError(i,.0001);
    }else {
      //FillerHist->SetBinContent(i,1.0); FillerHist->SetBinError(i,.0001);
      c3[0][0][0][0][Mb_PbPb]->SetBinContent(i,10.0);
      c3_Sys[0][0][0][0][Mb_PbPb]->SetBinContent(i,10.0);
    }
  }
  //FillerHist->SetMarkerStyle(24);
  //
  c3_Sys[2][0][0][0][Mb_pp]->Draw("E2 same");
  c3_Sys[1][0][0][0][Mb_pPb]->Draw("E2 same");
  c3_Sys[0][0][0][0][Mb_PbPb]->Draw("E2 same");
  c3[2][0][0][0][Mb_pp]->Draw("same");
  c3[1][0][0][0][Mb_pPb]->Draw("same");
  c3[0][0][0][0][Mb_PbPb]->Draw("same");
  //FillerHist->Draw("same");
  ExRangeTerm1->Draw("same");
  
  TLatex *Specif_Kt3_1 = new TLatex(0.5,0.59,"0.16<#font[12]{K}_{T,3}<1.0 GeV/#font[12]{c}");
  Specif_Kt3_1->SetNDC();
  Specif_Kt3_1->SetTextFont(TextFont);
  Specif_Kt3_1->SetTextSize(SizeSpecif*SF2);
  Specif_Kt3_1->Draw("same");
  TLatex *Specif_Kinematics_1 = new TLatex(0.5,0.52,"|#eta|<0.8, 0.16<#font[12]{p}_{T}<1.0 GeV/#font[12]{c}");
  Specif_Kinematics_1->SetNDC();
  Specif_Kinematics_1->SetTextFont(TextFont);
  Specif_Kinematics_1->SetTextSize(SizeSpecif*SF2);
  Specif_Kinematics_1->Draw("same");
  TLatex *Species_1 = new TLatex(0.5,0.45,"#pi^{+}#pi^{+}#pi^{+} & #pi^{-}#pi^{-}#pi^{-} combined");
  Species_1->SetNDC();
  Species_1->SetTextFont(TextFont);
  Species_1->SetTextSize(SizeSpecif*SF2);
  Species_1->Draw("same");

  legend3->Draw("same");
  DrawALICELogo(kTRUE, .75, .25, .9, .4);

  if(SaveFiles) can2->SaveAs("ThreePionCorrelation_Evolution.eps");
  */
  ////////////////////////////////////////////////////
  ////////////////////////////////////////////////////
  // Radii and lambda plot
  
  TCanvas *can3 = new TCanvas("can3", "can3",700,0,600,900);// 11,53,700,500
  can3->SetHighLightColor(2);
  gStyle->SetOptFit(0111);
  can3->SetFillColor(10);//10
  can3->SetBorderMode(0);
  can3->SetBorderSize(2);
  can3->SetFrameFillColor(0);
  can3->SetFrameBorderMode(0);
  can3->SetFrameBorderMode(0);
  can3->cd();
  TPad *pad3 = new TPad("pad3","pad3",0.0,0.0,1.,1.);
  gPad->SetTickx();
  gPad->SetTicky();
  pad3->SetTopMargin(0.02);//0.05
  pad3->SetRightMargin(0.01);//3e-2
  pad3->SetBottomMargin(0.07);//0.12
  pad3->Divide(1,2,0,0);
  pad3->Draw();
  pad3->cd();
  TLegend *legend4 = new TLegend(.15,.75, .55,.95,NULL,"brNDC");//.45 or .4 for x1
  legend4->SetBorderSize(0);
  legend4->SetFillColor(0);
  legend4->SetTextFont(TextFont);
  legend4->SetTextSize(SizeLegend*SF1);
  //
  pad3->cd(1);
  gPad->SetRightMargin(RightMargin); gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.1);
  gPad->SetTickx(); gPad->SetTicky();
  //gPad->SetGridx(); gPad->SetGridy();
  TH1D *RadiiPbPb=(TH1D*)Parameters_c3[0][KT3Bin][2]->Clone();
  TH1D *RadiipPb=(TH1D*)Parameters_c3[1][KT3Bin][2]->Clone();
  TH1D *Radiipp=(TH1D*)Parameters_c3[2][KT3Bin][2]->Clone();
  RadiiPbPb->GetXaxis()->SetLabelFont(TextFont); RadiiPbPb->GetYaxis()->SetLabelFont(TextFont); 
  RadiiPbPb->GetXaxis()->SetLabelSize(SizeLabel*SF1); RadiiPbPb->GetYaxis()->SetLabelSize(SizeLabel*SF1);
  RadiiPbPb->GetXaxis()->SetNdivisions(808);
  RadiiPbPb->GetXaxis()->SetTitleOffset(1.3);
  RadiiPbPb->GetYaxis()->SetTitleOffset(1.1);
  RadiiPbPb->GetYaxis()->SetTitleFont(TextFont); RadiiPbPb->GetYaxis()->SetTitleSize(SizeTitle*SF1);
  if(IncludeEW) RadiiPbPb->GetYaxis()->SetTitle("R_{inv,3}^{Ew} (fm)");
  RadiiPbPb->SetMinimum(0); RadiiPbPb->SetMaximum(11.5);
  
  RadiiPbPb->Draw();
  //
  double xAxis_e[20]={0};
  double yAxisPbPb[20]={0};
  double yAxisPbPb_e[20]={0};
  double yAxispPb[20]={0};
  double yAxispPb_e[20]={0};
  double yAxispp[20]={0};
  double yAxispp_e[20]={0};
 
  
  TF1 *QrangeSys_c3=new TF1("QrangeSys_c3","0.1 - .07*(x/3.3)",0,4);
  for(int cb=0; cb<20; cb++){
    int binPbPb = RadiiPbPb->GetXaxis()->FindBin(meanNchPbPb[cb]);
    int binpPb = RadiipPb->GetXaxis()->FindBin(meanNchpPb[cb]);
    int binpp = Radiipp->GetXaxis()->FindBin(meanNchpp[cb]);
    double RsysPbPb = sqrt(pow(QrangeSys_c3->Eval(meanNchPbPb[cb])*RadiiPbPb->GetBinContent(binPbPb),2) + pow(0.05,2));
    double RsyspPb = sqrt(pow(QrangeSys_c3->Eval(meanNchpPb[cb])*RadiipPb->GetBinContent(binpPb),2) + pow(0.05,2));
    double Rsyspp = sqrt(pow(QrangeSys_c3->Eval(meanNchpp[cb])*Radiipp->GetBinContent(binpp),2) + pow(0.05,2));
    if(RadiiPbPb->GetBinContent(binPbPb)==0) {yAxisPbPb[cb]=100; yAxisPbPb_e[cb]=100;}
    else {yAxisPbPb[cb]=RadiiPbPb->GetBinContent(binPbPb); yAxisPbPb_e[cb]=RsysPbPb;}
    if(RadiipPb->GetBinContent(binpPb)==0) {yAxispPb[cb]=100; yAxispPb_e[cb]=100;}
    else {yAxispPb[cb]=RadiipPb->GetBinContent(binpPb); yAxispPb_e[cb]=RsyspPb;}
    if(Radiipp->GetBinContent(binpp)==0) {yAxispp[cb]=100; yAxispp_e[cb]=100;}
    else {yAxispp[cb]=Radiipp->GetBinContent(binpp); yAxispp_e[cb]=Rsyspp;}
    // Nch systematics
    // old method: log10((0.93 + 0.05*(cb/19))*pow(10,<Nch>)), 7% discrepancy max for central PbPb, 2% for pp and pPb.
    // Now 5% for all bins
    if(cb<13) xAxis_e[cb]=fabs(meanNchPbPb[cb] - log10((0.95)*pow(10,meanNchPbPb[cb])));
    else xAxis_e[cb]=fabs(meanNchpPb[cb] - log10((0.95)*pow(10,meanNchpPb[cb])));
  }
  TGraphErrors *grRadiiSys_PbPb=new TGraphErrors(20,meanNchPbPb,yAxisPbPb,xAxis_e,yAxisPbPb_e);
  TGraphErrors *grRadiiSys_pPb=new TGraphErrors(20,meanNchpPb,yAxispPb,xAxis_e,yAxispPb_e);
  TGraphErrors *grRadiiSys_pp=new TGraphErrors(20,meanNchpp,yAxispp,xAxis_e,yAxispp_e);
  grRadiiSys_pp->SetMarkerSize(0); grRadiiSys_pp->SetFillStyle(1000); grRadiiSys_pp->SetFillColor(kBlue-10);
  grRadiiSys_pPb->SetMarkerSize(0); grRadiiSys_pPb->SetFillStyle(1000); grRadiiSys_pPb->SetFillColor(kRed-10);
  grRadiiSys_PbPb->SetMarkerSize(0); grRadiiSys_PbPb->SetFillStyle(1000); grRadiiSys_PbPb->SetFillColor(kGray);
  
  // C2 
  TH1D *RadiiC2PbPb=(TH1D*)Parameters_C2[0][KT3Bin][2]->Clone();
  TH1D *RadiiC2pPb=(TH1D*)Parameters_C2[1][KT3Bin][2]->Clone();
  TH1D *RadiiC2pp=(TH1D*)Parameters_C2[2][KT3Bin][2]->Clone();
  for(int cb=0; cb<20; cb++){
    int binPbPb = RadiiC2PbPb->GetXaxis()->FindBin(meanNchPbPb[cb]);
    int binpPb = RadiiC2pPb->GetXaxis()->FindBin(meanNchpPb[cb]);
    int binpp = RadiiC2pp->GetXaxis()->FindBin(meanNchpp[cb]);
    double RsysPbPb = sqrt(pow(0.08*RadiiC2PbPb->GetBinContent(binPbPb),2) + pow(0.02*RadiiC2PbPb->GetBinContent(binPbPb),2));// 8% fit range, 2% MRC
    double RsyspPb = sqrt(pow(0.08*RadiiC2pPb->GetBinContent(binpPb),2) + pow(0.02*RadiiC2pPb->GetBinContent(binpPb),2));
    double Rsyspp = sqrt(pow(0.08*RadiiC2pp->GetBinContent(binpp),2) + pow(0.02*RadiiC2pp->GetBinContent(binpp),2));
    if(RadiiC2PbPb->GetBinContent(binPbPb)==0) {yAxisPbPb[cb]=100; yAxisPbPb_e[cb]=100;}
    else {yAxisPbPb[cb]=RadiiC2PbPb->GetBinContent(binPbPb); yAxisPbPb_e[cb]=RsysPbPb;}
    if(RadiiC2pPb->GetBinContent(binpPb)==0) {yAxispPb[cb]=100; yAxispPb_e[cb]=100;}
    else {yAxispPb[cb]=RadiiC2pPb->GetBinContent(binpPb); yAxispPb_e[cb]=RsyspPb;}
    if(RadiiC2pp->GetBinContent(binpp)==0) {yAxispp[cb]=100; yAxispp_e[cb]=100;}
    else {yAxispp[cb]=RadiiC2pp->GetBinContent(binpp); yAxispp_e[cb]=Rsyspp;}
  }
  TGraphErrors *grRadiiC2Sys_PbPb=new TGraphErrors(20,meanNchPbPb,yAxisPbPb,xAxis_e,yAxisPbPb_e);
  TGraphErrors *grRadiiC2Sys_pPb=new TGraphErrors(20,meanNchpPb,yAxispPb,xAxis_e,yAxispPb_e);
  TGraphErrors *grRadiiC2Sys_pp=new TGraphErrors(20,meanNchpp,yAxispp,xAxis_e,yAxispp_e);
  grRadiiC2Sys_pp->SetMarkerSize(0); grRadiiC2Sys_pp->SetFillStyle(3001); grRadiiC2Sys_pp->SetFillColor(kBlue-10);
  grRadiiC2Sys_pPb->SetMarkerSize(0); grRadiiC2Sys_pPb->SetFillStyle(3001); grRadiiC2Sys_pPb->SetFillColor(kRed-10);
  grRadiiC2Sys_PbPb->SetMarkerSize(0); grRadiiC2Sys_PbPb->SetFillStyle(3001); grRadiiC2Sys_PbPb->SetFillColor(kGray);
  //
  //grRadiiC2Sys_pp->Draw("E2 p");
  //grRadiiC2Sys_pPb->Draw("E2 p");
  grRadiiC2Sys_PbPb->Draw("E2 p");
  grRadiiSys_pp->Draw("E2 p");
  grRadiiSys_pPb->Draw("E2 p");
  grRadiiSys_PbPb->Draw("E2 p");
  RadiiPbPb->Draw("same");
  RadiipPb->Draw("same");
  Radiipp->Draw("same");
  //
  TF1 *PbPbFit=new TF1("PbPbFit","pol1",2.1,3.6); PbPbFit->SetLineColor(1);
  //RadiiPbPb->Fit(PbPbFit,"IMENQ","",2.1,3.6);
  //PbPbFit->Draw("same");
  //cout<<"PbPb: "<<PbPbFit->GetParameter(0)<<", "<<PbPbFit->GetParameter(1)<<endl;
  TF1 *pPbFit=new TF1("pPbFit","pol1",0.8,1.8); pPbFit->SetLineColor(2);
  //RadiipPb->Fit(pPbFit,"IMENQ","",0.8,1.8);
  //pPbFit->Draw("same");
  //cout<<"pPb: "<<pPbFit->GetParameter(0)<<", "<<pPbFit->GetParameter(1)<<endl;
  TF1 *ppFit=new TF1("ppFit","pol1",0,1.5); ppFit->SetLineColor(4);
  //Radiipp->Fit(ppFit,"IMENQ","",0,1.5);
  //ppFit->Draw("same");
  //cout<<"pp: "<<ppFit->GetParameter(0)<<", "<<ppFit->GetParameter(1)<<endl;
  //
  double parsPbPb[3][2]={{-5.34719, 4.30537},{-4.52015, 3.95295},{-4.02046, 3.59879}};// FB7 values of pol1 fit (3 Kt3 bins)
  double parspPb[3][2]={{0.622312, 0.911277},{0.403686, 1.04216},{0.430842, 0.904828}};
  double parspp[3][2]={{0.717018, 0.777676},{0.754759, 0.691637},{0.756552, 0.611244}};
  PbPbFit->FixParameter(0,parsPbPb[KT3Bin][0]); PbPbFit->FixParameter(1,parsPbPb[KT3Bin][1]);
  pPbFit->FixParameter(0,parspPb[KT3Bin][0]); pPbFit->FixParameter(1,parspPb[KT3Bin][1]);
  ppFit->FixParameter(0,parspp[KT3Bin][0]); ppFit->FixParameter(1,parspp[KT3Bin][1]);
  //PbPbFit->Draw("same"); pPbFit->Draw("same"); ppFit->Draw("same");
  //
  RadiiC2PbPb->Draw("same");
  //RadiiC2pPb->Draw("same");
  //RadiiC2pp->Draw("same");
  
  legend4->AddEntry(Radiipp,"pp #sqrt{s}=7 TeV","p");
  legend4->AddEntry(RadiipPb,"p-Pb #sqrt{s_{NN}}=5.02 TeV","p");
  legend4->AddEntry(RadiiPbPb,"Pb-Pb #sqrt{s_{NN}}=2.76 TeV","p");

  TLatex *Specif_Marker_1 = new TLatex(0.57,0.23,"Hollow=R_{inv,2}  Solid=R_{inv,3}");
  Specif_Marker_1->SetNDC();
  Specif_Marker_1->SetTextFont(TextFont);
  Specif_Marker_1->SetTextSize(SizeSpecif*SF1);
  Specif_Marker_1->Draw("same");
  TLatex *Specif_kappas = new TLatex(0.57,0.16,"#kappa_{3}=0.12, #kappa_{4}=0.43");
  Specif_kappas->SetNDC();
  Specif_kappas->SetTextFont(TextFont);
  Specif_kappas->SetTextSize(SizeSpecif*SF1);
  if(IncludeEW) Specif_kappas->Draw("same");
  TLatex *Specif_Gauss = new TLatex(0.57,0.16,"Gaussian Radii");
  Specif_Gauss->SetNDC();
  Specif_Gauss->SetTextFont(TextFont);
  Specif_Gauss->SetTextSize(SizeSpecif*SF1);
  if(!IncludeEW) Specif_Gauss->Draw("same");

  TLatex *Specif_Kinematics_2 = new TLatex(0.57,0.16,"|#eta|<0.8, 0.16<p_{T}<1.0 GeV/c");
  Specif_Kinematics_2->SetNDC();
  Specif_Kinematics_2->SetTextFont(TextFont);
  Specif_Kinematics_2->SetTextSize(SizeSpecif*SF1);
  //Specif_Kinematics_2->Draw("same");

  legend4->SetTextFont(TextFont);
  legend4->SetTextSize(SizeLegend*SF1);
  legend4->Draw("same");
  
  ///////////////////////////////////////////////////////////////////
  pad3->cd(2);
  gPad->SetRightMargin(RightMargin); gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.1); gPad->SetTopMargin(0.02);
  gPad->SetTickx(); gPad->SetTicky();
  //gPad->SetGridx(); gPad->SetGridy();
  TH1D *LambdaPbPb=(TH1D*)Parameters_c3[0][KT3Bin][1]->Clone();
  TH1D *LambdapPb=(TH1D*)Parameters_c3[1][KT3Bin][1]->Clone();
  TH1D *Lambdapp=(TH1D*)Parameters_c3[2][KT3Bin][1]->Clone();
  
  LambdaPbPb->GetXaxis()->SetLabelFont(TextFont); LambdaPbPb->GetYaxis()->SetLabelFont(TextFont); 
  LambdaPbPb->GetXaxis()->SetLabelSize(SizeLabel*SF1); LambdaPbPb->GetYaxis()->SetLabelSize(SizeLabel*SF1);
  LambdaPbPb->GetXaxis()->SetNdivisions(808);
  LambdaPbPb->GetYaxis()->SetTitleFont(TextFont); LambdaPbPb->GetYaxis()->SetTitleSize(SizeTitle*SF1);
  LambdaPbPb->SetMaximum(2.0);
  LambdaPbPb->GetXaxis()->SetTitleOffset(1.3);
  LambdaPbPb->GetYaxis()->SetTitleOffset(1.1);
  
  LambdaPbPb->Draw();
  
  for(int cb=0; cb<20; cb++){
    int binPbPb = LambdaPbPb->GetXaxis()->FindBin(meanNchPbPb[cb]);
    int binpPb = LambdapPb->GetXaxis()->FindBin(meanNchpPb[cb]);
    int binpp = Lambdapp->GetXaxis()->FindBin(meanNchpp[cb]);
    double LambdasysPbPb = sqrt(pow(0.05,2)+pow(0.1,2));// MRC + Qrange
    double LambdasyspPb = sqrt(pow(0.03,2)+pow(0.05,2));// MRC + Qrange
    double Lambdasyspp = sqrt(pow(0.03,2)+pow(0.05,2));// MRC + Qrange
    if(LambdaPbPb->GetBinContent(binPbPb)==0) {yAxisPbPb[cb]=100; yAxisPbPb_e[cb]=100;}
    else {yAxisPbPb[cb]=LambdaPbPb->GetBinContent(binPbPb); yAxisPbPb_e[cb]=LambdasysPbPb;}
    if(LambdapPb->GetBinContent(binpPb)==0) {yAxispPb[cb]=100; yAxispPb_e[cb]=100;}
    else {yAxispPb[cb]=LambdapPb->GetBinContent(binpPb); yAxispPb_e[cb]=LambdasyspPb;}
    if(Lambdapp->GetBinContent(binpp)==0) {yAxispp[cb]=100; yAxispp_e[cb]=100;}
    else {yAxispp[cb]=Lambdapp->GetBinContent(binpp); yAxispp_e[cb]=Lambdasyspp;}
  }
  TGraphErrors *grLambdaSys_PbPb=new TGraphErrors(20,meanNchPbPb,yAxisPbPb,xAxis_e,yAxisPbPb_e);
  TGraphErrors *grLambdaSys_pPb=new TGraphErrors(20,meanNchpPb,yAxispPb,xAxis_e,yAxispPb_e);
  TGraphErrors *grLambdaSys_pp=new TGraphErrors(20,meanNchpp,yAxispp,xAxis_e,yAxispp_e);
  grLambdaSys_pp->SetMarkerSize(0); grLambdaSys_pp->SetFillStyle(1000); grLambdaSys_pp->SetFillColor(kBlue-10);
  grLambdaSys_pPb->SetMarkerSize(0); grLambdaSys_pPb->SetFillStyle(1000); grLambdaSys_pPb->SetFillColor(kRed-10);
  grLambdaSys_PbPb->SetMarkerSize(0); grLambdaSys_PbPb->SetFillStyle(1000); grLambdaSys_PbPb->SetFillColor(kGray);
  // C2 sys
  TH1D *LambdaC2PbPb=(TH1D*)Parameters_C2[0][KT3Bin][1]->Clone();
  TH1D *LambdaC2pPb=(TH1D*)Parameters_C2[1][KT3Bin][1]->Clone();
  TH1D *LambdaC2pp=(TH1D*)Parameters_C2[2][KT3Bin][1]->Clone();
  for(int cb=0; cb<20; cb++){
    int binPbPb = LambdaC2PbPb->GetXaxis()->FindBin(meanNchPbPb[cb]);
    int binpPb = LambdaC2pPb->GetXaxis()->FindBin(meanNchpPb[cb]);
    int binpp = LambdaC2pp->GetXaxis()->FindBin(meanNchpp[cb]);
    double LambdasysPbPb = sqrt(pow(0.03,2)+pow(0.04,2));// MRC + Qrange
    double LambdasyspPb = sqrt(pow(0.02,2)+pow(0.01,2));// MRC + Qrange
    double Lambdasyspp = sqrt(pow(0.02,2)+pow(0.01,2));// MRC + Qrange
    if(LambdaC2PbPb->GetBinContent(binPbPb)==0) {yAxisPbPb[cb]=100; yAxisPbPb_e[cb]=100;}
    else {yAxisPbPb[cb]=LambdaC2PbPb->GetBinContent(binPbPb); yAxisPbPb_e[cb]=LambdasysPbPb;}
    if(LambdaC2pPb->GetBinContent(binpPb)==0) {yAxispPb[cb]=100; yAxispPb_e[cb]=100;}
    else {yAxispPb[cb]=LambdaC2pPb->GetBinContent(binpPb); yAxispPb_e[cb]=LambdasyspPb;}
    if(LambdaC2pp->GetBinContent(binpp)==0) {yAxispp[cb]=100; yAxispp_e[cb]=100;}
    else {yAxispp[cb]=LambdaC2pp->GetBinContent(binpp); yAxispp_e[cb]=Lambdasyspp;}
  }
  TGraphErrors *grLambdaC2Sys_PbPb=new TGraphErrors(20,meanNchPbPb,yAxisPbPb,xAxis_e,yAxisPbPb_e);
  TGraphErrors *grLambdaC2Sys_pPb=new TGraphErrors(20,meanNchpPb,yAxispPb,xAxis_e,yAxispPb_e);
  TGraphErrors *grLambdaC2Sys_pp=new TGraphErrors(20,meanNchpp,yAxispp,xAxis_e,yAxispp_e);
  grLambdaC2Sys_pp->SetMarkerSize(0); grLambdaC2Sys_pp->SetFillStyle(3001); grLambdaC2Sys_pp->SetFillColor(kBlue-10);
  grLambdaC2Sys_pPb->SetMarkerSize(0); grLambdaC2Sys_pPb->SetFillStyle(3001); grLambdaC2Sys_pPb->SetFillColor(kRed-10);
  grLambdaC2Sys_PbPb->SetMarkerSize(0); grLambdaC2Sys_PbPb->SetFillStyle(3001); grLambdaC2Sys_PbPb->SetFillColor(kGray);
  //
  //grLambdaC2Sys_pp->Draw("E2 p");
  //grLambdaC2Sys_pPb->Draw("E2 p");
  grLambdaC2Sys_PbPb->Draw("E2 p");

  grLambdaSys_pp->Draw("E2 p");
  grLambdaSys_pPb->Draw("E2 p");
  grLambdaSys_PbPb->Draw("E2 p");
  
  //
  LambdaPbPb->Draw("same");
  LambdapPb->Draw("same");
  Lambdapp->Draw("same");
  //
  LambdaC2PbPb->Draw("same");
  //LambdaC2pPb->Draw("same");
  //LambdaC2pp->Draw("same");
  
  TLatex *Specif_Marker_2 = new TLatex(0.57,0.9,"Hollow=#lambda  Solid=#lambda_{3}");
  Specif_Marker_2->SetNDC();
  Specif_Marker_2->SetTextFont(TextFont);
  Specif_Marker_2->SetTextSize(SizeSpecif*SF1);
  Specif_Marker_2->Draw("same");
  TLatex *Specif_Kt3;
  TLatex *Specif_kt;
  if(KT3Bin==0) {Specif_Kt3 = new TLatex(0.57, 0.83,"0.16<K_{T,3}<0.3, <k_{T}>=0.3 GeV/c"); Specif_kt = new TLatex(0.57, 0.76,"0.3<k_{T}<0.4, <k_{T}>=0.32 GeV/c");}
  if(KT3Bin==1) {Specif_Kt3 = new TLatex(0.57, 0.83,"0.3<K_{T,3}<1.0, <k_{T}>=0.4 GeV/c"); Specif_kt = new TLatex(0.57, 0.76,"0.4<k_{T}<0.5, <k_{T}>=0.42 GeV/c");}
  //if(KT3Bin==0) {Specif_Kt3 = new TLatex(0.57, 0.83,"0.16<K_{T,3}<0.25 GeV/c"); Specif_kt = new TLatex(0.57, 0.76,"0.2<k_{T}<0.3 GeV/c");}
  //if(KT3Bin==1) {Specif_Kt3 = new TLatex(0.57, 0.83,"0.25<K_{T,3}<0.35 GeV/c"); Specif_kt = new TLatex(0.57, 0.76,"0.3<k_{T}<0.4 GeV/c");}
  //if(KT3Bin==2) {Specif_Kt3 = new TLatex(0.57, 0.83,"0.35<K_{T,3}<1.0 GeV/c"); Specif_kt = new TLatex(0.57, 0.76,"0.4<k_{T}<0.5 GeV/c");}
  Specif_Kt3->SetTextFont(TextFont); Specif_kt->SetTextFont(TextFont);
  Specif_Kt3->SetTextSize(SizeSpecif*SF1); Specif_kt->SetTextSize(SizeSpecif*SF1);
  Specif_Kt3->SetNDC(); Specif_kt->SetNDC();
  Specif_Kt3->Draw("same");
  Specif_kt->Draw("same");

  if(SaveFiles && !IncludeEW) can3->SaveAs("ThreePionFitParametersGauss.eps");
  if(SaveFiles && IncludeEW) can3->SaveAs("ThreePionFitParametersEW.eps");
  
  
  ////////////////////////////////////////////////////
  ////////////////////////////////////////////////////
  // Correlation functions and Monte Carlo

  TCanvas *can4 = new TCanvas("can4", "can4",1700,0,600,600);// 11,53,700,500
  can4->SetHighLightColor(2);
  gStyle->SetOptFit(0111);
  can4->SetFillColor(10);//10
  can4->SetBorderMode(0);
  can4->SetBorderSize(2);
  can4->SetFrameFillColor(0);
  can4->SetFrameBorderMode(0);
  can4->SetFrameBorderMode(0);
  can4->cd();
  TPad *pad4 = new TPad("pad4","pad4",0.0,0.0,1.,1.);
  gPad->SetGridx(0);
  gPad->SetGridy(0);
  gPad->SetTickx();
  gPad->SetTicky();
  pad4->SetTopMargin(0.02);//0.05
  pad4->SetRightMargin(0.01);//3e-2
  pad4->SetBottomMargin(0.07);//0.12
  pad4->Draw();
  pad4->cd();
  TLegend *legend5 = new TLegend(.55,.78, .97,.97,NULL,"brNDC");//.45 or .4 for x1
  legend5->SetBorderSize(0);
  legend5->SetFillColor(0);
  legend5->SetTextFont(TextFont);
  legend5->SetTextSize(SizeLegend*SF2);
  //
  gPad->SetRightMargin(0.03); gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.1); gPad->SetTopMargin(0.02);
  //
  int System_proof=0;
  int ChComb_proof=0;
  int Mb_proof=2;// 0 and 12 for PbPb
  
  C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->SetMinimum(0.9); C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->SetMaximum(3.0);
  if(System_proof==0 && Mb_proof<10) C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->SetMaximum(2.0);
  //
  //c3[2][0][0][KT3Bin][Mb_pp]->GetXaxis()->SetLabelSize(SizeLabel*SF2); c3[2][0][0][KT3Bin][Mb_pp]->GetYaxis()->SetLabelSize(SizeLabel*SF2);
  //c3[2][0][0][KT3Bin][Mb_pp]->GetXaxis()->SetNdivisions(808);
  //c3[2][0][0][KT3Bin][Mb_pp]->GetYaxis()->SetTitleSize(SizeTitle*SF1);
  C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetYaxis()->SetTitle("#font[12]{C}_{3} or #font[12]{#bf{c}}_{3}");
  C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetXaxis()->SetTitleOffset(1.2); C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetYaxis()->SetTitleOffset(1.2);
  C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetXaxis()->SetNdivisions(808);
  double Q3Limit;
  if(System_proof==0) Q3Limit = 0.1 + 0.1*Mb_proof/16.;
  else Q3Limit = 0.3 + 0.2*fabs(Mb_proof-10)/9.;
  C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetXaxis()->SetRangeUser(0,Q3Limit);
  C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->Draw();
  

  C3_Sys[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->Draw("E2 same");
  c3_Sys[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->Draw("E2 same");
  C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->Draw("same");
  c3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->Draw("same");
  if(System_proof!=0) c3[System_proof][1][ChComb_proof][KT3Bin][Mb_proof]->Draw("same");
  if(ChComb_proof==0){
    legend5->AddEntry(C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof],"#font[12]{C}_{3}^{#pm#pm#pm}","p");
    legend5->AddEntry(c3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof],"#font[12]{#bf{c}}_{3}^{#pm#pm#pm}","p");
    if(System_proof==1) legend5->AddEntry(c3[System_proof][1][ChComb_proof][KT3Bin][Mb_proof],"DPMJET #font[12]{#bf{c}}_{3}^{#pm#pm#pm}","p");
    if(System_proof==2) legend5->AddEntry(c3[System_proof][1][ChComb_proof][KT3Bin][Mb_proof],"Pythia Perugia 6 #font[12]{#bf{c}}_{3}^{#pm#pm#pm}","p");
  }else{
    legend5->AddEntry(C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof],"#font[12]{C}_{3}^{#pm#pm#mp}","p");
    legend5->AddEntry(c3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof],"#font[12]{#bf{c}}_{3}^{#pm#pm#mp}","p");
    if(System_proof==1) legend5->AddEntry(c3[System_proof][1][ChComb_proof][KT3Bin][Mb_proof],"DPMJET #font[12]{#bf{c}}_{3}^{#pm#pm#mp}","p");
    if(System_proof==2) legend5->AddEntry(c3[System_proof][1][ChComb_proof][KT3Bin][Mb_proof],"Pythia Perugia 6 #font[12]{#bf{c}}_{3}^{#pm#pm#mp}","p");
  } 
  
  if(ChComb_proof==0) {
    //c3_fit[System_proof][ChComb_proof][KT3Bin][Mb_proof]->Draw("c same");
    if(!IncludeEW) legend5->AddEntry(c3_fit[System_proof][ChComb_proof][KT3Bin][Mb_proof],"Gaussian Fit to #font[12]{#bf{c}}_{3}^{#pm#pm#pm}","l");
    if(IncludeEW) legend5->AddEntry(c3_fit[System_proof][ChComb_proof][KT3Bin][Mb_proof],"Edgeworth Fit to #font[12]{#bf{c}}_{3}^{#pm#pm#pm}","l");
  }
  
  TLatex *CTLabel;
  if(System_proof==0) CTLabel = new TLatex(0.55,0.73,"Pb-Pb #sqrt{s_{NN}}=2.76 TeV");// 0.003,.988
  if(System_proof==1) CTLabel = new TLatex(0.55,0.73,"p-Pb #sqrt{s_{NN}}=5.02 TeV");// 0.003,.988
  if(System_proof==2) CTLabel = new TLatex(0.55,0.73,"pp #sqrt{s}=7 TeV");// 0.003,.988
  CTLabel->SetNDC();
  CTLabel->SetTextFont(TextFont);
  CTLabel->SetTextSize(SizeTitle*SF2);
  CTLabel->Draw("same");

  TString *nameCh=new TString("N_{ch} = ");
  float Nch=1;
  if(System_proof==0) Nch = pow(10,meanNchPbPb[Mb_proof]);
  else if(System_proof==1) Nch = pow(10,meanNchpPb[Mb_proof]);
  else Nch = pow(10,meanNchpp[Mb_proof]);
  *nameCh += int(Nch);
  nameCh->Append(" #pm ");
  float SysPercent = 0.05;
  int SigFig=0;
  if(SysPercent*Nch < 1) {nameCh->Append("0."); SigFig=int(10*SysPercent*Nch + 0.5);}
  else {SigFig=int(SysPercent*Nch + 0.5);}
  
  *nameCh += SigFig;
  TLatex *MLabel = new TLatex(0.55,0.66,nameCh->Data());// (0.1 + 0.1*Mb_proof/16.) * 0.5,2.13
  MLabel->SetNDC();
  MLabel->SetTextFont(TextFont);
  MLabel->SetTextSize(SizeTitle*SF2);
  MLabel->Draw("same");
  

  TLatex *Specif_Kt3_3 = new TLatex(0.55,0.59,"0.16<K_{T,3}<0.3 GeV/c");
  Specif_Kt3_3->SetNDC();
  Specif_Kt3_3->SetTextFont(TextFont);
  Specif_Kt3_3->SetTextSize(SizeTitle*SF2);
  Specif_Kt3_3->Draw("same");
  TLatex *Specif_Kinematics_3 = new TLatex(0.55,0.52,"|#eta|<0.8, 0.16<p_{T}<1.0 GeV/c");
  Specif_Kinematics_3->SetNDC();
  Specif_Kinematics_3->SetTextFont(TextFont);
  Specif_Kinematics_3->SetTextSize(SizeSpecif*SF2);
  //Specif_Kinematics_3->Draw("same");

  //TString *nameCh=new TString("N_{rec}: ");
  //nameCh->Append(labels[Mb_proof]);
  //TLatex *MLabel;
  //if(System_proof==0) MLabel = new TLatex((0.1 + 0.1*Mb_proof/16.) * 0.5,2.0,nameCh->Data());
  //else MLabel = new TLatex( (0.3 + 0.2*fabs(Mb_proof-10)/9.)* 0.5,2.0,nameCh->Data());
  //MLabel->SetTextFont(TextFont);
  //MLabel->SetTextSize(SizeTitle*SF2);
  //MLabel->Draw("same");
  //TLatex *CTLabel;
  //if(System_proof==0) CTLabel = new TLatex((0.1 + 0.1*Mb_proof/16.) * 0.5,2.2,"Pb-Pb #sqrt{s_{NN}}=2.76 TeV");// 0.003,.988
  //if(System_proof==1) CTLabel = new TLatex((0.3 + 0.2*fabs(Mb_proof-10)/9.)* 0.5,2.2,"p-Pb #sqrt{s_{NN}}=5.02 TeV");// 0.003,.988
  //if(System_proof==2) CTLabel = new TLatex((0.3 + 0.2*fabs(Mb_proof-10)/9.)* 0.5,2.2,"pp #sqrt{s}=7 TeV");// 0.003,.988
  //CTLabel->SetTextFont(TextFont);
  //CTLabel->SetTextSize(SizeTitle*SF2);
  //CTLabel->Draw("same");
  
  legend5->Draw("same");

  if(SaveFiles) {
    TString *name = new TString("ThreePionCorrelations");
    if(ChComb_proof==0) name->Append("SC");
    else name->Append("MC");
    if(IncludeEW) name->Append("EW");
    if(System_proof==0) name->Append("_PbPb_M");
    else if(System_proof==1) name->Append("_pPb_M");
    else name->Append("_pp_M");
    *name += Mb_proof;
    name->Append(".eps");

    can4->SaveAs(name->Data());
    
  }
  


  
  ////////////////////////////////////////////////////
  ////////////////////////////////////////////////////
  // 3 system correlation function comparison
  TCanvas *can5 = new TCanvas("can5", "can5",1700,700,600,600);// 11,53,700,500
  can5->SetHighLightColor(2);
  gStyle->SetOptFit(0111);
  can5->SetFillColor(10);//10
  can5->SetBorderMode(0);
  can5->SetBorderSize(2);
  can5->SetFrameFillColor(0);
  can5->SetFrameBorderMode(0);
  can5->SetFrameBorderMode(0);
  can5->cd();
  TPad *pad5 = new TPad("pad5","pad5",0.0,0.0,1.,1.);
  gPad->SetGridx(0);
  gPad->SetGridy(0);
  gPad->SetTickx();
  gPad->SetTicky();
  pad5->SetTopMargin(0.02);//0.05
  pad5->SetRightMargin(0.01);//3e-2
  pad5->SetBottomMargin(0.07);//0.12
  pad5->Draw();
  pad5->cd();
  TLegend *legend6 = new TLegend(.35,.75, .97,.95,NULL,"brNDC");//.45 or .4 for x1
  legend6->SetBorderSize(0);
  legend6->SetFillColor(0);
  legend6->SetTextFont(TextFont);
  legend6->SetTextSize(SizeLegend*SF2);
  //
  gPad->SetRightMargin(0.03); gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.1); gPad->SetTopMargin(0.02);
  //
  int Mbin_SysCompPbPb=15;// 15
  int Mbin_SysComppPb=18;// 14, 16, 18
  int Mbin_SysComppp=18;// 13, 15, 18
  TH1D *c3_PbPb=(TH1D*)c3[0][0][0][KT3Bin][Mbin_SysCompPbPb]->Clone();
  TH1D *c3_pPb=(TH1D*)c3[1][0][0][KT3Bin][Mbin_SysComppPb]->Clone();
  TH1D *c3_pp=(TH1D*)c3[2][0][0][KT3Bin][Mbin_SysComppp]->Clone();

  c3_pp->GetYaxis()->SetTitle("#bf{c}_{3}^{#pm#pm#pm}");
  c3_pp->GetXaxis()->SetRangeUser(0,0.5);
  c3_pp->Draw();
  //
  c3_Sys[2][0][0][KT3Bin][Mbin_SysComppp]->Draw("E2 same");
  c3_Sys[1][0][0][KT3Bin][Mbin_SysComppPb]->Draw("E2 same");
  //c3_Sys[0][0][0][KT3Bin][Mbin_SysCompPbPb]->Draw("E2 same");
  c3_pp->Draw("same");
  c3_pPb->Draw("same");
  //c3_PbPb->Draw("same");
  // SysPercent = 0.07 - 0.05*(Mb_proof/19.);
  //legend6->AddEntry(c3_pp,"N_{ch} = 48#pm1.7 pp #sqrt{s}=7 TeV","p");// Mpp=13
  //legend6->AddEntry(c3_pPb,"N_{ch} = 46#pm1.5 p-Pb #sqrt{s_{NN}}=5.02 TeV","p");// MpPb=14
  //legend6->AddEntry(c3_PbPb,"N_{ch} = 39#pm1.2 Pb-Pb #sqrt{s_{NN}}=2.76 TeV","p");// MPbPb=15
  //legend6->AddEntry(c3_pp,"N_{ch} = 24#pm1.2 pp #sqrt{s}=7 TeV","p");// Mpp=15
  //legend6->AddEntry(c3_pPb,"N_{ch} = 20#pm1.0 p-Pb #sqrt{s_{NN}}=5.02 TeV","p");// MpPb=16
  legend6->AddEntry(c3_pp,"N_{ch} = 8.8#pm0.4 pp #sqrt{s}=7 TeV","p");// Mpp=18
  legend6->AddEntry(c3_pPb,"N_{ch} = 10.2#pm0.5 p-Pb #sqrt{s_{NN}}=5.02 TeV","p");// MpPb=18

  //c3_fit[0][0][KT3Bin][Mbin_SysCompPbPb]->Draw("c same");
  c3_fit[1][0][KT3Bin][Mbin_SysComppPb]->Draw("c same");
  c3_fit[2][0][KT3Bin][Mbin_SysComppp]->Draw("c same");
  legend6->Draw("same");

  
  TLatex *Specif_Kt3_4 = new TLatex(0.55,0.59,"0.16<K_{T,3}<0.3 GeV/c");
  Specif_Kt3_4->SetNDC();
  Specif_Kt3_4->SetTextFont(TextFont);
  Specif_Kt3_4->SetTextSize(SizeTitle*SF2);
  Specif_Kt3_4->Draw("same");
  TLatex *Specif_Kinematics_4 = new TLatex(0.55,0.52,"|#eta|<0.8, 0.16<p_{T}<1.0 GeV/c");
  Specif_Kinematics_4->SetNDC();
  Specif_Kinematics_4->SetTextFont(TextFont);
  Specif_Kinematics_4->SetTextSize(SizeSpecif*SF2);
  //Specif_Kinematics_4->Draw("same");

  //
  if(SaveFiles) {
    TString *name = new TString("c3SystemComp");
    if(IncludeEW) name->Append("EW");
    name->Append("_MpPb");
    *name += Mbin_SysComppPb;
    name->Append(".eps");
    //
    can5->SaveAs(name->Data());
  }
  

  /*
  ////////////////////////////////////////////////////
  ////////////////////////////////////////////////////
  // C2 correlation functions
  TCanvas *can6 = new TCanvas("can6", "can6",500,700,600,600);// 11,53,700,500
  can6->SetHighLightColor(2);
  gStyle->SetOptFit(0111);
  can6->SetFillColor(10);//10
  can6->SetBorderMode(0);
  can6->SetBorderSize(2);
  can6->SetFrameFillColor(0);
  can6->SetFrameBorderMode(0);
  can6->SetFrameBorderMode(0);
  can6->cd();
  TPad *pad6 = new TPad("pad6","pad6",0.0,0.0,1.,1.);
  gPad->SetGridx(0);
  gPad->SetGridy(0);
  gPad->SetTickx();
  gPad->SetTicky();
  pad6->SetTopMargin(0.02);//0.05
  pad6->SetRightMargin(0.01);//3e-2
  pad6->SetBottomMargin(0.07);//0.12
  pad6->Draw();
  pad6->cd();
  TLegend *legend7 = new TLegend(.35,.75, .97,.95,NULL,"brNDC");//.45 or .4 for x1
  legend7->SetBorderSize(0);
  legend7->SetFillColor(0);
  legend7->SetTextFont(TextFont);
  legend7->SetTextSize(SizeLegend*SF2);
  //
  gPad->SetRightMargin(0.03); gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.1); gPad->SetTopMargin(0.02);
  
  int Mbin_C2=19;
  //C2[2][Mbin_C2]->GetYaxis()->SetTitleOffset(1.2);
  C2[2][Mbin_C2]->Draw();
  //
  C2_Sys[2][Mbin_C2]->Draw("E2 same");
  C2_Sys[1][Mbin_C2]->Draw("E2 same");
  C2[2][Mbin_C2]->Draw("same");
  C2[1][Mbin_C2]->Draw("same");
  //C2[0][15]->Draw("same");
  Fit_h_C2[2][Mbin_C2]->Draw("C same");
  Fit_h_C2[1][Mbin_C2]->Draw("C same");
  
  //legend7->AddEntry(C2[2][Mbin_C2],"N_{ch} = 8.8 pp #sqrt{s}=7 TeV","p");
  //legend7->AddEntry(C2[1][Mbin_C2],"N_{ch} = 10.2 p-Pb #sqrt{s_{NN}}=5.02 TeV","p");
  legend7->AddEntry(C2[2][Mbin_C2],"N_{ch} = 4.6 pp #sqrt{s}=7 TeV","p");
  legend7->AddEntry(C2[1][Mbin_C2],"N_{ch} = 5.4 p-Pb #sqrt{s_{NN}}=5.02 TeV","p");
  //
  legend7->Draw("same");
  TLatex *Specif_Kt = new TLatex(0.22,1.35,"0.16<K_{t}<1.0 GeV/c");
  Specif_Kt->SetTextFont(TextFont);
  Specif_Kt->SetTextSize(SizeSpecif*SF2);
  Specif_Kt->Draw("same");
  TLatex *Specif_Kinematics_5 = new TLatex(0.172,1.3,"|#eta|<0.8, 0.16<p_{T}<1.0 GeV/c");
  Specif_Kinematics_5->SetTextFont(TextFont);
  Specif_Kinematics_5->SetTextSize(SizeSpecif*SF2);
  Specif_Kinematics_5->Draw("same");
  */















  /*
  TLegend *legend6 = new TLegend(.2,.70, .4,.95,NULL,"brNDC");//.45 or .4 for x1
  legend6->SetBorderSize(0);
  legend6->SetFillColor(0);
  legend6->SetTextFont(TextFont);
  legend6->SetTextSize(SizeLegend*SF2);
  Parameters_c3[0][4]->GetXaxis()->SetTitleOffset(1.2); Parameters_c3[0][3]->GetYaxis()->SetTitleOffset(1.4);
  Parameters_c3[0][4]->Draw();
  Parameters_c3[1][4]->Draw("same");
  Parameters_c3[2][4]->Draw("same");
  legend6->AddEntry(Parameters_c3[2][2],"pp #sqrt{s}=7 TeV","p");
  legend6->AddEntry(Parameters_c3[1][2],"p-Pb #sqrt{s_{NN}}=5.02 TeV","p");
  legend6->AddEntry(Parameters_c3[0][2],"Pb-Pb #sqrt{s_{NN}}=2.76 TeV","p");
  legend6->Draw("same");
  double sumk3=0, Enk3=0, sumk4=0, Enk4=0;
  for(int ct=0; ct<3; ct++){ 
    for(int cb=0; cb<20; cb++){
      if(Parameters_c3[ct][3]->GetBinError(cb+1) >0) {
	sumk3 += Parameters_c3[ct][3]->GetBinContent(cb+1)/Parameters_c3[ct][3]->GetBinError(cb+1); 
	Enk3 += 1/Parameters_c3[ct][3]->GetBinError(cb+1);}
      if(Parameters_c3[ct][4]->GetBinError(cb+1) >0) {
	sumk4 += Parameters_c3[ct][4]->GetBinContent(cb+1)/Parameters_c3[ct][4]->GetBinError(cb+1); 
	Enk4 += 1/Parameters_c3[ct][4]->GetBinError(cb+1);}
    }
  }
  cout<<"Avg k3 = "<<sumk3/Enk3<<endl;
  cout<<"Avg k4 = "<<sumk4/Enk4<<endl;
  */

 
  //DrawALICELogo(kTRUE, 0.82,0.24,1.12,0.44);// C3(+++)
  //DrawALICELogo(kTRUE, 0.32,0.47,0.65,0.7);// C2(+-)
  //DrawALICELogo(kTRUE, 0.3,0.18,0.63,0.36);// r3
  //can1->SaveAs("test.eps");

}
//____________________________________________________________________________________________________
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

