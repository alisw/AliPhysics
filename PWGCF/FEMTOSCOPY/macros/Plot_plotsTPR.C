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
const int KT3Bin=1;// Kt3 bin. 0-1
int FitType=1;// 0 (Gaussian), 1 (Edgeworth)
bool AddedCC=kTRUE;// Charge Conjugate already added?
//
int MbinMaxPbPb=15;// 15
int MbinMinpPb=11;// 13
int MbinMinpp=14;// 14 
int MbinMinPbPb=0;// 0
int MbinMaxpPb=18;// 18
//
//
//
const int MaxKT3Bins=2;
int TextFont=42;// 63, or 42
float SizeLabel=0.03;// 20(63 font), 0.08(42 font)
float SizeLegend=0.03;// 
float SizeTitle=0.03;// 
float SizeSpecif=0.03;// 




double RightMargin=0.002;// 0.002
//
double Chi2_C2global;
double NFitPoints_C2global;
TH1D *C2_ss[20];
TH1D *C2_os[20];



void DrawALICELogo(Bool_t, Float_t, Float_t, Float_t, Float_t);
TCanvas *make_canvas(const Char_t*,const Char_t*,Int_t,Int_t,Int_t,Int_t,Int_t);

void Plot_plotsTPR(){

  gStyle->SetOptStat(0);
  gStyle->SetOptDate(0);
  gStyle->SetOptFit(0111);

  ////////////////////////////////////
  // Get Nrec to Nch mapping
  double meanNchPbPb[20]={0};
  //double meanNchPbPb_e[20]={0};
  double meanNchpPb[20]={0};
  //double meanNchpPb_e[20]={0};
  double meanNchpp[20]={0};
  //double meanNchpp_e[20]={0};
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
  TF1 *Fit_C2[3][2][3][20];// CollType, Gauss/EW, EDbin, cb
  //TH1D *Fit_h_C2[3][3][20];// CollType, EDbin, cb
  TH1D *Parameters_C2[3][2][3][5];// CollType, Gauss/EW, EDbin, Parameter#
  TH1D *C2[3][3][20];// CollType, EDbin, cb
  TH1D *C2_Sys[3][3][20];// CollType, EDbin, cb
  TH1D *C3[3][2][2][3][20];// CollType, Real/MonteCarlo, SC/MC, EDbin, cb
  TH1D *c3[3][2][2][3][20];// CollType, Real/MonteCarlo, SC/MC, EDbin, cb
  TH1D *C3_Sys[3][2][2][3][20];// CollType, Real/MonteCarlo, SC/MC, EDbin, cb
  TH1D *c3_Sys[3][2][2][3][20];// CollType, Real/MonteCarlo, SC/MC, EDbin, cb
  TF1 *c3_fit[3][2][3][20];// CollType, Gauss/EW, EDbin, cb
  TMinuit *Min[3][2][3][20];// CollType, +/-, EDbin, cb
  TH1D *Parameters_c3[3][2][3][5];// CollType, Gaussian/EW, EDbin, Parameter#
  //char *labels[20]={"1050-2000","850-1050","700-850","600-700","500-600","400-500","320-400","260-320","200-260","150-200","100-150","70-100","50-70","40-50","30-40","20-30","15-20","10-15","5-10","0-5"};
  for(int ct=0; ct<3; ct++){
    for(int ft=0; ft<2; ft++){// Gaussian or EW
      for(int kt3=0; kt3<3; kt3++){
	for(int par=0; par<5; par++){
	  TString *name_C2=new TString("Parameters_C2_");
	  *name_C2 += ct;
	  *name_C2 += ft;
	  *name_C2 += kt3;
	  *name_C2 += par;
	  TString *name_c3=new TString("Parameters_c3_");
	  *name_c3 += ct;
	  *name_c3 += ft;
	  *name_c3 += kt3;
	  *name_c3 += par;
	  //Parameters_c3[ct][par] = new TH1D(name->Data(),"",20,0.5,20.5);
	  //for(int cb=0; cb<20; cb++) Parameters_c3[ct][par]->GetXaxis()->SetBinLabel(cb+1,labels[cb]);
	  Parameters_C2[ct][ft][kt3][par] = new TH1D(name_C2->Data(),"",2900,0.6,3.4);// 2900,0.6,3.4
	  Parameters_c3[ct][ft][kt3][par] = new TH1D(name_c3->Data(),"",2900,0.6,3.4);
	  
	}
      }
    }
  }
  
  TH1D *RadiiC2pp_Published = new TH1D("RadiiC2pp_Published","",2900,0.6,3.4);
  //
  double N_1 = 0, N_1_e=0;
  double lambda_1 = 0, lambda_1_e=0;
  double radius_1 = 0, radius_1_e=0;
  double EW1_1 = 0, EW1_1_e=0;
  double EW2_1 = 0, EW2_1_e=0;
  
  
  // Start File access
  for(int ct=0; ct<3; ct++){
    for(int dt=0; dt<2; dt++){// data type (Real or Monte Carlo)
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
	      //if(IncludeEW && dt==0 && ChComb==0) name3->Append("EW");
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
		    //Fit_h_C2[ct][KT3][cb]=(TH1D*)file->Get("fitC2ss_h");
		    Fit_C2[ct][0][KT3][cb]=(TF1*)file->Get("fitC2ss_Gauss");
		    Fit_C2[ct][1][KT3][cb]=(TF1*)file->Get("fitC2ss_EW");
		    C2_Sys[ct][KT3][cb]=(TH1D*)C2[ct][KT3][cb]->Clone();
		    if(ct==0){
		      C2[ct][KT3][cb]->SetMarkerColor(1); C2[ct][KT3][cb]->SetLineColor(1);
		      //Fit_h_C2[ct][KT3][cb]->SetLineColor(1);
		      C2_Sys[ct][KT3][cb]->SetFillColor(kGray);
		    }else if(ct==1){
		      C2[ct][KT3][cb]->SetMarkerColor(2); C2[ct][KT3][cb]->SetLineColor(2);
		      //Fit_h_C2[ct][KT3][cb]->SetLineColor(2);
		      C2_Sys[ct][KT3][cb]->SetFillColor(kRed-10);
		    }else{
		      C2[ct][KT3][cb]->SetMarkerColor(4); C2[ct][KT3][cb]->SetLineColor(4);
		      //Fit_h_C2[ct][KT3][cb]->SetLineColor(4);
		      C2_Sys[ct][KT3][cb]->SetFillColor(kBlue-10);
		    }
		    C2_Sys[ct][KT3][cb]->SetMarkerSize(0); C2_Sys[ct][KT3][cb]->SetFillStyle(1000);
		    // C2 Systematics
		    for(int bin=1; bin<=C2_Sys[ct][KT3][cb]->GetNbinsX(); bin++){
		      C2_Sys[ct][KT3][cb]->SetBinError(bin, 0.01*C2_Sys[ct][KT3][cb]->GetBinContent(bin));
		    }
		    C2[ct][KT3][cb]->SetDirectory(0); //Fit_h_C2[ct][KT3][cb]->SetDirectory(0);
		    C2_Sys[ct][KT3][cb]->SetDirectory(0);
		    C2[ct][KT3][cb]->SetMarkerStyle(24+ct);
		    C2[ct][KT3][cb]->GetXaxis()->SetRangeUser(0,0.33);
		    C2[ct][KT3][cb]->GetXaxis()->SetLabelFont(TextFont); C2[ct][KT3][cb]->GetYaxis()->SetLabelFont(TextFont); 
		    C2[ct][KT3][cb]->GetXaxis()->SetTitleOffset(1.2); C2[ct][KT3][cb]->GetYaxis()->SetTitleOffset(1.2);
		    C2[ct][KT3][cb]->GetXaxis()->SetTitleFont(TextFont); C2[ct][KT3][cb]->GetYaxis()->SetTitleFont(TextFont);
		    C2[ct][KT3][cb]->GetXaxis()->SetTitleSize(SizeTitle); //C2[ct][KT3][cb]->GetYaxis()->SetTitleFont(SizeTitle*SF2);
		    //C2[ct][KT3][cb]->GetYaxis()->SetTitle("C_{2}^{#pm#pm}");
		    //
		    c3_fit[ct][0][KT3][cb]=(TF1*)file->Get("c3Fit1D_Gauss");
		    //cout<<c3_fit[ct][0][ChComb][KT3][cb]->GetName()<<endl;
		    c3_fit[ct][1][KT3][cb]=(TF1*)file->Get("c3Fit1D_EW");
		    //cout<<c3_fit[ct][1][ChComb][KT3][cb]->GetName()<<endl;
		    c3_fit[ct][0][KT3][cb]->SetLineStyle(2);
		    c3_fit[ct][1][KT3][cb]->SetLineStyle(1);
		    //c3_fit[ct][KT3][cb]->SetDirectory(0);
		    if(ct==0) {c3_fit[ct][0][KT3][cb]->SetLineColor(1); c3_fit[ct][1][KT3][cb]->SetLineColor(1);}
		    if(ct==1) {c3_fit[ct][0][KT3][cb]->SetLineColor(2); c3_fit[ct][1][KT3][cb]->SetLineColor(2);}
		    if(ct==2) {c3_fit[ct][0][KT3][cb]->SetLineColor(4); c3_fit[ct][1][KT3][cb]->SetLineColor(4);}
		    
		    
		  }// ChComb==0
			  
		}
	      }else{
		//TH1D *tempC3=(TH1D*)file->Get("C3");
		//TH1D *tempc3=(TH1D*)file->Get("c3");
		//TH1D *tempc3_fit=(TH1D*)file->Get("dentimesFit_c3");
		//if(dt==0) c3_fit[ct][ChComb][KT3][cb]->Add(tempc3_fit);
		if(ChComb==0){
		  N_1 = 0; N_1_e=0;
		  lambda_1 = 0; lambda_1_e=0;
		  radius_1 = 0; radius_1_e=0;
		  EW1_1 = 0; EW1_1_e=0;
		  EW2_1 = 0; EW2_1_e=0;
		  //
		  if(!AddedCC){
		    cout<<"Not Supported!!"<<endl;
		  }
		}
		//
		if(!AddedCC){
		  cout<<"Not Supported!!"<<endl;
		}
	      }// ch==1
	   
	      
	      for(int bin=1; bin<10; bin++){// Remove large error bins
		if(C3[ct][dt][ChComb][KT3][cb]->GetBinError(bin) > 0.33*C3[ct][dt][ChComb][KT3][cb]->GetBinContent(bin)){
		  C3[ct][dt][ChComb][KT3][cb]->SetBinContent(bin,10); C3[ct][dt][ChComb][KT3][cb]->SetBinError(bin,10);
		}
		if(c3[ct][dt][ChComb][KT3][cb]->GetBinError(bin) > 0.33*c3[ct][dt][ChComb][KT3][cb]->GetBinContent(bin)){
		  c3[ct][dt][ChComb][KT3][cb]->SetBinContent(bin,10); c3[ct][dt][ChComb][KT3][cb]->SetBinError(bin,10);
		}
		//cout<<C3[ct][dt][ChComb][KT3][cb]->GetBinError(bin)<<"  "<<C3[ct][dt][ChComb][KT3][cb]->GetBinContent(bin)<<endl;
	      }
	     
	      if(AddedCC && dt==0){
		//Min[ct][0][KT3][cb]->GetParameter(0,N_1,N_1_e);
		//Min[ct][0][KT3][cb]->GetParameter(1,lambda_1,lambda_1_e);
		//Min[ct][0][KT3][cb]->GetParameter(2,radius_1,radius_1_e);
		//Min[ct][0][KT3][cb]->GetParameter(3,EW1_1,EW1_1_e);
		//Min[ct][0][KT3][cb]->GetParameter(4,EW2_1,EW2_1_e);
		//
		if(ChComb==0){
		  double logNch=0;
		  if(ct==0) logNch=meanNchPbPb[cb];
		  else if(ct==1) logNch=meanNchpPb[cb];
		  else logNch=meanNchpp[cb];
		  int logNchBin = Parameters_c3[ct][0][KT3][0]->GetXaxis()->FindBin(logNch);
		  for(int ft=0; ft<2; ft++){// Gaussian or EW
		    N_1 = c3_fit[ct][ft][KT3][cb]->GetParameter(0); N_1_e = c3_fit[ct][ft][KT3][cb]->GetParError(0);
		    lambda_1 = c3_fit[ct][ft][KT3][cb]->GetParameter(1); lambda_1_e = c3_fit[ct][ft][KT3][cb]->GetParError(1);
		    radius_1 = c3_fit[ct][ft][KT3][cb]->GetParameter(2); radius_1_e = c3_fit[ct][ft][KT3][cb]->GetParError(2);
		    EW1_1 = c3_fit[ct][ft][KT3][cb]->GetParameter(3); EW1_1_e = c3_fit[ct][ft][KT3][cb]->GetParError(3);
		    EW2_1 = c3_fit[ct][ft][KT3][cb]->GetParameter(4); EW2_1_e = c3_fit[ct][ft][KT3][cb]->GetParError(4);
		    //
		    Parameters_c3[ct][ft][KT3][0]->SetBinContent(logNchBin, N_1); Parameters_c3[ct][ft][KT3][0]->SetBinError(logNchBin, N_1_e);
		    Parameters_c3[ct][ft][KT3][1]->SetBinContent(logNchBin, lambda_1); Parameters_c3[ct][ft][KT3][1]->SetBinError(logNchBin, lambda_1_e);
		    Parameters_c3[ct][ft][KT3][2]->SetBinContent(logNchBin, radius_1); Parameters_c3[ct][ft][KT3][2]->SetBinError(logNchBin, radius_1_e);
		    Parameters_c3[ct][ft][KT3][3]->SetBinContent(logNchBin, EW1_1); Parameters_c3[ct][ft][KT3][3]->SetBinError(logNchBin, EW1_1_e);
		    Parameters_c3[ct][ft][KT3][4]->SetBinContent(logNchBin, EW2_1); Parameters_c3[ct][ft][KT3][4]->SetBinError(logNchBin, EW2_1_e);
		    // remove unstable c3 Fit points
		    bool badbin=kFALSE;
		    if(ct==0 && cb>12) badbin=kTRUE; 
		    if(ct==1 && cb<12) badbin=kTRUE; 
		    if(ct==2 && cb<14) badbin=kTRUE;
		    if(badbin){
		      Parameters_c3[ct][ft][KT3][0]->SetBinContent(logNchBin, 100); Parameters_c3[ct][ft][KT3][0]->SetBinError(logNchBin, 100);
		      Parameters_c3[ct][ft][KT3][1]->SetBinContent(logNchBin, 100); Parameters_c3[ct][ft][KT3][1]->SetBinError(logNchBin, 100);
		      Parameters_c3[ct][ft][KT3][2]->SetBinContent(logNchBin, 100); Parameters_c3[ct][ft][KT3][2]->SetBinError(logNchBin, 100);
		      Parameters_c3[ct][ft][KT3][3]->SetBinContent(logNchBin, 100); Parameters_c3[ct][ft][KT3][3]->SetBinError(logNchBin, 100);
		      Parameters_c3[ct][ft][KT3][4]->SetBinContent(logNchBin, 100); Parameters_c3[ct][ft][KT3][4]->SetBinError(logNchBin, 100);
		    }
		    //
		    Parameters_C2[ct][ft][KT3][0]->SetBinContent(logNchBin, Fit_C2[ct][ft][KT3][cb]->GetParameter(0));// N
		    Parameters_C2[ct][ft][KT3][0]->SetBinError(logNchBin, Fit_C2[ct][ft][KT3][cb]->GetParError(0));
		    Parameters_C2[ct][ft][KT3][1]->SetBinContent(logNchBin, Fit_C2[ct][ft][KT3][cb]->GetParameter(1));// lambda
		    Parameters_C2[ct][ft][KT3][1]->SetBinError(logNchBin, Fit_C2[ct][ft][KT3][cb]->GetParError(1));
		    Parameters_C2[ct][ft][KT3][2]->SetBinContent(logNchBin, Fit_C2[ct][ft][KT3][cb]->GetParameter(3));// R
		    Parameters_C2[ct][ft][KT3][2]->SetBinError(logNchBin, Fit_C2[ct][ft][KT3][cb]->GetParError(3));
		    Parameters_C2[ct][ft][KT3][3]->SetBinContent(logNchBin, Fit_C2[ct][ft][KT3][cb]->GetParameter(5));// kappa3
		    Parameters_C2[ct][ft][KT3][3]->SetBinError(logNchBin, Fit_C2[ct][ft][KT3][cb]->GetParError(5));
		    Parameters_C2[ct][ft][KT3][4]->SetBinContent(logNchBin, Fit_C2[ct][ft][KT3][cb]->GetParameter(6));// kappa4
		    Parameters_C2[ct][ft][KT3][4]->SetBinError(logNchBin, Fit_C2[ct][ft][KT3][cb]->GetParError(6));
		  }// ft
		}// ChComb==0
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
	for(int ft=0; ft<2; ft++){// Gaussian or EW
	  for(int KT3=0; KT3<3; KT3++){
	    for(int par=0; par<5; par++){
	      if(ct<2){
		Parameters_C2[ct][ft][KT3][par]->SetMarkerStyle(24+ct);
		Parameters_c3[ct][ft][KT3][par]->SetMarkerStyle(20+ct);
	      }else{
		Parameters_C2[ct][ft][KT3][par]->SetMarkerStyle(28);
		Parameters_c3[ct][ft][KT3][par]->SetMarkerStyle(34);
		RadiiC2pp_Published->SetMarkerStyle(28);
	      }
	      
	      if(ct==0){
		Parameters_C2[ct][ft][KT3][par]->SetMarkerColor(1); 
		Parameters_C2[ct][ft][KT3][par]->SetLineColor(1);
		Parameters_c3[ct][ft][KT3][par]->SetMarkerColor(1); 
		Parameters_c3[ct][ft][KT3][par]->SetLineColor(1);
	      }else if(ct==1){
		Parameters_C2[ct][ft][KT3][par]->SetMarkerColor(2); 
		Parameters_C2[ct][ft][KT3][par]->SetLineColor(2);
		Parameters_c3[ct][ft][KT3][par]->SetMarkerColor(2); 
		Parameters_c3[ct][ft][KT3][par]->SetLineColor(2);
	      }else {
		Parameters_C2[ct][ft][KT3][par]->SetMarkerColor(4); 
		Parameters_C2[ct][ft][KT3][par]->SetLineColor(4);
		Parameters_c3[ct][ft][KT3][par]->SetMarkerColor(4); 
		Parameters_c3[ct][ft][KT3][par]->SetLineColor(4);
		RadiiC2pp_Published->SetMarkerColor(4);
		RadiiC2pp_Published->SetLineColor(4);
	      }
	      
	      if(par==0) Parameters_c3[ct][ft][KT3][par]->GetYaxis()->SetTitle("N");
	      //if(par==1) Parameters_c3[ct][ft][KT3][par]->GetYaxis()->SetTitle("#lambda_{3}");
	      if(par==1) Parameters_c3[ct][ft][KT3][par]->GetYaxis()->SetTitle("#lambda or #lambda_{3}");
	      //if(par==2) Parameters_c3[ct][ft][KT3][par]->GetYaxis()->SetTitle("R_{inv,3} (fm)");
	      if(par==2) {
		if(FitType==0) Parameters_c3[ct][ft][KT3][par]->GetYaxis()->SetTitle("R_{inv,2} or R_{inv,3} (fm)");
		else Parameters_c3[ct][ft][KT3][par]->GetYaxis()->SetTitle("R^{Ew}_{inv,2} or R^{Ew}_{inv,3} (fm)");
	      }		
	      if(par==3) Parameters_c3[ct][ft][KT3][par]->GetYaxis()->SetTitle("#kappa_{3}");
	      if(par==4) Parameters_c3[ct][ft][KT3][par]->GetYaxis()->SetTitle("#kappa_{4}");
	      Parameters_c3[ct][ft][KT3][par]->GetXaxis()->SetTitle("log_{10}(N_{ch})");
	    
	      
	    }// par
	  }// KT3
	}// ft
      }// dt==0 
      
    }// dt
  }// ct
  
  cout<<"Done Getting Histograms"<<endl;
  
  TF1 *Unity = new TF1("Unity","1",0,100);
  Unity->SetLineStyle(2);
  Unity->SetLineColor(1);
  
 
  //return;
  
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
  legend3->SetTextSize(SizeLegend);

  gPad->SetRightMargin(0.03); gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.1); gPad->SetTopMargin(0.02);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  int Mb_pp=18;
  int Mb_pPb=13;
  int Mb_PbPb=0;
  c3[2][0][0][0][Mb_pp]->GetXaxis()->SetLabelSize(SizeLabel); c3[2][0][0][0][Mb_pp]->GetYaxis()->SetLabelSize(SizeLabel);
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
  Specif_Kt3_1->SetTextSize(SizeSpecif);
  Specif_Kt3_1->Draw("same");
  TLatex *Specif_Kinematics_1 = new TLatex(0.5,0.52,"|#eta|<0.8, 0.16<#font[12]{p}_{T}<1.0 GeV/#font[12]{c}");
  Specif_Kinematics_1->SetNDC();
  Specif_Kinematics_1->SetTextFont(TextFont);
  Specif_Kinematics_1->SetTextSize(SizeSpecif);
  Specif_Kinematics_1->Draw("same");
  TLatex *Species_1 = new TLatex(0.5,0.45,"#pi^{+}#pi^{+}#pi^{+} & #pi^{-}#pi^{-}#pi^{-} combined");
  Species_1->SetNDC();
  Species_1->SetTextFont(TextFont);
  Species_1->SetTextSize(SizeSpecif);
  Species_1->Draw("same");

  legend3->Draw("same");
  DrawALICELogo(kTRUE, .75, .25, .9, .4);

  if(SaveFiles) can2->SaveAs("ThreePionCorrelation_Evolution.eps");
  */
  ////////////////////////////////////////////////////
  ////////////////////////////////////////////////////
  // Radii and lambda plot
  float RadiiC2ppPubPoints[2][8]={{0.88,1.09,1.23,1.28,1.34,1.37,1.42,1.48},{0.86,1.06,1.16,1.23,1.23,1.28,1.32,1.38}};
  float RadiiC2ppPubPoints_e[2][8]={{0.058,0.064,0.071,0.071,0.078,0.078,0.086,0.098},{0.12,0.12,0.13,0.13,0.13,0.13,0.13,0.13}};
  float MeanPubNch[8]={6.98,14.95,19.68,24.68,29.38,33.95,38.46,42.66};
  float SF=1.3;
  TCanvas *can3 = new TCanvas("can3", "can3",1000,0,600,900);// 11,53,700,500
  can3->SetHighLightColor(2);
  gStyle->SetOptFit(0111);
  can3->SetFillColor(10);//10
  can3->SetBorderMode(0);
  can3->SetBorderSize(2);
  can3->SetFrameFillColor(0);
  can3->SetFrameBorderMode(0);
  can3->SetFrameBorderMode(0);
  can3->cd();
  float PadLeftMargin=0.06, PadBottomMargin=0.;
  TPad *pad3 = new TPad("pad3","pad3",PadLeftMargin,PadBottomMargin,1.,1.);
  gPad->SetTickx();
  gPad->SetTicky();
  pad3->SetTopMargin(0.0);//0.05
  pad3->SetRightMargin(0.0);//3e-2
  pad3->SetBottomMargin(0.0);//0.12
  pad3->SetLeftMargin(0.0);
  pad3->Divide(1,2,0,0);
  pad3->Draw();
  pad3->cd();
  TLegend *legend4 = new TLegend(.1,.65, .5,.95,NULL,"brNDC");//.45 or .4 for x1
  legend4->SetBorderSize(0);
  legend4->SetFillColor(0);
  legend4->SetTextFont(TextFont);
  legend4->SetTextSize(SizeLegend);
  //
  pad3->cd(1);
  gPad->SetLeftMargin(0.05);
  gPad->SetRightMargin(0.01);
  gPad->SetTopMargin(0.02);
  gPad->SetBottomMargin(0.1);
  //gPad->SetRightMargin(RightMargin); gPad->SetLeftMargin(0.12);
  //gPad->SetBottomMargin(0.1);
  gPad->SetTickx(); gPad->SetTicky();
  //gPad->SetGridx(); gPad->SetGridy();
  TH1D *RadiiPbPb=(TH1D*)Parameters_c3[0][FitType][KT3Bin][2]->Clone();
  TH1D *RadiipPb=(TH1D*)Parameters_c3[1][FitType][KT3Bin][2]->Clone();
  TH1D *Radiipp=(TH1D*)Parameters_c3[2][FitType][KT3Bin][2]->Clone();
  RadiiPbPb->GetXaxis()->SetLabelFont(TextFont); RadiiPbPb->GetYaxis()->SetLabelFont(TextFont); 
  RadiiPbPb->GetXaxis()->SetLabelSize(SizeLabel*SF); RadiiPbPb->GetYaxis()->SetLabelSize(SizeLabel*SF);
  RadiiPbPb->GetXaxis()->SetNdivisions(808);
  RadiiPbPb->GetXaxis()->SetTitleOffset(1.1);//1.3
  RadiiPbPb->GetYaxis()->SetTitleOffset(100);//1.1
  RadiiPbPb->GetXaxis()->SetTitleFont(TextFont); RadiiPbPb->GetXaxis()->SetTitleSize(SizeTitle*SF);
  RadiiPbPb->GetYaxis()->SetTitleFont(TextFont); RadiiPbPb->GetYaxis()->SetTitleSize(SizeTitle*SF);
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
  TH1D *RadiiC2PbPb=(TH1D*)Parameters_C2[0][FitType][KT3Bin][2]->Clone();
  TH1D *RadiiC2pPb=(TH1D*)Parameters_C2[1][FitType][KT3Bin][2]->Clone();
  TH1D *RadiiC2pp=(TH1D*)Parameters_C2[2][FitType][KT3Bin][2]->Clone();
  
  for(int mbin=0; mbin<8; mbin++){
    int bin = RadiiC2pp_Published->GetXaxis()->FindBin(log10(MeanPubNch[mbin])); 
    RadiiC2pp_Published->SetBinContent(bin, RadiiC2ppPubPoints[KT3Bin][mbin]); 
    RadiiC2pp_Published->SetBinError(bin, RadiiC2ppPubPoints_e[KT3Bin][mbin]);
  }  

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
  if(FitType==1){
    grRadiiC2Sys_pp->Draw("E2 p");
    grRadiiC2Sys_pPb->Draw("E2 p");
  }  
  grRadiiC2Sys_PbPb->Draw("E2 p");
  grRadiiSys_pp->Draw("E2 p");
  grRadiiSys_pPb->Draw("E2 p");
  grRadiiSys_PbPb->Draw("E2 p");
  RadiiPbPb->Draw("same");
  RadiipPb->Draw("same");
  Radiipp->Draw("same");
  if(FitType==0) RadiiC2pp_Published->Draw("same");
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
  if(FitType==1) {
    RadiiC2pPb->Draw("same");
    RadiiC2pp->Draw("same");
  }

  legend4->AddEntry(Radiipp,"pp #sqrt{s}=7 TeV","p");
  legend4->AddEntry(RadiipPb,"p-Pb #sqrt{s_{NN}}=5.02 TeV","p");
  legend4->AddEntry(RadiiPbPb,"Pb-Pb #sqrt{s_{NN}}=2.76 TeV","p");

  TLatex *Specif_Marker_1 = new TLatex(0.55,0.23,"Hollow=R_{inv,2}  Solid=R_{inv,3}");
  Specif_Marker_1->SetNDC();
  Specif_Marker_1->SetTextFont(TextFont);
  Specif_Marker_1->SetTextSize(SizeSpecif*SF);
  Specif_Marker_1->Draw("same");
  TLatex *Specif_kappas = new TLatex(0.55,0.16,"#kappa_{3}=0.16, #kappa_{4}=0.4");
  Specif_kappas->SetNDC();
  Specif_kappas->SetTextFont(TextFont);
  Specif_kappas->SetTextSize(SizeSpecif*SF);
  if(FitType==1) Specif_kappas->Draw("same");
  TLatex *Specif_Gauss = new TLatex(0.55,0.16,"Gaussian Radii");
  Specif_Gauss->SetNDC();
  Specif_Gauss->SetTextFont(TextFont);
  Specif_Gauss->SetTextSize(SizeSpecif*SF);
  if(FitType==0) Specif_Gauss->Draw("same");

  TLatex *Specif_Kt3;
  TLatex *Specif_kt;
  if(KT3Bin==0) {Specif_Kt3 = new TLatex(0.1, 0.6,"0.16<K_{T,3}<0.3, <k_{T}>=0.24 GeV/c"); Specif_kt = new TLatex(0.1, 0.5,"0.2<k_{T}<0.3, <k_{T}>=0.25 GeV/c");}
  if(KT3Bin==1) {Specif_Kt3 = new TLatex(0.1, 0.6,"0.3<K_{T,3}<1.0, <k_{T}>=0.39 GeV/c"); Specif_kt = new TLatex(0.1, 0.5,"0.3<k_{T}<0.4, <k_{T}>=0.35 GeV/c");}
  //if(KT3Bin==0) {Specif_Kt3 = new TLatex(0.57, 0.83,"0.16<K_{T,3}<0.25 GeV/c"); Specif_kt = new TLatex(0.57, 0.76,"0.2<k_{T}<0.3 GeV/c");}
  //if(KT3Bin==1) {Specif_Kt3 = new TLatex(0.57, 0.83,"0.25<K_{T,3}<0.35 GeV/c"); Specif_kt = new TLatex(0.57, 0.76,"0.3<k_{T}<0.4 GeV/c");}
  //if(KT3Bin==2) {Specif_Kt3 = new TLatex(0.57, 0.83,"0.35<K_{T,3}<1.0 GeV/c"); Specif_kt = new TLatex(0.57, 0.76,"0.4<k_{T}<0.5 GeV/c");}
  Specif_Kt3->SetTextFont(TextFont); Specif_kt->SetTextFont(TextFont);
  Specif_Kt3->SetTextSize(SizeSpecif*SF); Specif_kt->SetTextSize(SizeSpecif*SF);
  Specif_Kt3->SetNDC(); Specif_kt->SetNDC();
  Specif_Kt3->Draw("same");
  Specif_kt->Draw("same");
  
  legend4->SetTextFont(TextFont);
  legend4->SetTextSize(SizeLegend*SF);
  legend4->Draw("same");
  
  ///////////////////////////////////////////////////////////////////
  pad3->cd(2);
  gPad->SetLeftMargin(0.05);
  gPad->SetRightMargin(0.01);
  gPad->SetTopMargin(0.02);
  gPad->SetBottomMargin(0.1);
  gPad->SetTickx(); gPad->SetTicky();
  //gPad->SetGridx(); gPad->SetGridy();
  TH1D *LambdaPbPb=(TH1D*)Parameters_c3[0][FitType][KT3Bin][1]->Clone();
  TH1D *LambdapPb=(TH1D*)Parameters_c3[1][FitType][KT3Bin][1]->Clone();
  TH1D *Lambdapp=(TH1D*)Parameters_c3[2][FitType][KT3Bin][1]->Clone();
  
  LambdaPbPb->GetXaxis()->SetLabelFont(TextFont); LambdaPbPb->GetYaxis()->SetLabelFont(TextFont); 
  LambdaPbPb->GetXaxis()->SetLabelSize(SizeLabel*SF); LambdaPbPb->GetYaxis()->SetLabelSize(SizeLabel*SF);
  LambdaPbPb->GetXaxis()->SetNdivisions(808);
  LambdaPbPb->GetYaxis()->SetNdivisions(606);
  LambdaPbPb->GetXaxis()->SetTitleFont(TextFont); LambdaPbPb->GetXaxis()->SetTitleSize(SizeTitle*SF);
  LambdaPbPb->GetYaxis()->SetTitleFont(TextFont); LambdaPbPb->GetYaxis()->SetTitleSize(SizeTitle*SF);
  LambdaPbPb->SetMaximum(2.3);
  LambdaPbPb->GetXaxis()->SetTitleOffset(1.1);
  LambdaPbPb->GetYaxis()->SetTitleOffset(100);//1.1
  
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
  TH1D *LambdaC2PbPb=(TH1D*)Parameters_C2[0][FitType][KT3Bin][1]->Clone();
  TH1D *LambdaC2pPb=(TH1D*)Parameters_C2[1][FitType][KT3Bin][1]->Clone();
  TH1D *LambdaC2pp=(TH1D*)Parameters_C2[2][FitType][KT3Bin][1]->Clone();
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
  if(FitType==1){
    grLambdaC2Sys_pp->Draw("E2 p");
    grLambdaC2Sys_pPb->Draw("E2 p");
  }
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
  if(FitType==1){
    LambdaC2pPb->Draw("same");
    LambdaC2pp->Draw("same");
  }
  TLatex *Specif_Marker_2 = new TLatex(0.55,0.9,"Hollow=#lambda  Solid=#lambda_{3}");
  Specif_Marker_2->SetNDC();
  Specif_Marker_2->SetTextFont(TextFont);
  Specif_Marker_2->SetTextSize(SizeSpecif*SF);
  Specif_Marker_2->Draw("same");
  
  can3->cd();
  TLatex *RinvTitle;
  if(FitType==0) RinvTitle=new TLatex(0.04,0.8,"#font[12]{R}_{inv,2} or #font[12]{R}_{inv,3} (fm)");
  else RinvTitle=new TLatex(0.04,0.8,"#font[12]{R}^{#font[12]{Ew}}_{inv,2} or #font[12]{R}^{#font[12]{Ew}}_{inv,3} (fm)");
  RinvTitle->SetNDC();
  RinvTitle->SetTextFont(TextFont);
  RinvTitle->SetTextSize(SizeTitle*SF);
  RinvTitle->SetTextAngle(90);
  RinvTitle->Draw("same");
  TLatex *LambdaTitle;
  if(FitType==0) LambdaTitle=new TLatex(0.04,0.4,"#lambda or #lambda_{3}");
  else LambdaTitle=new TLatex(0.04,0.37,"#lambda^{#font[12]{Ew}} or #lambda^{#font[12]{Ew}}_{3}");
  LambdaTitle->SetNDC();
  LambdaTitle->SetTextFont(TextFont);
  LambdaTitle->SetTextSize(SizeTitle*SF);
  LambdaTitle->SetTextAngle(90);
  LambdaTitle->Draw("same");
  

  if(SaveFiles && FitType==0) can3->SaveAs("ThreePionFitParametersGauss.eps");
  if(SaveFiles && FitType==1) can3->SaveAs("ThreePionFitParametersEW.eps");
  
  

  if(KT3Bin>0) {cout<<"Only print Radii/Lambda Plot for this setting"<<endl; return;}

  ////////////////////////////////////////////////////
  ////////////////////////////////////////////////////
  // Correlation functions and Monte Carlo
  SF=2;
  /*
  TCanvas *can4 = new TCanvas("can4", "can4",10,0,900,600);// 11,53,700,500
  gStyle->SetOptFit(0111);
  can4->SetFillColor(10);//10
  can4->SetBorderMode(0);
  //can4->SetBorderSize(2);
  can4->SetFrameFillColor(0);
  can4->SetFrameBorderMode(0);
  can4->SetFrameBorderMode(0);
  can4->cd();
  PadLeftMargin=0.06; PadBottomMargin=0.05*3/2.;
  float CanTRMarg=0.005;
  TPad *pad4 = new TPad("pad4","pad4",PadLeftMargin,PadBottomMargin,1.0,1.0);
  gPad->SetGridx(0);
  gPad->SetGridy(0);
  gPad->SetTickx();
  gPad->SetTicky();
  pad4->SetTopMargin(0.);
  pad4->SetBottomMargin(0.0);//0.15
  pad4->SetRightMargin(0.0);
  pad4->SetLeftMargin(0.0);//0.01
  pad4->Divide(3,2,0,0);
  pad4->Draw();
  */
  TCanvas *can4 = (TCanvas*)make_canvas("can4","can4",3,2,0,900,600);
  can4->Draw();

  TLegend *legend5[6];
  legend5[0] = new TLegend(.4,.6, .97,.97,NULL,"brNDC");//.45 or .4 for x1
  legend5[0]->SetBorderSize(0);
  legend5[0]->SetFillColor(0);
  legend5[0]->SetTextFont(TextFont);
  legend5[1]=(TLegend*)legend5[0]->Clone();
  legend5[2]=(TLegend*)legend5[0]->Clone();
  legend5[3]=(TLegend*)legend5[0]->Clone();
  legend5[4]=(TLegend*)legend5[0]->Clone();
  legend5[5]=(TLegend*)legend5[0]->Clone();
  //
  //TGaxis *Xaxes[3];
  //TGaxis *Yaxes[2];
  //float AxesLimitsX[3][2]={{0}};
  //float AxesLimitsY[2][2]={{0}};
  double HIJING_c3_SC_K1_M3[15]={0, 0.851168, 0.979088, 1.0209, 0.976797, 1.01555, 0.992262, 1.00773, 0.991634, 0.991504, 0.997317, 0.993006, 0.99694, 0.999369, 0.998404};
  double HIJING_c3_SC_K1_M3_e[15]={0, 0.546937, 0.118551, 0.0436675, 0.0226652, 0.0139659, 0.00906562, 0.00649369, 0.00488794, 0.00380819, 0.00306916, 0.00255166, 0.00219781, 0.00235171, 0.00292962};
  double HIJING_c3_MC_K1_M3[15]={0, 0.886712, 1.02583, 0.985831, 1.00453, 1.01572, 1.00153, 0.991872, 0.997636, 0.997151, 0.996838, 0.999562, 0.998487, 0.996162, 1.001};
  double HIJING_c3_MC_K1_M3_e[15]={0, 0.190786, 0.0527107, 0.0220311, 0.0119954, 0.00756783, 0.00499146, 0.00360927, 0.0027377, 0.00214376, 0.00173608, 0.00144723, 0.00124891, 0.00133898, 0.0016771};
  TH1D *HIJING_c3_SC=(TH1D*)c3[1][1][0][0][17]->Clone();// choose any one that exists to copy histo attributes
  TH1D *HIJING_c3_MC=(TH1D*)c3[1][1][0][0][17]->Clone();// choose any one that exists to copy histo attributes
  for(int i=1; i<=15; i++){
    HIJING_c3_SC->SetBinContent(i, HIJING_c3_SC_K1_M3[i-1]);
    HIJING_c3_SC->SetBinError(i, HIJING_c3_SC_K1_M3_e[i-1]);
    HIJING_c3_MC->SetBinContent(i, HIJING_c3_MC_K1_M3[i-1]);
    HIJING_c3_MC->SetBinError(i, HIJING_c3_MC_K1_M3_e[i-1]);
  }
  //
  for(int padNum=1; padNum<=6; padNum++){
    
    /*pad4->cd(padNum);
    
    if(padNum==1 || padNum==4) {gPad->SetLeftMargin(0.0);}
    else{gPad->SetLeftMargin(0);}
    
    if(padNum>=4){gPad->SetBottomMargin(CanTRMarg); gPad->SetTopMargin(0);}
    else{gPad->SetBottomMargin(0); gPad->SetTopMargin(CanTRMarg);}
    
    if(padNum==3 || padNum==6) gPad->SetRightMargin(0.01);
    else gPad->SetRightMargin(0);
    */
    can4->cd(padNum);
    if(padNum==3 || padNum==6) gPad->SetRightMargin(0.005);
    float SF_6pannel=2;
    
    
    //
    int System_proof=0;
    int ChComb_proof=0;
    int Mb_proof=0;
    if(padNum==1) {System_proof=2; ChComb_proof=0; Mb_proof=18;}
    if(padNum==2) {System_proof=1; ChComb_proof=0; Mb_proof=14;}
    if(padNum==3) {System_proof=0; ChComb_proof=0; Mb_proof=3;}
    if(padNum==4) {System_proof=2; ChComb_proof=1; Mb_proof=18;}
    if(padNum==5) {System_proof=1; ChComb_proof=1; Mb_proof=14;}
    if(padNum==6) {System_proof=0; ChComb_proof=1; Mb_proof=3;}
    
    C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->SetMinimum(0.9); C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->SetMaximum(3.1);
    //
    //c3[2][0][0][KT3Bin][Mb_pp]->GetXaxis()->SetLabelSize(SizeLabel*SF2); c3[2][0][0][KT3Bin][Mb_pp]->GetYaxis()->SetLabelSize(SizeLabel*SF2);
    //c3[2][0][0][KT3Bin][Mb_pp]->GetXaxis()->SetNdivisions(808);
    //c3[2][0][0][KT3Bin][Mb_pp]->GetYaxis()->SetTitleSize(SizeTitle*SF1);
    C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetYaxis()->SetTitle("#font[12]{C}_{3} or #font[12]{#bf{c}}_{3} ");
    if(padNum<=5){
      C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetXaxis()->SetTitleOffset(10); 
    }else C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetXaxis()->SetTitleOffset(1.1);
    C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetYaxis()->SetTitleOffset(1.1);
    
    C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetXaxis()->SetLabelSize(SizeLabel*SF);
    C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetYaxis()->SetLabelSize(SizeLabel*SF);
    C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetXaxis()->SetTitleSize(SizeTitle*SF);
    C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetYaxis()->SetTitleSize(SizeTitle*SF);
    C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetXaxis()->SetNdivisions(606);
    C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetYaxis()->SetNdivisions(606);
    
    double Q3Limit;
    //if(System_proof==0) Q3Limit = 0.1 + 0.1*Mb_proof/16.;
    //else Q3Limit = 0.3 + 0.2*fabs(Mb_proof-10)/9.;
    if(System_proof==1 || System_proof==2) Q3Limit = 0.49;
    else Q3Limit = 0.1099;
    C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetXaxis()->SetRangeUser(0,Q3Limit);
    C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->DrawCopy();
    
    /*
    if(padNum==1) {
      AxesLimitsY[1][0]=C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetMinimum();
      AxesLimitsY[1][1]=C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetMaximum();
    }
    if(padNum==4) {
      AxesLimitsX[0][0]=0;
      AxesLimitsX[0][1]=Q3Limit;
      AxesLimitsY[0][0]=C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetMinimum();
      AxesLimitsY[0][1]=C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetMaximum();
    }
    if(padNum>4) {
      AxesLimitsX[padNum-4][0]=0;
      AxesLimitsX[padNum-4][1]=Q3Limit;
    }
    */
    //if(padNum>=4) Xaxes[padNum-4]=(TGaxis*)(C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetXaxis())->Clone();
    //if(padNum==1) Yaxes[1]=(TGaxis*)(C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetYaxis())->Clone();
    //if(padNum==4) Yaxes[0]=(TGaxis*)(C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetYaxis())->Clone();
    //TGaxis *axis = new TGaxis(gPad->GetUxmin(),gPad->GetUymax(),gPad->GetUxmax(),gPad->GetUymax(),0,10,510,"+L");
    
    C3_Sys[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->DrawCopy("E2 same");
    c3_Sys[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->DrawCopy("E2 same");
    C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->DrawCopy("same");
    c3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->DrawCopy("same");
    if(System_proof!=0) c3[System_proof][1][ChComb_proof][KT3Bin][Mb_proof]->DrawCopy("same");// MC cumulants
    else{
      if(ChComb_proof==0) HIJING_c3_SC->DrawCopy("same");
      else HIJING_c3_MC->DrawCopy("same");
    }
    if(padNum<=3){
      legend5[padNum-1]->AddEntry(C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof],"#font[12]{C}_{3}^{#pm#pm#pm}","p");
      legend5[padNum-1]->AddEntry(c3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof],"#font[12]{#bf{c}}_{3}^{#pm#pm#pm}","p");
      if(System_proof==0) legend5[padNum-1]->AddEntry(HIJING_c3_SC,"HIJING #font[12]{#bf{c}}_{3}^{#pm#pm#pm}","p");
      if(System_proof==1) legend5[padNum-1]->AddEntry(c3[System_proof][1][ChComb_proof][KT3Bin][Mb_proof],"DPMJET #font[12]{#bf{c}}_{3}^{#pm#pm#pm}","p");
      if(System_proof==2) legend5[padNum-1]->AddEntry(c3[System_proof][1][ChComb_proof][KT3Bin][Mb_proof],"Pythia #font[12]{#bf{c}}_{3}^{#pm#pm#pm}","p");
    }else{
      legend5[padNum-1]->AddEntry(C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof],"#font[12]{C}_{3}^{#pm#pm#mp}","p");
      legend5[padNum-1]->AddEntry(c3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof],"#font[12]{#bf{c}}_{3}^{#pm#pm#mp}","p");
      if(System_proof==0) legend5[padNum-1]->AddEntry(HIJING_c3_MC,"HIJING #font[12]{#bf{c}}_{3}^{#pm#pm#mp}","p");
      if(System_proof==1) legend5[padNum-1]->AddEntry(c3[System_proof][1][ChComb_proof][KT3Bin][Mb_proof],"DPMJET #font[12]{#bf{c}}_{3}^{#pm#pm#mp}","p");
      if(System_proof==2) legend5[padNum-1]->AddEntry(c3[System_proof][1][ChComb_proof][KT3Bin][Mb_proof],"Pythia #font[12]{#bf{c}}_{3}^{#pm#pm#mp}","p");
    }
   
    if(ChComb_proof==0) {
      c3_fit[System_proof][0][KT3Bin][Mb_proof]->Draw("c same");
      c3_fit[System_proof][1][KT3Bin][Mb_proof]->Draw("c same");
      legend5[padNum-1]->AddEntry(c3_fit[System_proof][0][KT3Bin][Mb_proof],"Gaussian Fit","l");
      legend5[padNum-1]->AddEntry(c3_fit[System_proof][1][KT3Bin][Mb_proof],"Edgeworth Fit","l");
    }
    double SF_correction=1.0;
    if(padNum==4 || padNum==5) SF_correction=0.95;
    if(padNum==6) SF_correction=0.92;


    TLatex *CTLabel;
    if(System_proof==0) CTLabel = new TLatex(0.45,0.52,"Pb-Pb #sqrt{s_{NN}}=2.76 TeV");// 0.003,.988
    if(System_proof==1) CTLabel = new TLatex(0.45,0.52,"p-Pb #sqrt{s_{NN}}=5.02 TeV");// 0.003,.988
    if(System_proof==2) CTLabel = new TLatex(0.65,0.52,"pp #sqrt{s}=7 TeV");// 0.003,.988
    CTLabel->SetNDC();
    CTLabel->SetTextFont(TextFont);
    CTLabel->SetTextSize(SizeSpecif*SF_6pannel*SF_correction);
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
    TLatex *MLabel;
    if(System_proof!=2) MLabel = new TLatex(0.45,0.4,nameCh->Data());
    else MLabel = new TLatex(0.65,0.4,nameCh->Data());
    MLabel->SetNDC();
    MLabel->SetTextFont(TextFont);
    MLabel->SetTextSize(SizeSpecif*SF_6pannel*SF_correction);
    MLabel->Draw("same");
    
    
    legend5[padNum-1]->SetTextSize(SizeLegend*SF_6pannel);
    if(padNum==1 || padNum==4) legend5[padNum-1]->SetX1(0.55);
    else legend5[padNum-1]->SetX1(0.45);
    legend5[padNum-1]->Draw("same");
  }
  
  

  can4->cd();
  
  TPad *pad4_2 = new TPad("pad4_2","pad4_2",0.0,0.0,1.,1.);
  pad4_2->SetFillStyle(0);
  pad4_2->Draw();
  pad4_2->cd();
  TBox *CoverUp1 = new TBox(0.35,0.05,0.38,.075);
  CoverUp1->SetFillColor(10);
  CoverUp1->Draw();
  TBox *CoverUp2 = new TBox(0.66,0.05,0.69,.075);
  CoverUp2->SetFillColor(10);
  CoverUp2->Draw();
  /*
  TLatex *Q3Title=new TLatex(0.9,0.02,"#font[12]{Q}_{3} (GeV/#font[12]{c})");
  Q3Title->SetNDC();
  Q3Title->SetTextFont(TextFont);
  Q3Title->SetTextSize(SizeTitle);
  Q3Title->Draw("same");
  
  TLatex *C3c3Title=new TLatex(0.02,0.9,"#font[12]{C}_{3} or #font[12]{#bf{c}}_{3}");
  C3c3Title->SetNDC();
  C3c3Title->SetTextFont(TextFont);
  C3c3Title->SetTextSize(SizeTitle);
  C3c3Title->SetTextAngle(90);
  C3c3Title->Draw("same");


  
  Xaxes[0]=new TGaxis(PadLeftMargin,PadBottomMargin+CanTRMarg, PadLeftMargin+.333*(1-PadLeftMargin),PadBottomMargin+CanTRMarg, AxesLimitsX[0][0],AxesLimitsX[0][1],606,"");
  Xaxes[0]->SetLabelFont(TextFont);
  Xaxes[0]->SetLabelSize(SizeLabel);
  Xaxes[0]->Draw();
  //
  Xaxes[1]=new TGaxis(PadLeftMargin+.3333*(1-PadLeftMargin),PadBottomMargin+CanTRMarg, PadLeftMargin+.666*(1-PadLeftMargin),PadBottomMargin+CanTRMarg, AxesLimitsX[1][0],AxesLimitsX[1][1],606,"+L");
  Xaxes[1]->SetLabelFont(TextFont);
  Xaxes[1]->SetLabelSize(SizeLabel);
  Xaxes[1]->Draw();
  //
  Xaxes[2]=new TGaxis(PadLeftMargin+.667*(1-PadLeftMargin),PadBottomMargin+CanTRMarg, PadLeftMargin+(1-PadLeftMargin),PadBottomMargin+CanTRMarg, AxesLimitsX[2][0],AxesLimitsX[2][1],606,"+L");
  Xaxes[2]->SetLabelFont(TextFont);
  Xaxes[2]->SetLabelSize(SizeLabel);
  Xaxes[2]->Draw();
  //
  Yaxes[0]=new TGaxis(PadLeftMargin,PadBottomMargin+CanTRMarg, PadLeftMargin,PadBottomMargin+CanTRMarg+0.5*(1-PadBottomMargin), AxesLimitsY[0][0],AxesLimitsY[0][1],606,"+L");
  Yaxes[0]->SetLabelFont(TextFont);
  Yaxes[0]->SetLabelSize(SizeLabel);
  Yaxes[0]->SetOption("-");
  Yaxes[0]->Draw();
  //
  Yaxes[1]=new TGaxis(PadLeftMargin,PadBottomMargin+0.5*(1-PadBottomMargin), PadLeftMargin,PadBottomMargin+(1-PadBottomMargin), AxesLimitsY[1][0],AxesLimitsY[1][1],606,"+L");
  Yaxes[1]->SetLabelFont(TextFont);
  Yaxes[1]->SetLabelSize(SizeLabel);
  Yaxes[1]->SetOption("-");
  Yaxes[1]->Draw();
  */
  

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
  


  

 
  
  ////////////////////////////////////////////////////
  ////////////////////////////////////////////////////
  // 2 system correlation function comparison
  SF=1.1;
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
  pad5->SetTopMargin(0.0);//0.05
  pad5->SetRightMargin(0.0);//3e-2
  pad5->SetBottomMargin(0.0);//0.12
  pad5->SetLeftMargin(0.0);//0.12
  pad5->Draw();
  pad5->cd();
  TLegend *legend6 = new TLegend(.3,.75, .8,.95,NULL,"brNDC");//.45 or .4 for x1
  legend6->SetBorderSize(0);
  legend6->SetFillColor(0);
  legend6->SetTextFont(TextFont);
  legend6->SetTextSize(SizeLegend*SF);
  //
  gPad->SetRightMargin(0.03); gPad->SetLeftMargin(0.09);
  gPad->SetBottomMargin(0.09); gPad->SetTopMargin(0.02);
  gPad->SetTickx(); gPad->SetTicky();
  //
  int KT3Bin_CorrComp=0;
  int Mbin_SysCompPbPb=12;// 12
  int Mbin_SysComppPb=16;// 12, 16,
  int Mbin_SysComppp=15;//   ,  15,
  TH1D *c3_PbPb=(TH1D*)c3[0][0][0][KT3Bin_CorrComp][Mbin_SysCompPbPb]->Clone();
  TH1D *c3_pPb=(TH1D*)c3[1][0][0][KT3Bin_CorrComp][Mbin_SysComppPb]->Clone();
  TH1D *c3_pp=(TH1D*)c3[2][0][0][KT3Bin_CorrComp][Mbin_SysComppp]->Clone();

  c3_pPb->GetYaxis()->SetTitle("#bf{c}_{3}^{#pm#pm#pm}");
  c3_pPb->GetXaxis()->SetRangeUser(0,0.5);// pp and pPb
  //c3_pPb->GetXaxis()->SetRangeUser(0,0.27);// PbPb and pPb
  c3_pPb->GetXaxis()->SetLabelSize(SizeLabel*SF); c3_pp->GetYaxis()->SetLabelSize(SizeLabel*SF);
  c3_pPb->GetXaxis()->SetTitleSize(SizeLabel*SF); c3_pp->GetYaxis()->SetTitleSize(SizeLabel*SF);
  c3_pPb->GetYaxis()->SetTitle("#font[12]{#bf{c}}_{3}^{#pm#pm#pm}");
  c3_pPb->GetXaxis()->SetTitleOffset(1.2); c3_pp->GetYaxis()->SetTitleOffset(1.3);
  c3_pPb->SetMinimum(0.9); c3_pPb->SetMaximum(3.1);
  c3_pPb->GetXaxis()->SetNdivisions(606);
  c3_pPb->Draw();
  //
  c3_Sys[2][0][0][KT3Bin_CorrComp][Mbin_SysComppp]->Draw("E2 same");
  c3_Sys[1][0][0][KT3Bin_CorrComp][Mbin_SysComppPb]->Draw("E2 same");
  //c3_Sys[0][0][0][KT3Bin][Mbin_SysCompPbPb]->Draw("E2 same");
  c3_pp->Draw("same");
  c3_pPb->Draw("same");
  //c3_PbPb->Draw("same");
  // SysPercent = 0.07 - 0.05*(Mb_proof/19.);
  
  legend6->AddEntry(c3_pPb,"p-Pb #sqrt{s_{NN}}=5.02 TeV, N_{ch} = 24#pm1.2","p");// MpPb=16
  legend6->AddEntry(c3_pp,"pp #sqrt{s}=7 TeV, N_{ch} = 27#pm1.4","p");// Mpp=15
  //legend6->AddEntry(c3_pPb,"p-Pb #sqrt{s_{NN}}=5.02 TeV, N_{ch} = 59#pm3.0","p");// MpPb=12
  //legend6->AddEntry(c3_PbPb,"Pb-Pb #sqrt{s}=7 TeV, N_{ch} = 69#pm3.5","p");// MPbPb=12

  //c3_fit[0][0][KT3Bin][Mbin_SysCompPbPb]->Draw("c same");
  c3_fit[1][0][KT3Bin][Mbin_SysComppPb]->Draw("c same");
  c3_fit[2][0][KT3Bin][Mbin_SysComppp]->Draw("c same");
  legend6->Draw("same");

  
  TLatex *Specif_Kt3_4 = new TLatex(0.425,0.67,"0.16<K_{T,3}<0.3 GeV/c");
  Specif_Kt3_4->SetNDC();
  Specif_Kt3_4->SetTextFont(TextFont);
  Specif_Kt3_4->SetTextSize(SizeTitle*SF);
  Specif_Kt3_4->Draw("same");
  

  //
  if(SaveFiles) {
    TString *name = new TString("c3SystemComp");
    if(FitType==1) name->Append("EW");
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

 
  //TF1 *Poly=new TF1("Poly","pol1",0,5);
  //TH1D *Combined_kappaPlot=(TH1D*)Parameters_c3[0][1][KT3Bin][3]->Clone();
  //Combined_kappaPlot->Add(Parameters_c3[1][1][KT3Bin][3]);
  //Combined_kappaPlot->Add(Parameters_c3[2][1][KT3Bin][3]);
  //Combined_kappaPlot->Fit(Poly,"IME","",0.6,3.6);
  //Parameters_c3[0][1][KT3Bin][3]->Draw("same");
  //Parameters_c3[1][1][KT3Bin][3]->Draw("same");
  //Parameters_c3[2][1][KT3Bin][3]->Draw("same");
  //Poly->Draw("same");


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
//____________________________________________________________________________________________________
TCanvas *make_canvas(const Char_t *name, const Char_t *title, Int_t n_x, Int_t n_y, Int_t share_axes, Int_t ww, Int_t wh){
  const Int_t width=350;
  const Int_t height=350;
  TCanvas *canvas;
 
  if (share_axes==0) {
    if (ww==0)
      ww=width*n_x;
    if (wh==0)
      wh=height*n_y;
    canvas=new TCanvas(name,title,10,0,ww,wh);
    canvas->SetTopMargin(0.01);
    canvas->SetRightMargin(0.01);
    canvas->SetLeftMargin(0.13); 
    canvas->SetBottomMargin(0.13); 
    
    canvas->Divide(n_x,n_y,0,0);
  }
  else {
    Float_t pix_width=(1-gStyle->GetPadLeftMargin()-gStyle->GetPadRightMargin())*width;
    Float_t pix_height=(1-gStyle->GetPadTopMargin()-gStyle->GetPadBottomMargin())*height;
    if (ww==0)
      ww=width+(n_x-1)*pix_width;
    if (wh==0)
      wh=height+(n_y-1)*pix_height;

    canvas=new TCanvas(name,title,ww,wh);
    canvas->SetTopMargin(0); canvas->SetBottomMargin(0); 
    canvas->SetLeftMargin(0); canvas->SetRightMargin(0); 
    Float_t tot_width;
    if (n_x>1) 
      tot_width=(n_x-2)+1./(1-gStyle->GetPadLeftMargin())
	+1./(1-gStyle->GetPadRightMargin());
    else
      tot_width=1./(1-gStyle->GetPadLeftMargin()-gStyle->GetPadRightMargin());
    Float_t tot_height;
    if (n_y>1) 
      tot_height=(n_y-2)+1./(1-gStyle->GetPadTopMargin())
	+1./(1-gStyle->GetPadBottomMargin());
    else
      tot_height=1./(1-gStyle->GetPadTopMargin()-gStyle->GetPadBottomMargin());

    //Int_t idx=n_x*n_y;
    
    for (Int_t j=0; j<n_y; j++) {
      for (Int_t i=0; i<n_x; i++) {
	//for (Int_t j=n_y-1; j>=0; j--) {
	//for (Int_t i=n_x-1; i>=0; i--) {
	Char_t tmp_str[256];
	Char_t p_title[256];
	Int_t idx=n_x*j+i+1;
	sprintf(tmp_str,"%s_%d",canvas->GetName(),idx);
	sprintf(p_title,"Pad %d",idx);
	Float_t x1=0,y1=0;
	Float_t x2=1,y2=1;
	if (n_x>1) {
	  if (i==0) 
	    x2=1./(1-gStyle->GetPadLeftMargin())/tot_width;
	  else {
	    x1=(1./(1-gStyle->GetPadLeftMargin())+i-1)/tot_width;
	    if (i<n_x-1)
	      x2=(1./(1-gStyle->GetPadLeftMargin())+i)/tot_width;
	  }
	}
	if (n_y>1) {
	  if (j==0) 
	    y1=1-1./(1-gStyle->GetPadTopMargin())/tot_height;
	  else {
	    y2=1-(1./(1-gStyle->GetPadTopMargin())+j-1)/tot_height;
	    if (j<n_y-1)
	      y1=1-(1./(1-gStyle->GetPadTopMargin())+j)/tot_height;
	  }
	}
	//cout << "x1 " << x1 << ", x2 " << x2 << endl;
	TPad *pad=new TPad(tmp_str,title,x1,y1,x2,y2);
	//pad->SetFillColor(idx+1);
	if (i>0)
	  pad->SetLeftMargin(0.001);
	if (i<n_x-1)
	  pad->SetRightMargin(0.001);
	if (j>0)
	  pad->SetTopMargin(0.001);
	if (j<n_y-1)
	  pad->SetBottomMargin(0.001);


	pad->SetNumber(idx);
	//pad->SetNumber(n_x*j+i+1);
	pad->Draw();
	//idx--;
	//idx++;
	//pad->SetP
      }
    }
  }

  return canvas;
}
