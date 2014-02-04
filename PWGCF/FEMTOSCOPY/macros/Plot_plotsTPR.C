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
const int KT3Bin=0;// Kt3 bin. 0-1
int FitType=1;// 0 (Gaussian), 1 (Edgeworth)
bool p_pPb_Comp=0;
bool AddedCC=kTRUE;// Charge Conjugate already added?
//
int MbinMaxPbPb=15;// 15
int MbinMinpPb=11;// 13
int MbinMinpp=13;// 14 
int MbinMinPbPb=0;// 0
int MbinMaxpPb=18;// 18
//
//
//
const int MaxKT3Bins=2;
int TextFont=42;// 63, or 42
float SizeLabel=0.06;// 20(63 font), 0.08(42 font)
float SizeLegend=0.05;// 
float SizeTitle=0.06;// 
float SizeSpecif=0.05;// 




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
  //gStyle->SetOptFit(0111);

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
  NrecMapFile = new TFile("Results/NrecMapping_12a17a.root","READ");// standard (with P < 1.0 cut)
  //NrecMapFile = new TFile("Results/Old_NrecMappingFiles/NrecMapping_12a17a.root","READ");// paper v1 file (without P < 1.0 cut)
  //NrecMapFile = new TFile("Results/NrecMapping_12a11a.root","READ");// MC variation
  //NrecMapFile = new TFile("Results/NrecMapping_12a17a_FB5and7overlap.root","READ");
  MyList=(TList*)NrecMapFile->Get("MyList");
  NrecMap = (TH2D*)MyList->FindObject("fNchTrueDist");
  for(int bin=1; bin<=2; bin++){// 1 to 2 (FB7),  1 to 1 (AMPT),  1 to 4 (FB5and7overlap)
    NrecMap->GetXaxis()->SetRangeUser(bin,bin);
    //if(NrecMap->GetMean(2)>0) meanNchPbPb[bin-1] = log10(NrecMap->GetMean(2));
    if(NrecMap->GetMean(2)>0) {
      meanNchPbPb[bin-1] = NrecMap->GetMean(2);
      //cout<<NrecMap->GetMean(2)<<"  "<<fabs(NrecMap->GetRMS(2))/NrecMap->GetMean(2)<<endl;
    }
    //if(NrecMap->GetMean(2)>0) meanNchPbPb[bin-1] = pow(NrecMap->GetMean(2),1/3.);
  }
  NrecMapFile->Close();
  //
  NrecMapFile = new TFile("Results/NrecMapping_12a17e.root","READ");// standard (with P < 1.0 cut)
  //NrecMapFile = new TFile("Results/Old_NrecMappingFiles/NrecMapping_12a17e.root","READ");// paper v1 file (without P < 1.0 cut)
  //NrecMapFile = new TFile("Results/NrecMapping_12a11b.root","READ");// MC variation
  //NrecMapFile = new TFile("Results/NrecMapping_12a17e_FB5and7overlap.root","READ");
  MyList=(TList*)NrecMapFile->Get("MyList");
  NrecMap = (TH2D*)MyList->FindObject("fNchTrueDist");
  for(int bin=3; bin<=10; bin++){// 3 to 10 (FB7),  2 to 3 (AMPT), 5 to 12 (FB5and7overlap)
    NrecMap->GetXaxis()->SetRangeUser(bin,bin);
    //if(NrecMap->GetMean(2)>0) meanNchPbPb[bin-1] = log10(NrecMap->GetMean(2));
    if(NrecMap->GetMean(2)>0) {
      meanNchPbPb[bin-1] = NrecMap->GetMean(2);
      //cout<<NrecMap->GetMean(2)<<"  "<<fabs(NrecMap->GetRMS(2))/NrecMap->GetMean(2)<<endl;
    }
    //if(NrecMap->GetMean(2)>0) meanNchPbPb[bin-1] = pow(NrecMap->GetMean(2),1/3.);
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
  NrecMapFile = new TFile("Results/NrecMapping_12a17c.root","READ");// standard (with P < 1.0 cut)
  //NrecMapFile = new TFile("Results/Old_NrecMappingFiles/NrecMapping_12a17c.root","READ");// paper v1 file (without P < 1.0 cut)
  //NrecMapFile = new TFile("Results/NrecMapping_12a11g.root","READ");// MC variation
  //NrecMapFile = new TFile("Results/NrecMapping_12a17c_FB5and7overlap.root","READ");
  MyList=(TList*)NrecMapFile->Get("MyList");
  NrecMap = (TH2D*)MyList->FindObject("fNchTrueDist");
  for(int bin=11; bin<=19; bin++){// 11 to 19 (FB7),  1 to 1 (AMPT), 13 to 19 (FB5and7overlap)
    NrecMap->GetXaxis()->SetRangeUser(bin,bin);
    //if(NrecMap->GetMean(2)>0) meanNchPbPb[bin-1] = log10(NrecMap->GetMean(2));
    if(NrecMap->GetMean(2)>0) {
      meanNchPbPb[bin-1] = NrecMap->GetMean(2);
      //cout<<NrecMap->GetMean(2)<<"  "<<fabs(NrecMap->GetRMS(2))/NrecMap->GetMean(2)<<endl;
    }    
    //if(NrecMap->GetMean(2)>0) meanNchPbPb[bin-1] = pow(NrecMap->GetMean(2),1/3.);
  }
  NrecMapFile->Close();
  //cout<<endl;
  //
  NrecMapFile = new TFile("Results/NrecMapping_13b2_efix_p1.root","READ");// standard (with P < 1.0 cut)
  //NrecMapFile = new TFile("Results/PDC_13b2_efix_p1_R2.root","READ");// paper v1 file (without P < 1.0 cut)
  //NrecMapFile = new TFile("Results/NrecMapping_13b3.root","READ");// MC variation
  //NrecMapFile = new TFile("Results/NrecMapping_13b2_efix_p1_FB5and7overlap.root","READ");
  MyList=(TList*)NrecMapFile->Get("MyList");
  NrecMap = (TH2D*)MyList->FindObject("fNchTrueDist");
  for(int bin=10; bin<=20; bin++){
    NrecMap->GetXaxis()->SetRangeUser(bin,bin);
    //if(NrecMap->GetMean(2)>0) meanNchpPb[bin-1] = log10(NrecMap->GetMean(2));
    if(NrecMap->GetMean(2)>0) {
      meanNchpPb[bin-1] = NrecMap->GetMean(2);
      //cout<<NrecMap->GetMean(2)<<"  "<<fabs(NrecMap->GetRMS(2))/NrecMap->GetMean(2)<<endl;
    }
    //if(NrecMap->GetMean(2)>0) meanNchpPb[bin-1] = pow(NrecMap->GetMean(2),1/3.);
  }
  NrecMapFile->Close();
  //
  NrecMapFile = new TFile("Results/NrecMapping_10f6a.root","READ");// standard (with P < 1.0 cut)
  //NrecMapFile = new TFile("Results/PDC_10f6a_R2.root","READ");// paper v1 file (without P < 1.0 cut)
  //NrecMapFile = new TFile("Results/NrecMapping_10f6.root","READ");// MC variation
  //NrecMapFile = new TFile("Results/NrecMapping_10f6a_FB5and7overlap.root","READ");
  MyList=(TList*)NrecMapFile->Get("MyList");
  NrecMap = (TH2D*)MyList->FindObject("fNchTrueDist");
  for(int bin=10; bin<=20; bin++){
    NrecMap->GetXaxis()->SetRangeUser(bin,bin);
    //if(NrecMap->GetMean(2)>0) meanNchpp[bin-1] = log10(NrecMap->GetMean(2));
    if(NrecMap->GetMean(2)>0) {
      meanNchpp[bin-1] = NrecMap->GetMean(2);
      //cout<<NrecMap->GetMean(2)<<"  "<<fabs(NrecMap->GetRMS(2))/NrecMap->GetMean(2)<<endl;
    }
    //if(NrecMap->GetMean(2)>0) meanNchpp[bin-1] = pow(NrecMap->GetMean(2),1/3.);
  }
  NrecMapFile->Close();
  
  //for(int i=0; i<20; i++) cout<<pow(10,meanNchPbPb[i])<<endl;
  //cout<<"+++++++++++++++"<<endl;
  //for(int i=0; i<20; i++) cout<<pow(10,meanNchpPb[i])<<endl;
  //cout<<"+++++++++++++++"<<endl;
  //for(int i=0; i<20; i++) cout<<meanNchpp[i]<<endl;

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
  //TMinuit *Min[3][2][3][20];// CollType, +/-, EDbin, cb
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
	  //Parameters_C2[ct][ft][kt3][par] = new TH1D(name_C2->Data(),"",2900,0.6,3.4);// log10(Nch)
	  //Parameters_c3[ct][ft][kt3][par] = new TH1D(name_c3->Data(),"",2900,0.6,3.4);// log10(Nch)
	  Parameters_C2[ct][ft][kt3][par] = new TH1D(name_C2->Data(),"",30000,1,3001);// Nch
	  Parameters_c3[ct][ft][kt3][par] = new TH1D(name_c3->Data(),"",30000,1,3001);// Nch
	  //Parameters_C2[ct][ft][kt3][par] = new TH1D(name_C2->Data(),"",3000,1,14.4);// Nch^(1/3)
	  //Parameters_c3[ct][ft][kt3][par] = new TH1D(name_c3->Data(),"",3000,1,14.4);// Nch^(1/3)
	}
      }
    }
  }
  
  TH1D *RadiiC2pp_Published = new TH1D("RadiiC2pp_Published","",3000,1.0,3001);
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
		//Min[ct][ch][KT3][cb]=(TMinuit*)file->Get("MyMinuit_c3");
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
		//C3[ct][dt][ChComb][KT3][cb]->GetXaxis()->SetTitleOffset(1.2);
		//C3[ct][dt][ChComb][KT3][cb]->GetYaxis()->SetTitleOffset(1.2);
		C3[ct][dt][ChComb][KT3][cb]->GetYaxis()->SetTitleFont(TextFont);
		c3[ct][dt][ChComb][KT3][cb]->GetXaxis()->SetRangeUser(0,0.5);
		c3[ct][dt][ChComb][KT3][cb]->GetXaxis()->SetLabelFont(TextFont); c3[ct][dt][ChComb][KT3][cb]->GetYaxis()->SetLabelFont(TextFont); 
		//c3[ct][dt][ChComb][KT3][cb]->GetXaxis()->SetTitleOffset(1.2);
		//c3[ct][dt][ChComb][KT3][cb]->GetYaxis()->SetTitleOffset(1.2);
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
		if(ct==0 || ct==1) c3[ct][dt][ChComb][KT3][cb]->SetMarkerSize(1.12*C3[ct][dt][ChComb][KT3][cb]->GetMarkerSize());
		else c3[ct][dt][ChComb][KT3][cb]->SetMarkerSize(1.2*C3[ct][dt][ChComb][KT3][cb]->GetMarkerSize());
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
		    C3_Sys[ct][0][0][KT3][cb]->SetBinError(i, 0.01 * C3_Sys[ct][0][0][KT3][cb]->GetBinContent(i));
		    c3_Sys[ct][0][0][KT3][cb]->SetBinError(i, sqrt(pow(MixedChargeSysFit->Eval(Q3)-1.0,2) + pow(0.1*(c3_Sys[ct][0][0][KT3][cb]->GetBinContent(i)-1.0),2)));// residue + lambda undilution variation (0.7 to 0.65)
		    // MixedCharge
		    C3_Sys[ct][0][1][KT3][cb]->SetBinError(i, 0.01 * C3_Sys[ct][0][1][KT3][cb]->GetBinContent(i));// correlation function uncertainty
		    c3_Sys[ct][0][1][KT3][cb]->SetBinError(i, sqrt(pow(0.01 * c3_Sys[ct][0][1][KT3][cb]->GetBinContent(i),2) + pow(0.1*(c3_Sys[ct][0][1][KT3][cb]->GetBinContent(i)-1.0),2)));// correlation function uncertainty + lambda undilution variation (0.7 to 0.65)
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
		if(FitType==0) Parameters_c3[ct][ft][KT3][par]->GetYaxis()->SetTitle("#font[12]{R}_{inv,2} or #font[12]{R}_{inv,3} (fm)");
		else Parameters_c3[ct][ft][KT3][par]->GetYaxis()->SetTitle("#font[12]{R^{E_{w}}}_{inv,2} or #font[12]{R^{E_{w}}}_{inv,3} (fm)");
	      }		
	      if(par==3) Parameters_c3[ct][ft][KT3][par]->GetYaxis()->SetTitle("#kappa_{3}");
	      if(par==4) Parameters_c3[ct][ft][KT3][par]->GetYaxis()->SetTitle("#kappa_{4}");
	      //Parameters_c3[ct][ft][KT3][par]->GetXaxis()->SetTitle("log_{10}(N_{ch})");
	      Parameters_c3[ct][ft][KT3][par]->GetXaxis()->SetTitle("N_{ch}");
	      //Parameters_c3[ct][ft][KT3][par]->GetXaxis()->SetTitle("N_{ch}^{1/3}");
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
  can2->SetFillColor(0);//10
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
  //legend3->AddEntry(c3[1][0][0][Mb_pPb],"N_{ch} = 59#pm2, p-Pb #sqrt{#font[12]{s}_{NN}}=5.02 TeV","p");
  //legend3->AddEntry(c3[0][0][0][Mb_PbPb],"N_{ch} = 1922#pm135, Pb-Pb #sqrt{#font[12]{s}_{NN}}=2.76 TeV","p");
  legend3->AddEntry(c3[2][0][0][0][Mb_pp],"Low N_{ch}, pp #sqrt{#font[12]{s}}=7 TeV","p");
  legend3->AddEntry(c3[1][0][0][0][Mb_pPb],"Mid N_{ch}, p-Pb #sqrt{#font[12]{#font[12]{s}_{NN}}}=5.02 TeV","p");
  legend3->AddEntry(c3[0][0][0][0][Mb_PbPb],"High N_{ch}, Pb-Pb #sqrt{#font[12]{#font[12]{s}_{NN}}}=2.76 TeV","p");
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
  // K0s removal plot
  TCanvas *can1 = new TCanvas("can1", "can1",10,700,600,600);// 11,53,700,500
  can1->SetHighLightColor(2);
  gStyle->SetOptFit(0111);
  can1->SetFillColor(0);//10
  can1->SetBorderMode(0);
  can1->SetBorderSize(2);
  can1->SetFrameFillColor(0);
  can1->SetFrameBorderMode(0);
  can1->SetFrameBorderMode(0);
  can1->cd();
  TPad *pad1 = new TPad("pad1","pad1",0.,0.,1.,1.);
  gPad->SetTickx();
  gPad->SetTicky();
  pad1->SetTopMargin(0.0);//0.05
  pad1->SetRightMargin(0.0);//3e-2
  pad1->SetBottomMargin(0.0);//0.12
  pad1->SetLeftMargin(0.0);
  pad1->Draw();
  pad1->cd();
  TLegend *legend1 = new TLegend(.2,.65, .55,.85,NULL,"brNDC");//.45 or .4 for x1
  legend1->SetBorderSize(0);
  legend1->SetFillColor(0);
  legend1->SetTextFont(TextFont);
  legend1->SetTextSize(SizeLegend);

  gPad->SetRightMargin(0.01); gPad->SetLeftMargin(0.14);
  gPad->SetBottomMargin(0.14); gPad->SetTopMargin(0.02);
  gPad->SetTickx(); gPad->SetTicky();
  // pp, M18
  double points_K0s_C3[50]={1.44542, 1.33366, 1.25518, 1.21051, 1.20484, 1.18742, 1.18246, 1.17024, 1.16383, 1.1555, 1.15536, 1.15272, 1.14571, 1.14822, 1.14708, 1.14066, 1.14133, 1.13866, 1.13755, 1.13308, 1.1386, 1.13821, 1.13429, 1.13707, 1.13294, 1.13545, 1.13929, 1.13374, 1.13546, 1.13205, 1.12846, 1.1281, 1.12947, 1.13162, 1.13822, 1.13626, 1.13912, 1.14379, 1.15594, 1.19387, 1.29685, 1.21203, 1.14452, 1.14288, 1.14161, 1.1408, 1.13796, 1.13661, 1.13565, 1.13312};
  double points_K0s_C3_e[50]={0.0270571, 0.00940436, 0.0058043, 0.00430829, 0.00355541, 0.00306969, 0.00275566, 0.00253068, 0.00236966, 0.00224498, 0.00215682, 0.00208189, 0.00201555, 0.00196896, 0.00192381, 0.0018782, 0.00184374, 0.0018119, 0.00178103, 0.00175109, 0.00173313, 0.00171119, 0.00168873, 0.00167507, 0.00165525, 0.00164399, 0.00163488, 0.00161686, 0.00160668, 0.00159163, 0.00157479, 0.00155996, 0.00154251, 0.00152405, 0.00150469, 0.00147687, 0.00144939, 0.00142223, 0.0013988, 0.00139214, 0.00142684, 0.0013476, 0.00128368, 0.00126468, 0.0012497, 0.00123827, 0.001228, 0.00122253, 0.00121921, 0.00121782};
  double points_K0s_c3[50]={0.965148, 1.05077, 0.983133, 0.972848, 0.993922, 0.994064, 1.00618, 0.993434, 1.00322, 0.999747, 0.982221, 1.00029, 0.992184, 0.993566, 1.00007, 0.995263, 0.998559, 0.987037, 0.987919, 0.991035, 0.991488, 0.99592, 0.99677, 0.990709, 0.990352, 0.994706, 0.99606, 0.998525, 0.994121, 0.986511, 0.991009, 0.990309, 1.00161, 0.997478, 1.00117, 0.989375, 0.996592, 0.990832, 0.998989, 0.995124, 0.996301, 1.00189, 0.995304, 0.992511, 0.994084, 0.994777, 1.00074, 0.996959, 0.998713, 0.996205};
  double points_K0s_c3_e[50]={0.0630244, 0.0223583, 0.0140589, 0.0105446, 0.00871274, 0.00755392, 0.00678715, 0.00625387, 0.00586333, 0.00556614, 0.00535136, 0.00516553, 0.00501148, 0.00489185, 0.00478054, 0.00467542, 0.00458851, 0.00451377, 0.00443821, 0.00436787, 0.00431762, 0.00426282, 0.0042107, 0.00417461, 0.00412901, 0.00409783, 0.00407137, 0.0040317, 0.00400527, 0.00397212, 0.00393349, 0.00389671, 0.00385012, 0.00380258, 0.00374821, 0.0036816, 0.0036096, 0.00353852, 0.00346903, 0.00342075, 0.00342477, 0.00329675, 0.0031932, 0.00314748, 0.00311129, 0.00308337, 0.00305968, 0.00304743, 0.00303989, 0.00303859};
  TH1D *K0s_C3=new TH1D("K0s_C3","",50,0,0.5);
  TH1D *K0s_c3=new TH1D("K0s_c3","",50,0,0.5);
  for(int i=0; i<50; i++){
    K0s_C3->SetBinContent(i+1, points_K0s_C3[i]);
    K0s_C3->SetBinError(i+1, points_K0s_C3_e[i]);
    K0s_c3->SetBinContent(i+1, points_K0s_c3[i]);
    K0s_c3->SetBinError(i+1, points_K0s_c3_e[i]);
  }
  K0s_c3->SetMarkerStyle(20);
  K0s_c3->SetMarkerColor(4);
  K0s_C3->SetMarkerStyle(24);
  K0s_C3->SetMarkerColor(4);
  K0s_c3->SetMarkerSize(K0s_C3->GetMarkerSize() * 1.12);
  K0s_C3->SetMinimum(0.92);
  K0s_C3->GetXaxis()->SetTitle("#font[12]{q_{31}^{#pm#mp}} (GeV/#font[12]{c})");
  K0s_C3->GetYaxis()->SetTitle("#font[12]{C_{3}} or #font[12]{#bf{c}_{3}}"); 
  K0s_C3->GetXaxis()->SetTitleOffset(1.05); K0s_C3->GetYaxis()->SetTitleOffset(1.12);
  K0s_C3->GetXaxis()->SetTitleSize(SizeTitle); K0s_C3->GetYaxis()->SetTitleSize(SizeTitle);
  K0s_C3->GetXaxis()->SetLabelSize(SizeTitle); K0s_C3->GetYaxis()->SetLabelSize(SizeTitle);
  K0s_C3->GetXaxis()->SetNdivisions(404);
  K0s_C3->GetYaxis()->SetNdivisions(404);
  //K0s_C3->GetXaxis()->SetRangeUser(-0.001,0.5);

  K0s_C3->Draw();
  // Sys errors
  TH1D *Sys_K0s_C3=(TH1D*)K0s_C3->Clone();
  TH1D *Sys_K0s_c3=(TH1D*)K0s_c3->Clone();
  Sys_K0s_C3->SetMarkerSize(0); Sys_K0s_c3->SetMarkerSize(0);
  Sys_K0s_C3->SetFillColor(kBlue-10); Sys_K0s_c3->SetFillColor(kBlue-10);
  Sys_K0s_C3->SetFillStyle(1000); Sys_K0s_c3->SetFillStyle(1000);
  for(int i=0; i<Sys_K0s_C3->GetNbinsX(); i++) { 
    Sys_K0s_C3->SetBinError(i+1, 0.01 * Sys_K0s_C3->GetBinContent(i+1));
    Sys_K0s_c3->SetBinError(i+1, sqrt(pow(0.01 * Sys_K0s_c3->GetBinContent(i+1),2) + pow(0.1*(Sys_K0s_c3->GetBinContent(i+1)-1),2)));
  }
  Sys_K0s_C3->Draw("E2 same");
  Sys_K0s_c3->Draw("E2 same");
  K0s_C3->Draw("same");
  K0s_c3->Draw("same");

  legend1->AddEntry(K0s_C3,"#font[12]{C_{3}^{#pm#pm#mp}} projection","p");
  legend1->AddEntry(K0s_c3,"#font[12]{#bf{c}_{3}^{#pm#pm#mp}} projection","p");
  legend1->Draw("same");
  
  TLatex *Specif_qRange = new TLatex(0.25,0.6,"0.2 < #font[12]{q_{12}^{#pm#pm},q_{23}^{#pm#mp}} < 0.5 GeV/#font[12]{c}");
  Specif_qRange->SetNDC();
  Specif_qRange->SetTextFont(42);
  Specif_qRange->SetTextSize(SizeSpecif);
  Specif_qRange->Draw("same");
  TLatex *Specif_System = new TLatex(0.25,0.9,"pp #sqrt{#font[12]{s}}=7 TeV, #LTN_{ch}#GT = 8.0 #pm 0.4");
  Specif_System->SetNDC();
  Specif_System->SetTextFont(42);
  Specif_System->SetTextSize(SizeSpecif);
  Specif_System->Draw("same");
  




  ////////////////////////////////////////////////////
  ////////////////////////////////////////////////////
  // Radii and lambda plot
  float SF_2panel=1.2;
  float RadiiC2ppPubPoints[2][8]={{0.88,1.09,1.23,1.28,1.34,1.37,1.42,1.48},{0.86,1.06,1.16,1.23,1.23,1.28,1.32,1.38}};
  float RadiiC2ppPubPoints_e[2][8]={{0.058,0.064,0.071,0.071,0.078,0.078,0.086,0.098},{0.12,0.12,0.13,0.13,0.13,0.13,0.13,0.13}};
  float MeanPubNch[8]={6.98,14.95,19.68,24.68,29.38,33.95,38.46,42.66};
  
  TCanvas *can3 = new TCanvas("can3", "can3",1000,0,600,900);// 11,53,700,500
  can3->SetHighLightColor(2);
  gStyle->SetOptFit(0111);
  can3->SetFillColor(0);//10
  can3->SetBorderMode(0);
  can3->SetBorderSize(2);
  can3->SetFrameFillColor(0);
  can3->SetFrameBorderMode(0);
  can3->SetFrameBorderMode(0);
  can3->cd();
  float PadLeftMargin=0.0, PadBottomMargin=0.;
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
  TLegend *legend4 = new TLegend(.15,.65, .55,.95,NULL,"brNDC");//.45 or .4 for x1
  legend4->SetBorderSize(0);
  legend4->SetFillColor(0);
  legend4->SetTextFont(TextFont);
  legend4->SetTextSize(SizeLegend);
  //
  pad3->cd(1);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.01);
  gPad->SetTopMargin(0.01);
  gPad->SetBottomMargin(0.16);
  gPad->SetTickx(); gPad->SetTicky();
  gPad->SetLogx();
  //gPad->SetGridx(); gPad->SetGridy();
  TH1D *RadiiPbPb=(TH1D*)Parameters_c3[0][FitType][KT3Bin][2]->Clone();
  TH1D *RadiipPb=(TH1D*)Parameters_c3[1][FitType][KT3Bin][2]->Clone();
  TH1D *Radiipp=(TH1D*)Parameters_c3[2][FitType][KT3Bin][2]->Clone();
  RadiiPbPb->GetXaxis()->SetLabelFont(TextFont); RadiiPbPb->GetYaxis()->SetLabelFont(TextFont); 
  RadiiPbPb->GetXaxis()->SetLabelSize(SizeLabel*SF_2panel); RadiiPbPb->GetYaxis()->SetLabelSize(SizeLabel*SF_2panel);
  RadiiPbPb->GetXaxis()->SetNdivisions(808);
  RadiiPbPb->GetXaxis()->SetTitleOffset(0.95);//1.3
  RadiiPbPb->GetYaxis()->SetTitleOffset(100);//1.1
  RadiiPbPb->GetXaxis()->SetTitleFont(TextFont); RadiiPbPb->GetXaxis()->SetTitleSize(SizeTitle*SF_2panel);
  RadiiPbPb->GetYaxis()->SetTitleFont(TextFont); RadiiPbPb->GetYaxis()->SetTitleSize(SizeTitle*SF_2panel);
  RadiiPbPb->SetMinimum(0); RadiiPbPb->SetMaximum(11.9);
  //gStyle->SetErrorX(0);
  RadiiPbPb->GetXaxis()->SetRangeUser(3,3000);
  RadiiPbPb->Draw();
  //
  double xAxis_e[20]={0};
  double yAxisPbPb[20]={0};
  double yAxisPbPb_e[20]={0};
  double yAxispPb[20]={0};
  double yAxispPb_e[20]={0};
  double yAxispp[20]={0};
  double yAxispp_e[20]={0};
  
  
  //TF1 *QrangeSys_c3=new TF1("QrangeSys_c3","0.1 - .07*(x/3.3)",0,4);// old method
  // Old method for Rsys was Sys = sqrt(pow(QrangeSys_c3->Eval(meanNchPbPb[cb])*RadiiPbPb->GetBinContent(binPbPb),2) + pow(0.05,2)); 
  // same for pp and pPb but with proper Nch
  double SysPercent_PbPb[2][4]={{12,10,11,16},{5,7,5,10}};// Gauss/EW, parameter#(Rinv2, Rinv3, lambda, lambda_3)
  double SysPercent_pPb[2][4]={{15,10,6,10},{6,2,1,4}};// Gauss/EW, parameter#(Rinv2, Rinv3, lambda, lambda_3)
  double SysPercent_pp[2][4]={{11,9,2,5},{12,5,2,5}};// Gauss/EW, parameter#(Rinv2, Rinv3, lambda, lambda_3)
  for(int cb=0; cb<20; cb++){
    int binPbPb = RadiiPbPb->GetXaxis()->FindBin(meanNchPbPb[cb]);
    int binpPb = RadiipPb->GetXaxis()->FindBin(meanNchpPb[cb]);
    int binpp = Radiipp->GetXaxis()->FindBin(meanNchpp[cb]);
    double RsysPbPb = 0.01*sqrt(pow(SysPercent_PbPb[FitType][1],2) + pow(1,2)) *  RadiiPbPb->GetBinContent(binPbPb);// fit Variation + MRC 
    double RsyspPb = 0.01*sqrt(pow(SysPercent_pPb[FitType][1],2) + pow(1,2)) *  RadiipPb->GetBinContent(binpPb);// fit Variation + MRC 
    double Rsyspp = 0.01*sqrt(pow(SysPercent_pp[FitType][1],2) + pow(1,2)) *  Radiipp->GetBinContent(binpp);// fit Variation + MRC 
    if(RadiiPbPb->GetBinContent(binPbPb)==0) {yAxisPbPb[cb]=100; yAxisPbPb_e[cb]=100;}
    else {yAxisPbPb[cb]=RadiiPbPb->GetBinContent(binPbPb); yAxisPbPb_e[cb]=RsysPbPb;}
    if(RadiipPb->GetBinContent(binpPb)==0) {yAxispPb[cb]=100; yAxispPb_e[cb]=100;}
    else {yAxispPb[cb]=RadiipPb->GetBinContent(binpPb); yAxispPb_e[cb]=RsyspPb;}
    if(Radiipp->GetBinContent(binpp)==0) {yAxispp[cb]=100; yAxispp_e[cb]=100;}
    else {yAxispp[cb]=Radiipp->GetBinContent(binpp); yAxispp_e[cb]=Rsyspp;}
    //cout<<RadiipPb->GetBinContent(binpPb)<<"  "<<Radiipp->GetBinContent(binpp)<<endl;
    // X-axis systematics 5% for all bins
    //if(cb<13) xAxis_e[cb]=fabs(meanNchPbPb[cb] - log10((0.95)*pow(10,meanNchPbPb[cb])));
    //else xAxis_e[cb]=fabs(meanNchpPb[cb] - log10((0.95)*pow(10,meanNchpPb[cb])));
    if(cb<13) xAxis_e[cb]=fabs(meanNchPbPb[cb] - (0.95)*meanNchPbPb[cb]);
    else xAxis_e[cb]=fabs(meanNchpPb[cb] - (0.95)*meanNchpPb[cb]);
  }
 
  TGraphErrors *grRadiiSys_PbPb=new TGraphErrors(20,meanNchPbPb,yAxisPbPb,xAxis_e,yAxisPbPb_e);
  TGraphErrors *grRadiiSys_pPb=new TGraphErrors(20,meanNchpPb,yAxispPb,xAxis_e,yAxispPb_e);
  TGraphErrors *grRadiiSys_pp=new TGraphErrors(20,meanNchpp,yAxispp,xAxis_e,yAxispp_e);
  grRadiiSys_pp->SetMarkerSize(0); grRadiiSys_pp->SetFillStyle(1000); grRadiiSys_pp->SetFillColor(kBlue-10);
  grRadiiSys_pPb->SetMarkerSize(0); grRadiiSys_pPb->SetFillStyle(1000); grRadiiSys_pPb->SetFillColor(kRed-10);
  grRadiiSys_PbPb->SetMarkerSize(0); grRadiiSys_PbPb->SetFillStyle(1000); grRadiiSys_PbPb->SetFillColor(kGray);
  grRadiiSys_pp->SetMarkerColor(kBlue-10); grRadiiSys_pp->SetMarkerColor(kRed-10); grRadiiSys_pp->SetMarkerColor(kGray);
  // C2 
  TH1D *RadiiC2PbPb=(TH1D*)Parameters_C2[0][FitType][KT3Bin][2]->Clone();
  TH1D *RadiiC2pPb=(TH1D*)Parameters_C2[1][FitType][KT3Bin][2]->Clone();
  TH1D *RadiiC2pp=(TH1D*)Parameters_C2[2][FitType][KT3Bin][2]->Clone();
  RadiiC2pp_Published->SetMarkerStyle(30);
  if(FitType==0) RadiiC2pp->SetMarkerStyle(30);

  for(int mbin=0; mbin<8; mbin++){
    //int bin = RadiiC2pp_Published->GetXaxis()->FindBin(log10(MeanPubNch[mbin]));
    int bin = RadiiC2pp_Published->GetXaxis()->FindBin(MeanPubNch[mbin]);
    RadiiC2pp_Published->SetBinContent(bin, RadiiC2ppPubPoints[KT3Bin][mbin]); 
    RadiiC2pp_Published->SetBinError(bin, RadiiC2ppPubPoints_e[KT3Bin][mbin]);
  }  

  for(int cb=0; cb<20; cb++){
    int binPbPb = RadiiC2PbPb->GetXaxis()->FindBin(meanNchPbPb[cb]);
    int binpPb = RadiiC2pPb->GetXaxis()->FindBin(meanNchpPb[cb]);
    int binpp = RadiiC2pp->GetXaxis()->FindBin(meanNchpp[cb]);
    double RsysPbPb = 0.01*sqrt(pow(SysPercent_PbPb[FitType][0],2) + pow(1,2)) *  RadiiC2PbPb->GetBinContent(binPbPb);// fit Variation + MRC
    double RsyspPb = 0.01*sqrt(pow(SysPercent_pPb[FitType][0],2) + pow(1,2)) *  RadiiC2pPb->GetBinContent(binpPb);// fit Variation + MRC
    double Rsyspp = 0.01*sqrt(pow(SysPercent_pp[FitType][0],2) + pow(1,2)) *  RadiiC2pp->GetBinContent(binpp);// fit Variation + MRC
    if(RadiiC2PbPb->GetBinContent(binPbPb)==0) {yAxisPbPb[cb]=100; yAxisPbPb_e[cb]=1000;}
    else {yAxisPbPb[cb]=RadiiC2PbPb->GetBinContent(binPbPb); yAxisPbPb_e[cb]=RsysPbPb;}
    if(RadiiC2pPb->GetBinContent(binpPb)==0) {yAxispPb[cb]=100; yAxispPb_e[cb]=1000;}
    else {yAxispPb[cb]=RadiiC2pPb->GetBinContent(binpPb); yAxispPb_e[cb]=RsyspPb;}
    if(RadiiC2pp->GetBinContent(binpp)==0) {yAxispp[cb]=100; yAxispp_e[cb]=1000;}
    else {yAxispp[cb]=RadiiC2pp->GetBinContent(binpp); yAxispp_e[cb]=Rsyspp;}
  }
  TGraphErrors *grRadiiC2Sys_PbPb=new TGraphErrors(20,meanNchPbPb,yAxisPbPb,xAxis_e,yAxisPbPb_e);
  TGraphErrors *grRadiiC2Sys_pPb=new TGraphErrors(20,meanNchpPb,yAxispPb,xAxis_e,yAxispPb_e);
  TGraphErrors *grRadiiC2Sys_pp=new TGraphErrors(20,meanNchpp,yAxispp,xAxis_e,yAxispp_e);
  grRadiiC2Sys_pp->SetMarkerSize(0); //grRadiiC2Sys_pp->SetFillStyle(3001); grRadiiC2Sys_pp->SetFillColor(kBlue-10);
  grRadiiC2Sys_pPb->SetMarkerSize(0); //grRadiiC2Sys_pPb->SetFillStyle(3001); grRadiiC2Sys_pPb->SetFillColor(kRed-10);
  grRadiiC2Sys_PbPb->SetMarkerSize(0); //grRadiiC2Sys_PbPb->SetFillStyle(3001); grRadiiC2Sys_PbPb->SetFillColor(kGray);
  //grRadiiC2Sys_PbPb->SetMarkerSize(3.3*RadiiC2PbPb->GetMarkerSize());
  //grRadiiC2Sys_PbPb->SetLineWidth(3.3*RadiiC2PbPb->GetLineWidth());
  //
  if(FitType==1){
    //grRadiiC2Sys_pp->Draw("|| p");
    //grRadiiC2Sys_pPb->Draw("|| p");
  }
  grRadiiSys_pp->Draw("E2 p");
  grRadiiSys_pPb->Draw("E2 p");
  grRadiiSys_PbPb->Draw("E2 p");
  RadiiPbPb->Draw("same");
  RadiipPb->Draw("same");
  Radiipp->Draw("same");
  grRadiiC2Sys_PbPb->Draw("|| p");
  if(FitType==0) RadiiC2pp_Published->Draw("same");
  //
  //for(int i=0; i<20; i++){
  //cout<<yAxispPb[i]<<"  "<<yAxispPb_e[i]<<endl;
  //}
  TF1 *PbPbFit=new TF1("PbPbFit","pol1",1,14); PbPbFit->SetLineColor(1);
  //grRadiiSys_PbPb->Fit(PbPbFit,"","",1,13);
  //PbPbFit->Draw("same");
  //cout<<"PbPb: "<<PbPbFit->GetParameter(0)<<", "<<PbPbFit->GetParameter(1)<<endl;
  TF1 *pPbFit=new TF1("pPbFit","pol1",1,14); pPbFit->SetLineColor(2);
  //grRadiiSys_pPb->Fit(pPbFit,"","",2,4.5);
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
  //RadiiC2PbPb->SetMarkerSize(1.3*RadiiC2PbPb->GetMarkerSize());
  
  RadiiC2PbPb->Draw("same");
  if(FitType==1) {
    RadiiC2pPb->Draw("same");
    RadiiC2pp->Draw("same");
  }

  legend4->AddEntry(Radiipp,"pp #sqrt{s}=7 TeV","p");
  legend4->AddEntry(RadiipPb,"p-Pb #sqrt{#font[12]{s}_{NN}}=5.02 TeV","p");
  legend4->AddEntry(RadiiPbPb,"Pb-Pb #sqrt{#font[12]{s}_{NN}}=2.76 TeV","p");

  

  TLatex *Specif_Kt3;
  TLatex *Specif_kt;
  if(KT3Bin==0) {
    Specif_Kt3 = new TLatex(0.17, 0.57,"0.16<#font[12]{K}_{T,3}<0.3 GeV/#font[12]{c}"); 
    Specif_kt = new TLatex(0.17, 0.47,"0.2<#font[12]{k}_{T}<0.3 GeV/#font[12]{c}");
    // KT3:  #LT#font[12]{k}_{T}#GT=0.24 GeV/#font[12]{c}
    // kT:  #LT#font[12]{k}_{T}#GT=0.25 GeV/#font[12]{c}
  }
  if(KT3Bin==1) {
    Specif_Kt3 = new TLatex(0.17, 0.57,"0.3<#font[12]{K}_{T,3}<1.0 GeV/#font[12]{c}"); 
    Specif_kt = new TLatex(0.17, 0.47,"0.3<#font[12]{k}_{T}<1.0 GeV/#font[12]{c}");
    // KT3:  #LT#font[12]{k}_{T}#GT=0.39 GeV/#font[12]{c}
    // kT:  #LT#font[12]{k}_{T}#GT=0.43 GeV/#font[12]{c}
  }
  //if(KT3Bin==0) {Specif_Kt3 = new TLatex(0.57, 0.83,"0.16<K_{T,3}<0.25 GeV/#font[12]{c}"); Specif_kt = new TLatex(0.57, 0.76,"0.2<k_{T}<0.3 GeV/#font[12]{c}");}
  //if(KT3Bin==1) {Specif_Kt3 = new TLatex(0.57, 0.83,"0.25<K_{T,3}<0.35 GeV/#font[12]{c}"); Specif_kt = new TLatex(0.57, 0.76,"0.3<k_{T}<0.4 GeV/#font[12]{c}");}
  //if(KT3Bin==2) {Specif_Kt3 = new TLatex(0.57, 0.83,"0.35<K_{T,3}<1.0 GeV/#font[12]{c}"); Specif_kt = new TLatex(0.57, 0.76,"0.4<k_{T}<0.5 GeV/#font[12]{c}");}
  Specif_Kt3->SetTextFont(TextFont); Specif_kt->SetTextFont(TextFont);
  Specif_Kt3->SetTextSize(SizeSpecif*SF_2panel); Specif_kt->SetTextSize(SizeSpecif*SF_2panel);
  Specif_Kt3->SetNDC(); Specif_kt->SetNDC();
  Specif_Kt3->Draw("same");
  Specif_kt->Draw("same");
  
  legend4->SetTextFont(TextFont);
  legend4->SetTextSize(SizeLegend*SF_2panel);
  legend4->Draw("same");
  
  TH1D *MarkerPbPb_3=(TH1D*)RadiiPbPb->Clone();
  TH1D *MarkerpPb_3=(TH1D*)RadiipPb->Clone();
  TH1D *Markerpp_3=(TH1D*)Radiipp->Clone();
  TH1D *MarkerPbPb_2=(TH1D*)RadiiC2PbPb->Clone();
  TH1D *MarkerpPb_2=(TH1D*)RadiiC2pPb->Clone();
  TH1D *Markerpp_2=(TH1D*)RadiiC2pp->Clone();
  for(int i=1; i<=MarkerPbPb_3->GetNbinsX(); i++){
    MarkerPbPb_3->SetBinContent(i,1000); MarkerpPb_3->SetBinContent(i,1000); Markerpp_3->SetBinContent(i,1000);
    MarkerPbPb_2->SetBinContent(i,1000); MarkerpPb_2->SetBinContent(i,1000); Markerpp_2->SetBinContent(i,1000);
    MarkerPbPb_3->SetBinError(i,0.001); MarkerpPb_3->SetBinError(i,0.001); Markerpp_3->SetBinError(i,0.001);
    MarkerPbPb_2->SetBinError(i,0.001); MarkerpPb_2->SetBinError(i,0.001); Markerpp_2->SetBinError(i,0.001);
  }
  MarkerPbPb_3->SetBinContent(MarkerPbPb_3->GetXaxis()->FindBin(450), 1);
  MarkerpPb_3->SetBinContent(MarkerPbPb_3->GetXaxis()->FindBin(600), 1);
  Markerpp_3->SetBinContent(MarkerPbPb_3->GetXaxis()->FindBin(800), 1);
  MarkerPbPb_2->SetBinContent(MarkerPbPb_3->GetXaxis()->FindBin(450), 3.1);
  MarkerpPb_2->SetBinContent(MarkerPbPb_3->GetXaxis()->FindBin(600), 3.1);
  Markerpp_2->SetBinContent(MarkerPbPb_3->GetXaxis()->FindBin(800), 3.1);
  if(FitType==0) MarkerpPb_2->SetBinContent(MarkerPbPb_3->GetXaxis()->FindBin(600), 100);
  MarkerPbPb_3->Draw("same"); MarkerpPb_3->Draw("same"); Markerpp_3->Draw("same");
  MarkerPbPb_2->Draw("same"); MarkerpPb_2->Draw("same"); Markerpp_2->Draw("same");

  TLatex *TwoPionText = new TLatex(0.74,0.41,"Two-Pions");
  TLatex *ThreePionText = new TLatex(0.74,0.265,"Three-Pions");
  TwoPionText->SetNDC(); ThreePionText->SetNDC(); 
  TwoPionText->SetTextFont(TextFont); ThreePionText->SetTextFont(TextFont);
  TwoPionText->SetTextSize(SizeSpecif*SF_2panel); ThreePionText->SetTextSize(SizeSpecif*SF_2panel);
  TwoPionText->Draw("same");
  ThreePionText->Draw("same");

  TLatex *Specif_kappas = new TLatex(0.4,0.2,"#kappa_{3}=0.16, #kappa_{4}=0.4");// 0.16 and 0.4
  Specif_kappas->SetNDC();
  Specif_kappas->SetTextFont(TextFont);
  Specif_kappas->SetTextSize(SizeSpecif*SF_2panel);
  if(FitType==1) Specif_kappas->Draw("same");

  ///////////////////////////////////////////////////////////////////
  pad3->cd(2);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.01);
  gPad->SetTopMargin(0.01);
  gPad->SetBottomMargin(0.16);
  gPad->SetTickx(); gPad->SetTicky();
  gPad->SetLogx();
  //gPad->SetGridx(); gPad->SetGridy();
  TH1D *LambdaPbPb=(TH1D*)Parameters_c3[0][FitType][KT3Bin][1]->Clone();
  TH1D *LambdapPb=(TH1D*)Parameters_c3[1][FitType][KT3Bin][1]->Clone();
  TH1D *Lambdapp=(TH1D*)Parameters_c3[2][FitType][KT3Bin][1]->Clone();
  
  LambdaPbPb->GetXaxis()->SetLabelFont(TextFont); LambdaPbPb->GetYaxis()->SetLabelFont(TextFont); 
  LambdaPbPb->GetXaxis()->SetLabelSize(SizeLabel*SF_2panel); LambdaPbPb->GetYaxis()->SetLabelSize(SizeLabel*SF_2panel);
  LambdaPbPb->GetXaxis()->SetNdivisions(808);
  LambdaPbPb->GetYaxis()->SetNdivisions(606);
  LambdaPbPb->GetXaxis()->SetTitleFont(TextFont); LambdaPbPb->GetXaxis()->SetTitleSize(SizeTitle*SF_2panel);
  LambdaPbPb->GetYaxis()->SetTitleFont(TextFont); LambdaPbPb->GetYaxis()->SetTitleSize(SizeTitle*SF_2panel);
  LambdaPbPb->SetMaximum(2.3);
  LambdaPbPb->GetXaxis()->SetTitleOffset(0.95);
  LambdaPbPb->GetYaxis()->SetTitleOffset(100);//1.1
  
  LambdaPbPb->GetXaxis()->SetRangeUser(3,3000);
  LambdaPbPb->Draw();
  
  // old method was PbPb Sys = sqrt(pow(0.05,2)+pow(0.1,2));// MRC + Qrange
  // for pp and pPb: sqrt(pow(0.03,2)+pow(0.05,2));
  for(int cb=0; cb<20; cb++){
    int binPbPb = LambdaPbPb->GetXaxis()->FindBin(meanNchPbPb[cb]);
    int binpPb = LambdapPb->GetXaxis()->FindBin(meanNchpPb[cb]);
    int binpp = Lambdapp->GetXaxis()->FindBin(meanNchpp[cb]);
    double LambdasysPbPb = 0.01*sqrt(pow(SysPercent_PbPb[FitType][3],2) + pow(1,2) + pow(5,2) + pow(10,2)) *  LambdaPbPb->GetBinContent(binPbPb);// fit Variation + MRC + TTC
    double LambdasyspPb = 0.01*sqrt(pow(SysPercent_pPb[FitType][3],2) + pow(1,2) + pow(10,2)) *  LambdapPb->GetBinContent(binpPb);// fit Variation + MRC
    double Lambdasyspp = 0.01*sqrt(pow(SysPercent_pp[FitType][3],2) + pow(1,2) + pow(10,2)) *  Lambdapp->GetBinContent(binpp);// fit Variation + MRC
    if(LambdaPbPb->GetBinContent(binPbPb)==0) {yAxisPbPb[cb]=100; yAxisPbPb_e[cb]=100;}
    else {yAxisPbPb[cb]=LambdaPbPb->GetBinContent(binPbPb); yAxisPbPb_e[cb]=LambdasysPbPb;}
    //
    if(LambdapPb->GetBinContent(binpPb)==0) {yAxispPb[cb]=100; yAxispPb_e[cb]=100;}
    else {yAxispPb[cb]=LambdapPb->GetBinContent(binpPb); yAxispPb_e[cb]=LambdasyspPb;}
    //
    if(Lambdapp->GetBinContent(binpp)==0) {yAxispp[cb]=100; yAxispp_e[cb]=100;}
    else {yAxispp[cb]=Lambdapp->GetBinContent(binpp); yAxispp_e[cb]=Lambdasyspp;}
  }
  TGraphErrors *grLambdaSys_PbPb=new TGraphErrors(20,meanNchPbPb,yAxisPbPb,xAxis_e,yAxisPbPb_e);
  TGraphErrors *grLambdaSys_pPb=new TGraphErrors(20,meanNchpPb,yAxispPb,xAxis_e,yAxispPb_e);
  TGraphErrors *grLambdaSys_pp=new TGraphErrors(20,meanNchpp,yAxispp,xAxis_e,yAxispp_e);
  grLambdaSys_pp->SetMarkerSize(0); grLambdaSys_pp->SetFillStyle(1000); grLambdaSys_pp->SetFillColor(kBlue-10);
  grLambdaSys_pPb->SetMarkerSize(0); grLambdaSys_pPb->SetFillStyle(1000); grLambdaSys_pPb->SetFillColor(kRed-10);
  grLambdaSys_PbPb->SetMarkerSize(0); grLambdaSys_PbPb->SetFillStyle(1000); grLambdaSys_PbPb->SetFillColor(kGray);
  grLambdaSys_pp->SetMarkerColor(kBlue-10); grLambdaSys_pPb->SetMarkerColor(kRed-10); grLambdaSys_PbPb->SetMarkerColor(kGray);
  // C2 sys
  TH1D *LambdaC2PbPb=(TH1D*)Parameters_C2[0][FitType][KT3Bin][1]->Clone();
  TH1D *LambdaC2pPb=(TH1D*)Parameters_C2[1][FitType][KT3Bin][1]->Clone();
  TH1D *LambdaC2pp=(TH1D*)Parameters_C2[2][FitType][KT3Bin][1]->Clone();
  for(int cb=0; cb<20; cb++){
    int binPbPb = LambdaC2PbPb->GetXaxis()->FindBin(meanNchPbPb[cb]);
    int binpPb = LambdaC2pPb->GetXaxis()->FindBin(meanNchpPb[cb]);
    int binpp = LambdaC2pp->GetXaxis()->FindBin(meanNchpp[cb]);
    double LambdasysPbPb = 0.01*sqrt(pow(SysPercent_PbPb[FitType][2],2) + pow(1,2) + pow(5,2)) *  LambdaC2PbPb->GetBinContent(binPbPb);// fit Variation + MRC + TTC
    double LambdasyspPb = 0.01*sqrt(pow(SysPercent_pPb[FitType][2],2) + pow(1,2)) *  LambdaC2pPb->GetBinContent(binpPb);// fit Variation + MRC
    double Lambdasyspp = 0.01*sqrt(pow(SysPercent_pp[FitType][2],2) + pow(1,2)) *  LambdaC2pp->GetBinContent(binpp);// fit Variation + MRC
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
  grLambdaC2Sys_pp->SetMarkerSize(0); grLambdaC2Sys_pp->SetFillStyle(3001); grLambdaC2Sys_pp->SetFillColor(0);
  grLambdaC2Sys_pPb->SetMarkerSize(0); grLambdaC2Sys_pPb->SetFillStyle(3001); grLambdaC2Sys_pPb->SetFillColor(0);
  grLambdaC2Sys_PbPb->SetMarkerSize(0); grLambdaC2Sys_PbPb->SetFillStyle(3001); grLambdaC2Sys_PbPb->SetFillColor(0);
  //
  
  if(FitType==1){
    //grLambdaC2Sys_pp->Draw("|| p");
    //grLambdaC2Sys_pPb->Draw("|| p");
  }
  grLambdaC2Sys_PbPb->Draw("|| p");

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
 
  
  
  /*TLatex *Specif_Marker_2;
  if(FitType==0) Specif_Marker_2 = new TLatex(0.55,0.9,"Hollow=#lambda^{#font[12]{G}}  Solid=#lambda^{#font[12]{G}}_{3}");
  else Specif_Marker_2 = new TLatex(0.3,0.9,"Hollow=#lambda^{#font[12]{E_{w}}}, #font[12]{R^{E_{w}}}_{inv}   Solid=#lambda^{#font[12]{E_{w}}}_{3}, #font[12]{R^{E_{w}}}_{inv,3}");
  Specif_Marker_2->SetNDC();
  Specif_Marker_2->SetTextFont(TextFont);
  Specif_Marker_2->SetTextSize(SizeSpecif*SF_2panel);
  Specif_Marker_2->Draw("same");*/
  
  can3->cd();
  TPad *pad3_2 = new TPad("pad3_2","pad3_2",0.0,0.0,1.,1.);
  pad3_2->SetFillStyle(0);
  pad3_2->Draw();
  pad3_2->cd();
  TLatex *RinvTitle;
  if(FitType==0) RinvTitle=new TLatex(0.062,0.72,"#font[12]{R}^{#font[12]{G}}_{inv} or #font[12]{R}^{#font[12]{G}}_{inv,3} (fm)");
  else RinvTitle=new TLatex(0.062,0.72,"#font[12]{R}^{#font[12]{E_{w}}}_{inv} or #font[12]{R}^{#font[12]{E_{w}}}_{inv,3} (fm)");
  RinvTitle->SetNDC();
  RinvTitle->SetTextFont(TextFont);
  RinvTitle->SetTextSize(SizeTitle);
  RinvTitle->SetTextAngle(90);
  RinvTitle->Draw("same");
  TLatex *LambdaTitle;
  if(FitType==0) LambdaTitle=new TLatex(0.064,0.35,"#lambda^{#font[12]{G}} or #lambda^{#font[12]{G}}_{3}");
  else LambdaTitle=new TLatex(0.064,0.33,"#lambda^{#font[12]{E_{w}}} or #lambda^{#font[12]{E_{w}}}_{3}");
  LambdaTitle->SetNDC();
  LambdaTitle->SetTextFont(TextFont);
  LambdaTitle->SetTextSize(SizeTitle);
  LambdaTitle->SetTextAngle(90);
  LambdaTitle->Draw("same");
  

  if(SaveFiles && FitType==0) can3->SaveAs("ThreePionFitParametersGauss.eps");
  if(SaveFiles && FitType==1) can3->SaveAs("ThreePionFitParametersEW.eps");
  
  

  if(KT3Bin>0) {cout<<"Only print Radii/Lambda Plot for this setting"<<endl; return;}

  ////////////////////////////////////////////////////
  ////////////////////////////////////////////////////
  // Correlation functions and Monte Carlo
  /*
  TCanvas *can4 = new TCanvas("can4", "can4",10,0,900,600);// 11,53,700,500
  gStyle->SetOptFit(0111);
  can4->SetFillColor(0);//10
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
  legend5[0] = new TLegend(.3,.52, .97,.99,NULL,"brNDC");//.45 or .4 for x1
  legend5[0]->SetBorderSize(0);
  legend5[0]->SetFillColor(0);
  legend5[0]->SetTextFont(TextFont);
  legend5[1]=(TLegend*)legend5[0]->Clone();
  legend5[2]=(TLegend*)legend5[0]->Clone();
  legend5[3]=(TLegend*)legend5[0]->Clone();
  legend5[4]=(TLegend*)legend5[0]->Clone();
  legend5[5]=(TLegend*)legend5[0]->Clone();
  TLegend *legendFitTypes = (TLegend*)legend5[0]->Clone();
  
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
    //gPad->SetTickx(); gPad->SetTicky();
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
    double SF_correction=1.0;
    if(padNum==3) SF_correction=1.1;
    if(padNum==4) SF_correction=0.8;
    if(padNum==5 || padNum==6) SF_correction=0.92;
    
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
    
    C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->SetMinimum(0.9); 
    C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->SetMaximum(3.3);
    //
    //c3[2][0][0][KT3Bin][Mb_pp]->GetXaxis()->SetLabelSize(SizeLabel*SF2); c3[2][0][0][KT3Bin][Mb_pp]->GetYaxis()->SetLabelSize(SizeLabel*SF2);
    //c3[2][0][0][KT3Bin][Mb_pp]->GetXaxis()->SetNdivisions(808);
    //c3[2][0][0][KT3Bin][Mb_pp]->GetYaxis()->SetTitleSize(SizeTitle*SF1);
    C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetYaxis()->SetTitle("#font[12]{C}_{3} or #font[12]{#bf{c}}_{3} ");
    if(padNum<=5){
      C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetXaxis()->SetTitleOffset(10); 
    }else C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetXaxis()->SetTitleOffset(0.88);
    if(padNum==1) C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetYaxis()->SetTitleOffset(0.55);
    else C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetYaxis()->SetTitleOffset(10.);
    if(padNum>=5) C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetXaxis()->SetLabelOffset(-.0);

    //C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetXaxis()->SetLabelSize(SizeLabel);
    //C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetYaxis()->SetLabelSize(SizeLabel);
    //C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetXaxis()->SetTitleSize(SizeTitle);
    //C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetYaxis()->SetTitleSize(SizeTitle);
    C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetXaxis()->SetNdivisions(504);
    C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetYaxis()->SetNdivisions(504);
    C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetXaxis()->SetTitleSize(SizeTitle*SF_6pannel*SF_correction);
    C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetYaxis()->SetTitleSize(SizeTitle*SF_6pannel*SF_correction);
    C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetXaxis()->SetLabelSize(SizeTitle*SF_6pannel*SF_correction);
    C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof]->GetYaxis()->SetLabelSize(SizeTitle*SF_6pannel*SF_correction);
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
      //if(System_proof==0) legend5[padNum-1]->AddEntry(HIJING_c3_SC,"HIJING #font[12]{#bf{c}}_{3}^{#pm#pm#pm}","p");
      //if(System_proof==1) legend5[padNum-1]->AddEntry(c3[System_proof][1][ChComb_proof][KT3Bin][Mb_proof],"DPMJET #font[12]{#bf{c}}_{3}^{#pm#pm#pm}","p");
      //if(System_proof==2) legend5[padNum-1]->AddEntry(c3[System_proof][1][ChComb_proof][KT3Bin][Mb_proof],"PYTHIA #font[12]{#bf{c}}_{3}^{#pm#pm#pm}","p");
    }else{
      legend5[padNum-1]->AddEntry(C3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof],"#font[12]{C}_{3}^{#pm#pm#mp}","p");
      legend5[padNum-1]->AddEntry(c3[System_proof][0][ChComb_proof][KT3Bin][Mb_proof],"#font[12]{#bf{c}}_{3}^{#pm#pm#mp}","p");
      if(System_proof==0) legend5[padNum-1]->AddEntry(HIJING_c3_MC,"HIJING #font[12]{#bf{c}}_{3}^{#pm#pm#mp}","p");
      if(System_proof==1) legend5[padNum-1]->AddEntry(c3[System_proof][1][ChComb_proof][KT3Bin][Mb_proof],"DPMJET #font[12]{#bf{c}}_{3}^{#pm#pm#mp}","p");
      if(System_proof==2) legend5[padNum-1]->AddEntry(c3[System_proof][1][ChComb_proof][KT3Bin][Mb_proof],"PYTHIA #font[12]{#bf{c}}_{3}^{#pm#pm#mp}","p");
    }
   
    if(ChComb_proof==0) {
      c3_fit[System_proof][0][KT3Bin][Mb_proof]->Draw("c same");
      c3_fit[System_proof][1][KT3Bin][Mb_proof]->Draw("c same");
      if(padNum==3){
	legendFitTypes->AddEntry(c3_fit[System_proof][0][KT3Bin][Mb_proof],"Gaussian","l");
	legendFitTypes->AddEntry(c3_fit[System_proof][1][KT3Bin][Mb_proof],"Edgeworth","l");
      }
    }
    
    TLatex *CTLabel;
    if(System_proof==0) CTLabel = new TLatex(0.12,0.9,"Pb-Pb #sqrt{#font[12]{s}_{NN}}=2.76 TeV");// 0.003,.988
    if(System_proof==1) CTLabel = new TLatex(0.15,0.9,"p-Pb #sqrt{#font[12]{s}_{NN}}=5.02 TeV");// 0.003,.988
    if(System_proof==2) CTLabel = new TLatex(0.4,0.9,"pp #sqrt{s}=7 TeV");// 0.003,.988
    CTLabel->SetNDC();
    CTLabel->SetTextFont(TextFont);
    CTLabel->SetTextSize(SizeSpecif*SF_6pannel*SF_correction);
    if(padNum>=4) CTLabel->Draw("same");
    
    TString *nameCh=new TString("#LTN_{ch}#GT = ");
    float Nch=1;
    //if(System_proof==0) Nch = pow(10,meanNchPbPb[Mb_proof]);
    //else if(System_proof==1) Nch = pow(10,meanNchpPb[Mb_proof]);
    //else Nch = pow(10,meanNchpp[Mb_proof]);
    if(System_proof==0) Nch = meanNchPbPb[Mb_proof];
    else if(System_proof==1) Nch = meanNchpPb[Mb_proof];
    else Nch = meanNchpp[Mb_proof];
    *nameCh += int(Nch + 0.5);
    nameCh->Append(" #pm ");
    float SysPercent = 0.05;
    int SigFig=0;
    if(SysPercent*Nch < 1) {nameCh->Append("0."); SigFig=int(10*SysPercent*Nch + 0.5);}
    else {SigFig=int(SysPercent*Nch + 0.5);}
    
    *nameCh += SigFig;
    TLatex *MLabel;
    if(padNum==1) MLabel = new TLatex(0.52,0.6,nameCh->Data());
    if(padNum==2) MLabel = new TLatex(0.4,0.6,nameCh->Data());
    if(padNum==3) MLabel = new TLatex(0.4,0.6,nameCh->Data());
   
    MLabel->SetNDC();
    MLabel->SetTextFont(TextFont);
    MLabel->SetTextSize(SizeSpecif*SF_6pannel*SF_correction);
    if(padNum<=3) MLabel->Draw("same");
    
    
    legend5[padNum-1]->SetTextSize(SizeLegend*SF_6pannel*SF_correction);
    if(padNum==1) {legend5[padNum-1]->SetX1(0.45); legend5[padNum-1]->SetY1(0.69);}
    if(padNum==2) {legend5[padNum-1]->SetX1(0.32); legend5[padNum-1]->SetY1(0.69);}
    if(padNum==3) {
      legend5[padNum-1]->SetX1(0.32); legend5[padNum-1]->SetY1(0.69);
      legendFitTypes->SetX1(0.44); 
      legendFitTypes->SetY1(0.22); legendFitTypes->SetY2(0.56);
      legendFitTypes->Draw("same");
    }
    if(padNum==4) {legend5[padNum-1]->SetX1(0.45); legend5[padNum-1]->SetY1(0.45); legend5[padNum-1]->SetY2(0.85);}
    if(padNum==5) {legend5[padNum-1]->SetX1(0.32); legend5[padNum-1]->SetY1(0.45); legend5[padNum-1]->SetY2(0.85);}
    if(padNum==6) {legend5[padNum-1]->SetX1(0.32); legend5[padNum-1]->SetY1(0.45); legend5[padNum-1]->SetY2(0.85);}
    legend5[padNum-1]->Draw("same");
  }
  
  

  can4->cd();
  
  TPad *pad4_2 = new TPad("pad4_2","pad4_2",0.0,0.0,1.,1.);
  pad4_2->SetFillStyle(0);
  pad4_2->Draw();
  pad4_2->cd();
  TBox *CoverUp1 = new TBox(0.35,0.05,0.42,.115);
  CoverUp1->SetFillColor(0);
  CoverUp1->Draw();
  TBox *CoverUp2 = new TBox(0.66,0.05,0.73,.115);
  CoverUp2->SetFillColor(0);
  CoverUp2->Draw();
  
  

 
  
  ////////////////////////////////////////////////////
  ////////////////////////////////////////////////////
  // 2 system correlation function comparison
  TCanvas *can5 = new TCanvas("can5", "can5",1700,700,600,600);// 11,53,700,500
  can5->SetHighLightColor(2);
  gStyle->SetOptFit(0);// 0111 to show fit stat box
  can5->SetFillColor(0);//10
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
  TLegend *legend6 = new TLegend(.4,.6, .9,.95,NULL,"brNDC");//.45 or .4 for x1
  legend6->SetBorderSize(0);
  legend6->SetFillColor(0);
  legend6->SetTextFont(TextFont);
  legend6->SetTextSize(SizeLegend);
  TLegend *legend7 = new TLegend(.65,.35, .97,.52,NULL,"brNDC");//.45 or .4 for x1
  legend7->SetBorderSize(0);
  legend7->SetFillColor(0);
  legend7->SetTextFont(TextFont);
  legend7->SetTextSize(SizeLegend);
  //
  gPad->SetRightMargin(0.01); gPad->SetLeftMargin(0.10);
  gPad->SetBottomMargin(0.14); gPad->SetTopMargin(0.02);
  gPad->SetTickx(); gPad->SetTicky();
  //
  int KT3Bin_CorrComp=0;
  int Mbin_SysComp_PbPb=12;// 12
  int Mbin_SysComp_pPb;// 12, 16,
  int Mbin_SysComp_pp=15;//   ,  15,
  if(p_pPb_Comp) Mbin_SysComp_pPb=16; 
  else Mbin_SysComp_pPb=12;

  TH1D *c3_PbPb=(TH1D*)c3[0][0][0][KT3Bin_CorrComp][Mbin_SysComp_PbPb]->Clone();
  TH1D *c3_pPb=(TH1D*)c3[1][0][0][KT3Bin_CorrComp][Mbin_SysComp_pPb]->Clone();
  TH1D *c3_pp=(TH1D*)c3[2][0][0][KT3Bin_CorrComp][Mbin_SysComp_pp]->Clone();

  c3_pPb->GetYaxis()->SetTitle("#bf{c}_{3}^{#pm#pm#pm}");
  c3_pPb->GetXaxis()->SetRangeUser(0,0.27);
  
  c3_pPb->GetXaxis()->SetLabelSize(SizeLabel); c3_pPb->GetYaxis()->SetLabelSize(SizeLabel);
  c3_pPb->GetXaxis()->SetTitleSize(SizeTitle); c3_pPb->GetYaxis()->SetTitleSize(SizeTitle);
  
  c3_pPb->GetYaxis()->SetTitle("#font[12]{#bf{c}}_{3}^{#pm#pm#pm}");
  c3_pPb->GetXaxis()->SetTitleOffset(1.0); 
  c3_pPb->GetYaxis()->SetTitleOffset(0.7); 
  c3_pPb->SetMinimum(0.9); c3_pPb->SetMaximum(3.3);
  c3_pPb->GetXaxis()->SetNdivisions(503);
  c3_pPb->GetYaxis()->SetNdivisions(503);

  c3_pPb->Draw();
  //
  if(p_pPb_Comp) c3_Sys[2][0][0][KT3Bin_CorrComp][Mbin_SysComp_pp]->Draw("E2 same");
  c3_Sys[1][0][0][KT3Bin_CorrComp][Mbin_SysComp_pPb]->Draw("E2 same");
  if(!p_pPb_Comp) c3_Sys[0][0][0][KT3Bin][Mbin_SysComp_PbPb]->Draw("E2 same");
  if(p_pPb_Comp) c3_pp->Draw("same");
  c3_pPb->Draw("same");
  if(!p_pPb_Comp) c3_PbPb->Draw("same");
  // SysPercent = 0.07 - 0.05*(Mb_proof/19.);
  
  if(p_pPb_Comp) {
    legend6->AddEntry(c3_pPb,"#splitline{p-Pb #sqrt{#font[12]{s}_{NN}}=5.02 TeV}{#LTN_{ch}#GT = 23 #pm 1.2}","p");// MpPb=16
    legend6->AddEntry(c3_pp,"#splitline{pp #sqrt{s}=7 TeV}{#LTN_{ch}#GT = 26 #pm 1.3}","p");// Mpp=15
  }else{  
    legend6->AddEntry(c3_pPb,"#splitline{p-Pb #sqrt{#font[12]{s}_{NN}}=5.02 TeV}{#LTN_{ch}#GT = 71 #pm 4}","p");// MpPb=12
    legend6->AddEntry(c3_PbPb,"#splitline{Pb-Pb #sqrt{#font[12]{s}_{NN}}=2.76 TeV}{#LTN_{ch}#GT = 87 #pm 4}","p");// MPbPb=12
  }
  
  //
  
  if(!p_pPb_Comp) c3_fit[0][0][KT3Bin][Mbin_SysComp_PbPb]->Draw("c same");
  c3_fit[1][0][KT3Bin][Mbin_SysComp_pPb]->Draw("c same");
  if(p_pPb_Comp) c3_fit[2][0][KT3Bin][Mbin_SysComp_pp]->Draw("c same");
  //
  TF1 *GaussFit_forLegend=(TF1*)c3_fit[1][0][KT3Bin][Mbin_SysComp_pPb]->Clone();
  TF1 *EwFit_forLegend=(TF1*)c3_fit[1][1][KT3Bin][Mbin_SysComp_pPb]->Clone();
  GaussFit_forLegend->SetLineColor(1);
  EwFit_forLegend->SetLineColor(1);
  legend7->AddEntry(GaussFit_forLegend,"Gaussian","l");
  legend7->AddEntry(EwFit_forLegend,"Edgeworth","l");
  //
  
  if(!p_pPb_Comp) c3_fit[0][1][KT3Bin][Mbin_SysComp_PbPb]->Draw("c same");
  c3_fit[1][1][KT3Bin][Mbin_SysComp_pPb]->Draw("c same");
  if(p_pPb_Comp) c3_fit[2][1][KT3Bin][Mbin_SysComp_pp]->Draw("c same");
  legend6->Draw("same");
  legend7->Draw("same");
  
  TLatex *Specif_Kt3_4 = new TLatex(0.52,0.55,"0.16<#font[12]{K}_{T,3}<0.3 GeV/#font[12]{c}");
  Specif_Kt3_4->SetNDC();
  Specif_Kt3_4->SetTextFont(TextFont);
  Specif_Kt3_4->SetTextSize(SizeSpecif);
  Specif_Kt3_4->Draw("same");
  

  //
  if(SaveFiles) {
    TString *name = new TString("c3SystemComp");
    if(FitType==1) name->Append("EW");
    name->Append("_MpPb");
    *name += Mbin_SysComp_pPb;
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
  can6->SetFillColor(0);//10
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
  //legend7->AddEntry(C2[1][Mbin_C2],"N_{ch} = 10.2 p-Pb #sqrt{#font[12]{s}_{NN}}=5.02 TeV","p");
  legend7->AddEntry(C2[2][Mbin_C2],"N_{ch} = 4.6 pp #sqrt{s}=7 TeV","p");
  legend7->AddEntry(C2[1][Mbin_C2],"N_{ch} = 5.4 p-Pb #sqrt{#font[12]{s}_{NN}}=5.02 TeV","p");
  //
  legend7->Draw("same");
  TLatex *Specif_Kt = new TLatex(0.22,1.35,"0.16<K_{t}<1.0 GeV/#font[12]{c}");
  Specif_Kt->SetTextFont(TextFont);
  Specif_Kt->SetTextSize(SizeSpecif*SF2);
  Specif_Kt->Draw("same");
  TLatex *Specif_Kinematics_5 = new TLatex(0.172,1.3,"|#eta|<0.8, 0.16<p_{T}<1.0 GeV/#font[12]{c}");
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
  legend6->AddEntry(Parameters_c3[1][2],"p-Pb #sqrt{#font[12]{s}_{NN}}=5.02 TeV","p");
  legend6->AddEntry(Parameters_c3[0][2],"Pb-Pb #sqrt{#font[12]{s}_{NN}}=2.76 TeV","p");
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
    canvas->SetLeftMargin(0.14); 
    canvas->SetBottomMargin(0.18); 
    
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
