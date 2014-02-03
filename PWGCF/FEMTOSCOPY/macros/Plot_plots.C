#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <Riostream.h>

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

#define BohrR 1963.6885
#define FmToGeV 0.19733 // conversion to fm
#define PI 3.1415926
#define masspiC 0.1395702 // pi+ mass (GeV/c^2)

using namespace std;

bool TherminatorC2=kFALSE;
const int BOI_1=0;// centrality bin (0-9)
const int BOI_2=9;// centrality of second bin for C2 fit parameter plot only
const int ChProdBOI=0;// 0=SameCharge, 1=MixedCharge
const int KT3Bin=0;// Kt3 bin. 0=low Kt3 bin.  1=high Kt3 bin
//
const int CoulChoice=0;// 0 for GRS (default), 1 for Omega0

double MIN[10]={0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97};// C2 y-axis min
double MAX[10]={1.39, 1.38, 1.38, 1.38, 1.38, 1.38, 1.38, 1.38, 1.38, 1.53};// C2 y-axis max
//
const int KTBINS = 6;
int KTINDEX;
bool ChargeConstraint=kFALSE;
bool LinkRadii=kFALSE;
//
int TextFont=42;// 63, or 42
float SizeLabel=0.1;// 20(63 font), 0.08(42 font)
float SizeLegend=0.1;// .08
float SizeTitle=0.12;// 
float SizeSpecif=0.075;// 
float SF1=2/3.*0.95;
float SF2=1/2.*0.95;

double RightMargin=0.004;// 0.002
//
double Chi2_C2global;
double NFitPoints_C2global;
TH1D *C2_ss[KTBINS][10];
TH1D *C2_os[KTBINS][10];
TH1D *C2_ss_MomSys[KTBINS][10];
TH1D *C2_ss_KSys[KTBINS][10];
TH1D *C2_os_MomSys[KTBINS][10];
TH1D *C2_os_KSys[KTBINS][10];
TH1D *hfitC2ss_noG[KTBINS][10];
TH1D *hfitC2os_noG[KTBINS][10];
TF1 *fitC2ss_noG[KTBINS][10];
TF1 *fitC2os_noG[KTBINS][10];
TF1 *fitC2ss_yesG[KTBINS][10];
TF1 *fitC2os_yesG[KTBINS][10];
TF1 *fitC2ss_noGEWfromTherm[KTBINS][10];
TF1 *fitC2os_noGEWfromTherm[KTBINS][10];
TF1 *fitC2ss_yesGEWfromTherm[KTBINS][10];
TF1 *fitC2os_yesGEWfromTherm[KTBINS][10];
TF1 *fitC2ss_noGEW[KTBINS][10];
TF1 *fitC2os_noGEW[KTBINS][10];
TH1D *hfitC2ss_noGEW[KTBINS][10];
TH1D *hfitC2os_noGEW[KTBINS][10];
TH1D *K2_ss[KTBINS][10];
TH1D *K2_os[KTBINS][10];
TH1D *C2Therm_ss[KTBINS][10];
TH1D *C2Therm_os[KTBINS][10];
TF1 *fitC2ss_Therm[KTBINS][10];
TF1 *fitC2os_Therm[KTBINS][10];
TF1 *fitC2ss_ThermGaus[KTBINS][10];
TF1 *fitC2os_ThermGaus[KTBINS][10];


void DrawALICELogo(Bool_t, Float_t, Float_t, Float_t, Float_t);

void Plot_plots(){

  gStyle->SetOptStat(0);
  gStyle->SetOptDate(0);
  

  TFile *files_2_noG[KTBINS][10];
  TFile *files_2_yesG[KTBINS][10];
  TFile *files_2_noGEWfromTherm[KTBINS][10];
  TFile *files_2_yesGEWfromTherm[KTBINS][10];
  TFile *files_2_noGEW[KTBINS][10];
  TFile *files_3[2][2][2][2][10];// SC/MC, +/-, GRS/Omega0, Therm/Gauss, MBINS
  //
  double intercept1[10]={0};
  double intercept1_e[10]={0};
  double intercept2[10]={0};
  double intercept2_e[10]={0};
  double ChiNDF1[10]={0};
  double ChiNDF2[10]={0};
  double cent_mean[10]={0};
  double cent_mean_e[10]={0};

  
  TH1D *Chi2NDF_yesG[KTBINS][10];
  TH1D *Chi2NDF_yesGEWfromTherm[KTBINS][10];
  TH1D *Chi2NDF_noG[KTBINS][10];
  TH1D *Chi2NDF_noGEWfromTherm[KTBINS][10];
  TH1D *Cterm1;
  TH1D *C3[2][2][2][10];// SC/MC, GRS/Omega0, Therm/Gauss, cb
  TH1D *c3[2][2][2][10];// SC/MC, GRS/Omega0, Therm/Gauss, cb
  TH1D *r3_Q3[2][2][10];// GRS/Omega0, Therm/Gauss, cb
  TH1D *C3EW[2][2][10];//  GRS/Omega0, Therm/Gauss, cb
  TH1D *K3[2][2][2][10];// SC/MC, GRS/Omega0, Therm/Gauss, cb
  //
  TH1D *ParHisto_ch[4][10];
  TH1D *ParHisto_coh[4][10];
  TH1D *ParHisto_chEWfromTherm[4][10];
  TH1D *ParHisto_cohEWfromTherm[4][10];
  TH1D *ParHisto_chEW[4][10];

  for(int ii=0; ii<10; ii++){
    for(int par=1; par<=4; par++){
      TString *name_ch = new TString("ParHisto_ch");
      *name_ch += ii; *name_ch += par;
      TString *name_coh = new TString("ParHisto_coh");
      *name_coh += ii; *name_coh += par;
      TString *name_EWfromThermch = new TString("ParHisto_EWfromThermch");
      *name_EWfromThermch += ii; *name_EWfromThermch += par;
      TString *name_EWch = new TString("ParHisto_EWch");
      *name_EWch += ii; *name_EWch += par;
      TString *name_EWcoh = new TString("ParHisto_EWcoh");
      *name_EWcoh += ii; *name_EWcoh += par;
      
      int MStyle_ch, MStyle_coh;
      if(ii==BOI_1) {MStyle_ch = 20; MStyle_coh = 24;}
      else {MStyle_ch = 22; MStyle_coh = 26;}
      
      ParHisto_ch[par-1][ii] = new TH1D(name_ch->Data(),"",10,0,1);
      ParHisto_ch[par-1][ii]->SetMarkerStyle(MStyle_ch);
      ParHisto_ch[par-1][ii]->SetMarkerColor(1);
      ParHisto_ch[par-1][ii]->SetLineColor(1);
      ParHisto_ch[par-1][ii]->SetMarkerSize(1.5);
      ParHisto_coh[par-1][ii] = new TH1D(name_coh->Data(),"",10,0,1);
      ParHisto_coh[par-1][ii]->SetMarkerStyle(MStyle_coh);
      ParHisto_coh[par-1][ii]->SetMarkerColor(4);
      ParHisto_coh[par-1][ii]->SetLineColor(4);
      ParHisto_coh[par-1][ii]->SetMarkerSize(1.5);
      ParHisto_chEWfromTherm[par-1][ii] = new TH1D(name_EWfromThermch->Data(),"",10,0,1);
      ParHisto_chEWfromTherm[par-1][ii]->SetMarkerStyle(MStyle_ch);
      ParHisto_chEWfromTherm[par-1][ii]->SetMarkerColor(1);
      ParHisto_chEWfromTherm[par-1][ii]->SetLineColor(1);
      ParHisto_chEWfromTherm[par-1][ii]->SetMarkerSize(1.5);
      ParHisto_chEW[par-1][ii] = new TH1D(name_EWch->Data(),"",10,0,1);
      ParHisto_chEW[par-1][ii]->SetMarkerStyle(MStyle_coh);
      ParHisto_chEW[par-1][ii]->SetMarkerColor(4);
      ParHisto_chEW[par-1][ii]->SetLineColor(4);
      ParHisto_chEW[par-1][ii]->SetMarkerSize(1.5);
      ParHisto_cohEWfromTherm[par-1][ii] = new TH1D(name_EWcoh->Data(),"",10,0,1);
      ParHisto_cohEWfromTherm[par-1][ii]->SetMarkerStyle(MStyle_coh);
      ParHisto_cohEWfromTherm[par-1][ii]->SetMarkerColor(4);
      ParHisto_cohEWfromTherm[par-1][ii]->SetLineColor(4);
      ParHisto_cohEWfromTherm[par-1][ii]->SetMarkerSize(1.5);
    }
  }
  //////////////////////////////

  // Start File access
  for(int cb=0; cb<10; cb++){
    for(int ChComb=0; ChComb<2; ChComb++) {// SC or MC
      for(int ch=0; ch<2; ch++) {// - or +
	for(int KT3=0; KT3<2; KT3++) {// Kt3 bin
	  TString *name3 = new TString("OutFiles/OutFile");
	  if(ChComb==0) name3->Append("SC");
	  else name3->Append("MC");
	  if(ch==0) name3->Append("Neg");
	  else name3->Append("Pos");
	  TString *name3_Omega0=new TString(name3->Data());
	  name3->Append("NoGEWGRS");
	  name3_Omega0->Append("NoGEWOmega0");
	  name3->Append("Kt3_"); name3_Omega0->Append("Kt3_"); 
	  *name3 += KT3+1; *name3_Omega0 += KT3+1;
	  
	  name3->Append("_Kt10_M");
	  name3_Omega0->Append("_Kt10_M");
	  if(cb<10) {*name3 += cb; *name3_Omega0 += cb;}
	  else {*name3 += 0; *name3_Omega0 += 0;}
	  name3->Append(".root");
	  name3_Omega0->Append(".root");
	  files_3[ChComb][ch][0][KT3][cb] = new TFile(name3->Data(),"READ");
	  files_3[ChComb][ch][1][KT3][cb] = new TFile(name3_Omega0->Data(),"READ");
	  ///////////////////////////////
	  for(int Coul=0; Coul<2; Coul++){
	    if(ch==0) {
	      C3[ChComb][Coul][KT3][cb]=(TH1D*)files_3[ChComb][ch][Coul][KT3][cb]->Get("C3");
	      C3[ChComb][Coul][KT3][cb]->SetDirectory(0);
	      c3[ChComb][Coul][KT3][cb]=(TH1D*)files_3[ChComb][ch][Coul][KT3][cb]->Get("c3");
	      c3[ChComb][Coul][KT3][cb]->SetDirectory(0);
	      if(ChComb==0) C3EW[Coul][KT3][cb]=(TH1D*)files_3[ChComb][ch][Coul][KT3][cb]->Get("C3_EWexpectation");
	      if(ChComb==0) r3_Q3[Coul][KT3][cb]=(TH1D*)files_3[ChComb][ch][Coul][KT3][cb]->Get("r3_Q3");
	      C3EW[Coul][KT3][cb]->SetDirectory(0);
	      r3_Q3[Coul][KT3][cb]->SetDirectory(0);
	      
	      if(Coul==1){
		K3[ChComb][0][KT3][cb]=(TH1D*)files_3[ChComb][ch][Coul][KT3][cb]->Get("Coul_GRiverside");
		K3[ChComb][0][KT3][cb]->SetDirectory(0);
		K3[ChComb][1][KT3][cb]=(TH1D*)files_3[ChComb][ch][Coul][KT3][cb]->Get("Coul_Omega0");
		K3[ChComb][1][KT3][cb]->SetDirectory(0);
	      }
	    }else{
	      TH1D *tempC3=(TH1D*)files_3[ChComb][ch][Coul][KT3][cb]->Get("C3");
	      TH1D *tempc3=(TH1D*)files_3[ChComb][ch][Coul][KT3][cb]->Get("c3");
	      TH1D *tempC3EW=(TH1D*)files_3[ChComb][ch][Coul][KT3][cb]->Get("C3_EWexpectation");
	      TH1D *tempr3_Q3;
	      if(ChComb==0) tempr3_Q3=(TH1D*)files_3[ChComb][ch][Coul][KT3][cb]->Get("r3_Q3");
	      for(int bin=1; bin<=tempC3->GetNbinsX(); bin++){
		double valueC3 = (C3[ChComb][Coul][KT3][cb]->GetBinContent(bin) + tempC3->GetBinContent(bin))/2.;
		double valueC3_e = sqrt(pow(C3[ChComb][Coul][KT3][cb]->GetBinError(bin),2)+pow(tempC3->GetBinError(bin),2))/sqrt(2.);
		double valuec3 = (c3[ChComb][Coul][KT3][cb]->GetBinContent(bin) + tempc3->GetBinContent(bin))/2.;
		double valuec3_e = sqrt(pow(c3[ChComb][Coul][KT3][cb]->GetBinError(bin),2)+pow(tempc3->GetBinError(bin),2))/sqrt(2.);
		C3[ChComb][Coul][KT3][cb]->SetBinContent(bin, valueC3);
		C3[ChComb][Coul][KT3][cb]->SetBinError(bin, valueC3_e);
		c3[ChComb][Coul][KT3][cb]->SetBinContent(bin, valuec3);
		c3[ChComb][Coul][KT3][cb]->SetBinError(bin, valuec3_e);
		if(Coul==0) {c3[ChComb][Coul][KT3][cb]->SetMarkerStyle(20); c3[ChComb][Coul][KT3][cb]->SetMarkerColor(4); c3[ChComb][Coul][KT3][cb]->SetLineColor(4);}
		else {c3[ChComb][Coul][KT3][cb]->SetMarkerStyle(22); c3[ChComb][Coul][KT3][cb]->SetMarkerColor(2); c3[ChComb][Coul][KT3][cb]->SetLineColor(2);}
		c3[ChComb][Coul][KT3][cb]->GetYaxis()->SetTitleOffset(1.2);
		if(ChComb==1){
		  c3[ChComb][Coul][KT3][cb]->SetMinimum(0.96);// 0.99
		  c3[ChComb][Coul][KT3][cb]->SetMaximum(1.02);// 1.02
		  c3[ChComb][Coul][KT3][cb]->GetXaxis()->SetRangeUser(0,0.13);
		  c3[ChComb][Coul][KT3][cb]->GetYaxis()->SetTitleOffset(1.2);
		  c3[ChComb][Coul][KT3][cb]->GetYaxis()->SetTitleSize(.05);
		}
		
		if(ChComb==0) {// only same-charge for r3
		  double valuer3_Q3 = (r3_Q3[Coul][KT3][cb]->GetBinContent(bin) + tempr3_Q3->GetBinContent(bin))/2.;
		  double valuer3_Q3_e = sqrt(pow(r3_Q3[Coul][KT3][cb]->GetBinError(bin),2) + pow(tempr3_Q3->GetBinError(bin),2))/sqrt(2.);
		  double valueC3EW = (C3EW[Coul][KT3][cb]->GetBinContent(bin) + tempC3EW->GetBinContent(bin))/2.;
		  double valueC3EW_e = sqrt(pow(C3EW[Coul][KT3][cb]->GetBinError(bin),2)+pow(tempC3EW->GetBinError(bin),2))/sqrt(2.);
		  C3EW[Coul][KT3][cb]->SetBinContent(bin, valueC3EW);
		  C3EW[Coul][KT3][cb]->SetBinError(bin, valueC3EW_e);
		  r3_Q3[Coul][KT3][cb]->SetBinContent(bin, valuer3_Q3);
		  r3_Q3[Coul][KT3][cb]->SetBinError(bin, valuer3_Q3_e);
		  r3_Q3[Coul][KT3][cb]->GetXaxis()->SetRangeUser(0,0.14);
		  if(cb<2) r3_Q3[Coul][KT3][cb]->SetMinimum(-3);
		  else r3_Q3[Coul][KT3][cb]->SetMinimum(-0.5);
		  r3_Q3[Coul][KT3][cb]->SetMaximum(2.6);
		  if(Coul==0) {r3_Q3[Coul][KT3][cb]->SetMarkerStyle(20); r3_Q3[Coul][KT3][cb]->SetMarkerColor(4); r3_Q3[Coul][KT3][cb]->SetLineColor(4);}
		  else {r3_Q3[Coul][KT3][cb]->SetMarkerStyle(22); r3_Q3[Coul][KT3][cb]->SetMarkerColor(2); r3_Q3[Coul][KT3][cb]->SetLineColor(2);}
		}
	      }
	    }
	    files_3[ChComb][ch][Coul][KT3][cb]->Close();
	  }// Coul
	}
      }
    }
  
  
    for(int kt=0; kt<KTBINS; kt++){
      
      TString *name = new TString("OutFiles/OutFileSCPosNoG"); 
      TString *nameEWfromTherm = new TString(name->Data()); nameEWfromTherm->Append("EWfromThermGRS");
      TString *nameEW = new TString(name->Data()); nameEW->Append("EWGRS");
      name->Append("GRS");
      name->Append("Kt3_1"); nameEWfromTherm->Append("Kt3_1"); nameEW->Append("Kt3_1");
      name->Append("_Kt"); nameEWfromTherm->Append("_Kt"); nameEW->Append("_Kt");
      *name += kt+1; *nameEWfromTherm += kt+1; *nameEW += kt+1;
      name->Append("_M"); nameEWfromTherm->Append("_M"); nameEW->Append("_M");
      if(cb<10) {*name += cb; *nameEWfromTherm += cb; *nameEW += cb;}
      else {*name += 0; *nameEWfromTherm += 0; *nameEW += 0;}
      name->Append(".root"); nameEWfromTherm->Append(".root"); nameEW->Append(".root");
      files_2_noG[kt][cb] = new TFile(name->Data(),"READ");
      files_2_noGEWfromTherm[kt][cb] = new TFile(nameEWfromTherm->Data(),"READ");
      files_2_noGEW[kt][cb] = new TFile(nameEW->Data(),"READ");
      //
      
      //
      TString *nameYesG = new TString("OutFiles/OutFileSCPosYesG"); 
      TString *nameYesGEWfromTherm = new TString(nameYesG->Data()); nameYesGEWfromTherm->Append("EWfromThermGRS");
      nameYesG->Append("GRS");
      nameYesG->Append("Kt3_1"); nameYesGEWfromTherm->Append("Kt3_1");
      nameYesG->Append("_Kt"); nameYesGEWfromTherm->Append("_Kt");
      *nameYesG += kt+1; *nameYesGEWfromTherm += kt+1;
      nameYesG->Append("_M"); nameYesGEWfromTherm->Append("_M");
      if(cb<10) {*nameYesG += cb; *nameYesGEWfromTherm += cb;}
      else {*nameYesG += 0; *nameYesGEWfromTherm += 0;}
      nameYesG->Append(".root"); nameYesGEWfromTherm->Append(".root");
      files_2_yesG[kt][cb] = new TFile(nameYesG->Data(),"READ");
      files_2_yesGEWfromTherm[kt][cb] = new TFile(nameYesGEWfromTherm->Data(),"READ");
      //////////////////////////////////////////

      TString *Centname = new TString("Centrality ");
      *Centname += (cb)*5;
      Centname->Append("-");
      *Centname += (cb+1)*5;
      Centname->Append("%");
      
      
      C2_ss[kt][cb] = (TH1D*)files_2_noG[kt][cb]->Get("C2_ss"); C2_ss[kt][cb]->SetDirectory(0);
      C2_os[kt][cb] = (TH1D*)files_2_noG[kt][cb]->Get("C2_os"); C2_os[kt][cb]->SetDirectory(0);
      C2_ss_MomSys[kt][cb] = (TH1D*)files_2_noG[kt][cb]->Get("C2_ss_Momsys"); C2_ss_MomSys[kt][cb]->SetDirectory(0);
      C2_os_MomSys[kt][cb] = (TH1D*)files_2_noG[kt][cb]->Get("C2_os_Momsys"); C2_os_MomSys[kt][cb]->SetDirectory(0);
      C2_ss_KSys[kt][cb] = (TH1D*)files_2_noG[kt][cb]->Get("C2_ss_Ksys"); C2_ss_KSys[kt][cb]->SetDirectory(0);
      C2_os_KSys[kt][cb] = (TH1D*)files_2_noG[kt][cb]->Get("C2_os_Ksys"); C2_os_KSys[kt][cb]->SetDirectory(0);
      fitC2ss_noG[kt][cb] = (TF1*)files_2_noG[kt][cb]->Get("fitC2ss");
      fitC2os_noG[kt][cb] = (TF1*)files_2_noG[kt][cb]->Get("fitC2os");
      hfitC2ss_noG[kt][cb] = (TH1D*)files_2_noG[kt][cb]->Get("fitC2ss_h"); hfitC2ss_noG[kt][cb]->SetDirectory(0);
      hfitC2os_noG[kt][cb] = (TH1D*)files_2_noG[kt][cb]->Get("fitC2os_h"); hfitC2os_noG[kt][cb]->SetDirectory(0);
      fitC2ss_yesG[kt][cb] = (TF1*)files_2_yesG[kt][cb]->Get("fitC2ss");
      fitC2os_yesG[kt][cb] = (TF1*)files_2_yesG[kt][cb]->Get("fitC2os");
      fitC2ss_noGEW[kt][cb] = (TF1*)files_2_noGEW[kt][cb]->Get("fitC2ss");
      fitC2os_noGEW[kt][cb] = (TF1*)files_2_noGEW[kt][cb]->Get("fitC2os");
      hfitC2ss_noGEW[kt][cb] = (TH1D*)files_2_noGEW[kt][cb]->Get("fitC2ss_h"); hfitC2ss_noGEW[kt][cb]->SetDirectory(0);
      hfitC2os_noGEW[kt][cb] = (TH1D*)files_2_noGEW[kt][cb]->Get("fitC2os_h"); hfitC2os_noGEW[kt][cb]->SetDirectory(0);
      fitC2ss_noGEWfromTherm[kt][cb] = (TF1*)files_2_noGEWfromTherm[kt][cb]->Get("fitC2ss");
      fitC2os_noGEWfromTherm[kt][cb] = (TF1*)files_2_noGEWfromTherm[kt][cb]->Get("fitC2os");
      fitC2ss_yesGEWfromTherm[kt][cb] = (TF1*)files_2_yesGEWfromTherm[kt][cb]->Get("fitC2ss");
      fitC2os_yesGEWfromTherm[kt][cb] = (TF1*)files_2_yesGEWfromTherm[kt][cb]->Get("fitC2os");
      fitC2ss_Therm[kt][cb] = (TF1*)files_2_noGEWfromTherm[kt][cb]->Get("fitC2ssTherm");
      fitC2os_Therm[kt][cb] = (TF1*)files_2_noGEWfromTherm[kt][cb]->Get("fitC2osTherm");
      fitC2ss_ThermGaus[kt][cb] = (TF1*)files_2_noGEWfromTherm[kt][cb]->Get("fitC2ssThermGaus");
      fitC2os_ThermGaus[kt][cb] = (TF1*)files_2_noGEWfromTherm[kt][cb]->Get("fitC2osThermGaus");
      K2_ss[kt][cb] = (TH1D*)files_2_noG[kt][cb]->Get("K2_ss"); K2_ss[kt][cb]->SetDirectory(0);
      K2_os[kt][cb] = (TH1D*)files_2_noG[kt][cb]->Get("K2_os"); K2_os[kt][cb]->SetDirectory(0);
      C2Therm_ss[kt][cb] = (TH1D*)files_2_noG[kt][cb]->Get("C2Therm_ss"); C2Therm_ss[kt][cb]->SetDirectory(0);
      C2Therm_os[kt][cb] = (TH1D*)files_2_noG[kt][cb]->Get("C2Therm_os"); C2Therm_os[kt][cb]->SetDirectory(0);
      Chi2NDF_noG[kt][cb] = (TH1D*)files_2_noG[kt][cb]->Get("ChiSquaredNDF"); Chi2NDF_noG[kt][cb]->SetDirectory(0);
      Chi2NDF_noGEWfromTherm[kt][cb] = (TH1D*)files_2_noGEWfromTherm[kt][cb]->Get("ChiSquaredNDF"); Chi2NDF_noGEWfromTherm[kt][cb]->SetDirectory(0);
      Chi2NDF_yesG[kt][cb] = (TH1D*)files_2_yesG[kt][cb]->Get("ChiSquaredNDF"); Chi2NDF_yesG[kt][cb]->SetDirectory(0);
      Chi2NDF_yesGEWfromTherm[kt][cb] = (TH1D*)files_2_yesGEWfromTherm[kt][cb]->Get("ChiSquaredNDF"); Chi2NDF_yesGEWfromTherm[kt][cb]->SetDirectory(0);
      
      //
      for(int par=1; par<=4; par++){
	//if(par!=4 && par!=2) 
	ParHisto_ch[par-1][cb]->SetBinContent(kt+3, fitC2ss_noG[kt][cb]->GetParameter(par));
	ParHisto_ch[par-1][cb]->SetBinError(kt+3, fitC2ss_noG[kt][cb]->GetParError(par));
	ParHisto_coh[par-1][cb]->SetBinContent(kt+3, fitC2ss_yesG[kt][cb]->GetParameter(par));
	ParHisto_coh[par-1][cb]->SetBinError(kt+3, fitC2ss_yesG[kt][cb]->GetParError(par));
	//if(par!=4 && par!=2) 
	ParHisto_chEWfromTherm[par-1][cb]->SetBinContent(kt+3, fitC2ss_noGEWfromTherm[kt][cb]->GetParameter(par));
	ParHisto_chEWfromTherm[par-1][cb]->SetBinError(kt+3, fitC2ss_noGEWfromTherm[kt][cb]->GetParError(par));
	ParHisto_cohEWfromTherm[par-1][cb]->SetBinContent(kt+3, fitC2ss_yesGEWfromTherm[kt][cb]->GetParameter(par));
	ParHisto_cohEWfromTherm[par-1][cb]->SetBinError(kt+3, fitC2ss_yesGEWfromTherm[kt][cb]->GetParError(par));
	ParHisto_chEW[par-1][cb]->SetBinContent(kt+3, fitC2ss_noGEW[kt][cb]->GetParameter(par));
	ParHisto_chEW[par-1][cb]->SetBinError(kt+3, fitC2ss_noGEW[kt][cb]->GetParError(par));
      }
      //////////////////////////////////////
      files_2_noG[kt][cb]->Close();
      files_2_noGEWfromTherm[kt][cb]->Close();
      files_2_noGEW[kt][cb]->Close();
      files_2_yesG[kt][cb]->Close();
      files_2_yesGEWfromTherm[kt][cb]->Close();
    }// kt loop
  }
  
 
  TF1 *Unity = new TF1("Unity","1",0,100);
  Unity->SetLineStyle(2);
  Unity->SetLineColor(1);

 
 

  // merge C3 histogram centralities
  TH1D *C3merged[2][2][2][6];// SC/MC, GRS/Omega0, Therm/Gauss, cb_merged
  TH1D *c3merged[2][2][2][6];// SC/MC, GRS/Omega0, Therm/Gauss, cb_merged
  TH1D *r3merged[2][2][6];// GRS/Omega0, Therm/Gauss, cb_merged
  for(int Coul=0; Coul<2; Coul++){// GRS or Omega0
    for(int ChComb=0; ChComb<2; ChComb++) {
      for(int KT3=0; KT3<2; KT3++) {// Therminator or Gaussian FSI source
	for(int cb=0; cb<6; cb++){
	  C3merged[ChComb][Coul][KT3][cb]=(TH1D*)C3[ChComb][Coul][KT3][cb]->Clone();
	  c3merged[ChComb][Coul][KT3][cb]=(TH1D*)c3[ChComb][Coul][KT3][cb]->Clone();
	  if(ChComb==0) r3merged[Coul][KT3][cb]=(TH1D*)r3_Q3[Coul][KT3][cb]->Clone();
	}      
      }
    }
  }

  for(int MergedBin=0; MergedBin<4; MergedBin++){
    for(int Coul=0; Coul<2; Coul++){// GRS or Omega0
      for(int ChComb=0; ChComb<2; ChComb++) {
	for(int KT3=0; KT3<2; KT3++) {// Therminator or Gaussian FSI source
	  for(int bin=1; bin<=C3[0][Coul][KT3][0]->GetNbinsX(); bin++){
	    int startbin = MergedBin*2 + 2;
	 
	    double valueC3 = C3[ChComb][Coul][KT3][startbin]->GetBinContent(bin) + C3[ChComb][Coul][KT3][startbin+1]->GetBinContent(bin);
	    valueC3 /= 2.;
	    double valueC3_e = pow(C3[ChComb][Coul][KT3][startbin]->GetBinError(bin),2) + pow(C3[ChComb][Coul][KT3][startbin+1]->GetBinError(bin),2);
	    valueC3_e = sqrt(valueC3_e);
	    valueC3_e /= 2.;
	    double valuec3 = c3[ChComb][Coul][KT3][startbin]->GetBinContent(bin) + c3[ChComb][Coul][KT3][startbin+1]->GetBinContent(bin);
	    valuec3 /= 2.;
	    double valuec3_e = pow(c3[ChComb][Coul][KT3][startbin]->GetBinError(bin),2) + pow(c3[ChComb][Coul][KT3][startbin+1]->GetBinError(bin),2);
	    valuec3_e = sqrt(valuec3_e);
	    valuec3_e /= 2.;
	    //
	    C3merged[ChComb][Coul][KT3][MergedBin+2]->SetBinContent(bin, valueC3);
	    C3merged[ChComb][Coul][KT3][MergedBin+2]->SetBinError(bin, valueC3_e);
	    c3merged[ChComb][Coul][KT3][MergedBin+2]->SetBinContent(bin, valuec3);
	    c3merged[ChComb][Coul][KT3][MergedBin+2]->SetBinError(bin, valuec3_e);
	    if(ChComb==0){
	      double valueC3EW = C3EW[Coul][KT3][startbin]->GetBinContent(bin) + C3EW[Coul][KT3][startbin+1]->GetBinContent(bin);
	      valueC3EW /= 2.;
	      double valueC3EW_e = pow(C3EW[Coul][KT3][startbin]->GetBinError(bin),2) + pow(C3EW[Coul][KT3][startbin+1]->GetBinError(bin),2);
	      valueC3EW_e = sqrt(valueC3EW_e);
	      valueC3EW_e /= 2.;
	      double valuer3_Q3 = r3_Q3[Coul][KT3][startbin]->GetBinContent(bin) + r3_Q3[Coul][KT3][startbin+1]->GetBinContent(bin);
	      valuer3_Q3 /= 2.;
	      double valuer3_Q3_e = pow(r3_Q3[Coul][KT3][startbin]->GetBinError(bin),2) + pow(r3_Q3[Coul][KT3][startbin+1]->GetBinError(bin),2);
	      valuer3_Q3_e = sqrt(valuer3_Q3_e);
	      valuer3_Q3_e /= 2.;
	      r3merged[Coul][KT3][MergedBin+2]->SetBinContent(bin, valuer3_Q3);
	      r3merged[Coul][KT3][MergedBin+2]->SetBinError(bin, valuer3_Q3_e);
	      C3EW[Coul][KT3][MergedBin+2]->SetBinContent(bin, valueC3EW);
	      C3EW[Coul][KT3][MergedBin+2]->SetBinError(bin, valueC3EW_e);
	    }
	    //
	  }
	}// KT3
      }
    }
  }// Merged Bin



  //////////////////////////////////////////////////////////////////////////
  // Draw K3 factors
  TCanvas *canK3 = new TCanvas("canK3", "canK3",10,0,600,900);// 11,53,700,500
  canK3->SetHighLightColor(2);
  //canK3->Range(-0.7838115,-1.033258,0.7838115,1.033258);
  gStyle->SetOptFit(0111);
  canK3->SetFillColor(0);//10
  canK3->SetBorderMode(0);
  canK3->SetBorderSize(2);
  canK3->SetFrameFillColor(0);
  canK3->SetFrameBorderMode(0);
  canK3->SetFrameBorderMode(0);
  
  TPad *padK3 = new TPad("padK3","padK3",0.0,0.0,1.,1.);
  gPad->SetGridx(0);
  gPad->SetGridy(0);
  gPad->SetTickx();
  gPad->SetTicky();
  padK3->SetTopMargin(0.02);//0.05
  padK3->SetRightMargin(0.01);//1e-2
  padK3->SetBottomMargin(0.07);//0.12
  padK3->Divide(1,2,0,0);
  padK3->Draw();
  padK3->cd(1);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.03);
  //double Dim1=gPad->GetAbsHNDC();
  //cout<<gPad->GetAbsHNDC()<<endl;
  //cout<<gPad->GetAspectRatio()<<endl;
  TLegend *legendK3 = new TLegend(.45,.05,.95,.3,NULL,"brNDC");
  legendK3->SetBorderSize(0);
  legendK3->SetTextSize(SizeLabel*SF1);// *2/3.*0.96
  legendK3->SetFillColor(0);

  K3[0][1][KT3Bin][BOI_1]->SetMinimum(0.48);
  K3[0][1][KT3Bin][BOI_1]->SetMaximum(1.33);

  K3[0][1][KT3Bin][BOI_1]->SetMarkerStyle(20);
  K3[0][1][KT3Bin][BOI_1]->SetMarkerColor(2);
  K3[0][0][KT3Bin][BOI_1]->SetMarkerStyle(24);
  K3[0][0][KT3Bin][BOI_1]->SetMarkerColor(2);
  K3[1][1][KT3Bin][BOI_1]->SetMarkerStyle(22);
  K3[1][1][KT3Bin][BOI_1]->SetMarkerColor(4);
  K3[1][0][KT3Bin][BOI_1]->SetMarkerStyle(26);
  K3[1][0][KT3Bin][BOI_1]->SetMarkerColor(4);
  //
  K3[0][1][KT3Bin][BOI_1]->GetXaxis()->SetLabelFont(TextFont); K3[0][1][KT3Bin][BOI_1]->GetYaxis()->SetLabelFont(TextFont);
  K3[0][1][KT3Bin][BOI_1]->GetXaxis()->SetLabelSize(SizeLabel*SF1); K3[0][1][KT3Bin][BOI_1]->GetYaxis()->SetLabelSize(SizeLabel*SF1);
  K3[0][1][KT3Bin][BOI_1]->GetXaxis()->SetNdivisions(808); K3[0][1][KT3Bin][BOI_1]->GetYaxis()->SetNdivisions(909);
  TGaxis::SetMaxDigits(3); 
  K3[0][1][KT3Bin][BOI_1]->GetYaxis()->SetTitleFont(TextFont);
  K3[0][1][KT3Bin][BOI_1]->GetYaxis()->SetTitleSize(SizeTitle*SF1);
  K3[0][1][KT3Bin][BOI_1]->GetYaxis()->SetTitle("#font[12]{K_{3}}");
  K3[0][1][KT3Bin][BOI_1]->GetYaxis()->SetTitleOffset(0.8);

  //
  K3[0][1][KT3Bin][BOI_1]->Draw();
  K3[0][0][KT3Bin][BOI_1]->Draw("same");
  K3[1][1][KT3Bin][BOI_1]->Draw("same");
  K3[1][0][KT3Bin][BOI_1]->Draw("same");
  Unity->Draw("same");
  legendK3->AddEntry(K3[0][1][KT3Bin][BOI_1],"same-charge, #Omega_{0}","p");
  legendK3->AddEntry(K3[0][0][KT3Bin][BOI_1],"same-charge, GRS","p");
  legendK3->AddEntry(K3[1][1][KT3Bin][BOI_1],"mixed-charge, #Omega_{0}","p");
  legendK3->AddEntry(K3[1][0][KT3Bin][BOI_1],"mixed-charge, GRS","p");
  legendK3->Draw("same");
  //TLatex *K3Label = new TLatex(0.02,1.25,"K_{3}");// -0.011,0.92 (ss), 1.26 (os)
  //K3Label->SetTextFont(TextFont);
  //K3Label->SetTextSize(SizeTitle*SF1);
  //K3Label->SetTextAngle(90);
  //K3Label->Draw();
  //
  padK3->cd(2);
  double SF_correction=0.97;
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.03);
  gPad->SetBottomMargin(0.16);
  //double Dim2=gPad->GetAbsHNDC();
  TLegend *legendK3comp = new TLegend(.45,.8,.95,.95,NULL,"brNDC");
  legendK3comp->SetBorderSize(0);
  legendK3comp->SetTextSize(SizeLabel*SF1*SF_correction);// small .03; large .06
  legendK3comp->SetFillColor(0);

  TH1D *K3sc_compOmega0 = (TH1D*)K3[0][1][KT3Bin][BOI_1]->Clone();
  TH1D *K3sc_compGRS = (TH1D*)K3[0][0][KT3Bin][BOI_1]->Clone();
  TH1D *K3mc_compOmega0 = (TH1D*)K3[1][1][KT3Bin][BOI_1]->Clone();
  TH1D *K3mc_compGRS = (TH1D*)K3[1][0][KT3Bin][BOI_1]->Clone();
  for(int ii=0; ii<K3sc_compOmega0->GetNbinsX(); ii++){ 
    K3sc_compOmega0->SetBinContent(ii+1, K3sc_compOmega0->GetBinContent(ii+1)-1.0);
    K3sc_compGRS->SetBinContent(ii+1, K3sc_compGRS->GetBinContent(ii+1)-1.0);
    K3mc_compOmega0->SetBinContent(ii+1, K3mc_compOmega0->GetBinContent(ii+1)-1.0);
    K3mc_compGRS->SetBinContent(ii+1, K3mc_compGRS->GetBinContent(ii+1)-1.0);
  }
  K3sc_compOmega0->Add(K3sc_compGRS,-1);
  K3sc_compOmega0->Divide(K3sc_compGRS);
  K3mc_compOmega0->Add(K3mc_compGRS,-1);
  K3mc_compOmega0->Divide(K3mc_compGRS);
  K3sc_compOmega0->SetMarkerStyle(20);
  K3sc_compOmega0->SetMarkerColor(2);
  K3sc_compOmega0->SetLineColor(2);
  K3mc_compOmega0->SetMarkerStyle(24);
  K3mc_compOmega0->SetMarkerColor(4);
  K3mc_compOmega0->SetLineColor(4);
  K3sc_compOmega0->SetMinimum(-0.12);// -0.021
  K3sc_compOmega0->SetMaximum(0.12);// 0.021
  K3sc_compOmega0->SetBinContent(1,-100); K3mc_compOmega0->SetBinContent(1,-100);
  K3sc_compOmega0->GetYaxis()->SetNdivisions(404);
  K3sc_compOmega0->GetYaxis()->SetTitleFont(TextFont); K3sc_compOmega0->GetYaxis()->SetTitleSize(SizeTitle*SF1*SF_correction);
  //K3sc_compOmega0->GetXaxis()->SetTitleFont(TextFont); K3sc_compOmega0->GetXaxis()->SetTitleSize(SizeTitle*2/3.*.96);
  //K3sc_compOmega0->GetXaxis()->SetTitle("Q_{3} (GeV/c)");
  K3sc_compOmega0->GetYaxis()->SetTitle("#Delta#font[12]{K_{3} / (K_{3}}(#Omega_{0})-1)");
  K3sc_compOmega0->GetXaxis()->SetTitleOffset(1.2); K3sc_compOmega0->GetYaxis()->SetTitleOffset(0.8);
  K3sc_compOmega0->Draw();
  K3mc_compOmega0->Draw("same");
  TF1 *Zero=new TF1("Zero","0",0,1);
  Zero->SetLineStyle(2); Zero->SetLineColor(1);
  Zero->Draw("same");
  legendK3comp->AddEntry(K3sc_compOmega0,"same-charge","p");
  legendK3comp->AddEntry(K3mc_compOmega0,"mixed-charge","p");
  legendK3comp->Draw("same");
  //TLatex *RatioLabel = new TLatex(-.011,0.02,"#Delta K_{3}/(K_{3}(#Omega_{0})-1)");// -0.011,0.04
  //RatioLabel->SetTextFont(TextFont);
  //RatioLabel->SetTextSize(SizeTitle*SF1);
  //RatioLabel->SetTextAngle(90);
  //RatioLabel->Draw();

  TLatex *Q3Label = new TLatex(.071,-0.155,"#font[12]{Q}_{3} (GeV/#font[12]{c})");//.065,-0.147
  Q3Label->SetTextFont(TextFont);
  Q3Label->SetTextSize(SizeTitle*SF1);// 0.08
  Q3Label->Draw();



 
  
  //////////////////////////////////////////
  // print global fit values

  // kappa3 and kappa4
  cout<<"ALICE EW"<<endl;
  cout<<"kappa3"<<endl;
  for(int cb=0; cb<10; cb++){
    if(cb!=0 && cb!=9) continue;
    TString *Cname = new TString("$");
    *Cname += cb*5; Cname->Append("-");
    *Cname += (cb+1)*5; Cname->Append("\\%$");
    cout<<Cname->Data()<<" & ";
    cout.precision(2);
    for(int kt=0; kt<KTBINS; kt++) {
      cout<<fitC2ss_noGEW[kt][cb]->GetParameter(5)<<" ";
      if((kt+1)!=KTBINS) cout<<"& ";
      else cout<<" \\\\\ \\hline";
    }
    cout<<endl;
  }
  //
  cout<<"kappa4"<<endl;
  for(int cb=0; cb<10; cb++){
    if(cb!=0 && cb!=9) continue;
    TString *Cname = new TString("$");
    *Cname += cb*5; Cname->Append("-");
    *Cname += (cb+1)*5; Cname->Append("\\%$");
    cout<<Cname->Data()<<" & ";
    cout.precision(2);
    for(int kt=0; kt<KTBINS; kt++) {
      cout<<fitC2ss_noGEW[kt][cb]->GetParameter(6)<<" ";
      if((kt+1)!=KTBINS) cout<<"& ";
      else cout<<" \\\\\ \\hline";
    }
    cout<<endl;
  }
  // Therminator
  cout<<"Therminator EW"<<endl;
  cout<<"kappa3"<<endl;
  for(int cb=0; cb<10; cb++){
    if(cb!=0 && cb!=9) continue;
    TString *Cname = new TString("$");
    *Cname += cb*5; Cname->Append("-");
    *Cname += (cb+1)*5; Cname->Append("\\%$");
    cout<<Cname->Data()<<" & ";
    cout.precision(2);
    for(int kt=0; kt<KTBINS; kt++) {
      cout<<fitC2ss_noGEWfromTherm[kt][cb]->GetParameter(5)<<" ";
      if((kt+1)!=KTBINS) cout<<"& ";
      else cout<<" \\\\\ \\hline";
    }
    cout<<endl;
  }
  //
  cout<<"kappa4"<<endl;
  for(int cb=0; cb<10; cb++){
    if(cb!=0 && cb!=9) continue;
    TString *Cname = new TString("$");
    *Cname += cb*5; Cname->Append("-");
    *Cname += (cb+1)*5; Cname->Append("\\%$");
    cout<<Cname->Data()<<" & ";
    cout.precision(2);
    for(int kt=0; kt<KTBINS; kt++) {
      cout<<fitC2ss_noGEWfromTherm[kt][cb]->GetParameter(6)<<" ";
      if((kt+1)!=KTBINS) cout<<"& ";
      else cout<<" \\\\\ \\hline";
    }
    cout<<endl;
  }




  // Chi2/NDF
  /*cout<<"Chi2/NDF"<<endl;
  cout<<"fits with Coherence"<<endl;
  for(int cb=0; cb<10; cb++){
    if(cb!=0 && cb!=9) continue;
    TString *Cname = new TString("$");
    *Cname += cb*5; Cname->Append("-");
    *Cname += (cb+1)*5; Cname->Append("\\%$");
    cout<<Cname->Data()<<" & ";
    cout.precision(2);
    for(int kt=0; kt<KTBINS; kt++) {
      //cout<<fitC2ss_yesGEWfromTherm[kt][cb]->GetParameter(4)<<" ";
      //cout<<fitC2ss_yesG[kt][cb]->GetParameter(4)<<" ";
      //cout<<fitC2ss_noGEWfromTherm[kt][cb]->GetParameter(3)<<" ";
      cout<<Chi2NDF_yesG[kt][cb]->GetBinContent(1)<<" ";// Chi^2/NDF
      if((kt+1)!=KTBINS) cout<<"& ";
      else cout<<" \\\\\ \\hline";
    }
    cout<<endl;
  }
  //
  cout<<"fits with Therminator EW"<<endl;
  for(int cb=0; cb<10; cb++){
    if(cb!=0 && cb!=9) continue;
    TString *Cname = new TString("$");
    *Cname += cb*5; Cname->Append("-");
    *Cname += (cb+1)*5; Cname->Append("\\%$");
    cout<<Cname->Data()<<" & ";
    cout.precision(2);
    for(int kt=0; kt<KTBINS; kt++) {
      cout<<Chi2NDF_noGEWfromTherm[kt][cb]->GetBinContent(1)<<" ";// Chi^2/NDF
      if((kt+1)!=KTBINS) cout<<"& ";
      else cout<<" \\\\\ \\hline";
    }
    cout<<endl;
  }
  */


  
  /////////////////////////////////////////////
  // Draw C2 in 6 kt bins
  TCanvas *canC2 = new TCanvas("canC2", "canC2",10,0,600,900);// 11,53,700,500
  canC2->SetHighLightColor(2);
  //canC2->Range(-0.7838115,-1.033258,0.7838115,1.033258);
  gStyle->SetOptFit(0111);
  canC2->SetFillColor(0);//10
  canC2->SetBorderMode(0);
  canC2->SetBorderSize(2);
  canC2->SetFrameFillColor(0);
  canC2->SetFrameBorderMode(0);
  canC2->SetFrameBorderMode(0);
  
  
  TPad *pad1 = new TPad("pad1","pad1",0.06,0.06,1.,1.);//0.05,0.05,1.,1.
  gPad->SetGridx(0);
  gPad->SetGridy(0);
  gPad->SetTickx();
  gPad->SetTicky();
  pad1->SetTopMargin(0.02);//0.05
  pad1->SetLeftMargin(0.12);//
  pad1->SetRightMargin(0.01);//3e-2
  pad1->SetBottomMargin(0.08);//0.12
  pad1->Divide(2,3,0,0);
  pad1->Draw();
 

  TLegend *legend1 = new TLegend(.2,.75,.99,.95,NULL,"brNDC");//.45 or .4 for x1
  legend1->SetBorderSize(0);
  //legend1->SetTextSize(.08);// small .03; large .036 
  legend1->SetTextFont(TextFont);
  legend1->SetTextSize(SizeLegend);
  legend1->SetFillColor(0);
  //
  TLegend *legend2 = new TLegend(.35,.74, .99,.98,NULL,"brNDC");//.45 or .4 for x1
  legend2->SetBorderSize(0);
  //legend2->SetTextSize(.08);// small .03; large .036 
  legend2->SetTextFont(TextFont);
  legend2->SetTextSize(SizeLegend);
  legend2->SetFillColor(0);


  TGaxis::SetMaxDigits(4);
  /////////////////////
  for(int kt=0; kt<KTBINS; kt++) {
    pad1->cd(kt+1);
    gPad->SetRightMargin(RightMargin);
    if(kt%2==0) gPad->SetRightMargin(0.001);
    if(kt==4) SF_correction=0.92;
    else if(kt==5) SF_correction=1.015;
    else SF_correction=1.0;
    //if(kt==0) Dim1=gPad->GetAbsWNDC()*gPad->GetAbsHNDC();
    //Dim2=gPad->GetAbsWNDC()*gPad->GetAbsHNDC();
    C2_ss[kt][BOI_1]->SetMinimum(MIN[BOI_1]);
    C2_ss[kt][BOI_1]->SetMaximum(MAX[BOI_1]);
    C2_ss[kt][BOI_1]->GetXaxis()->SetRangeUser(0,0.061);
    C2_ss[kt][BOI_1]->GetXaxis()->SetTitleOffset(0.85);
    C2_ss[kt][BOI_1]->GetYaxis()->SetTitleOffset(0.85);// 0.85
    C2_ss[kt][BOI_1]->GetXaxis()->SetTitle(""); C2_ss[kt][BOI_1]->GetYaxis()->SetTitle("");
    C2_ss[kt][BOI_1]->GetXaxis()->SetTitleFont(TextFont); C2_ss[kt][BOI_1]->GetYaxis()->SetTitleFont(TextFont);
    C2_ss[kt][BOI_1]->GetXaxis()->SetLabelFont(TextFont); C2_ss[kt][BOI_1]->GetYaxis()->SetLabelFont(TextFont);
    C2_ss[kt][BOI_1]->GetXaxis()->SetLabelSize(SizeLabel*SF_correction); C2_ss[kt][BOI_1]->GetYaxis()->SetLabelSize(SizeLabel*SF_correction);
    if(kt%2==0) C2_ss[kt][BOI_1]->GetYaxis()->SetLabelSize(SizeLabel*SF_correction);
    else {C2_ss[kt][BOI_1]->GetYaxis()->SetLabelSize(.0);}
    if(kt<KTBINS-2) C2_ss[kt][BOI_1]->GetXaxis()->SetLabelSize(.0);
    C2_ss[kt][BOI_1]->GetXaxis()->SetNdivisions(603); C2_ss[kt][BOI_1]->GetYaxis()->SetNdivisions(606);
    
    //
    C2Therm_ss[kt][BOI_1]->SetMinimum(MIN[BOI_1]);
    C2Therm_ss[kt][BOI_1]->SetMaximum(MAX[BOI_1]);
    C2Therm_ss[kt][BOI_1]->GetXaxis()->SetRangeUser(0,0.061);
    C2Therm_ss[kt][BOI_1]->GetXaxis()->SetTitleOffset(0.85);
    C2Therm_ss[kt][BOI_1]->GetYaxis()->SetTitleOffset(0.85);
    C2Therm_ss[kt][BOI_1]->GetXaxis()->SetTitle(""); C2Therm_ss[kt][BOI_1]->GetYaxis()->SetTitle("");
    C2Therm_ss[kt][BOI_1]->GetXaxis()->SetLabelFont(TextFont); C2Therm_ss[kt][BOI_1]->GetYaxis()->SetLabelFont(TextFont); 
    if(kt%2==0) C2Therm_ss[kt][BOI_1]->GetYaxis()->SetLabelSize(SizeLabel*SF_correction);
    else {C2Therm_ss[kt][BOI_1]->GetYaxis()->SetLabelSize(.0);}
    if(kt<KTBINS-2) C2Therm_ss[kt][BOI_1]->GetXaxis()->SetLabelSize(.0);
    else {C2Therm_ss[kt][BOI_1]->GetXaxis()->SetLabelSize(SizeLabel*SF_correction);}
    C2Therm_ss[kt][BOI_1]->GetXaxis()->SetNdivisions(606); C2Therm_ss[kt][BOI_1]->GetYaxis()->SetNdivisions(606);

    if(BOI_1==9 && kt==5) {// Point with large statistical error
      C2_ss[kt][BOI_1]->SetBinContent(3,-100);
      C2_ss_MomSys[kt][BOI_1]->SetBinContent(3,-100);
    }

    if(TherminatorC2==kFALSE){
      C2_ss[kt][BOI_1]->DrawCopy();
      C2_ss_MomSys[kt][BOI_1]->SetFillColor(kRed-9); C2_os_MomSys[kt][BOI_1]->SetFillColor(kBlue-9);
      C2_ss_MomSys[kt][BOI_1]->DrawCopy("E2 same");
      C2_os_MomSys[kt][BOI_1]->DrawCopy("E2 same");
      C2_ss[kt][BOI_1]->DrawCopy("same");
      C2_os[kt][BOI_1]->DrawCopy("same");
      hfitC2ss_noG[kt][BOI_1]->SetLineStyle(3); fitC2ss_noG[kt][BOI_1]->SetLineStyle(3); 
      hfitC2os_noG[kt][BOI_1]->SetLineStyle(3); fitC2os_noG[kt][BOI_1]->SetLineStyle(3);
      hfitC2ss_noG[kt][BOI_1]->DrawCopy("C same");
      hfitC2os_noG[kt][BOI_1]->DrawCopy("C same");
      hfitC2ss_noGEW[kt][BOI_1]->DrawCopy("C same");
      hfitC2os_noGEW[kt][BOI_1]->DrawCopy("C same");
      //fitC2ss_yesG[kt][BOI_1]->DrawCopy("same");
      //fitC2os_yesG[kt][BOI_1]->DrawCopy("same");

    }else{// Therminator
      C2Therm_ss[kt][BOI_1]->DrawCopy();
      C2Therm_os[kt][BOI_1]->DrawCopy("same");
      fitC2ss_Therm[kt][BOI_1]->DrawCopy("same");
      fitC2os_Therm[kt][BOI_1]->DrawCopy("same");
      fitC2ss_ThermGaus[kt][BOI_1]->SetLineStyle(3);
      fitC2os_ThermGaus[kt][BOI_1]->SetLineStyle(3);
      fitC2ss_ThermGaus[kt][BOI_1]->DrawCopy("same");
      fitC2os_ThermGaus[kt][BOI_1]->DrawCopy("same");
    }
    //
    if(kt==0){
      TLatex *Specif;
      if(!TherminatorC2) Specif = new TLatex(0.007,MIN[BOI_1]+(MAX[BOI_1]-MIN[BOI_1])*0.91,"ALICE Pb-Pb #\sqrt{#font[12]{s}_{NN}}=2.76 TeV");// ALICE specifications
      else Specif = new TLatex(0.0015,MIN[BOI_1]+(MAX[BOI_1]-MIN[BOI_1])*0.91,"Therminator");// Therminator
      Specif->SetTextFont(TextFont);
      Specif->SetTextSize(SizeSpecif*.9);
      TString *CentName;
      if(!TherminatorC2) {CentName = new TString(); *CentName += BOI_1*5; CentName->Append("-"); *CentName += (BOI_1+1)*5; CentName->Append("%");}
      else {
	if(BOI_1==0) CentName = new TString("0-5%");
	if(BOI_1==9) CentName = new TString("40-50%");
      }
      TLatex *Specif2 = new TLatex(0.027,MIN[BOI_1]+(MAX[BOI_1]-MIN[BOI_1])*0.8,CentName->Data());
      Specif2->SetTextFont(TextFont);
      Specif2->SetTextSize(SizeSpecif);
      Specif->Draw("same");
      Specif2->Draw("same");
    }
    if(kt==1){
      legend1->AddEntry(C2_ss[kt][BOI_1],"same-charge","p");
      legend1->AddEntry(C2_os[kt][BOI_1],"mixed-charge","p");
      legend1->Draw("same");
    }
    if(kt==2){
      fitC2ss_noG[kt][BOI_1]->SetLineColor(1);
      fitC2ss_yesG[kt][BOI_1]->SetLineColor(1);
      fitC2ss_noGEW[kt][BOI_1]->SetLineColor(1);
      fitC2ss_Therm[kt][BOI_1]->SetLineColor(1);
      fitC2ss_ThermGaus[kt][BOI_1]->SetLineColor(1);
      if(!TherminatorC2) {
	legend2->AddEntry(fitC2ss_noG[kt][BOI_1],"Gauss","l");
	legend2->AddEntry(fitC2ss_noGEW[kt][BOI_1],"Edgeworth","l");
	//legend2->AddEntry(fitC2ss_yesG[kt][BOI_1],"Gauss, G#neq0","l");
      }else {
	legend2->AddEntry(fitC2ss_ThermGaus[kt][BOI_1],"Gauss","l");// Therminator
	legend2->AddEntry(fitC2ss_Therm[kt][BOI_1],"Edgeworth","l");// Therminator
      }      
      legend2->Draw("same");
    }
    TLatex *KtLabel;
    if(kt==0) KtLabel = new TLatex(0.028,MIN[BOI_1]+(MAX[BOI_1]-MIN[BOI_1])*0.68,"0.2<#font[12]{k}_{T}<0.3 GeV/#font[12]{c}");// 0.003,.988
    else if(kt==1) KtLabel = new TLatex(0.028,MIN[BOI_1]+(MAX[BOI_1]-MIN[BOI_1])*0.68,"0.3<#font[12]{k}_{T}<0.4 GeV/#font[12]{c}");
    else if(kt==2) KtLabel = new TLatex(0.028,MIN[BOI_1]+(MAX[BOI_1]-MIN[BOI_1])*0.68,"0.4<#font[12]{k}_{T}<0.5 GeV/#font[12]{c}");
    else if(kt==3) KtLabel = new TLatex(0.028,MIN[BOI_1]+(MAX[BOI_1]-MIN[BOI_1])*0.68,"0.5<#font[12]{k}_{T}<0.6 GeV/#font[12]{c}");
    else if(kt==4) KtLabel = new TLatex(0.028,MIN[BOI_1]+(MAX[BOI_1]-MIN[BOI_1])*0.68,"0.6<#font[12]{k}_{T}<0.7 GeV/#font[12]{c}");
    else KtLabel = new TLatex(0.028,MIN[BOI_1]+(MAX[BOI_1]-MIN[BOI_1])*0.68,"0.7<#font[12]{k}_{T}<0.8 GeV/#font[12]{c}");
    KtLabel->SetTextFont(TextFont);// 23
    KtLabel->SetTextSize(SizeSpecif*SF_correction);// 24
    KtLabel->Draw("same");
  }
 
  canC2->cd();
    
  TPad *pad1_2 = new TPad("pad1_2","pad1_2",0.0,0.0,1.,1.);
  pad1_2->SetFillStyle(0);
  pad1_2->Draw();
  pad1_2->cd();
  TLatex *C2Label = new TLatex(.04,.94,"#font[12]{C}_{2}");// 0.04,.91
  C2Label->SetTextFont(TextFont);
  C2Label->SetTextSize(SizeTitle*SF2);
  C2Label->SetTextAngle(90);
  C2Label->Draw();
  TLatex *Q2Label = new TLatex(.75,.015,"#font[12]{q} (GeV/#font[12]{c})");// .7,.02
  Q2Label->SetTextFont(TextFont);
  Q2Label->SetTextSize(SizeTitle*SF2);
  Q2Label->Draw();
  TBox *CoverUp = new TBox(0.554,0.05,0.58,.083);
  CoverUp->SetFillColor(0);
  CoverUp->Draw();
  



  ////////////////////////////////////////////////////
  ////////////////////////////////////////////////////
  // C2 fit parameters
  TCanvas *can2 = new TCanvas("can2", "can2",680,0,600,900);// 680,0,600,900 with G
  can2->SetHighLightColor(2);
  gStyle->SetOptFit(0111);
  can2->SetFillColor(0);//10
  can2->SetBorderMode(0);
  can2->SetBorderSize(2);
  can2->SetFrameFillColor(0);
  can2->SetFrameBorderMode(0);
  can2->SetFrameBorderMode(0);
  can2->cd();

  TPad *pad2 = new TPad("pad2","pad2",0.,0.,1.,1.);
  gPad->SetGridx(0);
  gPad->SetGridy(0);
  gPad->SetTickx();
  gPad->SetTicky();
  pad2->SetTopMargin(0.02);
  pad2->SetLeftMargin(0.0);
  pad2->SetRightMargin(0.0);
  pad2->SetBottomMargin(0.07);
  pad2->Divide(1,2,0,0);// 1,3,0,0 with G
  pad2->Draw();
  pad2->cd();
  TLegend *legend3 = new TLegend(.44,.66, .96,.98,NULL,"brNDC");// .56,.6, .96,.98
  legend3->SetBorderSize(0);
  legend3->SetTextFont(TextFont);
  legend3->SetTextSize(SizeLegend*SF1);
  legend3->SetFillColor(0);
    
  double LowerLimits[4]={0.46, 0.32, 4.8, 0.08};// 0.46, 0.32, 4.8, 0.08
  double UpperLimits[4]={0.87, 0.63, 13.2, 1.28};// 0.84, 0.63, 13.2, 1.28
  for(int par=1; par<=3; par++){
    //pad2->cd(par);
    if(par!=1 && par!=3) continue;
    if(par==1) pad2->cd(1);
    else pad2->cd(2);
    gPad->SetRightMargin(0.03); gPad->SetLeftMargin(0.14);
    if(par==3) gPad->SetBottomMargin(0.16);// 0.22 with G
    // settings applied to ParHisto_coh[par-1][BOI_1] to include G.  
    ParHisto_ch[par-1][BOI_1]->SetMinimum(LowerLimits[par-1]);
    ParHisto_ch[par-1][BOI_1]->SetMaximum(UpperLimits[par-1]);
    ParHisto_ch[par-1][BOI_1]->GetXaxis()->SetTitleFont(TextFont); ParHisto_ch[par-1][BOI_1]->GetYaxis()->SetTitleFont(TextFont);
    ParHisto_ch[par-1][BOI_1]->GetXaxis()->SetLabelFont(TextFont); ParHisto_ch[par-1][BOI_1]->GetYaxis()->SetLabelFont(TextFont);
    if(par==3){
      ParHisto_ch[par-1][BOI_1]->GetXaxis()->SetLabelSize(SizeLabel*SF1);
      ParHisto_ch[par-1][BOI_1]->GetXaxis()->SetTitleSize(SizeTitle*SF1);
      ParHisto_ch[par-1][BOI_1]->GetXaxis()->SetLabelSize(SizeLabel*SF1);
    }else{
      ParHisto_ch[par-1][BOI_1]->GetXaxis()->SetLabelSize(0);
      ParHisto_ch[par-1][BOI_1]->GetXaxis()->SetTitleSize(0);
      ParHisto_ch[par-1][BOI_1]->GetXaxis()->SetLabelSize(0);
    }
    if(par==1) SF_correction=1.0;
    else SF_correction=0.97;
    ParHisto_ch[par-1][BOI_1]->GetYaxis()->SetTitleSize(SizeTitle*SF1*SF_correction);
    ParHisto_ch[par-1][BOI_1]->GetYaxis()->SetLabelSize(SizeLabel*SF1*SF_correction);
    ParHisto_ch[par-1][BOI_1]->GetXaxis()->SetTitle("#font[12]{k}_{T} (GeV/#font[12]{c})");
    ParHisto_ch[par-1][BOI_1]->GetXaxis()->SetRangeUser(0.2,0.79);
    if(par==1) ParHisto_ch[par-1][BOI_1]->GetYaxis()->SetTitle("#lambda");
    if(par==2) ParHisto_ch[par-1][BOI_1]->GetYaxis()->SetTitle("#font[12]{G}");
    if(par==3) ParHisto_ch[par-1][BOI_1]->GetYaxis()->SetTitle("#font[12]{R_{ch}} (fm)");
    if(par==4) ParHisto_ch[par-1][BOI_1]->GetYaxis()->SetTitle("#font[12]{R_{coh}} (fm)");
    ParHisto_ch[par-1][BOI_1]->GetXaxis()->SetTitleOffset(0.9);// 2.0
    ParHisto_ch[par-1][BOI_1]->GetYaxis()->SetTitleOffset(0.8);// 1.4
    ParHisto_ch[par-1][BOI_1]->GetYaxis()->SetNdivisions(505);
   
    // ParN=2 (lambda), 3(G), 4(Rch), 5(Rcoh)
    ParHisto_ch[par-1][BOI_1]->Draw();
    
    
    TH1D *ParHisto1_Sys=(TH1D*)ParHisto_ch[par-1][BOI_1]->Clone();
    TH1D *ParHisto2_Sys=(TH1D*)ParHisto_ch[par-1][BOI_2]->Clone();
    TH1D *ParHisto3_Sys=(TH1D*)ParHisto_chEW[par-1][BOI_1]->Clone();
    TH1D *ParHisto4_Sys=(TH1D*)ParHisto_chEW[par-1][BOI_2]->Clone();
    TF1 *RSys1=new TF1("RSys1","0.3",0,1);// M0, R Gauss Fit range syst in fm.
    TF1 *RSys2=new TF1("RSys2","0.2",0,1);// M9, R Gauss Fit range syst in fm.
    TF1 *RSys3=new TF1("RSys3","0.2",0,1);// M0, R EW Fit range syst in fm.
    TF1 *RSys4=new TF1("RSys4","0.2",0,1);// M9, R EW Fit range syst in fm.
    TF1 *LamSys1=new TF1("LamSys1","0.015",0,1);// M0, Lam Gauss Fit range syst in fm.
    TF1 *LamSys2=new TF1("LamSys2","0.015",0,1);// M9, Lam Gauss Fit range syst in fm.
    TF1 *LamSys3=new TF1("LamSys3","0.01",0,1);// M0, Lam EW Fit range syst in fm.
    TF1 *LamSys4=new TF1("LamSys4","0.01",0,1);// M9, Lam EW Fit range syst in fm.

    for(int bin=1; bin<=ParHisto1_Sys->GetNbinsX(); bin++){
      double ktvalue=ParHisto1_Sys->GetXaxis()->GetBinCenter(bin);
      if(par==1){// lambda Sys
	ParHisto1_Sys->SetBinError(bin, LamSys1->Eval(ktvalue));
	ParHisto2_Sys->SetBinError(bin, LamSys2->Eval(ktvalue));
	ParHisto3_Sys->SetBinError(bin, LamSys3->Eval(ktvalue));
	ParHisto4_Sys->SetBinError(bin, LamSys4->Eval(ktvalue));
      }
      if(par==3){// R sys
	ParHisto1_Sys->SetBinError(bin, RSys1->Eval(ktvalue));
	ParHisto2_Sys->SetBinError(bin, RSys2->Eval(ktvalue));
	ParHisto3_Sys->SetBinError(bin, RSys3->Eval(ktvalue));
	ParHisto4_Sys->SetBinError(bin, RSys4->Eval(ktvalue));
      }
    }
    ParHisto1_Sys->SetMarkerSize(0); ParHisto1_Sys->SetFillStyle(1000); ParHisto1_Sys->SetFillColor(kGray);
    ParHisto2_Sys->SetMarkerSize(0); ParHisto2_Sys->SetFillStyle(1000); ParHisto2_Sys->SetFillColor(kGray);
    ParHisto3_Sys->SetMarkerSize(0); ParHisto3_Sys->SetFillStyle(1000); ParHisto3_Sys->SetFillColor(kBlue-10);
    ParHisto4_Sys->SetMarkerSize(0); ParHisto4_Sys->SetFillStyle(1000); ParHisto4_Sys->SetFillColor(kBlue-10);

    ParHisto1_Sys->Draw("E2 same"); ParHisto2_Sys->Draw("E2 same"); ParHisto3_Sys->Draw("E2 same"); ParHisto4_Sys->Draw("E2 same");
    ParHisto_ch[par-1][BOI_1]->Draw("same");
    ParHisto_ch[par-1][BOI_2]->Draw("same");
    ParHisto_chEW[par-1][BOI_1]->Draw("same"); 
    ParHisto_chEW[par-1][BOI_2]->Draw("same");
    
    if(par==1){
      //legend3->AddEntry(ParHisto_chEWfromTherm[par-1][BOI_1],"E_{w}, 0-5%","p"); legend3->AddEntry(ParHisto_chEWfromTherm[par-1][BOI_2],"E_{w}, 45-50%","p");
      //legend3->AddEntry(ParHisto_coh[par-1][BOI_1],"G#neq0, 0-5%","p"); legend3->AddEntry(ParHisto_coh[par-1][BOI_2],"G#neq0, 45-50%","p");
      legend3->AddEntry(ParHisto_ch[par-1][BOI_1],"Gauss, 0-5%","p"); legend3->AddEntry(ParHisto_ch[par-1][BOI_2],"Gauss, 45-50%","p");
      legend3->AddEntry(ParHisto_chEW[par-1][BOI_1],"Edgeworth, 0-5%","p"); legend3->AddEntry(ParHisto_chEW[par-1][BOI_2],"Edgeworth, 45-50%","p");
    }
    
    if(par==3) legend3->Draw("same");// par==2 with G
    if(par==1){
      TLatex *Specif_2 = new TLatex(0.35,LowerLimits[par-1]+(UpperLimits[par-1]-LowerLimits[par-1])*.94,"ALICE Pb-Pb #\sqrt{#font[12]{s}_{NN}}=2.76 TeV");
      Specif_2->SetTextFont(TextFont);
      Specif_2->SetTextSize(SizeSpecif*SF1);
      Specif_2->Draw();
    }
  }

  ////////////////////////////////////////////////////
  ////////////////////////////////////////////////////
  // C2(+-) Therminator comparisons
  
  TCanvas *can3 = new TCanvas("can3", "can3",1800,0,600,600);// 11,53,700,500
  can3->SetHighLightColor(2);
  gStyle->SetOptFit(0111);
  can3->SetFillColor(0);//10
  can3->SetBorderMode(0);
  can3->SetBorderSize(2);
  can3->SetFrameFillColor(0);
  can3->SetFrameBorderMode(0);
  can3->SetFrameBorderMode(0);
  can3->cd();
  TPad *pad3 = new TPad("pad3","pad3",0.0,0.0,1.,1.);
  gPad->SetGridx(0);
  gPad->SetGridy(0);
  gPad->SetTickx();
  gPad->SetTicky();
  pad3->SetTopMargin(0.02);//0.05
  pad3->SetRightMargin(0.01);//3e-2
  pad3->SetBottomMargin(0.07);//0.12
  //pad3->Divide(1,2,0,0);
  pad3->Draw();
  pad3->cd();
  TLegend *legend4 = new TLegend(.24,.78, .85,.97,NULL,"brNDC");//.45 or .4 for x1
  legend4->SetBorderSize(0);
  legend4->SetFillColor(0);
  
  gPad->SetRightMargin(0.01); gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.11); gPad->SetTopMargin(0.01);

  TH1D *C2_os_ktcompare[2][KTBINS];
  TH1D *C2Therm_os_ktcompare[2][KTBINS];
  float MIN_ktcomp=0.985, MAX_ktcomp=1.045;// MIN_ktcomp=0.92, MAX_ktcomp=2.04 (cumulant comparison)
  for(int type=0; type<2; type++){
    //pad3->cd(type+1);
    //gPad->SetRightMargin(RightMargin); gPad->SetLeftMargin(0.12);
    //gPad->SetBottomMargin(0.145);
    for(int kt=0; kt<KTBINS; kt++){
      C2_os_ktcompare[type][kt] = (TH1D*)C2_os[kt][BOI_1]->Clone();
      TH1D *FSIbase;
      if(type==0) {
	FSIbase = (TH1D*)C2Therm_os[0][BOI_1]->Clone();// diluted
	C2Therm_os_ktcompare[type][kt] = (TH1D*)C2Therm_os[kt][BOI_1]->Clone();// diluted
      }else{
	FSIbase = (TH1D*)K2_os[0][BOI_1]->Clone();// undiluted
	C2Therm_os_ktcompare[type][kt] = (TH1D*)K2_os[kt][BOI_1]->Clone();// undiluted
      }
      //

      for(int qbin=2; qbin<=20; qbin++){
	if(C2_os[0][BOI_1]->GetBinContent(qbin)<=0) continue;
	if(C2Therm_os[0][BOI_1]->GetBinContent(qbin)<=0) continue;
	
	double value = C2_os_ktcompare[type][kt]->GetBinContent(qbin)/C2_os[0][BOI_1]->GetBinContent(qbin);
	double value_e = pow(C2_os_ktcompare[type][kt]->GetBinError(qbin)/C2_os[0][BOI_1]->GetBinContent(qbin),2);
	value_e += pow(C2_os[0][BOI_1]->GetBinError(qbin) * C2_os_ktcompare[type][kt]->GetBinContent(qbin)/pow(C2_os[0][BOI_1]->GetBinContent(qbin),2),2);
	value_e = sqrt(value_e);
	double valueC2Therm = C2Therm_os_ktcompare[type][kt]->GetBinContent(qbin)/FSIbase->GetBinContent(qbin);
	double valueC2Therm_e = pow(C2Therm_os_ktcompare[type][kt]->GetBinError(qbin)/FSIbase->GetBinContent(qbin),2);
	valueC2Therm_e += pow(FSIbase->GetBinError(qbin) * C2Therm_os_ktcompare[type][kt]->GetBinContent(qbin)/pow(FSIbase->GetBinContent(qbin),2),2);
	valueC2Therm_e = sqrt(valueC2Therm_e);
	// cumulant comparison
	//double value = (C2_os_ktcompare[type][kt]->GetBinContent(qbin)-1.)/(C2_os[0][BOI_1]->GetBinContent(qbin)-1.);
	//double value_e = pow(C2_os_ktcompare[type][kt]->GetBinError(qbin)/(C2_os[0][BOI_1]->GetBinContent(qbin)-1.),2);
	//value_e += pow(C2_os[0][BOI_1]->GetBinError(qbin) * (C2_os_ktcompare[type][kt]->GetBinContent(qbin)-1.)/pow((C2_os[0][BOI_1]->GetBinContent(qbin)-1.),2),2);
	//value_e = sqrt(value_e);
	//double valueC2Therm = (C2Therm_os_ktcompare[type][kt]->GetBinContent(qbin)-1.)/(FSIbase->GetBinContent(qbin)-1.);
	//double valueC2Therm_e = pow(C2Therm_os_ktcompare[type][kt]->GetBinError(qbin)/(FSIbase->GetBinContent(qbin)-1.),2);
	//valueC2Therm_e += pow(FSIbase->GetBinError(qbin) * (C2Therm_os_ktcompare[type][kt]->GetBinContent(qbin)-1.)/pow((FSIbase->GetBinContent(qbin)-1.),2),2);
	//valueC2Therm_e = sqrt(valueC2Therm_e);
	
	
	//
	C2_os_ktcompare[type][kt]->SetBinContent(qbin, value);
	C2_os_ktcompare[type][kt]->SetBinError(qbin, value_e);
	C2Therm_os_ktcompare[type][kt]->SetBinContent(qbin, valueC2Therm);
	if(type==0) C2Therm_os_ktcompare[type][kt]->SetBinError(qbin, valueC2Therm_e);
	else C2Therm_os_ktcompare[type][kt]->SetBinError(qbin, C2Therm_os_ktcompare[0][kt]->GetBinError(qbin));
	if(kt==0) {C2_os_ktcompare[type][kt]->SetBinError(qbin, 0);}
      }
      C2_os_ktcompare[type][kt]->SetBinContent(1,-10); C2Therm_os_ktcompare[type][kt]->SetBinContent(1,-10);
      C2_os_ktcompare[type][kt]->SetMinimum(MIN_ktcomp);
      C2_os_ktcompare[type][kt]->SetMaximum(MAX_ktcomp);
      C2_os_ktcompare[type][kt]->GetXaxis()->SetRangeUser(0,0.09);
      int MarkerStyle=20;
      if(kt==1) MarkerStyle=24;
      if(kt==5) MarkerStyle=21;
      C2_os_ktcompare[type][kt]->SetMarkerStyle(MarkerStyle);
      C2Therm_os_ktcompare[type][kt]->SetMarkerStyle(24);// MarkerStyle
      C2Therm_os_ktcompare[type][kt]->SetMarkerColor(1); C2Therm_os_ktcompare[type][kt]->SetLineColor(1); 
      C2Therm_os_ktcompare[type][kt]->SetMarkerSize(C2_os_ktcompare[type][kt]->GetMarkerSize());
      C2_os_ktcompare[type][kt]->GetXaxis()->SetTitleOffset(20); C2_os_ktcompare[type][kt]->GetYaxis()->SetTitleOffset(20);
      
      C2_os_ktcompare[type][kt]->GetXaxis()->SetNdivisions(606);

      //if(kt!=5 && kt!=1 ) continue;// old way
      if(kt!=5) continue;
      if(kt==5) {
	C2_os_ktcompare[type][kt]->GetXaxis()->SetTitleFont(TextFont); C2_os_ktcompare[type][kt]->GetYaxis()->SetTitleFont(TextFont);
	C2_os_ktcompare[type][kt]->GetXaxis()->SetTitleSize(SizeTitle*SF2); 
	C2_os_ktcompare[type][kt]->GetYaxis()->SetTitleSize(SizeTitle*SF2);
	C2_os_ktcompare[type][kt]->GetXaxis()->SetTitle("#font[12]{q} (GeV/#font[12]{c})");
	C2_os_ktcompare[type][kt]->GetXaxis()->SetTitleOffset(0.8);
	C2_os_ktcompare[type][kt]->GetXaxis()->SetLabelFont(TextFont); C2_os_ktcompare[type][kt]->GetYaxis()->SetLabelFont(TextFont);
	C2_os_ktcompare[type][kt]->GetXaxis()->SetLabelSize(SizeLabel*SF2); 
	C2_os_ktcompare[type][kt]->GetYaxis()->SetLabelSize(SizeLabel*SF2); 
	if(type==1) C2Therm_os_ktcompare[type][kt]->SetMarkerStyle(25);

	if(type==0) C2_os_ktcompare[type][kt]->DrawCopy(); 
	C2Therm_os_ktcompare[type][kt]->DrawCopy("same");
      }
      //else {
      //C2_os_ktcompare[type][kt]->DrawCopy("same");
      //}
      TString *ktname = new TString("#font[12]{k}_{T,");
      *ktname += kt+1;
      ktname->Append("}/#font[12]{k}_{T,1}");
      TString *nameReal=new TString(ktname->Data());
      nameReal->Append(", ALICE");
      TString *nameTherm=new TString(ktname->Data());
      nameTherm->Append(", Therminator");
      if(type==0) {
	nameTherm->Append(" (diluted)");
	legend4->AddEntry(C2_os_ktcompare[type][kt],nameReal->Data(),"p");
	legend4->AddEntry(C2Therm_os_ktcompare[type][kt],nameTherm->Data(),"p");
      }else {
	TString *nameTherm2=new TString(nameTherm->Data());
	//nameTherm2->Append(", r*<80 fm");
	nameTherm2->Append(" (undiluted)");
	legend4->AddEntry(C2Therm_os_ktcompare[type][kt],nameTherm2->Data(),"p");
      }
      
    }
    
    if(type==1) {
      TLatex *Specif_3 = new TLatex(0.4,0.15,"ALICE Pb-Pb #\sqrt{#font[12]{s}_{NN}}=2.76 TeV");
      Specif_3->SetNDC();
      Specif_3->SetTextFont(TextFont);
      Specif_3->SetTextSize(SizeSpecif*SF2);
      TLatex *Specif2_3 = new TLatex(0.55,0.22,"0-5%");
      Specif2_3->SetNDC();
      Specif2_3->SetTextFont(TextFont);
      Specif2_3->SetTextSize(SizeLegend*SF2);
      Specif_3->Draw();
      Specif2_3->Draw();
      legend4->SetTextFont(TextFont);
      legend4->SetTextSize(SizeLegend*SF2);
      legend4->Draw("same");
    }

    TF1 *Unity = new TF1("Unity","1",0,100);
    Unity->SetLineStyle(2);
    Unity->Draw("same");
  
  }
  can3->cd();
  TPad *pad3_2 = new TPad("pad3_2","pad3_2",0.0,0.0,1.,1.);
  pad3_2->SetFillStyle(0);
  pad3_2->Draw();
  pad3_2->cd();
  TLatex *DthLabel = new TLatex(.044,.58,"#font[12]{C}_{2}^{+-}(#font[12]{k}_{T,6}) / #font[12]{C}_{2}^{+-}(#font[12]{k}_{T,1})");// .04,.96, "D^{+-}_{Th}"
  DthLabel->SetTextFont(TextFont);
  DthLabel->SetTextSize(SizeTitle*SF2);
  DthLabel->SetTextAngle(90);
  DthLabel->Draw();
  

  ////////////////////////////////////////////////////
  ////////////////////////////////////////////////////
  // Three-pion correlations

  
  TCanvas *can4 = new TCanvas("can4", "can4",10,300,600,900);// 11,53,700,500
  can4->SetHighLightColor(2);
  gStyle->SetOptFit(0111);
  can4->SetFillColor(0);//10
  can4->SetBorderMode(0);
  can4->SetBorderSize(2);
  can4->SetFrameFillColor(0);
  can4->SetFrameBorderMode(0);
  can4->SetFrameBorderMode(0);
  
  
  TPad *pad4 = new TPad("pad4","pad4",0.06,0.06,1.,1.);
  gPad->SetGridx(0);
  gPad->SetGridy(0);
  gPad->SetTickx();
  gPad->SetTicky();
  pad4->SetTopMargin(0.02);//0.05
  pad4->SetLeftMargin(0.12);
  pad4->SetRightMargin(0.01);//3e-2
  pad4->SetBottomMargin(0.08);//0.12
  pad4->Divide(2,3,0,0);
  pad4->Draw();


  TLegend *legend5 = new TLegend(.55,.6,.99,.99,NULL,"brNDC");//.55,.6,.99,.95
  legend5->SetBorderSize(0);
  legend5->SetTextFont(TextFont);
  legend5->SetTextSize(SizeLegend);
  legend5->SetFillColor(0);


  TF1 *BaseLine = new TF1("BaseLine","1.0",0,1);
  BaseLine->SetLineStyle(3);
  BaseLine->SetLineColor(1);


  TF1 *GaussFit[6];
  for(int ii=0; ii<6; ii++) {
    TString *nameFit=new TString("GaussFit_");
    *nameFit += ii;
    GaussFit[ii] = new TF1(nameFit->Data(),"[0]+[1]*exp(-pow([2]*x,2))",0,0.2);
    GaussFit[ii]->SetParameter(0, 1.);
    GaussFit[ii]->SetParameter(1, 2.);
    GaussFit[ii]->SetParameter(2, 10./FmToGeV);
  }
  
  TF1 *C3_SysFit=new TF1("C3_SysFit","sqrt(pow(pol3,2)+pow(0.003,2))",0,0.15);// MRC + PID
  C3_SysFit->FixParameter(0, -0.01354);
  C3_SysFit->FixParameter(1, 0.2948);
  C3_SysFit->FixParameter(2, -2.237);
  C3_SysFit->FixParameter(3, 5.697);
  TH1D *C3Sys[6];// cb
  TH1D *c3Sys[6];// cb
  TF1 *c3_Coul_SysFit=new TF1("c3_Coul_SysFit","pol1",0,0.15);
  TH1D *c3Sys_Coul[6];
  float C3MAX=2.3;
  if(ChProdBOI==1) C3MAX=2.1;

  for(int cb=0; cb<6; cb++) {
    pad4->cd(cb+1);
    gPad->SetRightMargin(RightMargin);
    if(cb%2==0) gPad->SetRightMargin(0.001);
    if(cb==4) SF_correction=0.92;
    else if(cb==5) SF_correction=1.015;
    else SF_correction=1.0;
    //if(cb==0) Dim1=gPad->GetAbsWNDC()*gPad->GetAbsHNDC()*gPad->PixeltoX(10)*gPad->PixeltoY(10);
    //Dim2=gPad->GetAbsWNDC()*gPad->GetAbsHNDC()*gPad->PixeltoX(10)*gPad->PixeltoY(10);
    
    C3merged[0][CoulChoice][KT3Bin][cb]->SetMarkerSize(1.5); C3merged[1][CoulChoice][KT3Bin][cb]->SetMarkerSize(1.5); c3merged[0][CoulChoice][KT3Bin][cb]->SetMarkerSize(1.5); c3merged[1][CoulChoice][KT3Bin][cb]->SetMarkerSize(1.5);
    C3merged[0][CoulChoice][KT3Bin][cb]->SetMarkerStyle(20); C3merged[1][CoulChoice][KT3Bin][cb]->SetMarkerStyle(21);
    c3merged[0][CoulChoice][KT3Bin][cb]->SetMarkerStyle(24); c3merged[1][CoulChoice][KT3Bin][cb]->SetMarkerStyle(25);
    C3EW[CoulChoice][KT3Bin][cb]->SetMarkerStyle(22);
    C3EW[CoulChoice][KT3Bin][cb]->SetMarkerColor(1);
    C3EW[CoulChoice][KT3Bin][cb]->SetLineColor(1);
    C3merged[0][CoulChoice][KT3Bin][cb]->SetMarkerColor(2); C3merged[0][CoulChoice][KT3Bin][cb]->SetLineColor(2); c3merged[0][CoulChoice][KT3Bin][cb]->SetMarkerColor(2); c3merged[0][CoulChoice][KT3Bin][cb]->SetLineColor(2); 
    C3merged[1][CoulChoice][KT3Bin][cb]->SetMarkerColor(4); C3merged[1][CoulChoice][KT3Bin][cb]->SetLineColor(4); c3merged[1][CoulChoice][KT3Bin][cb]->SetMarkerColor(4); c3merged[1][CoulChoice][KT3Bin][cb]->SetLineColor(4); 
    C3merged[ChProdBOI][CoulChoice][KT3Bin][cb]->GetXaxis()->SetTitle(""); C3merged[ChProdBOI][CoulChoice][KT3Bin][cb]->GetYaxis()->SetTitle("");
    C3merged[ChProdBOI][CoulChoice][KT3Bin][cb]->GetXaxis()->SetLabelFont(TextFont);
    C3merged[ChProdBOI][CoulChoice][KT3Bin][cb]->GetYaxis()->SetLabelFont(TextFont);
    C3merged[ChProdBOI][CoulChoice][KT3Bin][cb]->GetXaxis()->SetLabelSize(SizeLabel*SF_correction);
    if(cb%2==0) C3merged[ChProdBOI][CoulChoice][KT3Bin][cb]->GetYaxis()->SetLabelSize(SizeLabel*SF_correction);
    else {C3merged[ChProdBOI][CoulChoice][KT3Bin][cb]->GetYaxis()->SetLabelSize(.0);}
    if(cb<4) C3merged[ChProdBOI][CoulChoice][KT3Bin][cb]->GetXaxis()->SetLabelSize(.0);
    C3merged[ChProdBOI][CoulChoice][KT3Bin][cb]->GetXaxis()->SetRangeUser(0,0.11);
    C3merged[ChProdBOI][CoulChoice][KT3Bin][cb]->SetMinimum(0.9);//0.9
    C3merged[ChProdBOI][CoulChoice][KT3Bin][cb]->SetMaximum(C3MAX);
    if(ChProdBOI==0) c3merged[ChProdBOI][CoulChoice][KT3Bin][cb]->Fit(GaussFit[cb],"IMEQN","",0.02,0.1);
    //
    if(cb>2 && ChProdBOI==0) {
      C3merged[ChProdBOI][CoulChoice][KT3Bin][cb]->SetBinContent(2,-100);
      c3merged[ChProdBOI][CoulChoice][KT3Bin][cb]->SetBinContent(2,-100);
    }
    C3merged[ChProdBOI][CoulChoice][KT3Bin][cb]->GetXaxis()->SetNdivisions(404);
    C3merged[ChProdBOI][CoulChoice][KT3Bin][cb]->GetYaxis()->SetNdivisions(505);
    //TGaxis::SetMaxDigits(3);
    C3merged[ChProdBOI][CoulChoice][KT3Bin][cb]->DrawCopy();
    c3merged[ChProdBOI][CoulChoice][KT3Bin][cb]->DrawCopy("same");
        
    C3Sys[cb] = (TH1D*)C3merged[ChProdBOI][CoulChoice][KT3Bin][cb]->Clone();
    c3Sys[cb] = (TH1D*)c3merged[ChProdBOI][CoulChoice][KT3Bin][cb]->Clone();
    c3Sys_Coul[cb] = (TH1D*)c3merged[ChProdBOI][CoulChoice][KT3Bin][cb]->Clone();// Coulomb choice
    c3Sys_Coul[cb]->Add(c3merged[ChProdBOI][1-int(CoulChoice)][KT3Bin][cb],-1);// Alternate Coulomb
    c3Sys_Coul[cb]->Fit(c3_Coul_SysFit,"MENQ","",0,0.1);
    for(int ii=1; ii<=11; ii++){
      double Q3=(ii-0.5)*0.01;
      C3Sys[cb]->SetBinError(ii, C3_SysFit->Eval(Q3) * C3merged[ChProdBOI][CoulChoice][KT3Bin][cb]->GetBinContent(ii));// MRC + PID
      //
      double err = pow(0.2 * C3_SysFit->Eval(Q3) * c3merged[ChProdBOI][CoulChoice][KT3Bin][cb]->GetBinContent(ii),2);// MRC + PID (1/5 that of C3)
      err += pow(c3_Coul_SysFit->Eval(Q3),2);// Coulomb
      c3Sys[cb]->SetBinError(ii, sqrt(err));
    }
    C3Sys[cb]->SetMarkerSize(0); c3Sys[cb]->SetMarkerSize(0);
    C3Sys[cb]->SetFillStyle(1000); c3Sys[cb]->SetFillStyle(1000);
    if(ChProdBOI==0) { C3Sys[cb]->SetFillColor(kRed-9); c3Sys[cb]->SetFillColor(kRed-9);}
    else {C3Sys[cb]->SetFillColor(kBlue-9); c3Sys[cb]->SetFillColor(kBlue-9);}
    C3Sys[cb]->Draw("E2 same");
    c3Sys[cb]->Draw("E2 same");
    C3merged[ChProdBOI][CoulChoice][KT3Bin][cb]->DrawCopy("same");// redraw
    c3merged[ChProdBOI][CoulChoice][KT3Bin][cb]->DrawCopy("same");// redraw
    BaseLine->Draw("same");
    //ResidueFit->Draw("same");
   

    if(cb==0){
      if(ChProdBOI==0){
	legend5->AddEntry(C3merged[ChProdBOI][CoulChoice][KT3Bin][cb],"#font[12]{C}_{3}^{#pm#pm#pm}","p");
	legend5->AddEntry(c3merged[ChProdBOI][CoulChoice][KT3Bin][cb],"#font[12]{#bf{c}}_{3}^{#pm#pm#pm}","p");
      }else{
	legend5->AddEntry(C3merged[ChProdBOI][CoulChoice][KT3Bin][cb],"#font[12]{C}_{3}^{#pm#pm#mp}","p");
	legend5->AddEntry(c3merged[ChProdBOI][CoulChoice][KT3Bin][cb],"#font[12]{#bf{c}}_{3}^{#pm#pm#mp}","p");
      }
      legend5->Draw("same");
    }
    TLatex *Specif_3 = new TLatex(0.007,0.9+(C3MAX-0.9)*0.91,"ALICE Pb-Pb #\sqrt{#font[12]{s}_{NN}}=2.76 TeV");// ALICE specifications
    Specif_3->SetTextFont(TextFont);
    Specif_3->SetTextSize(SizeSpecif*SF_correction);
    if(cb==3) Specif_3->Draw("same");

    TString *KTString = new TString("");
    if(KT3Bin==0) KTString->Append("0.16<#font[12]{K}_{T,3}<0.3 GeV/#font[12]{c}");
    else KTString->Append("0.3<#font[12]{K}_{T,3}<1.0 GeV/#font[12]{c}");
    TLatex *Kt3Name_1 = new TLatex(0.016,0.9+(C3MAX-0.9)*0.91,KTString->Data());
    Kt3Name_1->SetTextFont(TextFont);
    Kt3Name_1->SetTextSize(SizeLabel*SF_correction);
    if(cb==1) Kt3Name_1->Draw("same");

    TLatex *CentLabel_1;
    if(cb==0) CentLabel_1 = new TLatex(0.07,C3MAX-0.75,"0-5%");// ALICE specifications
    else if(cb==1) CentLabel_1 = new TLatex(0.07,C3MAX-0.75,"5-10%");
    else if(cb==2) CentLabel_1 = new TLatex(0.07,C3MAX-0.75,"10-20%");
    else if(cb==3) CentLabel_1 = new TLatex(0.07,C3MAX-0.75,"20-30%");
    else if(cb==4) CentLabel_1 = new TLatex(0.07,C3MAX-0.75,"30-40%");
    else CentLabel_1 = new TLatex(0.07,C3MAX-0.75,"40-50%");
    
    CentLabel_1->SetTextFont(TextFont);
    CentLabel_1->SetTextSize(SizeLabel*SF_correction);//.1
    CentLabel_1->Draw("same");

    if(ChProdBOI==1){
      TF1 *ResidueFit=new TF1("ResidueFit","pol0",0,0.1);
      c3merged[ChProdBOI][CoulChoice][KT3Bin][cb]->Fit(ResidueFit,"NME","",0,0.1);
      double Chi2=0;
      for(int q3bin=2; q3bin<=10; q3bin++){
	Chi2 += pow(c3merged[ChProdBOI][CoulChoice][KT3Bin][cb]->GetBinContent(q3bin)-1.0,2)/pow(c3merged[ChProdBOI][CoulChoice][KT3Bin][cb]->GetBinError(q3bin),2);
      }
      cout<<"Chi2/NDF MC from baseline: "<<Chi2/9.<<endl;
    }
    
  }
  
  can4->cd();
  TPad *pad4_2 = new TPad("pad4_2","pad4_2",0.0,0.0,1.,1.);
  pad4_2->SetFillStyle(0);
  pad4_2->Draw();
  pad4_2->cd();
  TLatex *C3Label = new TLatex(.04,.86,"#font[12]{C}_{3} or #font[12]{#bf{c}}_{3}");// 0.05,0.9
  C3Label->SetTextFont(TextFont);
  C3Label->SetTextSize(SizeTitle*SF2);
  C3Label->SetTextAngle(90);
  C3Label->Draw();
  TLatex *Q3Label = new TLatex(.72,.02,"#font[12]{Q}_{3} (GeV/#font[12]{c})");// 0.05,0.9
  Q3Label->SetTextFont(TextFont);
  Q3Label->SetTextSize(SizeTitle*SF2);
  Q3Label->Draw();
  CoverUp->Draw();

  if(ChProdBOI==1) return;

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  // r3

  TF1 *ChaoticLimit = new TF1("ChaoticLimit","2.0",0,1);
  ChaoticLimit->SetLineStyle(3);
  
  TF1 *QuarticFit[4][6];// Stat/Det/Coul1/Coul2, cb
  TF1 *QuadraticFit[4][6];// Stat/Det/Coul1/Coul2, cb
  TF1 *QuadQuarticFit[4][6];// Stat/Det/Coul1/Coul2, cb
  for(int ii=0; ii<6; ii++){
    for(int jj=0; jj<4; jj++){
      TString *nameFit=new TString("Quartic_");
      *nameFit += ii;
      nameFit->Append("_E");
      *nameFit += jj;
      QuarticFit[jj][ii] = new TF1(nameFit->Data(),"[0]*(1-[1]*pow(x,4))",0,0.2);
      QuarticFit[jj][ii]->SetParameter(0,2.0);
      QuarticFit[jj][ii]->SetParameter(1,0);
      QuarticFit[jj][ii]->SetLineColor(4);
      TString *nameFit2=new TString("Quadratic_");
      *nameFit2 += ii;
      nameFit2->Append("_E");
      *nameFit2 += jj;
      QuadraticFit[jj][ii] = new TF1(nameFit2->Data(),"[0]*(1-[1]*pow(x,2))",0,0.2);
      QuadraticFit[jj][ii]->SetParameter(0,2.0);
      QuadraticFit[jj][ii]->SetParameter(1,0);
      QuadraticFit[jj][ii]->SetLineColor(1);
      QuadraticFit[jj][ii]->SetLineStyle(4);
      TString *nameFit3=new TString("QuadQuartic_");
      *nameFit3 += ii;
      nameFit3->Append("_E");
      *nameFit3 += jj;
      QuadQuarticFit[jj][ii] = new TF1(nameFit3->Data(),"[0]*(1-[1]*pow(x,2)-[2]*pow(x,4))",0,0.2);
      QuadQuarticFit[jj][ii]->SetParameter(0,2.0);
      QuadQuarticFit[jj][ii]->SetParameter(1,50);
      QuadQuarticFit[jj][ii]->SetParameter(2,5000);
      QuadQuarticFit[jj][ii]->SetLineColor(6);
      QuadQuarticFit[jj][ii]->SetLineStyle(3);

    }
  }
  
  TCanvas *can5 = new TCanvas("can5", "can5",680,300,600,900);// 11,53,700,500
  can5->SetHighLightColor(2);
  gStyle->SetOptFit(0111);
  can5->SetFillColor(0);//10
  can5->SetBorderMode(0);
  can5->SetBorderSize(2);
  can5->SetFrameFillColor(0);
  can5->SetFrameBorderMode(0);
  can5->SetFrameBorderMode(0);
  
  
  TPad *pad5 = new TPad("pad5","pad5",0.06,0.06,1.,1.);
  gPad->SetGridx(0);
  gPad->SetGridy(0);
  gPad->SetTickx();
  gPad->SetTicky();
  pad5->SetTopMargin(0.02);//0.05
  pad5->SetLeftMargin(0.12);//.14 for wide title, .10 for narrow title, 0.08 for DeltaQ
  pad5->SetRightMargin(0.01);//3e-2
  pad5->SetBottomMargin(0.08);//0.12
  pad5->Divide(2,3,0,0);
  pad5->Draw();
  
  TLegend *legend6 = new TLegend(.18,.72, .53,.98,NULL,"brNDC");// .25,.16, .6,.46
  legend6->SetBorderSize(0);
  legend6->SetTextFont(TextFont);
  legend6->SetTextSize(SizeLegend);
  legend6->SetFillColor(0);
  
  
  double LambdaSysPar0[6]={0.0143, 0.014, 0.01294, 0.00204, 0.001227, 0.00999};
  double LambdaSysPar1[6]={-0.857, -0.741, -0.6655, -0.3516, -0.2327, -0.3325};
  double r3MIN = 1.45;// 0.3
  double r3MAX = 2.45;// 2.68
  r3merged[CoulChoice][1][5]->SetBinContent(3,100);  r3merged[CoulChoice][1][5]->SetBinError(3,1000);// remove high stat error point
 
  for(int cb=0; cb<6; cb++) {
    pad5->cd(cb+1);
    gPad->SetRightMargin(RightMargin);
    if(cb%2==0) gPad->SetRightMargin(0.001);
    if(cb==4) SF_correction=0.92;
    else if(cb==5) SF_correction=1.015;
    else SF_correction=1.0;
    r3merged[CoulChoice][KT3Bin][cb]->SetMarkerSize(1.5);
    r3merged[CoulChoice][KT3Bin][cb]->SetMarkerStyle(20);
    r3merged[CoulChoice][KT3Bin][cb]->SetMarkerColor(2); r3merged[CoulChoice][KT3Bin][cb]->SetLineColor(2); 
    r3merged[CoulChoice][KT3Bin][cb]->GetXaxis()->SetTitle(""); 
    r3merged[CoulChoice][KT3Bin][cb]->GetYaxis()->SetTitle(""); 
    r3merged[CoulChoice][KT3Bin][cb]->GetXaxis()->SetLabelFont(TextFont);
    r3merged[CoulChoice][KT3Bin][cb]->GetYaxis()->SetLabelFont(TextFont);
    r3merged[CoulChoice][KT3Bin][cb]->GetXaxis()->SetLabelSize(SizeLabel*SF_correction);
    if(cb%2==0) r3merged[CoulChoice][KT3Bin][cb]->GetYaxis()->SetLabelSize(SizeLabel*SF_correction);
    else {r3merged[CoulChoice][KT3Bin][cb]->GetYaxis()->SetLabelSize(.0);}
    if(cb<4) r3merged[CoulChoice][KT3Bin][cb]->GetXaxis()->SetLabelSize(.0);
    r3merged[CoulChoice][KT3Bin][cb]->GetXaxis()->SetRangeUser(0,0.071);// was 0,0.11
    r3merged[CoulChoice][KT3Bin][cb]->SetMinimum(r3MIN);
    r3merged[CoulChoice][KT3Bin][cb]->SetMaximum(r3MAX);
    if(cb>2 || KT3Bin==1) r3merged[CoulChoice][KT3Bin][cb]->SetBinContent(2,-100);// 10-20 MeV bin is insignificant
    r3merged[CoulChoice][KT3Bin][cb]->GetXaxis()->SetNdivisions(502);// was 404
    r3merged[CoulChoice][KT3Bin][cb]->GetYaxis()->SetNdivisions(505);// was absent
    r3merged[CoulChoice][KT3Bin][cb]->DrawCopy();
    //ChaoticLimit->Draw("same");

    ///////////////
    // Systematics
    TF1 *ResidueFit=new TF1("ResidueFit","[0]+[1]*exp(-[2]*x)",0,0.2);
    ResidueFit->SetParameter(0,1);
    ResidueFit->SetParameter(1,0.001);
    ResidueFit->SetParameter(2,1);
    ResidueFit->SetParLimits(0,.9,1.1); ResidueFit->SetParLimits(1,-.2,.2); ResidueFit->SetParLimits(2,-40,200);
    c3merged[1][CoulChoice][KT3Bin][cb]->Fit(ResidueFit,"IMNQ","",0.01,0.14);// Remove "N" to see residue fits
    //cout<<setprecision(5)<<ResidueFit->Eval(0.08)<<endl;
    TH1D *r3MethSys = (TH1D*)r3merged[CoulChoice][KT3Bin][cb]->Clone();
    r3MethSys->SetBinContent(1,-10);
    r3MethSys->SetMarkerSize(0);
    r3MethSys->SetFillStyle(1000);
    r3MethSys->SetFillColor(kRed-10);
    TH1D *r3DetSys = (TH1D*)r3merged[CoulChoice][KT3Bin][cb]->Clone();
    r3DetSys->SetBinContent(1,-10);
    r3DetSys->SetBinContent(1,-10);
    r3DetSys->SetMarkerSize(0);
    r3DetSys->SetFillStyle(1000);
    r3DetSys->SetFillColor(kGray+2);
    TH1D *r3MethDetSys = (TH1D*)r3merged[CoulChoice][KT3Bin][cb]->Clone();// for fitting systematics
    TH1D *r3Coul1 = (TH1D*)r3merged[CoulChoice][KT3Bin][cb]->Clone();
    r3Coul1->SetMarkerSize(0);
    TH1D *r3Coul2 = (TH1D*)r3merged[1-int(CoulChoice)][KT3Bin][cb]->Clone();
    
    TF1 *r3LambdaSysFit = new TF1("r3LambdaSysFit","pol1",0,0.15);
    r3LambdaSysFit->FixParameter(0,LambdaSysPar0[cb]);
    r3LambdaSysFit->FixParameter(1,LambdaSysPar1[cb]);
    const int Nb = 20;
    double x[Nb]={0};
    double y[Nb]={0};
    double exL[Nb]={0};
    double exH[Nb]={0};
    double eyL[Nb]={0};
    double eyH[Nb]={0};
    for(int ii=0; ii<Nb; ii++){
      x[ii] = r3merged[CoulChoice][KT3Bin][cb]->GetXaxis()->GetBinCenter(ii+1);
      y[ii] = r3merged[CoulChoice][KT3Bin][cb]->GetBinContent(ii+1);
      double ey = r3merged[CoulChoice][KT3Bin][cb]->GetBinContent(ii+1)-r3merged[1-int(CoulChoice)][KT3Bin][cb]->GetBinContent(ii+1);
      if(ey>0) eyL[ii] = fabs(ey);
      else eyH[ii] = fabs(ey);
    }
    
    TGraphAsymmErrors *r3CoulSys = new TGraphAsymmErrors(Nb,x,y,exL,exH,eyL,eyH);
    r3CoulSys->SetLineStyle(2);
    for(int bin=1; bin<=11; bin++){
      double Q3 = (bin-0.5)*0.01;
      double SysMeth=0;
      if(C3merged[1][CoulChoice][KT3Bin][cb]->GetBinContent(bin)>0) {
	//SysMeth = pow(2 * (fabs(ResidueFit->Eval(Q3)-1.0) + 0.002)/(c3merged[0][CoulChoice][KT3Bin][cb]->GetBinContent(bin)-1.0),2);// mixed-charge baseline (new).  Includes 0.002 as allowed variation through lambda powers
	SysMeth = pow((fabs(ResidueFit->Eval(Q3)-1.0)+0.002)/(c3merged[0][CoulChoice][KT3Bin][cb]->GetBinContent(bin)-1.0) * r3merged[CoulChoice][KT3Bin][cb]->GetBinContent(bin),2);// factional error determined from MC residue fit divided by SC cumulant.  Error = Fraction * r3
      }
      SysMeth += pow(r3LambdaSysFit->Eval(Q3,2)*r3merged[CoulChoice][KT3Bin][cb]->GetBinContent(bin),2);// Lambda 0.7 vs 0.6
      double SysDet = pow(0.01*r3merged[CoulChoice][KT3Bin][cb]->GetBinContent(bin),2);// MRC
      SysDet += pow(0.01*r3merged[CoulChoice][KT3Bin][cb]->GetBinContent(bin),2);// PID
      
      r3MethSys->SetBinError(bin, sqrt(SysMeth));
      r3DetSys->SetBinError(bin, sqrt(SysDet));
      r3MethDetSys->SetBinError(bin, sqrt(pow(r3merged[CoulChoice][KT3Bin][cb]->GetBinError(bin),2) + SysMeth + SysDet));// was SysDet + SysDet by mistake
      r3Coul1->SetBinError(bin, sqrt(pow(r3merged[CoulChoice][KT3Bin][cb]->GetBinError(bin),2) + SysMeth + SysDet));
      r3Coul2->SetBinError(bin, sqrt(pow(r3merged[CoulChoice][KT3Bin][cb]->GetBinError(bin),2) + SysMeth + SysDet));
    }

    r3MethSys->Draw("E2 same");
    r3DetSys->Draw("E2 same");
    r3CoulSys->Draw("P");
    //
    float MinQ3Fit=0.01;//0.01
    float MaxQ3Fit=0.08;//0.1
    if(r3merged[CoulChoice][KT3Bin][cb]->GetBinContent(2)<-10) MinQ3Fit=0.02;

    r3merged[CoulChoice][KT3Bin][cb]->Fit(QuarticFit[0][cb],"IMENQ","",MinQ3Fit,MaxQ3Fit);// Stat only
    r3MethDetSys->Fit(QuarticFit[1][cb],"IMENQ","",MinQ3Fit,MaxQ3Fit);// Meth + Det Sys
    r3Coul1->Fit(QuarticFit[2][cb],"IMENQ","",MinQ3Fit,MaxQ3Fit);// GRS
    r3Coul2->Fit(QuarticFit[3][cb],"IMENQ","",MinQ3Fit,MaxQ3Fit);// Omega0
    //
    r3merged[CoulChoice][KT3Bin][cb]->Fit(QuadraticFit[0][cb],"IMENQ","",MinQ3Fit,MaxQ3Fit);// Stat only
    r3MethDetSys->Fit(QuadraticFit[1][cb],"IMENQ","",MinQ3Fit,MaxQ3Fit);// Meth + Det Sys
    r3Coul1->Fit(QuadraticFit[2][cb],"IMENQ","",MinQ3Fit,MaxQ3Fit);// GRS
    r3Coul2->Fit(QuadraticFit[3][cb],"IMENQ","",MinQ3Fit,MaxQ3Fit);// Omega0
    //
    //if(cb==1){
    //for(int i=0; i<10; i++){
    //	cout<<r3MethDetSys->GetBinError(i+1)<<endl;
    //}
    //}
    //
    QuarticFit[2][cb]->Draw("same");
    QuadraticFit[2][cb]->Draw("same");
    r3merged[CoulChoice][KT3Bin][cb]->Draw("same");
    //QuadQuarticFit[1][cb]->Draw("same");
    //cout<<"Quartic Chi^2/NDF = "<<QuarticFit[cb]->GetChisquare()/QuarticFit[cb]->GetNDF()<<endl;
    //cout<<"Quadratic Chi^2/NDF = "<<QuadraticFit[cb]->GetChisquare()/QuadraticFit[cb]->GetNDF()<<endl;
    
    TLatex *Specif_4 = new TLatex(0.007,r3MIN+(r3MAX-r3MIN)*0.91,"ALICE Pb-Pb #\sqrt{#font[12]{s}_{NN}}=2.76 TeV");// 0.005,2.1
    Specif_4->SetTextFont(TextFont);
    Specif_4->SetTextSize(SizeSpecif*SF_correction);
    if(cb==3) Specif_4->Draw("same");
   
    
    TString *KTString = new TString("");
    if(KT3Bin==0) KTString->Append("0.16<#font[12]{K}_{T,3}<0.3 GeV/#font[12]{c}");
    else KTString->Append("0.3<#font[12]{K}_{T,3}<1.0 GeV/#font[12]{c}");
    TLatex *Kt3Name_2 = new TLatex(0.002,r3MIN+(r3MAX-r3MIN)*0.91,KTString->Data());//
    Kt3Name_2->SetTextFont(TextFont);
    Kt3Name_2->SetTextSize(SizeLabel*SF_correction);
    if(cb==1) Kt3Name_2->Draw("same");
    
    TLatex *CentLabel_2;
    //TString *CentLabel2;
    if(cb==0) {CentLabel_2 = new TLatex(0.03,1.52,"0-5%"); }//CentLabel2 = new TString("0-5\\%");}// was 0.04,0.4
    else if(cb==1) {CentLabel_2 = new TLatex(0.03,1.52,"5-10%");}// CentLabel2 = new TString("5-10\\%");}
    else if(cb==2) {CentLabel_2 = new TLatex(0.03,1.52,"10-20%");}// CentLabel2 = new TString("10-20\\%");}
    else if(cb==3) {CentLabel_2 = new TLatex(0.03,1.52,"20-30%");}// CentLabel2 = new TString("20-30\\%");}
    else if(cb==4) {CentLabel_2 = new TLatex(0.03,1.52,"30-40%");}// CentLabel2 = new TString("30-40\\%");}
    else {CentLabel_2 = new TLatex(0.03,1.52,"40-50%");}// CentLabel2 = new TString("40-50\\%");}
    
    CentLabel_2->SetTextFont(TextFont);
    CentLabel_2->SetTextSize(SizeLabel*SF_correction);
    CentLabel_2->Draw("same");
    
    if(cb==0){
      legend6->AddEntry(QuarticFit[1][cb],"Quartic","l");
      legend6->AddEntry(QuadraticFit[1][cb],"Quadratic","l");
      legend6->Draw("same");
    }

    ChaoticLimit->Draw("same");
  }// cb bin
  can5->cd();
  TPad *pad5_2 = new TPad("pad5_2","pad5_2",0.0,0.0,1.,1.);
  pad5_2->SetFillStyle(0);
  pad5_2->Draw();
  pad5_2->cd();
  TLatex *r3Label = new TLatex(.035,.96,"#font[12]{r}_{3}");// 0.035,0.88
  r3Label->SetTextFont(TextFont);
  r3Label->SetTextSize(SizeTitle*SF2);
  r3Label->SetTextAngle(90);
  r3Label->Draw();
  CoverUp->Draw();
  Q3Label->Draw();
  
  
  ////////////////////////////
  // Print Quartic r3 fit values
  double I_quartic_Avg=0, I_quadratic_Avg=0;
  double I_quartic_AvgErr1=0, I_quadratic_AvgErr1=0;
  double I_quartic_AvgErr2=0, I_quadratic_AvgErr2=0;
  double a_quartic_Avg=0, a_quadratic_Avg=0;
  double a_quartic_AvgErr1=0, a_quadratic_AvgErr1=0;
  double a_quartic_AvgErr2=0, a_quadratic_AvgErr2=0;
  
  cout<<"Quartic"<<endl;
  for(int cb=0; cb<6; cb++) {
    TString *CentLabel2;
    if(cb==0) {CentLabel2 = new TString("0-5\\%");}
    else if(cb==1) {CentLabel2 = new TString("5-10\\%");}
    else if(cb==2) {CentLabel2 = new TString("10-20\\%");}
    else if(cb==3) {CentLabel2 = new TString("20-30\\%");}
    else if(cb==4) {CentLabel2 = new TString("30-40\\%");}
    else {CentLabel2 = new TString("40-50\\%");}

    double SF=1000.;
    double I_Sys1=fabs(QuarticFit[1][cb]->GetParameter(0)-QuarticFit[0][cb]->GetParameter(0));// Meth+Det Sys
    double I_Sys2=fabs(QuarticFit[2][cb]->GetParameter(0)-QuarticFit[3][cb]->GetParameter(0));// Coulomb Sys
    double I_SysTotal = sqrt(pow(I_Sys1,2)+pow(I_Sys2,2));
    double a_Sys1=fabs(QuarticFit[1][cb]->GetParameter(1)-QuarticFit[0][cb]->GetParameter(1));
    double a_Sys2=fabs(QuarticFit[2][cb]->GetParameter(1)-QuarticFit[3][cb]->GetParameter(1));
    a_Sys1 /= SF; a_Sys2 /= SF;
    double a_SysTotal = sqrt(pow(a_Sys1,2)+pow(a_Sys2,2));
    I_quartic_Avg += QuarticFit[1][cb]->GetParameter(0);// [1] has both Stat+Sys errors
    I_quartic_AvgErr1 += pow(QuarticFit[0][cb]->GetParError(0),2);// [0] has only Stat errors
    I_quartic_AvgErr2 += I_SysTotal;
    a_quartic_Avg += QuarticFit[1][cb]->GetParameter(1)/SF;
    a_quartic_AvgErr1 += pow(QuarticFit[0][cb]->GetParError(1)/SF,2);
    a_quartic_AvgErr2 += sqrt(pow(a_Sys1,2) + pow(a_Sys2,2));

    cout.precision(3);
    cout<<"$"<<setprecision(2)<<fixed<<CentLabel2->Data()<<"$ & $"<<QuarticFit[1][cb]->GetParameter(0)<<" \\pm "<<QuarticFit[0][cb]->GetParError(0)<<" \\pm "<<I_SysTotal<<"$ & $"<<setprecision(1)<<fixed<<QuarticFit[1][cb]->GetParameter(1)/SF<<" \\pm "<<QuarticFit[0][cb]->GetParError(1)/SF<<" \\pm "<<a_SysTotal<<"$ \\\\ \\hline"<<endl;
  }
  I_quartic_Avg /= 6.; I_quartic_AvgErr1 = sqrt(I_quartic_AvgErr1)/6.; I_quartic_AvgErr2 = I_quartic_AvgErr2/6.;
  a_quartic_Avg /= 6.; a_quartic_AvgErr1 = sqrt(a_quartic_AvgErr1)/6.; a_quartic_AvgErr2 = a_quartic_AvgErr2/6.;
  cout<<"$"<<setprecision(2)<<fixed<<"0-50\\%"<<"$ & $"<<I_quartic_Avg<<" \\pm "<<I_quartic_AvgErr1<<" \\pm "<<I_quartic_AvgErr2<<"$ & $"<<setprecision(1)<<fixed<<a_quartic_Avg<<" \\pm "<<a_quartic_AvgErr1<<" \\pm "<<a_quartic_AvgErr2<<"$ \\\\ \\hline"<<endl;
  
  /*TCanvas *can5_5 = new TCanvas("can5_5", "can5_5",680,300,600,600);// 11,53,700,500
  can5_5->SetHighLightColor(2);
  gStyle->SetOptFit(0111);
  can5_5->SetFillColor(0);//10
  can5_5->SetBorderMode(0);
  can5_5->SetBorderSize(2);
  can5_5->SetFrameFillColor(0);
  can5_5->SetFrameBorderMode(0);
  can5_5->SetFrameBorderMode(0);
  can5_5->Draw();
  double Iold_K1[6]={1.84,1.84,1.82,1.84,1.79,1.73};
  double Iold_K2[6]={1.94,1.92,2.02,1.98,2.01,2.03};
  double aold_K1[6]={13.8,10.6,6.9,6.7,3.7,4.9};
  double aold_K2[6]={8.7,4.6,4.4,1.8,-0.2,-1.8};

  TH1D *ratioMuonCorr=new TH1D("ratioMuonCorr","",6,0.5,6.5);
  ratioMuonCorr->GetXaxis()->SetBinLabel(1,"0-5%");
  ratioMuonCorr->GetXaxis()->SetBinLabel(2,"5-10%");
  ratioMuonCorr->GetXaxis()->SetBinLabel(3,"10-20%");
  ratioMuonCorr->GetXaxis()->SetBinLabel(4,"20-30%");
  ratioMuonCorr->GetXaxis()->SetBinLabel(5,"30-40%");
  ratioMuonCorr->GetXaxis()->SetBinLabel(6,"40-50%");
  ratioMuonCorr->GetYaxis()->SetTitle("New/Old"); ratioMuonCorr->GetYaxis()->SetTitleOffset(1.8);
  ratioMuonCorr->SetTitle("Quartic coefficient");
  for(int i=0; i<6; i++){
    ratioMuonCorr->SetBinContent(i+1, QuarticFit[1][i]->GetParameter(1)/SF/aold_K1[i]);
    ratioMuonCorr->SetBinError(i+1, 1000./QuarticFit[1][i]->GetParameter(1));// 0.015/1.8 for I_K1, 0.03/1.8 for I_K2, 
  }
  ratioMuonCorr->Draw();
  */
  
  ////////////////////////////
  cout<<"Quadratic"<<endl;
  for(int cb=0; cb<6; cb++) {
    TString *CentLabel2;
    if(cb==0) {CentLabel2 = new TString("0-5\\%");}
    else if(cb==1) {CentLabel2 = new TString("5-10\\%");}
    else if(cb==2) {CentLabel2 = new TString("10-20\\%");}
    else if(cb==3) {CentLabel2 = new TString("20-30\\%");}
    else if(cb==4) {CentLabel2 = new TString("30-40\\%");}
    else {CentLabel2 = new TString("40-50\\%");}
    double SF=10.;
    cout.precision(4);
    double I_Sys1 = fabs(QuadraticFit[1][cb]->GetParameter(0)-QuadraticFit[0][cb]->GetParameter(0));// Meth+Det Sys
    double I_Sys2 = fabs(QuadraticFit[2][cb]->GetParameter(0)-QuadraticFit[3][cb]->GetParameter(0));// Coulomb Sys
    double I_SysTotal = sqrt(pow(I_Sys1,2)+pow(I_Sys2,2));
    double a_Sys1=fabs(QuadraticFit[1][cb]->GetParameter(1)-QuadraticFit[0][cb]->GetParameter(1));
    double a_Sys2 = fabs(QuadraticFit[2][cb]->GetParameter(1)-QuadraticFit[3][cb]->GetParameter(1));
    a_Sys1 /= SF; a_Sys2 /= SF;
    double a_SysTotal = sqrt(pow(a_Sys1,2)+pow(a_Sys2,2));
    I_quadratic_Avg += QuadraticFit[1][cb]->GetParameter(0);
    I_quadratic_AvgErr1 += pow(QuadraticFit[0][cb]->GetParError(0),2);
    I_quadratic_AvgErr2 += I_SysTotal;
    a_quadratic_Avg += QuadraticFit[1][cb]->GetParameter(1)/SF;
    a_quadratic_AvgErr1 += pow(QuadraticFit[0][cb]->GetParError(1)/SF,2);
    a_quadratic_AvgErr2 += sqrt(pow(a_Sys1,2) + pow(a_Sys2,2));
    cout.precision(3);
    cout<<"$"<<setprecision(2)<<fixed<<CentLabel2->Data()<<"$ & $"<<QuadraticFit[1][cb]->GetParameter(0)<<" \\pm "<<QuadraticFit[0][cb]->GetParError(0)<<" \\pm "<<I_SysTotal<<"$ & $"<<setprecision(1)<<fixed<<QuadraticFit[1][cb]->GetParameter(1)/SF<<" \\pm "<<QuadraticFit[0][cb]->GetParError(1)/SF<<" \\pm "<<a_SysTotal<<"$ \\\\ \\hline"<<endl;
  }
  I_quadratic_Avg /= 6.; I_quadratic_AvgErr1 = sqrt(I_quadratic_AvgErr1)/6.; I_quadratic_AvgErr2 = I_quadratic_AvgErr2/6.;
  a_quadratic_Avg /= 6.; a_quadratic_AvgErr1 = sqrt(a_quadratic_AvgErr1)/6.; a_quadratic_AvgErr2 = a_quadratic_AvgErr2/6.;
  cout<<"$"<<setprecision(2)<<fixed<<"0-50\\%"<<"$ & $"<<I_quadratic_Avg<<" \\pm "<<I_quadratic_AvgErr1<<" \\pm "<<I_quadratic_AvgErr2<<"$ & $"<<setprecision(1)<<fixed<<a_quadratic_Avg<<" \\pm "<<a_quadratic_AvgErr1<<" \\pm "<<a_quadratic_AvgErr2<<"$ \\\\ \\hline"<<endl;
  
 


  
  ////////////////////////////////////////////////////////////////////////////////
  // Coulomb Variation
  
  int CBOI_SS=1;// 1=5-10%
  TCanvas *can6 = new TCanvas("can6", "can6",1800,300,600,900);// 680,300,600,600
  can6->SetHighLightColor(2);
  gStyle->SetOptFit(0111);
  can6->SetFillColor(0);//10
  can6->SetBorderMode(0);
  can6->SetBorderSize(2);
  can6->SetFrameFillColor(0);
  can6->SetFrameBorderMode(0);
  can6->SetFrameBorderMode(0);
  
  TPad *pad6 = new TPad("pad6","pad6",0.0,0.0,1.,1.);
  gPad->SetGridx(0);
  gPad->SetGridy(0);
  gPad->SetTickx();
  gPad->SetTicky();
  pad6->SetBottomMargin(0.07);//0.12
  pad6->SetTopMargin(0.02);//0.05
  pad6->SetRightMargin(0.01);//1e-2
  pad6->Divide(1,2,0,0);
  pad6->Draw();
  pad6->cd(1);
  TLegend *legend7 = new TLegend(.8,.05, .98,.2,NULL,"brNDC");//c3
  legend7->SetBorderSize(0);
  legend7->SetTextFont(TextFont);
  legend7->SetTextSize(SizeLegend*SF1);
  legend7->SetFillColor(0);
  
  //
  r3merged[1][KT3Bin][CBOI_SS]->SetMarkerStyle(24); r3merged[1][KT3Bin][CBOI_SS]->SetMarkerColor(2); r3merged[1][KT3Bin][CBOI_SS]->SetLineColor(2);// Omega0, Kt3_1
  r3merged[0][KT3Bin][CBOI_SS]->SetMinimum(1.1); r3merged[0][KT3Bin][CBOI_SS]->SetMaximum(2.35); 
  c3merged[1][1][KT3Bin][CBOI_SS]->SetMarkerStyle(21); c3merged[1][1][KT3Bin][CBOI_SS]->SetMarkerColor(4); c3merged[1][1][KT3Bin][CBOI_SS]->SetLineColor(4);// Omega0, Kt3_1
  c3merged[1][0][KT3Bin][CBOI_SS]->SetMinimum(0.995); c3merged[1][0][KT3Bin][CBOI_SS]->SetMaximum(1.065); 
  //
  gPad->SetRightMargin(RightMargin); gPad->SetLeftMargin(0.13);
  r3merged[0][KT3Bin][CBOI_SS]->GetXaxis()->SetLabelFont(TextFont); r3merged[0][KT3Bin][CBOI_SS]->GetYaxis()->SetLabelFont(TextFont); 
  r3merged[0][KT3Bin][CBOI_SS]->GetXaxis()->SetLabelSize(0); r3merged[0][KT3Bin][CBOI_SS]->GetYaxis()->SetLabelSize(SizeLabel*SF1);
  //
  r3merged[0][KT3Bin][CBOI_SS]->GetYaxis()->SetTitleOffset(0.8);
  r3merged[0][KT3Bin][CBOI_SS]->GetYaxis()->SetTitleFont(TextFont);
  r3merged[0][KT3Bin][CBOI_SS]->GetYaxis()->SetTitleSize(SizeTitle*SF1);
  r3merged[0][KT3Bin][CBOI_SS]->GetYaxis()->SetTitle("#font[12]{r}_{3}");
  r3merged[0][KT3Bin][CBOI_SS]->GetXaxis()->SetRangeUser(0,0.08);
  //
  r3merged[0][KT3Bin][CBOI_SS]->SetMinimum(0.9); r3merged[0][KT3Bin][CBOI_SS]->SetMaximum(2.4);
  r3merged[0][KT3Bin][CBOI_SS]->SetMarkerSize(1.5);
  r3merged[1][KT3Bin][CBOI_SS]->SetMarkerSize(1.7);
  r3merged[0][KT3Bin][CBOI_SS]->DrawCopy();
  r3merged[1][KT3Bin][CBOI_SS]->DrawCopy("same");
  ChaoticLimit->Draw("same");
  legend7->AddEntry(r3merged[0][KT3Bin][CBOI_SS],"GRS","p");
  legend7->AddEntry(r3merged[1][KT3Bin][CBOI_SS],"#Omega_{0}","p");
  legend7->Draw("same");
  TLatex *Specif_5 = new TLatex(0.36,0.94,"ALICE Pb-Pb #\sqrt{#font[12]{s}_{NN}}=2.76 TeV");
  Specif_5->SetNDC();
  Specif_5->SetTextFont(TextFont);
  Specif_5->SetTextSize(SizeSpecif*SF1);
  TLatex *Specif2_5 = new TLatex(0.5,0.1,"5-10%");
  Specif2_5->SetNDC();
  Specif2_5->SetTextFont(TextFont);
  Specif2_5->SetTextSize(SizeLabel*SF1);
  Specif_5->Draw();
  Specif2_5->Draw();
  TString *KTString = new TString("");
  if(KT3Bin==0) KTString->Append("0.16<#font[12]{K}_{T,3}<0.3 GeV/#font[12]{c}");
  else KTString->Append("0.3<#font[12]{K}_{T,3}<1.0 GeV/#font[12]{c}");
  TLatex *Kt3Name_3 = new TLatex(0.36,0.22,KTString->Data());
  Kt3Name_3->SetNDC();
  Kt3Name_3->SetTextFont(TextFont);
  Kt3Name_3->SetTextSize(SizeLabel*SF1);
  Kt3Name_3->Draw("same");
  
  //
  //
  pad6->cd(2);
  SF_correction=0.97;
  gPad->SetRightMargin(RightMargin); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.16);
  c3merged[1][0][KT3Bin][CBOI_SS]->GetXaxis()->SetLabelFont(TextFont); c3merged[1][0][KT3Bin][CBOI_SS]->GetYaxis()->SetLabelFont(TextFont);
  c3merged[1][0][KT3Bin][CBOI_SS]->GetXaxis()->SetLabelSize(SizeLabel*SF1*SF_correction); c3merged[1][0][KT3Bin][CBOI_SS]->GetYaxis()->SetLabelSize(SizeLabel*SF1*SF_correction);
  c3merged[1][0][KT3Bin][CBOI_SS]->SetMarkerStyle(21);
  c3merged[1][1][KT3Bin][CBOI_SS]->SetMarkerStyle(25);
  c3merged[1][0][KT3Bin][CBOI_SS]->SetMarkerSize(1.5);
  c3merged[1][1][KT3Bin][CBOI_SS]->SetMarkerSize(1.7);
  //
  c3merged[1][0][KT3Bin][CBOI_SS]->GetXaxis()->SetTitleOffset(0.9); c3merged[1][0][KT3Bin][CBOI_SS]->GetYaxis()->SetTitleOffset(0.8);
  c3merged[1][0][KT3Bin][CBOI_SS]->GetXaxis()->SetTitleFont(TextFont); c3merged[1][0][KT3Bin][CBOI_SS]->GetYaxis()->SetTitleFont(TextFont);
  c3merged[1][0][KT3Bin][CBOI_SS]->GetXaxis()->SetTitleSize(SizeTitle*SF1); c3merged[1][0][KT3Bin][CBOI_SS]->GetYaxis()->SetTitleSize(SizeTitle*SF1);
  c3merged[1][0][KT3Bin][CBOI_SS]->GetXaxis()->SetTitle("#font[12]{Q}_{3} (GeV/#font[12]{c})"); 
  c3merged[1][0][KT3Bin][CBOI_SS]->GetYaxis()->SetTitle("#font[12]{#bf{c}}_{3}^{#pm#pm#mp}");
  c3merged[1][0][KT3Bin][CBOI_SS]->GetYaxis()->SetNdivisions(504);
  c3merged[1][0][KT3Bin][CBOI_SS]->GetXaxis()->SetNdivisions(503);
  //
  c3merged[1][0][KT3Bin][CBOI_SS]->GetXaxis()->SetRangeUser(0,0.08);
  c3merged[1][0][KT3Bin][CBOI_SS]->SetMaximum(1.08);
  TLegend *legend8 = new TLegend(.8,.8, .98,.95,NULL,"brNDC");//c3
  legend8->SetBorderSize(0);
  legend8->SetTextFont(TextFont);
  legend8->SetTextSize(SizeLegend*SF1);
  legend8->SetFillColor(0);
  c3merged[1][0][KT3Bin][CBOI_SS]->DrawCopy();
  c3merged[1][1][KT3Bin][CBOI_SS]->DrawCopy("same");
  legend8->AddEntry(c3merged[1][0][KT3Bin][CBOI_SS],"GRS","p");
  legend8->AddEntry(c3merged[1][1][KT3Bin][CBOI_SS],"#Omega_{0}","p");
  legend8->Draw("same");
  BaseLine->Draw("same");
  
  
  //can6->SaveAs("SaveFigs/Fig_r3c3_CoulVar.eps");
  ///////////////////////////////////////////////////
  // Lambda variation
  
  TCanvas *can7 = new TCanvas("can7", "can7",1800,300,600,600);// 680,300,600,600
  can7->SetHighLightColor(2);
  gStyle->SetOptFit(0111);
  can7->SetFillColor(0);//10
  can7->SetBorderMode(0);
  can7->SetBorderSize(2);
  can7->SetFrameFillColor(0);
  can7->SetFrameBorderMode(0);
  can7->SetFrameBorderMode(0);
  
  TPad *pad7 = new TPad("pad7","pad7",0.0,0.0,1.,1.);
  gPad->SetGridx(0);
  gPad->SetGridy(0);
  //gPad->SetTickx();
  //gPad->SetTicky();
  pad7->SetBottomMargin(0.07);//0.12
  pad7->SetTopMargin(0.02);//0.05
  pad7->SetRightMargin(0.01);//1e-2
  pad7->Draw();
  TLegend *legend9 = new TLegend(.75,.7, .95,.85,NULL,"brNDC");//c3
  legend9->SetBorderSize(0);
  legend9->SetTextFont(TextFont);
  legend9->SetTextSize(SizeLegend*SF2);
  legend9->SetFillColor(0);
  
  gPad->SetRightMargin(0.01); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.01);
  r3merged[0][KT3Bin][CBOI_SS]->SetMinimum(r3MIN); r3merged[0][KT3Bin][CBOI_SS]->SetMaximum(r3MAX);
  r3merged[0][KT3Bin][CBOI_SS]->GetXaxis()->SetLabelFont(TextFont); r3merged[0][KT3Bin][CBOI_SS]->GetYaxis()->SetLabelFont(TextFont);
  r3merged[0][KT3Bin][CBOI_SS]->GetXaxis()->SetLabelSize(SizeLabel*SF2); r3merged[0][KT3Bin][CBOI_SS]->GetYaxis()->SetLabelSize(SizeLabel*SF2);
  //
  r3merged[0][KT3Bin][CBOI_SS]->GetXaxis()->SetTitleFont(TextFont); r3merged[0][KT3Bin][CBOI_SS]->GetYaxis()->SetTitleFont(TextFont); 
  r3merged[0][KT3Bin][CBOI_SS]->GetXaxis()->SetTitleSize(SizeTitle*SF2); r3merged[0][KT3Bin][CBOI_SS]->GetYaxis()->SetTitleSize(SizeTitle*SF2);
  r3merged[0][KT3Bin][CBOI_SS]->GetXaxis()->SetTitle("#font[12]{Q}_{3} (GeV/#font[12]{c})"); r3merged[0][KT3Bin][CBOI_SS]->GetYaxis()->SetTitle("#font[12]{r}_{3}");
  r3merged[0][KT3Bin][CBOI_SS]->GetXaxis()->SetTitleOffset(0.9); r3merged[0][KT3Bin][CBOI_SS]->GetYaxis()->SetTitleOffset(1.05); 
  TH1D *r3LambdaVaried = (TH1D*)r3merged[0][KT3Bin][CBOI_SS]->Clone();
  r3LambdaVaried->SetMarkerStyle(25);
  TF1 *r3LambdaSysFit = new TF1("r3LambdaSysFit","pol1",0,0.1);
  r3LambdaSysFit->FixParameter(0,LambdaSysPar0[CBOI_SS]);
  r3LambdaSysFit->FixParameter(1,LambdaSysPar1[CBOI_SS]);
  for(int bin=1; bin<=11; bin++){
    double Q3 = (bin-0.5)*0.01;
    r3LambdaVaried->SetBinContent(bin, (1-r3LambdaSysFit->Eval(Q3)) * r3LambdaVaried->GetBinContent(bin));
  }
  r3merged[0][KT3Bin][CBOI_SS]->GetXaxis()->SetNdivisions(503);
  r3merged[0][KT3Bin][CBOI_SS]->GetXaxis()->SetRangeUser(0,0.08);

  r3merged[0][KT3Bin][CBOI_SS]->DrawCopy();
  r3LambdaVaried->DrawCopy("same");
  ChaoticLimit->Draw("same");
  legend9->AddEntry(r3merged[0][KT3Bin][CBOI_SS],"#lambda=0.7","p");
  legend9->AddEntry(r3LambdaVaried,"#lambda=0.6","p");
  legend9->Draw("same");
  
  TLatex *Specif_6 = new TLatex(0.36,0.94,"ALICE Pb-Pb #\sqrt{#font[12]{s}_{NN}}=2.76 TeV");
  Specif_6->SetNDC();
  Specif_6->SetTextFont(TextFont);
  Specif_6->SetTextSize(SizeSpecif*SF2);
  TLatex *Specif2_6 = new TLatex(0.5,0.2,"5-10%");
  Specif2_6->SetNDC();
  Specif2_6->SetTextFont(TextFont);
  Specif2_6->SetTextSize(SizeLegend*SF2);
  Specif_6->Draw();
  Specif2_6->Draw();
  TString *KTString = new TString("");
  if(KT3Bin==0) KTString->Append("0.16<#font[12]{K}_{T,3}<0.3 GeV/#font[12]{c}");
  else KTString->Append("0.3<#font[12]{K}_{T,3}<1.0 GeV/#font[12]{c}");
  TLatex *Kt3Name_4 = new TLatex(0.36,0.28,KTString->Data());
  Kt3Name_4->SetNDC();
  Kt3Name_4->SetTextFont(TextFont);
  Kt3Name_4->SetTextSize(SizeLabel*SF2);
  Kt3Name_4->Draw("same");

  //can7->SaveAs("SaveFigs/Fig_r3_LamVar.eps");

  //TString *OutName=new TString("r3_M");
  //TString *OutName=new TString("c3mc_M");
  //TString *OutName=new TString("c3mc_GaussK3_M");
  //*OutName += BOI_1;
  //OutName->Append(".png");
  //can6->SaveAs(OutName->Data());
  

  
 
  /*
  // digitized from Wilson's Thesis
  double Wilson_Q3[18]={0};
  double Wilson_Q3_e[18]={0};
  double Wilson_C3[18]={2.02, 1.345, 1.51, 1.514, 1.51, 1.453, 1.394, 1.329, 1.277, 1.228, 1.188, 1.16, 1.132, 1.116, 1.092, 1.076, 1.064, 1.052};
  double Wilson_C3_eh[18]={2.97, 1.445, 1.54, 1.53, 1.518, 0,0,0,0,0,0,0,0,0,0,0,0,0};
  double Wilson_C3_e[18]={0};
  for(int i=0; i<18; i++){
    Wilson_Q3[i] = (i+0.5)*0.005 + .010;
    Wilson_C3_e[i] = Wilson_C3_eh[i]-Wilson_C3[i];
    if(Wilson_C3_e[i] > 1) Wilson_C3_e[i]=0;
  }
  TGraphErrors *gr_Wilson = new TGraphErrors(18, Wilson_Q3, Wilson_C3, Wilson_Q3_e, Wilson_C3_e);
  gr_Wilson->SetMarkerStyle(20);
  gr_Wilson->SetMarkerColor(2);
  gr_Wilson->SetLineColor(2);
  //gr_Wilson->Draw("P");
  //legend1->AddEntry(gr_Wilson,"STAR 0-12%: 130 GeV","p");
  
  // digitized from Phys. Rev. Lett. 87 (2001) 82301
  double STAR_130GeV_C2[20]={1.21, 1.211, 1.182, 1.135, 1.095, 1.064, 1.0389, 1.023, 1.014, 1.0047, 1.0016, 1.001, 1,1,1,1,1,1,1,1};
  double STAR_130GeV_C2_e[20]={0};
  double STAR_130GeV_C2_qinv[20]={0};
  for(int i=0; i<20; i++) STAR_130GeV_C2_qinv[i] = (i+0.5)*0.005;
  TGraphErrors *gr_STAR_C2_130 = new TGraphErrors(20,STAR_130GeV_C2_qinv, STAR_130GeV_C2, 0,0);
  gr_STAR_C2_130->SetMarkerStyle(20);
  gr_STAR_C2_130->SetMarkerColor(2);
  //gr_STAR_C2_130->Draw("P");
  //legend1->AddEntry(gr_STAR_C2_130,"STAR 0-12%: 130 GeV","p");
  
  TF1 *line=new TF1("line","3",0,100);
  line->SetLineStyle(1);
  //line->Draw("same");
  //legend1->AddEntry(line,"Chaotic limit #font[12]{c}_{3}","l");
  TF1 *line2=new TF1("line2","6",0,100);
  line2->SetLineStyle(1);
  line2->SetLineColor(4);
  //line2->Draw("same");
  //legend1->AddEntry(line2,"Chaotic limit C_{3}","l");
  TF1 *line3=new TF1("line3","1",0,100);
  line3->SetLineStyle(1);
  //line3->Draw("same");
  //legend1->AddEntry(line3,"Chaotic limit","l");
  */
  /*
  //TLatex *tex = new TLatex(0.001,0.09,"q_{ij} > 0.01 GeV/c");//1.54,2300
  //TLatex *tex = new TLatex(0.05,4,"#\lambda=0.4#pm0.04, R=7#pm0.5 fm");// for C3
  //TLatex *tex = new TLatex(0.05,1.17,"#\lambda=0.4#pm0.04, R=7#pm0.5 fm");// for C2(+-) 0-5%
  //TLatex *tex = new TLatex(0.05,1.17,"#\lambda=0.4#pm0.04, R=4#pm0.5 fm");// for C2(+-) 45-50%
  TLatex *tex = new TLatex(0.0005,-.9,"#\lambda=0.4#pm0.04, R=7#pm0.5 fm");// for r3
  
  tex->SetTextSize(.04);
  //tex->Draw("same");

  TLatex *tex2 = new TLatex(0.002,5.75,"Chaotic Limit C_{3}");
  tex2->SetTextSize(.04);
  tex2->SetTextColor(4);
  //tex2->Draw("same");

  TLatex *tex3 = new TLatex(0.002,2.75,"Chaotic Limit #font[12]{c}_{3}");
  tex3->SetTextSize(.04);
  tex3->SetTextColor(2);
  //tex3->Draw("same");
  
  TLatex *tex4 = new TLatex(0.0005,-0.5,"Spherical Gaussian Coulomb WF");
  tex4->SetTextSize(.05);
  tex4->SetTextColor();
  //tex4->Draw("same");

  TLatex *tex5 = new TLatex(0.0005,.5,"Pb-Pb #sqrt{#font[12]{s}_{NN}}=2.76 TeV: 0-10% Centrality");
  tex5->SetTextSize(.04);
  tex5->SetTextColor();
  //tex5->Draw("same");
  */
  //
  /*
  TF1 *C2noConstraint=new TF1("C2noConstraint","2-pow(x,2)",0,1);
  C2noConstraint->SetLineColor(2);
  TF1 *C2Constraint=new TF1("C2Constraint","2-4/5.*pow(x,2)",0,1);
  C2Constraint->SetLineColor(4);
  C2Constraint->GetXaxis()->SetTitle("G");
  C2Constraint->GetYaxis()->SetTitle("C_{2}^{++} Intercept");
  C2Constraint->GetXaxis()->SetTitleSize(0.06);
  C2Constraint->GetYaxis()->SetTitleSize(0.06);
  C2Constraint->GetXaxis()->SetTitleOffset(.8);
  C2Constraint->GetYaxis()->SetTitleOffset(.7);
  C2Constraint->SetTitle("");
  TF1 *C3noConstraint=new TF1("C3noConstraint","6-9*pow(x,2)+4*pow(x,3)",0,1);
  C3noConstraint->SetLineColor(2);
  //C3noConstraint->SetLineStyle(2);
  TF1 *C3Constraint=new TF1("C3Constraint","6-36/5.*pow(x,2)+96/35.*pow(x,3)",0,1);
  C3Constraint->SetLineColor(4);
  //C3Constraint->SetLineStyle(2);
  C3Constraint->GetXaxis()->SetTitle("G");
  C3Constraint->GetYaxis()->SetTitle("C_{3}^{+++} Intercept");
  C3Constraint->GetXaxis()->SetTitleSize(0.06);
  C3Constraint->GetYaxis()->SetTitleSize(0.06);
  C3Constraint->GetXaxis()->SetTitleOffset(.8);
  C3Constraint->GetYaxis()->SetTitleOffset(.7);
  C3Constraint->SetTitle("");
  C3Constraint->Draw();
  legend1->AddEntry(C3Constraint,"C_{3}^{+++} Charge-Constrained","l");
  C3noConstraint->Draw("same");
  legend1->AddEntry(C3noConstraint,"C_{3}^{+++} Charge-Unconstrained","l");
  //C2Constraint->Draw();
  //legend1->AddEntry(C2Constraint,"C_{2}^{++} Charge-Constrained","l");
  //C2noConstraint->Draw("same");
  //legend1->AddEntry(C2noConstraint,"C_{2}^{++} Charge-Unconstrained","l");
  legend1->SetTextSize(.03);
  */
  
  /*TF1 *r3noConstraint = new TF1("r3noConstraint","(2-6*pow(x,2)+4*pow(x,3))/pow(1-pow(x,2),1.5)",0,1);
  r3noConstraint->SetLineColor(2);
  r3noConstraint->GetXaxis()->SetTitle("G");
  r3noConstraint->GetYaxis()->SetTitle("r_{3} Intercept");
  r3noConstraint->GetXaxis()->SetTitleSize(0.06);
  r3noConstraint->GetYaxis()->SetTitleSize(0.06);
  r3noConstraint->GetXaxis()->SetTitleOffset(.8);
  r3noConstraint->GetYaxis()->SetTitleOffset(.7);
  r3noConstraint->SetTitle("");
  TF1 *r3Constraint = new TF1("r3Constraint","(2-24/5.*pow(x,2)+96/35.*pow(x,3))/pow(1-4/5.*pow(x,2),1.5)",0,1);
  r3Constraint->SetLineColor(4);
  r3Constraint->GetXaxis()->SetTitle("G");
  r3Constraint->GetYaxis()->SetTitle("r_{3} Intercept");
  r3Constraint->GetXaxis()->SetTitleSize(0.06);
  r3Constraint->GetYaxis()->SetTitleSize(0.06);
  r3Constraint->GetXaxis()->SetTitleOffset(.8);
  r3Constraint->GetYaxis()->SetTitleOffset(.7);
  r3Constraint->SetTitle("");
  r3Constraint->Draw();
  //legend1->AddEntry(r3Constraint,"Charge-Constrained","l");
  r3noConstraint->Draw("same");
  //legend1->AddEntry(r3noConstraint,"Charge-Unconstrained","l");
  //cout<<setprecision(5)<<r3noConstraint->Eval(0.3)<<"  "<<r3Constraint->Eval(0.3)<<endl;
  */
  //legend1->Draw("same");

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
  TASImage *myAliceLogo = new TASImage((prel) ? "~/Pictures/alice_logo_preliminary.eps" : "~/Pictures/alice_logo_performance.eps");
  myAliceLogo->Draw();
}

