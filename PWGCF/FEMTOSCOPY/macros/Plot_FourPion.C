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
// 0(standard), 1(fcSq=0.65), 2(fcSq=0.75), 3(B minus), 4(B plus), 
// 5(Linear instead of Cubic,also fc^2=0.75), 6(10h), 7(MRC 10% increase), 8(MuonCorr 92% pion-pair purity)
int FileSetting=0;
//
bool MCcase_def=0;// MC data?
int CollisionType=1;// PbPb, pPb, pp
//
int Mbin_def=0;// 0-9: centrality bin in widths of 5%
bool SameCharge_def=1;// same-charge?
int CHARGE_def=-1;// -1 or +1: + or - pions for same-charge case, --+ or -++,  ---+ or -+++
int MixedCharge4pionType_def = 2;// 1(---+) or 2(--++)
//
int EDbin_def=0;// 0 or 1: Kt3 bin
int TPNbin=0;// TPN bin for r3 and r4
int Gbin = int( (0) /2. ) + 55;// +5 (Rcoh=0), +25 (Rcoh=Rch) or +55 for extended G range 
int c3FitType = 2;// EW(1), LG(2)
int Ktbin_def=1;// 1(0.2-0.3),..., 6(0.7-0.8), 10 = Full Range
//
bool MRC=1;// Momentum Resolution Corrections?
bool MuonCorrection=1;// correct for Muons?
bool FSICorrection=1;// correct For Final-State-Interactions
bool InterpCorrection=1;// correct for finite bin interpolation
//
int f_choice=0;// 0(Core/Halo), 1(40fm), 2(70fm), 3(100fm)
//
//
bool SaveToFile_def=kFALSE;// Save outputs to file?
bool GeneratedSignal=kFALSE;// Was the QS+FSI signal generated in MC? 
//
//
//
//

int fFSIindex=0;
int Ktlowbin;
int Kthighbin;
float TwoFrac;// Lambda parameter
float OneFrac;// Lambda^{1/2}
float ThreeFrac;// Lambda^{3/2}
float FourFrac;// lambda^{4/2}
double Qcut_pp = 0.6;// 0.6
double Qcut_PbPb = 0.1;
double NormQcutLow_pp = 1.0;
double NormQcutHigh_pp = 1.5;
double NormQcutLow_PbPb = .15;
double NormQcutHigh_PbPb = .2;// was .175

int ThreeParticleRebin;
int FourParticleRebin;

int RbinMRC;

TH1D *fFSIss[12];
TH1D *fFSIos[12];

void SetFSICorrelations();
void SetFSIindex(Float_t);
Float_t FSICorrelation(Int_t, Int_t, Float_t);
void SetMuonCorrections();
void SetMomResCorrections();
double Gamov(int, double);
double cubicInterpolate(double[4], double);
double bicubicInterpolate(double[4][4], double, double);
double tricubicInterpolate(double[4][4][4], double, double, double);
//
float fact(float);


TH1D *MRC_SC_2[2];
TH1D *MRC_MC_2[2];
TH1D *MRC_SC_3[3];
TH1D *MRC_MC_3[4];
TH1D *MRC_SC_4[5];
TH1D *MRC_MC1_4[6];
TH1D *MRC_MC2_4[6];


double AvgQinvSS[30];
double AvgQinvOS[30];
double BinCentersQ4[20];

//
TH1D *C2muonCorrectionSC;
TH1D *C2muonCorrectionMC;
TH1D *WeightmuonCorrection;
TH1D *C3muonCorrectionSC[2];
TH1D *C3muonCorrectionMC[3];
TH1D *C4muonCorrectionSC[4];
TH1D *C4muonCorrectionMC1[5];
TH1D *C4muonCorrectionMC2[5];


void Plot_FourPion(bool SaveToFile=SaveToFile_def, bool MCcase=MCcase_def, bool SameCharge=SameCharge_def, int MixedCharge4pionType=MixedCharge4pionType_def, int EDbin=EDbin_def, int CHARGE=CHARGE_def, int Mbin=Mbin_def, int Ktbin=Ktbin_def){
  
  EDbin_def=EDbin;
  SaveToFile_def=SaveToFile;
  MCcase_def=MCcase;
  CHARGE_def=CHARGE;
  SameCharge_def=SameCharge;// 3-pion same-charge
  MixedCharge4pionType_def=MixedCharge4pionType;
  Mbin_def=Mbin;
  Ktbin_def=Ktbin;
  //
  Ktlowbin=(Ktbin)*2+3;// kt bins are 0.5 GeV/c wide (0-0.5, 0.5-1.0 ...)
  Kthighbin=(Ktbin)*2+4;
  
  //
  if(FileSetting==1) TwoFrac=0.65;
  else if(FileSetting==2 || FileSetting==5) TwoFrac=0.75;
  else TwoFrac=0.7;

  OneFrac = sqrt(TwoFrac);
  ThreeFrac = pow(TwoFrac, 1.5);
  FourFrac = pow(TwoFrac, 2.);

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

  //// Core/Halo, 40fm, 70fm, 100fm
  float ThermShift_f44[4]={0., 0.08741, 0.005284, -0.01064};
  float ThermShift_f43[4]={0., -0.01126, -0.001486, 0.001686};
  float ThermShift_f42[4]={0., -0.006466, -7.683e-05, 0.0004572};
  float ThermShift_f41[4]={0., -0.003556, 0.00112, 0.00115};
  float f44prime = FourFrac;
  float f43prime = ThreeFrac*(1-OneFrac);
  float f42prime = TwoFrac*pow(1-OneFrac,2);
  float f41prime = pow(1-OneFrac,4) + 4*OneFrac*pow(1-OneFrac,3);
  f44prime += ThermShift_f44[f_choice];
  f43prime += ThermShift_f43[f_choice];
  f42prime += ThermShift_f42[f_choice];
  f41prime += ThermShift_f41[f_choice];
  float f44 = f44prime;
  float f43 = f43prime/f33prime;
  float f42 = f42prime/TwoFrac;
  f42 -= 2*f43prime*f32prime/f33prime/TwoFrac;
  float f41 = f41prime;
  f41 -= 4*f43prime*f31prime/f33prime;
  f41 -= 6*f42prime*(1-TwoFrac)/TwoFrac;
  f41 += 12*f43prime/f33prime*f32prime/TwoFrac*(1-TwoFrac);
  //cout<<f44 + 4*f43 + 6*f42 + f41<<endl;
  //cout<<f44<<"  "<<f43<<"  "<<f42<<"  "<<f41<<endl;
  //
  // Core/Halo, 40fm, 70fm, 100fm
  //float TherminatorMod_f33[4]={1., 1.123, 1.022, 1.008};// 1.008
  //float TherminatorMod_f32[4]={1., 0.844, 0.933, 0.967};// .967
  //float TherminatorMod_f31[4]={1., 0.828, 0.982, 1.028};// 1.028
  /*float TherminatorMod_f44[4]={1., 1.188, 1.008, 0.985};// .985
  float TherminatorMod_f43[4]={1., 0.885, 0.979, 1.026};// 1.026
  float TherminatorMod_f42[4]={1., 0.687, 0.99, 1.08};// 1.08
  float TherminatorMod_f41[4]={1., 0.806, 1.33, 1.548};//1.548
  float f44 = FourFrac;
  float f43 = (1-OneFrac);
  float f42 = -pow(1-OneFrac,2);
  float f41 = -3*pow(1-OneFrac,4) - 8*OneFrac*pow(1-OneFrac,3) + 6*(1-TwoFrac)*pow(1-OneFrac,2);
  f44 *= TherminatorMod_f44[f_choice];
  f43 *= TherminatorMod_f43[f_choice];
  f42 *= TherminatorMod_f42[f_choice];
  f41 *= TherminatorMod_f41[f_choice];*/
  //f44prime *= TherminatorMod_f44[f_choice];
  //f43prime *= TherminatorMod_f43[f_choice];
  //f42prime *= TherminatorMod_f42[f_choice];
  //f41prime *= TherminatorMod_f41[f_choice];
  //float f44 = f44prime;
  //float f43 = f43prime/ThreeFrac;
  //float f42 = f42prime/TwoFrac - 2*pow(1-OneFrac,2);
  //float f41 = f41prime - 4*pow(1-OneFrac,4) - 12*OneFrac*pow(1-OneFrac,3) + 6*(1-TwoFrac)*pow(1-OneFrac,2);

  cout<<"Mbin = "<<Mbin<<"   KT3 = "<<EDbin<<"   SameCharge = "<<SameCharge<<endl;
  if(!SameCharge) cout<<"4-pion MixedCharge type = "<<MixedCharge4pionType<<endl;
  
  if(CollisionType==0){
    if(Mbin==0) {RbinMRC=10;}
    else if(Mbin==1) {RbinMRC=9;}
    else if(Mbin<=3) {RbinMRC=8;}
    else if(Mbin<=5) {RbinMRC=7;}
    else {RbinMRC=6;}
  }else{
    RbinMRC=2;
  }
  
  if(MCcase) MuonCorrection=kFALSE;

  if(CollisionType==0){
    ThreeParticleRebin=2;
    FourParticleRebin=3;
  }else{
    ThreeParticleRebin=2;
    FourParticleRebin=6;
  }
  float Cutoff_FullWeight_Q3[10]={0.065, 0.065, 0.075, 0.075, 0.075, 0.075, 0.085, 0.085, 0.085, 0.085};
  float Cutoff_FullWeight_Q4[10]={0.115, 0.115, 0.115, 0.130, 0.130, 0.130, 0.145, 0.145, 0.145, 0.145};

 
  //
  // avg Qinv within the 5 MeV bins (0-5, 5-10,...) for true bin mean values.  From unsmeared HIJING 0-5% with input signal
  double AvgQinvSS_temp[30]={0.00421164, 0.00796583, 0.0128473, 0.0178262, 0.0228416, 0.0276507, 0.0325368, 0.0376114, 0.0424707, 0.0476274, 0.0526015, 0.0575645, 0.0625478, 0.0675416, 0.0725359, 0.0775219, 0.0825521, 0.0875211, 0.0925303, 0.0975319, 0.102544, 0.10753, 0.112499, 0.117545, 0.122522, 0.127522, 0.132499, 0.137514, 0.142495, 0.147521};
  double AvgQinvOS_temp[30]={0.00352789, 0.00797596, 0.0128895, 0.0177198, 0.0226397, 0.0276331, 0.0326309, 0.0376471, 0.0426217, 0.047567, 0.0525659, 0.0575472, 0.0625886, 0.0675261, 0.0725543, 0.077564, 0.0825617, 0.0875465, 0.092539, 0.0975348, 0.102529, 0.107527, 0.112506, 0.117531, 0.122536, 0.1275, 0.132508, 0.137508, 0.14251, 0.147511};
  for(int i=0; i<30; i++){
    AvgQinvSS[i] = AvgQinvSS_temp[i];
    AvgQinvOS[i] = AvgQinvOS_temp[i];
  }
 

  
  TH2D *TwoParticle_2d[2][2][2];// ch1,ch2,term
  TH1D *TwoParticle[2][2][2];// ch1,ch2,term
  double norm_2[2]={0};
  //
  TH1D *ThreeParticle[2][2][2][5];// ch1,ch2,ch3,term
  TProfile *K3avg[2][2][2][4];
  TH2D *TPFullWeight_ThreeParticle_2D[2];// charge
  TH2D *TPNegFullWeight_ThreeParticle_2D[2];// charge
  TH1D *TPN_ThreeParticle[2];// charge
  TH1D *TPFullWeight_ThreeParticle[2];// charge
  TH1D *TPNegFullWeight_ThreeParticle[2];// charge
  double norm_3[5]={0};
  //
  TH1D *FourParticle[2][2][2][2][13];// ch1,ch2,ch3,ch4,term
  TProfile *K4avg[2][2][2][2][12];
  TH2D *TPFullWeight_FourParticle_2D[2];// charge
  TH2D *TPNegFullWeight_FourParticle_2D[2];// charge
  TH1D *TPFullWeight_FourParticle[2];// charge
  TH1D *TPNegFullWeight_FourParticle[2];// charge
  TH1D *TPN_FourParticle[2];// charge
  TH3D *FullBuildFromFits_3D;// charge
  TH1D *FullBuildFromFits;// charge
  TH3D *PartialBuildFromFits_3D;// charge
  TH1D *PartialBuildFromFits;// charge
  double norm_4[13]={0};
  

  gStyle->SetOptStat(0);
  gStyle->SetOptDate(0);
  gStyle->SetOptFit(0111);

  TFile *_file0;
  if(CollisionType==0){// PbPb
    if(MCcase) {
      if(Mbin<=1){
	//_file0 = new TFile("Results/PDC_12a17a_GeneratedSignal.root","READ");
	//_file0 = new TFile("Results/PDC_12a17a_pTSpectrumWeight.root","READ");
	_file0 = new TFile("Results/PDC_12a17a_Qweights.root","READ");
      }
    }else{
      //if(FileSetting==0) _file0 = new TFile("Results/PDC_11h_pT_0p2to1p0_FullRunWrongWeightsNoPadRowTTCandMCTTC_RawWeightFile.root","READ");
      //if(FileSetting==0) _file0 = new TFile("Results/PDC_11h_pT_0p16to0p25_0p25to1p0.root","READ");
      //if(FileSetting==0) _file0 = new TFile("Results/PDC_11h_pT_0p2to1p0.root","READ");
      //if(FileSetting==0) _file0 = new TFile("Results/PDC_11h_pT0p16to0p25_KT0p35.root","READ");
      //if(FileSetting==0) _file0 = new TFile("Results/PDC_11h_1percentCentEM.root","READ");
      //if(FileSetting==0) _file0 = new TFile("Results/PDC_11h_0p02eta0p045phi_0p03eta0p067phi.root","READ");
      if(FileSetting==0) _file0 = new TFile("Results/PDC_11h_c3FitBuild.root","READ");
      //if(FileSetting==0) _file0 = new TFile("Results/PDC_11h_extendedGweights.root","READ");// Preliminary results
      if(FileSetting==5) _file0 = new TFile("Results/PDC_11h_Cubic_Linear.root","READ");
      if(FileSetting==1 || FileSetting==2) _file0 = new TFile("Results/PDC_11h_Lam0p65_Lam0p75.root","READ");
      if(FileSetting==3) _file0 = new TFile("Results/PDC_11h_Cubic_Linear_Bminus.root","READ");
      if(FileSetting==4) _file0 = new TFile("Results/PDC_11h_Cubic_Linear_Bplus.root","READ");
      if(FileSetting==6) _file0 = new TFile("Results/PDC_10h_Cubic_Linear.root","READ");
      if(FileSetting==7 || FileSetting==8) _file0 = new TFile("Results/PDC_11h_MRC10percIncrease_Muon92percent.root","READ");
    }
  }else if(CollisionType==1){// pPb
    _file0 = new TFile("Results/PDC_13bc_kINT7.root","READ");
  }else{// pp
    _file0 = new TFile("Results/PDC_10bcde_kMB.root","READ");
  }
  
  SetFSIindex(10.);
  SetFSICorrelations();
  SetMomResCorrections();
  SetMuonCorrections();
  //
  /////////////////////////////////////////////////////////
  

  double NormQcutLow;
  double NormQcutHigh;
  if(CollisionType==0) {
    NormQcutLow = NormQcutLow_PbPb;
    NormQcutHigh = NormQcutHigh_PbPb;
  }else {
    NormQcutLow = NormQcutLow_pp;
    NormQcutHigh = NormQcutHigh_pp;
  }

 
  TList *MyList;
  if(!MCcase){
    TDirectoryFile *tdir = (TDirectoryFile*)_file0->Get("PWGCF.outputFourPionAnalysis.root");
    if(FileSetting==2 || FileSetting==5 || FileSetting==8) MyList=(TList*)tdir->Get("FourPionOutput_2");
    else MyList=(TList*)tdir->Get("FourPionOutput_1");
    //MyList=(TList*)_file0->Get("MyList");
  }else{
    MyList=(TList*)_file0->Get("MyList");
  }
  _file0->Close();

 
  TH1D *Events = (TH1D*)MyList->FindObject("fEvents2");
  //

  cout<<"#Events = "<<Events->Integral(Mbin+1,Mbin+1)<<endl;

  
  
  ///////////////////////////////////
  // Get Histograms
  
  // 2-particle
  for(int c1=0; c1<2; c1++){
    for(int c2=0; c2<2; c2++){
      for(int term=0; term<2; term++){
	
	TString *name2 = new TString("TwoParticle_Charge1_");
	*name2 += c1;
	name2->Append("_Charge2_");
	*name2 += c2;
	name2->Append("_M_");
	*name2 += Mbin;
	name2->Append("_ED_");
	*name2 += 0;// EDbin
	name2->Append("_Term_");
	*name2 += term+1;
	
	if( (c1+c2)==1 ) {if(c1!=0) continue;}// skip degenerate histogram
	
	TwoParticle_2d[c1][c2][term] = (TH2D*)MyList->FindObject(name2->Data());
	TwoParticle_2d[c1][c2][term]->Sumw2();
	TString *proname = new TString(name2->Data());
	proname->Append("_pro");
	
	if(Ktbin==10) {Ktlowbin=1; Kthighbin=TwoParticle_2d[c1][c2][term]->GetNbinsX();}
	TwoParticle[c1][c2][term] = (TH1D*)TwoParticle_2d[c1][c2][term]->ProjectionY(proname->Data(),Ktlowbin,Kthighbin);// bins:5-6 (kt:0.2-0.3)
	
	norm_2[term] = TwoParticle[c1][c2][term]->Integral(TwoParticle[c1][c2][term]->GetXaxis()->FindBin(NormQcutLow),TwoParticle[c1][c2][term]->GetXaxis()->FindBin(NormQcutHigh));
	//cout<<"2-pion norms  "<<norm_2[term]<<endl;
	TwoParticle[c1][c2][term]->Scale(norm_2[0]/norm_2[term]);
	
	TwoParticle[c1][c2][term]->SetMarkerStyle(20);
	TwoParticle[c1][c2][term]->GetXaxis()->SetTitle("q_{inv} (GeV/c)");
	TwoParticle[c1][c2][term]->GetYaxis()->SetTitle("C_{2}");
	TwoParticle[c1][c2][term]->SetTitle("");
	
      }// term
      
    
      // 3-particle
      for(int c3=0; c3<2; c3++){
	for(int term=0; term<5; term++){
	   
	  TString *name3 = new TString("ThreeParticle_Charge1_");
	  *name3 += c1;
	  name3->Append("_Charge2_");
	  *name3 += c2;
	  name3->Append("_Charge3_");
	  *name3 += c3;
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
	  nameK3->Append("_Kfactor");
	  //
	  TString *nameTPN3=new TString(name3->Data());
	  nameTPN3->Append("_TwoPartNorm");
	  //
	  TString *nameNegTPN3=new TString(name3->Data());
	  nameNegTPN3->Append("_TwoPartNegNorm");
	  //
	  name3->Append("_1D");
	  
	  ///////////////////////////////////////
	  // skip degenerate histograms
	  if( (c1+c2+c3)==1) {if(c1!=0 || c2!=0 || c3!=1) continue;}
	  if( (c1+c2+c3)==2) {if(c1!=0) continue;}
	  ////////////////////////////////////////
	  
	  norm_3[term] = ((TH1D*)MyList->FindObject(nameNorm3->Data()))->Integral();
	  ThreeParticle[c1][c2][c3][term] = (TH1D*)MyList->FindObject(name3->Data());
	  ThreeParticle[c1][c2][c3][term]->Sumw2();
	  //if(c1==c2 && c1==c3) cout<<"3-pion norms  "<<norm_3[term]<<endl;
	  ThreeParticle[c1][c2][c3][term]->Scale(norm_3[0]/norm_3[term]);
	  ThreeParticle[c1][c2][c3][term]->GetXaxis()->SetTitle("Q_{3} (GeV/c)");
	  ThreeParticle[c1][c2][c3][term]->GetYaxis()->SetTitle("C_{3}");
	  ThreeParticle[c1][c2][c3][term]->SetMarkerStyle(20);
	  ThreeParticle[c1][c2][c3][term]->SetTitle("");
	  //
	  ThreeParticle[c1][c2][c3][term]->Rebin(ThreeParticleRebin);
	  
	  //
	  if(term<4){
	    K3avg[c1][c2][c3][term] = (TProfile*)MyList->FindObject(nameK3->Data());
	    K3avg[c1][c2][c3][term]->Rebin(ThreeParticleRebin);
	    if(MCcase || !FSICorrection) {
	      K3avg[c1][c2][c3][term]->Reset();
	      for(int ii=1; ii<=K3avg[c1][c2][c3][term]->GetNbinsX(); ii++) {
		K3avg[c1][c2][c3][term]->Fill(K3avg[c1][c2][c3][term]->GetXaxis()->GetBinCenter(ii), 1);
	      }	
	    }
	  }
	  //
	  if(term==4 && c1==c2 && c1==c3){
	    TPFullWeight_ThreeParticle_2D[c1] = (TH2D*)MyList->FindObject(nameTPN3->Data());
	    TPFullWeight_ThreeParticle_2D[c1]->Scale(norm_3[0]/norm_3[term]);
	    TPFullWeight_ThreeParticle_2D[c1]->RebinY(ThreeParticleRebin);
	    //
	    TPNegFullWeight_ThreeParticle_2D[c1] = (TH2D*)MyList->FindObject(nameNegTPN3->Data());
	    TPNegFullWeight_ThreeParticle_2D[c1]->Scale(norm_3[0]/norm_3[term]);
	    TPNegFullWeight_ThreeParticle_2D[c1]->RebinY(ThreeParticleRebin);
	    //
	    for(int binx=1; binx<=TPFullWeight_ThreeParticle_2D[c1]->GetNbinsX(); binx++) {
	      for(int biny=1; biny<=TPFullWeight_ThreeParticle_2D[c1]->GetNbinsY(); biny++) {
		TPFullWeight_ThreeParticle_2D[c1]->SetBinError(binx, biny, 0);
		if(binx!=4){
		  TPFullWeight_ThreeParticle_2D[c1]->SetBinContent(binx, biny, TPFullWeight_ThreeParticle_2D[c1]->GetBinContent(binx, biny) + TPFullWeight_ThreeParticle_2D[c1]->GetBinContent(4, biny));
		  TPNegFullWeight_ThreeParticle_2D[c1]->SetBinContent(binx, biny, TPNegFullWeight_ThreeParticle_2D[c1]->GetBinContent(binx, biny) + TPNegFullWeight_ThreeParticle_2D[c1]->GetBinContent(4, biny));
		  if(InterpCorrection){
		    double q3 = TPFullWeight_ThreeParticle_2D[c1]->GetYaxis()->GetBinCenter(biny);
		    double InterCorr = 1.024 - 0.2765*q3;
		    //if(q3<0.1) cout<<binx<<"  "<<biny<<"  "<<q3<<" "<<InterCorr<<endl;
		    TPFullWeight_ThreeParticle_2D[c1]->SetBinContent(binx, biny, InterCorr * TPFullWeight_ThreeParticle_2D[c1]->GetBinContent(binx, biny));
		    TPNegFullWeight_ThreeParticle_2D[c1]->SetBinContent(binx, biny, InterCorr * TPNegFullWeight_ThreeParticle_2D[c1]->GetBinContent(binx, biny));
		  }
		}
	      }
	    }
	    
	    TString *proName=new TString(nameTPN3->Data()); TString *proNameNeg=new TString(nameNegTPN3->Data());
	    proName->Append("_pro");  proNameNeg->Append("_pro");
	    TPN_ThreeParticle[c1] = (TH1D*)TPFullWeight_ThreeParticle_2D[c1]->ProjectionY(proName->Data(), TPNbin, TPNbin);
	    //
	    proName->Append("_FullWeight"); proNameNeg->Append("_FullWeight");
	    TPFullWeight_ThreeParticle[c1] = (TH1D*)TPFullWeight_ThreeParticle_2D[c1]->ProjectionY(proName->Data(), Gbin, Gbin);
	    TPNegFullWeight_ThreeParticle[c1] = (TH1D*)TPNegFullWeight_ThreeParticle_2D[c1]->ProjectionY(proNameNeg->Data(), Gbin, Gbin);
	    proName->Append("_FullWeightDen"); proNameNeg->Append("_FullWeightDen");
	    TH1D *tempDen = (TH1D*)TPFullWeight_ThreeParticle_2D[c1]->ProjectionY(proName->Data(), 4, 4);
	    TH1D *tempDenNeg = (TH1D*)TPNegFullWeight_ThreeParticle_2D[c1]->ProjectionY(proNameNeg->Data(), 4, 4);
	    // Add Pos with Neg weights
	    tempDen->Add(tempDenNeg);
	    TPFullWeight_ThreeParticle[c1]->Add(TPNegFullWeight_ThreeParticle[c1]);
	    //
	    //TPFullWeight_ThreeParticle[c1]->Add(tempDen);// now added above in Interp section
	    TPFullWeight_ThreeParticle[c1]->Divide(tempDen);
	    TPFullWeight_ThreeParticle[c1]->SetLineColor(1);
	    //
	  }

	}// term 3-particle

	
	// 4-particle
	for(int c4=0; c4<2; c4++){
	  for(int term=0; term<13; term++){
	    
	    TString *name4 = new TString("FourParticle_Charge1_");
	    *name4 += c1;
	    name4->Append("_Charge2_");
	    *name4 += c2;
	    name4->Append("_Charge3_");
	    *name4 += c3;
	    name4->Append("_Charge4_");
	    *name4 += c4;
	    name4->Append("_M_");
	    *name4 += Mbin;
	    name4->Append("_ED_");
	    *name4 += EDbin;
	    name4->Append("_Term_");
	    *name4 += term+1;
	    
	    TString *nameNorm4=new TString(name4->Data());
	    nameNorm4->Append("_Norm");
	    //
	    TString *nameK4=new TString(name4->Data());
	    nameK4->Append("_Kfactor");// or Weighted
	    //
	    TString *nameTPN4=new TString(name4->Data());
	    nameTPN4->Append("_TwoPartNorm");
	    //
	    TString *nameNegTPN4=new TString(name4->Data());
	    nameNegTPN4->Append("_TwoPartNegNorm");
	    //
	    TString *nameFitBuild=new TString(name4->Data());
	    nameFitBuild->Append("_FullBuildFromFits");
	    //
	    TString *namePartialFitBuild=new TString(name4->Data());
	    namePartialFitBuild->Append("_PartialBuildFromFits");
	    //
	    name4->Append("_1D");
	    ///////////////////////////////////////
	    // skip degenerate histograms
	    if( (c1+c2+c3+c4)==1) {if(c4!=1) continue;}
	    if( (c1+c2+c3+c4)==2) {if(c3+c4!=2) continue;}
	    if( (c1+c2+c3+c4)==3) {if(c1!=0) continue;}
	    /////////////////////////////////////////
	    norm_4[term] = ((TH1D*)MyList->FindObject(nameNorm4->Data()))->Integral();
	    //if( (c1+c2+c3+c4)==4) cout<<"4-pion norms  "<<norm_4[term]<<endl;
	    if(norm_4[term] > 0){
	      FourParticle[c1][c2][c3][c4][term] = (TH1D*)MyList->FindObject(name4->Data());
	      FourParticle[c1][c2][c3][c4][term]->Sumw2();
	      FourParticle[c1][c2][c3][c4][term]->Scale(norm_4[0]/norm_4[term]);
	      FourParticle[c1][c2][c3][c4][term]->GetXaxis()->SetTitle("Q_{4} (GeV/c)");
	      FourParticle[c1][c2][c3][c4][term]->GetYaxis()->SetTitle("C_{4}");
	      FourParticle[c1][c2][c3][c4][term]->SetMarkerStyle(20);
	      FourParticle[c1][c2][c3][c4][term]->SetTitle("");
	      //
	      FourParticle[c1][c2][c3][c4][term]->Rebin(FourParticleRebin);
	     
	    }else{
	      cout<<"4-particle normalization = 0."<<endl;
	    }	      
	    if(term<12){
	      K4avg[c1][c2][c3][c4][term] = (TProfile*)MyList->FindObject(nameK4->Data());
	      K4avg[c1][c2][c3][c4][term]->Rebin(FourParticleRebin);
	      if(MCcase || !FSICorrection) {
		K4avg[c1][c2][c3][c4][term]->Reset();
		for(int ii=1; ii<=K4avg[c1][c2][c3][c4][term]->GetNbinsX(); ii++) {
		  K4avg[c1][c2][c3][c4][term]->Fill(K4avg[c1][c2][c3][c4][term]->GetXaxis()->GetBinCenter(ii), 1);
		}	
	      }
	    }
	    if(term==12 && c1==c2 && c1==c3 && c1==c4){
	     
	      TPFullWeight_FourParticle_2D[c1] = (TH2D*)MyList->FindObject(nameTPN4->Data());
	      TPFullWeight_FourParticle_2D[c1]->Scale(norm_4[0]/norm_4[term]);
	      TPFullWeight_FourParticle_2D[c1]->RebinY(FourParticleRebin);
	      //
	      TPNegFullWeight_FourParticle_2D[c1] = (TH2D*)MyList->FindObject(nameNegTPN4->Data());
	      TPNegFullWeight_FourParticle_2D[c1]->Scale(norm_4[0]/norm_4[term]);
	      TPNegFullWeight_FourParticle_2D[c1]->RebinY(FourParticleRebin);
	      //
	      if(c1==0){
		FullBuildFromFits_3D = (TH3D*)MyList->FindObject(nameFitBuild->Data());
		FullBuildFromFits_3D->Scale(norm_4[0]/norm_4[term]);
		FullBuildFromFits_3D->RebinZ(FourParticleRebin);
		//
		PartialBuildFromFits_3D = (TH3D*)MyList->FindObject(namePartialFitBuild->Data());
		PartialBuildFromFits_3D->Scale(norm_4[0]/norm_4[term]);
		PartialBuildFromFits_3D->RebinZ(FourParticleRebin);
	      }
	      //
	      for(int binx=1; binx<=TPFullWeight_FourParticle_2D[c1]->GetNbinsX(); binx++) {
		for(int biny=1; biny<=TPFullWeight_FourParticle_2D[c1]->GetNbinsY(); biny++) {
		  TPFullWeight_FourParticle_2D[c1]->SetBinError(binx, biny, 0);
		  if(binx!=4){
		    TPFullWeight_FourParticle_2D[c1]->SetBinContent(binx, biny, TPFullWeight_FourParticle_2D[c1]->GetBinContent(binx, biny) + TPFullWeight_FourParticle_2D[c1]->GetBinContent(4, biny));
		    TPNegFullWeight_FourParticle_2D[c1]->SetBinContent(binx, biny, TPNegFullWeight_FourParticle_2D[c1]->GetBinContent(binx, biny) + TPNegFullWeight_FourParticle_2D[c1]->GetBinContent(4, biny));
		    
		    if(InterpCorrection){// Interpolator correction
		      double q4 = TPFullWeight_FourParticle_2D[c1]->GetYaxis()->GetBinCenter(biny);
		      double InterCorr = 1.032 - 0.2452*q4;
		      TPFullWeight_FourParticle_2D[c1]->SetBinContent(binx, biny, InterCorr * TPFullWeight_FourParticle_2D[c1]->GetBinContent(binx, biny));
		      TPNegFullWeight_FourParticle_2D[c1]->SetBinContent(binx, biny, InterCorr * TPNegFullWeight_FourParticle_2D[c1]->GetBinContent(binx, biny));
		    }
		  }
		}
	      }
	      TString *proName=new TString(nameTPN4->Data()); TString *proNegName=new TString(nameNegTPN4->Data());
	      TString *proNameFitBuild=new TString(nameFitBuild->Data()); TString *proNamePartialFitBuild=new TString(namePartialFitBuild->Data());
	      proName->Append("_pro"); proNegName->Append("_pro"); proNameFitBuild->Append("_pro"); proNamePartialFitBuild->Append("_pro");
	      TPN_FourParticle[c1] = (TH1D*)TPFullWeight_FourParticle_2D[c1]->ProjectionY(proName->Data(), TPNbin, TPNbin);
	      //
	      proName->Append("_FullWeight"); proNegName->Append("_FullWeight");
	      TPFullWeight_FourParticle[c1] = (TH1D*)TPFullWeight_FourParticle_2D[c1]->ProjectionY(proName->Data(), Gbin, Gbin);
	      TPNegFullWeight_FourParticle[c1] = (TH1D*)TPNegFullWeight_FourParticle_2D[c1]->ProjectionY(proNegName->Data(), Gbin, Gbin);
	      proName->Append("_FullWeightDen"); proNegName->Append("_FullWeightDen");
	      TH1D *tempDen = (TH1D*)TPFullWeight_FourParticle_2D[c1]->ProjectionY(proName->Data(), 4, 4);
	      TH1D *tempDenNeg = (TH1D*)TPNegFullWeight_FourParticle_2D[c1]->ProjectionY(proNegName->Data(), 4, 4);
	      //
	      // Add Pos and Neg Weights
	      tempDen->Add(tempDenNeg);
	      TPFullWeight_FourParticle[c1]->Add(TPNegFullWeight_FourParticle[c1]);
	      //
	      //TPFullWeight_FourParticle[c1]->Add(tempDen);// now added above in Interp section
	      TPFullWeight_FourParticle[c1]->Divide(tempDen);
	      TPFullWeight_FourParticle[c1]->SetLineColor(1);
	      /*TString *ErrName=new TString(nameTPN4->Data());
	      ErrName->Append("Err");
	      TH2D *temperr2D = (TH2D*)MyList->FindObject(ErrName->Data());
	      TH1D *temperr = (TH1D*)temperr2D->ProjectionY("tesst",4,4);
	      temperr->Rebin(FourParticleRebin);
	      cout.precision(8);
	      cout<<temperr->GetBinContent(3)<<endl;
	      cout<<(temperr->GetBinContent(5) / tempDen->GetBinContent(5))<<"  "<<TPFullWeight_FourParticle[c1]->GetBinContent(5)<<endl;
	      */
	      if(c1==0){
		FullBuildFromFits = (TH1D*)FullBuildFromFits_3D->ProjectionZ(proNameFitBuild->Data(), c3FitType, c3FitType, Gbin, Gbin);
		TH1D *tempDen2 = (TH1D*)FullBuildFromFits_3D->ProjectionZ("tempDen2", c3FitType, c3FitType, 4, 4);
		tempDen2->Scale(1/100.);// It was filled 100 times with the same value
		FullBuildFromFits->Add(tempDen2);
		FullBuildFromFits->Divide(tempDen2);
		FullBuildFromFits->SetLineColor(1);
		//
		PartialBuildFromFits = (TH1D*)PartialBuildFromFits_3D->ProjectionZ(proNamePartialFitBuild->Data(), c3FitType, c3FitType, Gbin, Gbin);
		PartialBuildFromFits->Add(tempDen2);
		PartialBuildFromFits->Divide(tempDen2);
		PartialBuildFromFits->SetLineColor(2);
	      }
	    }
	  }// term 4-particle
	  
	}// c4
      }// c3
    }// c2
  }// c1
    
  TH1D *fTPNRejects3pion = (TH1D*)MyList->FindObject("fTPNRejects3pion2");
  TH1D *fTPNRejects4pion1 = (TH1D*)MyList->FindObject("fTPNRejects4pion1");

  //TH2D *c4QSFitNum2D = (TH2D*)MyList->FindObject("fc4QSFitNum");
  //TH2D *c4QSFitDen2D = (TH2D*)MyList->FindObject("fc4QSFitDen");
  //c4QSFitNum2D->RebinY(3);
  //c4QSFitDen2D->RebinY(3);

  cout<<"Done getting Histograms"<<endl;
  
  TF1 *Unity=new TF1("Unity","1",0,100);
  Unity->SetLineStyle(2);


  int ch1_2=0,ch2_2=0;
  int ch1_3=0,ch2_3=0,ch3_3=0;
  int ch1_4=0,ch2_4=0,ch3_4=0,ch4_4=0;
  
  if(SameCharge && CHARGE==+1) {ch1_2=1; ch2_2=1;   ch1_3=1; ch2_3=1; ch3_3=1;   ch1_4=1; ch2_4=1; ch3_4=1; ch4_4=1;}
  if(!SameCharge){
    ch1_2=0; ch2_2=1;
    //
    if(CHARGE==-1) { ch1_3=0; ch2_3=0; ch3_3=1; }
    else { ch1_3=0; ch2_3=1; ch3_3=1; }
    //
    if(CHARGE==-1 && MixedCharge4pionType==1) {ch1_4=0; ch2_4=0; ch3_4=0; ch4_4=1;}
    if(CHARGE==+1 && MixedCharge4pionType==1) {ch1_4=0; ch2_4=1; ch3_4=1; ch4_4=1;}
    if(MixedCharge4pionType==2) {ch1_4=0; ch2_4=0; ch3_4=1; ch4_4=1;}
  }
  
  ///////////////////////////////////////////////////////////
  // 2-pion section
  cout<<"2-pion section"<<endl;
  
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
  gPad->SetRightMargin(0.02); gPad->SetTopMargin(0.02);

  TLegend *legend1 = new TLegend(.6,.67,.87,.87,NULL,"brNDC");
  legend1->SetBorderSize(1);
  legend1->SetTextSize(.04);// small .03; large .036 
  legend1->SetFillColor(0);
  
  
  TH1D *TERM1_2=(TH1D*)TwoParticle[ch1_2][ch2_2][0]->Clone();
  TH1D *TERM2_2=(TH1D*)TwoParticle[ch1_2][ch2_2][1]->Clone();
  //
  if(SameCharge){
    TERM1_2->Multiply(MRC_SC_2[0]);
    TERM2_2->Multiply(MRC_SC_2[1]);
  }else {
    TERM1_2->Multiply(MRC_MC_2[0]);
    TERM2_2->Multiply(MRC_MC_2[1]);
  }
  
  TH1D *C2=(TH1D*)TERM1_2->Clone();
  C2->Divide(TERM2_2);
  C2->GetXaxis()->SetRangeUser(0,0.1);
  C2->SetMinimum(0.98); C2->SetMaximum(2.5);
  /*C2->GetXaxis()->SetTitleSize(0.06);
  C2->GetXaxis()->SetLabelSize(0.06);
  C2->GetYaxis()->SetTitleSize(0.06);
  C2->GetYaxis()->SetLabelSize(0.06);
  C2->GetXaxis()->SetTitleOffset(1.05);
  C2->GetYaxis()->SetTitleOffset(1.1);
  C2->GetXaxis()->SetNdivisions(606);
  C2->GetYaxis()->SetNdivisions(505);*/
  C2->Draw();

  TH1D *C2QS=(TH1D*)TERM1_2->Clone();
  C2QS->Add(TERM2_2, -(1-TwoFrac));
  C2QS->Scale(1/TwoFrac);
  for(int bin=1; bin<=C2QS->GetNbinsX(); bin++) {
    Float_t K2 = 1.0;
    if(FSICorrection) K2 = FSICorrelation(ch1_2, ch2_2, C2QS->GetXaxis()->GetBinCenter(bin));
    C2QS->SetBinContent(bin, C2QS->GetBinContent(bin)/K2);
    C2QS->SetBinError(bin, C2QS->GetBinError(bin)/K2);
  }
  C2QS->Divide(TERM2_2);
  if(SameCharge && MuonCorrection) C2QS->Multiply(C2muonCorrectionSC);
  if(!SameCharge && MuonCorrection) C2QS->Multiply(C2muonCorrectionMC);
  C2QS->SetMarkerColor(2); C2QS->SetLineColor(2);
  if(!MCcase) C2QS->Draw("same");

  // Print 2-pion values
  /*for(int bin=1; bin<=20; bin++){
    cout<<C2QS->GetBinContent(bin)<<", ";
  }
  cout<<endl;
  for(int bin=1; bin<=20; bin++){
    cout<<C2QS->GetBinError(bin)<<", ";
  }
  cout<<endl;*/
  double C2refPoints_0p65[20]={-0.218603, 1.00087, 1.00215, 1.00209, 1.00202, 1.00215, 1.00162, 1.00129, 1.00118, 1.00076, 1.00049, 1.00021, 0.999996, 0.999733, 0.999622, 0.999554, 0.999764, 0.999679, 0.999677, 0.999641};
  double C2refPoints_e_0p65[20]={0.309152, 0.000473757, 0.000305599, 0.00022603, 0.000179596, 0.000149295, 0.000128009, 0.000112371, 0.00010047, 9.11363e-05, 8.37061e-05, 7.76871e-05, 7.27445e-05, 6.86217e-05, 6.5155e-05, 6.22148e-05, 5.9718e-05, 5.75613e-05, 5.57066e-05, 5.41028e-05};
  double C2refPoints_0p75[20]={-0.135326, 0.968257, 0.984681, 0.991465, 0.99503, 0.997249, 0.998093, 0.998645, 0.999135, 0.999172, 0.999237, 0.999213, 0.999197, 0.999095, 0.999105, 0.999123, 0.999385, 0.999366, 0.999408, 0.999411};
  double C2refPoints_e_0p75[20]={0.19138, 0.000422359, 0.000272193, 0.000201183, 0.000159788, 0.000132797, 0.000113845, 9.99271e-05, 8.93377e-05, 8.10334e-05, 7.44239e-05, 6.90702e-05, 6.46743e-05, 6.10077e-05, 5.79248e-05, 5.53103e-05, 5.30903e-05, 5.11725e-05, 4.95234e-05, 4.80974e-05};
  TH1D *C2ref_0p65=(TH1D*)C2QS->Clone();
  TH1D *C2ref_0p75=(TH1D*)C2QS->Clone();
  for(int bin=1; bin<=20; bin++){
    C2ref_0p65->SetBinContent(bin, C2refPoints_0p65[bin-1]);
    C2ref_0p65->SetBinError(bin, C2refPoints_e_0p65[bin-1]);
    C2ref_0p75->SetBinContent(bin, C2refPoints_0p75[bin-1]);
    C2ref_0p75->SetBinError(bin, C2refPoints_e_0p75[bin-1]);
  }
  C2ref_0p65->SetMarkerColor(4); C2ref_0p65->SetLineColor(4);
  C2ref_0p75->SetMarkerColor(6); C2ref_0p75->SetLineColor(6);
  //C2ref_0p65->Draw("same");
  //C2ref_0p75->Draw("same");




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
  TLegend *legend2 = (TLegend*)legend1->Clone();
  legend2->SetX1(0.45); legend2->SetX2(0.95);  legend2->SetY1(0.6); legend2->SetY2(0.95);

  TH1D *TERM1_3=(TH1D*)ThreeParticle[ch1_3][ch2_3][ch3_3][0]->Clone();
  TH1D *TERM2_3=(TH1D*)ThreeParticle[ch1_3][ch2_3][ch3_3][1]->Clone();
  TH1D *TERM3_3=(TH1D*)ThreeParticle[ch1_3][ch2_3][ch3_3][2]->Clone();
  TH1D *TERM4_3=(TH1D*)ThreeParticle[ch1_3][ch2_3][ch3_3][3]->Clone();
  TH1D *TERM5_3=(TH1D*)ThreeParticle[ch1_3][ch2_3][ch3_3][4]->Clone();
  //
  
  if(SameCharge){
    TERM1_3->Multiply(MRC_SC_3[0]);
    TERM2_3->Multiply(MRC_SC_3[1]);
    TERM3_3->Multiply(MRC_SC_3[1]);
    TERM4_3->Multiply(MRC_SC_3[1]);
    TERM5_3->Multiply(MRC_SC_3[2]);
    //cout<<MRC_SC_3[2]->GetBinContent(2)<<"  "<<MRC_SC_3[2]->GetBinContent(3)<<"  "<<MRC_SC_3[2]->GetBinContent(4)<<endl;
  }else{
    if(CHARGE==-1){// --+ but MRC is stored as -++
      TERM1_3->Multiply(MRC_MC_3[0]);
      TERM2_3->Multiply(MRC_MC_3[2]);
      TERM3_3->Multiply(MRC_MC_3[1]);
      TERM4_3->Multiply(MRC_MC_3[1]);
      TERM5_3->Multiply(MRC_MC_3[3]);
    }else{
      TERM1_3->Multiply(MRC_MC_3[0]);
      TERM2_3->Multiply(MRC_MC_3[1]);
      TERM3_3->Multiply(MRC_MC_3[1]);
      TERM4_3->Multiply(MRC_MC_3[2]);
      TERM5_3->Multiply(MRC_MC_3[3]);
    }
  }
 
  
  TH1D *N3QS = (TH1D*)TERM1_3->Clone();
  N3QS->Add(TERM5_3, -f31);
  N3QS->Add(TERM2_3, -f32);
  N3QS->Add(TERM3_3, -f32);
  N3QS->Add(TERM4_3, -f32);
  N3QS->Scale(1/f33);
  N3QS->Multiply(K3avg[ch1_3][ch2_3][ch3_3][0]);

  
  TH1D *C3QS = (TH1D*)N3QS->Clone();
  C3QS->Divide(TERM5_3);
  if(SameCharge) C3QS->Multiply(C3muonCorrectionSC[0]);
  else C3QS->Multiply(C3muonCorrectionMC[0]);
  
  C3QS->SetMarkerStyle(20);
  C3QS->SetMarkerColor(4);
  C3QS->GetXaxis()->SetRangeUser(0,0.1);
  C3QS->SetMinimum(0.9); C3QS->SetMaximum(5.0);
  
  TH1D *n3QS = (TH1D*)N3QS->Clone();
  TH1D *c3QS = (TH1D*)N3QS->Clone();
  for(int bin=1; bin<=n3QS->GetNbinsX(); bin++){
    if(n3QS->GetBinContent(bin)==0) continue;
    double MuonCorr1=1.0, MuonCorr2=1.0, MuonCorr3=1.0, MuonCorr4=1.0;
    if(SameCharge){
      MuonCorr1 = C3muonCorrectionSC[0]->GetBinContent(bin);
      MuonCorr2 = C3muonCorrectionSC[1]->GetBinContent(bin);
      MuonCorr3 = MuonCorr2;
      MuonCorr4 = MuonCorr2;
    }else if(CHARGE==-1){
      MuonCorr1 = C3muonCorrectionMC[0]->GetBinContent(bin);
      MuonCorr2 = C3muonCorrectionMC[2]->GetBinContent(bin);// --
      MuonCorr3 = C3muonCorrectionMC[1]->GetBinContent(bin);// -+
      MuonCorr4 = MuonCorr3;// -+
    }else{
      MuonCorr1 = C3muonCorrectionMC[0]->GetBinContent(bin);
      MuonCorr2 = C3muonCorrectionMC[1]->GetBinContent(bin);// -+
      MuonCorr3 = MuonCorr2;
      MuonCorr4 = C3muonCorrectionMC[2]->GetBinContent(bin);
    }
    double value = n3QS->GetBinContent(bin) * MuonCorr1;
    value -= ( TERM2_3->GetBinContent(bin) - (1-TwoFrac)*TERM5_3->GetBinContent(bin)) * K3avg[ch1_3][ch2_3][ch3_3][1]->GetBinContent(bin) / TwoFrac * MuonCorr2;
    value -= ( TERM3_3->GetBinContent(bin) - (1-TwoFrac)*TERM5_3->GetBinContent(bin)) * K3avg[ch1_3][ch2_3][ch3_3][2]->GetBinContent(bin) / TwoFrac * MuonCorr3;
    value -= ( TERM4_3->GetBinContent(bin) - (1-TwoFrac)*TERM5_3->GetBinContent(bin)) * K3avg[ch1_3][ch2_3][ch3_3][3]->GetBinContent(bin) / TwoFrac * MuonCorr4;
    //value -= TERM5_3->GetBinContent(bin)*(MuonCorr1 - MuonCorr2 - MuonCorr3 - MuonCorr4);
    value += 2*TERM5_3->GetBinContent(bin);
    //
    n3QS->SetBinContent(bin, value);
    c3QS->SetBinContent(bin, value + TERM5_3->GetBinContent(bin));
  }
  c3QS->Divide(TERM5_3);

  C3QS->Draw();
  c3QS->Draw("same");
  //int lowBin = C3QS->GetXaxis()->FindBin(Cutoff_FullWeight_Q3[Mbin]);
  //int highBin = C3QS->GetXaxis()->FindBin(Cutoff_FullWeight_Q3[Mbin]);
  //double SF=C3QS->Integral(lowBin, highBin); SF /= TPFullWeight_ThreeParticle[ch1_3]->Integral(lowBin, highBin);
  //TPFullWeight_ThreeParticle[ch1_3]->Scale(SF);
  //
  if(SameCharge) TPFullWeight_ThreeParticle[ch1_3]->Draw("same");
  
  //cout<<TPFullWeight_ThreeParticle[ch1_3]->GetBinContent(2)<<endl;
  //TH1D *C3raw = (TH1D*)TERM1_3->Clone();
  //C3raw->Divide(TERM5_3);
  //C3raw->SetMarkerColor(4);
  //C3raw->Draw("same");

  legend2->AddEntry(C3QS,"C_{3}^{QS}","p");
  legend2->AddEntry(c3QS,"c_{3}^{QS}{2-pion removal}","p");
  if(SameCharge) legend2->AddEntry(TPFullWeight_ThreeParticle[ch1_3],"C_{3}^{QS} built","l");
  legend2->Draw("same");
  
  
  ///////////////////////////////////////
  // C3QS built ratio
  TCanvas *can2_2 = new TCanvas("can2_2", "can2_2",1200,53,700,500);
  can2_2->SetHighLightColor(2);
  can2_2->Range(-0.7838115,-1.033258,0.7838115,1.033258);
  gStyle->SetOptFit(0111);
  can2_2->SetFillColor(10);
  can2_2->SetBorderMode(0);
  can2_2->SetBorderSize(2);
  can2_2->SetGridx();
  can2_2->SetGridy();
  can2_2->SetFrameFillColor(0);
  can2_2->SetFrameBorderMode(0);
  can2_2->SetFrameBorderMode(0);
  gPad->SetRightMargin(0.02); gPad->SetTopMargin(0.02);
  TLegend *legend2_2 = (TLegend*)legend1->Clone();
  legend2_2->SetX1(0.25); legend2_2->SetX2(0.6);  legend2_2->SetY1(0.8); legend2_2->SetY2(0.95);

  TH1D *C3QSratio = (TH1D*)C3QS->Clone();
 
  //TH1D *Chi2NDF_PointSize_3 = new TH1D("Chi2NDF_PointSize_3","",40,-0.5,39.5);
  //TH1D *Chi2NDF_FullSize_3 = new TH1D("Chi2NDF_FullSize_3","",40,-0.5,39.5);
  TH1D *Chi2NDF_PointSize_3 = new TH1D("Chi2NDF_PointSize_3","",100,-0.5,99.5);
  TH1D *Chi2NDF_FullSize_3 = new TH1D("Chi2NDF_FullSize_3","",100,-0.5,99.5);
  Chi2NDF_PointSize_3->SetLineColor(4); Chi2NDF_FullSize_3->SetLineColor(2);
  Chi2NDF_PointSize_3->SetMarkerColor(4); Chi2NDF_FullSize_3->SetMarkerColor(2);
  Chi2NDF_PointSize_3->GetXaxis()->SetTitle("Coherent fraction (%)"); Chi2NDF_PointSize_3->GetYaxis()->SetTitle("#sqrt{#chi^{2}}");
  
  if(SameCharge && CollisionType==0){
    
    TH1D *tempDen = (TH1D*)TPFullWeight_ThreeParticle_2D[ch1_3]->ProjectionY("TPFullWeight3_Den", 4, 4);
    TH1D *tempDenNeg = (TH1D*)TPNegFullWeight_ThreeParticle_2D[ch1_3]->ProjectionY("TPNegFullWeight3_Den", 4, 4);
    tempDen->Add(tempDenNeg);// Add Pos and Neg Den

    for(int binG=5; binG<=104; binG++){// was 44
      TString *proName=new TString("TPFullWeight3_");
      *proName += binG;
      TH1D *tempNum = (TH1D*)TPFullWeight_ThreeParticle_2D[ch1_3]->ProjectionY(proName->Data(), binG, binG);
      proName->Append("_Neg");
      TH1D *tempNumNeg = (TH1D*)TPNegFullWeight_ThreeParticle_2D[ch1_3]->ProjectionY(proName->Data(), binG, binG);
      //
      // Add Pos and Neg Num
      tempNum->Add(tempNumNeg);
      //
      //tempNum->Add(tempDen);// now added in histogram read section
      tempNum->Divide(tempDen);
      //lowBin = C3QS->GetXaxis()->FindBin(Cutoff_FullWeight_Q3[Mbin]);
      //highBin = C3QS->GetXaxis()->FindBin(Cutoff_FullWeight_Q3[Mbin]);
      //SF=C3QS->Integral(lowBin, highBin);
      //SF /= tempNum->Integral(lowBin, highBin);
      //tempNum->Scale(SF);
      
      double Chi2=0;
      double NDF=0;
      for(int binQ3=3; binQ3<=3; binQ3++){
	if(tempNum->GetBinContent(binQ3) <=0) continue;
	double value = C3QS->GetBinContent(binQ3) - tempNum->GetBinContent(binQ3);
	double err = C3QS->GetBinError(binQ3);
	if(err <=0) continue;

	Chi2 += pow(value / err,2);
	//
	NDF += 1;
      }
      //if(binG<25) {Chi2NDF_PointSize_3->SetBinContent(1 + 2*(binG-5), sqrt(Chi2)/NDF); Chi2NDF_PointSize_3->SetBinError(1 + 2*(binG-5), 0.001);}
      //else {Chi2NDF_FullSize_3->SetBinContent(1 + 2*(binG-25), sqrt(Chi2)/NDF); Chi2NDF_FullSize_3->SetBinError(1 + 2*(binG-25), 0.001);}
      if(binG<55) {Chi2NDF_PointSize_3->SetBinContent(1 + 2*(binG-5), sqrt(Chi2)/NDF); Chi2NDF_PointSize_3->SetBinError(1 + 2*(binG-5), 0.001);}
      else {Chi2NDF_FullSize_3->SetBinContent(1 + 2*(binG-55), sqrt(Chi2)/NDF); Chi2NDF_FullSize_3->SetBinError(1 + 2*(binG-55), 0.001);}
      //if(NDF>0) cout<<binG<<"  "<<sqrt(Chi2)/NDF<<endl;
    }
    Chi2NDF_PointSize_3->SetMarkerStyle(20); Chi2NDF_FullSize_3->SetMarkerStyle(20);
    
    C3QSratio->Divide(TPFullWeight_ThreeParticle[ch1_3]);
    //
    C3QSratio->SetMinimum(0.94); C3QSratio->SetMaximum(1.02);
    C3QSratio->GetYaxis()->SetTitleOffset(1.2);
    C3QSratio->GetYaxis()->SetTitle("C_{3}^{QS} / C_{3}(built)");
    C3QSratio->Draw();
    Unity->Draw("same");
    
    Chi2NDF_PointSize_3->Draw();
    Chi2NDF_FullSize_3->Draw("same");
    legend2_2->AddEntry(Chi2NDF_PointSize_3,"R_{coh}=0 (Point Source)","p");
    legend2_2->AddEntry(Chi2NDF_FullSize_3,"R_{coh}=R_{ch}","p");
    legend2_2->Draw("same");
  }
  
  // r3
  TH1D *r3;
  if(SameCharge && CollisionType==0){
    r3 = (TH1D*)n3QS->Clone();
    TPN_ThreeParticle[ch1_3]->Multiply(MRC_SC_3[2]);
    r3->Divide(TPN_ThreeParticle[ch1_3]);
    r3->GetXaxis()->SetRangeUser(0,0.08);
    r3->SetMinimum(0.5); r3->SetMaximum(2.5);
    r3->GetYaxis()->SetTitle("r_{3}");
    //
    r3->Draw();
    //TPN_ThreeParticle[ch1_3]->Draw();
    //fTPNRejects3pion->SetLineColor(2);
    //fTPNRejects3pion->Draw("same");
  }
  

  // Print 3-pion points
  for(int bin=1; bin<=10; bin++){
    //cout<<C3QS->GetBinContent(bin)<<", ";
    //cout<<c3QS->GetBinContent(bin)<<", ";
    //cout<<TPFullWeight_ThreeParticle[ch1_3]->GetBinContent(bin)<<", ";
  }
  //cout<<endl;
  for(int bin=1; bin<=10; bin++){
    //cout<<C3QS->GetBinError(bin)<<", ";
    //cout<<c3QS->GetBinError(bin)<<", ";
  }
  //cout<<endl;


  ////////////////////////////////////////////////////////////////
  // 4-pion 
  cout<<"4-pion section"<<endl;
  
  TCanvas *can3 = new TCanvas("can3", "can3",11,600,700,500);// 60 was 600
  can3->SetHighLightColor(2);
  can3->Range(-0.7838115,-1.033258,0.7838115,1.033258);
  gStyle->SetOptFit(0111);
  can3->SetFillColor(10);
  can3->SetBorderMode(0);
  can3->SetBorderSize(2);
  can3->SetGridx();
  can3->SetGridy();
  can3->SetFrameFillColor(0);
  can3->SetFrameBorderMode(0);
  can3->SetFrameBorderMode(0);
  gPad->SetRightMargin(0.02); gPad->SetTopMargin(0.02);
  TLegend *legend3 = (TLegend*)legend1->Clone();
  legend3->SetX1(0.45); legend3->SetX2(0.98);  legend3->SetY1(0.6); legend3->SetY2(0.95);

  TH1D *TERMS_4[13];
  for(int index=0; index<13; index++){
    TERMS_4[index] = (TH1D*)FourParticle[ch1_4][ch2_4][ch3_4][ch4_4][index]->Clone();
    //
    TH1D *MRC_temp; 
    if(SameCharge){
      if(index==0) MRC_temp = (TH1D*)MRC_SC_4[0]->Clone();
      else if(index<=4) MRC_temp = (TH1D*)MRC_SC_4[1]->Clone();
      else if(index<=10) MRC_temp = (TH1D*)MRC_SC_4[2]->Clone();
      else if(index==11) MRC_temp = (TH1D*)MRC_SC_4[3]->Clone();
      else MRC_temp = (TH1D*)MRC_SC_4[4]->Clone();    
    }else if(CHARGE==-1 && MixedCharge4pionType==1){
      if(index==0) MRC_temp = (TH1D*)MRC_MC1_4[0]->Clone();//0
      else if(index!=1 && index<=4) MRC_temp = (TH1D*)MRC_MC1_4[1]->Clone();//1
      else if(index==1) MRC_temp = (TH1D*)MRC_MC1_4[2]->Clone();//2
      else if(index==7 || index==9 || index==10) MRC_temp = (TH1D*)MRC_MC1_4[3]->Clone();//3
      else if(index!=12) MRC_temp = (TH1D*)MRC_MC1_4[4]->Clone();//4
      else MRC_temp = (TH1D*)MRC_MC1_4[5]->Clone();//5
    }else if(CHARGE==+1 && MixedCharge4pionType==1){
      if(index==0) MRC_temp = (TH1D*)MRC_MC1_4[0]->Clone();
      else if(index<=3) MRC_temp = (TH1D*)MRC_MC1_4[1]->Clone();
      else if(index==4) MRC_temp = (TH1D*)MRC_MC1_4[2]->Clone();
      else if(index==5 || index==6 || index==7) MRC_temp = (TH1D*)MRC_MC1_4[3]->Clone();
      else if(index!=12) MRC_temp = (TH1D*)MRC_MC1_4[4]->Clone();
      else MRC_temp = (TH1D*)MRC_MC1_4[5]->Clone();
    }else{// --++
      if(index==0) MRC_temp = (TH1D*)MRC_MC2_4[0]->Clone();
      else if(index<=4) MRC_temp = (TH1D*)MRC_MC2_4[1]->Clone();
      else if(index==5 || index==10) MRC_temp = (TH1D*)MRC_MC2_4[2]->Clone();//2
      else if(index<=9) MRC_temp = (TH1D*)MRC_MC2_4[3]->Clone();// 3
      else if(index==11) MRC_temp = (TH1D*)MRC_MC2_4[4]->Clone();// 4
      else MRC_temp = (TH1D*)MRC_MC2_4[5]->Clone();
    }
    //
    TERMS_4[index]->Multiply(MRC_temp);
  }
  
  TH1D *C4raw=(TH1D*)TERMS_4[2]->Clone();
  C4raw->Divide(TERMS_4[12]);
  C4raw->GetXaxis()->SetRangeUser(0,0.2);
  

  TH1D *N4QS=(TH1D*)TERMS_4[0]->Clone();
  //
  N4QS->Add(TERMS_4[1], -f43); 
  N4QS->Add(TERMS_4[2], -f43);
  N4QS->Add(TERMS_4[3], -f43);
  N4QS->Add(TERMS_4[4], -f43);
  //
  N4QS->Add(TERMS_4[5], -f42);
  N4QS->Add(TERMS_4[6], -f42);
  N4QS->Add(TERMS_4[7], -f42);
  N4QS->Add(TERMS_4[8], -f42);
  N4QS->Add(TERMS_4[9], -f42);
  N4QS->Add(TERMS_4[10], -f42);
  //
  N4QS->Add(TERMS_4[12], -f41);
  //
  N4QS->Scale(1/f44);
  //K4avg[ch1_4][ch2_4][ch3_4][ch4_4][0]->Scale(1.005);// K scale variation
  N4QS->Multiply(K4avg[ch1_4][ch2_4][ch3_4][ch4_4][0]);
  //N4QS->Divide(K4avg[ch1_4][ch2_4][ch3_4][ch4_4][1]);// K factorization test (MC1)
  //N4QS->Divide(K4avg[ch1_4][ch2_4][ch3_4][ch4_4][11]);// K factorization test (MC2)
  
  TH1D *C4QS=(TH1D*)N4QS->Clone();
  C4QS->GetXaxis()->SetRangeUser(0,0.179);
  C4QS->SetMinimum(0.5);
  C4QS->SetMarkerColor(4);
 
  
  //
  /*TH1D *C4QS_basic=(TH1D*)TERMS_4[0]->Clone();
  for(int ii=1; ii<=C4QS_basic->GetNbinsX(); ii++){
    double value = C4QS_basic->GetBinContent(ii) - TERMS_4[12]->GetBinContent(ii);
    value /= f44;
    value += TERMS_4[12]->GetBinContent(ii);
    C4QS_basic->SetBinContent(ii, value);
    }
  C4QS_basic->Divide(TERMS_4[12]);
  //C4QS_basic->Multiply(K4avg[ch1_4][ch2_4][ch3_4][ch4_4][0]);
  C4QS_basic->SetMarkerColor(1);
  //C4QS_basic->Draw("same");
  */
  //
  
  TH1D *C4_pieces[12];
  for(int i=0; i<12; i++){
    C4_pieces[i]=(TH1D*)TERMS_4[i]->Clone();
    C4_pieces[i]->Divide(TERMS_4[12]);
    
    if(i==0) C4_pieces[i]->SetMarkerColor(1);
    else if(i<5) C4_pieces[i]->SetMarkerColor(2);
    else if(i<10) C4_pieces[i]->SetMarkerColor(3);
    else C4_pieces[i]->SetMarkerColor(4);
  }
  C4_pieces[0]->GetXaxis()->SetRangeUser(0,0.2);
  C4_pieces[0]->SetMinimum(0.5); C4_pieces[0]->SetMaximum(3); 
  //C4_pieces[0]->Draw();
  //C4_pieces[1]->Draw("same");
  //C4_pieces[5]->Draw("same");
  //C4_pieces[11]->Draw("same");
  //
  
  TH1D *C4_2pairs_FromSquares = (TH1D*)TERMS_4[5]->Clone();
  C4_2pairs_FromSquares->Divide(TERMS_4[12]);
  for(int bin=1; bin<=C4_2pairs_FromSquares->GetNbinsX(); bin++){
    double value = C4_2pairs_FromSquares->GetBinContent(bin);
    value -= 1;
    value *= value;
    value += 1;
    C4_2pairs_FromSquares->SetBinContent(bin, value);
  }
  C4_2pairs_FromSquares->SetMarkerColor(6);
  //C4_2pairs_FromSquares->Draw("same");
  //
  TH1D *C4QS_2pairs=(TH1D*)TERMS_4[11]->Clone();
  C4QS_2pairs->Add(TERMS_4[12], -pow(1-TwoFrac,2));// N1^4
  C4QS_2pairs->Add(TERMS_4[5], -2*TwoFrac*(1-TwoFrac));// N2 * N1^2
  C4QS_2pairs->Scale(1/pow(TwoFrac,2));
  C4QS_2pairs->Multiply(K4avg[ch1_4][ch2_4][ch3_4][ch4_4][11]);
  C4QS_2pairs->Divide(TERMS_4[12]);
  C4QS_2pairs->SetMarkerColor(6);
  //C4QS_2pairs->Draw("same");
  //
  

  
  TH1D *n4QS = (TH1D*)N4QS->Clone();
  TH1D *c4QS = (TH1D*)N4QS->Clone();
  TH1D *c4QS_RemovalStage1 = (TH1D*)N4QS->Clone();
  TH1D *c4QS_RemovalStage2 = (TH1D*)N4QS->Clone();
  TH1D *c4QS_RemovalStage3 = (TH1D*)N4QS->Clone();


  for(int bin=1; bin<=n4QS->GetNbinsX(); bin++){
    if(n4QS->GetBinContent(bin)==0) continue;
    bool exitCode=kFALSE;
    for(int ii=0; ii<13; ii++) {if(TERMS_4[ii]->GetBinContent(bin) < 1) exitCode=kTRUE;}
    if(exitCode) continue;
    //cout<<TERMS_4[0]->GetBinContent(bin)<<endl;
    //
    double N4QSvalue = N4QS->GetBinContent(bin);
    double SubtractionTerm[11] = {0};
    double ErrorTerm[12] = {0};
    ErrorTerm[0] = n4QS->GetBinError(bin);
    //
    // 3-pion terms
    if(SameCharge || MixedCharge4pionType==1){
      SubtractionTerm[0] = (TERMS_4[1]->GetBinContent(bin) - f32*TERMS_4[5]->GetBinContent(bin)  - f32*TERMS_4[6]->GetBinContent(bin)  - f32*TERMS_4[8]->GetBinContent(bin) - f31*TERMS_4[12]->GetBinContent(bin)) / f33;// 5 6 8
      SubtractionTerm[0] *= K4avg[ch1_4][ch2_4][ch3_4][ch4_4][1]->GetBinContent(bin);
      ErrorTerm[1] = sqrt(TERMS_4[1]->GetBinContent(bin) + f32*f32*TERMS_4[5]->GetBinContent(bin)  + f32*f32*TERMS_4[6]->GetBinContent(bin) + f32*f32*TERMS_4[8]->GetBinContent(bin) + f31*f31*TERMS_4[12]->GetBinContent(bin)) / f33;
      ErrorTerm[1] *= K4avg[ch1_4][ch2_4][ch3_4][ch4_4][1]->GetBinContent(bin);
      //
      SubtractionTerm[1] = (TERMS_4[2]->GetBinContent(bin) - f32*TERMS_4[5]->GetBinContent(bin)  - f32*TERMS_4[7]->GetBinContent(bin)  - f32*TERMS_4[9]->GetBinContent(bin) - f31*TERMS_4[12]->GetBinContent(bin)) / f33;
      SubtractionTerm[1] *= K4avg[ch1_4][ch2_4][ch3_4][ch4_4][2]->GetBinContent(bin);
      ErrorTerm[2] = sqrt(TERMS_4[2]->GetBinContent(bin) + f32*f32*TERMS_4[5]->GetBinContent(bin)  + f32*f32*TERMS_4[7]->GetBinContent(bin) + f32*f32*TERMS_4[9]->GetBinContent(bin) + f31*f31*TERMS_4[12]->GetBinContent(bin)) / f33;
      ErrorTerm[2] *= K4avg[ch1_4][ch2_4][ch3_4][ch4_4][2]->GetBinContent(bin);
      //
      SubtractionTerm[2] = (TERMS_4[3]->GetBinContent(bin) - f32*TERMS_4[6]->GetBinContent(bin)  - f32*TERMS_4[7]->GetBinContent(bin)  - f32*TERMS_4[10]->GetBinContent(bin) - f31*TERMS_4[12]->GetBinContent(bin)) / f33;
      SubtractionTerm[2] *= K4avg[ch1_4][ch2_4][ch3_4][ch4_4][3]->GetBinContent(bin);
      ErrorTerm[3] = sqrt(TERMS_4[3]->GetBinContent(bin) + f32*f32*TERMS_4[6]->GetBinContent(bin)  + f32*f32*TERMS_4[7]->GetBinContent(bin) + f32*f32*TERMS_4[10]->GetBinContent(bin) + f31*f31*TERMS_4[12]->GetBinContent(bin)) / f33;
      ErrorTerm[3] *= K4avg[ch1_4][ch2_4][ch3_4][ch4_4][3]->GetBinContent(bin);
      //
      SubtractionTerm[3] = (TERMS_4[4]->GetBinContent(bin) - f32*TERMS_4[8]->GetBinContent(bin)  - f32*TERMS_4[9]->GetBinContent(bin)  - f32*TERMS_4[10]->GetBinContent(bin) - f31*TERMS_4[12]->GetBinContent(bin)) / f33;
      SubtractionTerm[3] *= K4avg[ch1_4][ch2_4][ch3_4][ch4_4][4]->GetBinContent(bin);
      ErrorTerm[4] = sqrt(TERMS_4[4]->GetBinContent(bin) + f32*f32*TERMS_4[8]->GetBinContent(bin)  + f32*f32*TERMS_4[9]->GetBinContent(bin) + f32*f32*TERMS_4[10]->GetBinContent(bin) + f31*f31*TERMS_4[12]->GetBinContent(bin)) / f33;
      ErrorTerm[4] *= K4avg[ch1_4][ch2_4][ch3_4][ch4_4][4]->GetBinContent(bin);
      //
    }
    // 2-pion terms
    SubtractionTerm[4] = ( TERMS_4[5]->GetBinContent(bin) - (1-TwoFrac)*TERMS_4[12]->GetBinContent(bin)) / TwoFrac * K4avg[ch1_4][ch2_4][ch3_4][ch4_4][5]->GetBinContent(bin);
    ErrorTerm[5] = sqrt( TERMS_4[5]->GetBinContent(bin) + pow(1-TwoFrac,2)*TERMS_4[12]->GetBinContent(bin)) * K4avg[ch1_4][ch2_4][ch3_4][ch4_4][5]->GetBinContent(bin) / TwoFrac;
    //
    SubtractionTerm[5] = ( TERMS_4[6]->GetBinContent(bin) - (1-TwoFrac)*TERMS_4[12]->GetBinContent(bin)) * K4avg[ch1_4][ch2_4][ch3_4][ch4_4][6]->GetBinContent(bin) / TwoFrac;
    ErrorTerm[6] = sqrt( TERMS_4[6]->GetBinContent(bin) + pow(1-TwoFrac,2)*TERMS_4[12]->GetBinContent(bin)) * K4avg[ch1_4][ch2_4][ch3_4][ch4_4][6]->GetBinContent(bin) / TwoFrac;
    //
    SubtractionTerm[6] = ( TERMS_4[7]->GetBinContent(bin) - (1-TwoFrac)*TERMS_4[12]->GetBinContent(bin)) * K4avg[ch1_4][ch2_4][ch3_4][ch4_4][7]->GetBinContent(bin) / TwoFrac;
    ErrorTerm[7] = sqrt( TERMS_4[7]->GetBinContent(bin) + pow(1-TwoFrac,2)*TERMS_4[12]->GetBinContent(bin)) * K4avg[ch1_4][ch2_4][ch3_4][ch4_4][7]->GetBinContent(bin) / TwoFrac;
    //
    SubtractionTerm[7] = ( TERMS_4[8]->GetBinContent(bin) - (1-TwoFrac)*TERMS_4[12]->GetBinContent(bin)) * K4avg[ch1_4][ch2_4][ch3_4][ch4_4][8]->GetBinContent(bin) / TwoFrac;
    ErrorTerm[8] = sqrt( TERMS_4[8]->GetBinContent(bin) + pow(1-TwoFrac,2)*TERMS_4[12]->GetBinContent(bin)) * K4avg[ch1_4][ch2_4][ch3_4][ch4_4][8]->GetBinContent(bin) / TwoFrac;
    //
    SubtractionTerm[8] = ( TERMS_4[9]->GetBinContent(bin) - (1-TwoFrac)*TERMS_4[12]->GetBinContent(bin)) * K4avg[ch1_4][ch2_4][ch3_4][ch4_4][9]->GetBinContent(bin) / TwoFrac;
    ErrorTerm[9] = sqrt( TERMS_4[9]->GetBinContent(bin) + pow(1-TwoFrac,2)*TERMS_4[12]->GetBinContent(bin)) * K4avg[ch1_4][ch2_4][ch3_4][ch4_4][9]->GetBinContent(bin) / TwoFrac;
    //
    SubtractionTerm[9] = ( TERMS_4[10]->GetBinContent(bin) - (1-TwoFrac)*TERMS_4[12]->GetBinContent(bin)) * K4avg[ch1_4][ch2_4][ch3_4][ch4_4][10]->GetBinContent(bin) / TwoFrac;
    ErrorTerm[10] = sqrt( TERMS_4[10]->GetBinContent(bin) + pow(1-TwoFrac,2)*TERMS_4[12]->GetBinContent(bin)) * K4avg[ch1_4][ch2_4][ch3_4][ch4_4][10]->GetBinContent(bin) / TwoFrac;
    //
    
    //
    // 2 2-pion terms: cos(q12*x12)*cos(q34*x34)
    SubtractionTerm[10] = TERMS_4[11]->GetBinContent(bin);// N22
    SubtractionTerm[10] -= 2*(1-TwoFrac) * TERMS_4[5]->GetBinContent(bin);
    SubtractionTerm[10] += pow(1-TwoFrac,2) * TERMS_4[12]->GetBinContent(bin);
    SubtractionTerm[10] /= pow(TwoFrac,2);
    SubtractionTerm[10] *= K4avg[ch1_4][ch2_4][ch3_4][ch4_4][11]->GetBinContent(bin);
    SubtractionTerm[10] *= C4muonCorrectionSC[3]->GetBinContent(bin);
    double intermed = ( TERMS_4[5]->GetBinContent(bin) - (1-TwoFrac)*TERMS_4[12]->GetBinContent(bin)) / TwoFrac;
    intermed *= K4avg[ch1_4][ch2_4][ch3_4][ch4_4][5]->GetBinContent(bin);
    intermed = intermed * C4muonCorrectionSC[2]->GetBinContent(bin);
    SubtractionTerm[10] -= 2*intermed;
    SubtractionTerm[10] += 2*TERMS_4[12]->GetBinContent(bin);
    //
    ErrorTerm[11] = TERMS_4[11]->GetBinContent(bin);
    ErrorTerm[11] += pow(1-TwoFrac,4) * TERMS_4[12]->GetBinContent(bin);
    ErrorTerm[11] += pow(2*(1-TwoFrac)*TwoFrac,2) * TERMS_4[5]->GetBinContent(bin);
    ErrorTerm[11] = sqrt(ErrorTerm[11]);
    ErrorTerm[11] /= pow(TwoFrac,2);
    ErrorTerm[11] *= K4avg[ch1_4][ch2_4][ch3_4][ch4_4][11]->GetBinContent(bin);
    ErrorTerm[11] = pow(ErrorTerm[11],2);
    ErrorTerm[11] += pow( 2*sqrt(TERMS_4[5]->GetBinContent(bin) + pow(1-TwoFrac,2)*TERMS_4[12]->GetBinContent(bin)) * K4avg[ch1_4][ch2_4][ch3_4][ch4_4][5]->GetBinContent(bin) / TwoFrac, 2);
    ErrorTerm[11] += 4*TERMS_4[12]->GetBinContent(bin);
    ErrorTerm[11] = sqrt(ErrorTerm[11]);

    //
    double FinalValue[4]={0};
    double FinalValue_e[4]={0};
    if(SameCharge) {
      //N4QSvalue = (N4QSvalue - TERMS_4[12]->GetBinContent(bin)) * C4muonCorrectionSC[0]->GetBinContent(bin) + TERMS_4[12]->GetBinContent(bin);
      N4QSvalue *= C4muonCorrectionSC[0]->GetBinContent(bin); 
      ErrorTerm[0] *= C4muonCorrectionSC[0]->GetBinContent(bin);
    }else if(MixedCharge4pionType==1) {
      //N4QSvalue = (N4QSvalue - TERMS_4[12]->GetBinContent(bin)) * C4muonCorrectionMC1[0]->GetBinContent(bin) + TERMS_4[12]->GetBinContent(bin); 
      N4QSvalue *= C4muonCorrectionMC1[0]->GetBinContent(bin);
      ErrorTerm[0] *= C4muonCorrectionMC1[0]->GetBinContent(bin);
    }else {
      //N4QSvalue = (N4QSvalue - TERMS_4[12]->GetBinContent(bin)) * C4muonCorrectionMC2[0]->GetBinContent(bin) + TERMS_4[12]->GetBinContent(bin);
      N4QSvalue *= C4muonCorrectionMC2[0]->GetBinContent(bin); 
      ErrorTerm[0] *= C4muonCorrectionMC2[0]->GetBinContent(bin);
    }
    for(int ii=0; ii<4; ii++) {
      FinalValue[ii] = N4QSvalue; 
      FinalValue_e[ii] = pow(ErrorTerm[0],2);
    }
    
    if(SameCharge) {
      //FinalValue[0] -= 4*(SubtractionTerm[0] - TERMS_4[12]->GetBinContent(bin)) * C4muonCorrectionSC[1]->GetBinContent(bin);
      //FinalValue[0] += 6*(SubtractionTerm[4] - TERMS_4[12]->GetBinContent(bin)) * C4muonCorrectionSC[2]->GetBinContent(bin);
      FinalValue[0] -= 4*(SubtractionTerm[0] * C4muonCorrectionSC[1]->GetBinContent(bin) - TERMS_4[12]->GetBinContent(bin));
      FinalValue[0] += 6*(SubtractionTerm[4] * C4muonCorrectionSC[2]->GetBinContent(bin) - TERMS_4[12]->GetBinContent(bin));
      FinalValue[0] -= 3*(SubtractionTerm[10] - TERMS_4[12]->GetBinContent(bin));

      FinalValue_e[0] += 4*pow(ErrorTerm[1]*C4muonCorrectionSC[1]->GetBinContent(bin),2);
      FinalValue_e[0] += 6*pow(ErrorTerm[5]*C4muonCorrectionSC[2]->GetBinContent(bin),2);
      FinalValue_e[0] += 3*pow(ErrorTerm[11]*C4muonCorrectionSC[3]->GetBinContent(bin),2);
      FinalValue_e[0] += (4*C4muonCorrectionSC[1]->GetBinContent(bin) + 6*C4muonCorrectionSC[2]->GetBinContent(bin) + 3*C4muonCorrectionSC[3]->GetBinContent(bin)) * TERMS_4[12]->GetBinContent(bin);
      //
      FinalValue[1] -= 6*SubtractionTerm[4] - 6*TERMS_4[12]->GetBinContent(bin);// 2-pion removal
      FinalValue_e[1] += 6*pow(ErrorTerm[5],2) + 6*TERMS_4[12]->GetBinContent(bin);
      FinalValue[2] = FinalValue[1] - 3*SubtractionTerm[10] + 3*TERMS_4[12]->GetBinContent(bin);// further 2-pair removal
      FinalValue_e[2] += FinalValue_e[1] + 3*pow(ErrorTerm[11],2) + 3*TERMS_4[12]->GetBinContent(bin);
      FinalValue[3] = SubtractionTerm[10];
    }else{
      if(MixedCharge4pionType==1 && CHARGE==-1){
	//FinalValue[0] -= (SubtractionTerm[0] - TERMS_4[12]->GetBinContent(bin)) * C4muonCorrectionMC1[2]->GetBinContent(bin);// [2] is +++ MuonCorr
	FinalValue[0] -= SubtractionTerm[0] * C4muonCorrectionMC1[2]->GetBinContent(bin) - TERMS_4[12]->GetBinContent(bin);// [2] is +++ MuonCorr
	FinalValue[0] -= 3*(SubtractionTerm[6] * C4muonCorrectionMC1[3]->GetBinContent(bin) - TERMS_4[12]->GetBinContent(bin));
	//
	FinalValue_e[0] += pow(ErrorTerm[1]*C4muonCorrectionMC1[2]->GetBinContent(bin),2);
	FinalValue[1] -= 3*SubtractionTerm[4] - 3*TERMS_4[12]->GetBinContent(bin);
	FinalValue_e[1] += 3*pow(ErrorTerm[5],2) + 3*TERMS_4[12]->GetBinContent(bin);
      }else if(MixedCharge4pionType==1 && CHARGE==+1){
	//FinalValue[0] -= (SubtractionTerm[3] - TERMS_4[12]->GetBinContent(bin)) * C4muonCorrectionMC1[2]->GetBinContent(bin);
	FinalValue[0] -= SubtractionTerm[3] * C4muonCorrectionMC1[2]->GetBinContent(bin) - TERMS_4[12]->GetBinContent(bin);
	FinalValue[0] -= 3*(SubtractionTerm[6] * C4muonCorrectionMC1[3]->GetBinContent(bin) - TERMS_4[12]->GetBinContent(bin));
	//
	FinalValue_e[0] += pow(ErrorTerm[4]*C4muonCorrectionMC1[2]->GetBinContent(bin),2);
	FinalValue[1] -= 3*SubtractionTerm[9] - 3*TERMS_4[12]->GetBinContent(bin);
	FinalValue_e[1] += 3*pow(ErrorTerm[10],2) + 3*TERMS_4[12]->GetBinContent(bin);
      }else{// --++ case
	//FinalValue[0] -= 2*(SubtractionTerm[4] - TERMS_4[12]->GetBinContent(bin)) * C4muonCorrectionMC2[2]->GetBinContent(bin);// old way
	FinalValue[0] -= 2*(SubtractionTerm[4] * C4muonCorrectionMC2[2]->GetBinContent(bin) - TERMS_4[12]->GetBinContent(bin));
	FinalValue[0] -= (SubtractionTerm[10] - TERMS_4[12]->GetBinContent(bin));
	FinalValue[0] -= 4*(SubtractionTerm[6] * C4muonCorrectionMC2[3]->GetBinContent(bin) - TERMS_4[12]->GetBinContent(bin));
	FinalValue_e[0] += 2*pow(ErrorTerm[5],2) + pow(ErrorTerm[11],2) + 2*TERMS_4[12]->GetBinContent(bin);
	//
	FinalValue[1] -= 2*(SubtractionTerm[4] - TERMS_4[12]->GetBinContent(bin)) * C4muonCorrectionMC2[2]->GetBinContent(bin);
	FinalValue_e[1] += 2*pow(ErrorTerm[5],2) + 2*TERMS_4[12]->GetBinContent(bin);
      }
    }
    
    for(int ii=0; ii<4; ii++) FinalValue_e[ii] = sqrt(FinalValue_e[ii]);
    //
    
    n4QS->SetBinContent(bin, FinalValue[0] - TERMS_4[12]->GetBinContent(bin)); 
    n4QS->SetBinError(bin, sqrt( pow(FinalValue_e[0],2) + TERMS_4[12]->GetBinContent(bin)));
    c4QS->SetBinContent(bin, FinalValue[0]);
    c4QS->SetBinError(bin, FinalValue_e[0]);
    C4QS->SetBinContent(bin, N4QSvalue);
    C4QS->SetBinError(bin, ErrorTerm[0]);
    //
    c4QS_RemovalStage1->SetBinContent(bin, FinalValue[1]);
    c4QS_RemovalStage1->SetBinError(bin, FinalValue_e[1]);
    //
    c4QS_RemovalStage2->SetBinContent(bin, FinalValue[2]);
    c4QS_RemovalStage2->SetBinError(bin, FinalValue_e[2]);
    //
    c4QS_RemovalStage3->SetBinContent(bin, FinalValue[3]);
    c4QS_RemovalStage3->SetBinError(bin, FinalValue_e[3]);
        
  }

  C4QS->Divide(TERMS_4[12]);
  c4QS->Divide(TERMS_4[12]);
  c4QS->SetMarkerColor(1); c4QS->SetLineColor(1);
  //
  c4QS_RemovalStage1->Divide(TERMS_4[12]);
  c4QS_RemovalStage1->SetMarkerColor(2); c4QS_RemovalStage1->SetLineColor(2); 
  c4QS_RemovalStage2->Divide(TERMS_4[12]);
  c4QS_RemovalStage2->SetMarkerColor(6); c4QS_RemovalStage2->SetLineColor(6); 
  c4QS_RemovalStage3->Divide(TERMS_4[12]);
  c4QS_RemovalStage3->SetMarkerColor(7); c4QS_RemovalStage3->SetLineColor(7); 
  //
  //
  


  if(SameCharge) C4QS->SetMaximum(9);
  else if(MixedCharge4pionType==1) C4QS->SetMaximum(4);
  else C4QS->SetMaximum(3);
  
  if(CollisionType!=0) C4QS->GetXaxis()->SetRangeUser(0,0.6);

  C4QS->Draw();
  c4QS_RemovalStage1->Draw("same");
  if(SameCharge) c4QS_RemovalStage2->Draw("same");
  //c4QS_RemovalStage3->Draw("same");
  c4QS->Draw("same");

  
  //cout<<TPFullWeight_FourParticle[ch1_4]->GetBinContent(9)<<endl;

  if(SameCharge) {
    //TPFullWeight_FourParticle[ch1_4]->Draw("same");
    FullBuildFromFits->Draw("same");
    PartialBuildFromFits->Draw("same");
  }
  legend3->AddEntry(C4QS,"C_{4}^{QS}","p");
  legend3->AddEntry(c4QS_RemovalStage1,"c_{4}^{QS}{2-pion removal}","p");
  if(SameCharge) legend3->AddEntry(c4QS_RemovalStage2,"c_{4}^{QS}{2-pion+2-pair removal}","p");
  if(SameCharge) legend3->AddEntry(c4QS,"c_{4}^{QS}{2-pion+2-pair+3-pion removal}","p");
  if(!SameCharge && MixedCharge4pionType==1) legend3->AddEntry(c4QS,"c_{4}^{QS}{2-pion+3-pion removal}","p");
  if(!SameCharge && MixedCharge4pionType==2) legend3->AddEntry(c4QS,"c_{4}^{QS}{2-pion+2-pair removal}","p");
  //legend3->AddEntry(c4QS,"c_{4}^{QS}{2-pion+2-pair removal}","p");

  if(SameCharge) legend3->AddEntry(TPFullWeight_FourParticle[ch1_4],"C_{4}^{QS} built","l");
  legend3->Draw("same");

  //K4avg[ch1_4][ch2_4][ch3_4][ch4_4][0]->Draw();

  

  //TH1D *c4QSFitNum = (TH1D*)c4QSFitNum2D->ProjectionY("c4QSFitNum",5,5);
  //TH1D *c4QSFitDen = (TH1D*)c4QSFitDen2D->ProjectionY("c4QSFitDen",5,5);
  //c4QSFitNum->Divide(c4QSFitDen);
  //for(int bin=1; bin<20; bin++){ 
    //double fc4 = 
    //c4QSFitNum->SetBinContent(bin, (c4QSFitNum->GetBinContent(bin) - 1.0)*
  //}
  //c4QSFitNum->Draw("same");

  //C4QS_Syst->Draw("E2 same");
  //C4QS->Draw("same");
  //C4raw->Draw();
  //TH1D *testTerm = (TH1D*)TERMS_4[11]->Clone();
  //testTerm->Divide(TERMS_4[12]);
  //testTerm->SetMinimum(1); testTerm->SetMaximum(1.4); testTerm->GetXaxis()->SetRangeUser(0,0.14);
  //testTerm->Draw();

  ////////////////////////////////////////////////////////////////
  if(SameCharge && CollisionType==0){
    TCanvas *can4 = new TCanvas("can4", "can4",600,600,700,500);
    can4->SetHighLightColor(2);
    can4->Range(-0.7838115,-1.033258,0.7838115,1.033258);
    gStyle->SetOptFit(0111);
    can4->SetFillColor(10);
    can4->SetBorderMode(0);
    can4->SetBorderSize(2);
    can4->SetGridx();
    can4->SetGridy();
    can4->SetFrameFillColor(0);
    can4->SetFrameBorderMode(0);
    can4->SetFrameBorderMode(0);
    gPad->SetRightMargin(0.02); gPad->SetTopMargin(0.02);
    TLegend *legend3_2 = (TLegend*)legend1->Clone();
    legend3_2->SetX1(0.45); legend3_2->SetX2(0.98);  legend3_2->SetY1(0.6); legend3_2->SetY2(0.95);
    //
    //TH1D *Chi2NDF_PointSize_4 = new TH1D("Chi2NDF_PointSize_4","",40,-0.5,39.5);
    //TH1D *Chi2NDF_FullSize_4 = new TH1D("Chi2NDF_FullSize_4","",40,-0.5,39.5);
    TH1D *Chi2NDF_PointSize_4 = new TH1D("Chi2NDF_PointSize_4","",100,-0.5,99.5);
    TH1D *Chi2NDF_FullSize_4 = new TH1D("Chi2NDF_FullSize_4","",100,-0.5,99.5);
    Chi2NDF_PointSize_4->SetLineColor(4); Chi2NDF_FullSize_4->SetLineColor(2);
    Chi2NDF_PointSize_4->SetMarkerColor(4); Chi2NDF_FullSize_4->SetMarkerColor(2);
    Chi2NDF_PointSize_4->GetXaxis()->SetTitle("Coherent fraction (%)"); Chi2NDF_PointSize_4->GetYaxis()->SetTitle("#sqrt{#chi^{2}}");
    
    
    TH1D *tempDen = (TH1D*)TPFullWeight_FourParticle_2D[ch1_4]->ProjectionY("TPFullWeight4_Den", 4, 4);
    TH1D *tempDenNeg = (TH1D*)TPNegFullWeight_FourParticle_2D[ch1_4]->ProjectionY("TPNegFullWeight4_Den", 4, 4);
    tempDen->Add(tempDenNeg);// Add Pos and Neg Weight

    for(int binG=5; binG<=104; binG++){// 44
      TString *proName=new TString("TPFullWeight4_");
      *proName += binG;
      TH1D *tempNum = (TH1D*)TPFullWeight_FourParticle_2D[ch1_4]->ProjectionY(proName->Data(), binG, binG);
      proName->Append("_Neg");
      TH1D *tempNumNeg = (TH1D*)TPNegFullWeight_FourParticle_2D[ch1_4]->ProjectionY(proName->Data(), binG, binG);
      //
      // Add Pos and Neg Weights
      tempNum->Add(tempNumNeg);
      //
      //tempNum->Add(tempDen);// now added in InterpCorr section
      tempNum->Divide(tempDen);
      
      
      double Chi2=0;
      double NDF=0;
      for(int binQ4=4; binQ4<=4; binQ4++){
	if(tempNum->GetBinContent(binQ4) <=0) continue;
	double value = C4QS->GetBinContent(binQ4) - tempNum->GetBinContent(binQ4);
	double err = pow(C4QS->GetBinError(binQ4),2);

	err = sqrt(err);
	if(err<=0) continue;
		
	Chi2 += pow(value / err,2);
	//
	NDF += 1;
      }
      //if(binG<25) {Chi2NDF_PointSize_4->SetBinContent(1 + 2*(binG-5), sqrt(fabs(Chi2))/NDF); Chi2NDF_PointSize_4->SetBinError(1 + 2*(binG-5), 0.001);}
      //else {Chi2NDF_FullSize_4->SetBinContent(1 + 2*(binG-25), sqrt(fabs(Chi2))/NDF); Chi2NDF_FullSize_4->SetBinError(1 + 2*(binG-25), 0.001);}
      if(binG<55) {Chi2NDF_PointSize_4->SetBinContent(1 + 2*(binG-5), sqrt(fabs(Chi2))/NDF); Chi2NDF_PointSize_4->SetBinError(1 + 2*(binG-5), 0.001);}
      else {Chi2NDF_FullSize_4->SetBinContent(1 + 2*(binG-55), sqrt(fabs(Chi2))/NDF); Chi2NDF_FullSize_4->SetBinError(1 + 2*(binG-55), 0.001);}
    }
    Chi2NDF_PointSize_4->SetMarkerStyle(20); Chi2NDF_FullSize_4->SetMarkerStyle(20);
    

    Chi2NDF_PointSize_4->Draw();
    Chi2NDF_FullSize_4->Draw("same");
    legend3_2->AddEntry(Chi2NDF_PointSize_4,"R_{coh}=0 (Point Source)","p");
    legend3_2->AddEntry(Chi2NDF_FullSize_4,"R_{coh}=R_{ch}","p");
    legend3_2->Draw("same");


    
  }

  // Print 4-pion points
  for(int bin=1; bin<=12; bin++){
    //cout<<C4QS->GetBinContent(bin)<<", ";
    //cout<<c4QS->GetBinContent(bin)<<", ";
    //cout<<TPFullWeight_FourParticle[ch1_4]->GetBinContent(bin)<<", ";
    //cout<<C4raw->GetBinContent(bin)<<", ";
    //cout<<K4avg[ch1_4][ch2_4][ch3_4][ch4_4][0]->GetBinContent(bin)<<", ";
  }
  cout<<endl;
  for(int bin=1; bin<=12; bin++){
    //cout<<c4QS->GetBinContent(bin)<<", ";
    //cout<<C4QS->GetBinError(bin)<<", ";
    //cout<<c4QS->GetBinError(bin)<<", ";
    //cout<<C4raw->GetBinError(bin)<<", ";
  }
  //cout<<endl;
  ////////////////////////////////////////////////////////////////
  // r4
  
  TF1 *ChaoticLimit_r42 = new TF1("ChaoticLimit_r42","6",0,10);
  ChaoticLimit_r42->SetLineStyle(2);
  TH1D *r42;
  if(SameCharge && CollisionType==0){
    /*TCanvas *can5 = new TCanvas("can5", "can5",1200,600,700,500);
    can5->SetHighLightColor(2);
    can5->Range(-0.7838115,-1.033258,0.7838115,1.033258);
    gStyle->SetOptFit(0111);
    can5->SetFillColor(10);
    can5->SetBorderMode(0);
    can5->SetBorderSize(2);
    can5->SetGridx();
    can5->SetGridy();
    can5->SetFrameFillColor(0);
    can5->SetFrameBorderMode(0);
    can5->SetFrameBorderMode(0);
    gPad->SetRightMargin(0.02); gPad->SetTopMargin(0.02);*/

    r42 = (TH1D*)n4QS->Clone();
    TPN_FourParticle[ch1_4]->Multiply(MRC_SC_4[4]);
    r42->Divide(TPN_FourParticle[ch1_4]);
    r42->GetXaxis()->SetRangeUser(0,0.08);
    r42->SetMinimum(0.5); r42->SetMaximum(14);
    r42->GetYaxis()->SetTitle("r_{4}{2}");
    //
    //r42->Draw();
    //ChaoticLimit_r42->Draw("same");
    //TPN_FourParticle[ch1_4]->Draw();
    //fTPNRejects4pion1->SetLineColor(2);
    //fTPNRejects4pion1->Draw("same");
  }
  
  //double EA = exp(-pow(0.005 * 10/FmToGeV,2)/2.);
  //TF1 *C2mag = new TF1("C2mag","pow(1-x,2)*pow([0],2) + 2*x*(1-x)*[0]",0,1);
  //C2mag->FixParameter(0, EA);
  //TF1 *C4norm = new TF1("C4norm","( 6*pow(1-x,4)*pow([0],4) + 24*x*pow(1-x,3)*pow([0],3) + 4*(2*pow(1-x,3)*pow([0],3) + 6*x*pow(1-x,2)*pow([0],2)) + 3*pow(C2mag,2) + 6*C2mag + 1) / (6*pow(C2mag,2) + 8*pow(C2mag,1.5) + 3*pow(C2mag,2) + 6*C2mag + 1)",0,1);
  //C4norm->FixParameter(0, EA);
  //C4norm->Draw();

  /*TF1 *C2mag = new TF1("C2mag","pow(1-[0],2)*exp(-pow(x*[1],2)/6.) + 2*[0]*(1-[0])*exp(-pow(x*[1],2)/12.)",0,0.2);
  C2mag->FixParameter(0, 0.25);
  C2mag->FixParameter(1, 10./FmToGeV);
  TF1 *C4norm = new TF1("C4norm","( 6*pow(1-[0],4)*exp(-pow(x*[1],2)/3.) + 24*[0]*pow(1-[0],3)*exp(-pow(x*[1],2)/4.) + 4*(2*pow(1-[0],3)*exp(-pow(x*[1],2)/4.) + 6*[0]*pow(1-[0],2)*exp(-pow(x*[1],2)/6.)) + 3*pow(C2mag,2) + 6*C2mag + 1) / (6*pow(C2mag,2) + 8*pow(C2mag,1.5) + 3*pow(C2mag,2) + 6*C2mag + 1)",0,0.2);
  C4norm->FixParameter(0, 0.25);
  C4norm->FixParameter(1, 10./FmToGeV);*/
  //C4norm->Draw();


  if(SaveToFile){
    TString *savefilename = new TString("OutFiles/OutFile");
    if(MCcase) savefilename->Append("MonteCarlo");
    //
    if(SameCharge) savefilename->Append("SC");
    else savefilename->Append("MC");
    //
    if(!SameCharge) *savefilename += MixedCharge4pionType_def;
    //
    if(CHARGE==1) savefilename->Append("_Pos_");
    else savefilename->Append("_Neg_");
    //
    
    savefilename->Append("KT_");
    *savefilename += EDbin+1;
    
    savefilename->Append("_M");
    *savefilename += Mbin;
    savefilename->Append(".root");
    TFile *savefile = new TFile(savefilename->Data(),"RECREATE");
    //
    C2->Write("C2");
    C2QS->Write("C2QS");
    C3QS->Write("C3QS");
    c3QS->Write("c3QS");
    C4QS->Write("C4QS");
    c4QS->Write("c4QS");
    c4QS_RemovalStage1->Write("c4QS_RemovalStage1");
    if(SameCharge) {
      r3->Write("r3");
      r42->Write("r42");
      c4QS_RemovalStage2->Write("c4QS_RemovalStage2");
      TPFullWeight_ThreeParticle[ch1_3]->Write("C3QS_built");
      TPFullWeight_FourParticle[ch1_4]->Write("C4QS_built");
      //
      TPFullWeight_ThreeParticle_2D[ch1_3]->Write("C3QS_built2D");
      TPFullWeight_FourParticle_2D[ch1_4]->Write("C4QS_built2D");
      TPNegFullWeight_ThreeParticle_2D[ch1_3]->Write("C3QS_Negbuilt2D");
      TPNegFullWeight_FourParticle_2D[ch1_4]->Write("C4QS_Negbuilt2D");
    }
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
  if(FileSetting==7) momresfilename->Append("_10percentIncrease");
  if(CollisionType!=0) momresfilename->Append("_ppAndpPb");
  momresfilename->Append(".root");
  
  TFile *MomResFile = new TFile(momresfilename->Data(),"READ");
  TString *proName[28];
  for(int ii=0; ii<28; ii++){
    proName[ii] = new TString("MRC_pro_");
    *proName[ii] += ii;
  }

  MRC_SC_2[0] = (TH1D*)(((TH2D*)MomResFile->Get("MRC_2_SC_term1"))->ProjectionY(proName[0]->Data(), RbinMRC, RbinMRC));
  MRC_SC_2[1] = (TH1D*)(((TH2D*)MomResFile->Get("MRC_2_SC_term2"))->ProjectionY(proName[1]->Data(), RbinMRC, RbinMRC));
  MRC_SC_2[0]->SetDirectory(0); MRC_SC_2[1]->SetDirectory(0); 
  //
  MRC_MC_2[0] = (TH1D*)(((TH2D*)MomResFile->Get("MRC_2_MC_term1"))->ProjectionY(proName[2]->Data(), RbinMRC, RbinMRC));
  MRC_MC_2[1] = (TH1D*)(((TH2D*)MomResFile->Get("MRC_2_MC_term2"))->ProjectionY(proName[3]->Data(), RbinMRC, RbinMRC));
  MRC_MC_2[0]->SetDirectory(0); MRC_MC_2[1]->SetDirectory(0); 
  //
  MRC_SC_3[0] = (TH1D*)(((TH2D*)MomResFile->Get("MRC_3_SC_term1"))->ProjectionY(proName[4]->Data(), RbinMRC, RbinMRC));
  MRC_SC_3[1] = (TH1D*)(((TH2D*)MomResFile->Get("MRC_3_SC_term2"))->ProjectionY(proName[5]->Data(), RbinMRC, RbinMRC));
  MRC_SC_3[2] = (TH1D*)(((TH2D*)MomResFile->Get("MRC_3_SC_term3"))->ProjectionY(proName[6]->Data(), RbinMRC, RbinMRC));
  MRC_SC_3[0]->SetDirectory(0); MRC_SC_3[1]->SetDirectory(0); MRC_SC_3[2]->SetDirectory(0);
  //
  MRC_MC_3[0] = (TH1D*)(((TH2D*)MomResFile->Get("MRC_3_MC_term1"))->ProjectionY(proName[7]->Data(), RbinMRC, RbinMRC));
  MRC_MC_3[1] = (TH1D*)(((TH2D*)MomResFile->Get("MRC_3_MC_term2"))->ProjectionY(proName[8]->Data(), RbinMRC, RbinMRC));
  MRC_MC_3[2] = (TH1D*)(((TH2D*)MomResFile->Get("MRC_3_MC_term3"))->ProjectionY(proName[9]->Data(), RbinMRC, RbinMRC));
  MRC_MC_3[3] = (TH1D*)(((TH2D*)MomResFile->Get("MRC_3_MC_term4"))->ProjectionY(proName[10]->Data(), RbinMRC, RbinMRC));
  MRC_MC_3[0]->SetDirectory(0); MRC_MC_3[1]->SetDirectory(0); MRC_MC_3[2]->SetDirectory(0); MRC_MC_3[3]->SetDirectory(0);
  //
  MRC_SC_4[0] = (TH1D*)(((TH2D*)MomResFile->Get("MRC_4_SC_term1"))->ProjectionY(proName[11]->Data(), RbinMRC, RbinMRC));
  MRC_SC_4[1] = (TH1D*)(((TH2D*)MomResFile->Get("MRC_4_SC_term2"))->ProjectionY(proName[12]->Data(), RbinMRC, RbinMRC));
  MRC_SC_4[2] = (TH1D*)(((TH2D*)MomResFile->Get("MRC_4_SC_term3"))->ProjectionY(proName[13]->Data(), RbinMRC, RbinMRC));
  MRC_SC_4[3] = (TH1D*)(((TH2D*)MomResFile->Get("MRC_4_SC_term4"))->ProjectionY(proName[14]->Data(), RbinMRC, RbinMRC));
  MRC_SC_4[4] = (TH1D*)(((TH2D*)MomResFile->Get("MRC_4_SC_term5"))->ProjectionY(proName[15]->Data(), RbinMRC, RbinMRC));
  MRC_SC_4[0]->SetDirectory(0); MRC_SC_4[1]->SetDirectory(0); MRC_SC_4[2]->SetDirectory(0); MRC_SC_4[3]->SetDirectory(0);
  MRC_SC_4[4]->SetDirectory(0);
  //
  MRC_MC1_4[0] = (TH1D*)(((TH2D*)MomResFile->Get("MRC_4_MC1_term1"))->ProjectionY(proName[16]->Data(), RbinMRC, RbinMRC));
  MRC_MC1_4[1] = (TH1D*)(((TH2D*)MomResFile->Get("MRC_4_MC1_term2"))->ProjectionY(proName[17]->Data(), RbinMRC, RbinMRC));
  MRC_MC1_4[2] = (TH1D*)(((TH2D*)MomResFile->Get("MRC_4_MC1_term3"))->ProjectionY(proName[18]->Data(), RbinMRC, RbinMRC));
  MRC_MC1_4[3] = (TH1D*)(((TH2D*)MomResFile->Get("MRC_4_MC1_term4"))->ProjectionY(proName[19]->Data(), RbinMRC, RbinMRC));
  MRC_MC1_4[4] = (TH1D*)(((TH2D*)MomResFile->Get("MRC_4_MC1_term5"))->ProjectionY(proName[20]->Data(), RbinMRC, RbinMRC));
  MRC_MC1_4[5] = (TH1D*)(((TH2D*)MomResFile->Get("MRC_4_MC1_term6"))->ProjectionY(proName[21]->Data(), RbinMRC, RbinMRC));
  MRC_MC1_4[0]->SetDirectory(0); MRC_MC1_4[1]->SetDirectory(0); MRC_MC1_4[2]->SetDirectory(0); MRC_MC1_4[3]->SetDirectory(0);
  MRC_MC1_4[4]->SetDirectory(0); MRC_MC1_4[5]->SetDirectory(0);
  //
  MRC_MC2_4[0] = (TH1D*)(((TH2D*)MomResFile->Get("MRC_4_MC2_term1"))->ProjectionY(proName[22]->Data(), RbinMRC, RbinMRC));
  MRC_MC2_4[1] = (TH1D*)(((TH2D*)MomResFile->Get("MRC_4_MC2_term2"))->ProjectionY(proName[23]->Data(), RbinMRC, RbinMRC));
  MRC_MC2_4[2] = (TH1D*)(((TH2D*)MomResFile->Get("MRC_4_MC2_term3"))->ProjectionY(proName[24]->Data(), RbinMRC, RbinMRC));
  MRC_MC2_4[3] = (TH1D*)(((TH2D*)MomResFile->Get("MRC_4_MC2_term4"))->ProjectionY(proName[25]->Data(), RbinMRC, RbinMRC));
  MRC_MC2_4[4] = (TH1D*)(((TH2D*)MomResFile->Get("MRC_4_MC2_term5"))->ProjectionY(proName[26]->Data(), RbinMRC, RbinMRC));
  MRC_MC2_4[5] = (TH1D*)(((TH2D*)MomResFile->Get("MRC_4_MC2_term6"))->ProjectionY(proName[27]->Data(), RbinMRC, RbinMRC));
  MRC_MC2_4[0]->SetDirectory(0); MRC_MC2_4[1]->SetDirectory(0); MRC_MC2_4[2]->SetDirectory(0); MRC_MC2_4[3]->SetDirectory(0);
  MRC_MC2_4[4]->SetDirectory(0); MRC_MC2_4[5]->SetDirectory(0);

  if(!MRC || MCcase_def){
    for(int bin=1; bin<=MRC_SC_2[0]->GetNbinsX(); bin++){
      MRC_SC_2[0]->SetBinContent(bin, 1.0); MRC_SC_2[1]->SetBinContent(bin, 1.0);
      MRC_MC_2[0]->SetBinContent(bin, 1.0); MRC_MC_2[1]->SetBinContent(bin, 1.0);
    }
    for(int bin=1; bin<=MRC_SC_3[0]->GetNbinsX(); bin++){
      MRC_SC_3[0]->SetBinContent(bin, 1.0); MRC_SC_3[1]->SetBinContent(bin, 1.0); MRC_SC_3[2]->SetBinContent(bin, 1.0);
      MRC_MC_3[0]->SetBinContent(bin, 1.0); MRC_MC_3[1]->SetBinContent(bin, 1.0); MRC_MC_3[2]->SetBinContent(bin, 1.0);
      MRC_MC_3[3]->SetBinContent(bin, 1.0);
    }
    for(int bin=1; bin<=MRC_SC_4[0]->GetNbinsX(); bin++){
      MRC_SC_4[0]->SetBinContent(bin, 1.0); MRC_SC_4[1]->SetBinContent(bin, 1.0); MRC_SC_4[2]->SetBinContent(bin, 1.0);
      MRC_SC_4[3]->SetBinContent(bin, 1.0); MRC_SC_4[4]->SetBinContent(bin, 1.0);
      MRC_MC1_4[0]->SetBinContent(bin, 1.0); MRC_MC1_4[1]->SetBinContent(bin, 1.0); MRC_MC1_4[2]->SetBinContent(bin, 1.0);
      MRC_MC1_4[3]->SetBinContent(bin, 1.0); MRC_MC1_4[4]->SetBinContent(bin, 1.0); MRC_MC1_4[5]->SetBinContent(bin, 1.0);
      //
      MRC_MC2_4[0]->SetBinContent(bin, 1.0); MRC_MC2_4[1]->SetBinContent(bin, 1.0); MRC_MC2_4[2]->SetBinContent(bin, 1.0);
      MRC_MC2_4[3]->SetBinContent(bin, 1.0); MRC_MC2_4[4]->SetBinContent(bin, 1.0); MRC_MC2_4[5]->SetBinContent(bin, 1.0);
    }
  }
  MomResFile->Close();
  
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
//________________________________________________________________________________________
void SetMuonCorrections(){
  TString *name = new TString("MuonCorrection");
  if(FileSetting==8) name->Append("_92percent");
  if(CollisionType!=0) name->Append("_ppAndpPb");
  name->Append(".root");
  TFile *MuonFile=new TFile(name->Data(),"READ");
  TString *proName[22];
  for(int ii=0; ii<22; ii++){
    proName[ii] = new TString("MuonCorr_pro_");
    *proName[ii] += ii;
  }
  C2muonCorrectionSC = (TH1D*)(((TH2D*)MuonFile->Get("C2muonCorrection_SC"))->ProjectionY(proName[0]->Data(), RbinMRC, RbinMRC));
  C2muonCorrectionMC = (TH1D*)(((TH2D*)MuonFile->Get("C2muonCorrection_MC"))->ProjectionY(proName[1]->Data(), RbinMRC, RbinMRC));
  WeightmuonCorrection = (TH1D*)(((TH2D*)MuonFile->Get("WeightmuonCorrection"))->ProjectionY(proName[2]->Data(), RbinMRC, RbinMRC));
  C2muonCorrectionSC->SetDirectory(0); C2muonCorrectionMC->SetDirectory(0); WeightmuonCorrection->SetDirectory(0);
  //
  C3muonCorrectionSC[0] = (TH1D*)(((TH2D*)MuonFile->Get("C3muonCorrection_SC_term1"))->ProjectionY(proName[3]->Data(), RbinMRC, RbinMRC));
  C3muonCorrectionSC[1] = (TH1D*)(((TH2D*)MuonFile->Get("C3muonCorrection_SC_term2"))->ProjectionY(proName[4]->Data(), RbinMRC, RbinMRC));
  C3muonCorrectionMC[0] = (TH1D*)(((TH2D*)MuonFile->Get("C3muonCorrection_MC_term1"))->ProjectionY(proName[5]->Data(), RbinMRC, RbinMRC));
  C3muonCorrectionMC[1] = (TH1D*)(((TH2D*)MuonFile->Get("C3muonCorrection_MC_term2"))->ProjectionY(proName[6]->Data(), RbinMRC, RbinMRC));
  C3muonCorrectionMC[2] = (TH1D*)(((TH2D*)MuonFile->Get("C3muonCorrection_MC_term3"))->ProjectionY(proName[7]->Data(), RbinMRC, RbinMRC));
  C3muonCorrectionSC[0]->SetDirectory(0); C3muonCorrectionSC[1]->SetDirectory(0);
  C3muonCorrectionMC[0]->SetDirectory(0); C3muonCorrectionMC[1]->SetDirectory(0); C3muonCorrectionMC[2]->SetDirectory(0);
  //
  C4muonCorrectionSC[0] = (TH1D*)(((TH2D*)MuonFile->Get("C4muonCorrection_SC_term1"))->ProjectionY(proName[8]->Data(), RbinMRC, RbinMRC));
  C4muonCorrectionSC[1] = (TH1D*)(((TH2D*)MuonFile->Get("C4muonCorrection_SC_term2"))->ProjectionY(proName[9]->Data(), RbinMRC, RbinMRC));
  C4muonCorrectionSC[2] = (TH1D*)(((TH2D*)MuonFile->Get("C4muonCorrection_SC_term3"))->ProjectionY(proName[10]->Data(), RbinMRC, RbinMRC));
  C4muonCorrectionSC[3] = (TH1D*)(((TH2D*)MuonFile->Get("C4muonCorrection_SC_term4"))->ProjectionY(proName[11]->Data(), RbinMRC, RbinMRC));
  C4muonCorrectionSC[0]->SetDirectory(0); C4muonCorrectionSC[1]->SetDirectory(0);
  C4muonCorrectionSC[2]->SetDirectory(0); C4muonCorrectionSC[3]->SetDirectory(0);
  //
  C4muonCorrectionMC1[0] = (TH1D*)(((TH2D*)MuonFile->Get("C4muonCorrection_MC1_term1"))->ProjectionY(proName[12]->Data(), RbinMRC, RbinMRC));
  C4muonCorrectionMC1[1] = (TH1D*)(((TH2D*)MuonFile->Get("C4muonCorrection_MC1_term2"))->ProjectionY(proName[13]->Data(), RbinMRC, RbinMRC));
  C4muonCorrectionMC1[2] = (TH1D*)(((TH2D*)MuonFile->Get("C4muonCorrection_MC1_term3"))->ProjectionY(proName[14]->Data(), RbinMRC, RbinMRC));
  C4muonCorrectionMC1[3] = (TH1D*)(((TH2D*)MuonFile->Get("C4muonCorrection_MC1_term4"))->ProjectionY(proName[15]->Data(), RbinMRC, RbinMRC));
  C4muonCorrectionMC1[4] = (TH1D*)(((TH2D*)MuonFile->Get("C4muonCorrection_MC1_term5"))->ProjectionY(proName[16]->Data(), RbinMRC, RbinMRC));
  C4muonCorrectionMC1[0]->SetDirectory(0); C4muonCorrectionMC1[1]->SetDirectory(0);
  C4muonCorrectionMC1[2]->SetDirectory(0); C4muonCorrectionMC1[3]->SetDirectory(0); C4muonCorrectionMC1[4]->SetDirectory(0);
  //
  C4muonCorrectionMC2[0] = (TH1D*)(((TH2D*)MuonFile->Get("C4muonCorrection_MC2_term1"))->ProjectionY(proName[17]->Data(), RbinMRC, RbinMRC));
  C4muonCorrectionMC2[1] = (TH1D*)(((TH2D*)MuonFile->Get("C4muonCorrection_MC2_term2"))->ProjectionY(proName[18]->Data(), RbinMRC, RbinMRC));
  C4muonCorrectionMC2[2] = (TH1D*)(((TH2D*)MuonFile->Get("C4muonCorrection_MC2_term3"))->ProjectionY(proName[19]->Data(), RbinMRC, RbinMRC));
  C4muonCorrectionMC2[3] = (TH1D*)(((TH2D*)MuonFile->Get("C4muonCorrection_MC2_term4"))->ProjectionY(proName[20]->Data(), RbinMRC, RbinMRC));
  C4muonCorrectionMC2[4] = (TH1D*)C4muonCorrectionSC[3]->Clone();
  C4muonCorrectionMC2[0]->SetDirectory(0); C4muonCorrectionMC2[1]->SetDirectory(0);
  C4muonCorrectionMC2[2]->SetDirectory(0); C4muonCorrectionMC2[3]->SetDirectory(0); C4muonCorrectionMC2[4]->SetDirectory(0);
  //
  //
  if(!MuonCorrection || MCcase_def){
    for(int bin=1; bin<=C2muonCorrectionSC->GetNbinsX(); bin++){
      C2muonCorrectionSC->SetBinContent(bin, 1.0);
      C2muonCorrectionMC->SetBinContent(bin, 1.0);
      WeightmuonCorrection->SetBinContent(bin, 1.0);
    }
    for(int bin=1; bin<=C3muonCorrectionSC[0]->GetNbinsX(); bin++){
      C3muonCorrectionSC[0]->SetBinContent(bin, 1.0); C3muonCorrectionSC[1]->SetBinContent(bin, 1.0);
      C3muonCorrectionMC[0]->SetBinContent(bin, 1.0); C3muonCorrectionMC[1]->SetBinContent(bin, 1.0); C3muonCorrectionMC[2]->SetBinContent(bin, 1.0); 
    }
    for(int bin=1; bin<=C4muonCorrectionSC[0]->GetNbinsX(); bin++){
      C4muonCorrectionSC[0]->SetBinContent(bin, 1.0); C4muonCorrectionSC[1]->SetBinContent(bin, 1.0); C4muonCorrectionSC[2]->SetBinContent(bin, 1.0);
      C4muonCorrectionSC[3]->SetBinContent(bin, 1.0);
      C4muonCorrectionMC1[0]->SetBinContent(bin, 1.0); C4muonCorrectionMC1[1]->SetBinContent(bin, 1.0); C4muonCorrectionMC1[2]->SetBinContent(bin, 1.0); 
      C4muonCorrectionMC1[3]->SetBinContent(bin, 1.0); C4muonCorrectionMC1[4]->SetBinContent(bin, 1.0);
      //
      C4muonCorrectionMC2[0]->SetBinContent(bin, 1.0); C4muonCorrectionMC2[1]->SetBinContent(bin, 1.0); C4muonCorrectionMC2[2]->SetBinContent(bin, 1.0); 
      C4muonCorrectionMC2[3]->SetBinContent(bin, 1.0); C4muonCorrectionMC2[4]->SetBinContent(bin, 1.0);
    }
  }
  
  MuonFile->Close();
}
//________________________________________________________________________
void SetFSIindex(Float_t R){
  if(!MCcase_def){
    if(CollisionType==0){
      if(Mbin_def==0) fFSIindex = 0;//0-5%
      else if(Mbin_def==1) fFSIindex = 1;//5-10%
      else if(Mbin_def<=3) fFSIindex = 2;//10-20%
      else if(Mbin_def<=5) fFSIindex = 3;//20-30%
      else if(Mbin_def<=7) fFSIindex = 4;//30-40%
      else if(Mbin_def<=9) fFSIindex = 5;//40-50%
      else if(Mbin_def<=12) fFSIindex = 6;//40-50%
      else if(Mbin_def<=15) fFSIindex = 7;//40-50%
      else if(Mbin_def<=18) fFSIindex = 8;//40-50%
      else fFSIindex = 8;//90-100%
    }else fFSIindex = 9;// pp and pPb
  }else{// FSI binning for MC 
    if(R>=10.) fFSIindex = 0;
    else if(R>=9.) fFSIindex = 1;
    else if(R>=8.) fFSIindex = 2;
    else if(R>=7.) fFSIindex = 3;
    else if(R>=6.) fFSIindex = 4;
    else if(R>=5.) fFSIindex = 5;
    else if(R>=4.) fFSIindex = 6;
    else if(R>=3.) fFSIindex = 7;
    else if(R>=2.) fFSIindex = 8;
    else fFSIindex = 9;
  }
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
