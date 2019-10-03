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

float TwoFrac=0.7;
//
bool MCcase_def=0;// MC data?
int CollisionType_def=1;// PbPb, pPb, pp
//
int Mbin_def=0;// 0-9: centrality bin in widths of 5%
bool SameCharge_def=0;// same-charge?
int CHARGE_def=-1;// -1 or +1: + or - pions for same-charge case, --+ or -++,  ---+ or -+++
int MixedCharge4pionType_def = 2;// 1(---+) or 2(--++)
//
int EDbin_def=0;// 0 or 1: Kt3 bin
int TPNbin=3;// TPN bin for r3 and r4
int Gbin = int( (0) /2. ) + 155;// + 5 + 25*RcohIndex: 5, 30, 55, 80, 105, 130, 155 (t=T).  215 for 35Gbin File
int Ktbin_def=1;// 1(0.2-0.3),..., 6(0.7-0.8), 10 = Full Range
bool FitBuild=0;
//
bool ReNormBuiltBaseline=1;// Re-normalize Built Baseline
bool MRC=1;// Momentum Resolution Corrections?
bool MuonCorrection=1;// correct for Muons?
bool FSICorrection=1;// correct For Final-State-Interactions
//
int f_choice=0;// 0(Core/Halo), 1(60fm), 2(80fm), 3(100fm)
//
//
bool SaveToFile_def=0;// Save outputs to file?
bool GeneratedSignal=kFALSE;// Was the QS+FSI signal generated in MC? 
//
//
//
//

int fFSIindex=0;
int Ktlowbin;
int Kthighbin;
float OneFrac;// Lambda^{1/2}
float ThreeFrac;// Lambda^{3/2}
float FourFrac;// lambda^{4/2}
double Qcut_pp = 0.6;// 0.6
double Qcut_PbPb = 0.1;
double NormQcutLow_pp = 0.27;// 0.9, 0.45
double NormQcutHigh_pp = 0.33;// 1.2, 0.55
double NormQcutLow_PbPb = .15;// .15
double NormQcutHigh_PbPb = .2;// .2
double ReNormL_3;
double ReNormH_3;
double ReNormL_4;
double ReNormH_4;

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


void Plot_FourPion(bool SaveToFile=SaveToFile_def, bool MCcase=MCcase_def, bool SameCharge=SameCharge_def, int CollisionType=CollisionType_def, int MixedCharge4pionType=MixedCharge4pionType_def, int EDbin=EDbin_def, int CHARGE=CHARGE_def, int Mbin=Mbin_def, int Ktbin=Ktbin_def){
  
  EDbin_def=EDbin;
  SaveToFile_def=SaveToFile;
  MCcase_def=MCcase;
  CHARGE_def=CHARGE;
  SameCharge_def=SameCharge;// 3-pion same-charge
  CollisionType_def=CollisionType;
  MixedCharge4pionType_def=MixedCharge4pionType;
  Mbin_def=Mbin;
  Ktbin_def=Ktbin;
  //
  Ktlowbin=(Ktbin)*2+3;// kt bins are 0.5 GeV/c wide (0-0.05, 0.05-0.1 ...)
  Kthighbin=(Ktbin)*2+4;
  //
  //Ktlowbin=5;// 1, 5,...
  //Kthighbin=5;// 4, 5,... 
  //
  //Ktlowbin=(Ktbin)*10+15;// kt bins are 0.5 GeV/c wide (0-0.05, 0.05-0.1 ...)
  //Kthighbin=(Ktbin)*10+20;
  //

 

  OneFrac = sqrt(TwoFrac);
  ThreeFrac = pow(TwoFrac, 1.5);
  FourFrac = pow(TwoFrac, 2.);

  int Ysize=600;
  if(SameCharge==1) Ysize=800;
  double Ystart=0.3;
  if(SameCharge==0) Ystart=0.005;
  double Xmarge=0.01;
  if(SameCharge==0) Xmarge=0.14;
  
  //// Core/Halo, 60fm, 80fm, 100fm
  float ThermShift_f33[4]={0., 0.05884, 0.04487, 0.04132};
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

  //// Core/Halo, 60fm, 80fm, 100fm
  float ThermShift_f44[4]={0., 0.05816, 0.03786, 0.03221};
  float ThermShift_f43[4]={0., -0.009794, -0.006847, -0.005939};
  float ThermShift_f42[4]={0., -0.002967, -0.00167 , -0.001371};
  float ThermShift_f41[4]={0., -0.001182, -0.0004481, -0.0002327};
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
  //
  /*float shift4=0.08;
  f44 -= shift4;
  f43 += 1/3.*shift4/4.;
  f42 += 1/3.*shift4/6.;
  f41 += 1/3.*shift4;*/
  //
  /*float shift3=0.04;
  f33 -= shift3;
  f32 += 1/2.*shift3/3.;
  f31 += 1/2.*shift3;*/
  


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

  cout<<"CollisionType = "<<CollisionType<<"   Mbin = "<<Mbin<<"   EDbin = "<<EDbin<<"   SameCharge = "<<SameCharge<<endl;
  if(!SameCharge) cout<<"4-pion MixedCharge type = "<<MixedCharge4pionType<<endl;
  cout<<"TwoFrac = "<<TwoFrac<<endl;
  cout<<"Gbin = "<<Gbin<<endl;

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

  if(CollisionType==0){
    ReNormL_3=0.095 + 0.03*(Mbin/9.);//.095 (default), 0.075
    ReNormH_3=0.105 + 0.03*(Mbin/9.);//.105 (default), 0.085
    ReNormL_4=0.125 + 0.04*(Mbin/9.);//.14 (default), 0.13
    ReNormH_4=0.145 + 0.04*(Mbin/9.);//.16 (default), 0.14
    // older values
    //ReNormL_3=0.095;
    //ReNormH_3=0.105;
    //ReNormL_4=0.16;
    //ReNormH_4=0.17;
  }else{
    ReNormL_3=0.3;
    ReNormH_3=0.32;
    ReNormL_4=0.46;//.46
    ReNormH_4=0.49;//.49
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
 
  TF1 *C3ratioFit = new TF1("C3ratioFit","[0] + [1]*exp(-pow(x*[2],2))",0,1);
  C3ratioFit->FixParameter(0,9.98651e-01); C3ratioFit->FixParameter(1,2.48247e-02); C3ratioFit->FixParameter(2,3.13184e+01);
  TF1 *C4ratioFit = new TF1("C4ratioFit","[0] + [1]*exp(-pow(x*[2],2))",0,1);
  C4ratioFit->FixParameter(0,9.98432e-01); C4ratioFit->FixParameter(1,9.43508e-02); C4ratioFit->FixParameter(2,2.89368e+01);


  TH2D *TwoParticle_2d[2][2][2];// ch1,ch2,term
  TH1D *TwoParticle[2][2][2];// ch1,ch2,term
  TH3D *Build_TwoParticle_2D[2];// charge
  TH1D *Build_TwoParticle[2];// charge
  TH2D *UnitMult_2d[2][2][2];// ch1,ch2,term
  TH1D *UnitMult[2][2][2];// ch1,ch2,term
  double norm_2[2]={0};
  double norm_2_UM[2]={0};
  //
  TH1D *ThreeParticle[2][2][2][5];// ch1,ch2,ch3,term
  TProfile *K3avg[2][2][2][4];
  TH2D *Build_ThreeParticle_2D[2];// charge
  TH2D *BuildNeg_ThreeParticle_2D[2];// charge
  TH2D *CumulantBuild_ThreeParticle_2D[2];// charge
  TH2D *CumulantBuildNeg_ThreeParticle_2D[2];// charge
  TH1D *TPN_ThreeParticle[2];// charge
  TH1D *Build_ThreeParticle[2];// charge
  TH1D *BuildNeg_ThreeParticle[2];// charge
  TH1D *CumulantBuild_ThreeParticle[2];// charge
  TH1D *CumulantBuildNeg_ThreeParticle[2];// charge
  double norm_3[5]={0};
  //
  TH1D *FourParticle[2][2][2][2][13];// ch1,ch2,ch3,ch4,term
  TProfile *K4avg[2][2][2][2][12];
  TH2D *Build_FourParticle_2D[2];// charge
  TH2D *PrimeBuild_FourParticle_2D[2];// charge
  TH2D *PrimePrimeBuild_FourParticle_2D[2];// charge
  TH2D *CumulantBuild_FourParticle_2D[2];// charge
  TH2D *BuildNeg_FourParticle_2D[2];// charge
  TH2D *PrimeBuildNeg_FourParticle_2D[2];// charge
  TH2D *PrimePrimeBuildNeg_FourParticle_2D[2];// charge
  TH2D *CumulantBuildNeg_FourParticle_2D[2];// charge
  TH1D *Build_FourParticle[2];// charge
  TH1D *PrimeBuild_FourParticle[2];// charge
  TH1D *PrimePrimeBuild_FourParticle[2];// charge
  TH1D *CumulantBuild_FourParticle[2];// charge
  TH1D *BuildNeg_FourParticle[2];// charge
  TH1D *PrimeBuildNeg_FourParticle[2];// charge
  TH1D *PrimePrimeBuildNeg_FourParticle[2];// charge
  TH1D *CumulantBuildNeg_FourParticle[2];// charge
  
  TH1D *TPN_FourParticle[2];// charge
  TH3D *BuildFromFits_3D;// FitType, Gbin, Q4
  TH3D *PrimeBuildFromFits_3D;
  TH3D *PrimePrimeBuildFromFits_3D;
  TH3D *CumulantBuildFromFits_3D;
  TH1D *BuildFromFits1;
  TH1D *PrimeBuildFromFits1;
  TH1D *PrimePrimeBuildFromFits1;
  TH1D *CumulantBuildFromFits1;
  TH1D *BuildFromFits2;
  TH1D *PrimeBuildFromFits2;
  TH1D *PrimePrimeBuildFromFits2;
  TH1D *CumulantBuildFromFits2;

  TH1D *BuildFromFits_M;
  TH1D *PrimeBuildFromFits_M;
  TH1D *PrimePrimeBuildFromFits_M;
  TH1D *CumulantBuildFromFits_M;
  
  double norm_4[13]={0};
  

  gStyle->SetOptStat(0);
  gStyle->SetOptDate(0);
  gStyle->SetOptFit(0111);

  TFile *_file0;
  if(CollisionType==0){// PbPb
    if(MCcase) {
      if(Mbin<=5){
	//_file0 = new TFile("Results/PDC_12a17a_GeneratedSignal.root","READ");
	//_file0 = new TFile("Results/PDC_12a17a_pTSpectrumWeight.root","READ");
	//_file0 = new TFile("Results/PDC_12a17a_Qweights.root","READ");
	_file0 = new TFile("Results/PDC_12a11a_GeneratorLevel.root","READ");
      }
    }else{
      if(Mbin<=1){
	//_file0 = new TFile("Results/PDC_11h_extendedGweights.root","READ");// Preliminary results
	//_file0 = new TFile("Results/PDC_11h_LowPtMultBinning_M0to1.root","READ");
	//_file0 = new TFile("Results/PDC_11h_LowPtMultBinningHighPtConstraint_M0to1.root","READ");
	//_file0 = new TFile("Results/PDC_11h_noMRCnoMuon.root","READ");
	//_file0 = new TFile("Results/PDC_11h_fc2_0p4_0p85.root","READ");
       	//if(FileSetting==0) _file0 = new TFile("Results/PDC_11h_fc2_0p65andK5percentIncrease_WeightFileFSIandMRCcorrected.root","READ");
	//if(FileSetting==0) _file0 = new TFile("Results/PDC_11h_4kappas.root","READ");
	//if(FileSetting==5) _file0 = new TFile("Results/PDC_11h_Cubic_Linear.root","READ");
	//if(FileSetting==6) _file0 = new TFile("Results/PDC_10h_but11hWeights.root","READ");
	//if(FileSetting==7 || FileSetting==8) _file0 = new TFile("Results/PDC_11h_MRC10percIncrease_Muon92percent.root","READ");
	//if(FileSetting==9) _file0 = new TFile("Results/PDC_11h_6kT.root","READ");
	//_file0 = new TFile("Results/PDC_11h_7kT.root","READ");
	//_file0 = new TFile("Results/PDC_11h_2pionSelfBuild_OfflineCorrection_OnlineCorrection.root","READ");
	//_file0 = new TFile("Results/PDC_11h_Corrected2pionSelfBuild_OfflineCorrection_OnlineCorrection.root","READ");
	//_file0 = new TFile("Results/PDC_11h_QinvInterp_Cubic.root","READ");
	//_file0 = new TFile("Results/PDC_11h_LinearInterp_6kT.root","READ");
	//_file0 = new TFile("Results/PDC_11h_35Gsteps.root","READ");
	//_file0 = new TFile("Results/PDC_11h_35Gsteps_Bfield1.root","READ");
	//_file0 = new TFile("Results/PDC_11h_35Gsteps_Bfield2.root","READ");
	_file0 = new TFile("Results/PDC_11h_InterpCorrected_M0to1.root","READ");// new one
      }else{
	//_file0 = new TFile("Results/PDC_11h_M3to10_ZvtxEM.root","READ");
	//_file0 = new TFile("Results/PDC_11h_M2to10_6kT.root","READ");// problematic muon corrections
	//_file0 = new TFile("Results/PDC_11h_6kT_M2to10.root","READ");
	//_file0 = new TFile("Results/PDC_11h_M2to8.root","READ");// new one
	//_file0 = new TFile("Results/PDC_11h_InterpCorrections_Offline_M2to7.root","READ");
	//_file0 = new TFile("Results/PDC_11h_M8to10.root","READ");
	//_file0 = new TFile("Results/PDC_11h_InterpCorrections_M8to9.root","READ");
	_file0 = new TFile("Results/PDC_11h_InterpCorrected_M2to9.root","READ");// new one
      }
    }
  }else if(CollisionType==1){// pPb
    //_file0 = new TFile("Results/PDC_13bc_kINT7_2FitBuilds.root","READ");// used for paper proposal
    //_file0 = new TFile("Results/PDC_13bc_kINT7_BothBuildTypes.root","READ");
    //_file0 = new TFile("Results/PDC_13c_InterpCorrected.root","READ");// new one
    _file0 = new TFile("Results/PDC_13bc_kINT7_BothBuildTypes.root","READ");// new one
  }else{// pp
    //_file0 = new TFile("Results/PDC_10bcde_kMB_2FitBuilds.root","READ");
    //_file0 = new TFile("Results/PDC_10bcde_kMB_Norm0p8to1p5.root","READ");// used for paper proposal
    //_file0 = new TFile("Results/PDC_10d_InterpCorrected.root","READ");
    //_file0 = new TFile("Results/PDC_10d_InterpCorrected_0p18to0p22Norm.root","READ");
    //_file0 = new TFile("Results/PDC_10bcde_kMB_1DBuilds.root","READ");
    //_file0 = new TFile("Results/PDC_10d_InterpCorrected.root","READ");
    _file0 = new TFile("Results/PDC_10bcde_kMB_BothBuildTypes.root","READ");// new one
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
    if(CollisionType==0 && Mbin>1) MyList=(TList*)tdir->Get("FourPionOutput_2");
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
	//cout<<"Low edge = "<<TwoParticle_2d[c1][c2][term]->GetXaxis()->GetBinLowEdge(Ktlowbin)<<"  High edge = "<<TwoParticle_2d[c1][c2][term]->GetXaxis()->GetBinUpEdge(Kthighbin)<<endl;
	norm_2[term] = TwoParticle[c1][c2][term]->Integral(TwoParticle[c1][c2][term]->GetXaxis()->FindBin(NormQcutLow),TwoParticle[c1][c2][term]->GetXaxis()->FindBin(NormQcutHigh));
	//cout<<"2-pion norms  "<<norm_2[term]<<endl;
	TwoParticle[c1][c2][term]->Scale(norm_2[0]/norm_2[term]);
	
	TwoParticle[c1][c2][term]->SetMarkerStyle(20);
	TwoParticle[c1][c2][term]->GetXaxis()->SetTitle("q_{inv} (GeV/c)");
	TwoParticle[c1][c2][term]->GetYaxis()->SetTitle("C_{2}");
	TwoParticle[c1][c2][term]->SetTitle("");
	//
	// 2-pion self builds
	/*if(term==1 && c1==c2){
	  TString *twoBuild = new TString(name2->Data());
	  twoBuild->Append("_Build");
	  Build_TwoParticle_2D[c1] = (TH3D*)MyList->FindObject(twoBuild->Data());
	  twoBuild->Append("_pro");
	  Build_TwoParticle[c1] = (TH1D*)Build_TwoParticle_2D[c1]->ProjectionZ(twoBuild->Data(),5,5, Ktlowbin,Kthighbin);
	  twoBuild->Append("Den");
	  TH1D *tempDenC2 = (TH1D*)Build_TwoParticle_2D[c1]->ProjectionZ(twoBuild->Data(),4,4, Ktlowbin,Kthighbin);
	  Build_TwoParticle[c1]->Add(tempDenC2);
	  Build_TwoParticle[c1]->Divide(tempDenC2);
	  for(int i=0; i<Build_TwoParticle[c1]->GetNbinsX(); i++){
	    Build_TwoParticle[c1]->SetBinError(i+1,0.000001);
	  }
	  }*/
	
	if(Mbin==0 && CollisionType==0 && !MCcase){
	  TString *nameUM=new TString(name2->Data());
	  nameUM->Append("_UnitMult");
	  UnitMult_2d[c1][c2][term] = (TH2D*)MyList->FindObject(nameUM->Data());
	  UnitMult_2d[c1][c2][term]->Sumw2();
	  TString *pronameUM = new TString(nameUM->Data());
	  pronameUM->Append("_pro");
	  UnitMult[c1][c2][term] = (TH1D*)UnitMult_2d[c1][c2][term]->ProjectionY(pronameUM->Data(),13,13);// 11 means 1000 pions
	  norm_2_UM[term] = UnitMult[c1][c2][term]->Integral(UnitMult[c1][c2][term]->GetXaxis()->FindBin(1.0),UnitMult[c1][c2][term]->GetXaxis()->FindBin(1.2));
	  if(norm_2_UM[term] > 0) UnitMult[c1][c2][term]->Scale(norm_2_UM[0]/norm_2_UM[term]);
	  //
	  UnitMult[c1][c2][term]->SetMarkerStyle(20);
	  UnitMult[c1][c2][term]->GetXaxis()->SetTitle("q_{inv} (GeV/c)");
	  UnitMult[c1][c2][term]->GetYaxis()->SetTitle("C_{2}");
	  UnitMult[c1][c2][term]->SetTitle("");
	}
      }// term
      //continue;
      
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
	  nameTPN3->Append("_Build");
	  //
	  TString *nameNegTPN3=new TString(name3->Data());
	  nameNegTPN3->Append("_BuildNeg");
	  //
	  TString *nameBuildCumulant3=new TString(name3->Data());
	  nameBuildCumulant3->Append("_CumulantBuild");
	  //
	  TString *nameBuildNegCumulant3=new TString(name3->Data());
	  nameBuildNegCumulant3->Append("_CumulantBuildNeg");
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
	  //cout<<"3-pion norms  "<<norm_3[term]<<endl;
	  ThreeParticle[c1][c2][c3][term]->Scale(norm_3[0]/norm_3[term]);
	  ThreeParticle[c1][c2][c3][term]->GetXaxis()->SetTitle("Q_{3} (GeV/c)");
	  ThreeParticle[c1][c2][c3][term]->GetYaxis()->SetTitle("Three Pion");
	  ThreeParticle[c1][c2][c3][term]->SetMarkerStyle(20);
	  ThreeParticle[c1][c2][c3][term]->SetTitle("");
	  //
	  ThreeParticle[c1][c2][c3][term]->Rebin(ThreeParticleRebin);
	  
	  //
	  if(term<4){
	    K3avg[c1][c2][c3][term] = (TProfile*)MyList->FindObject(nameK3->Data());
	    K3avg[c1][c2][c3][term]->Rebin(ThreeParticleRebin);
	    K3avg[c1][c2][c3][term]->GetXaxis()->SetTitle("Q_{3} (GeV/c)"); K3avg[c1][c2][c3][term]->GetYaxis()->SetTitle("1/K_{3}");
	    if(MCcase || !FSICorrection) {
	      K3avg[c1][c2][c3][term]->Reset();
	      for(int ii=1; ii<=K3avg[c1][c2][c3][term]->GetNbinsX(); ii++) {
		K3avg[c1][c2][c3][term]->Fill(K3avg[c1][c2][c3][term]->GetXaxis()->GetBinCenter(ii), 1);
	      }	
	    }
	  }
	  //
	  if(term==4 && c1==c2 && c1==c3 && !MCcase){
	    Build_ThreeParticle_2D[c1] = (TH2D*)MyList->FindObject(nameTPN3->Data());
	    Build_ThreeParticle_2D[c1]->Scale(norm_3[0]/norm_3[term]);
	    Build_ThreeParticle_2D[c1]->RebinY(ThreeParticleRebin);
	    //
	    BuildNeg_ThreeParticle_2D[c1] = (TH2D*)MyList->FindObject(nameNegTPN3->Data());
	    BuildNeg_ThreeParticle_2D[c1]->Scale(norm_3[0]/norm_3[term]);
	    BuildNeg_ThreeParticle_2D[c1]->RebinY(ThreeParticleRebin);
	    //
	    CumulantBuild_ThreeParticle_2D[c1] = (TH2D*)MyList->FindObject(nameBuildCumulant3->Data());
	    CumulantBuild_ThreeParticle_2D[c1]->Scale(norm_3[0]/norm_3[term]);
	    CumulantBuild_ThreeParticle_2D[c1]->RebinY(ThreeParticleRebin);
	    //
	    CumulantBuildNeg_ThreeParticle_2D[c1] = (TH2D*)MyList->FindObject(nameBuildNegCumulant3->Data());
	    CumulantBuildNeg_ThreeParticle_2D[c1]->Scale(norm_3[0]/norm_3[term]);
	    CumulantBuildNeg_ThreeParticle_2D[c1]->RebinY(ThreeParticleRebin);
	    //
	    TString *proName=new TString(nameTPN3->Data()); TString *proNameNeg=new TString(nameNegTPN3->Data());
	    proName->Append("_pro");  proNameNeg->Append("_pro");
	    TPN_ThreeParticle[c1] = (TH1D*)Build_ThreeParticle_2D[c1]->ProjectionY(proName->Data(), TPNbin, TPNbin);
	    //
	    for(int binx=1; binx<=Build_ThreeParticle_2D[c1]->GetNbinsX(); binx++) {
	      for(int biny=1; biny<=Build_ThreeParticle_2D[c1]->GetNbinsY(); biny++) {
		Build_ThreeParticle_2D[c1]->SetBinError(binx, biny, 0);
		if(binx!=4){
		  Build_ThreeParticle_2D[c1]->SetBinContent(binx, biny, Build_ThreeParticle_2D[c1]->GetBinContent(binx, biny) + Build_ThreeParticle_2D[c1]->GetBinContent(4, biny));
		  BuildNeg_ThreeParticle_2D[c1]->SetBinContent(binx, biny, BuildNeg_ThreeParticle_2D[c1]->GetBinContent(binx, biny) + BuildNeg_ThreeParticle_2D[c1]->GetBinContent(4, biny));
		  //
		  CumulantBuild_ThreeParticle_2D[c1]->SetBinContent(binx, biny, CumulantBuild_ThreeParticle_2D[c1]->GetBinContent(binx, biny) + Build_ThreeParticle_2D[c1]->GetBinContent(4, biny));
		  CumulantBuildNeg_ThreeParticle_2D[c1]->SetBinContent(binx, biny, CumulantBuildNeg_ThreeParticle_2D[c1]->GetBinContent(binx, biny) + BuildNeg_ThreeParticle_2D[c1]->GetBinContent(4, biny));
		}
	      }
	    }
	    
	    //
	    proName->Append("_FullWeight"); proNameNeg->Append("_FullWeight");
	    Build_ThreeParticle[c1] = (TH1D*)Build_ThreeParticle_2D[c1]->ProjectionY(proName->Data(), Gbin, Gbin);
	    BuildNeg_ThreeParticle[c1] = (TH1D*)BuildNeg_ThreeParticle_2D[c1]->ProjectionY(proNameNeg->Data(), Gbin, Gbin);
	    proName->Append("_cumulant"); proNameNeg->Append("_cumulant");
	    CumulantBuild_ThreeParticle[c1] = (TH1D*)CumulantBuild_ThreeParticle_2D[c1]->ProjectionY(proName->Data(), Gbin, Gbin);
	    CumulantBuildNeg_ThreeParticle[c1] = (TH1D*)CumulantBuildNeg_ThreeParticle_2D[c1]->ProjectionY(proNameNeg->Data(), Gbin, Gbin);
	    proName->Append("_FullWeightDen"); proNameNeg->Append("_FullWeightDen");
	    TH1D *tempDen = (TH1D*)Build_ThreeParticle_2D[c1]->ProjectionY(proName->Data(), 4, 4);
	    TH1D *tempDenNeg = (TH1D*)BuildNeg_ThreeParticle_2D[c1]->ProjectionY(proNameNeg->Data(), 4, 4);
	    // Add Pos with Neg weights
	    tempDen->Add(tempDenNeg);
	    Build_ThreeParticle[c1]->Add(BuildNeg_ThreeParticle[c1]);
	    CumulantBuild_ThreeParticle[c1]->Add(CumulantBuildNeg_ThreeParticle[c1]);
	    //
	    Build_ThreeParticle[c1]->Divide(tempDen);
	    CumulantBuild_ThreeParticle[c1]->Divide(tempDen);
	    Build_ThreeParticle[c1]->SetLineColor(4);
	    CumulantBuild_ThreeParticle[c1]->SetLineColor(1);
	    //
	   
	    /*proName->Append("_2"); proNameNeg->Append("_2"); 
	    CumulantBuild_ThreeParticle[c1] = (TH1D*)Build_ThreeParticle_2D[c1]->ProjectionY(proName->Data(), TPNbin, TPNbin);
	    CumulantBuildNeg_ThreeParticle[c1] = (TH1D*)BuildNeg_ThreeParticle_2D[c1]->ProjectionY(proNameNeg->Data(), TPNbin, TPNbin);
	    CumulantBuild_ThreeParticle[c1]->Add(CumulantBuildNeg_ThreeParticle[c1]);
	    CumulantBuild_ThreeParticle[c1]->Add(tempDen,-1);
	    CumulantBuild_ThreeParticle[c1]->Scale(2);
	    CumulantBuild_ThreeParticle[c1]->Add(tempDen,+1);
	    CumulantBuild_ThreeParticle[c1]->Divide(tempDen);
	    CumulantBuild_ThreeParticle[c1]->SetLineColor(1);
	    */
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
	    nameTPN4->Append("_Build");
	    //
	    TString *nameTPN4p=new TString(name4->Data());
	    nameTPN4p->Append("_PrimeBuild");
	    //
	    TString *nameTPN4pp=new TString(name4->Data());
	    nameTPN4pp->Append("_PrimePrimeBuild");
	    //
	    TString *nameTPN4c=new TString(name4->Data());
	    nameTPN4c->Append("_CumulantBuild");
	    //
	    TString *nameNegTPN4=new TString(name4->Data());
	    nameNegTPN4->Append("_BuildNeg");
	    //
	    TString *nameNegTPN4p=new TString(name4->Data());
	    nameNegTPN4p->Append("_PrimeBuildNeg");
	    //
	    TString *nameNegTPN4pp=new TString(name4->Data());
	    nameNegTPN4pp->Append("_PrimePrimeBuildNeg");
	    //
	    TString *nameNegTPN4c=new TString(name4->Data());
	    nameNegTPN4c->Append("_CumulantBuildNeg");
	    //
	    TString *nameFitBuild=new TString(name4->Data());
	    nameFitBuild->Append("_BuildFromFits");
	    //
	    TString *namePrimeFitBuild=new TString(name4->Data());
	    namePrimeFitBuild->Append("_PrimeBuildFromFits");
	    //
	    TString *namePrimePrimeFitBuild=new TString(name4->Data());
	    namePrimePrimeFitBuild->Append("_PrimePrimeBuildFromFits");
	    //
	    TString *nameCumulantFitBuild=new TString(name4->Data());
	    nameCumulantFitBuild->Append("_CumulantBuildFromFits");
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
	      //if(c1==c2 && c1==c3 && c1==c4) cout<<term<<"  "<<norm_4[0]/norm_4[term]<<endl;
	      FourParticle[c1][c2][c3][c4][term] = (TH1D*)MyList->FindObject(name4->Data());
	      FourParticle[c1][c2][c3][c4][term]->Sumw2();
	      FourParticle[c1][c2][c3][c4][term]->Scale(norm_4[0]/norm_4[term]);
	      FourParticle[c1][c2][c3][c4][term]->GetXaxis()->SetTitle("Q_{4} (GeV/c)");
	      FourParticle[c1][c2][c3][c4][term]->GetYaxis()->SetTitle("Four Pion");
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
	      K4avg[c1][c2][c3][c4][term]->GetXaxis()->SetTitle("Q_{4} (GeV/c)"); K4avg[c1][c2][c3][c4][term]->GetYaxis()->SetTitle("1/K_{4}");
	      if(MCcase || !FSICorrection) {
		K4avg[c1][c2][c3][c4][term]->Reset();
		for(int ii=1; ii<=K4avg[c1][c2][c3][c4][term]->GetNbinsX(); ii++) {
		  K4avg[c1][c2][c3][c4][term]->Fill(K4avg[c1][c2][c3][c4][term]->GetXaxis()->GetBinCenter(ii), 1);
		}	
	      }
	    }
	    if(term==12 && c1==c2 && c1==c3 && c1==c4 && !MCcase){
	      
	      Build_FourParticle_2D[c1] = (TH2D*)MyList->FindObject(nameTPN4->Data());
	      Build_FourParticle_2D[c1]->Scale(norm_4[0]/norm_4[term]);
	      Build_FourParticle_2D[c1]->RebinY(FourParticleRebin);
	      PrimeBuild_FourParticle_2D[c1] = (TH2D*)MyList->FindObject(nameTPN4p->Data());
	      PrimeBuild_FourParticle_2D[c1]->Scale(norm_4[0]/norm_4[term]);
	      PrimeBuild_FourParticle_2D[c1]->RebinY(FourParticleRebin);
	      PrimePrimeBuild_FourParticle_2D[c1] = (TH2D*)MyList->FindObject(nameTPN4pp->Data());
	      PrimePrimeBuild_FourParticle_2D[c1]->Scale(norm_4[0]/norm_4[term]);
	      PrimePrimeBuild_FourParticle_2D[c1]->RebinY(FourParticleRebin);
	      CumulantBuild_FourParticle_2D[c1] = (TH2D*)MyList->FindObject(nameTPN4c->Data());
	      CumulantBuild_FourParticle_2D[c1]->Scale(norm_4[0]/norm_4[term]);
	      CumulantBuild_FourParticle_2D[c1]->RebinY(FourParticleRebin);
	      //
	      BuildNeg_FourParticle_2D[c1] = (TH2D*)MyList->FindObject(nameNegTPN4->Data());
	      BuildNeg_FourParticle_2D[c1]->Scale(norm_4[0]/norm_4[term]);
	      BuildNeg_FourParticle_2D[c1]->RebinY(FourParticleRebin);
	      PrimeBuildNeg_FourParticle_2D[c1] = (TH2D*)MyList->FindObject(nameNegTPN4p->Data());
	      PrimeBuildNeg_FourParticle_2D[c1]->Scale(norm_4[0]/norm_4[term]);
	      PrimeBuildNeg_FourParticle_2D[c1]->RebinY(FourParticleRebin);
	      PrimePrimeBuildNeg_FourParticle_2D[c1] = (TH2D*)MyList->FindObject(nameNegTPN4pp->Data());
	      PrimePrimeBuildNeg_FourParticle_2D[c1]->Scale(norm_4[0]/norm_4[term]);
	      PrimePrimeBuildNeg_FourParticle_2D[c1]->RebinY(FourParticleRebin);
	      CumulantBuildNeg_FourParticle_2D[c1] = (TH2D*)MyList->FindObject(nameNegTPN4c->Data());
	      CumulantBuildNeg_FourParticle_2D[c1]->Scale(norm_4[0]/norm_4[term]);
	      CumulantBuildNeg_FourParticle_2D[c1]->RebinY(FourParticleRebin);
	      //
	      if(c1==0 && !MCcase){
		BuildFromFits_3D = (TH3D*)MyList->FindObject(nameFitBuild->Data());
		BuildFromFits_3D->Scale(norm_4[0]/norm_4[term]);
		BuildFromFits_3D->RebinZ(FourParticleRebin);
		//
		PrimeBuildFromFits_3D = (TH3D*)MyList->FindObject(namePrimeFitBuild->Data());
		PrimeBuildFromFits_3D->Scale(norm_4[0]/norm_4[term]);
		PrimeBuildFromFits_3D->RebinZ(FourParticleRebin);
		//
		PrimePrimeBuildFromFits_3D = (TH3D*)MyList->FindObject(namePrimePrimeFitBuild->Data());
		PrimePrimeBuildFromFits_3D->Scale(norm_4[0]/norm_4[term]);
		PrimePrimeBuildFromFits_3D->RebinZ(FourParticleRebin);
		//
		CumulantBuildFromFits_3D = (TH3D*)MyList->FindObject(nameCumulantFitBuild->Data());
		CumulantBuildFromFits_3D->Scale(norm_4[0]/norm_4[term]);
		CumulantBuildFromFits_3D->RebinZ(FourParticleRebin);
	      }
	      //
	      for(int binx=1; binx<=Build_FourParticle_2D[c1]->GetNbinsX(); binx++) {
		for(int biny=1; biny<=Build_FourParticle_2D[c1]->GetNbinsY(); biny++) {
		  Build_FourParticle_2D[c1]->SetBinError(binx, biny, 0);
		  PrimeBuild_FourParticle_2D[c1]->SetBinError(binx, biny, 0);
		  PrimePrimeBuild_FourParticle_2D[c1]->SetBinError(binx, biny, 0);
		  CumulantBuild_FourParticle_2D[c1]->SetBinError(binx, biny, 0);
		  if(binx!=4){
		    Build_FourParticle_2D[c1]->SetBinContent(binx, biny, Build_FourParticle_2D[c1]->GetBinContent(binx, biny) + Build_FourParticle_2D[c1]->GetBinContent(4, biny));
		    PrimeBuild_FourParticle_2D[c1]->SetBinContent(binx, biny, PrimeBuild_FourParticle_2D[c1]->GetBinContent(binx, biny) + PrimeBuild_FourParticle_2D[c1]->GetBinContent(4, biny));
		    PrimePrimeBuild_FourParticle_2D[c1]->SetBinContent(binx, biny, PrimePrimeBuild_FourParticle_2D[c1]->GetBinContent(binx, biny) + PrimePrimeBuild_FourParticle_2D[c1]->GetBinContent(4, biny));
		    CumulantBuild_FourParticle_2D[c1]->SetBinContent(binx, biny, CumulantBuild_FourParticle_2D[c1]->GetBinContent(binx, biny) + CumulantBuild_FourParticle_2D[c1]->GetBinContent(4, biny));
		    //
		    BuildNeg_FourParticle_2D[c1]->SetBinContent(binx, biny, BuildNeg_FourParticle_2D[c1]->GetBinContent(binx, biny) + BuildNeg_FourParticle_2D[c1]->GetBinContent(4, biny));
		    PrimeBuildNeg_FourParticle_2D[c1]->SetBinContent(binx, biny, PrimeBuildNeg_FourParticle_2D[c1]->GetBinContent(binx, biny) + PrimeBuildNeg_FourParticle_2D[c1]->GetBinContent(4, biny));
		    PrimePrimeBuildNeg_FourParticle_2D[c1]->SetBinContent(binx, biny, PrimePrimeBuildNeg_FourParticle_2D[c1]->GetBinContent(binx, biny) + PrimePrimeBuildNeg_FourParticle_2D[c1]->GetBinContent(4, biny));
		    CumulantBuildNeg_FourParticle_2D[c1]->SetBinContent(binx, biny, CumulantBuildNeg_FourParticle_2D[c1]->GetBinContent(binx, biny) + CumulantBuildNeg_FourParticle_2D[c1]->GetBinContent(4, biny));
		    
		  }
		}
	      }
	      TString *proName=new TString(nameTPN4->Data()); 
	      TString *proNamePrime=new TString(nameTPN4p->Data()); 
	      TString *proNamePrimePrime=new TString(nameTPN4pp->Data()); 
	      TString *proNameCumulant=new TString(nameTPN4c->Data()); 
	      TString *proNegName=new TString(nameNegTPN4->Data());
	      TString *proNegNamePrime=new TString(nameNegTPN4p->Data());
	      TString *proNegNamePrimePrime=new TString(nameNegTPN4pp->Data());
	      TString *proNegNameCumulant=new TString(nameNegTPN4c->Data());
	      TString *proNameFitBuild=new TString(nameFitBuild->Data()); 
	      TString *proNamePrimeFitBuild=new TString(namePrimeFitBuild->Data());
	      TString *proNamePrimePrimeFitBuild=new TString(namePrimePrimeFitBuild->Data());
	      TString *proNameCumulantFitBuild=new TString(nameCumulantFitBuild->Data());
	      proName->Append("_pro"); proNamePrime->Append("_pro"); proNamePrimePrime->Append("_pro"); proNameCumulant->Append("_pro");
	      proNegName->Append("_pro"); proNegNamePrime->Append("_pro"); proNegNamePrimePrime->Append("_pro"); proNegNameCumulant->Append("_pro"); 
	      proNameFitBuild->Append("_pro"); 
	      proNamePrimeFitBuild->Append("_pro");
	      proNamePrimePrimeFitBuild->Append("_pro");
	      proNameCumulantFitBuild->Append("_pro"); 

	      TPN_FourParticle[c1] = (TH1D*)Build_FourParticle_2D[c1]->ProjectionY("r4DenPro", TPNbin, TPNbin);
	      //
	      Build_FourParticle[c1] = (TH1D*)Build_FourParticle_2D[c1]->ProjectionY(proName->Data(), Gbin, Gbin);
	      PrimeBuild_FourParticle[c1] = (TH1D*)PrimeBuild_FourParticle_2D[c1]->ProjectionY(proNamePrime->Data(), Gbin, Gbin);
	      PrimePrimeBuild_FourParticle[c1] = (TH1D*)PrimePrimeBuild_FourParticle_2D[c1]->ProjectionY(proNamePrimePrime->Data(), Gbin, Gbin);
	      CumulantBuild_FourParticle[c1] = (TH1D*)CumulantBuild_FourParticle_2D[c1]->ProjectionY(proNameCumulant->Data(), Gbin, Gbin);
	      BuildNeg_FourParticle[c1] = (TH1D*)BuildNeg_FourParticle_2D[c1]->ProjectionY(proNegName->Data(), Gbin, Gbin);
	      PrimeBuildNeg_FourParticle[c1] = (TH1D*)PrimeBuildNeg_FourParticle_2D[c1]->ProjectionY(proNegNamePrime->Data(), Gbin, Gbin);
	      PrimePrimeBuildNeg_FourParticle[c1] = (TH1D*)PrimePrimeBuildNeg_FourParticle_2D[c1]->ProjectionY(proNegNamePrimePrime->Data(), Gbin, Gbin);
	      CumulantBuildNeg_FourParticle[c1] = (TH1D*)CumulantBuildNeg_FourParticle_2D[c1]->ProjectionY(proNegNameCumulant->Data(), Gbin, Gbin);
	      //
	      proName->Append("_FullWeightDen"); proNegName->Append("_FullWeightDen");
	      TH1D *tempDen = (TH1D*)Build_FourParticle_2D[c1]->ProjectionY(proName->Data(), 4, 4);
	      TH1D *tempDenNeg = (TH1D*)BuildNeg_FourParticle_2D[c1]->ProjectionY(proNegName->Data(), 4, 4);
	      //
	      // Add Pos and Neg Weights
	      tempDen->Add(tempDenNeg);
	      Build_FourParticle[c1]->Add(BuildNeg_FourParticle[c1]);
	      PrimeBuild_FourParticle[c1]->Add(PrimeBuildNeg_FourParticle[c1]);
	      PrimePrimeBuild_FourParticle[c1]->Add(PrimePrimeBuildNeg_FourParticle[c1]);
	      CumulantBuild_FourParticle[c1]->Add(CumulantBuildNeg_FourParticle[c1]);
	      //
	      Build_FourParticle[c1]->Divide(tempDen);
	      PrimeBuild_FourParticle[c1]->Divide(tempDen);
	      PrimePrimeBuild_FourParticle[c1]->Divide(tempDen);
	      CumulantBuild_FourParticle[c1]->Divide(tempDen);
	      //
	      Build_FourParticle[c1]->SetLineColor(4); Build_FourParticle[c1]->SetLineStyle(2);
	      PrimeBuild_FourParticle[c1]->SetLineColor(2); PrimeBuild_FourParticle[c1]->SetLineStyle(2);
	      PrimePrimeBuild_FourParticle[c1]->SetLineColor(6); PrimePrimeBuild_FourParticle[c1]->SetLineStyle(2);
	      CumulantBuild_FourParticle[c1]->SetLineColor(1); CumulantBuild_FourParticle[c1]->SetLineStyle(2);
	      
	      if(c1==0 && !MCcase){
		BuildFromFits1 = (TH1D*)BuildFromFits_3D->ProjectionZ(proNameFitBuild->Data(), 1,1, Gbin, Gbin);
		TH1D *tempDen2 = (TH1D*)BuildFromFits_3D->ProjectionZ("tempDen2", 1,1, 4, 4);
		BuildFromFits1->Add(tempDen2);
		BuildFromFits1->Divide(tempDen2);
		BuildFromFits1->SetLineColor(4); BuildFromFits1->SetLineStyle(2);
		//
		PrimeBuildFromFits1 = (TH1D*)PrimeBuildFromFits_3D->ProjectionZ(proNamePrimeFitBuild->Data(), 1,1, Gbin, Gbin);
		PrimeBuildFromFits1->Add(tempDen2);
		PrimeBuildFromFits1->Divide(tempDen2);
		PrimeBuildFromFits1->SetLineColor(2); PrimeBuildFromFits1->SetLineStyle(2);
		//
		PrimePrimeBuildFromFits1 = (TH1D*)PrimePrimeBuildFromFits_3D->ProjectionZ(proNamePrimePrimeFitBuild->Data(), 1,1, Gbin, Gbin);
		PrimePrimeBuildFromFits1->Add(tempDen2);
		PrimePrimeBuildFromFits1->Divide(tempDen2);
		PrimePrimeBuildFromFits1->SetLineColor(6); PrimePrimeBuildFromFits1->SetLineStyle(2);
		//
		CumulantBuildFromFits1 = (TH1D*)CumulantBuildFromFits_3D->ProjectionZ(proNameCumulantFitBuild->Data(), 1,1, Gbin, Gbin);
		CumulantBuildFromFits1->Add(tempDen2);
		CumulantBuildFromFits1->Divide(tempDen2);
		CumulantBuildFromFits1->SetLineColor(1); CumulantBuildFromFits1->SetLineStyle(2);
		//////
		proNameFitBuild->Append("_2"); proNamePrimeFitBuild->Append("_2"); proNamePrimePrimeFitBuild->Append("_2"); proNameCumulantFitBuild->Append("_2");
		BuildFromFits2 = (TH1D*)BuildFromFits_3D->ProjectionZ(proNameFitBuild->Data(), 2,2, Gbin, Gbin);
		BuildFromFits2->Add(tempDen2);
		BuildFromFits2->Divide(tempDen2);
		BuildFromFits2->SetLineColor(4);
		//
		PrimeBuildFromFits2 = (TH1D*)PrimeBuildFromFits_3D->ProjectionZ(proNamePrimeFitBuild->Data(), 2,2, Gbin, Gbin);
		PrimeBuildFromFits2->Add(tempDen2);
		PrimeBuildFromFits2->Divide(tempDen2);
		PrimeBuildFromFits2->SetLineColor(2);
		//
		PrimePrimeBuildFromFits2 = (TH1D*)PrimePrimeBuildFromFits_3D->ProjectionZ(proNamePrimePrimeFitBuild->Data(), 2,2, Gbin, Gbin);
		PrimePrimeBuildFromFits2->Add(tempDen2);
		PrimePrimeBuildFromFits2->Divide(tempDen2);
		PrimePrimeBuildFromFits2->SetLineColor(6);
		//
		CumulantBuildFromFits2 = (TH1D*)CumulantBuildFromFits_3D->ProjectionZ(proNameCumulantFitBuild->Data(), 2,2, Gbin, Gbin);
		CumulantBuildFromFits2->Add(tempDen2);
		CumulantBuildFromFits2->Divide(tempDen2);
		CumulantBuildFromFits2->SetLineColor(1);
		//
		//
		//
		BuildFromFits_M = (TH1D*)BuildFromFits1->Clone();
		for(int i=1; i<=BuildFromFits_M->GetNbinsX(); i++){
		  BuildFromFits_M->SetBinContent(i, (BuildFromFits1->GetBinContent(i) + BuildFromFits2->GetBinContent(i))/2.);
		  BuildFromFits_M->SetBinError(i, (BuildFromFits1->GetBinContent(i) - BuildFromFits2->GetBinContent(i))/2.);
		}
		BuildFromFits_M->SetMarkerSize(0);
		BuildFromFits_M->SetFillColor(kBlue-10); BuildFromFits_M->SetMarkerColor(kBlue-10);
		//
		PrimeBuildFromFits_M = (TH1D*)PrimeBuildFromFits1->Clone();
		for(int i=1; i<=PrimeBuildFromFits_M->GetNbinsX(); i++){
		  PrimeBuildFromFits_M->SetBinContent(i, (PrimeBuildFromFits1->GetBinContent(i) + PrimeBuildFromFits2->GetBinContent(i))/2.);
		  PrimeBuildFromFits_M->SetBinError(i, (PrimeBuildFromFits1->GetBinContent(i) - PrimeBuildFromFits2->GetBinContent(i))/2.);
		}
		PrimeBuildFromFits_M->SetMarkerSize(0);
		PrimeBuildFromFits_M->SetFillColor(kRed-10); PrimeBuildFromFits_M->SetMarkerColor(kRed-10);
		//
		PrimePrimeBuildFromFits_M = (TH1D*)PrimePrimeBuildFromFits1->Clone();
		for(int i=1; i<=PrimePrimeBuildFromFits_M->GetNbinsX(); i++){
		  PrimePrimeBuildFromFits_M->SetBinContent(i, (PrimePrimeBuildFromFits1->GetBinContent(i) + PrimePrimeBuildFromFits2->GetBinContent(i))/2.);
		  PrimePrimeBuildFromFits_M->SetBinError(i, (PrimePrimeBuildFromFits1->GetBinContent(i) - PrimePrimeBuildFromFits2->GetBinContent(i))/2.);
		}
		PrimePrimeBuildFromFits_M->SetMarkerSize(0);
		PrimePrimeBuildFromFits_M->SetFillColor(kMagenta-10); PrimePrimeBuildFromFits_M->SetMarkerColor(kMagenta-10);
		//
		CumulantBuildFromFits_M = (TH1D*)CumulantBuildFromFits1->Clone();
		for(int i=1; i<=CumulantBuildFromFits_M->GetNbinsX(); i++){
		  CumulantBuildFromFits_M->SetBinContent(i, (CumulantBuildFromFits1->GetBinContent(i) + CumulantBuildFromFits2->GetBinContent(i))/2.);
		  CumulantBuildFromFits_M->SetBinError(i, (CumulantBuildFromFits1->GetBinContent(i) - CumulantBuildFromFits2->GetBinContent(i))/2.);
		}
		CumulantBuildFromFits_M->SetMarkerSize(0);
		CumulantBuildFromFits_M->SetFillColor(kGray); CumulantBuildFromFits_M->SetMarkerColor(kGray);
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
  C2->GetXaxis()->SetRangeUser(0,0.06);//0 to 0.1
  C2->SetMinimum(0.95); C2->SetMaximum(2.1);// 0.98 to 2.5
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
    // <qinv> is estimated as the bin center.  Should be somewhat shifted to larger values. The shift is approximately 0.0005 
    if(FSICorrection) K2 = FSICorrelation(ch1_2, ch2_2, C2QS->GetXaxis()->GetBinCenter(bin) + 0.0005);
    //K2 *= 0.95;
    C2QS->SetBinContent(bin, C2QS->GetBinContent(bin)/K2);
    C2QS->SetBinError(bin, C2QS->GetBinError(bin)/K2);
  }
  C2QS->Divide(TERM2_2);
  if(SameCharge && MuonCorrection) C2QS->Multiply(C2muonCorrectionSC);
  if(!SameCharge && MuonCorrection) C2QS->Multiply(C2muonCorrectionMC);
  C2QS->SetMarkerColor(2); C2QS->SetLineColor(2);
  if(!MCcase) C2QS->Draw("same");

  //for(int i=1; i<14; i++) cout<<C2QS->GetBinContent(i)<<", ";
  //cout<<endl;
  float C2QSmm_fc2_0p8[13]={0, 1.74456, 1.48544, 1.30077, 1.17806, 1.10706, 1.06533, 1.04121, 1.02598, 1.01656, 1.01048, 1.00706, 1.00516};
  float C2QSmm_fc2_0p75[13]={0, 1.76929, 1.50671, 1.31458, 1.18612, 1.11172, 1.06797, 1.0427, 1.02677, 1.01694, 1.0106, 1.00706, 1.00511};
  float C2QSmm_fc2_0p7[13]={0, 1.79756, 1.53103, 1.33036, 1.19534, 1.11703, 1.07098, 1.0444, 1.02768, 1.01736, 1.01073, 1.00706, 1.00504};
  float C2QSmm_fc2_0p65[13]={0, 1.83017, 1.55908, 1.34856, 1.20597, 1.12317, 1.07445, 1.04637, 1.02872, 1.01785, 1.01089, 1.00705, 1.00497};
  float C2QSmm_fc2_0p6[13]={0, 1.86823, 1.59181, 1.3698, 1.21838, 1.13033, 1.0785, 1.04867, 1.02993, 1.01843, 1.01107, 1.00705, 1.00489};
  //
  float C2QSmp_fc2_0p8[13]={0, 0.970557, 0.986022, 0.990979, 0.994373, 0.996436, 0.997731, 0.998119, 0.998702, 0.998714, 0.998899, 0.998943, 0.99898};
  float C2QSmp_fc2_0p75[13]={0, 0.983776, 0.993128, 0.995315, 0.997216, 0.998421, 0.999181, 0.999196, 0.999532, 0.999353, 0.999404, 0.999346, 0.999304};
  float C2QSmp_fc2_0p7[13]={0, 0.998884, 1.00125, 1.00027, 1.00046, 1.00069, 1.00084, 1.00043, 1.00048, 1.00008, 0.99998, 0.999807, 0.999675};
  float C2QSmp_fc2_0p65[13]={0, 1.01632, 1.01062, 1.00599, 1.00421, 1.00331, 1.00275, 1.00185, 1.00157, 1.00093, 1.00065, 1.00034, 1.0001};
  float C2QSmp_fc2_0p6[13]={0, 1.03665, 1.02155, 1.01266, 1.00859, 1.00636, 1.00498, 1.00351, 1.00285, 1.00191, 1.00142, 1.00096, 1.0006};
  //
  //
  C2QS->GetYaxis()->SetTitle("C_{2}^{QS}");
  TH1D *C2mm_0p8 = (TH1D*)C2QS->Clone(); TH1D *C2mm_0p75 = (TH1D*)C2QS->Clone(); TH1D *C2mm_0p7 = (TH1D*)C2QS->Clone();
  TH1D *C2mm_0p65 = (TH1D*)C2QS->Clone(); TH1D *C2mm_0p6 = (TH1D*)C2QS->Clone();
  TH1D *C2mp_0p8 = (TH1D*)C2QS->Clone(); TH1D *C2mp_0p75 = (TH1D*)C2QS->Clone(); TH1D *C2mp_0p7 = (TH1D*)C2QS->Clone();
  TH1D *C2mp_0p65 = (TH1D*)C2QS->Clone(); TH1D *C2mp_0p6 = (TH1D*)C2QS->Clone();
  //
  for(int i=1; i<14; i++){
    C2mm_0p8->SetBinContent(i,C2QSmm_fc2_0p8[i-1]); C2mm_0p8->SetBinError(i, C2QS->GetBinError(i));
    C2mm_0p75->SetBinContent(i,C2QSmm_fc2_0p75[i-1]); C2mm_0p75->SetBinError(i, C2QS->GetBinError(i));
    C2mm_0p7->SetBinContent(i,C2QSmm_fc2_0p7[i-1]); C2mm_0p7->SetBinError(i, C2QS->GetBinError(i));
    C2mm_0p65->SetBinContent(i,C2QSmm_fc2_0p65[i-1]); C2mm_0p65->SetBinError(i, C2QS->GetBinError(i));
    C2mm_0p6->SetBinContent(i,C2QSmm_fc2_0p6[i-1]); C2mm_0p6->SetBinError(i, C2QS->GetBinError(i));
    //
    C2mp_0p8->SetBinContent(i,C2QSmp_fc2_0p8[i-1]); C2mp_0p8->SetBinError(i, C2QS->GetBinError(i));
    C2mp_0p75->SetBinContent(i,C2QSmp_fc2_0p75[i-1]); C2mp_0p75->SetBinError(i, C2QS->GetBinError(i));
    C2mp_0p7->SetBinContent(i,C2QSmp_fc2_0p7[i-1]); C2mp_0p7->SetBinError(i, C2QS->GetBinError(i));
    C2mp_0p65->SetBinContent(i,C2QSmp_fc2_0p65[i-1]); C2mp_0p65->SetBinError(i, C2QS->GetBinError(i));
    C2mp_0p6->SetBinContent(i,C2QSmp_fc2_0p6[i-1]); C2mp_0p6->SetBinError(i, C2QS->GetBinError(i));
  }
  C2mm_0p8->SetMarkerColor(1); C2mm_0p8->SetLineColor(1); C2mm_0p75->SetMarkerColor(2); C2mm_0p75->SetLineColor(2); 
  C2mm_0p7->SetMarkerColor(4); C2mm_0p7->SetLineColor(4); C2mm_0p65->SetMarkerColor(6); C2mm_0p65->SetLineColor(6);
  C2mm_0p6->SetMarkerColor(7); C2mm_0p6->SetLineColor(7); 
  C2mp_0p8->SetMarkerColor(1); C2mp_0p8->SetLineColor(1); C2mp_0p75->SetMarkerColor(2); C2mp_0p75->SetLineColor(2);
  C2mp_0p7->SetMarkerColor(4); C2mp_0p7->SetLineColor(4); C2mp_0p65->SetMarkerColor(6); C2mp_0p65->SetLineColor(6);
  C2mp_0p6->SetMarkerColor(7); C2mp_0p6->SetLineColor(7);
  //C2mm_0p8->Draw(); C2mm_0p75->Draw("same"); C2mm_0p7->Draw("same"); C2mm_0p65->Draw("same"); C2mm_0p6->Draw("same");
  //C2mp_0p8->Draw(); C2mp_0p75->Draw("same"); C2mp_0p7->Draw("same"); C2mp_0p65->Draw("same"); C2mp_0p6->Draw("same");
  //legend1->AddEntry(C2mm_0p8,"f_{c}^{2}=0.8","p"); legend1->AddEntry(C2mm_0p75,"f_{c}^{2}=0.75","p"); 
  //legend1->AddEntry(C2mm_0p7,"f_{c}^{2}=0.7","p"); legend1->AddEntry(C2mm_0p65,"f_{c}^{2}=0.65","p"); 
  //legend1->AddEntry(C2mm_0p6,"f_{c}^{2}=0.6","p"); 
  //
  //legend1->Draw("same");
  
  //cout.precision(3);
  /*Build_TwoParticle[ch1_2]->SetLineColor(4); 
  Build_TwoParticle[ch1_2]->Draw("same");
  TH1D *Ratio_C2=(TH1D*)C2QS->Clone();
  for(int i=2; i<=60; i++){
    //cout<<(C2QS->GetBinContent(i)-1) / (Build_TwoParticle[ch1_2]->GetBinContent(i)-1)<<endl;
    Ratio_C2->SetBinContent(i, (Ratio_C2->GetBinContent(i)-1) / (Build_TwoParticle[ch1_2]->GetBinContent(i)-1));
    //cout<<Ratio_C2->GetBinContent(i)<<", ";
  }
  cout<<endl;
  //Ratio_C2->Divide(Build_TwoParticle[ch1_2]);
  Ratio_C2->SetMinimum(0.9); Ratio_C2->SetMaximum(1.35); Ratio_C2->GetXaxis()->SetRangeUser(0,0.07);
  Ratio_C2->GetYaxis()->SetTitle("(Measured-1) / (Built-1)"); Ratio_C2->GetYaxis()->SetTitleOffset(1.4);
  //Ratio_C2->Draw();
  //TF1 *GaussResidue=new TF1("GaussResidue","[0]+[1]*exp(-pow([2]*x/0.19733,2))",0,0.5);
  //GaussResidue->SetParameter(0,1); GaussResidue->SetParameter(1,0.1); GaussResidue->SetParameter(2,10);
  //Ratio_C2->Fit(GaussResidue,"","",0,0.05);
  */

  //return;
  
  //TF1 *ExpFit = new TF1("ExpFit","[0]*(1+[1]*exp(-pow([2]*x/0.19733,2)))",0,0.5);
  //ExpFit->SetParameter(0,1);
  //ExpFit->SetParameter(1,1);
  //ExpFit->SetParameter(2,12);
  //C2QS->Fit(ExpFit,"","",0.005,0.06);
  //ExpFit->Draw("same");

  // To visualize the Qcut and Norm Regions
  //TH1D *QcutRegion = new TH1D("QcutRegion","",400,0,2);
  //TH1D *NormRegion1 = new TH1D("NormRegion1","",400,0,2);  
  //TH1D *NormRegion2 = new TH1D("NormRegion2","",400,0,2);  
  //TERM1_2->Scale(1/TERM1_2->Integral());
  //for(int bin=1; bin<=16; bin++) QcutRegion->SetBinContent(bin,TERM1_2->GetBinContent(bin));
  //for(int bin=213; bin<=220; bin++) NormRegion1->SetBinContent(bin,Two_ex[0][0][0][0]->GetBinContent(bin));
  //for(int bin=31; bin<=35; bin++) NormRegion2->SetBinContent(bin,Two_ex[0][0][0][0]->GetBinContent(bin));
  //QcutRegion->SetFillColor(4);
  //NormRegion1->SetFillColor(2);
  //NormRegion2->SetFillColor(3);
  //TERM1_2->GetYaxis()->SetTitle("Pair probability");
  //TERM1_2->GetYaxis()->SetTitleOffset(1.3);
  //TERM1_2->GetXaxis()->SetRangeUser(0,1.6);
  //TERM1_2->Draw();
  //QcutRegion->Draw("same");
  //cout<<TERM1_2->Integral(1,16) / TERM1_2->Integral()<<endl;

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

  Unity->Draw("same");

  TF1 *ChaoticLimitC2=new TF1("ChaoticLimitC2","2",0,1);
  ChaoticLimitC2->Draw("same");

  /*if(Mbin==0 && CollisionType==0){
    TH1D *UnitMultC2=(TH1D*)UnitMult[0][0][0]->Clone();
    TH1D *UnitMultC2den=(TH1D*)UnitMult[0][0][1]->Clone();
    UnitMultC2->Divide(UnitMultC2den);
    UnitMultC2->GetXaxis()->SetRangeUser(0,0.3);
    UnitMultC2->SetMinimum(0.97);
    UnitMultC2->SetMaximum(1.25);
    UnitMultC2->Draw();
    
    TH1D *UnitMultC2QS=(TH1D*)UnitMult[0][0][0]->Clone();
    UnitMultC2QS->Add(UnitMult[0][0][1], -(1-TwoFrac));
    UnitMultC2QS->Scale(1/TwoFrac);
    for(int bin=1; bin<=UnitMultC2QS->GetNbinsX(); bin++) {
      Float_t K2 = 1.0;
      if(FSICorrection) K2 = FSICorrelation(ch1_2, ch2_2, UnitMultC2QS->GetXaxis()->GetBinCenter(bin));
      UnitMultC2QS->SetBinContent(bin, UnitMultC2QS->GetBinContent(bin)/K2);
      UnitMultC2QS->SetBinError(bin, UnitMultC2QS->GetBinError(bin)/K2);
    }
    UnitMultC2QS->Divide(UnitMultC2den);
    UnitMultC2QS->SetMarkerColor(2); UnitMultC2QS->SetLineColor(2);
    UnitMultC2QS->Draw("same");
    }*/
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
  TProfile *K3Test=new TProfile("K3Test","",TERM1_3->GetNbinsX(), 0, TERM1_3->GetXaxis()->GetBinUpEdge(TERM1_3->GetNbinsX()),0,100.,"");
  K3Test->SetLineColor(2);
  for(int ii=0; ii<=K3avg[ch1_3][ch2_3][ch3_3][0]->GetNbinsX(); ii++){
    double newValue = K3avg[ch1_3][ch2_3][ch3_3][0]->GetBinContent(ii+1) - 1.0;
    newValue *= 1.2;// value used to modify K3
    newValue += 1.;
    double q3 = K3Test->GetXaxis()->GetBinCenter(ii+1);
    K3Test->Fill(q3, newValue);
  }
  
      
  TH1D *N3QS = (TH1D*)TERM1_3->Clone();
  N3QS->Add(TERM5_3, -f31);
  N3QS->Add(TERM2_3, -f32);
  N3QS->Add(TERM3_3, -f32);
  N3QS->Add(TERM4_3, -f32);
  N3QS->Scale(1/f33);
  
  N3QS->Multiply(K3avg[ch1_3][ch2_3][ch3_3][0]);
  //N3QS->Multiply(K3Test);
  
  
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
  TH1D *c3QS_partial = (TH1D*)N3QS->Clone();
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
    value += 2*TERM5_3->GetBinContent(bin);
    //
    n3QS->SetBinContent(bin, value);
    c3QS->SetBinContent(bin, value + TERM5_3->GetBinContent(bin));
    c3QS_partial->SetBinContent(bin, (TERM3_3->GetBinContent(bin) - (1-TwoFrac)*TERM5_3->GetBinContent(bin)) * K3avg[ch1_3][ch2_3][ch3_3][2]->GetBinContent(bin) / TwoFrac * MuonCorr3);
    //c3QS_partial->SetBinContent(bin, (TERM3_3->GetBinContent(bin)* MuonCorr3));
  }
  c3QS->Divide(TERM5_3);
  c3QS_partial->Divide(TERM5_3);
  
  C3QS->Draw();
  if(!MCcase) c3QS->Draw("same");
  //
  if(SameCharge && !MCcase) {
    if(ReNormBuiltBaseline){
      int BinL = Build_ThreeParticle[ch1_3]->GetXaxis()->FindBin(ReNormL_3);
      int BinH = Build_ThreeParticle[ch1_3]->GetXaxis()->FindBin(ReNormH_3);
      cout<<"C3 Build Norm = "<<C3QS->Integral(BinL,BinH) / Build_ThreeParticle[ch1_3]->Integral(BinL,BinH)<<endl;
      Build_ThreeParticle[ch1_3]->Scale(C3QS->Integral(BinL,BinH) / Build_ThreeParticle[ch1_3]->Integral(BinL,BinH));
      CumulantBuild_ThreeParticle[ch1_3]->Scale(c3QS->Integral(BinL,BinH) / CumulantBuild_ThreeParticle[ch1_3]->Integral(BinL,BinH));
    }
    Build_ThreeParticle[ch1_3]->Draw("same");
    CumulantBuild_ThreeParticle[ch1_3]->Draw("same");
  }
  
  //c3QS_partial->Draw();
  //cout<<Build_ThreeParticle[ch1_3]->GetBinContent(2)<<endl;
  //TH1D *C3raw = (TH1D*)TERM3_3->Clone();
  //C3raw->Divide(TERM5_3);
  //C3raw->SetMarkerColor(4);
  //C3raw->Draw();

  legend2->AddEntry(C3QS,"C_{3}^{QS}","p");
  if(!MCcase){
    legend2->AddEntry(c3QS,"c_{3}^{QS}{2-pion removal}","p");
    if(SameCharge) legend2->AddEntry(Build_ThreeParticle[ch1_3],"C_{3}^{QS} built","l");
    legend2->Draw("same");
  }
  //if(!SameCharge) Chi2_dist_c3->Draw("colz");

  //K3avg[ch1_3][ch2_3][ch3_3][0]->Draw();
  //K3Test->Draw("same");


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
 
  TH1D *Chi2NDF_PointSize_3 = new TH1D("Chi2NDF_PointSize_3","",100,-0.5,99.5);
  TH1D *Chi2NDF_FullSize_3 = new TH1D("Chi2NDF_FullSize_3","",100,-0.5,99.5);
  Chi2NDF_PointSize_3->SetLineColor(4); Chi2NDF_FullSize_3->SetLineColor(2);
  Chi2NDF_PointSize_3->SetMarkerColor(4); Chi2NDF_FullSize_3->SetMarkerColor(2);
  Chi2NDF_PointSize_3->GetXaxis()->SetTitle("Coherent fraction (%)"); Chi2NDF_PointSize_3->GetYaxis()->SetTitle("#sqrt{#chi^{2}}");
  
  if(SameCharge && CollisionType==0 && !MCcase){
    
    TH1D *tempDen = (TH1D*)Build_ThreeParticle_2D[ch1_3]->ProjectionY("Build3_Den", 4, 4);
    TH1D *tempDenNeg = (TH1D*)BuildNeg_ThreeParticle_2D[ch1_3]->ProjectionY("BuildNeg3_Den", 4, 4);
    tempDen->Add(tempDenNeg);// Add Pos and Neg Den

    for(int binG=5; binG<=179; binG++){// was 104
      TString *proName=new TString("Build3_");
      *proName += binG;
      TH1D *tempNum = (TH1D*)Build_ThreeParticle_2D[ch1_3]->ProjectionY(proName->Data(), binG, binG);
      proName->Append("_Neg");
      TH1D *tempNumNeg = (TH1D*)BuildNeg_ThreeParticle_2D[ch1_3]->ProjectionY(proName->Data(), binG, binG);
      //
      // Add Pos and Neg Num
      tempNum->Add(tempNumNeg);
      //
      tempNum->Divide(tempDen);
      
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
      if(binG<55) {Chi2NDF_PointSize_3->SetBinContent(1 + 2*(binG-5), sqrt(Chi2)/NDF); Chi2NDF_PointSize_3->SetBinError(1 + 2*(binG-5), 0.001);}
      else {Chi2NDF_FullSize_3->SetBinContent(1 + 2*(binG-55), sqrt(Chi2)/NDF); Chi2NDF_FullSize_3->SetBinError(1 + 2*(binG-55), 0.001);}
    }
    Chi2NDF_PointSize_3->SetMarkerStyle(20); Chi2NDF_FullSize_3->SetMarkerStyle(20);
   
    C3QSratio->Divide(Build_ThreeParticle[ch1_3]);
    //
    C3QSratio->SetMinimum(0.94); C3QSratio->SetMaximum(1.02);
    C3QSratio->GetYaxis()->SetTitleOffset(1.2);
    C3QSratio->GetYaxis()->SetTitle("C_{3}^{QS} / C_{3}(built)");
    C3QSratio->Draw();
    Unity->Draw("same");
    
  }
 
  // r3
  TH1D *r3;
  if(SameCharge && CollisionType==0 && !MCcase){
    r3 = (TH1D*)n3QS->Clone();
    TPN_ThreeParticle[ch1_3]->Multiply(MRC_SC_3[2]);
    r3->Divide(TPN_ThreeParticle[ch1_3]);
    //r3->Divide(CumulantBuild_ThreeParticle[ch1_3]);
    r3->GetXaxis()->SetRangeUser(0,0.08);
    r3->SetMinimum(0.5); r3->SetMaximum(2.5);
    r3->GetYaxis()->SetTitle("r_{3}");
    //
    //r3->Draw();
    //TPN_ThreeParticle[ch1_3]->Draw();
    //fTPNRejects3pion->SetLineColor(2);
    //fTPNRejects3pion->Draw("same");
  }
  
  cout.precision(4);
  // Print 3-pion points
  for(int bin=1; bin<=10; bin++){
    //cout<<C3QS->GetBinContent(bin)<<", ";
    //cout<<c3QS->GetBinContent(bin)<<", ";
    //cout<<C3QSratio->GetBinContent(bin)<<", ";
    //cout<<Build_ThreeParticle[ch1_3]->GetBinContent(bin)<<", ";
  }
  //cout<<endl;
  for(int bin=1; bin<=10; bin++){
    //cout<<C3QS->GetBinError(bin)<<", ";
    //cout<<c3QS->GetBinError(bin)<<", ";
  }
  //cout<<endl;

  double Chi2_K3resolve=0;
  //for(int bin=2; bin<=6; bin++) Chi2_K3resolve += pow((C3QSratio->GetBinContent(bin)-1) / (C3QSratio->GetBinError(bin)),2);
  for(int bin=2; bin<=4; bin++) {
    //double TotError = sqrt(pow(C3QSratio->GetBinError(bin),2) + pow(0.002,2));
    //double TotError = sqrt(pow(c3QS->GetBinError(bin),2) + pow(0.002,2));
    //Chi2_K3resolve += pow((C3QSratio->GetBinContent(bin)-1) / (TotError),2);
    //Chi2_K3resolve += pow((c3QS->GetBinContent(bin)-1) / (TotError),2);
  }  
  //cout<<"Chi2/NDF from K3 resolution = "<<Chi2_K3resolve/3.<<endl;
  

  ////////////////////////////////////////////////////////////////
  // 4-pion 
  cout<<"4-pion section"<<endl;
  
  TCanvas *can3 = new TCanvas("can3", "can3",11,600,700,800);// 60 was 600
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
  TPad *pad1 = new TPad("pad1","pad1",0.,Ystart,1.,1.);
  TPad *pad1_2 = new TPad("pad1_2","pad1_2",0.,0.,1.,Ystart);
  pad1->Draw();
  if(SameCharge==1) pad1_2->Draw();
  pad1->cd();
  gPad->SetRightMargin(0.02); gPad->SetTopMargin(0.02);
  TLegend *legend3 = (TLegend*)legend1->Clone();
  legend3->SetX1(0.47); legend3->SetX2(0.98);  legend3->SetY1(0.7); legend3->SetY2(0.98);

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
  // K4 scale variation
  TProfile *K4Test=new TProfile("K4Test","",TERMS_4[0]->GetNbinsX(), 0, TERMS_4[0]->GetXaxis()->GetBinUpEdge(TERMS_4[0]->GetNbinsX()),0,100.,"");
  K4Test->SetLineColor(2);
  for(int ii=0; ii<=K4avg[ch1_4][ch2_4][ch3_4][ch4_4][0]->GetNbinsX(); ii++){
    double newValue = K4avg[ch1_4][ch2_4][ch3_4][ch4_4][0]->GetBinContent(ii+1) - 1.0;
    newValue *= .89;// value used to increase K4
    newValue += 1;
    double q4 = K4Test->GetXaxis()->GetBinCenter(ii+1);
    K4Test->Fill(q4, newValue);
  }
  //N4QS->Multiply(K4Test);// K4 scale variation
  N4QS->Multiply(K4avg[ch1_4][ch2_4][ch3_4][ch4_4][0]);

  //K4avg[ch1_4][ch2_4][ch3_4][ch4_4][0]->Draw();
  //K4Test->Draw("same");
  
  //
  
  TH1D *C4QS=(TH1D*)N4QS->Clone();
  C4QS->GetXaxis()->SetRangeUser(0,0.16);
  C4QS->SetMinimum(0.5);
  C4QS->SetMarkerColor(4);
  
 
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
  TH1D *a4QS = (TH1D*)N4QS->Clone();
  TH1D *b4QS = (TH1D*)N4QS->Clone();
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
      N4QSvalue *= C4muonCorrectionSC[0]->GetBinContent(bin); 
      ErrorTerm[0] *= C4muonCorrectionSC[0]->GetBinContent(bin);
    }else if(MixedCharge4pionType==1) {
      N4QSvalue *= C4muonCorrectionMC1[0]->GetBinContent(bin);
      ErrorTerm[0] *= C4muonCorrectionMC1[0]->GetBinContent(bin);
    }else {
      N4QSvalue *= C4muonCorrectionMC2[0]->GetBinContent(bin); 
      ErrorTerm[0] *= C4muonCorrectionMC2[0]->GetBinContent(bin);
    }
    for(int ii=0; ii<4; ii++) {
      FinalValue[ii] = N4QSvalue; 
      FinalValue_e[ii] = pow(ErrorTerm[0],2);
    }
    
    if(SameCharge) {
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
	FinalValue[0] -= SubtractionTerm[0] * C4muonCorrectionMC1[2]->GetBinContent(bin) - TERMS_4[12]->GetBinContent(bin);// [2] is +++ MuonCorr
	FinalValue[0] -= 3*(SubtractionTerm[6] * C4muonCorrectionMC1[3]->GetBinContent(bin) - TERMS_4[12]->GetBinContent(bin));
	//
	FinalValue_e[0] += pow(ErrorTerm[1]*C4muonCorrectionMC1[2]->GetBinContent(bin),2);
	FinalValue[1] -= 3*SubtractionTerm[4] - 3*TERMS_4[12]->GetBinContent(bin);
	FinalValue_e[1] += 3*pow(ErrorTerm[5],2) + 3*TERMS_4[12]->GetBinContent(bin);
      }else if(MixedCharge4pionType==1 && CHARGE==+1){
	FinalValue[0] -= SubtractionTerm[3] * C4muonCorrectionMC1[2]->GetBinContent(bin) - TERMS_4[12]->GetBinContent(bin);
	FinalValue[0] -= 3*(SubtractionTerm[6] * C4muonCorrectionMC1[3]->GetBinContent(bin) - TERMS_4[12]->GetBinContent(bin));
	//
	FinalValue_e[0] += pow(ErrorTerm[4]*C4muonCorrectionMC1[2]->GetBinContent(bin),2);
	FinalValue[1] -= 3*SubtractionTerm[9] - 3*TERMS_4[12]->GetBinContent(bin);
	FinalValue_e[1] += 3*pow(ErrorTerm[10],2) + 3*TERMS_4[12]->GetBinContent(bin);
      }else{// --++ case
	FinalValue[0] -= 2*(SubtractionTerm[4] * C4muonCorrectionMC2[2]->GetBinContent(bin) - TERMS_4[12]->GetBinContent(bin));
	FinalValue[0] -= (SubtractionTerm[10] - TERMS_4[12]->GetBinContent(bin));
	FinalValue[0] -= 4*(SubtractionTerm[6] * C4muonCorrectionMC2[3]->GetBinContent(bin) - TERMS_4[12]->GetBinContent(bin));
	FinalValue_e[0] += 2*pow(ErrorTerm[5],2) + pow(ErrorTerm[11],2) + 2*TERMS_4[12]->GetBinContent(bin);
	//
	FinalValue[1] -= 2*(SubtractionTerm[4] * C4muonCorrectionMC2[2]->GetBinContent(bin) - TERMS_4[12]->GetBinContent(bin));
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
    a4QS->SetBinContent(bin, FinalValue[1]);
    a4QS->SetBinError(bin, FinalValue_e[1]);
    //
    b4QS->SetBinContent(bin, FinalValue[2]);
    b4QS->SetBinError(bin, FinalValue_e[2]);
    //
    c4QS_RemovalStage3->SetBinContent(bin, FinalValue[3]);
    c4QS_RemovalStage3->SetBinError(bin, FinalValue_e[3]);
    C4_pieces[5]->SetBinContent(bin, SubtractionTerm[6] / TERMS_4[12]->GetBinContent(bin));
  }

  C4QS->Divide(TERMS_4[12]);
  c4QS->Divide(TERMS_4[12]);
  c4QS->SetMarkerColor(1); c4QS->SetLineColor(1);
  //
  a4QS->Divide(TERMS_4[12]);
  a4QS->SetMarkerColor(2); a4QS->SetLineColor(2); 
  b4QS->Divide(TERMS_4[12]);
  b4QS->SetMarkerColor(6); b4QS->SetLineColor(6); 
  c4QS_RemovalStage3->Divide(TERMS_4[12]);
  c4QS_RemovalStage3->SetMarkerColor(7); c4QS_RemovalStage3->SetLineColor(7); 
  //
  //
  
  
  
  if(SameCharge) {
    if(CollisionType==0) C4QS->SetMaximum(7.5);
    else C4QS->SetMaximum(12);
  }else if(MixedCharge4pionType==1) C4QS->SetMaximum(4);
  else C4QS->SetMaximum(3);
  
  
  if(CollisionType!=0) C4QS->GetXaxis()->SetRangeUser(0,0.5);
  C4QS->Draw();
  if(!MCcase) a4QS->Draw("same");
  if(SameCharge && !MCcase) b4QS->Draw("same");
  if(!MCcase) c4QS->Draw("same");
  
 
  
  if(SameCharge && !MCcase) {
    //BuildFromFits_M->Draw("E2 same");
    //
    // renormalize
    if(ReNormBuiltBaseline){
      int BinL = Build_FourParticle[ch1_4]->GetXaxis()->FindBin(ReNormL_4);
      int BinH = Build_FourParticle[ch1_4]->GetXaxis()->FindBin(ReNormH_4);
      //if(CollisionType==0){
	cout<<"C2 Build ReNorm = "<<C4QS->Integral(BinL,BinH) / Build_FourParticle[ch1_4]->Integral(BinL,BinH)<<endl;
	Build_FourParticle[ch1_4]->Scale(C4QS->Integral(BinL,BinH) / Build_FourParticle[ch1_4]->Integral(BinL,BinH));
	PrimeBuild_FourParticle[ch1_4]->Scale(a4QS->Integral(BinL,BinH) / PrimeBuild_FourParticle[ch1_4]->Integral(BinL,BinH));
	PrimePrimeBuild_FourParticle[ch1_4]->Scale(b4QS->Integral(BinL,BinH) / PrimePrimeBuild_FourParticle[ch1_4]->Integral(BinL,BinH));
	CumulantBuild_FourParticle[ch1_4]->Scale(c4QS->Integral(BinL,BinH) / CumulantBuild_FourParticle[ch1_4]->Integral(BinL,BinH));
	//}
      //
      cout<<"c3 Fit Build ReNorm = "<<C4QS->Integral(BinL,BinH) / BuildFromFits1->Integral(BinL,BinH)<<endl;
      cout<<"C3 Fit Build ReNorm = "<<C4QS->Integral(BinL,BinH) / BuildFromFits2->Integral(BinL,BinH)<<endl;
      BuildFromFits1->Scale(C4QS->Integral(BinL,BinH) / BuildFromFits1->Integral(BinL,BinH));
      PrimeBuildFromFits1->Scale(a4QS->Integral(BinL,BinH) / PrimeBuildFromFits1->Integral(BinL,BinH));
      PrimePrimeBuildFromFits1->Scale(b4QS->Integral(BinL,BinH) / PrimePrimeBuildFromFits1->Integral(BinL,BinH));
      CumulantBuildFromFits1->Scale(c4QS->Integral(BinL,BinH) / CumulantBuildFromFits1->Integral(BinL,BinH));
      BuildFromFits2->Scale(C4QS->Integral(BinL,BinH) / BuildFromFits2->Integral(BinL,BinH));
      PrimeBuildFromFits2->Scale(a4QS->Integral(BinL,BinH) / PrimeBuildFromFits2->Integral(BinL,BinH));
      PrimePrimeBuildFromFits2->Scale(b4QS->Integral(BinL,BinH) / PrimePrimeBuildFromFits2->Integral(BinL,BinH));
      CumulantBuildFromFits2->Scale(c4QS->Integral(BinL,BinH) / CumulantBuildFromFits2->Integral(BinL,BinH));
    }
    //if(CollisionType==0){
    if(!FitBuild){
      Build_FourParticle[ch1_4]->Draw("same");
      PrimeBuild_FourParticle[ch1_4]->Draw("same");
      PrimePrimeBuild_FourParticle[ch1_4]->Draw("same");
      CumulantBuild_FourParticle[ch1_4]->Draw("same");
    }else{
      BuildFromFits1->Draw("same");
      PrimeBuildFromFits1->Draw("same");
      PrimePrimeBuildFromFits1->Draw("same");
      CumulantBuildFromFits1->Draw("same");
      BuildFromFits2->Draw("same");
      PrimeBuildFromFits2->Draw("same");
      PrimePrimeBuildFromFits2->Draw("same");
      CumulantBuildFromFits2->Draw("same");
    }
  }
  
  legend3->AddEntry(C4QS,"C_{4}^{QS}","p");
  
  if(!MCcase){
    legend3->AddEntry(a4QS,"a_{4}^{QS}","p");
    if(SameCharge) legend3->AddEntry(b4QS,"b_{4}^{QS}","p");
    if(SameCharge) legend3->AddEntry(c4QS,"c_{4}^{QS}","p");
    if(!SameCharge && MixedCharge4pionType==1) legend3->AddEntry(c4QS,"c_{4}^{QS}","p");
    if(!SameCharge && MixedCharge4pionType==2) legend3->AddEntry(c4QS,"c_{4}^{QS}","p");
    
    if(SameCharge && !FitBuild) legend3->AddEntry(Build_FourParticle[ch1_4],"Built from C_{2}","l");
    if(SameCharge && FitBuild) legend3->AddEntry(BuildFromFits1,"C_{4}^{QS} built from c_{3}","l");
    if(SameCharge && FitBuild) legend3->AddEntry(BuildFromFits2,"C_{4}^{QS} built from C_{3}","l");
    legend3->Draw("same");
  }
  Unity->Draw("same");

 

  //C4_pieces[5]->Draw();

  //int BOI=4;
  //double Perc=1.07;
  //double BuildMod = 1 + Perc*(Build_FourParticle[ch1_4]->GetBinContent(BOI)-PrimeBuild_FourParticle[ch1_4]->GetBinContent(BOI));
  //BuildMod += pow(Perc,2)*(PrimeBuild_FourParticle[ch1_4]->GetBinContent(BOI)-PrimePrimeBuild_FourParticle[ch1_4]->GetBinContent(BOI));
  //BuildMod += pow(Perc,1.5)*(PrimePrimeBuild_FourParticle[ch1_4]->GetBinContent(BOI)-CumulantBuild_FourParticle[ch1_4]->GetBinContent(BOI));
  //BuildMod += pow(Perc,2)*(CumulantBuild_FourParticle[ch1_4]->GetBinContent(BOI)-1);
  //cout<<" BuildMod ratio = "<<Build_FourParticle[ch1_4]->GetBinContent(BOI) / BuildMod<<endl;
  

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
  
  TH1D* C4QSratio = (TH1D*)C4QS->Clone();
  if(SameCharge && !FitBuild){
     
    pad1_2->cd();
    /*TCanvas *can4 = new TCanvas("can4", "can4",600,600,700,500);
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
    can4->SetFrameBorderMode(0);*/

    gPad->SetRightMargin(0.02); gPad->SetTopMargin(0.02);
    gPad->SetBottomMargin(0.15);
    TLegend *legend3_2 = (TLegend*)legend1->Clone();
    legend3_2->SetX1(0.45); legend3_2->SetX2(0.98);  legend3_2->SetY1(0.6); legend3_2->SetY2(0.95);
    //
    TH1D *Chi2NDF_PointSize_4 = new TH1D("Chi2NDF_PointSize_4","",100,-0.5,99.5);
    TH1D *Chi2NDF_FullSize_4 = new TH1D("Chi2NDF_FullSize_4","",100,-0.5,99.5);
    Chi2NDF_PointSize_4->SetLineColor(4); Chi2NDF_FullSize_4->SetLineColor(2);
    Chi2NDF_PointSize_4->SetMarkerColor(4); Chi2NDF_FullSize_4->SetMarkerColor(2);
    Chi2NDF_PointSize_4->GetXaxis()->SetTitle("Coherent fraction (%)"); Chi2NDF_PointSize_4->GetYaxis()->SetTitle("#sqrt{#chi^{2}}");
    
    C4QSratio->Divide(Build_FourParticle[ch1_4]);
    C4QSratio->GetYaxis()->SetTitle(" C_{4} Measured / Built");
    C4QSratio->SetMinimum(0.86); C4QSratio->SetMaximum(1.04);
    C4QSratio->GetXaxis()->SetTitleSize(2.*C4QSratio->GetXaxis()->GetTitleSize());
    C4QSratio->GetYaxis()->SetTitleSize(2.*C4QSratio->GetYaxis()->GetTitleSize());
    C4QSratio->GetYaxis()->SetTitleOffset(0.7); C4QSratio->GetYaxis()->SetNdivisions(404);
    C4QSratio->GetXaxis()->SetLabelSize(2.*C4QSratio->GetXaxis()->GetLabelSize());
    C4QSratio->GetYaxis()->SetLabelSize(2.*C4QSratio->GetYaxis()->GetLabelSize());
    C4QSratio->Draw();
    Unity->Draw("same");
     
    
  }

  // Print 4-pion points
  int MaxPrintBins=12;
  if(CollisionType!=0) MaxPrintBins=17;
  for(int bin=1; bin<=MaxPrintBins; bin++){
    //cout<<C4QS->GetBinContent(bin)<<", ";
    //cout<<c4QS->GetBinContent(bin)<<", ";
    //cout<<Build_FourParticle[ch1_4]->GetBinContent(bin)<<", ";
    cout<<C4QSratio->GetBinContent(bin)<<", ";
    //cout<<C4raw->GetBinContent(bin)<<", ";
    //cout<<K4avg[ch1_4][ch2_4][ch3_4][ch4_4][0]->GetBinContent(bin)<<", ";
  }
  cout<<endl;
  for(int bin=1; bin<=MaxPrintBins; bin++){
    //cout<<c4QS->GetBinContent(bin)<<", ";
    //cout<<C4QS->GetBinError(bin)<<", ";
    //cout<<c4QS->GetBinError(bin)<<", ";
    //cout<<C4raw->GetBinError(bin)<<", ";
    //cout<<C4QSratio->GetBinError(bin)<<", ";
  }
  //cout<<endl;

  double Chi2_K4resolve=0;
  for(int bin=3; bin<=5; bin++) {// was 3 to 7
    //double TotError= sqrt(pow(C4QSratio->GetBinError(bin),2) + pow(0.001,2));//approximate addition of systematic error difference between low and high KT4 (same-charge) 
    //double TotError= sqrt(pow(c4QS->GetBinError(bin),2) + pow(0.001,2));//approximate addition of systematic error (mixed-charge)
    //
    //Chi2_K4resolve += pow((C4QSratio->GetBinContent(bin)-1) / (TotError),2);
    //Chi2_K4resolve += pow((c4QS->GetBinContent(bin)-1) / (TotError),2);
  }
  //cout<<"Chi2/NDF from K4 resolution = "<<Chi2_K4resolve/3.<<endl;

  ////////////////////////////////////////////////////////////////
  // r4

  
  TCanvas *can5 = new TCanvas("can5", "can5",1200,53,700,700);
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
  gPad->SetRightMargin(0.02); gPad->SetTopMargin(0.02);
 
  //K4avg[ch1_4][ch2_4][ch3_4][ch4_4][0]->Draw();
  //K4Test->Draw("same");

  TF1 *ChaoticLimit_r42 = new TF1("ChaoticLimit_r42","6",0,10);
  ChaoticLimit_r42->SetLineStyle(2);
  TH1D *r42;
  if(SameCharge && CollisionType==0 && !MCcase){
    
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
    savefilename->Append("_CT");
    *savefilename += CollisionType;
    //
    if(SameCharge) savefilename->Append("_SC");
    else savefilename->Append("_MC");
    //
    if(!SameCharge) *savefilename += MixedCharge4pionType_def;
    //
    if(CHARGE==1) savefilename->Append("_Pos_");
    else savefilename->Append("_Neg_");
    //
    
    savefilename->Append("ED");
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
    a4QS->Write("a4QS");
    if(SameCharge) {
      b4QS->Write("b4QS");
      if(!MCcase){
	if(CollisionType==0){
	  r3->Write("r3");
	  r42->Write("r42");
	}
	Build_ThreeParticle[ch1_3]->Write("C3QS_built");
	CumulantBuild_ThreeParticle[ch1_3]->Write("c3QS_built");
	Build_FourParticle[ch1_4]->Write("C4QS_built");
	PrimeBuild_FourParticle[ch1_4]->Write("a4QS_built");
	PrimePrimeBuild_FourParticle[ch1_4]->Write("b4QS_built");
	CumulantBuild_FourParticle[ch1_4]->Write("c4QS_built");
	//
	Build_ThreeParticle_2D[ch1_3]->Write("C3QS_built2D");
	CumulantBuild_ThreeParticle_2D[ch1_3]->Write("c3QS_built2D");
	Build_FourParticle_2D[ch1_4]->Write("C4QS_built2D");
	PrimeBuild_FourParticle_2D[ch1_4]->Write("a4QS_built2D");
	PrimePrimeBuild_FourParticle_2D[ch1_4]->Write("b4QS_built2D");
	CumulantBuild_FourParticle_2D[ch1_4]->Write("c4QS_built2D");
	BuildNeg_ThreeParticle_2D[ch1_3]->Write("C3QS_Negbuilt2D");
	CumulantBuildNeg_ThreeParticle_2D[ch1_3]->Write("c3QS_Negbuilt2D");
	BuildNeg_FourParticle_2D[ch1_4]->Write("C4QS_Negbuilt2D");
	PrimeBuildNeg_FourParticle_2D[ch1_4]->Write("a4QS_Negbuilt2D");
	PrimePrimeBuildNeg_FourParticle_2D[ch1_4]->Write("b4QS_Negbuilt2D");
	CumulantBuildNeg_FourParticle_2D[ch1_4]->Write("c4QS_Negbuilt2D");
      }
      //
      BuildFromFits_3D->Write("C4QS_BuiltFromFits3D");
      PrimeBuildFromFits_3D->Write("a4QS_BuiltFromFits3D");
      PrimePrimeBuildFromFits_3D->Write("b4QS_BuiltFromFits3D");
      CumulantBuildFromFits_3D->Write("c4QS_BuiltFromFits3D");
      BuildFromFits1->Write("C4QS_BuiltFromFits1");
      PrimeBuildFromFits1->Write("a4QS_BuiltFromFits1");
      PrimePrimeBuildFromFits1->Write("b4QS_BuiltFromFits1");
      CumulantBuildFromFits1->Write("c4QS_BuiltFromFits1");
      BuildFromFits2->Write("C4QS_BuiltFromFits2");
      PrimeBuildFromFits2->Write("a4QS_BuiltFromFits2");
      PrimePrimeBuildFromFits2->Write("b4QS_BuiltFromFits2");
      CumulantBuildFromFits2->Write("c4QS_BuiltFromFits2");
      //
      BuildFromFits_M->Write("C4QS_BuiltFromFits_M");
      PrimeBuildFromFits_M->Write("a4QS_BuiltFromFits_M");
      PrimePrimeBuildFromFits_M->Write("b4QS_BuiltFromFits_M");
      CumulantBuildFromFits_M->Write("c4QS_BuiltFromFits_M");
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
  if(CollisionType_def!=0) momresfilename->Append("_ppAndpPb");
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
  if(CollisionType_def!=0) name->Append("_ppAndpPb");
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
    if(CollisionType_def==0){
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
