#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <Riostream.h>
#include <complex>

#include "TObject.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TVector2.h"
#include "TVector3.h"
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
#include "TChain.h"
#include "TLorentzVector.h"

#define PI 3.1415926
#define BohrR 1963.6885 // Mate's value(1963.6885) ~ 387.5/0.197327(1963.7455)
#define FmToGeV 0.197327 // 0.197327
#define fzero_aa 0.942597 // 0.186fm/FmToGeV (scattering length of pi+pi-) = 0.942597
#define fzero_ab -0.89192 // -0.176fm/FmtoGeV (scattering length (pi+pi-)-->(pi0pi0) = -0.89192
#define fzero_bb 0.3119 // fzero_aa + 1/sqrt(2)*fzero_ab = 0.0615/FmToGeV = 0.3119
#define dzero -50.6773 // -10/FmToGeV = -50.6773 (effective range)
#define EulerC 0.5772156649 // Euler constant
#define masspi0 0.1349766 // pi0 mass
#define masspiC 0.1395702 // pi+ mass
#define massKch 0.493677 // K+- mass
#define massK0 0.497614 // K0 mass

#define SF 1.0 // Scale Factor for spatial coordinates
#define SF_Mom 1.0 // Scale Factor to reduce momenta
#define RstarMax 405.4/SF // 405.4 (80 fm / FmToGeV), tried 253.4 (50 fm / FmToGeV) as a variation
#define RstarMin 0.507/SF // 0.507 (0.1 fm / FmToGeV)
using namespace std;

// g and p used in Gamma function
const int g_gamma = 7;
double p_gamma[9] = {0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};
// h coefficients for very small kstar (eta > 0.2)
double hcoef[15]={.8333333333e-1, .8333333333e-2, .396825396825e-2, .4166666667e-2, .7575757576e-2, .2109279609e-1, .8333333333e-1, .4432598039e0 , .305395433e1, .2645621212e2, .2814601449e3, .3607510546e4, .5482758333e5, .9749368235e6, .200526958e8};
// Scattering length (in GeV) coefficients in isospin representation from Lednicky's code from Colangelo (2001)
double fzero[2][6]={{.220, .268,-.0139,-.00139, 36.77, 15*0}, {-.0444,-.0857,-.00221,-.000129,-21.62, 15*0}};

double massPOI;
int PIDPOI;
double qCutOff;

float Gamov(float);
float fact(float);
double hFunction(double);
double KinverseFunction(double, int);
void Shuffle(int*, int, int);
void BoostPRF(float [4], float [4], float [4], TVector3*);


complex<double> gamma(complex<double> z)
{// Lanczos approximation
   
  if(floor(abs(z)) == z)
    return fact(abs(z)-1.0);
      
  if(real(z)<(double)(0.5)) return PI/(sin(PI*z)*gamma((1.0-z)));
  z -= 1.0;
  complex<double> x = p_gamma[0];
  for(int i=1;i<g_gamma+2;i++) x = x+p_gamma[i]/(z+(double)i);
  complex<double> t = z+(double)g_gamma+0.5;
  return sqrt(2*PI)*pow(t,z+0.5)*exp(-t)*x;
}

// confluent hypergeometric function
complex<double> conf_hyp(complex<double> a, complex<double> b, complex<double> z)
{// b is assumed to be 1.0 here
  complex<double> S(1.0, 0.0);
  complex<double> F=S;
  complex<double> alf1 = a - 1.0;
  double J=0;
  double errorCutOff=0.00001;// series truncation point (0.00001, 1e-6 for strict case)
  if(abs(z) > 40.) {// high k*r* case (was 20., now 40. Improved accuracy for qinv>60 MeV)
    double eta = -a.imag();
    double RKS = z.imag();
    double D1 = log(RKS)*eta;
    double D2 = pow(eta,2)/RKS;
    complex<double> arg(1.0, eta);
    complex<double> EIDC = gamma(arg);
    EIDC /= abs(EIDC);
    complex<double> ZZ1(cos(D1), sin(D1));
    ZZ1 /= EIDC;
    complex<double> value(1.0, -D2);
    value *= ZZ1;
    complex<double> value2(cos(RKS), sin(RKS));
    F = value - value2*eta/RKS/ZZ1;
    return F;
  }else {// low k*r* case
    while(fabs(S.real())+fabs(S.imag()) > errorCutOff){
      J++;
      complex<double> A = (alf1 + J)/(pow(J,2));
      S *= A*z;
      F += S;
    }
    return F;
  }
}

void Therm(){

  TRandom3 *RandomNumber = new TRandom3();
  RandomNumber->SetSeed(0);
  
  
  bool StrongFSI=kTRUE;
  bool LambdaDirect=kFALSE;
  bool RemoveXP=kFALSE;
  bool Gaussian=kFALSE;
  bool ThreeParticle=kFALSE;
  bool BoostPairsIntoPRF=kTRUE;
  massPOI = masspiC;
  PIDPOI = 211;// 211 for pi+-
  qCutOff = 0.1;// 0.1 or 0.4 for pp
  int qBins = int(qCutOff/0.002);
  //
  int Nfiles=1000;// 1000 for b2, 500 b3-b5, 1500 b7-b8, 3000 b9
  int bValue=2;

  double RValue=8;// Gaussian radius
  double NsigmaGauss=2;// number of sigma to integrate over for a Gaussian source
  if(bValue==2) RValue=8;//8
  if(bValue==3) RValue=8;//8
  if(bValue==5) RValue=7;//7
  if(bValue==7) RValue=6;//6
  if(bValue==8) RValue=5;//5
  if(bValue==9) RValue=5;//5
  TF1 *SingleParticleGaussian = new TF1("SingleParticleGaussian","exp(-pow(x/(sqrt(2.)*[0]),2))",0,1000);
  SingleParticleGaussian->SetParameter(0,RValue);

  TFile *outfile = new TFile("Therm_dist_temp.root","RECREATE");
  TH1D *r_dist=new TH1D("r_dist","",400,0,100);
  TH1D *x_dist=new TH1D("x_dist","",400,0,100);
  TH1D *Pt_dist=new TH1D("Pt_dist","",200,0,2);
  TH1D *P_dist=new TH1D("P_dist","",900,0.1,1);
  TH1D *P_dist_lowq=new TH1D("P_dist_lowq","",900,0.1,1);
  TProfile *AvgMult = new TProfile("AvgMult","",1,0.5,1.5, 0,2000,"");
  //
  TH2D *rstarPRF_dist_ss=new TH2D("rstarPRF_dist_ss","",20,0,1.0, 400,0,100);
  TH2D *tstar_dist_ss=new TH2D("tstar_dist_ss","",20,0,1.0, 400,0,100);
  //
  TH2D *rstarPRF_dist_os=new TH2D("rstarPRF_dist_os","",20,0,1.0, 400,0,100);
  TH2D *tstar_dist_os=new TH2D("tstar_dist_os","",20,0,1.0, 400,0,100);

  TH2D *Num_Cos_ss=new TH2D("Num_Cos_ss","",20,0,1, 40,0,0.2);
  TH2D *Num_Cos_os=new TH2D("Num_Cos_os","",20,0,1, 40,0,0.2);
  TH2D *Num_CosFSI_ss=new TH2D("Num_CosFSI_ss","",20,0,1, 40,0,0.2);
  TH2D *Num_CosFSI_os=new TH2D("Num_CosFSI_os","",20,0,1, 40,0,0.2);
  TH2D *NumSq_CosFSI_ss=new TH2D("NumSq_CosFSI_ss","",20,0,1, 40,0,0.2);
  TH2D *NumSq_CosFSI_os=new TH2D("NumSq_CosFSI_os","",20,0,1, 40,0,0.2);
  TH2D *Num_PrimCosFSI_ss=new TH2D("Num_PrimCosFSI_ss","",20,0,1, 40,0,0.2);
  TH2D *Num_PrimCosFSI_os=new TH2D("Num_PrimCosFSI_os","",20,0,1, 40,0,0.2);
  //
  TH2D *Den_ss=new TH2D("Den_ss","",20,0,1, 40,0,0.2);
  TH2D *Den_os=new TH2D("Den_os","",20,0,1, 40,0,0.2);
  TH2D *DenEM_ss=new TH2D("DenEM_ss","",20,0,1, 40,0,0.2);
  TH2D *DenEM_os=new TH2D("DenEM_os","",20,0,1, 40,0,0.2);
  TH2D *LargeRpairs_ss=new TH2D("LargeRpairs_ss","",20,0,1, 40,0,0.2);
  TH2D *LargeRpairs_os=new TH2D("LargeRpairs_os","",20,0,1, 40,0,0.2);
  //
  TH2D *NumLambda1=new TH2D("NumLambda1","",40,0,0.2, 100,0,100);
  TH1D *DenLambda1=new TH1D("DenLambda1","",40,0,0.2);
  TH3D *NumLambda2_ss=new TH3D("NumLambda2_ss","",20,0,1, 40,0,0.2, 100,0,100);
  TH3D *NumLambda2_os=new TH3D("NumLambda2_os","",20,0,1, 40,0,0.2, 100,0,100);
  TH2D *DenLambda2_ss=new TH2D("DenLambda2_ss","",20,0,1, 40,0,0.2);
  TH2D *DenLambda2_os=new TH2D("DenLambda2_os","",20,0,1, 40,0,0.2);
  TH3D *NumLambda3_ss=new TH3D("NumLambda3_ss","",40,0,.2, 40,0,.2, 40,0,.2);
  TH3D *DenLambda3_ss=new TH3D("DenLambda3_ss","",40,0,.2, 40,0,.2, 40,0,.2);
  TH3D *NumLambda3_os=new TH3D("NumLambda3_os","",40,0,.2, 40,0,.2, 40,0,.2);
  TH3D *DenLambda3_os=new TH3D("DenLambda3_os","",40,0,.2, 40,0,.2, 40,0,.2);
  //
  TH3D *rstar3_dist_ss=new TH3D("rstar3_dist_ss","",200,0,100, 200,0,100, 200,0,100);
  TH3D *tstar3_dist_ss=new TH3D("tstar3_dist_ss","",200,0,100, 200,0,100, 200,0,100);
  //
  TH3D *r3num=new TH3D("r3num","",qBins,0,qCutOff, qBins,0,qCutOff, qBins,0,qCutOff);
  TH3D *r3den1=new TH3D("r3den1","",qBins,0,qCutOff, qBins,0,qCutOff, qBins,0,qCutOff);
  TH3D *r3den2=new TH3D("r3den2","",qBins,0,qCutOff, qBins,0,qCutOff, qBins,0,qCutOff);
  TH3D *r3den3=new TH3D("r3den3","",qBins,0,qCutOff, qBins,0,qCutOff, qBins,0,qCutOff);
  TH3D *r3numSq=new TH3D("r3numSq","",qBins,0,qCutOff, qBins,0,qCutOff, qBins,0,qCutOff);
  TH3D *r3den1Sq=new TH3D("r3den1Sq","",qBins,0,qCutOff, qBins,0,qCutOff, qBins,0,qCutOff);
  TH3D *r3den2Sq=new TH3D("r3den2Sq","",qBins,0,qCutOff, qBins,0,qCutOff, qBins,0,qCutOff);
  TH3D *r3den3Sq=new TH3D("r3den3Sq","",qBins,0,qCutOff, qBins,0,qCutOff, qBins,0,qCutOff);
  TH3D *r3numEn=new TH3D("r3numEn","",qBins,0,qCutOff, qBins,0,qCutOff, qBins,0,qCutOff);
  TH3D *r3den1En=new TH3D("r3den1En","",qBins,0,qCutOff, qBins,0,qCutOff, qBins,0,qCutOff);
  TH3D *r3den2En=new TH3D("r3den2En","",qBins,0,qCutOff, qBins,0,qCutOff, qBins,0,qCutOff);
  TH3D *r3den3En=new TH3D("r3den3En","",qBins,0,qCutOff, qBins,0,qCutOff, qBins,0,qCutOff);
  TH3D *r3numME=new TH3D("r3numME","",qBins,0,qCutOff, qBins,0,qCutOff, qBins,0,qCutOff);
  TH3D *r3den1ME=new TH3D("r3den1ME","",qBins,0,qCutOff, qBins,0,qCutOff, qBins,0,qCutOff);
  TH3D *r3den2ME=new TH3D("r3den2ME","",qBins,0,qCutOff, qBins,0,qCutOff, qBins,0,qCutOff);
  TH3D *r3den3ME=new TH3D("r3den3ME","",qBins,0,qCutOff, qBins,0,qCutOff, qBins,0,qCutOff);
  //
  TH1D *Num_Cos12=new TH1D("Num_Cos12","",40,0,0.2);
  TH1D *Num_Cos23=new TH1D("Num_Cos23","",40,0,0.2);
  TH1D *Num_Cos31=new TH1D("Num_Cos31","",40,0,0.2);
  TH1D *Den_Cos12=new TH1D("Den_Cos12","",40,0,0.2);
  TH1D *Den_Cos23=new TH1D("Den_Cos23","",40,0,0.2);
  TH1D *Den_Cos31=new TH1D("Den_Cos31","",40,0,0.2);
  //
  //
  TH2D *K2_ss = new TH2D("K2_ss","",20,0,1, qBins,0,qCutOff);
  TH2D *PlaneWF_ss = new TH2D("PlaneWF_ss","",20,0,1, qBins,0,qCutOff);
  TH2D *K2_os = new TH2D("K2_os","",20,0,1, qBins,0,qCutOff);
  TH2D *PlaneWF_os = new TH2D("PlaneWF_os","",20,0,1, qBins,0,qCutOff);
  //
  TH1D *PlaneWF3ss = new TH1D("PlaneWF3ss","",qBins,0,qCutOff);
  TH1D *PlaneWF3os = new TH1D("PlaneWF3os","",qBins,0,qCutOff);
  TH1D *K3ss = new TH1D("K3ss","",qBins,0,qCutOff);
  K3ss->GetXaxis()->SetTitle("Q3 (GeV/c)");
  K3ss->GetYaxis()->SetTitle("K_{3}");
  TH1D *K3os = new TH1D("K3os","",qBins,0,qCutOff);
  K3os->GetXaxis()->SetTitle("Q3 (GeV/c)");
  K3os->GetYaxis()->SetTitle("K_{3}");
  TH3D *K3ss_3D = new TH3D("K3ss_3D","", qBins,0,qCutOff, qBins,0,qCutOff, qBins,0,qCutOff);
  K3ss_3D->GetXaxis()->SetTitle("q12");
  K3ss_3D->GetYaxis()->SetTitle("q23");
  K3ss_3D->GetZaxis()->SetTitle("q31");
  TH3D *PlaneWF3ss_3D = new TH3D("PlaneWF3ss_3D","", qBins,0,qCutOff, qBins,0,qCutOff, qBins,0,qCutOff);
  TH3D *K3os_3D = new TH3D("K3os_3D","", qBins,0,qCutOff, qBins,0,qCutOff, qBins,0,qCutOff);
  K3os_3D->GetXaxis()->SetTitle("q12");
  K3os_3D->GetYaxis()->SetTitle("q23");
  K3os_3D->GetZaxis()->SetTitle("q31");
  TH3D *PlaneWF3os_3D = new TH3D("PlaneWF3os_3D","", qBins,0,qCutOff, qBins,0,qCutOff, qBins,0,qCutOff);

  TH1D *test=new TH1D("test","",200,-2*PI,2*PI);

  ////////////////////////
 
  
  complex<double> binput(1.0,0.0);
  complex<double> CHG(0.,0.);// for Strong FSI in 3-body WF
  double kstar=0;
  double rstar=0;
  double kstar_1=0, kstar_2=0;
  double rstar_1=0, rstar_2=0;
  double phase_1=0, phase_2=0, phase_3=0, phase_4=0;
  double qinvMixedC_ss=0, qinvMixedC_os1=0, qinvMixedC_os2=0;
  double WeightWF_aa=1.0, WeightWF_bb=1.0;
  complex<double> CHG_1(0.,0.);
  complex<double> CHG_2(0.,0.);
  complex<double> CHG_3(0.,0.);
  complex<double> CHG_4(0.,0.);
  complex<double> AlphaWF[4];
  complex<double> BetaWF[4];
  double ThreePhase_1=0, ThreePhase_2=0, ThreePhase_3=0;
  double ThreePhase_4=0, ThreePhase_5=0, ThreePhase_6=0;
  double TwoPhase_12=0, TwoPhase_23=0, TwoPhase_31=0;

  const int MAXPIONS = 2000;
  const int NevtBuff = 1;
  float P[NevtBuff][MAXPIONS][4]={{{0}}};
  float X[NevtBuff][MAXPIONS][4]={{{0}}};
  int PID[NevtBuff][MAXPIONS]={{0}};
  bool Primary[NevtBuff][MAXPIONS]={{0}};
  int Nparticles[NevtBuff]={0};
  Int_t ShuffledIndex[MAXPIONS]={0};// x-p decorrelation
  TVector3 P1_PRF;
  TVector3 P2_PRF;
  TVector3 P3_PRF;
  TVector3 X1_PRF;
  TVector3 X2_PRF;
  TVector3 X3_PRF;
  TVector3 X1_RF;
  TVector3 X2_RF;
  float P_dummy[4]={0};
	    


  // First Create B(rho,eta) and P(rho,eta) and h(eta) functions (for G hypergeometric functions (strong+Coulomb FSI))
  const int Btermslimit=500;// 500
  double rhoLimit = .05*100/FmToGeV;// 100 fm "limit"
  int rhoBins = rhoLimit/(0.001/FmToGeV);// this binning provides small changes between nearby bins
  double rhoStep = rhoLimit/rhoBins;
  double etaLimit = 1/(0.0005*BohrR);
  int etaBins = 2.0*etaLimit/(0.001);// this binning provides small changes between nearby bins
  double etaStep = (2*etaLimit)/etaBins;
  //cout<<rhoBins<<"  "<<etaBins<<"  "<<rhoLimit/rhoBins<<"  "<<etaLimit/etaBins<<endl;
  TH2D *B_histogram = new TH2D("B_histogram","",rhoBins,0,rhoLimit, etaBins,-etaLimit,etaLimit);
  TH2D *P_histogram = new TH2D("P_histogram","",rhoBins,0,rhoLimit, etaBins,-etaLimit,etaLimit);
  TH1D *h_histogram = new TH1D("h_histogram","",etaBins,-etaLimit,etaLimit);
  
  for(int i=1; i<=rhoBins; i++){
    for(int j=1; j<=etaBins; j++){
      double rho = (i+0.5)*rhoStep;
      double eta = (j+0.5)*etaStep - etaLimit;// starts from negative values
      if(fabs(eta) < 0.0001) continue;
      
      double Bfunc[Btermslimit]={0};
      double Pfunc[Btermslimit]={0};
      int StopPoint=Btermslimit;
      Bfunc[0]=1; Bfunc[1]=eta*rho;
      Pfunc[0]=1; Pfunc[1]=0;
      for(int ii=2; ii<Btermslimit; ii++) { 
	Bfunc[ii] = (2*eta*rho*Bfunc[ii-1] - rho*rho*Bfunc[ii-2])/(double(ii*(ii+1.)));
	Pfunc[ii] = (2*eta*rho*Pfunc[ii-1] - rho*rho*Pfunc[ii-2] - (2*(ii-1)+1)*2*eta*rho*Bfunc[ii-1])/(double((ii-1)*ii));
	if(fabs(Bfunc[ii])+fabs(Bfunc[ii-1])+fabs(Pfunc[ii])+fabs(Pfunc[ii-1]) < 1.0e-7) {StopPoint=ii; break;}
      }
      //cout<<StopPoint+1<<endl;
      double Bsum=0, Psum=0;
      for(int ii=0; ii<StopPoint+1; ii++) {Bsum += Bfunc[ii]; Psum += Pfunc[ii];}
      B_histogram->SetBinContent(i,j, Bsum);
      P_histogram->SetBinContent(i,j, Psum);
    }
  }
  cout<<"Done preparing B and P histograms"<<endl;
  for(int j=1; j<=etaBins; j++){
    double eta = (j+0.5)*etaStep - etaLimit;// starts from negative values
    if(fabs(eta) < 0.0001) continue;
    h_histogram->SetBinContent(j, hFunction(eta));
  }
  cout<<"Done preparing h histogram"<<endl;
  
  
  int seconds = time(0);
  
  for(int FileN=0; FileN<Nfiles; FileN++){
    cout<<"File #"<<FileN+1<<endl;
    //if(FileN<1260) continue;
    //if(FileN>230) continue;
    TString *fname=new TString("b");
    *fname += bValue;
    fname->Append("/event");
    *fname += FileN+1;
    fname->Append(".root");
    TFile *infile=new TFile(fname->Data(),"READ");
    if(!infile->IsOpen()) {cout<<fname->Data()<<" does not exist"<<endl; continue;}
    TTree *event_tree=(TTree*)infile->Get("events");
    TTree *particles_tree=(TTree*)infile->Get("particles");
    

    int Nevents = event_tree->GetEntries();
    int NpartList=0;
    int NpartList_past = 0;
    /////////////////////////////////////////
    // Create Pion list first
    for(int Evt=0; Evt<Nevents; Evt++){
      
      //cout<<"Event # "<<Evt+1<<endl;
      event_tree->GetEntry(Evt);
      
      TBranch *br=(TBranch*)event_tree->GetBranch("event");
      TLeaf *lf=(TLeaf*)br->GetLeaf("entries");
      NpartList_past += NpartList;
      NpartList = lf->GetValue();
      
      /////////////////////////////////////
      // past event buffer
      for(int pastEvt=NevtBuff-1; pastEvt>0; pastEvt--){
	
	for(int particle=0; particle<Nparticles[pastEvt-1]; particle++){
	  P[pastEvt][particle][0] = P[pastEvt-1][particle][0];
	  P[pastEvt][particle][1] = P[pastEvt-1][particle][1];
	  P[pastEvt][particle][2] = P[pastEvt-1][particle][2];
	  P[pastEvt][particle][3] = P[pastEvt-1][particle][3];
	  //
	  X[pastEvt][particle][0] = X[pastEvt-1][particle][0];
	  X[pastEvt][particle][1] = X[pastEvt-1][particle][1];
	  X[pastEvt][particle][2] = X[pastEvt-1][particle][2];
	  X[pastEvt][particle][3] = X[pastEvt-1][particle][3];
	  //
	  PID[pastEvt][particle] = PID[pastEvt-1][particle];
	  Primary[pastEvt][particle] = Primary[pastEvt-1][particle];
	}
	
	Nparticles[pastEvt]=Nparticles[pastEvt-1];
      }
      Nparticles[0]=0;
      /////////////////////////////////////
      
     
      // Create Pion list first
      int Npions=0;
      
      for(int index1=NpartList_past; index1<NpartList_past + NpartList; index1++){
	
	if(Npions >= MAXPIONS) {cout<<"Too Many Pions!!!"<<endl; break;}
	particles_tree->GetEntry(index1);
	//
	
	TLeaf *lf_pid=particles_tree->GetLeaf("pid");
	int pid = int(lf_pid->GetValue());
	if(abs(pid) != PIDPOI) continue;// pions only
	TLeaf *lf_decayed=(TLeaf*)particles_tree->GetLeaf("decayed");
	if(lf_decayed->GetValue() != 0) continue;
	
	
	TLeaf *lf_px=(TLeaf*)particles_tree->GetLeaf("px");
	TLeaf *lf_py=(TLeaf*)particles_tree->GetLeaf("py");
	TLeaf *lf_pz=(TLeaf*)particles_tree->GetLeaf("pz");
	TLeaf *lf_x=(TLeaf*)particles_tree->GetLeaf("x");
	TLeaf *lf_y=(TLeaf*)particles_tree->GetLeaf("y");
	TLeaf *lf_z=(TLeaf*)particles_tree->GetLeaf("z");
	TLeaf *lf_t=(TLeaf*)particles_tree->GetLeaf("t");
	
	float px = lf_px->GetValue(); float py = lf_py->GetValue(); float pz = lf_pz->GetValue();
	float x = lf_x->GetValue(); float y = lf_y->GetValue(); float z = lf_z->GetValue();
	double t = lf_t->GetValue();
	
	if(Gaussian){// Gaussian randomization
	  bool goodchoice=kFALSE;
	  while(goodchoice==kFALSE){
	    x = NsigmaGauss*RValue*RandomNumber->Rndm();
	    float height = RandomNumber->Rndm();
	    if(SingleParticleGaussian->Eval(x) >= height) goodchoice=kTRUE;
	  }
	  if(RandomNumber->Rndm()<0.5) x = -x;
	  //
	  goodchoice=kFALSE;
	  while(goodchoice==kFALSE){
	    y = NsigmaGauss*RValue*RandomNumber->Rndm();
	    float height = RandomNumber->Rndm();
	    if(SingleParticleGaussian->Eval(y) >= height) goodchoice=kTRUE;
	  }
	  if(RandomNumber->Rndm()<0.5) y = -y;
	  //
	  goodchoice=kFALSE;
	  while(goodchoice==kFALSE){
	    z = NsigmaGauss*RValue*RandomNumber->Rndm();
	    float height = RandomNumber->Rndm();
	    if(SingleParticleGaussian->Eval(z) >= height) goodchoice=kTRUE;
	  }
	  if(RandomNumber->Rndm()<0.5) z = -z;
	  //
	  goodchoice=kFALSE;
	  while(goodchoice==kFALSE){
	    t = NsigmaGauss*RValue*RandomNumber->Rndm();
	    float height = RandomNumber->Rndm();
	    if(SingleParticleGaussian->Eval(t) >= height) goodchoice=kTRUE;
	  }
	  if(RandomNumber->Rndm()<0.5) t = -t;
	  //t=0;
	}

	float pt = sqrt(pow(px,2)+pow(py,2));
	float eta = asinh(pz/pt);
	float E = sqrt(pow(px,2)+pow(py,2)+pow(pz,2)+pow(massPOI,2));
	
	if(pt<0.16 || pt>1.0) continue;
	if(fabs(eta)>0.8) continue;
	if(sqrt(pow(px,2)+pow(py,2)+pow(pz,2)) > 1.0) continue;
	//
	P[0][Npions][0] = E/SF_Mom; P[0][Npions][1] = px/SF_Mom; P[0][Npions][2] = py/SF_Mom; P[0][Npions][3] = pz/SF_Mom;
	X[0][Npions][0] = t/FmToGeV/SF; X[0][Npions][1] = x/FmToGeV/SF; X[0][Npions][2] = y/FmToGeV/SF; X[0][Npions][3] = z/FmToGeV/SF;
	PID[0][Npions] = pid;
	float r = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
	if(t < 1.e6) Primary[0][Npions]=kTRUE;// 1e6
	else {Primary[0][Npions]=kFALSE;}

	
	if(RemoveXP){
	  if(r > 100./SF) continue;
	}

	if(Primary[0][Npions]) P_dist->Fill(sqrt( pow(px,2) + pow(py,2) + pow(pz,2)));
	
	if(Primary[0][Npions]){
	  r_dist->Fill(r);
	  x_dist->Fill(fabs(X[0][Npions][1])*FmToGeV);
	  Pt_dist->Fill(pt);
	}
	
	//if(Primary[0][Npions]==kFALSE) continue;
	
	//
	
	
	Npions++;
      }
      
      Nparticles[0]=Npions;
      AvgMult->Fill(1, Npions);

      // Shuffle
      for (Int_t i = 0; i < Npions; i++) {
	// generate random index call
	int j = (Int_t) (gRandom->Rndm() * Npions);
	// temporarily store values from random index
	float t = X[0][j][0], x = X[0][j][1], y = X[0][j][2], z = X[0][j][3];
	float E = P[0][j][0], px = P[0][j][1], py = P[0][j][2], pz = P[0][j][3];
	int pid = PID[0][j];
	bool prim = Primary[0][j];
	//
	// Swap value locations
	X[0][j][0] = X[0][i][0]; X[0][j][1] = X[0][i][1]; X[0][j][2] = X[0][i][2]; X[0][j][3] = X[0][i][3];
	P[0][j][0] = P[0][i][0]; P[0][j][1] = P[0][i][1]; P[0][j][2] = P[0][i][2]; P[0][j][3] = P[0][i][3];
	PID[0][j] = PID[0][i];
	Primary[0][j] = Primary[0][i];
	//
	X[0][i][0] = t; X[0][i][1] = x; X[0][i][2] = y; X[0][i][3] = z;
	P[0][i][0] = E; P[0][i][1] = px; P[0][i][2] = py; P[0][i][3] = pz;
	PID[0][i] = pid;
	Primary[0][i] = prim;
      }
      
      ///////////////////////////////////////////
      // de-correlate x-p
      if(RemoveXP){
	for (Int_t i = 0; i < Npions; i++) ShuffledIndex[i] = i;
	Shuffle(ShuffledIndex, 0, Npions-1);
	for (Int_t i = 0; i < Npions; i++) {// decorrelate x-p
	  float t = X[0][i][0], x = X[0][i][1], y = X[0][i][2], z = X[0][i][3];
	  X[0][i][0] = X[0][ShuffledIndex[i]][0]; 
	  X[0][i][1] = X[0][ShuffledIndex[i]][1]; 
	  X[0][i][2] = X[0][ShuffledIndex[i]][2]; 
	  X[0][i][3] = X[0][ShuffledIndex[i]][3];
	  //
	  X[0][ShuffledIndex[i]][0] = t;
	  X[0][ShuffledIndex[i]][1] = x;
	  X[0][ShuffledIndex[i]][2] = y;
	  X[0][ShuffledIndex[i]][3] = z;
	}
      }
      
	
      //////////////////////////////////////////
      

      
      //////////////////////////////////////////////////
      // Start Correlation Analysis

      for(int EN=0; EN<1; EN++){// current or 1st past event
	
	for(int index1=0; index1<Nparticles[EN]; index1++){// current or past event pion loop
	  
	  int start=index1+1;
	  if(EN>0) start=0;
	  for(int index2=start; index2<Npions; index2++){// current event pion loop (EN=0)
	    
	    float kt = sqrt(pow(P[EN][index1][1]+P[0][index2][1],2)+pow(P[EN][index1][2]+P[0][index2][2],2))/2.;
	    float qinv12 = sqrt(pow(P[EN][index1][1]-P[0][index2][1],2)+pow(P[EN][index1][2]-P[0][index2][2],2)+pow(P[EN][index1][3]-P[0][index2][3],2)-pow(P[EN][index1][0]-P[0][index2][0],2));
	    float rstar12_CMS = sqrt(pow(X[EN][index1][1]-X[0][index2][1],2)+pow(X[EN][index1][2]-X[0][index2][2],2)+pow(X[EN][index1][3]-X[0][index2][3],2));
	    
	    if(qinv12 < 0.001) continue;
	    if(qinv12 > qCutOff) continue;

	    if(EN>0){
	      if(PID[EN][index1] == PID[0][index2]){
		DenEM_ss->Fill(kt, qinv12);
	      }else{
		DenEM_os->Fill(kt, qinv12);
	      }
	      if(LambdaDirect) continue;
	    }
	    
	    BoostPRF(P[EN][index1], P_dummy, X[EN][index1], &X1_RF);
	    DenLambda1->Fill(qinv12);
	    if(X1_RF.Mag()<RstarMax) NumLambda1->Fill(qinv12, X1_RF.Mag()*FmToGeV);
	    BoostPRF(P[0][index2], P_dummy, X[0][index2], &X2_RF);
	    DenLambda1->Fill(qinv12);
	    if(X2_RF.Mag()<RstarMax) NumLambda1->Fill(qinv12, X2_RF.Mag()*FmToGeV);
	    
	    
	    P_dist_lowq->Fill(sqrt(pow(P[EN][index1][1],2)+pow(P[EN][index1][2],2)+pow(P[EN][index1][3],2)));
	    P_dist_lowq->Fill(sqrt(pow(P[0][index2][1],2)+pow(P[0][index2][2],2)+pow(P[0][index2][3],2)));
	    
	    // boost P into PRF
	    BoostPRF(P[EN][index1], P[0][index2], P[EN][index1], &P1_PRF);
	    BoostPRF(P[EN][index1], P[0][index2], P[0][index2], &P2_PRF);
	    if(BoostPairsIntoPRF){// boost X into PRF
	      BoostPRF(P[EN][index1], P[0][index2], X[EN][index1], &X1_PRF);
	      BoostPRF(P[EN][index1], P[0][index2], X[0][index2], &X2_PRF);
	    }else{
	      X1_PRF[0]=X[EN][index1][1]; X1_PRF[1]=X[EN][index1][2]; X1_PRF[2]=X[EN][index1][3];
	      X2_PRF[0]=X[0][index2][1]; X2_PRF[1]=X[0][index2][2]; X2_PRF[2]=X[0][index2][3];
	    }
	    //cout<<"out: "<<P1_PRF[0]<<"  "<<P2_PRF[0]<<endl;
	    //cout<<"side: "<<P1_PRF[1]<<"  "<<P2_PRF[1]<<endl;
	    //cout<<"long: "<<P1_PRF[2]<<"  "<<P2_PRF[2]<<" ++++++"<<endl;
	    //cout<<qinv12<<"  "<<sqrt(pow(P1_PRF[0]-P2_PRF[0],2)+pow(P1_PRF[1]-P2_PRF[1],2)+pow(P1_PRF[2]-P2_PRF[2],2))<<endl;
	    //
	    double rstar12_PRF = (X1_PRF-X2_PRF).Mag();
	    
	    if(PID[EN][index1] == PID[0][index2]) DenLambda2_ss->Fill(kt, qinv12);
	    else DenLambda2_os->Fill(kt, qinv12);
	    	    
	    if(rstar12_PRF < RstarMin) continue;// removes most pairs from the same decay
	    
	    if(rstar12_PRF < RstarMax){
	      if(PID[EN][index1] == PID[0][index2]) NumLambda2_ss->Fill(kt, qinv12, rstar12_PRF*FmToGeV);
	      else NumLambda2_os->Fill(kt, qinv12, rstar12_PRF*FmToGeV);
	    }

	    if(rstar12_PRF > RstarMax) {
	      if(PID[EN][index1]==PID[0][index2]) LargeRpairs_ss->Fill(kt, qinv12);
	      else LargeRpairs_os->Fill(kt, qinv12);
	      if(!LambdaDirect) continue;
	    }
	    float FVP12 = (P1_PRF-P2_PRF).Mag();
	    FVP12 -= pow(P[EN][index1][0]-P[0][index2][0],2);
	    double tstar12_PRF = sqrt(fabs(pow(rstar12_PRF,2)-FVP12));// fabs biases only 0.00001 % of the pairs (likely from same decay)
	    
	    	    
	 	    
	    //
	    double rho = rstar12_PRF * qinv12/2.;
	    double phase = (P1_PRF-P2_PRF)*(X1_PRF-X2_PRF);
	    double Cosphase = phase/(rstar12_PRF)/(qinv12);
	    
	    double eta = 1/(qinv12/2.*BohrR);
	    if(PID[EN][index1] != PID[0][index2]) eta = -eta;
	    complex<double> ainput(0.0, -eta);
	    complex<double> zinput12(0.0, rho + rho*Cosphase);
	    complex<double> zinput21(0.0, rho - rho*Cosphase);
	    complex<double> CHG12(1.0, 0);
	    complex<double> CHG21(1.0, 0);
	    double GamovFactor=1.0;
	    double BEE = 0.;
	    double WeightWF = 1.0;
	    double WeightPlaneWF = 1.0;
	    // Strong FSI
	    complex<double> Psi_alpha(0., 0.);
	    if(EN==0 && !LambdaDirect && qinv12 < qCutOff && rstar12_PRF < RstarMax){
	      if(Primary[EN][index1]==kTRUE && Primary[0][index2]==kTRUE){
		if(StrongFSI && PID[EN][index1]!=PID[0][index2]){
		  double kb = sqrt(pow(massPOI,2)+pow(qinv12/2.,2) - pow(masspi0,2));
		  double G2 = KinverseFunction(pow(qinv12/2.,2),1);
		  double G1 = KinverseFunction(pow(qinv12/2.,2),0);
		  double RK11 = 2/3.*G1 + 1/3.*G2;
		  double RK22 = 2/3.*G2 + 1/3.*G1;
		  double RK12 = -sqrt(2/9.)*(G1-G2);
		  complex<double> Chi(h_histogram->GetBinContent(h_histogram->GetXaxis()->FindBin(eta)), Gamov(eta)/(2*eta));
		  complex<double> C3 = RK11 - 2.0*Chi/BohrR;
		  complex<double> C5(RK22, -kb);
		  complex<double> C10 = C3*C5 - pow(RK12,2);
		  complex<double> fc_aa = C5/C10;
		  complex<double> G = P_histogram->GetBinContent(P_histogram->GetXaxis()->FindBin(rho), P_histogram->GetYaxis()->FindBin(eta));
		  G += 2.0*eta*rho*( log(fabs(2.0*eta*rho)) + 2.0*EulerC - 1.0 + Chi )*B_histogram->GetBinContent(B_histogram->GetXaxis()->FindBin(rho), B_histogram->GetYaxis()->FindBin(eta));
		  Psi_alpha = fc_aa*G/(rstar12_PRF);
		}
		
		CHG12= conf_hyp(ainput,binput,zinput12);
		CHG21= conf_hyp(ainput,binput,zinput21);
		GamovFactor = Gamov(eta);
		BEE = cos(phase);
	
		complex<double> planewave12(double(TMath::Cos(-rho*Cosphase)), double(TMath::Sin(-rho*Cosphase)));
		complex<double> planewave21(planewave12.real(), -planewave12.imag());
		//
		complex<double> Full12(0.0,0.0);
		Full12 += planewave12;
		Full12 *= CHG12;
		if(PID[EN][index1]!=PID[0][index2]) Full12 += Psi_alpha;
		//
		complex<double> Full21(0.0,0.0);
		Full21 += planewave21;
		Full21 *= CHG21;
		//
		complex<double> FullWF = Full12;
		complex<double> FullPlaneWF = planewave12;
		if(PID[EN][index1] == PID[0][index2]) {
		  FullWF += Full21;
		  FullWF /= sqrt(2.);
		  FullPlaneWF += planewave21;
		  FullPlaneWF /= sqrt(2.);
		}
		
		WeightWF = GamovFactor*pow(abs(FullWF),2);
		WeightPlaneWF = pow(abs(FullPlaneWF),2);
		
		// Fill undiluted histos
		if(PID[EN][index1] == PID[0][index2]){
		  Num_PrimCosFSI_ss->Fill(kt, qinv12, WeightWF);
		}else {
		  Num_PrimCosFSI_os->Fill(kt, qinv12, WeightWF);
		}
		
	      }// primary
	    }// same-event, LambdaDirect, rstar and qinv cut
	    
	    
	    if(PID[EN][index1] == PID[0][index2]){// same-charge
	      rstarPRF_dist_ss->Fill(kt, rstar12_PRF*FmToGeV);
	      tstar_dist_ss->Fill(kt, fabs(X[EN][index1][0]-X[0][index2][0])*FmToGeV);
	      Num_Cos_ss->Fill(kt, qinv12, BEE);
	      Num_CosFSI_ss->Fill(kt, qinv12, WeightWF);
	      NumSq_CosFSI_ss->Fill(kt, qinv12, pow(WeightWF,2));
	      Den_ss->Fill(kt, qinv12);
	      //
	      K2_ss->Fill(kt, qinv12, WeightWF);
	      PlaneWF_ss->Fill(kt, qinv12, WeightPlaneWF);
	    }else {// mixed-charge
	      rstarPRF_dist_os->Fill(kt, rstar12_PRF*FmToGeV);
	      tstar_dist_os->Fill(kt, fabs(X[EN][index1][0]-X[0][index2][0])*FmToGeV);
	      Num_Cos_os->Fill(kt, qinv12, BEE);
	      Num_CosFSI_os->Fill(kt, qinv12, WeightWF);
	      NumSq_CosFSI_os->Fill(kt, qinv12, pow(WeightWF,2));
	      Den_os->Fill(kt, qinv12);
	      //
	      K2_os->Fill(kt, qinv12, WeightWF);
	      PlaneWF_os->Fill(kt, qinv12, WeightPlaneWF);
	    }


	    if(!ThreeParticle) continue;


	    TLorentzVector LP1(P[EN][index1][1], P[EN][index1][2], P[EN][index1][3], P[EN][index1][0]);
	    TLorentzVector LX1(X[EN][index1][1], X[EN][index1][2], X[EN][index1][3], X[EN][index1][0]);
	    TLorentzVector LP2(P[0][index2][1], P[0][index2][2], P[0][index2][3], P[0][index2][0]);
	    TLorentzVector LX2(X[0][index2][1], X[0][index2][2], X[0][index2][3], X[0][index2][0]);
	    
	    
	    if(!LambdaDirect) {if(!Primary[EN][index1] || !Primary[0][index2]) continue;}
	    if(qinv12 > qCutOff) continue;
	    //
	    int EN2=0;
	    int start3=index2+1;
	    if(EN==1) {EN2=2; start3=0;}
	    
	    for(int index3=start3; index3<Nparticles[EN2]; index3++){
	     
	      if(!LambdaDirect && !Primary[EN2][index3]) continue;
	      
	      double qinv23 = sqrt(pow(P[0][index2][1]-P[EN2][index3][1],2)+pow(P[0][index2][2]-P[EN2][index3][2],2)+pow(P[0][index2][3]-P[EN2][index3][3],2)-pow(P[0][index2][0]-P[EN2][index3][0],2));
	      double qinv31 = sqrt(pow(P[EN2][index3][1]-P[EN][index1][1],2)+pow(P[EN2][index3][2]-P[EN][index1][2],2)+pow(P[EN2][index3][3]-P[EN][index1][3],2)-pow(P[EN2][index3][0]-P[EN][index1][0],2));
	      if(qinv23 > qCutOff || qinv31 > qCutOff) continue;
	      if(qinv23 < 0.001 || qinv31 < 0.001) continue;
	      
	      if(PID[EN][index1] == PID[0][index2] && PID[EN][index1] == PID[EN2][index3]){
		DenLambda3_ss->Fill(qinv12, qinv23, qinv31);
	      }else{
		DenLambda3_os->Fill(qinv12, qinv23, qinv31);
	      }
	      
	      	      
	      bool MixedCharge=kFALSE;
	      if(PID[EN][index1] != PID[0][index2] || PID[EN][index1] != PID[EN2][index3] || PID[0][index2] != PID[EN2][index3]) MixedCharge=kTRUE;
	      
	      if(EN2==2){
		if(MixedCharge==kFALSE){
		  r3numME->Fill(qinv12, qinv23, qinv31);
		  r3den1ME->Fill(qinv12, qinv23, qinv31);
		  r3den2ME->Fill(qinv12, qinv23, qinv31);
		  r3den3ME->Fill(qinv12, qinv23, qinv31); 
		}
		continue;
	      }
	      
	      

	      // 12 boost
	      BoostPRF(P[EN][index1], P[0][index2], P[EN2][index3], &P3_PRF);
	      BoostPRF(P[EN][index1], P[0][index2], P[0][index2], &P2_PRF);
	      BoostPRF(P[EN][index1], P[0][index2], P[EN][index1], &P1_PRF);
	      X1_PRF[0]=X[EN][index1][1]; X1_PRF[1]=X[EN][index1][2]; X1_PRF[2]=X[EN][index1][3];// default CMS
	      X2_PRF[0]=X[0][index2][1]; X2_PRF[1]=X[0][index2][2]; X2_PRF[2]=X[0][index2][3];// default CMS
	      X3_PRF[0]=X[EN2][index3][1]; X3_PRF[1]=X[EN2][index3][2]; X3_PRF[2]=X[EN2][index3][3];// default CMS
	      
	      //
	      if(BoostPairsIntoPRF){
		BoostPRF(P[EN][index1], P[0][index2], X[EN2][index3], &X3_PRF);
		BoostPRF(P[EN][index1], P[0][index2], X[0][index2], &X2_PRF);
		BoostPRF(P[EN][index1], P[0][index2], X[EN][index1], &X1_PRF);
	      }
	

	      rstar12_PRF = (X1_PRF-X2_PRF).Mag();
	      if(rstar12_PRF > RstarMax) continue;
	      if(rstar12_PRF < RstarMin) continue;
	      double phase_p12_x12 = (P1_PRF-P2_PRF)*(X1_PRF-X2_PRF)/2.;
	      double phase_p12_x32 = (P1_PRF-P2_PRF)*(X3_PRF-X2_PRF)/2.;
	      double phase_p12_x13 = (P1_PRF-P2_PRF)*(X1_PRF-X3_PRF)/2.;
	      complex<double> z_p12_x12(0.0, (P1_PRF-P2_PRF)*(X1_PRF-X2_PRF)/2. + qinv12/2. * (X1_PRF-X2_PRF).Mag());
	      complex<double> z_p12_x21(0.0, (P1_PRF-P2_PRF)*(X2_PRF-X1_PRF)/2. + qinv12/2. * (X2_PRF-X1_PRF).Mag());
	      complex<double> z_p12_x32(0.0, (P1_PRF-P2_PRF)*(X3_PRF-X2_PRF)/2. + qinv12/2. * (X3_PRF-X2_PRF).Mag());
	      complex<double> z_p12_x13(0.0, (P1_PRF-P2_PRF)*(X1_PRF-X3_PRF)/2. + qinv12/2. * (X1_PRF-X3_PRF).Mag());
	      complex<double> z_p12_x23(0.0, (P1_PRF-P2_PRF)*(X2_PRF-X3_PRF)/2. + qinv12/2. * (X2_PRF-X3_PRF).Mag());
	      complex<double> z_p12_x31(0.0, (P1_PRF-P2_PRF)*(X3_PRF-X1_PRF)/2. + qinv12/2. * (X3_PRF-X1_PRF).Mag());
	      
	      
	      // 23 boost
	      BoostPRF(P[0][index2], P[EN2][index3], P[EN2][index3], &P3_PRF);
	      BoostPRF(P[0][index2], P[EN2][index3], P[0][index2], &P2_PRF);
	      BoostPRF(P[0][index2], P[EN2][index3], P[EN][index1], &P1_PRF);
	      if(BoostPairsIntoPRF){
		BoostPRF(P[0][index2], P[EN2][index3], X[EN2][index3], &X3_PRF);
		BoostPRF(P[0][index2], P[EN2][index3], X[0][index2], &X2_PRF);
		BoostPRF(P[0][index2], P[EN2][index3], X[EN][index1], &X1_PRF);
	      }
	      double rstar23_PRF = (X2_PRF-X3_PRF).Mag();
	      if(rstar23_PRF > RstarMax) continue;
	      if(rstar23_PRF < RstarMin) continue;
	      double phase_p23_x13 = (P2_PRF-P3_PRF)*(X1_PRF-X3_PRF)/2.;
	      double phase_p23_x23 = (P2_PRF-P3_PRF)*(X2_PRF-X3_PRF)/2.;
	      double phase_p23_x21 = (P2_PRF-P3_PRF)*(X2_PRF-X1_PRF)/2.;
	      double FVP23 = (P2_PRF-P3_PRF).Mag();
	      FVP23 -= pow(P[0][index2][0]-P[EN2][index3][0],2);
	      double tstar23_PRF = sqrt(fabs(pow(rstar23_PRF,2)-FVP23));
	      //double rstar23_CMS = sqrt(pow(X[0][index2][1]-X[EN2][index3][1],2)+pow(X[0][index2][2]-X[EN2][index3][2],2)+pow(X[0][index2][3]-X[EN2][index3][3],2));
	      complex<double> z_p23_x23(0.0, (P2_PRF-P3_PRF)*(X2_PRF-X3_PRF)/2. + qinv23/2. * (X2_PRF-X3_PRF).Mag());
	      complex<double> z_p23_x13(0.0, (P2_PRF-P3_PRF)*(X1_PRF-X3_PRF)/2. + qinv23/2. * (X1_PRF-X3_PRF).Mag());
	      complex<double> z_p23_x21(0.0, (P2_PRF-P3_PRF)*(X2_PRF-X1_PRF)/2. + qinv23/2. * (X2_PRF-X1_PRF).Mag());
	      complex<double> z_p23_x32(0.0, (P2_PRF-P3_PRF)*(X3_PRF-X2_PRF)/2. + qinv23/2. * (X3_PRF-X2_PRF).Mag());
	      complex<double> z_p23_x31(0.0, (P2_PRF-P3_PRF)*(X3_PRF-X1_PRF)/2. + qinv23/2. * (X3_PRF-X1_PRF).Mag());
	      complex<double> z_p23_x12(0.0, (P2_PRF-P3_PRF)*(X1_PRF-X2_PRF)/2. + qinv23/2. * (X1_PRF-X2_PRF).Mag());
	      

	      // 31 boost
	      BoostPRF(P[EN2][index3], P[EN][index1], P[EN2][index3], &P3_PRF);
	      BoostPRF(P[EN2][index3], P[EN][index1], P[0][index2], &P2_PRF);
	      BoostPRF(P[EN2][index3], P[EN][index1], P[EN][index1], &P1_PRF);
	      if(BoostPairsIntoPRF){
		BoostPRF(P[EN2][index3], P[EN][index1], X[EN2][index3], &X3_PRF);
		BoostPRF(P[EN2][index3], P[EN][index1], X[0][index2], &X2_PRF);
		BoostPRF(P[EN2][index3], P[EN][index1], X[EN][index1], &X1_PRF);
	      }
	      double rstar31_PRF = (X3_PRF-X1_PRF).Mag();
	      if(rstar31_PRF > RstarMax) continue;
	      if(rstar31_PRF < RstarMin) continue;
	      double phase_p31_x31 = (P3_PRF-P1_PRF)*(X3_PRF-X1_PRF)/2.;
	      double phase_p31_x32 = (P3_PRF-P1_PRF)*(X3_PRF-X2_PRF)/2.;
	      double phase_p31_x21 = (P3_PRF-P1_PRF)*(X2_PRF-X1_PRF)/2.;
	      double FVP31 = (P3_PRF-P1_PRF).Mag();
	      FVP31 -= pow(P[EN2][index3][0]-P[EN][index1][0],2);
	      double tstar31_PRF = sqrt(fabs(pow(rstar31_PRF,2)-FVP31));
	      //double rstar31_CMS = sqrt(pow(X[EN2][index3][1]-X[EN][index1][1],2)+pow(X[EN2][index3][2]-X[EN][index1][2],2)+pow(X[EN2][index3][3]-X[EN][index1][3],2));
	      complex<double> z_p31_x31(0.0, (P3_PRF-P1_PRF)*(X3_PRF-X1_PRF)/2. + qinv31/2. * (X3_PRF-X1_PRF).Mag());
	      complex<double> z_p31_x32(0.0, (P3_PRF-P1_PRF)*(X3_PRF-X2_PRF)/2. + qinv31/2. * (X3_PRF-X2_PRF).Mag());
	      complex<double> z_p31_x13(0.0, (P3_PRF-P1_PRF)*(X1_PRF-X3_PRF)/2. + qinv31/2. * (X1_PRF-X3_PRF).Mag());
	      complex<double> z_p31_x21(0.0, (P3_PRF-P1_PRF)*(X2_PRF-X1_PRF)/2. + qinv31/2. * (X2_PRF-X1_PRF).Mag());
	      complex<double> z_p31_x12(0.0, (P3_PRF-P1_PRF)*(X1_PRF-X2_PRF)/2. + qinv31/2. * (X1_PRF-X2_PRF).Mag());
	      complex<double> z_p31_x23(0.0, (P3_PRF-P1_PRF)*(X2_PRF-X3_PRF)/2. + qinv31/2. * (X2_PRF-X3_PRF).Mag());
	      
	      

	      ////////////////////////////////////////////
	      ////////////////////////////////////////////
	      double Q3 = sqrt(pow(qinv12,2) + pow(qinv23,2) + pow(qinv31,2));
	      
	      if(!LambdaDirect){
		int CP12=+1, CP23=+1, CP31=+1;
		if(PID[EN][index1] != PID[0][index2]) CP12 = -1;
		if(PID[0][index2] != PID[EN2][index3]) CP23 = -1;
		if(PID[EN2][index3] != PID[EN][index1]) CP31 = -1;
		
		double Ac_12 = Gamov(CP12/(qinv12/2.*BohrR));
		double Ac_23 = Gamov(CP23/(qinv23/2.*BohrR));
		double Ac_31 = Gamov(CP31/(qinv31/2.*BohrR));
	
		//
		complex<double> a_12(0.0, -CP12/(qinv12/2.*BohrR));
		complex<double> a_23(0.0, -CP23/(qinv23/2.*BohrR));
		complex<double> a_31(0.0, -CP31/(qinv31/2.*BohrR));
		//
		//
		//
	      TLorentzVector LP3(P[EN2][index3][1], P[EN2][index3][2], P[EN2][index3][3], P[EN2][index3][0]);
	      TLorentzVector LX3(X[EN2][index3][1], X[EN2][index3][2], X[EN2][index3][3], X[EN2][index3][0]);
	      
	      double phasePW = ( (LP1-LP2)*(LX1-LX2) + (LP1-LP3)*(LX1-LX3) + (LP2-LP3)*(LX2-LX3) )/3.;// base
	      complex<double> planewave1(cos(phasePW), sin(phasePW));
	      //
	      phasePW = ( (LP1-LP2)*(LX2-LX1) + (LP1-LP3)*(LX2-LX3) + (LP2-LP3)*(LX1-LX3) )/3.;// 12 swap
	      complex<double> planewave2(cos(phasePW), sin(phasePW));
	      //
	      phasePW = ( (LP1-LP2)*(LX3-LX2) + (LP1-LP3)*(LX3-LX1) + (LP2-LP3)*(LX2-LX1) )/3.;// 13 swap
	      complex<double> planewave3(cos(phasePW), sin(phasePW));
	      //
	      phasePW = ( (LP1-LP2)*(LX1-LX3) + (LP1-LP3)*(LX1-LX2) + (LP2-LP3)*(LX3-LX2) )/3.;// 23 swap
	      complex<double> planewave4(cos(phasePW), sin(phasePW));
	      //
	      phasePW = ( (LP1-LP2)*(LX2-LX3) + (LP1-LP3)*(LX2-LX1) + (LP2-LP3)*(LX3-LX1) )/3.;// 123 swap (1st type)
	      complex<double> planewave5(cos(phasePW), sin(phasePW));
	      //
	      phasePW = ( (LP1-LP2)*(LX3-LX1) + (LP1-LP3)*(LX3-LX2) + (LP2-LP3)*(LX1-LX2) )/3.;// 123 swap (2nd type)
	      complex<double> planewave6(cos(phasePW), sin(phasePW));
	      //
	      //
	      //
	      complex<double> CHG_p12_x12= conf_hyp(a_12, binput, z_p12_x12);
	      complex<double> CHG_p12_x21= conf_hyp(a_12, binput, z_p12_x21);
	      complex<double> CHG_p12_x32= conf_hyp(a_12, binput, z_p12_x32);
	      complex<double> CHG_p12_x13= conf_hyp(a_12, binput, z_p12_x13);
	      complex<double> CHG_p12_x23= conf_hyp(a_12, binput, z_p12_x23);
	      complex<double> CHG_p12_x31= conf_hyp(a_12, binput, z_p12_x31);
	      //
	      complex<double> CHG_p31_x31= conf_hyp(a_31, binput, z_p31_x31);
	      complex<double> CHG_p31_x32= conf_hyp(a_31, binput, z_p31_x32);
	      complex<double> CHG_p31_x13= conf_hyp(a_31, binput, z_p31_x13);
	      complex<double> CHG_p31_x21= conf_hyp(a_31, binput, z_p31_x21);
	      complex<double> CHG_p31_x12= conf_hyp(a_31, binput, z_p31_x12);
	      complex<double> CHG_p31_x23= conf_hyp(a_31, binput, z_p31_x23);
	      //
	      complex<double> CHG_p23_x23= conf_hyp(a_23, binput, z_p23_x23);
	      complex<double> CHG_p23_x13= conf_hyp(a_23, binput, z_p23_x13);
	      complex<double> CHG_p23_x21= conf_hyp(a_23, binput, z_p23_x21);
	      complex<double> CHG_p23_x32= conf_hyp(a_23, binput, z_p23_x32);
	      complex<double> CHG_p23_x31= conf_hyp(a_23, binput, z_p23_x31);
	      complex<double> CHG_p23_x12= conf_hyp(a_23, binput, z_p23_x12);

	      ThreePhase_1 = X1_PRF*(P1_PRF-P2_PRF) + X2_PRF*(P2_PRF-P3_PRF) + X3_PRF*(P3_PRF-P1_PRF);
	      ThreePhase_2 = X1_PRF*(P1_PRF-P2_PRF) + X3_PRF*(P2_PRF-P3_PRF) + X2_PRF*(P3_PRF-P1_PRF);
	      ThreePhase_3 = X2_PRF*(P1_PRF-P2_PRF) + X1_PRF*(P2_PRF-P3_PRF) + X3_PRF*(P3_PRF-P1_PRF);
	      ThreePhase_4 = X2_PRF*(P1_PRF-P2_PRF) + X3_PRF*(P2_PRF-P3_PRF) + X1_PRF*(P3_PRF-P1_PRF);
	      ThreePhase_5 = X3_PRF*(P1_PRF-P2_PRF) + X1_PRF*(P2_PRF-P3_PRF) + X2_PRF*(P3_PRF-P1_PRF);
	      ThreePhase_6 = X3_PRF*(P1_PRF-P2_PRF) + X2_PRF*(P2_PRF-P3_PRF) + X1_PRF*(P3_PRF-P1_PRF);
	      TwoPhase_12 = (P1_PRF-P2_PRF)*(X1_PRF-X2_PRF);
	      TwoPhase_23 = (P2_PRF-P3_PRF)*(X2_PRF-X3_PRF);
	      TwoPhase_31 = (P3_PRF-P1_PRF)*(X3_PRF-X1_PRF);
	      //
	      //
	      
	      // Strong FSI
	      if(MixedCharge){
		if(CP12==+1){// 1 and 2 are same-charge
		  kstar_1 = qinv31/2.; kstar_2 = qinv23/2.;
		  rstar_1 = rstar31_PRF; rstar_2 = rstar23_PRF;
		  phase_1 = phase_p31_x31; phase_2 = phase_p31_x32; phase_3 = phase_p23_x13; phase_4 = phase_p23_x23;
		  CHG_1 = CHG_p31_x31;
		  CHG_2 = CHG_p31_x32;
		  CHG_3 = CHG_p23_x13;
		  CHG_4 = CHG_p23_x23;
		  qinvMixedC_ss=qinv12; qinvMixedC_os1=qinv31; qinvMixedC_os2=qinv23; 
		}else if(CP31==+1){// 1 and 3 are same-charge
		  kstar_1 = qinv12/2.; kstar_2 = qinv23/2.;
		  rstar_1 = rstar12_PRF; rstar_2 = rstar23_PRF;
		  phase_1 = phase_p12_x12; phase_2 = phase_p12_x32; phase_3 = phase_p23_x21; phase_4 = phase_p23_x23;
		  CHG_1 = CHG_p12_x12;
		  CHG_2 = CHG_p12_x32;
		  CHG_3 = CHG_p23_x21;
		  CHG_4 = CHG_p23_x23;
		  qinvMixedC_ss=qinv31; qinvMixedC_os1=qinv12; qinvMixedC_os2=qinv23; 
		}else {// 2 and 3 are same-charge
		  kstar_1 = qinv12/2.; kstar_2 = qinv31/2.;
		  rstar_1 = rstar12_PRF; rstar_2 = rstar31_PRF;
		  phase_1 = phase_p12_x12; phase_2 = phase_p12_x13; phase_3 = phase_p31_x21; phase_4 = phase_p31_x31;
		  CHG_1 = CHG_p12_x12;
		  CHG_2 = CHG_p12_x13;
		  CHG_3 = CHG_p31_x21;
		  CHG_4 = CHG_p31_x31;
		  qinvMixedC_ss=qinv23; qinvMixedC_os1=qinv31; qinvMixedC_os2=qinv31; 
		}


		for(int piece=0; piece<4; piece++){
		  if(piece==0) {// P31X31
		    kstar=kstar_1;
		    rstar=rstar_1; 
		    phase=phase_1;
		    CHG=CHG_1;
		  }else if(piece==1) {// P31X32
		    kstar=kstar_1;
		    rstar=rstar_2; 
		    phase=phase_2;
		    CHG=CHG_2;
		  }else if(piece==2) {// P23X13
		    kstar=kstar_2;
		    rstar=rstar_1; 
		    phase=phase_3;
		    CHG=CHG_3;
		  }else {// P23X23
		    kstar=kstar_2;
		    rstar=rstar_2; 
		    phase=phase_4;
		    CHG=CHG_4;
		  }
		  
		  rho = kstar*rstar;
		  eta = -1./(kstar*BohrR);
		  double kb = sqrt(pow(massPOI,2)+pow(kstar,2) - pow(masspi0,2));
		  double G2 = KinverseFunction(pow(kstar,2),1);
		  double G1 = KinverseFunction(pow(kstar,2),0);
		  double RK11 = 2/3.*G1 + 1/3.*G2;
		  double RK22 = 2/3.*G2 + 1/3.*G1;
		  double RK12 = -sqrt(2/9.)*(G1-G2);
		  complex<double> Chi(h_histogram->GetBinContent(h_histogram->GetXaxis()->FindBin(eta)), Gamov(eta)/(2*eta));
		  complex<double> C3 = RK11 - 2.0*Chi/BohrR;
		  complex<double> C5(RK22, -kb);
		  complex<double> C10 = C3*C5 - pow(RK12,2);
		  complex<double> fc_aa = C5/C10;
		  complex<double> fc_ab = -RK12/C10;
		  complex<double> G = P_histogram->GetBinContent(P_histogram->GetXaxis()->FindBin(rho), P_histogram->GetYaxis()->FindBin(eta));
		  G += 2.0*eta*rho*( log(fabs(2.0*eta*rho)) + 2.0*EulerC - 1.0 + Chi )*B_histogram->GetBinContent(B_histogram->GetXaxis()->FindBin(rho), B_histogram->GetYaxis()->FindBin(eta));
		  complex<double> PW(cos(phase), sin(phase));
		  AlphaWF[piece] = CHG + PW*fc_aa*G/rstar;
		  complex<double> spherewave(cos(rstar*kb), sin(rstar*kb));
		  BetaWF[piece] = PW*fc_ab*sqrt(masspi0/massPOI)*spherewave/rstar; 
		  
		}// piece
	      
	      }// Mixed-charge case with Strong FSI
	      
	      complex<double> FullFSI_WF[2];// aa, bb
	      if(MixedCharge==kTRUE){
		if(StrongFSI==kFALSE) {
		  FullFSI_WF[0] = planewave1*CHG_p12_x12*CHG_p31_x31*CHG_p23_x23;
		  if(CP12==+1){// 1 and 2 are same-charge
		    FullFSI_WF[0] += planewave2*CHG_p12_x21*CHG_p31_x32*CHG_p23_x13;// 12 symmetrization
		  }else if(CP31==+1){// 1 and 3 are same-charge
		    FullFSI_WF[0] += planewave3*CHG_p12_x32*CHG_p31_x13*CHG_p23_x21;// 13 symmetrization
		  }else{
		    FullFSI_WF[0] += planewave4*CHG_p12_x13*CHG_p31_x21*CHG_p23_x32;// 23 symmetrization
		  }
		  FullFSI_WF[1] = 0.;
		}else {// strong FSI 
		  if(CP12==+1){// 1 and 2 are same-charge
		    FullFSI_WF[0] = planewave1*CHG_p12_x12*AlphaWF[0]*AlphaWF[3];
		    FullFSI_WF[1] = planewave1*CHG_p12_x12*BetaWF[0]*BetaWF[3];
		    FullFSI_WF[0] += planewave2*CHG_p12_x21*AlphaWF[1]*AlphaWF[2];// 12 symmetrization
		    FullFSI_WF[1] += planewave2*CHG_p12_x21*BetaWF[1]*BetaWF[2];// 12 symmetrization
		  }else if(CP31==+1){// 1 and 3 are same-charge
		    FullFSI_WF[0] = planewave1*CHG_p31_x31*AlphaWF[0]*AlphaWF[3];
		    FullFSI_WF[1] = planewave1*CHG_p31_x31*BetaWF[0]*BetaWF[3];
		    FullFSI_WF[0] += planewave3*CHG_p31_x13*AlphaWF[1]*AlphaWF[2];// 13 symmetrization
		    FullFSI_WF[1] += planewave3*CHG_p31_x13*BetaWF[1]*BetaWF[2];// 13 symmetrization
		  }else{
		    FullFSI_WF[0] = planewave1*CHG_p23_x23*AlphaWF[0]*AlphaWF[3];
		    FullFSI_WF[1] = planewave1*CHG_p23_x23*BetaWF[0]*BetaWF[3];
		    FullFSI_WF[0] += planewave4*CHG_p23_x32*AlphaWF[1]*AlphaWF[2];// 13 symmetrization
		    FullFSI_WF[1] += planewave4*CHG_p23_x32*BetaWF[1]*BetaWF[2];// 13 symmetrization
		  }
		}
		FullFSI_WF[0] /= sqrt(2.); FullFSI_WF[1] /= sqrt(2.);
	      }else {
		 FullFSI_WF[0] = planewave1*CHG_p12_x12*CHG_p31_x31*CHG_p23_x23;// base
		 FullFSI_WF[0] += planewave2*CHG_p12_x21*CHG_p31_x32*CHG_p23_x13;// 12 symmetrization
		 FullFSI_WF[0] += planewave3*CHG_p12_x32*CHG_p31_x13*CHG_p23_x21;// 13
		 FullFSI_WF[0] += planewave4*CHG_p12_x13*CHG_p31_x21*CHG_p23_x32;// 23
		 FullFSI_WF[0] += planewave5*CHG_p12_x23*CHG_p31_x12*CHG_p23_x31;// 123
		 FullFSI_WF[0] += planewave6*CHG_p12_x31*CHG_p31_x23*CHG_p23_x12;// 123
		 FullFSI_WF[1] = 0.;
		 FullFSI_WF[0] /= sqrt(6.); FullFSI_WF[1] /= sqrt(6.);
	      }

	      complex<double> FullPlaneWF = planewave1;// base
	      if(MixedCharge && CP12==+1) FullPlaneWF += planewave2;
	      if(MixedCharge && CP31==+1) FullPlaneWF += planewave3;
	      if(MixedCharge && CP23==+1) FullPlaneWF += planewave4;
	      if(!MixedCharge) {
		FullPlaneWF += planewave2 + planewave3 + planewave4 + planewave5 + planewave6;
		FullPlaneWF /= sqrt(6.);
	      }else FullPlaneWF /= sqrt(2.);
	      
	      WeightWF_aa = pow(abs(FullFSI_WF[0]),2);
	      WeightWF_bb = pow(abs(FullFSI_WF[1]),2);
	      WeightWF_aa *= Ac_12*Ac_23*Ac_31;
	      WeightWF_bb *= Ac_12*Ac_23*Ac_31;
	      WeightPlaneWF = pow(abs(FullPlaneWF),2);
	      
	      }// !LambdaDirect
	      ////////////////////////////////////////////
	      ////////////////////////////////////////////
	      //
	     	     
	      
	      if(MixedCharge){
		K3os->Fill(Q3, WeightWF_aa + WeightWF_bb);
		PlaneWF3os->Fill(Q3, WeightPlaneWF);
		K3os_3D->Fill(qinvMixedC_ss, qinvMixedC_os1, qinvMixedC_os2, WeightWF_aa + WeightWF_bb);
		PlaneWF3os_3D->Fill(qinvMixedC_ss, qinvMixedC_os1, qinvMixedC_os2, WeightPlaneWF);
		NumLambda3_os->Fill(qinv12, qinv23, qinv31);
	      }else{
		K3ss->Fill(Q3, WeightWF_aa + WeightWF_bb);
		PlaneWF3ss->Fill(Q3, WeightPlaneWF);
		K3ss_3D->Fill(qinv12, qinv23, qinv31, WeightWF_aa + WeightWF_bb);
		PlaneWF3ss_3D->Fill(qinv12, qinv23, qinv31, WeightPlaneWF);
		NumLambda3_ss->Fill(qinv12, qinv23, qinv31);
		rstar3_dist_ss->Fill(rstar12_PRF, rstar23_PRF, rstar31_PRF);
		//rstar3_dist_ss->Fill(rstar12_CMS, rstar23_CMS, rstar31_CMS);
		tstar3_dist_ss->Fill(tstar12_PRF, tstar23_PRF, tstar31_PRF);
		//
		double cosThreePhase = 1/6. * (cos(ThreePhase_1)+cos(ThreePhase_2)+cos(ThreePhase_3)+cos(ThreePhase_4)+cos(ThreePhase_5)+cos(ThreePhase_6));
		r3num->Fill(qinv12, qinv23, qinv31, 2*cosThreePhase);
		r3den1->Fill(qinv12, qinv23, qinv31, cos(TwoPhase_12));
		r3den2->Fill(qinv12, qinv23, qinv31, cos(TwoPhase_23));
		r3den3->Fill(qinv12, qinv23, qinv31, cos(TwoPhase_31));
		r3numSq->Fill(qinv12, qinv23, qinv31, pow(2*cosThreePhase,2));
		r3den1Sq->Fill(qinv12, qinv23, qinv31, pow(cos(TwoPhase_12),2));
		r3den2Sq->Fill(qinv12, qinv23, qinv31, pow(cos(TwoPhase_23),2));
		r3den3Sq->Fill(qinv12, qinv23, qinv31, pow(cos(TwoPhase_31),2));
		r3numEn->Fill(qinv12, qinv23, qinv31);
		r3den1En->Fill(qinv12, qinv23, qinv31);
		r3den2En->Fill(qinv12, qinv23, qinv31);
		r3den3En->Fill(qinv12, qinv23, qinv31); 
		//
		Num_Cos12->Fill(qinv12,cos(TwoPhase_12));
		Num_Cos23->Fill(qinv23,cos(TwoPhase_23));
		Num_Cos31->Fill(qinv31,cos(TwoPhase_31));
		Den_Cos12->Fill(qinv12);
		Den_Cos23->Fill(qinv23);
		Den_Cos31->Fill(qinv31);
	      }
	      
	      
	    }// index3
	    
	  }// index2
	  
	}// index1
	
      }// EN buffer      
      
    }// Nevents
    
    infile->Close();
  }// Nfiles
  
  outfile->cd();
  //
  r_dist->Write();
  x_dist->Write();
  Pt_dist->Write();
  P_dist->Write();
  P_dist_lowq->Write();
  AvgMult->Write();

  rstarPRF_dist_ss->Write();
  tstar_dist_ss->Write();
  rstarPRF_dist_os->Write();
  tstar_dist_os->Write();
  NumLambda1->Write();
  DenLambda1->Write();
  NumLambda2_ss->Write();
  NumLambda2_os->Write();
  DenLambda2_ss->Write();
  DenLambda2_os->Write();
  NumLambda3_ss->Write();
  DenLambda3_ss->Write();
  NumLambda3_os->Write();
  DenLambda3_os->Write();
  Den_ss->Write();
  Den_os->Write();
  DenEM_ss->Write();
  DenEM_os->Write();
  Num_Cos_ss->Write();
  Num_Cos_os->Write();
  Num_CosFSI_ss->Write();
  Num_CosFSI_os->Write();
  NumSq_CosFSI_ss->Write();
  NumSq_CosFSI_os->Write();
  LargeRpairs_ss->Write();
  LargeRpairs_os->Write();
  //
  Num_PrimCosFSI_ss->Write();
  Num_PrimCosFSI_os->Write();
  //
  rstar3_dist_ss->Write();
  tstar3_dist_ss->Write();
  //
  K2_ss->Write();
  PlaneWF_ss->Write();
  K2_os->Write();
  PlaneWF_os->Write();
  //
  K3os->Write();
  PlaneWF3os->Write();
  K3os_3D->Write();
  PlaneWF3os_3D->Write();
  K3ss->Write();
  PlaneWF3ss->Write();
  K3ss_3D->Write();
  PlaneWF3ss_3D->Write();
  //
  r3num->Write();
  r3den1->Write();
  r3den2->Write();
  r3den3->Write();
  r3numSq->Write();
  r3den1Sq->Write();
  r3den2Sq->Write();
  r3den3Sq->Write();
  r3numEn->Write();
  r3den1En->Write();
  r3den2En->Write();
  r3den3En->Write();
  r3numME->Write();
  r3den1ME->Write();
  r3den2ME->Write();
  r3den3ME->Write();
  //
  Num_Cos12->Write();
  Num_Cos23->Write();
  Num_Cos31->Write();
  Den_Cos12->Write();
  Den_Cos23->Write();
  Den_Cos31->Write();
  //
  test->Write();
  //
  outfile->Close();

  seconds = time(0) - seconds;
  cout<<"Minutes = "<<seconds/60.<<endl;


}
float Gamov(float eta){
  return (2*PI*eta/(exp(2*PI*eta)-1));
}
float fact(float n){
  return (n < 1.00001) ? 1 : fact(n - 1) * n;
}
void Shuffle(Int_t *iarr, Int_t i1, Int_t i2)
{
  Int_t j, k;
  Int_t a = i2 - i1;
  for (Int_t i = i1; i < i2+1; i++) {
    j = (Int_t) (gRandom->Rndm() * a);
    k = iarr[j];
    iarr[j] = iarr[i];
    iarr[i] = k;
  }
}
void BoostPRF(float P1[4], float P2[4], float V[4], TVector3 *T){
  // P1 is particle 1 momentum in reaction CMS (0=E,1=px,2=py,3=pz)
  // P2 is particle 2 momentum in reaction CMS (0=E,1=px,2=py,3=pz)
  // V is the vector to be transformed (V is in CMS)
  // T is the PRF vector

  float MT=sqrt(pow(P1[0]+P2[0],2)-pow(P1[3]+P2[3],2));
  float PT=sqrt(pow(P1[1]+P2[1],2)+pow(P1[2]+P2[2],2));
  float MINV=sqrt(pow(P1[0]+P2[0],2) - pow(P1[1]+P2[1],2) - pow(P1[2]+P2[2],2) - pow(P1[3]+P2[3],2));
  // LCMS
  T->SetX( ((P1[1]+P2[1])*V[1] + (P1[2]+P2[2])*V[2])/PT );
  T->SetY( ((P1[1]+P2[1])*V[2] - (P1[2]+P2[2])*V[1])/PT );
  T->SetZ( ((P1[0]+P2[0])*V[3] - (P1[3]+P2[3])*V[0])/MT );
  // PRF
  T->SetX( T->X()*MINV/MT - PT/(MT*MINV)*((P1[0]+P2[0])*V[0] - (P1[1]+P2[1])*V[1] - (P1[2]+P2[2])*V[2] - (P1[3]+P2[3])*V[3]) );
}
// h function from Lednicky, phys. of part. and nuclei 40:307 (2009)
double hFunction(double eta){
  double DS=1.0;
  double N=0;
  double S=0;
  double returnValue=0;
  if(eta < 3.0){
    while(DS > 1.0e-13){
      N++;
      DS = 1/(N*(pow(N/eta,2)+1.0));
      S += DS;
    }
    returnValue = (S - EulerC - log(fabs(eta)));
  }else{// small kstar, high eta
    double etaSquared=pow(eta,2);
    double etaPowered=etaSquared;
    for(int ii=0; ii<9; ii++){// 9 was maximum value for Lednicky's code
      returnValue += 1/etaPowered * hcoef[ii];
      etaPowered *= etaSquared;
    }
  }
  return returnValue; 
}
double KinverseFunction(double kstarSq, int J){
  double ESq = kstarSq + pow(massPOI,2);
  double E = sqrt(ESq);
  double xSq = kstarSq/pow(massPOI,2);
  double Kinverse = E*(4*ESq-fzero[J][4])/(4*pow(massPOI,2)-fzero[J][4]);
  Kinverse /= fzero[J][0] + fzero[J][1]*xSq + fzero[J][2]*pow(xSq,2) + fzero[J][3]*pow(xSq,3);
  return Kinverse;
}

