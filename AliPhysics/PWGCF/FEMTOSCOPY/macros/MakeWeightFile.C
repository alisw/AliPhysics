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

#define BohrR 1963.6885 // Bohr Radius for pions
#define FmToGeV 0.19733 // conversion to fm
#define PI 3.1415926
#define masspiC 0.1395702 // pi+ mass (GeV/c^2)

using namespace std;

void SetFSICorrelations();
void SetMomResCorrections();
double Gamov(int, double);
float FSICorrelation(float);
void SetMuonCorrections();

TH1D *fFSIss[12];
TH2D *fMomResC2SC;
TH2D *fWeightmuonCorrection;

int fFSIindex=0;
short CollisionType=2;// 0(PbPb), 1(pPb), 2(pp)
bool FSIandMRandMuonCorrect=1;
double TwoFrac=0.7;
const int NumEDbins=2;
bool EDBinning=kFALSE;

void MakeWeightFile()
{
  
  gStyle->SetOptStat(0);
  gStyle->SetOptDate(0);
  gStyle->SetOptFit(1);
  
  SetFSICorrelations();
  SetMomResCorrections();
  SetMuonCorrections();

  //
  TFile *InputFile;
  if(CollisionType==0){
    //TFile *InputFile = new TFile("Results/RawWeightFile_11h.root","READ");
    //TFile *InputFile = new TFile("Results/PDC_11h_pT_0p2to1p0_FullRunWrongWeightsNoPadRowTTCandMCTTC_RawWeightFile.root","READ");
    //TFile *InputFile = new TFile("Results/RawWeightFile_11h_0p02eta0p045phi_0p03eta0p067phi.root","READ");
    InputFile = new TFile("Results/RawWeightFile_11h_7kT.root","READ");
    //InputFile = new TFile("Results/PDC_11h_0p06to0p08Norm_RawWeightFile100kTbins.root","READ");
    //TFile *InputFile = new TFile("Results/RawWeightFile_11h_q2Binning_LowPtMultBinning.root","READ");
    //TFile *InputFile = new TFile("Results/PDC_11h_Norm0p06to0p08_RawWeightFileLowpTBinningWithHighpTConstraint.root","READ");
    //TFile *InputFile = new TFile("Results/RawWeightFile_11h_8MeVbins.root","READ");
    //TFile *InputFile = new TFile("Results/RawWeightFile_11h_81EMbins.root","READ");
  }else if(CollisionType==1){
    InputFile = new TFile("Results/RawWeightFile_13bc_6kT.root","READ");
  }else{
    InputFile = new TFile("Results/RawWeightFile_10bcde_6kT.root","READ");
  }

  double NormL=0.135;// 0.06 or 0.135
  double NormH=0.2;// 0.08 or 0.2
  double NormL_1D=0.15;
  double NormH_1D=0.2;
  if(CollisionType!=0){
    NormL=0.3;// 0.3, 0.18
    NormH=0.35;// 0.35, 0.22
    NormL_1D=0.7;
    NormH_1D=0.8;
  }
  

  TDirectoryFile *tdir = (TDirectoryFile*)InputFile->Get("PWGCF.outputFourPionAnalysis.root");
  TList *MyList=(TList*)tdir->Get("FourPionOutput_1");
  //TList *MyList=(TList*)InputFile->Get("MyList");
  InputFile->Close();
  
  TH1D *Events = (TH1D*)MyList->FindObject("fEvents2");
  cout<<Events->GetBinContent(1)<<endl;
  
  
  const int KtBins=6;//4
  const int KtBinsOneD=28;
  const int KyBins=1;
  const int MBins=10;
  const int MBLimit=10;// 0-10. 1 for GenSignal, 10 for the rest!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  double fKstepT[KtBins]={0};
  double fKmeanT[KtBins]={0};
  double fKmiddleT[KtBins]={0};

  // 7x1 (Kt: 0-0.2, 0.2-0.25, 0.25-0.3, 0.3-0.35, 0.35-0.4, 0.4-0.45, 0.45-1.0)
  if(KtBins==7){
    fKstepT[0] = 0.2; fKstepT[1] = 0.05; fKstepT[2] = 0.05; fKstepT[3] = 0.05; fKstepT[4] = 0.05; fKstepT[5] = 0.05; fKstepT[6] = 0.55;
    fKmeanT[0] = 0.188; fKmeanT[1] = 0.227; fKmeanT[2] = 0.275; fKmeanT[3] = 0.324; fKmeanT[4] = 0.374; fKmeanT[5] = 0.424; fKmeanT[6] = 0.552; 
    fKmiddleT[0] = 0.1; fKmiddleT[1] = 0.225; fKmiddleT[2] = 0.275; fKmiddleT[3] = 0.325; fKmiddleT[4] = 0.375; fKmiddleT[5] = 0.425; fKmiddleT[6] = 0.725;
  }
  // 6x1 (Kt: 0-0.2, 0.2-0.24, 0.24-0.3, 0.3-0.35, 0.35-0.45, 0.45-1.0)
  if(KtBins==6){
    fKstepT[0] = 0.2; fKstepT[1] = 0.04; fKstepT[2] = 0.06; fKstepT[3] = 0.05; fKstepT[4] = 0.1; fKstepT[5] = 0.55;
    fKmeanT[0] = 0.188; fKmeanT[1] = 0.222; fKmeanT[2] = 0.270; fKmeanT[3] = 0.324; fKmeanT[4] = 0.395; fKmeanT[5] = 0.551; 
    fKmiddleT[0] = 0.1; fKmiddleT[1] = 0.22; fKmiddleT[2] = 0.27; fKmiddleT[3] = 0.325; fKmiddleT[4] = 0.4; fKmiddleT[5] = 0.725;
  }
  // 4x1 (Kt: 0-0.25, 0.25-0.35, 0.35-0.45, 0.45-1.0)
  if(KtBins==4){
    fKstepT[0] = 0.25; fKstepT[1] = 0.1; fKstepT[2] = 0.1; fKstepT[3] = 0.55;
    fKmeanT[0] = 0.212; fKmeanT[1] = 0.299; fKmeanT[2] = 0.398; fKmeanT[3] = 0.576;
    fKmiddleT[0] = 0.125; fKmiddleT[1] = 0.3; fKmiddleT[2] = 0.4; fKmiddleT[3] = 0.725;
  }
  // OneD kT bin edges
  int BinEdgesL[KtBinsOneD]={1,5,6,7,8,10,  10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10};
  int BinEdgesH[KtBinsOneD]={4,5,6,7,9,20,  20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20};
  /*int BinEdgesL[KtBinsOneD]={0};
  int BinEdgesH[KtBinsOneD]={0};
  for(int i=0; i<KtBinsOneD; i++){
    BinEdgesL[i] = 3*i + 17;
    BinEdgesH[i] = BinEdgesL[i] + 2;
    //cout<<BinEdgesL[i]<<"  "<<BinEdgesH[i]<<endl;
    }*/
  
  TFile *OutFile = new TFile("WeightFile_temp.root","RECREATE");
  TH3F *WeightHistos[KtBins][MBins][NumEDbins];
  TH2F *WeightHistosOneD[MBins];
 
  TH3D *QinvMean = (TH3D*)MyList->FindObject("TwoParticle_Charge1_0_Charge2_0_M_0_ED_0_Term_1_osl_b1_QW");
  TH3D *QinvMeanDen = (TH3D*)MyList->FindObject("TwoParticle_Charge1_0_Charge2_0_M_0_ED_0_Term_1_osl_b1");
  QinvMean->Divide(QinvMeanDen);

 
  for(int ktB=1; ktB<=KtBinsOneD; ktB++){
    cout<<"kT Bin "<<ktB<<endl;
    for(int MB=1; MB<=MBins; MB++){
      if(CollisionType!=0 && MB>1) continue;
      //
      if(CollisionType==0){
	if(MB==1) fFSIindex = 0;//0-5%
	else if(MB==2) fFSIindex = 1;//5-10%
	else if(MB<=3) fFSIindex = 2;//10-20%
	else if(MB<=6) fFSIindex = 3;//20-30%
	else if(MB<=8) fFSIindex = 4;//30-40%
	else fFSIindex = 5;//40-50%
      }else if(CollisionType==1) fFSIindex=9;
      else fFSIindex=11;
      
      //
      Int_t rBinForTPNMomRes = 10;
      if(CollisionType==0){
	if(MB==1) {rBinForTPNMomRes=10;}// 10 fm with EW (fRMax should be 11 for normal running)
	else if(MB==2) {rBinForTPNMomRes=9;}
	else if(MB<=4) {rBinForTPNMomRes=8;}
	else if(MB<=6) {rBinForTPNMomRes=7;}
	else {rBinForTPNMomRes=6;}
      }else rBinForTPNMomRes=2;
      //
      
      for(int q2B=1; q2B<=NumEDbins; q2B++){
	TString *InNameNum = new TString("TPN_num_Kt_");
	*InNameNum += ktB-1;
	InNameNum->Append("_Ky_0_M_");
	if(MB<=MBLimit) *InNameNum += MB-1;// make sure MB assignment is correct!!!!!!!!!!!!!!!
	else *InNameNum += 1;
	InNameNum->Append("_ED_");
	if(EDBinning) *InNameNum += q2B-1;
	else *InNameNum += 0;// for non q2 binning
	//
	TString *InNameDen = new TString("TPN_den_Kt_");
	*InNameDen += ktB-1;
	InNameDen->Append("_Ky_0_M_");
	if(MB<=MBLimit) *InNameDen += MB-1;// make sure MB assignment is correct!!!!!!!!!!!!!!!
	else *InNameDen += 1;
	InNameDen->Append("_ED_");
	if(EDBinning) *InNameDen += q2B-1;
	else *InNameDen += 0;// for non q2 binning
	//
	TString *OneDNameNum = new TString("TwoParticle_Charge1_0_Charge2_0_M_");//TwoParticle_Charge1_0_Charge2_0_M_0_ED_0_Term_1
	if(MB<=MBLimit) *OneDNameNum += MB-1;// make sure MB assignment is correct!!!!!!!!!!!!!!!
	else *OneDNameNum += 1;
	OneDNameNum->Append("_ED_0_Term_1");
	//
	TString *OneDNameDen = new TString("TwoParticle_Charge1_0_Charge2_0_M_");//TwoParticle_Charge1_0_Charge2_0_M_0_ED_0_Term_1
	if(MB<=MBLimit) *OneDNameDen += MB-1;// make sure MB assignment is correct!!!!!!!!!!!!!!!
	else *OneDNameDen += 1;
	OneDNameDen->Append("_ED_0_Term_2");
	//
	//
	TH2D *tempNum2D = (TH2D*)MyList->FindObject(OneDNameNum->Data());
	TH2D *tempDen2D = (TH2D*)MyList->FindObject(OneDNameDen->Data());
	TH1D *tempNum1D = (TH1D*)tempNum2D->ProjectionY("tempNum1D",BinEdgesL[ktB-1], BinEdgesH[ktB-1]);
	TH1D *tempDen1D = (TH1D*)tempDen2D->ProjectionY("tempDen1D",BinEdgesL[ktB-1], BinEdgesH[ktB-1]);
	//
	int NormBinStart1D = tempNum1D->GetXaxis()->FindBin(NormL_1D);
	int NormBinEnd1D = tempNum1D->GetXaxis()->FindBin(NormH_1D);
	double Norm1D = tempNum1D->Integral(NormBinStart1D, NormBinEnd1D);
	if(tempDen1D->Integral(NormBinStart1D, NormBinEnd1D) > 0){
	  Norm1D /= tempDen1D->Integral(NormBinStart1D, NormBinEnd1D);
	}else Norm1D=0;
	cout<<"1D Norm = "<<Norm1D<<endl;
	//for(int i=0; i<KtBinsOneD; i++){
	//cout<<tempNum2D->GetXaxis()->GetBinLowEdge(BinEdgesL[i])<<"  "<<tempNum2D->GetXaxis()->GetBinUpEdge(BinEdgesH[i])<<endl;
	//}	
	//
	//
	if(ktB<=KtBins){
	  TH3D *tempNum = (TH3D*)MyList->FindObject(InNameNum->Data());
	  TH3D *tempDen = (TH3D*)MyList->FindObject(InNameDen->Data());
	  int NormBinStartOut = tempNum->GetXaxis()->FindBin(NormL);
	  int NormBinStartSideLong = tempNum->GetXaxis()->FindBin(NormL);
	  int NormBinEnd = tempNum->GetXaxis()->FindBin(NormH);
	  double Norm = tempNum->Integral(NormBinStartOut,NormBinEnd, NormBinStartSideLong,NormBinEnd, NormBinStartSideLong,NormBinEnd);
	  if(tempDen->Integral(NormBinStartOut,NormBinEnd, NormBinStartSideLong,NormBinEnd, NormBinStartSideLong,NormBinEnd) > 0){
	    Norm /= tempDen->Integral(NormBinStartOut,NormBinEnd, NormBinStartSideLong,NormBinEnd, NormBinStartSideLong,NormBinEnd);
	  }else Norm=0;
	  cout<<"Mbin = "<<MB<<"  Normalization = "<<Norm<<"  Num Count = "<<tempNum->Integral(NormBinStartOut,NormBinEnd, NormBinStartSideLong,NormBinEnd, NormBinStartSideLong,NormBinEnd)<<"  Den Count = "<<tempDen->Integral(NormBinStartOut,NormBinEnd, NormBinStartSideLong,NormBinEnd, NormBinStartSideLong,NormBinEnd)<<endl;
	}
	//
	TString *OutNameWeight = new TString("Weight_Kt_");
	*OutNameWeight += ktB-1;
	OutNameWeight->Append("_Ky_0_M_");
	*OutNameWeight += MB-1;
	OutNameWeight->Append("_ED_");
	*OutNameWeight += q2B-1;
	TString *OutNameWeight1D = new TString("Weight_M_");
	*OutNameWeight1D += MB-1;
	OutNameWeight1D->Append("_1D");
	//
	int Nbins=tempNum->GetNbinsX();
	double QLimit = tempNum->GetXaxis()->GetBinUpEdge(Nbins);
	if(ktB<=KtBins){
	  WeightHistos[ktB-1][MB-1][q2B-1] = new TH3F(OutNameWeight->Data(),"r3 Weights", Nbins,0,QLimit, Nbins,0,QLimit, Nbins,0,QLimit);
	  WeightHistos[ktB-1][MB-1][q2B-1]->GetXaxis()->SetTitle("out");
	  WeightHistos[ktB-1][MB-1][q2B-1]->GetYaxis()->SetTitle("side");
	  WeightHistos[ktB-1][MB-1][q2B-1]->GetZaxis()->SetTitle("long");
	  WeightHistos[ktB-1][MB-1][q2B-1]->GetXaxis()->SetTitleOffset(1.8);
	  WeightHistos[ktB-1][MB-1][q2B-1]->GetYaxis()->SetTitleOffset(1.8);
	  WeightHistos[ktB-1][MB-1][q2B-1]->GetZaxis()->SetTitleOffset(1.8);
	}
	//
	// 2D weights
	//
	if(q2B==1) {
	  if(ktB==1) WeightHistosOneD[MB-1] = new TH2F(OutNameWeight1D->Data(),"1D Weights", KtBinsOneD,0.16,1.0, 100, 0, 0.5);
	  for( int invB=1; invB<=100; invB++){
	    double weight=1, weight_e=0;
	    if(tempDen1D->GetBinContent(invB) > 0 && tempNum1D->GetBinContent(invB) > 0) {
	      weight = double(tempNum1D->GetBinContent(invB))/double(tempDen1D->GetBinContent(invB)) / Norm1D;
	      if(FSIandMRandMuonCorrect){
		float qinv = tempNum1D->GetBinCenter(invB);
		int momBin = fMomResC2SC->GetYaxis()->FindBin(qinv);
		double FSICorr = FSICorrelation(qinv);
		double MomResCorr = fMomResC2SC->GetBinContent(rBinForTPNMomRes, momBin);
		double MuonCorr = fWeightmuonCorrection->GetBinContent(rBinForTPNMomRes, momBin);
		weight = weight*MomResCorr - TwoFrac*FSICorr - (1-TwoFrac);
		weight /= FSICorr*TwoFrac;
		weight *= MuonCorr;
		weight += 1;
	      }
	      weight_e = pow(sqrt(double(tempNum1D->GetBinContent(invB)))/double(tempDen1D->GetBinContent(invB)) / Norm1D,2);
	      weight_e += pow(sqrt(double(tempDen1D->GetBinContent(invB)))*double(tempNum1D->GetBinContent(invB))/pow(double(tempDen1D->GetBinContent(invB)),2) / Norm1D,2);
	      weight_e = sqrt(weight_e);
	      WeightHistosOneD[MB-1]->SetBinContent(ktB, invB, weight-1);// difference from unity
	      WeightHistosOneD[MB-1]->SetBinError(ktB, invB, weight_e);
	    }
	  }
	  // Set 0 entry low q bins to nearest non-zero neighbor
	  for(int b=5; b>0; b--) {
	    if( WeightHistosOneD[MB-1]->GetBinContent(ktB, b) == 0){
	      WeightHistosOneD[MB-1]->SetBinContent(ktB, b, WeightHistosOneD[MB-1]->GetBinContent(ktB, b+1));
	    }
	  }
	  //
	  if(ktB==KtBinsOneD) WeightHistosOneD[MB-1]->Write();
	}
	//
	//  4D weights
	//
	if(ktB>KtBins) continue;
	//
	double MaxQout = 2*(fKmiddleT[ktB-1]+fKstepT[ktB-1]/2.) - 2*0.16;
	int SaturationBin = WeightHistos[ktB-1][MB-1][q2B-1]->GetXaxis()->FindBin(MaxQout) - 1;
	//
	double LowQcount=0, HighQcount=0;
	for(int outB=1; outB<=Nbins; outB++){
	  for(int sideB=1; sideB<=Nbins; sideB++){
	    for(int longB=1; longB<=Nbins; longB++){
	      if(Norm==0) continue;
	      double weight=1, weight_e=0;
	      if(tempDen->GetBinContent(outB,sideB,longB) > 0 && tempNum->GetBinContent(outB,sideB,longB) > 0) {
		
		weight = double(tempNum->GetBinContent(outB,sideB,longB))/double(tempDen->GetBinContent(outB,sideB,longB)) / Norm;
		if(FSIandMRandMuonCorrect){
		  float qinv = QinvMean->GetBinContent(outB,sideB,longB);
		  if(qinv==0){// never actually used in Pb-Pb
		    float qout = QinvMean->GetXaxis()->GetBinCenter(outB);
		    float gamma = sqrt(pow(0.139,2) + pow(fKmeanT[ktB-1],2)) / 0.139;
		    qout /= gamma;
		    qinv = sqrt(pow(qout,2) + pow(QinvMean->GetYaxis()->GetBinCenter(sideB),2) + pow(QinvMean->GetZaxis()->GetBinCenter(longB),2));
		  }
		  
		  int momBin = fMomResC2SC->GetYaxis()->FindBin(qinv);
		  double FSICorr = FSICorrelation(qinv);
		  double MomResCorr = fMomResC2SC->GetBinContent(rBinForTPNMomRes, momBin);
		  double MuonCorr = fWeightmuonCorrection->GetBinContent(rBinForTPNMomRes, momBin);
		  weight = weight*MomResCorr - TwoFrac*FSICorr - (1-TwoFrac);
		  weight /= FSICorr*TwoFrac;
		  weight *= MuonCorr;
		  weight += 1;
		  //
		  //if(MB==1 && ktB==5 && longB>20 && weight==1) cout<<outB<<"  "<<sideB<<"  "<<longB<<"  "<<qinv<<"    "<<weight<<"       "<<FSICorr<<"  "<<MuonCorr<<"  "<<MomResCorr<<endl;
		}
		weight_e = pow(sqrt(double(tempNum->GetBinContent(outB,sideB,longB)))/double(tempDen->GetBinContent(outB,sideB,longB)) / Norm,2);
		weight_e += pow(sqrt(double(tempDen->GetBinContent(outB,sideB,longB)))*double(tempNum->GetBinContent(outB,sideB,longB))/pow(double(tempDen->GetBinContent(outB,sideB,longB)),2) / Norm,2);
		weight_e = sqrt(weight_e);
		
		double qo = tempNum->GetXaxis()->GetBinCenter(outB);
		double qs = tempNum->GetYaxis()->GetBinCenter(sideB);
		double ql = tempNum->GetZaxis()->GetBinCenter(longB);
		double qmag = sqrt(pow(qo,2) + pow(qs,2) + pow(ql,2));
		if(qmag > 0.01 && qmag < 0.04) LowQcount++;
		if(qmag > 0.04 && qmag < 0.07) HighQcount++;
	      }
	      if(outB > SaturationBin){
		WeightHistos[ktB-1][MB-1][q2B-1]->SetBinContent(outB,sideB,longB, WeightHistos[ktB-1][MB-1][q2B-1]->GetBinContent(SaturationBin,sideB,longB));
		//WeightHistos[ktB-1][MB-1][q2B-1]->SetBinContent(outB,sideB,longB, weight-1);// testing
		WeightHistos[ktB-1][MB-1][q2B-1]->SetBinError(outB,sideB,longB, 0);
	      }else if(weight==0){
		WeightHistos[ktB-1][MB-1][q2B-1]->SetBinContent(outB,sideB,longB, 0);
		WeightHistos[ktB-1][MB-1][q2B-1]->SetBinError(outB,sideB,longB, 0);
	      }else {
		WeightHistos[ktB-1][MB-1][q2B-1]->SetBinContent(outB,sideB,longB, weight-1);// difference from unity
		WeightHistos[ktB-1][MB-1][q2B-1]->SetBinError(outB,sideB,longB, weight_e);
	      }
	      
	      
	    }
	  }
	}
	//cout<<"PairCount   "<<LowQcount<<"  "<<HighQcount<<endl;
	WeightHistos[ktB-1][MB-1][q2B-1]->Write();
	
      }// q2B
    }// MB
  }// ktB

 

  int BOI1=2;
  int BOI2=2;
  //
  int KTBOI=0;
  //
  TH3D *histoOI = WeightHistos[KTBOI][0][0]->Clone();
  histoOI->SetDirectory(0);
    
  TH1D *pro = WeightHistos[KTBOI][0][0]->ProjectionX("pro",BOI1,BOI1,BOI2,BOI2);
  pro->Draw();

  TH1D *pro2 = WeightHistos[KTBOI][0][0]->ProjectionY("pro2",BOI1,BOI1,BOI2,BOI2);
  pro2->SetLineColor(2);
  //pro2->Draw("same");

  TH1D *pro3 = WeightHistos[KTBOI][0][0]->ProjectionZ("pro3",BOI1,BOI1,BOI2,BOI2);
  pro3->SetLineColor(4);
  //pro3->Draw("same");

  //for(int bin=1; bin<=40; bin++) cout<<pro->GetBinContent(bin)<<", ";
  //cout<<endl;
  //cout<<endl;
  //for(int bin=1; bin<=40; bin++) cout<<pro->GetBinError(bin)<<", ";
  //cout<<endl;

  // qout ,2--2 other bin integration (lowest kT, M0)
  //double Ref[40]={0.186581, 0.180099, 0.16335, 0.137551, 0.103902, 0.0754772, 0.0381106, 0.0250045, 0.0194101, 0.00606144, 0.00175905, -0.00335014, 0.00347716, -0.0020038, 0.00276897, -0.0114744, -0.00412381, 0.00880619, -0.00533615, 0.000337426, -0.00603902, 0.00860179, -0.00154849, 0.00605317, -0.00328388, -0.0257134, 0.00816, 0.00588766, 0.00238239, -0.0208381, 0.00262856, 0.000531408, -9.19566e-05, -0.011538, -0.0216935, -0.0439466, 0, 0, 0, 0};
  //double Ref_e[40]={0.00258484, 0.00262131, 0.00265747, 0.00268686, 0.00271356, 0.00276084, 0.00277321, 0.00281418, 0.00285723, 0.00284373, 0.00286302, 0.00289162, 0.00297097, 0.00304969, 0.00316452, 0.00323752, 0.00339239, 0.0035809, 0.00369192, 0.00389189, 0.00406313, 0.00436902, 0.00457948, 0.00493212, 0.00524766, 0.00552003, 0.00621051, 0.0068126, 0.00754255, 0.00825506, 0.00960699, 0.0111344, 0.0131076, 0.0162349, 0.021695, 0.0372195, 0, 0, 0, 0};
  // qside, 2--2 other bin integration
  //double Ref[40]={0.180498, 0.180099, 0.185806, 0.163894, 0.132356, 0.098442, 0.0732086, 0.0419, 0.0273973, 0.0140072, 0.00491044, 0.00828187, -0.00346015, 0.00556346, 0.00280363, 0.000493745, 0.00483145, -0.000884939, 0.00825155, -0.00551758, -0.00208245, -0.00369946, -0.000581864, -0.000484325, -0.000560278, 0.00125169, 0.00402387, 0.00369215, 0.000496307, 0.00343727, 0.00382731, -0.00107187, -0.00187021, -0.00498415, 0.00312242, -0.000784691, 0.00153594, 0.00514852, -0.00480193, 0.00576707};
  //double Ref_e[40]={0.00263533, 0.00262131, 0.00262328, 0.00255935, 0.00250353, 0.00244688, 0.00240138, 0.00234001, 0.00231244, 0.00228984, 0.00226491, 0.00227351, 0.00224481, 0.00226501, 0.00224979, 0.00223847, 0.00224863, 0.0022289, 0.00224206, 0.00220784, 0.0022124, 0.00219928, 0.00219753, 0.0021926, 0.00218304, 0.00218404, 0.00217921, 0.00217676, 0.00216391, 0.00216305, 0.00215278, 0.00213398, 0.00212274, 0.00211043, 0.0021198, 0.00210246, 0.00210134, 0.00209683, 0.00207121, 0.00208596};
  // qlong, 2--2 other bin integration
  //double Ref[40]={0.19068, 0.180099, 0.16035, 0.118137, 0.0800825, 0.0484194, 0.0314626, 0.0140836, 0.015328, 0.00456526, 0.00461427, 0.00444591, -0.00385123, 0.0076453, 0.00418319, 0.00385833, 0.00192732, 0.0015057, 0.00151401, 0.00441839, 0.00255238, -0.00562173, -0.00434487, 0.00327595, -0.000923301, 0.00340199, -0.00116685, -0.00639669, -0.000535647, 0.000632789, -0.00431929, -0.00590313, -0.00334011, -0.00415945, 0.0113907, -0.00252384, 0.00152043, 0.000260006, -0.00476303, -0.00232199};
  //double Ref_e[40]={0.00394138, 0.00262131, 0.00257243, 0.00251757, 0.00247518, 0.00244153, 0.00243551, 0.00242434, 0.00245618, 0.00245905, 0.00248088, 0.00250852, 0.00251633, 0.00257305, 0.0025939, 0.00261776, 0.0026447, 0.00267184, 0.00269996, 0.00273682, 0.00275901, 0.00277727, 0.00280127, 0.00286128, 0.00288497, 0.00292533, 0.00294705, 0.00296652, 0.00302051, 0.00306191, 0.00308495, 0.00311484, 0.00315736, 0.00320054, 0.00328945, 0.00328307, 0.00333872, 0.00337698, 0.0034125, 0.00346363};
  //
  // qout ,2--2 other bin integration (highest kT, M0)
  //double Ref[40]={0.347024, 0.340981, 0.207019, 0.204354, 0.235698, 0.113875, 0.103287, 0.0390906, 0.049902, 0.0632777, 0.0196594, 0.0168798, -0.0258711, -0.0173084, -0.0255809, -0.03879, -0.0467173, -0.0285905, -0.0280937, -0.0144413, -0.0242922, -0.0183046, -0.0232596, -0.0242507, -0.0322753, -0.0337653, -0.0270012, -0.00878507, -0.00345743, -0.0201966, -0.023139, -0.0138009, -0.0236573, -0.0212273, -0.0261684, -0.014746, 0.00113992, -0.00296493, -0.00872464, -0.0149615};
  //double Ref_e[40]={0.0291587, 0.028769, 0.0247413, 0.0233308, 0.0228364, 0.0193632, 0.0182532, 0.0163304, 0.0156345, 0.0151442, 0.01384, 0.0133412, 0.0123078, 0.0116285, 0.0108021, 0.00981288, 0.00906128, 0.00856431, 0.00804651, 0.00772649, 0.00733884, 0.00705491, 0.0068498, 0.00661994, 0.00638918, 0.00626731, 0.00617582, 0.00618696, 0.00611707, 0.00593496, 0.00581311, 0.00576991, 0.00562799, 0.00553088, 0.00539531, 0.00533839, 0.00533781, 0.00523377, 0.00512349, 0.00503239};
  // qout ,2--2 other bin integration (kT1, M9)
  //double Ref[40]={0.274214, 0.266342, 0.262087, 0.258269, 0.342804, 0.271642, 0.262439, 0.191731, 0.161266, 0.0480234, 0.0678067, 0.0298723, 0.140997, 0.0225363, 0.0369282, 0.0474108, 0.00666558, -0.00719739, 0.0397449, -0.0154514, 0.0563656, -0.0283267, -0.0233019, 0.00502379, -0.0142704, 0.00171802, 0.00224826, -0.0103424, -0.00215456, -0.00290923, 0.0294002, 0.00502061, -0.0311593, -0.00107628, 0.0157102, 0.0242639, -0.00358242, 0.0273829, -0.00498812, -0.0191332};
  //double Ref_e[40]={0.0219786, 0.0219724, 0.0219296, 0.021582, 0.023235, 0.0219217, 0.0219606, 0.0209754, 0.0203662, 0.0185691, 0.0189119, 0.0183129, 0.0204752, 0.017965, 0.0179198, 0.0177246, 0.0170545, 0.0166412, 0.0172614, 0.0163279, 0.0174978, 0.016035, 0.0160211, 0.0165804, 0.0164245, 0.016433, 0.0167745, 0.0163847, 0.0166272, 0.0166435, 0.0172967, 0.0171336, 0.0165754, 0.0172471, 0.0177293, 0.0181375, 0.018193, 0.0189333, 0.0187681, 0.0189596};
  
  /*TH1D *Ratio = (TH1D*)pro->Clone();
  Ratio->SetTitle("C_{2}^{#pm#pm} ratio of nanoAOD160 to AOD160");
  Ratio->GetXaxis()->SetTitleOffset(0.9);
  Ratio->GetXaxis()->SetTitle("q_{out} GeV/c");
  for(int bin=1; bin<=40; bin++) {
    double value = (pro->GetBinContent(bin)+1) / (Ref[bin-1]+1);
    double value_e = sqrt( fabs(pow(pro->GetBinError(bin),2) - pow(Ref_e[bin-1],2)));
    Ratio->SetBinContent(bin, value);
    Ratio->SetBinError(bin, value_e);
  }
  Ratio->Draw();
  Ratio->Fit("pol0","","",0,0.2);
  */

  /*TFile *projections=new TFile("projections.root","RECREATE");
  histoOI->GetXaxis()->SetRange(BOI1,BOI1);
  TH2D *proYZ = (TH2D*)histoOI->Project3D("yz");
  histoOI->GetXaxis()->SetRange(1,40);
  histoOI->GetYaxis()->SetRange(BOI1,BOI1);
  TH2D *proXZ = (TH2D*)histoOI->Project3D("xz");
  histoOI->GetYaxis()->SetRange(1,40);
  histoOI->GetZaxis()->SetRange(BOI1,BOI1);
  TH2D *proXY = (TH2D*)histoOI->Project3D("xy");
  proYZ->GetXaxis()->SetRangeUser(0,0.1);
  proYZ->GetYaxis()->SetRangeUser(0,0.1);
  proXZ->GetXaxis()->SetRangeUser(0,0.1);
  proXZ->GetYaxis()->SetRangeUser(0,0.1);
  proXY->GetXaxis()->SetRangeUser(0,0.1);
  proXY->GetYaxis()->SetRangeUser(0,0.1);
  proXY->Draw("lego");
  
  proYZ->Write();
  proXZ->Write();
  proXY->Write();
  */

}
//________________________________________________________________________
void SetFSICorrelations(){
  // read in 2-particle and 3-particle FSI correlations = K2 & K3
  // 2-particle input histo from file is binned in qinv.  3-particle in qinv of each pair
  TFile *fsifile = new TFile("KFile.root","READ");
  if(!fsifile->IsOpen()) {
    cout<<"No FSI file found!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
  }else {cout<<"Good FSI File Found!"<<endl;}
  
  TH1D *temphistoSS[12];
  TH1D *temphistoOS[12];
  for(Int_t MB=0; MB<12; MB++) {
    TString *nameK2SS = new TString("K2ss_");
    *nameK2SS += MB;
    temphistoSS[MB] = (TH1D*)fsifile->Get(nameK2SS->Data());
    //
    fFSIss[MB] = (TH1D*)temphistoSS[MB]->Clone();
    fFSIss[MB]->SetDirectory(0);
  }
  //
  
  fsifile->Close();
  
}
//________________________________________________________________________
void SetMomResCorrections(){
 
  TFile *momResFile;
  if(CollisionType==0) momResFile = new TFile("MomResFile_FourPion.root","READ");
  else momResFile = new TFile("MomResFile_ppAndpPb.root","READ");
  if(!momResFile->IsOpen()) {
    cout<<"No momentum resolution file found"<<endl;
    AliFatal("No momentum resolution file found.  Kill process.");
  }else {cout<<"Good Momentum Resolution File Found!"<<endl;}
  
  TH2D *temp2DSC2 = (TH2D*)momResFile->Get("MRC_C2_SC");
  fMomResC2SC = (TH2D*)temp2DSC2->Clone();
  fMomResC2SC->SetDirectory(0);
  //
  momResFile->Close();
  
}
//________________________________________________________________________
void SetMuonCorrections(){
  
  TFile *MuonFile;
  if(CollisionType==0) MuonFile = new TFile("MuonCorrection_FourPion.root","READ");
  else MuonFile = new TFile("MuonCorrection_ppAndpPb.root","READ");
  if(!MuonFile->IsOpen()) {
    cout<<"No Muon file found"<<endl;
    AliFatal("No Muon file found.  Kill process.");
  }else {cout<<"Good Muon File Found!"<<endl;}
  
  fWeightmuonCorrection = (TH2D*)MuonFile->Get("WeightmuonCorrection");
  fWeightmuonCorrection->SetDirectory(0);
  //
  MuonFile->Close();
  
}
//________________________________________________________________________
double Gamov(int chargeProduct, double qinv){
  
  double arg = chargeProduct*2.*PI/(BohrR*qinv/2.);
  
  return arg/(exp(arg)-1);
}
//________________________________________________________________________
float FSICorrelation(float qinv){
  // returns 2-particle Coulomb correlations = K2
  Int_t qbinL = fFSIss[fFSIindex]->GetXaxis()->FindBin(qinv-fFSIss[fFSIindex]->GetXaxis()->GetBinWidth(1)/2.);
  Int_t qbinH = qbinL+1;
  if(qbinL <= 0) return 1.0;
  if(qbinH > fFSIss[fFSIindex]->GetNbinsX()) {
    Float_t ScaleFac = (fFSIss[fFSIindex]->GetBinContent(fFSIss[fFSIindex]->GetNbinsX()-1) - 1);
    ScaleFac /= (Gamov(1, fFSIss[fFSIindex]->GetXaxis()->GetBinCenter(fFSIss[fFSIindex]->GetNbinsX()-1)) - 1);
    return ( (Gamov(1, qinv)-1)*ScaleFac + 1);
  }else{
    Float_t slope=0;
    slope = fFSIss[fFSIindex]->GetBinContent(qbinL) - fFSIss[fFSIindex]->GetBinContent(qbinH);
    slope /= fFSIss[fFSIindex]->GetXaxis()->GetBinCenter(qbinL) - fFSIss[fFSIindex]->GetXaxis()->GetBinCenter(qbinH);
    return (slope*(qinv - fFSIss[fFSIindex]->GetXaxis()->GetBinCenter(qbinL)) + fFSIss[fFSIindex]->GetBinContent(qbinL));
  }
}
