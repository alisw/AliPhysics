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
#include "TCanvas.h"
#include "TPad.h"

#define BohrR 1963.6885
#define FmToGeV 0.19733 // conversion to fm
#define PI 3.1415926
#define masspiC 0.1395702 // pi+ mass (GeV/c^2)

using namespace std;

const int ChProdBOI=0;// 0=SameCharge, 1=MixedCharge1, 2=MixedCharge2
const int KTBin=0;// Kt3 bin. 0=low Kt3 bin.  1=high Kt3 bin
const int MBOI=0;// Centrality bin: 0-9
const int GbinPlot=int( (34) /2. ) + 55;// +5 (Rcoh=0), +25 (Rcoh=Rch) or +55 for extended G range
const int Q3binChi2= 4;// 2-5
const int Q4binChi2= 7;// 3-7

//
//
//
int TextFont=42;// 63, or 42
float SizeLabel=0.06;// 20(63 font), 0.08(42 font)
float SizeLegend=0.04;// .08
float SizeTitle=0.06;// 
float SizeSpecif=0.045;// 
float SF1=2/3.*0.95;
float SF2=1/2.*0.95;

double RightMargin=0.004;// 0.002
//

void Plot_plotsFourPion(){

  gStyle->SetOptStat(0);
  gStyle->SetOptDate(0);
  
  int Q3_meanpT[2][5]={{0,213,228,233,237},{0,316,336,354,366}};// low then high KT3
  int Q4_meanpT[2][7]={{0,0,221,229,234,238,241},{0,0,325,335,344,351,355}};// low then high KT4

  TFile *files[3][2][2][10];// SC/MC1/MC2, +/-, KT, MBINS
  //
 

  TH1D *C3QS[2][2][2][10];// SC/MC, +/-, KT, Mbin
  TH1D *c3QS[2][2][2][10];// SC/MC, +/-, KT, Mbin
  TH1D *C3QSBuilt[2][2][10];// +/-, KT, Mbin
  TH2D *C3QSBuilt2D[2][2][10];// +/-, KT, Mbin
  TH2D *C3QSNegBuilt2D[2][2][10];// +/-, KT, Mbin
  //
  TH1D *C4QS[3][2][2][10];// SC/MC, +/-, KT, Mbin
  TH1D *c4QS[3][2][2][10];// SC/MC, +/-, KT, Mbin
  TH1D *c4QSstage1[3][2][2][10];// SC/MC, +/-, KT, Mbin
  TH1D *c4QSstage2[2][2][10];// +/-, KT, Mbin
  TH1D *C4QSBuilt[2][2][10];// +/-, KT, Mbin
  TH2D *C4QSBuilt2D[2][2][10];// +/-, KT, Mbin
  TH2D *C4QSNegBuilt2D[2][2][10];// +/-, KT, Mbin
  //
  TH1D *r3[2][2][10];// +/-, KT, Mbin
  TH1D *r42[2][2][10];// +/-, KT, Mbin
  
  //
  //////////////////////////////

  // Start File access
  for(int mb=0; mb<10; mb++){
    for(int ChComb=0; ChComb<3; ChComb++) {// SC or MC1 or MC2
      for(int ch=0; ch<2; ch++) {// - or +
	for(int KT=0; KT<2; KT++) {// KT bin
	  if(ChComb==2 && ch!=0) continue;

	  TString *name = new TString("OutFiles/OutFile");
	  if(ChComb==0) name->Append("SC");
	  else if(ChComb==1) name->Append("MC1");
	  else name->Append("MC2");
	  //
	  if(ch==0) name->Append("_Neg_");
	  else name->Append("_Pos_");
	  name->Append("KT_");
	  *name += KT+1;
	  name->Append("_M");
	  if(mb<10) {*name += mb;}
	  else {*name += 0;}
	  name->Append(".root");
	  files[ChComb][ch][KT][mb] = new TFile(name->Data(),"READ");
	  ///////////////////////////////
	  if(ChComb!=2){
	    C3QS[ChComb][ch][KT][mb]=(TH1D*)files[ChComb][ch][KT][mb]->Get("C3QS");
	    c3QS[ChComb][ch][KT][mb]=(TH1D*)files[ChComb][ch][KT][mb]->Get("c3QS");
	    C3QS[ChComb][ch][KT][mb]->SetDirectory(0); c3QS[ChComb][ch][KT][mb]->SetDirectory(0); 
	    C3QS[ChComb][ch][KT][mb]->GetXaxis()->SetLabelFont(TextFont);  C3QS[ChComb][ch][KT][mb]->GetYaxis()->SetLabelFont(TextFont);
	    c3QS[ChComb][ch][KT][mb]->GetXaxis()->SetLabelFont(TextFont);  c3QS[ChComb][ch][KT][mb]->GetYaxis()->SetLabelFont(TextFont);
	    C3QS[ChComb][ch][KT][mb]->GetXaxis()->SetTitleFont(TextFont);  C3QS[ChComb][ch][KT][mb]->GetYaxis()->SetTitleFont(TextFont);
	    c3QS[ChComb][ch][KT][mb]->GetXaxis()->SetTitleFont(TextFont);  c3QS[ChComb][ch][KT][mb]->GetYaxis()->SetTitleFont(TextFont);
	    C3QS[ChComb][ch][KT][mb]->GetXaxis()->SetTitle("#font[12]{Q}_{3} (GeV/#font[12]{c})");
	    c3QS[ChComb][ch][KT][mb]->GetXaxis()->SetTitle("#font[12]{Q}_{3} (GeV/#font[12]{c})");
	    C3QS[ChComb][ch][KT][mb]->GetYaxis()->SetTitle("#font[12]{C}_{3} or #font[12]{#bf{c}}_{3}");
	    c3QS[ChComb][ch][KT][mb]->GetYaxis()->SetTitle("#font[12]{#bf{c}}_{3}");
	  }
	  C4QS[ChComb][ch][KT][mb]=(TH1D*)files[ChComb][ch][KT][mb]->Get("C4QS");
	  c4QS[ChComb][ch][KT][mb]=(TH1D*)files[ChComb][ch][KT][mb]->Get("c4QS");
	  c4QSstage1[ChComb][ch][KT][mb]=(TH1D*)files[ChComb][ch][KT][mb]->Get("c4QS_RemovalStage1");
	  C4QS[ChComb][ch][KT][mb]->SetDirectory(0); c4QS[ChComb][ch][KT][mb]->SetDirectory(0); 
	  c4QSstage1[ChComb][ch][KT][mb]->SetDirectory(0);
	  C4QS[ChComb][ch][KT][mb]->GetXaxis()->SetTitle("#font[12]{Q}_{4} (GeV/#font[12]{c})");
	  c4QS[ChComb][ch][KT][mb]->GetXaxis()->SetTitle("#font[12]{Q}_{4} (GeV/#font[12]{c})");
	  C4QS[ChComb][ch][KT][mb]->GetYaxis()->SetTitle("#font[12]{C}_{4} or #font[12]{#bf{c}}_{4}");
	  c4QS[ChComb][ch][KT][mb]->GetYaxis()->SetTitle("#font[12]{#bf{c}}_{4}");
	  
	  if(ChComb==0){
	    c4QSstage2[ch][KT][mb]=(TH1D*)files[ChComb][ch][KT][mb]->Get("c4QS_RemovalStage2");
	    C3QSBuilt[ch][KT][mb]=(TH1D*)files[ChComb][ch][KT][mb]->Get("C3QS_built");
	    C4QSBuilt[ch][KT][mb]=(TH1D*)files[ChComb][ch][KT][mb]->Get("C4QS_built");
	    C3QSBuilt2D[ch][KT][mb]=(TH2D*)files[ChComb][ch][KT][mb]->Get("C3QS_built2D");
	    C4QSBuilt2D[ch][KT][mb]=(TH2D*)files[ChComb][ch][KT][mb]->Get("C4QS_built2D");
	    C3QSNegBuilt2D[ch][KT][mb]=(TH2D*)files[ChComb][ch][KT][mb]->Get("C3QS_Negbuilt2D");
	    C4QSNegBuilt2D[ch][KT][mb]=(TH2D*)files[ChComb][ch][KT][mb]->Get("C4QS_Negbuilt2D");
	    //
	    c4QSstage2[ch][KT][mb]->SetDirectory(0);  C3QSBuilt[ch][KT][mb]->SetDirectory(0);  C4QSBuilt[ch][KT][mb]->SetDirectory(0); 
	    C3QSBuilt2D[ch][KT][mb]->SetDirectory(0);  C4QSBuilt2D[ch][KT][mb]->SetDirectory(0); 
	    C3QSNegBuilt2D[ch][KT][mb]->SetDirectory(0);  C4QSNegBuilt2D[ch][KT][mb]->SetDirectory(0);
	    //
	    r3[ch][KT][mb]=(TH1D*)files[ChComb][ch][KT][mb]->Get("r3");
	    r3[ch][KT][mb]->SetDirectory(0);
	    //
	    r42[ch][KT][mb]=(TH1D*)files[ChComb][ch][KT][mb]->Get("r42");
	    r42[ch][KT][mb]->SetDirectory(0);
	    //
	    C4QS[ChComb][ch][KT][mb]->SetBinContent(2,100); C4QS[ChComb][ch][KT][mb]->SetBinError(2,100);
	    c4QS[ChComb][ch][KT][mb]->SetBinContent(2,100); c4QS[ChComb][ch][KT][mb]->SetBinError(2,100);
	    c4QSstage1[ChComb][ch][KT][mb]->SetBinContent(2,100); c4QSstage1[ChComb][ch][KT][mb]->SetBinError(2,100);
	    c4QSstage2[ch][KT][mb]->SetBinContent(2,100); c4QSstage2[ch][KT][mb]->SetBinError(2,100);
	  }
	  files[ChComb][ch][KT][mb]->Close();
	}// KT
      }// ch
    }// ChComb
  }// mb
  
 
 
  TF1 *Unity = new TF1("Unity","1",0,100);
  Unity->SetLineStyle(2);
  Unity->SetLineColor(1);

  
  
  TH1D *C3QSmerged[2][2][10];// SC/MC, KT, Mbin
  TH1D *c3QSmerged[2][2][10];// SC/MC, KT, Mbin
  TH1D *C3QSBuiltmerged[2][10];// KT, Mbin
  TH2D *C3QSBuiltmerged2D[2][10];// KT, Mbin
  TH2D *C3QSNegBuiltmerged2D[2][10];// KT, Mbin
  //
  TH1D *C4QSmerged[3][2][10];// SC/MC, KT, Mbin
  TH1D *c4QSmerged[3][2][10];// SC/MC, KT, Mbin
  TH1D *c4QSstage1merged[3][2][10];// SC/MC, KT, Mbin
  TH1D *c4QSstage2merged[2][10];// KT, Mbin
  TH1D *C4QSBuiltmerged[2][10];// KT, Mbin
  TH2D *C4QSBuiltmerged2D[2][10];// KT, Mbin
  TH2D *C4QSNegBuiltmerged2D[2][10];// KT, Mbin
  
 for(int mb=0; mb<10; mb++){
   for(int ChComb=0; ChComb<3; ChComb++) {
     for(int KT=0; KT<2; KT++) {
       if(ChComb!=2){
	 C3QSmerged[ChComb][KT][mb] = (TH1D*)C3QS[ChComb][0][KT][mb]->Clone();
	 c3QSmerged[ChComb][KT][mb] = (TH1D*)c3QS[ChComb][0][KT][mb]->Clone();
       }
       //
       C4QSmerged[ChComb][KT][mb] = (TH1D*)C4QS[ChComb][0][KT][mb]->Clone();
       c4QSmerged[ChComb][KT][mb] = (TH1D*)c4QS[ChComb][0][KT][mb]->Clone();
       c4QSstage1merged[ChComb][KT][mb] = (TH1D*)c4QSstage1[ChComb][0][KT][mb]->Clone();
       //
       if(ChComb==0) {
	 c4QSstage2merged[KT][mb] = (TH1D*)c4QSstage2[0][KT][mb]->Clone();
	 C3QSBuiltmerged[KT][mb] = (TH1D*)C3QSBuilt[0][KT][mb]->Clone();
	 C4QSBuiltmerged[KT][mb] = (TH1D*)C4QSBuilt[0][KT][mb]->Clone();
	 C3QSBuiltmerged2D[KT][mb] = (TH2D*)C3QSBuilt2D[0][KT][mb]->Clone();
	 C4QSBuiltmerged2D[KT][mb] = (TH2D*)C4QSBuilt2D[0][KT][mb]->Clone();
	 C3QSNegBuiltmerged2D[KT][mb] = (TH2D*)C3QSNegBuilt2D[0][KT][mb]->Clone();
	 C4QSNegBuiltmerged2D[KT][mb] = (TH2D*)C4QSNegBuilt2D[0][KT][mb]->Clone();
       }
       
       if(ChComb!=2){
	 for(int bin=1; bin<=C3QSmerged[ChComb][KT][mb]->GetNbinsX(); bin++){
	   double value=0, value_e=0;
	   value = (C3QS[ChComb][0][KT][mb]->GetBinContent(bin) + C3QS[ChComb][1][KT][mb]->GetBinContent(bin)) / 2.;
	   value_e = sqrt(pow(C3QS[ChComb][0][KT][mb]->GetBinError(bin),2) + pow(C3QS[ChComb][1][KT][mb]->GetBinError(bin),2)) / sqrt(2.);
	   C3QSmerged[ChComb][KT][mb]->SetBinContent(bin, value);  C3QSmerged[ChComb][KT][mb]->SetBinError(bin, value_e);
	   //
	   value = (c3QS[ChComb][0][KT][mb]->GetBinContent(bin) + c3QS[ChComb][1][KT][mb]->GetBinContent(bin)) / 2.;
	   value_e = sqrt(pow(c3QS[ChComb][0][KT][mb]->GetBinError(bin),2) + pow(c3QS[ChComb][1][KT][mb]->GetBinError(bin),2)) / sqrt(2.);
	   c3QSmerged[ChComb][KT][mb]->SetBinContent(bin, value);  c3QSmerged[ChComb][KT][mb]->SetBinError(bin, value_e);
	   //
	   if(ChComb==0){
	     value = (C3QSBuilt[0][KT][mb]->GetBinContent(bin) + C3QSBuilt[1][KT][mb]->GetBinContent(bin)) / 2.;
	     value_e = 0;
	     C3QSBuiltmerged[KT][mb]->SetBinContent(bin, value);  C3QSBuiltmerged[KT][mb]->SetBinError(bin, value_e);
	     //
	     for(int binG=1; binG<=C3QSBuiltmerged2D[KT][mb]->GetNbinsX(); binG++){
	       value = (C3QSBuilt2D[0][KT][mb]->GetBinContent(binG, bin) + C3QSBuilt2D[1][KT][mb]->GetBinContent(binG, bin)) / 2.;
	       value_e = 0;
	       C3QSBuiltmerged2D[KT][mb]->SetBinContent(binG, bin, value);  C3QSBuiltmerged2D[KT][mb]->SetBinError(binG, bin, value_e);
	       //
	       value = (C3QSNegBuilt2D[0][KT][mb]->GetBinContent(binG, bin) + C3QSNegBuilt2D[1][KT][mb]->GetBinContent(binG, bin)) / 2.;
	       value_e = 0;
	       C3QSNegBuiltmerged2D[KT][mb]->SetBinContent(binG, bin, value);  C3QSNegBuiltmerged2D[KT][mb]->SetBinError(binG, bin, value_e);
	     }
	   }
	 }
       }
       //cout<<ChComb<<"  "<<KT<<"  "<<mb<<endl;
       //cout<<C4QS[ChComb][1][KT][mb]->GetBinContent(4)<<endl;
       //
       if(ChComb==2) continue;
       for(int bin=1; bin<=C4QSmerged[ChComb][KT][mb]->GetNbinsX(); bin++){
	 double value = (C4QS[ChComb][0][KT][mb]->GetBinContent(bin) + C4QS[ChComb][1][KT][mb]->GetBinContent(bin)) / 2.;
	 double value_e = sqrt(pow(C4QS[ChComb][0][KT][mb]->GetBinError(bin),2) + pow(C4QS[ChComb][1][KT][mb]->GetBinError(bin),2)) / sqrt(2.);
	 C4QSmerged[ChComb][KT][mb]->SetBinContent(bin, value);  C4QSmerged[ChComb][KT][mb]->SetBinError(bin, value_e);
	 //
	 value = (c4QS[ChComb][0][KT][mb]->GetBinContent(bin) + c4QS[ChComb][1][KT][mb]->GetBinContent(bin)) / 2.;
	 value_e = sqrt(pow(c4QS[ChComb][0][KT][mb]->GetBinError(bin),2) + pow(c4QS[ChComb][1][KT][mb]->GetBinError(bin),2)) / sqrt(2.);
	 c4QSmerged[ChComb][KT][mb]->SetBinContent(bin, value);  c4QSmerged[ChComb][KT][mb]->SetBinError(bin, value_e);
	 //
	 value = (c4QSstage1[ChComb][0][KT][mb]->GetBinContent(bin) + c4QSstage1[ChComb][1][KT][mb]->GetBinContent(bin)) / 2.;
	 value_e = sqrt(pow(c4QSstage1[ChComb][0][KT][mb]->GetBinError(bin),2) + pow(c4QSstage1[ChComb][1][KT][mb]->GetBinError(bin),2)) / sqrt(2.);
	 c4QSstage1merged[ChComb][KT][mb]->SetBinContent(bin, value);  c4QSstage1merged[ChComb][KT][mb]->SetBinError(bin, value_e);
	 //
	 if(ChComb==0){
	   value = (c4QSstage2[0][KT][mb]->GetBinContent(bin) + c4QSstage2[1][KT][mb]->GetBinContent(bin)) / 2.;
	   value_e = sqrt(pow(c4QSstage2[0][KT][mb]->GetBinError(bin),2) + pow(c4QSstage2[1][KT][mb]->GetBinError(bin),2)) / sqrt(2.);
	   c4QSstage2merged[KT][mb]->SetBinContent(bin, value);  c4QSstage2merged[KT][mb]->SetBinError(bin, value_e);
	   //
	   value = (C4QSBuilt[0][KT][mb]->GetBinContent(bin) + C4QSBuilt[1][KT][mb]->GetBinContent(bin)) / 2.;
	   value_e = 0;
	   C4QSBuiltmerged[KT][mb]->SetBinContent(bin, value);  C4QSBuiltmerged[KT][mb]->SetBinError(bin, value_e);
	   //
	   for(int binG=1; binG<=C4QSBuiltmerged2D[KT][mb]->GetNbinsX(); binG++){
	     value = (C4QSBuilt2D[0][KT][mb]->GetBinContent(binG, bin) + C4QSBuilt2D[1][KT][mb]->GetBinContent(binG, bin)) / 2.;
	     value_e = 0;
	     C4QSBuiltmerged2D[KT][mb]->SetBinContent(binG, bin, value);  C4QSBuiltmerged2D[KT][mb]->SetBinError(binG, bin, value_e);
	     //
	     value = (C4QSNegBuilt2D[0][KT][mb]->GetBinContent(binG, bin) + C4QSNegBuilt2D[1][KT][mb]->GetBinContent(binG, bin)) / 2.;
	     value_e = 0;
	     C4QSNegBuiltmerged2D[KT][mb]->SetBinContent(binG, bin, value);  C4QSNegBuiltmerged2D[KT][mb]->SetBinError(binG, bin, value_e);
	   }
	 }
       }
       
     }// KT
   }// ChComb
 }// mb
  
  // merge r3 histogram centralities
 TH1D *r3merged[2];// KT
 TH1D *r4merged[2];// KT
  for(int KT=0; KT<2; KT++) {
    r3merged[KT]=(TH1D*)r3[0][KT][0]->Clone();
    r4merged[KT]=(TH1D*)r42[0][KT][0]->Clone();
  }

  double mergedValue_r3[2][20]={{0}};// KT
  double mergedError_r3[2][20]={{0}};// KT
  double ErrorWeightSum_r3[2][20]={{0}};// KT
  double EnSum_r3[2][20]={{0}};// KT
  //
  double mergedValue_r4[2][20]={{0}};// KT
  double mergedError_r4[2][20]={{0}};// KT
  double ErrorWeightSum_r4[2][20]={{0}};// KT
  double EnSum_r4[2][20]={{0}};// KT
  for(int mb=0; mb<10; mb++){
    for(int ch=0; ch<2; ch++) {
      for(int KT=0; KT<2; KT++) {
	for(int bin=1; bin<=20; bin++){
	  if(r3[ch][KT][mb]->GetBinError(bin) == 0) continue;
	  mergedValue_r3[KT][bin] += r3[ch][KT][mb]->GetBinContent(bin) / pow(r3[ch][KT][mb]->GetBinError(bin),2);
	  mergedError_r3[KT][bin] += pow(r3[ch][KT][mb]->GetBinError(bin),2) / pow(r3[ch][KT][mb]->GetBinError(bin),2);
	  ErrorWeightSum_r3[KT][bin] += 1.0 / pow(r3[ch][KT][mb]->GetBinError(bin),2);
	  EnSum_r3[KT][bin]++;
	}// bin
      }// KT
    }// ch
  }// mb
  
  for(int mb=0; mb<10; mb++){
    for(int ch=0; ch<2; ch++) {
      for(int KT=0; KT<2; KT++) {
	for(int bin=1; bin<=20; bin++){
	  if(r42[ch][KT][mb]->GetBinError(bin) == 0) continue;
	  mergedValue_r4[KT][bin] += r42[ch][KT][mb]->GetBinContent(bin) / pow(r42[ch][KT][mb]->GetBinError(bin),2);
	  mergedError_r4[KT][bin] += pow(r42[ch][KT][mb]->GetBinError(bin),2) / pow(r42[ch][KT][mb]->GetBinError(bin),2);
	  ErrorWeightSum_r4[KT][bin] += 1.0 / pow(r42[ch][KT][mb]->GetBinError(bin),2);
	  EnSum_r4[KT][bin]++;
	}// bin
      }// KT
    }// ch
  }// mb


  for(int bin=1; bin<=20; bin++){
    for(int KT=0; KT<2; KT++) {
      if(ErrorWeightSum_r3[KT][bin] ==0) continue;
      if(EnSum_r3[KT][bin] == 0) continue;
      r3merged[KT]->SetBinContent(bin, mergedValue_r3[KT][bin] / ErrorWeightSum_r3[KT][bin]);
      r3merged[KT]->SetBinError(bin, sqrt(mergedError_r3[KT][bin] / ErrorWeightSum_r3[KT][bin]));
    }
  }
  
  for(int bin=1; bin<=20; bin++){
    for(int KT=0; KT<2; KT++) {
      if(ErrorWeightSum_r4[KT][bin] ==0) continue;
      if(EnSum_r4[KT][bin] == 0) continue;
      r4merged[KT]->SetBinContent(bin, mergedValue_r4[KT][bin] / ErrorWeightSum_r4[KT][bin]);
      r4merged[KT]->SetBinError(bin, sqrt(mergedError_r4[KT][bin] / ErrorWeightSum_r4[KT][bin]) / sqrt(EnSum_r4[KT][bin]));
    }
  }
  

  //////////////////////////////////////////////////////////////////////////
  // 3-pion
  if(ChProdBOI!=2){
    
    TCanvas *can1 = new TCanvas("can1", "can1",10,0,700,600);// 11,53,700,500
    can1->SetHighLightColor(2);
    gStyle->SetOptFit(0111);
    can1->SetFillColor(0);//10
    can1->SetBorderMode(0);
    can1->SetBorderSize(2);
    can1->SetFrameFillColor(0);
    can1->SetFrameBorderMode(0);
    can1->SetFrameBorderMode(0);
  
  TPad *pad1 = new TPad("pad1","pad1",0.0,0.0,1.,1.);
  //gPad->SetGridx(1);
  //gPad->SetGridy(1);
  gPad->SetTickx();
  gPad->SetTicky();
  pad1->SetTopMargin(0.0);//0.05
  pad1->SetRightMargin(0.0);//1e-2
  pad1->SetBottomMargin(0.0);//0.12
  pad1->Draw();
  pad1->cd(1);
  gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.03); gPad->SetBottomMargin(0.14);
  TLegend *legend1 = new TLegend(.42,.5, .92,.8,NULL,"brNDC");//.45 or .4 for x1
  legend1->SetBorderSize(0);
  legend1->SetFillColor(0);
  legend1->SetTextFont(TextFont);
  legend1->SetTextSize(SizeLegend);

  C3QSmerged[ChProdBOI][KTBin][MBOI]->GetXaxis()->SetTitleSize(SizeTitle);
  C3QSmerged[ChProdBOI][KTBin][MBOI]->GetXaxis()->SetLabelSize(SizeLabel);
  C3QSmerged[ChProdBOI][KTBin][MBOI]->GetYaxis()->SetTitleSize(SizeTitle);
  C3QSmerged[ChProdBOI][KTBin][MBOI]->GetYaxis()->SetLabelSize(SizeLabel);
  C3QSmerged[ChProdBOI][KTBin][MBOI]->GetXaxis()->SetTitleOffset(1.05);
  C3QSmerged[ChProdBOI][KTBin][MBOI]->GetYaxis()->SetTitleOffset(1.1);
  C3QSmerged[ChProdBOI][KTBin][MBOI]->GetXaxis()->SetNdivisions(606);
  C3QSmerged[ChProdBOI][KTBin][MBOI]->GetYaxis()->SetNdivisions(505);
  //
  //TH1D *C3QS_Syst = new TH1D("C3QS_Syst","",200,0.,0.2);
  //TH1D *C3QSBuilt_Syst = new TH1D("C3QSBuilt_Syst","",200,0.,0.2);
  TH1D *C3QS_Syst = (TH1D*)C3QSmerged[ChProdBOI][KTBin][MBOI]->Clone();
  TH1D *C3QSBuilt_Syst = (TH1D*)C3QSBuiltmerged[KTBin][MBOI]->Clone();

  for(int bin=1; bin<=C3QS_Syst->GetNbinsX(); bin++){
    double q3 = C3QS_Syst->GetXaxis()->GetBinCenter(bin);
    C3QS_Syst->SetBinContent(bin, 4.7);
    double syst1 = pow(0.001,2);// cc
    syst1 += pow(0.002 - 0.002*q3/0.1,2);// 11h to 10h
    syst1 += pow(0.9913 - 0.2231*q3 - 1,2);// f coefficients, r*<70
    syst1 += pow(0.9847 + 0.358*q3 - 2.133*q3*q3 - 1,2);// MRC
    syst1 += pow(0.975 + 0.4189*q3 - 2.055*q3*q3 - 1,2);// Muon, 92%
    syst1 += pow(0.936 + 1.194*q3 - 5.912*q3*q3 - 1,2);// fc2 scale
    syst1 += pow(0.125*exp(-61.38*q3),2);// K factorization
    syst1 = sqrt(syst1);
    syst1 *= C3QSmerged[ChProdBOI][KTBin][MBOI]->GetBinContent(bin);
    C3QS_Syst->SetBinError(bin, syst1);
    // Built
    C3QSBuilt_Syst->SetBinContent(bin, 4.7);
    double syst2 = pow(0.002 - 0.002*q3/0.1,2);// 11h to 10h
    syst2 += pow(0.9856 + 0.3285*q3 - 1.897*q3*q3 - 1,2);// MRC
    syst2 += pow(0.9786 + 0.421*q3 - 2.108*q3*q3 - 1,2);// Muon, 92%
    syst2 += pow(0.946 + 0.849*q3 - 3.316*q3*q3 - 1,2);// fc2 scale
    syst2 += pow(0.0103*exp(-41.68*q3),2);// Interpolator
    syst2 = sqrt(syst2);
    syst2 *= C3QSBuiltmerged[KTBin][MBOI]->GetBinContent(bin);
    C3QSBuilt_Syst->SetBinError(bin, syst2);
  }
  double Syst_forChi2_3[15]={0};
  for(int bin=1; bin<=15; bin++){
    double q3 = C3QSmerged[ChProdBOI][KTBin][MBOI]->GetXaxis()->GetBinCenter(bin);
    int SystBin = C3QSBuilt_Syst->GetXaxis()->FindBin(q3);
    //Syst_forChi2_3[bin-1] = fabs(C3QS_Syst->GetBinError(SystBin) - C3QSBuilt_Syst->GetBinError(SystBin));
    double SystPercent_Diff = sqrt(pow(0.125*exp(-61.38*q3*sqrt(2.)),2) + pow(0.9913 - 0.2231*q3 - 1,2) + pow(0.0103*exp(-41.68*q3),2));// K, f coefficients, Interpolator
    Syst_forChi2_3[bin-1] = SystPercent_Diff * C3QSmerged[ChProdBOI][KTBin][MBOI]->GetBinContent(bin);
    //cout<<Syst_forChi2_3[bin-1]<<endl;
  }
  C3QS_Syst->SetBinContent(1,100); C3QSBuilt_Syst->SetBinContent(1,100); 
  C3QS_Syst->SetMarkerSize(0); C3QS_Syst->SetFillColor(kBlue-10);
  C3QS_Syst->SetMarkerColor(kBlue-10);
  C3QSBuilt_Syst->SetMarkerSize(0); C3QSBuilt_Syst->SetFillColor(1); //C3QSBuilt_Syst->SetFillStyle(3004);
  C3QSBuilt_Syst->SetMarkerColor(1);
  C3QS_Syst->GetXaxis()->SetRangeUser(0.01,0.2); C3QSBuilt_Syst->GetXaxis()->SetRangeUser(0.01,0.2);
  //
  C3QSBuiltmerged[KTBin][MBOI]->GetXaxis()->SetRange(2,15);
  C3QSmerged[ChProdBOI][KTBin][MBOI]->SetMaximum(5.1);
  C3QSmerged[ChProdBOI][KTBin][MBOI]->Draw();
  C3QS_Syst->Draw("E2 same");
  C3QSBuilt_Syst->Draw("E1 same");

  C3QSmerged[ChProdBOI][KTBin][MBOI]->Draw("same");
  c3QSmerged[ChProdBOI][KTBin][MBOI]->Draw("same");
  C3QSBuiltmerged[KTBin][MBOI]->SetLineWidth(1.2);
  if(ChProdBOI==0) C3QSBuiltmerged[KTBin][MBOI]->Draw("same");
  //
  TString *proName=new TString("C3QSbuilt_G"); TString *proNameNeg=new TString("C3QSNegbuilt_G");
  TH1D *C3QSbuilt_G = (TH1D*)C3QSBuiltmerged2D[KTBin][MBOI]->ProjectionY(proName->Data(), GbinPlot, GbinPlot);
  TH1D *C3QSNegbuilt_G = (TH1D*)C3QSNegBuiltmerged2D[KTBin][MBOI]->ProjectionY(proNameNeg->Data(), GbinPlot, GbinPlot);
  proName->Append("_FullWeightDen"); proNameNeg->Append("_FullWeightDen");
  TH1D *tempDen = (TH1D*)C3QSBuiltmerged2D[KTBin][MBOI]->ProjectionY(proName->Data(), 4, 4);
  TH1D *tempDenNeg = (TH1D*)C3QSNegBuiltmerged2D[KTBin][MBOI]->ProjectionY(proNameNeg->Data(), 4, 4);
  // Add Pos with Neg weights
  tempDen->Add(tempDenNeg);
  C3QSbuilt_G->Add(C3QSNegbuilt_G);
  //  
  C3QSbuilt_G->Add(tempDen);
  C3QSbuilt_G->Divide(tempDen);
  C3QSbuilt_G->SetLineColor(2);
  C3QSbuilt_G->GetXaxis()->SetRange(2,15);
  if(ChProdBOI==0) C3QSbuilt_G->Draw("same");

  legend1->AddEntry(C3QSmerged[ChProdBOI][KTBin][MBOI],"#font[12]{C}_{3}^{QS}","p");
  legend1->AddEntry(c3QSmerged[ChProdBOI][KTBin][MBOI],"#font[12]{#bf{c}}_{3}^{QS}","p");
  if(ChProdBOI==0) legend1->AddEntry(C3QSBuiltmerged[KTBin][MBOI],"Built #font[12]{C}_{3}^{QS} (G=0%)","l");
  if(ChProdBOI==0) legend1->AddEntry(C3QSbuilt_G,"Built #font[12]{C}_{3}^{QS} (G=34%, R_{coh}=R_{ch})","l");
  legend1->Draw("same");

  Unity->Draw("same");



  

  if(ChProdBOI==0){
    TCanvas *can2 = new TCanvas("can2", "can2",800,0,700,600);// 11,53,700,500
    can2->SetHighLightColor(2);
    gStyle->SetOptFit(0111);
    can2->SetFillColor(0);//10
    can2->SetBorderMode(0);
    can2->SetBorderSize(2);
    can2->SetFrameFillColor(0);
    can2->SetFrameBorderMode(0);
    can2->SetFrameBorderMode(0);
    
    TPad *pad2 = new TPad("pad2","pad2",0.0,0.0,1.,1.);
    gPad->SetTickx();
    gPad->SetTicky();
    pad2->SetTopMargin(0.0);//0.05
    pad2->SetRightMargin(0.0);//1e-2
    pad2->SetBottomMargin(0.0);//0.12
    pad2->Draw();
    pad2->cd(1);
    gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.03); gPad->SetBottomMargin(0.14);
    TLegend *legend2 = new TLegend(.55,.75, .95,.95,NULL,"brNDC");//.45 or .4 for x1
    legend2->SetBorderSize(0);
    legend2->SetFillColor(0);
    legend2->SetTextFont(TextFont);
    legend2->SetTextSize(SizeLegend);
 

    //TH1D *chi2_PointSize_3 = new TH1D("chi2_PointSize_3","",40,-0.5,39.5);
    //TH1D *chi2_FullSize_3 = new TH1D("chi2_FullSize_3","",40,-0.5,39.5);
    TH1D *chi2_PointSize_3 = new TH1D("chi2_PointSize_3","",100,-0.5,99.5);
    TH1D *chi2_FullSize_3 = new TH1D("chi2_FullSize_3","",100,-0.5,99.5);
    chi2_PointSize_3->SetLineColor(4); chi2_FullSize_3->SetLineColor(2);
    chi2_PointSize_3->SetMarkerColor(4); chi2_FullSize_3->SetMarkerColor(2);
    chi2_PointSize_3->GetXaxis()->SetTitle("Coherent fraction (%)"); chi2_PointSize_3->GetYaxis()->SetTitle("#sqrt{#chi^{2}}");
    chi2_PointSize_3->GetXaxis()->SetTitleSize(SizeTitle);  chi2_PointSize_3->GetYaxis()->SetTitleSize(SizeTitle);
    chi2_PointSize_3->GetXaxis()->SetLabelSize(SizeLabel);  chi2_PointSize_3->GetYaxis()->SetLabelSize(SizeLabel);
    TH2D *chi2_2D_3 = new TH2D("chi2_2D_3","",5,0.5,5.5, 100,-0.5,99.5);
       
    TH1D *tempDen = (TH1D*)C3QSBuiltmerged2D[KTBin][MBOI]->ProjectionY("TPFullWeight3_Den", 4, 4);
    TH1D *tempDenNeg = (TH1D*)C3QSNegBuiltmerged2D[KTBin][MBOI]->ProjectionY("TPNegFullWeight3_Den", 4, 4);
    tempDen->Add(tempDenNeg);// Add Pos and Neg Den

    for(int binG=5; binG<=104; binG++){// 44
      TString *proName=new TString("TPFullWeight3_");
      *proName += binG;
      TH1D *tempNum = (TH1D*)C3QSBuiltmerged2D[KTBin][MBOI]->ProjectionY(proName->Data(), binG, binG);
      proName->Append("_Neg");
      TH1D *tempNumNeg = (TH1D*)C3QSNegBuiltmerged2D[KTBin][MBOI]->ProjectionY(proName->Data(), binG, binG);
      //
      // Add Pos and Neg Num
      tempNum->Add(tempNumNeg);
      //
      tempNum->Add(tempDen);
      tempNum->Divide(tempDen);
      //lowBin = C3QS->GetXaxis()->FindBin(Cutoff_FullWeight_Q3[Mbin]);
      //highBin = C3QS->GetXaxis()->FindBin(Cutoff_FullWeight_Q3[Mbin]);
      //SF=C3QS->Integral(lowBin, highBin);
      //SF /= tempNum->Integral(lowBin, highBin);
      //tempNum->Scale(SF);
      
      double value = C3QSmerged[ChProdBOI][KTBin][MBOI]->GetBinContent(Q3binChi2) - tempNum->GetBinContent(Q3binChi2);
      double err = pow(C3QSmerged[ChProdBOI][KTBin][MBOI]->GetBinError(Q3binChi2),2);// stat
      err += pow(Syst_forChi2_3[Q3binChi2-1],2);// syst
      err = sqrt(err);
      if(err <=0) continue;
      double Chi2 = pow(value / err,2);
      //
      
      //if(binG<25) {chi2_PointSize_3->SetBinContent(1 + 2*(binG-5), sqrt(Chi2)); chi2_PointSize_3->SetBinError(1 + 2*(binG-5), 0.001);}
      //else {chi2_FullSize_3->SetBinContent(1 + 2*(binG-25), sqrt(Chi2)); chi2_FullSize_3->SetBinError(1 + 2*(binG-25), 0.001);}
      if(binG<55) {chi2_PointSize_3->SetBinContent(1 + 2*(binG-5), sqrt(Chi2)); chi2_PointSize_3->SetBinError(1 + 2*(binG-5), 0.001);}
      else {chi2_FullSize_3->SetBinContent(1 + 2*(binG-55), sqrt(Chi2)); chi2_FullSize_3->SetBinError(1 + 2*(binG-55), 0.001);}
      //
      Chi2=0;
      double NDF=0;
      for(int binQ3=2; binQ3<=5; binQ3++){
	if(tempNum->GetBinContent(binQ3) <=0) continue;
	double value = C3QSmerged[ChProdBOI][KTBin][MBOI]->GetBinContent(binQ3) - tempNum->GetBinContent(binQ3);
	double err = pow(C3QSmerged[ChProdBOI][KTBin][MBOI]->GetBinError(binQ3),2);// stat
	err += pow(Syst_forChi2_3[binQ3-1],2);// syst
	err = sqrt(err);
	if(err <=0) continue;
	Chi2 = pow(value / err,2);
	//
	chi2_2D_3->SetBinContent(binQ3, binG-4, sqrt(Chi2));
      }


    }
    chi2_PointSize_3->SetMarkerStyle(20); chi2_FullSize_3->SetMarkerStyle(20);
    chi2_PointSize_3->SetMinimum(0); chi2_PointSize_3->SetMaximum(13); 
    chi2_PointSize_3->Draw();
    chi2_FullSize_3->Draw("same");
    TString *Q3binName = new TString("0.0");
    *Q3binName += Q3binChi2-1;
    Q3binName->Append(" < #font[12]{Q_{3}} < 0.0");
    *Q3binName += Q3binChi2;
    Q3binName->Append(" GeV/#font[12]{c}");
    legend2->SetHeader(Q3binName->Data());
    legend2->AddEntry(chi2_PointSize_3,"R_{coh}=0","p");
    legend2->AddEntry(chi2_FullSize_3,"R_{coh}=R_{ch}","p");
    legend2->Draw("same");

    TString *meanpTName3 = new TString("#LT #font[12]{p}_{T} #GT = 0.");
    *meanpTName3 += Q3_meanpT[KTBin][Q3binChi2-1];
    meanpTName3->Append(" GeV/#font[12]{c}");
    TLatex *Specif_pT3 = new TLatex(0.15,0.9,meanpTName3->Data());
    Specif_pT3->SetNDC();
    Specif_pT3->SetTextFont(TextFont);
    Specif_pT3->SetTextSize(SizeSpecif);
    Specif_pT3->Draw("same");

    TString *SaveNameChi2_3 = new TString("ChiSq_C3_bin");
    *SaveNameChi2_3 += Q3binChi2;
    SaveNameChi2_3->Append("_K");
    *SaveNameChi2_3 += KTBin;
    SaveNameChi2_3->Append("_M");
    *SaveNameChi2_3 += MBOI;
    SaveNameChi2_3->Append(".eps");
    //can2->SaveAs(SaveNameChi2_3->Data());




    ///////////////////////////////////////////////////////////////////////////
    // G versus Q3

    TCanvas *can2_2 = new TCanvas("can2_2", "can2_2",1300,0,700,600);// 11,53,700,500
    can2_2->SetHighLightColor(2);
    gStyle->SetOptFit(0111);
    can2_2->SetFillColor(0);//10
    can2_2->SetBorderMode(0);
    can2_2->SetBorderSize(2);
    can2_2->SetFrameFillColor(0);
    can2_2->SetFrameBorderMode(0);
    can2_2->SetFrameBorderMode(0);
    
    TPad *pad2_2 = new TPad("pad2_2","pad2_2",0.0,0.0,1.,1.);
    gPad->SetTickx();
    gPad->SetTicky();
    pad2_2->SetTopMargin(0.0);//0.05
    pad2_2->SetRightMargin(0.0);//1e-2
    pad2_2->SetBottomMargin(0.0);//0.12
    pad2_2->Draw();
    pad2_2->cd(1);
    gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.03); gPad->SetBottomMargin(0.14);
    TLegend *legend2_2 = new TLegend(.15,.15, .35,.35,NULL,"brNDC");//.45 or .4 for x1
    legend2_2->SetBorderSize(0);
    legend2_2->SetFillColor(0);
    legend2_2->SetTextFont(TextFont);
    legend2_2->SetTextSize(SizeLegend);

    TH1D *GversusQ3_Point = new TH1D("GversusQ3_Point","",5,0,0.05);
    TH1D *GversusQ3_Full = new TH1D("GversusQ3_Full","",5,0,0.05);
    for(int binQ3=2; binQ3<=5; binQ3++){
      double minG = 0;
      double minG_e1=0, minG_e2=0;
      double minChi=100;
      // Point Source
      for(int binG=1; binG<=50; binG++){// min
	if(minChi > chi2_2D_3->GetBinContent(binQ3, binG)) {
	  minChi = chi2_2D_3->GetBinContent(binQ3, binG);
	  minG = 2*(binG-1);
	}
      }
      //cout<<binQ3<<"  "<<minChi<<"  "<<minG<<endl;
      for(int binG=1; binG<=50; binG++){// error
	if(minG > 0) {
	  if(fabs(minChi - chi2_2D_3->GetBinContent(binQ3, binG)) < 1.) {
	    if(minG>2*(binG-1)) minG_e1 = fabs(minG - 2*(binG-1)); 
	    else minG_e2 = fabs(minG - 2*(binG-1));
	  }
	}else{
	  if(fabs(minChi - chi2_2D_3->GetBinContent(binQ3, binG)) < 1.) {
	    minG_e1 = fabs(minG - 2*(binG-1)); 
	  }
	}
      }
      GversusQ3_Point->SetBinContent(binQ3, minG);
      if(minG_e1 > minG_e2) GversusQ3_Point->SetBinError(binQ3, minG_e1);
      else GversusQ3_Point->SetBinError(binQ3, minG_e2);
      //
      // Full Source
      minG = 0;
      minG_e1 = 0, minG_e2=0;
      minChi=100;
      for(int binG=51; binG<=100; binG++){// min
	if(minChi > chi2_2D_3->GetBinContent(binQ3, binG)) {
	  minChi = chi2_2D_3->GetBinContent(binQ3, binG);
	  minG = 2*(binG-51);
	}
      }
      for(int binG=51; binG<=100; binG++){// error
	if(minG > 0) {
	  if(fabs(minChi - chi2_2D_3->GetBinContent(binQ3, binG)) < 1.) {
	    if(minG>2*(binG-51)) minG_e1 = fabs(minG - 2*(binG-51)); 
	    else minG_e2 = fabs(minG - 2*(binG-51));
	  }
	}else{
	  if(fabs(minChi - chi2_2D_3->GetBinContent(binQ3, binG)) < 1.) {
	    minG_e1 = fabs(minG - 2*(binG-51)); 
	  }
	}
      }
      //cout<<binQ3<<"  "<<minG<<"  "<<minG_e<<endl;
      GversusQ3_Full->SetBinContent(binQ3, minG);
      if(minG_e1 > minG_e2) GversusQ3_Full->SetBinError(binQ3, minG_e1);
      else GversusQ3_Full->SetBinError(binQ3, minG_e2);
    }
    //
    GversusQ3_Point->SetMarkerStyle(20); GversusQ3_Point->SetMarkerColor(4); GversusQ3_Point->SetLineColor(4);
    GversusQ3_Full->SetMarkerStyle(20); GversusQ3_Full->SetMarkerColor(2); GversusQ3_Full->SetLineColor(2);
    GversusQ3_Point->SetMinimum(0); GversusQ3_Point->SetMaximum(80); 
    GversusQ3_Point->GetXaxis()->SetTitle("#font[12]{Q_{3}} (GeV/#font[12]{c})"); GversusQ3_Point->GetYaxis()->SetTitle("Coherent fraction (%)");
    GversusQ3_Point->GetXaxis()->SetTitleSize(SizeTitle);  GversusQ3_Point->GetYaxis()->SetTitleSize(SizeTitle);
    GversusQ3_Point->GetXaxis()->SetLabelSize(SizeLabel);  GversusQ3_Point->GetYaxis()->SetLabelSize(SizeLabel);
    GversusQ3_Point->GetXaxis()->SetNdivisions(404); GversusQ3_Point->GetYaxis()->SetNdivisions(505);
    GversusQ3_Point->Draw();
    GversusQ3_Full->Draw("same");
    //
    legend2_2->AddEntry(GversusQ3_Point,"R_{coh}=0","p");
    legend2_2->AddEntry(GversusQ3_Full,"R_{coh}=R_{ch}","p");
    legend2_2->Draw("same");


  }
  }// ChProdBOI!=2  




  //////////////////////////////////////////////////////////////////////////
  // 4-pion
  TCanvas *can3 = new TCanvas("can3", "can3",10,700,700,600);// 11,53,700,500
  can3->SetHighLightColor(2);
  gStyle->SetOptFit(0111);
  can3->SetFillColor(0);//10
  can3->SetBorderMode(0);
  can3->SetBorderSize(2);
  can3->SetFrameFillColor(0);
  can3->SetFrameBorderMode(0);
  can3->SetFrameBorderMode(0);
  
  TPad *pad3 = new TPad("pad3","pad3",0.0,0.0,1.,1.);
  //gPad->SetGridx(1);
  //gPad->SetGridy(1);
  gPad->SetTickx();
  gPad->SetTicky();
  pad3->SetTopMargin(0.0);//0.05
  pad3->SetRightMargin(0.0);//1e-2
  pad3->SetBottomMargin(0.0);//0.12
  pad3->Draw();
  pad3->cd(1);
  gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.03); gPad->SetBottomMargin(0.14);
  
  TLegend *legend3 = new TLegend(.45,.4, .85,.8,NULL,"brNDC");//.45 or .4 for x1
  legend3->SetBorderSize(0);
  legend3->SetFillColor(0);
  legend3->SetTextFont(TextFont);
  legend3->SetTextSize(SizeLegend);

  
  C4QSmerged[ChProdBOI][KTBin][MBOI]->GetXaxis()->SetTitleSize(SizeTitle);
  C4QSmerged[ChProdBOI][KTBin][MBOI]->GetXaxis()->SetLabelSize(SizeLabel);
  C4QSmerged[ChProdBOI][KTBin][MBOI]->GetYaxis()->SetTitleSize(SizeTitle);
  C4QSmerged[ChProdBOI][KTBin][MBOI]->GetYaxis()->SetLabelSize(SizeLabel);
  C4QSmerged[ChProdBOI][KTBin][MBOI]->GetXaxis()->SetTitleOffset(1.05);
  C4QSmerged[ChProdBOI][KTBin][MBOI]->GetYaxis()->SetTitleOffset(1.1);
  C4QSmerged[ChProdBOI][KTBin][MBOI]->GetXaxis()->SetNdivisions(606);
  C4QSmerged[ChProdBOI][KTBin][MBOI]->GetYaxis()->SetNdivisions(505);
  //
  TString *proName4=new TString("C4QSbuilt_G"); TString *proNameNeg4=new TString("C4QSNegbuilt_G");
  TH1D *C4QSbuilt_G = (TH1D*)C4QSBuiltmerged2D[KTBin][MBOI]->ProjectionY(proName4->Data(), GbinPlot, GbinPlot);
  TH1D *C4QSNegbuilt_G = (TH1D*)C4QSNegBuiltmerged2D[KTBin][MBOI]->ProjectionY(proNameNeg4->Data(), GbinPlot, GbinPlot);
  proName4->Append("_FullWeightDen"); proNameNeg4->Append("_FullWeightDen");
  TH1D *tempDen4 = (TH1D*)C4QSBuiltmerged2D[KTBin][MBOI]->ProjectionY(proName4->Data(), 4, 4);
  TH1D *tempDenNeg4 = (TH1D*)C4QSNegBuiltmerged2D[KTBin][MBOI]->ProjectionY(proNameNeg4->Data(), 4, 4);
  // Add Pos with Neg weights
  tempDen4->Add(tempDenNeg4);
  C4QSbuilt_G->Add(C4QSNegbuilt_G);
  //  
  C4QSbuilt_G->Add(tempDen4);
  C4QSbuilt_G->Divide(tempDen4);
  C4QSbuilt_G->SetLineColor(2);
  //
 
  C4QSmerged[ChProdBOI][KTBin][MBOI]->SetMaximum(8.8);
  C4QSBuiltmerged[KTBin][MBOI]->GetXaxis()->SetRange(3,15);
  C4QSbuilt_G->GetXaxis()->SetRange(3,15);
  //
  //TH1D *C4QS_Syst = new TH1D("C4QS_Syst","",200,0.,0.2);
  //TH1D *C4QSBuilt_Syst = new TH1D("C4QSBuilt_Syst","",200,0.,0.2);
  TH1D *C4QS_Syst = (TH1D*)C4QSmerged[ChProdBOI][KTBin][MBOI]->Clone();
  TH1D *C4QSBuilt_Syst = (TH1D*)C4QSBuiltmerged[KTBin][MBOI]->Clone();

  for(int bin=1; bin<=C4QS_Syst->GetNbinsX(); bin++){
    double q4 = C4QS_Syst->GetXaxis()->GetBinCenter(bin);
    C4QS_Syst->SetBinContent(bin, 8.);
    double syst1 = pow(0.001,2);// cc
    syst1 += pow(0.004 - 0.004*q4/0.18,2);// 11h to 10h
    syst1 += pow(0.9975 - 0.09*q4 - 1,2);// f coefficients, r*<70
    syst1 += pow(0.9814 + 0.2471*q4 - 0.8312*q4*q4 - 1,2);// MRC
    syst1 += pow(0.9635 + 0.3475*q4 - 0.9729*q4*q4 - 1,2);// Muon, 92%
    syst1 += pow(0.900 + 1.126*q4 - 3.354*q4*q4 - 1,2);// fc2 scale
    syst1 += pow(0.125*exp(-61.38*q4),2);// K factorization
    syst1 = sqrt(syst1);
    syst1 *= C4QSmerged[ChProdBOI][KTBin][MBOI]->GetBinContent(bin);
    C4QS_Syst->SetBinError(bin, syst1);
    // Built
    C4QSBuilt_Syst->SetBinContent(bin, 8.);
    double syst2 = pow(0.004 - 0.004*q4/0.18,2);// 11h to 10h
    syst2 += pow(0.9793 + 0.2857*q4 - 0.9888*q4*q4 - 1,2);// MRC
    syst2 += pow(0.9725 + 0.2991*q4 - 0.8058*q4*q4 - 1,2);// Muon, 92%
    syst2 += pow(0.905 + 1.03*q4 - 2.977*q4*q4 - 1,2);// fc2 scale
    syst2 += pow(0.0379*exp(-42.82*q4),2);// Interpolator
    syst2 = sqrt(syst2);
    syst2 *= C4QSBuiltmerged[KTBin][MBOI]->GetBinContent(bin);
    C4QSBuilt_Syst->SetBinError(bin, syst2);
  }
  C4QS_Syst->SetBinContent(2,100); C4QSBuilt_Syst->SetBinContent(2,100); 
  C4QS_Syst->SetMarkerSize(0); C4QS_Syst->SetFillColor(kBlue-10);
  C4QS_Syst->SetMarkerColor(kBlue-10);
  C4QSBuilt_Syst->SetMarkerSize(0); C4QSBuilt_Syst->SetFillColor(1); //C4QSBuilt_Syst->SetFillStyle(3004);
  C4QSBuilt_Syst->SetMarkerColor(1);
  C4QS_Syst->GetXaxis()->SetRangeUser(0.03,0.2); C4QSBuilt_Syst->GetXaxis()->SetRangeUser(0.03,0.2);
  double Syst_forChi2_4[15]={0};
  for(int bin=1; bin<=15; bin++){
    double q4 = C4QSmerged[ChProdBOI][KTBin][MBOI]->GetXaxis()->GetBinCenter(bin);
    int SystBin = C4QSBuilt_Syst->GetXaxis()->FindBin(q4);
    //Syst_forChi2_4[bin-1] = fabs(C4QS_Syst->GetBinError(SystBin) - C4QSBuilt_Syst->GetBinError(SystBin));
    double SystPercent_Diff = sqrt(pow(0.125*exp(-61.38*q4),2) + pow(0.9975 - 0.09*q4 - 1,2) + pow(0.0379*exp(-42.82*q4),2));// K, f coefficients, Interpolator
    Syst_forChi2_4[bin-1] = SystPercent_Diff * C4QSmerged[ChProdBOI][KTBin][MBOI]->GetBinContent(bin);
  }
  //
  /*for(int bin=1; bin<=C4QSmerged[ChProdBOI][KTBin][MBOI]->GetNbinsX(); bin++){
    C4QSmerged[ChProdBOI][KTBin][MBOI]->SetBinContent(bin, fabs(C4QSmerged[ChProdBOI][KTBin][MBOI]->GetBinContent(bin) - 1));
    C4QSBuiltmerged[KTBin][MBOI]->SetBinContent(bin, fabs(C4QSBuiltmerged[KTBin][MBOI]->GetBinContent(bin) - 1));
    }*/
  C4QSmerged[ChProdBOI][KTBin][MBOI]->SetBinContent(1,100); C4QSmerged[ChProdBOI][KTBin][MBOI]->SetBinError(1,100);
  C4QSmerged[ChProdBOI][KTBin][MBOI]->Draw();
  C4QS_Syst->Draw("E2 same");
  C4QSBuilt_Syst->Draw("E1 same");
  C4QSmerged[ChProdBOI][KTBin][MBOI]->Draw("same");
  
  c4QSstage1merged[ChProdBOI][KTBin][MBOI]->Draw("same");
  if(ChProdBOI==0) c4QSstage2merged[KTBin][MBOI]->Draw("same");
  c4QSmerged[ChProdBOI][KTBin][MBOI]->Draw("same");
  C4QSBuiltmerged[KTBin][MBOI]->SetLineWidth(1.2);
  if(ChProdBOI==0) C4QSBuiltmerged[KTBin][MBOI]->Draw("same");
  
  if(ChProdBOI==0) C4QSbuilt_G->Draw("same");

  legend3->AddEntry(C4QSmerged[ChProdBOI][KTBin][MBOI],"#font[12]{C}_{4}^{QS}","p");
  legend3->AddEntry(c4QSstage1merged[ChProdBOI][KTBin][MBOI],"#font[12]{#bf{c}}_{4}^{QS} 2-pion removal","p");
  if(ChProdBOI==0) legend3->AddEntry(c4QSstage2merged[KTBin][MBOI],"#font[12]{#bf{c}}_{4}^{QS} 2-pion + 2-pair removal","p");
  legend3->AddEntry(c4QSmerged[ChProdBOI][KTBin][MBOI],"#font[12]{#bf{c}}_{4}^{QS}","p");
  if(ChProdBOI==0) legend3->AddEntry(C4QSBuiltmerged[KTBin][MBOI],"Built #font[12]{C}_{4}^{QS} (G=0%)","l");
  if(ChProdBOI==0) legend3->AddEntry(C4QSbuilt_G,"Built #font[12]{C}_{4}^{QS} (G=34%, R_{coh}=R_{ch})","l");
  legend3->Draw("same");
  
  /*TF1 *Gauss_c4Fit=new TF1("Gauss_c4Fit","[0]*(1+[1]*exp(-pow(x*[2]/0.19733,2)/3.))",0,1);
  Gauss_c4Fit->SetParameter(0,1); Gauss_c4Fit->SetParameter(1,3); Gauss_c4Fit->SetParameter(2,8); 
  Gauss_c4Fit->SetParName(0,"N");  Gauss_c4Fit->SetParName(1,"#lambda_{4}"); Gauss_c4Fit->SetParName(2,"R");
  c4QSmerged[ChProdBOI][KTBin][MBOI]->Fit(Gauss_c4Fit,"IME","",0.03,0.14);
  Gauss_c4Fit->Draw("same");*/

  // hight KT4 reference
  double y_ref[12]={0, 0, 1.00133, 0.980848, 0.988251, 0.994434, 0.999677, 1.00269, 1.00642, 1.00881, 1.01082, 1.01554};
  double y_ref_e[12]={0, 0, 0.054465, 0.00678447, 0.00194947, 0.000799564, 0.00039767, 0.000222628, 0.000135335, 8.75305e-05, 6.31392e-05, 5.53329e-05};

  TH1D *Ratio=(TH1D*)C4QSmerged[ChProdBOI][KTBin][MBOI]->Clone();
  Ratio->Divide(C4QSBuiltmerged[KTBin][MBOI]);
  Ratio->GetYaxis()->SetTitle("#font[12]{C_{4}^{QS}} / #font[12]{C_{4}^{QS}}(built)");
  Ratio->SetMinimum(0.85); Ratio->SetMaximum(1.05);
  TH1D *DoubleRatio =(TH1D*)Ratio->Clone();
  DoubleRatio->GetYaxis()->SetTitle("Low K_{T,4} ratio / High K_{T,4} ratio");

  for(int bin=1; bin<=12; bin++){
    if(y_ref[bin-1]==0) continue;
    double value = Ratio->GetBinContent(bin) / y_ref[bin-1];
    double value_e = sqrt(pow(Ratio->GetBinError(bin) / y_ref[bin-1],2) + pow(y_ref_e[bin-1]*Ratio->GetBinContent(bin) /y_ref[bin-1]/y_ref[bin-1],2));
    DoubleRatio->SetBinContent(bin, value);
    DoubleRatio->SetBinError(bin, value_e);
  }
  //Ratio->Draw();
  //DoubleRatio->Draw();

  //for(int bin=1; bin<=12; bin++) cout<<Ratio->GetBinContent(bin)<<", ";
  //cout<<endl;
  //for(int bin=1; bin<=12; bin++) cout<<Ratio->GetBinError(bin)<<", ";
  //cout<<endl;
  
  Unity->Draw("same");

  
  if(ChProdBOI==0){// chi2
    TCanvas *can4 = new TCanvas("can4", "can4",800,700,700,600);// 11,53,700,500
    can4->SetHighLightColor(2);
    gStyle->SetOptFit(0111);
    can4->SetFillColor(0);//10
    can4->SetBorderMode(0);
    can4->SetBorderSize(2);
    can4->SetFrameFillColor(0);
    can4->SetFrameBorderMode(0);
    can4->SetFrameBorderMode(0);
    
    TPad *pad4 = new TPad("pad4","pad4",0.0,0.0,1.,1.);
    gPad->SetTickx();
    gPad->SetTicky();
    pad4->SetTopMargin(0.0);//0.05
    pad4->SetRightMargin(0.0);//1e-2
    pad4->SetBottomMargin(0.0);//0.12
    pad4->Draw();
    pad4->cd(1);
    gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.03); gPad->SetBottomMargin(0.14);

    TLegend *legend4 = new TLegend(.15,.65, .4,.85,NULL,"brNDC");//.45 or .4 for x1
    legend4->SetBorderSize(0);
    legend4->SetFillColor(0);
    legend4->SetTextFont(TextFont);
    legend4->SetTextSize(SizeLegend);

    //TH1D *chi2_PointSize_4 = new TH1D("chi2_PointSize_4","",40,-0.5,39.5);
    //TH1D *chi2_FullSize_4 = new TH1D("chi2_FullSize_4","",40,-0.5,39.5);
    TH1D *chi2_PointSize_4 = new TH1D("chi2_PointSize_4","",100,-0.5,99.5);
    TH1D *chi2_FullSize_4 = new TH1D("chi2_FullSize_4","",100,-0.5,99.5);
    chi2_PointSize_4->SetLineColor(4); chi2_FullSize_4->SetLineColor(2);
    chi2_PointSize_4->SetMarkerColor(4); chi2_FullSize_4->SetMarkerColor(2);
    chi2_PointSize_4->GetXaxis()->SetTitle("Coherent fraction (%)"); chi2_PointSize_4->GetYaxis()->SetTitle("#sqrt{#chi^{2}}");
    chi2_PointSize_4->GetXaxis()->SetTitleSize(SizeTitle);  chi2_PointSize_4->GetYaxis()->SetTitleSize(SizeTitle);
    chi2_PointSize_4->GetXaxis()->SetLabelSize(SizeLabel);  chi2_PointSize_4->GetYaxis()->SetLabelSize(SizeLabel);
    chi2_PointSize_4->GetYaxis()->SetNdivisions(505);
    TH2D *chi2_2D_4 = new TH2D("chi2_2D_4","",7,0.5,7.5, 100,-0.5,99.5);
       

    TH1D *tempDen = (TH1D*)C4QSBuiltmerged2D[KTBin][MBOI]->ProjectionY("TPFullWeight4_Den", 4, 4);
    TH1D *tempDenNeg = (TH1D*)C4QSNegBuiltmerged2D[KTBin][MBOI]->ProjectionY("TPNegFullWeight4_Den", 4, 4);
    tempDen->Add(tempDenNeg);// Add Pos and Neg Weight

    for(int binG=5; binG<=104; binG++){// 44
      TString *proName=new TString("TPFullWeight4_");
      *proName += binG;
      TH1D *tempNum = (TH1D*)C4QSBuiltmerged2D[KTBin][MBOI]->ProjectionY(proName->Data(), binG, binG);
      proName->Append("_Neg");
      TH1D *tempNumNeg = (TH1D*)C4QSNegBuiltmerged2D[KTBin][MBOI]->ProjectionY(proName->Data(), binG, binG);
      //
      // Add Pos and Neg Weights
      tempNum->Add(tempNumNeg);
      //
      tempNum->Add(tempDen);
      tempNum->Divide(tempDen);
      //lowBin = C4QS->GetXaxis()->FindBin(Cutoff_FullWeight_Q4[Mbin]);
      //highBin = C4QS->GetXaxis()->FindBin(Cutoff_FullWeight_Q4[Mbin]);
      //SF=C4QS->Integral(lowBin, highBin);
      //SF /= tempNum->Integral(lowBin, highBin);
      //tempNum->Scale(SF);
      

      if(tempNum->GetBinContent(Q4binChi2) <=0) continue;
      double value = C4QSmerged[ChProdBOI][KTBin][MBOI]->GetBinContent(Q4binChi2) - tempNum->GetBinContent(Q4binChi2);
      double err = pow(C4QSmerged[ChProdBOI][KTBin][MBOI]->GetBinError(Q4binChi2),2);
      err += pow(Syst_forChi2_4[Q4binChi2-1],2);
      err = sqrt(err);
      if(err<=0) continue;
      double Chi2 = pow(value / err,2);
      //
      
      //if(binG<25) {chi2_PointSize_4->SetBinContent(1 + 2*(binG-5), sqrt(fabs(Chi2))); chi2_PointSize_4->SetBinError(1 + 2*(binG-5), 0.001);}
      //else {chi2_FullSize_4->SetBinContent(1 + 2*(binG-25), sqrt(fabs(Chi2))); chi2_FullSize_4->SetBinError(1 + 2*(binG-25), 0.001);}
      if(binG<55) {chi2_PointSize_4->SetBinContent(1 + 2*(binG-5), sqrt(fabs(Chi2))); chi2_PointSize_4->SetBinError(1 + 2*(binG-5), 0.001);}
      else {chi2_FullSize_4->SetBinContent(1 + 2*(binG-55), sqrt(fabs(Chi2))); chi2_FullSize_4->SetBinError(1 + 2*(binG-55), 0.001);}
      //
      
      for(int binQ4=3; binQ4<=7; binQ4++){
	if(tempNum->GetBinContent(binQ4) <=0) continue;
	double value = C4QSmerged[ChProdBOI][KTBin][MBOI]->GetBinContent(binQ4) - tempNum->GetBinContent(binQ4);
	double err = pow(C4QSmerged[ChProdBOI][KTBin][MBOI]->GetBinError(binQ4),2);
	err += pow(Syst_forChi2_4[binQ4-1],2);
	err = sqrt(err);
	if(err<=0) continue;
	Chi2 = pow(value / err,2);
	//
	chi2_2D_4->SetBinContent(binQ4, binG-4, sqrt(fabs(Chi2)));
      }
      
    }
    chi2_PointSize_4->SetMarkerStyle(20); chi2_FullSize_4->SetMarkerStyle(20);
    chi2_PointSize_4->SetMinimum(0); chi2_PointSize_4->SetMaximum(13);
    
    TString *Q4binName = new TString("0.0");
    
    if(int((Q4binChi2-1)*1.5*10)%10 == 0) *Q4binName += int((Q4binChi2-1)*1.5);
    else {*Q4binName += int((Q4binChi2-1)*1.5); *Q4binName += 5;}
    Q4binName->Append(" < #font[12]{Q_{4}} < 0.0");
    if(int((Q4binChi2)*1.5*10)%10 == 0) *Q4binName += int((Q4binChi2)*1.5);
    else {*Q4binName += int((Q4binChi2)*1.5); *Q4binName += 5;}
    Q4binName->Append(" GeV/#font[12]{c}");
    legend4->SetHeader(Q4binName->Data());
    chi2_PointSize_4->Draw();
    chi2_FullSize_4->Draw("same");
    legend4->AddEntry(chi2_PointSize_4,"R_{coh}=0","p");
    legend4->AddEntry(chi2_FullSize_4,"R_{coh}=R_{ch}","p");
    legend4->Draw("same");

    TString *meanpTName = new TString("#LT #font[12]{p}_{T} #GT = 0.");
    *meanpTName += Q4_meanpT[KTBin][Q4binChi2-1];
    meanpTName->Append(" GeV/#font[12]{c}");
    TLatex *Specif_pT = new TLatex(0.15,0.9,meanpTName->Data());
    Specif_pT->SetNDC();
    Specif_pT->SetTextFont(TextFont);
    Specif_pT->SetTextSize(SizeSpecif);
    Specif_pT->Draw("same");
   

    TString *SaveNameChi2_4 = new TString("ChiSq_C4_bin");
    *SaveNameChi2_4 += Q4binChi2;
    SaveNameChi2_4->Append("_K");
    *SaveNameChi2_4 += KTBin;
    SaveNameChi2_4->Append("_M");
    *SaveNameChi2_4 += MBOI;
    SaveNameChi2_4->Append(".eps");
    //can4->SaveAs(SaveNameChi2_4->Data());
    



    ///////////////////////////////////////////////////////////////////////////
    // G versus Q4

    TCanvas *can5 = new TCanvas("can5", "can5",1300,700,700,600);// 11,53,700,500
    can5->SetHighLightColor(2);
    gStyle->SetOptFit(0111);
    can5->SetFillColor(0);//10
    can5->SetBorderMode(0);
    can5->SetBorderSize(2);
    can5->SetFrameFillColor(0);
    can5->SetFrameBorderMode(0);
    can5->SetFrameBorderMode(0);
    
    TPad *pad5 = new TPad("pad5","pad5",0.0,0.0,1.,1.);
    gPad->SetTickx();
    gPad->SetTicky();
    pad5->SetTopMargin(0.0);//0.05
    pad5->SetRightMargin(0.0);//1e-2
    pad5->SetBottomMargin(0.0);//0.12
    pad5->Draw();
    pad5->cd(1);
    gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.03); gPad->SetBottomMargin(0.14);

    TLegend *legend5 = new TLegend(.15,.75, .3,.95,NULL,"brNDC");//.45 or .4 for x1
    legend5->SetBorderSize(0);
    legend5->SetFillColor(0);
    legend5->SetTextFont(TextFont);
    legend5->SetTextSize(SizeLegend);

    TH1D *GversusQ4_Point = new TH1D("GversusQ4_Point","",7,0,0.105);
    TH1D *GversusQ4_Full = new TH1D("GversusQ4_Full","",7,0,0.105);
    for(int binQ4=3; binQ4<=7; binQ4++){
      double minG = 0;
      double minG_e1 = 0, minG_e2=0;
      double minChi=100;
      // Point Source
      for(int binG=1; binG<=50; binG++){// min
	if(minChi > chi2_2D_4->GetBinContent(binQ4, binG)) {
	  minChi = chi2_2D_4->GetBinContent(binQ4, binG);
	  minG = 2*(binG-1);
	}
      }
      //cout<<binQ4<<"  "<<minChi<<endl;
      for(int binG=1; binG<=50; binG++){// error
	if(minG > 0) {
	  if(fabs(minChi - chi2_2D_4->GetBinContent(binQ4, binG)) < 1.) {
	    if(minG>2*(binG-1)) minG_e1 = fabs(minG - 2*(binG-1)); 
	    else minG_e2 = fabs(minG - 2*(binG-1));
	  }
	}else{
	  if(fabs(minChi - chi2_2D_4->GetBinContent(binQ4, binG)) < 1.) {
	    minG_e1 = fabs(minG - 2*(binG-1)); 
	  }
	}
      }
      GversusQ4_Point->SetBinContent(binQ4, minG);
      if(minG_e1>minG_e2) GversusQ4_Point->SetBinError(binQ4, minG_e1);
      else GversusQ4_Point->SetBinError(binQ4, minG_e2);
      //
      // Full Source
      minG = 0;
      minG_e1 = 0, minG_e2=0;
      minChi=100;
      for(int binG=51; binG<=100; binG++){// min
	if(minChi > chi2_2D_4->GetBinContent(binQ4, binG)) {
	  minChi = chi2_2D_4->GetBinContent(binQ4, binG);
	  minG = 2*(binG-51);
	}
      }
      for(int binG=51; binG<=100; binG++){// error
	if(minG > 0) {
	  if(fabs(minChi - chi2_2D_4->GetBinContent(binQ4, binG)) < 1.) {
	    if(minG>2*(binG-51)) minG_e1 = fabs(minG - 2*(binG-51)); 
	    else minG_e2 = fabs(minG - 2*(binG-51));
	  }
	}else{
	  if(fabs(minChi - chi2_2D_4->GetBinContent(binQ4, binG)) < 1.) {
	    minG_e1 = fabs(minG - 2*(binG-51)); 
	  }
	}
      }
      //cout<<binQ4<<"  "<<minG<<"  "<<minG_e<<endl;
      GversusQ4_Full->SetBinContent(binQ4, minG);
      if(minG_e1>minG_e2) GversusQ4_Full->SetBinError(binQ4, minG_e1);
      else GversusQ4_Full->SetBinError(binQ4, minG_e2);
    }
    //
    GversusQ4_Point->SetMarkerStyle(20); GversusQ4_Point->SetMarkerColor(4); GversusQ4_Point->SetLineColor(4);
    GversusQ4_Full->SetMarkerStyle(20); GversusQ4_Full->SetMarkerColor(2); GversusQ4_Full->SetLineColor(2);
    GversusQ4_Point->SetMinimum(0); GversusQ4_Point->SetMaximum(40); 
    GversusQ4_Point->GetXaxis()->SetTitle("#font[12]{Q_{4}} (GeV/#font[12]{c})"); GversusQ4_Point->GetYaxis()->SetTitle("Coherent fraction (%)");
    GversusQ4_Point->GetXaxis()->SetTitleSize(SizeTitle);  GversusQ4_Point->GetYaxis()->SetTitleSize(SizeTitle);
    GversusQ4_Point->GetXaxis()->SetLabelSize(SizeLabel);  GversusQ4_Point->GetYaxis()->SetLabelSize(SizeLabel);
    GversusQ4_Point->GetYaxis()->SetNdivisions(505);
    GversusQ4_Point->Draw();
    GversusQ4_Full->Draw("same");
    //
    legend5->AddEntry(GversusQ4_Point,"R_{coh}=0","p");
    legend5->AddEntry(GversusQ4_Full,"R_{coh}=R_{ch}","p");
    legend5->Draw("same");

  }

  


  //////////////////////////////////////////////////////////////////////////
  // r3
  /* TCanvas *can1 = new TCanvas("can1", "can1",10,0,600,600);// 11,53,700,500
  can1->SetHighLightColor(2);
  gStyle->SetOptFit(0111);
  can1->SetFillColor(0);//10
  can1->SetBorderMode(0);
  can1->SetBorderSize(2);
  can1->SetFrameFillColor(0);
  can1->SetFrameBorderMode(0);
  can1->SetFrameBorderMode(0);
  
  TPad *pad1 = new TPad("pad1","pad1",0.0,0.0,1.,1.);
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  gPad->SetTickx();
  gPad->SetTicky();
  pad1->SetTopMargin(0.02);//0.05
  pad1->SetRightMargin(0.01);//1e-2
  pad1->SetBottomMargin(0.07);//0.12
  pad1->Draw();
  pad1->cd(1);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.03);
  
  r3merged[0]->SetMinimum(1.45); r3merged[0]->SetMaximum(2.45); 
  r3merged[1]->SetMarkerColor(2); r3merged[1]->SetLineColor(2); 
  //r3merged[0]->Draw();
  //r3merged[1]->Draw("same");
  //
  r4merged[0]->SetMinimum(1.45); r4merged[0]->SetMaximum(10.45); 
  r4merged[1]->SetMarkerColor(2); r4merged[1]->SetLineColor(2);
  r4merged[0]->Draw();
  r4merged[1]->Draw("same");
  */

}

