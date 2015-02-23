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

const int CollisionType=2;// Pb-Pb(0), p-Pb(1), pp(2)
const int ChProdBOI=0;// 0=SameCharge, 1=MixedCharge1, 2=MixedCharge2
const int EDBin=0;// KT3,4 bin. 0=low KT bin.  1=high KT bin
const int MBOI=0;// Centrality bin: 0-9
const int Q3binChi2= 5;// 2-5
const int Q4binChi2= 3;// 3-7
const int MBins=1;
const int EDBins=2;
const int RcohIndex=6;// 0(point source), 6(Full Size)
const int GValue = 0;
const int Gbinstart= 5 + RcohIndex*25;
const int GbinPlot=int( (GValue) /2. ) + 25*RcohIndex + 5;// +5 (Rcoh=0), +30 (Rcoh=1),....
bool FitBuild=1;
bool ReNormBuiltBaseline=1;
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
  
  int Ysize=600;
  if(ChProdBOI==0) Ysize=800;
  double Ystart=0.3;
  if(ChProdBOI!=0) Ystart=0.005;
  double Xmarge=0.01;
  if(ChProdBOI!=0) Xmarge=0.14;
  
  TString *System=new TString("");
  if(CollisionType==0) System->Append("Pb-Pb #\sqrt{#font[12]{s}_{NN}}=2.76 TeV");
  else if(CollisionType==1) System->Append("p-Pb #\sqrt{#font[12]{s}_{NN}}=5.02 TeV");
  else System->Append("pp #\sqrt{#font[12]{s}}=7 TeV");
  TLatex *ALICEspecif = new TLatex(0.62,.8,System->Data());// ALICE specifications
  ALICEspecif->SetNDC(1);
  ALICEspecif->SetTextFont(TextFont);
  ALICEspecif->SetTextSize(SizeSpecif);
  //
  TString *KT=new TString("");
  if(EDBin==0) KT->Append("0.16<#font[12]{K}_{T4}<0.3 GeV/#font[12]{c}");
  else KT->Append("0.3<#font[12]{K}_{T4}<1.0 GeV/#font[12]{c}");
  TLatex *KTspecif = new TLatex(0.64,0.93,KT->Data());
  KTspecif->SetNDC(1);
  KTspecif->SetTextFont(TextFont);
  KTspecif->SetTextSize(SizeSpecif);
    
  TString *Centname=new TString("0-5%");
  TLatex *Centspecif = new TLatex(0.52,0.8,Centname->Data());
  Centspecif->SetNDC(1);
  Centspecif->SetTextFont(TextFont);
  Centspecif->SetTextSize(SizeSpecif);

  TString *BuiltNameC2 = new TString("Built from #font[12]{C}_{2} (G=");
  TString *BuiltNameC3 = new TString("Built from #font[12]{C}_{3} (G=");
  TString *BuiltNamec3 = new TString("Built from #font[12]{#bf{c}}_{3} (G=");
  *BuiltNameC2 += GValue;
  *BuiltNameC3 += GValue;
  *BuiltNamec3 += GValue;
  BuiltNameC2->Append("%)");
  BuiltNameC3->Append("%)");
  BuiltNamec3->Append("%)");
  
  int Q3_meanpT[2][5]={{0,213,228,233,237},{0,316,336,354,366}};// low then high KT3
  int Q4_meanpT[2][7]={{0,0,221,229,234,238,241},{0,0,325,335,344,351,355}};// low then high KT4

  double ReNormL_3;
  double ReNormH_3;
  double ReNormL_4;
  double ReNormH_4;
  if(CollisionType==0){
    ReNormL_3=0.085;
    ReNormH_3=0.095;
    ReNormL_4=0.14;
    ReNormH_4=0.16;
  }else{
    ReNormL_4=0.46;
    ReNormH_4=0.49;
  }


  TFile *files[3][2][2][10][3];// SC/MC1/MC2, +/-, KT, MBINS, CT
  //
 

  TH1D *C3QS[2][2][EDBins][10][3];// SC/MC, +/-, KT, Mbin, CT
  TH1D *c3QS[2][2][EDBins][10][3];// SC/MC, +/-, KT, Mbin, CT
  TH1D *C3QSBuilt[2][EDBins][10][3];// +/-, KT, Mbin, CT
  TH1D *c3QSBuilt[2][EDBins][10][3];// +/-, KT, Mbin, CT
  TH2D *C3QSBuilt2D[2][EDBins][10][3];// +/-, KT, Mbin, CT
  TH2D *C3QSNegBuilt2D[2][EDBins][10][3];// +/-, KT, Mbin, CT
  //TH2D *c3QSBuilt2D[2][EDBins][10][3];// +/-, KT, Mbin, CT
  //TH2D *c3QSNegBuilt2D[2][EDBins][10][3];// +/-, KT, Mbin, CT
  //
  TH1D *C4QS[3][2][EDBins][10][3];// SC/MC, +/-, KT, Mbin, CT
  TH1D *c4QS[3][2][EDBins][10][3];// SC/MC, +/-, KT, Mbin, CT
  TH1D *a4QS[3][2][EDBins][10][3];// SC/MC, +/-, KT, Mbin, CT
  TH1D *b4QS[2][EDBins][10][3];// +/-, KT, Mbin, CT
  TH1D *C4QSBuilt[2][EDBins][10][3];// +/-, KT, Mbin, CT
  TH1D *a4QSBuilt[2][EDBins][10][3];// +/-, KT, Mbin, CT
  TH1D *b4QSBuilt[2][EDBins][10][3];// +/-, KT, Mbin, CT
  TH1D *c4QSBuilt[2][EDBins][10][3];// +/-, KT, Mbin, CT
  TH2D *C4QSBuilt2D[2][EDBins][10][3];// +/-, KT, Mbin, CT
  TH2D *C4QSNegBuilt2D[2][EDBins][10][3];// +/-, KT, Mbin, CT
  TH2D *a4QSBuilt2D[2][EDBins][10][3];// +/-, KT, Mbin, CT
  TH2D *a4QSNegBuilt2D[2][EDBins][10][3];// +/-, KT, Mbin, CT
  TH2D *b4QSBuilt2D[2][EDBins][10][3];// +/-, KT, Mbin, CT
  TH2D *b4QSNegBuilt2D[2][EDBins][10][3];// +/-, KT, Mbin, CT
  TH2D *c4QSBuilt2D[2][EDBins][10][3];// +/-, KT, Mbin, CT
  TH2D *c4QSNegBuilt2D[2][EDBins][10][3];// +/-, KT, Mbin, CT
  //
  TH3D *C4QSBuiltFromFits3D[2][EDBins][10][3];// +/-, ED, Mbin, CT
  TH3D *a4QSBuiltFromFits3D[2][EDBins][10][3];// +/-, ED, Mbin, CT
  TH3D *b4QSBuiltFromFits3D[2][EDBins][10][3];// +/-, ED, Mbin, CT
  TH3D *c4QSBuiltFromFits3D[2][EDBins][10][3];// +/-, ED, Mbin, CT
  //
  TH1D *r3[2][EDBins][10];// +/-, KT, Mbin
  TH1D *r42[2][EDBins][10];// +/-, KT, Mbin
  
  //
  //////////////////////////////

  // Start File access
  for(int CT=0; CT<3; CT++){
    for(int mb=0; mb<MBins; mb++){
      if(CT>0 && mb>0) continue;
      for(int ChComb=0; ChComb<3; ChComb++) {// SC or MC1 or MC2
	for(int ch=0; ch<2; ch++) {// - or +
	  for(int ED=0; ED<EDBins; ED++) {// KT3 or q2 bin
	    if(ChComb==2 && ch!=0) continue;
	    
	    TString *name = new TString("OutFiles/OutFile");
	    //
	    name->Append("_CT");
	    *name += CT;
	    if(ChComb==0) name->Append("_SC");
	    else if(ChComb==1) name->Append("_MC1");
	    else name->Append("_MC2");
	    //
	    if(ch==0) name->Append("_Neg_");
	    else name->Append("_Pos_");
	    name->Append("ED");
	    *name += ED+1;
	    name->Append("_M");
	    if(mb<10) {*name += mb;}
	    else {*name += 0;}
	    name->Append(".root");
	    files[ChComb][ch][ED][mb][CT] = new TFile(name->Data(),"READ");
	    ///////////////////////////////
	    if(ChComb!=2){
	      C3QS[ChComb][ch][ED][mb][CT]=(TH1D*)files[ChComb][ch][ED][mb][CT]->Get("C3QS");
	      c3QS[ChComb][ch][ED][mb][CT]=(TH1D*)files[ChComb][ch][ED][mb][CT]->Get("c3QS");
	      C3QS[ChComb][ch][ED][mb][CT]->SetDirectory(0); c3QS[ChComb][ch][ED][mb][CT]->SetDirectory(0); 
	      C3QS[ChComb][ch][ED][mb][CT]->GetXaxis()->SetLabelFont(TextFont);  C3QS[ChComb][ch][ED][mb][CT]->GetYaxis()->SetLabelFont(TextFont);
	      c3QS[ChComb][ch][ED][mb][CT]->GetXaxis()->SetLabelFont(TextFont);  c3QS[ChComb][ch][ED][mb][CT]->GetYaxis()->SetLabelFont(TextFont);
	      C3QS[ChComb][ch][ED][mb][CT]->GetXaxis()->SetTitleFont(TextFont);  C3QS[ChComb][ch][ED][mb][CT]->GetYaxis()->SetTitleFont(TextFont);
	      c3QS[ChComb][ch][ED][mb][CT]->GetXaxis()->SetTitleFont(TextFont);  c3QS[ChComb][ch][ED][mb][CT]->GetYaxis()->SetTitleFont(TextFont);
	      C3QS[ChComb][ch][ED][mb][CT]->GetXaxis()->SetTitle("#font[12]{Q}_{3} (GeV/#font[12]{c})");
	      c3QS[ChComb][ch][ED][mb][CT]->GetXaxis()->SetTitle("#font[12]{Q}_{3} (GeV/#font[12]{c})");
	      C3QS[ChComb][ch][ED][mb][CT]->GetYaxis()->SetTitle("Three Pion");
	      c3QS[ChComb][ch][ED][mb][CT]->GetYaxis()->SetTitle("#font[12]{#bf{c}}_{3}");
	    }
	    C4QS[ChComb][ch][ED][mb][CT]=(TH1D*)files[ChComb][ch][ED][mb][CT]->Get("C4QS");
	    c4QS[ChComb][ch][ED][mb][CT]=(TH1D*)files[ChComb][ch][ED][mb][CT]->Get("c4QS");
	    a4QS[ChComb][ch][ED][mb][CT]=(TH1D*)files[ChComb][ch][ED][mb][CT]->Get("a4QS");
	    C4QS[ChComb][ch][ED][mb][CT]->SetDirectory(0); c4QS[ChComb][ch][ED][mb][CT]->SetDirectory(0); 
	    a4QS[ChComb][ch][ED][mb][CT]->SetDirectory(0);
	    C4QS[ChComb][ch][ED][mb][CT]->GetXaxis()->SetTitle("#font[12]{Q}_{4} (GeV/#font[12]{c})");
	    c4QS[ChComb][ch][ED][mb][CT]->GetXaxis()->SetTitle("#font[12]{Q}_{4} (GeV/#font[12]{c})");
	    C4QS[ChComb][ch][ED][mb][CT]->GetYaxis()->SetTitle("Four Pion");
	    c4QS[ChComb][ch][ED][mb][CT]->GetYaxis()->SetTitle("#font[12]{#bf{c}}_{4}");
	    
	    if(ChComb==0){
	      b4QS[ch][ED][mb][CT]=(TH1D*)files[ChComb][ch][ED][mb][CT]->Get("b4QS");
	      b4QS[ch][ED][mb][CT]->SetDirectory(0);
	      if(CT==0){
		C3QSBuilt[ch][ED][mb][CT]=(TH1D*)files[ChComb][ch][ED][mb][CT]->Get("C3QS_built");
		c3QSBuilt[ch][ED][mb][CT]=(TH1D*)files[ChComb][ch][ED][mb][CT]->Get("c3QS_built");
		C4QSBuilt[ch][ED][mb][CT]=(TH1D*)files[ChComb][ch][ED][mb][CT]->Get("C4QS_built");
		a4QSBuilt[ch][ED][mb][CT]=(TH1D*)files[ChComb][ch][ED][mb][CT]->Get("a4QS_built");
		b4QSBuilt[ch][ED][mb][CT]=(TH1D*)files[ChComb][ch][ED][mb][CT]->Get("b4QS_built");
		c4QSBuilt[ch][ED][mb][CT]=(TH1D*)files[ChComb][ch][ED][mb][CT]->Get("c4QS_built");

		C3QSBuilt2D[ch][ED][mb][CT]=(TH2D*)files[ChComb][ch][ED][mb][CT]->Get("C3QS_built2D");
		//c3QSBuilt2D[ch][ED][mb][CT]=(TH2D*)files[ChComb][ch][ED][mb][CT]->Get("c3QS_built2D");
		C4QSBuilt2D[ch][ED][mb][CT]=(TH2D*)files[ChComb][ch][ED][mb][CT]->Get("C4QS_built2D");
		C3QSNegBuilt2D[ch][ED][mb][CT]=(TH2D*)files[ChComb][ch][ED][mb][CT]->Get("C3QS_Negbuilt2D");
		//c3QSNegBuilt2D[ch][ED][mb][CT]=(TH2D*)files[ChComb][ch][ED][mb][CT]->Get("c3QS_Negbuilt2D");
		C4QSNegBuilt2D[ch][ED][mb][CT]=(TH2D*)files[ChComb][ch][ED][mb][CT]->Get("C4QS_Negbuilt2D");
		a4QSBuilt2D[ch][ED][mb][CT]=(TH2D*)files[ChComb][ch][ED][mb][CT]->Get("a4QS_built2D");
		a4QSNegBuilt2D[ch][ED][mb][CT]=(TH2D*)files[ChComb][ch][ED][mb][CT]->Get("a4QS_Negbuilt2D");
		b4QSBuilt2D[ch][ED][mb][CT]=(TH2D*)files[ChComb][ch][ED][mb][CT]->Get("b4QS_built2D");
		b4QSNegBuilt2D[ch][ED][mb][CT]=(TH2D*)files[ChComb][ch][ED][mb][CT]->Get("b4QS_Negbuilt2D");
		c4QSBuilt2D[ch][ED][mb][CT]=(TH2D*)files[ChComb][ch][ED][mb][CT]->Get("c4QS_built2D");
		c4QSNegBuilt2D[ch][ED][mb][CT]=(TH2D*)files[ChComb][ch][ED][mb][CT]->Get("c4QS_Negbuilt2D");
		//
		C3QSBuilt[ch][ED][mb][CT]->SetDirectory(0); c3QSBuilt[ch][ED][mb][CT]->SetDirectory(0);
		C3QSBuilt2D[ch][ED][mb][CT]->SetDirectory(0); C3QSNegBuilt2D[ch][ED][mb][CT]->SetDirectory(0);
		C4QSBuilt[ch][ED][mb][CT]->SetDirectory(0); C4QSBuilt2D[ch][ED][mb][CT]->SetDirectory(0); 
		C4QSNegBuilt2D[ch][ED][mb][CT]->SetDirectory(0);
		a4QSBuilt2D[ch][ED][mb][CT]->SetDirectory(0); a4QSNegBuilt2D[ch][ED][mb][CT]->SetDirectory(0);
		b4QSBuilt2D[ch][ED][mb][CT]->SetDirectory(0); b4QSNegBuilt2D[ch][ED][mb][CT]->SetDirectory(0);
		c4QSBuilt2D[ch][ED][mb][CT]->SetDirectory(0); c4QSNegBuilt2D[ch][ED][mb][CT]->SetDirectory(0);
		a4QSBuilt[ch][ED][mb][CT]->SetDirectory(0); b4QSBuilt[ch][ED][mb][CT]->SetDirectory(0); c4QSBuilt[ch][ED][mb][CT]->SetDirectory(0);
		//
		r3[ch][ED][mb]=(TH1D*)files[ChComb][ch][ED][mb][CT]->Get("r3");
		r3[ch][ED][mb]->SetDirectory(0);
		//
		r42[ch][ED][mb]=(TH1D*)files[ChComb][ch][ED][mb][CT]->Get("r42");
		r42[ch][ED][mb]->SetDirectory(0);
	      }
	      C4QSBuiltFromFits3D[ch][ED][mb][CT]=(TH3D*)files[ChComb][ch][ED][mb][CT]->Get("C4QS_BuiltFromFits3D");
	      a4QSBuiltFromFits3D[ch][ED][mb][CT]=(TH3D*)files[ChComb][ch][ED][mb][CT]->Get("a4QS_BuiltFromFits3D");
	      b4QSBuiltFromFits3D[ch][ED][mb][CT]=(TH3D*)files[ChComb][ch][ED][mb][CT]->Get("b4QS_BuiltFromFits3D");
	      c4QSBuiltFromFits3D[ch][ED][mb][CT]=(TH3D*)files[ChComb][ch][ED][mb][CT]->Get("c4QS_BuiltFromFits3D");
	      //
	      C4QSBuiltFromFits3D[ch][ED][mb][CT]->SetDirectory(0);
	      a4QSBuiltFromFits3D[ch][ED][mb][CT]->SetDirectory(0);
	      b4QSBuiltFromFits3D[ch][ED][mb][CT]->SetDirectory(0);
	      c4QSBuiltFromFits3D[ch][ED][mb][CT]->SetDirectory(0);
	      //
	      C4QS[ChComb][ch][ED][mb][CT]->SetBinContent(2,100); C4QS[ChComb][ch][ED][mb][CT]->SetBinError(2,100);
	      c4QS[ChComb][ch][ED][mb][CT]->SetBinContent(2,100); c4QS[ChComb][ch][ED][mb][CT]->SetBinError(2,100);
	      a4QS[ChComb][ch][ED][mb][CT]->SetBinContent(2,100); a4QS[ChComb][ch][ED][mb][CT]->SetBinError(2,100);
	      b4QS[ch][ED][mb][CT]->SetBinContent(2,100); b4QS[ch][ED][mb][CT]->SetBinError(2,100);
	      
	    }
	    files[ChComb][ch][ED][mb][CT]->Close();
	  }// ED
	}// ch
      }// ChComb
    }// mb
  }// CT
    
  
  TF1 *Unity = new TF1("Unity","1",0,100);
  Unity->SetLineStyle(2);
  Unity->SetLineColor(1);

  
  
  TH1D *C3QSmerged[2][EDBins][10][3];// SC/MC, ED, Mbin, CT
  TH1D *c3QSmerged[2][EDBins][10][3];// SC/MC, ED, Mbin, CT
  TH1D *C3QSBuiltmerged[EDBins][10][3];// ED, Mbin, CT
  TH1D *c3QSBuiltmerged[EDBins][10][3];// ED, Mbin, CT
  TH2D *C3QSBuiltmerged2D[EDBins][10][3];// ED, Mbin, CT
  TH2D *C3QSNegBuiltmerged2D[EDBins][10][3];// ED, Mbin, CT
  //TH2D *c3QSBuiltmerged2D[EDBins][10][3];// ED, Mbin, CT
  //TH2D *c3QSNegBuiltmerged2D[EDBins][10][3];// ED, Mbin, CT
  //
  TH1D *C4QSmerged[3][EDBins][10][3];// SC/MC, ED, Mbin, CT
  TH1D *c4QSmerged[3][EDBins][10][3];// SC/MC, ED, Mbin, CT
  TH1D *a4QSmerged[3][EDBins][10][3];// SC/MC, ED, Mbin, CT
  TH1D *b4QSmerged[EDBins][10][3];// ED, Mbin, CT
  TH1D *C4QSBuiltmerged[EDBins][10][3];// ED, Mbin, CT
  TH1D *a4QSBuiltmerged[EDBins][10][3];// ED, Mbin, CT
  TH1D *b4QSBuiltmerged[EDBins][10][3];// ED, Mbin, CT
  TH1D *c4QSBuiltmerged[EDBins][10][3];// ED, Mbin, CT
  TH2D *C4QSBuiltmerged2D[EDBins][10][3];// ED, Mbin, CT
  TH2D *C4QSNegBuiltmerged2D[EDBins][10][3];// ED, Mbin, CT
  TH2D *a4QSBuiltmerged2D[EDBins][10][3];// ED, Mbin, CT
  TH2D *a4QSNegBuiltmerged2D[EDBins][10][3];// ED, Mbin, CT
  TH2D *b4QSBuiltmerged2D[EDBins][10][3];// ED, Mbin, CT
  TH2D *b4QSNegBuiltmerged2D[EDBins][10][3];// ED, Mbin, CT
  TH2D *c4QSBuiltmerged2D[EDBins][10][3];// ED, Mbin, CT
  TH2D *c4QSNegBuiltmerged2D[EDBins][10][3];// ED, Mbin, CT
  //
  TH3D *C4QSBuiltFromFitsmerged3D[EDBins][10][3];// ED, Mbin, CT
  TH3D *a4QSBuiltFromFitsmerged3D[EDBins][10][3];// ED, Mbin, CT
  TH3D *b4QSBuiltFromFitsmerged3D[EDBins][10][3];// ED, Mbin, CT
  TH3D *c4QSBuiltFromFitsmerged3D[EDBins][10][3];// ED, Mbin, CT
  //
  TH1D *C4QSBuiltFromFitsmerged1[EDBins][10][3];// ED, Mbin, CT
  TH1D *C4QSBuiltFromFitsmerged2[EDBins][10][3];// ED, Mbin, CT
  TH1D *a4QSBuiltFromFitsmerged1[EDBins][10][3];// ED, Mbin, CT
  TH1D *a4QSBuiltFromFitsmerged2[EDBins][10][3];// ED, Mbin, CT
  TH1D *b4QSBuiltFromFitsmerged1[EDBins][10][3];// ED, Mbin, CT
  TH1D *b4QSBuiltFromFitsmerged2[EDBins][10][3];// ED, Mbin, CT
  TH1D *c4QSBuiltFromFitsmerged1[EDBins][10][3];// ED, Mbin, CT
  TH1D *c4QSBuiltFromFitsmerged2[EDBins][10][3];// ED, Mbin, CT
  
  

  for(int CT=0; CT<3; CT++){
    for(int mb=0; mb<MBins; mb++){
      for(int ChComb=0; ChComb<3; ChComb++) {
	for(int ED=0; ED<EDBins; ED++) {
	  if(ChComb!=2){
	    C3QSmerged[ChComb][ED][mb][CT] = (TH1D*)C3QS[ChComb][0][ED][mb][CT]->Clone();
	    c3QSmerged[ChComb][ED][mb][CT] = (TH1D*)c3QS[ChComb][0][ED][mb][CT]->Clone();
	  }
	  //
	  C4QSmerged[ChComb][ED][mb][CT] = (TH1D*)C4QS[ChComb][0][ED][mb][CT]->Clone();
	  c4QSmerged[ChComb][ED][mb][CT] = (TH1D*)c4QS[ChComb][0][ED][mb][CT]->Clone();
	  a4QSmerged[ChComb][ED][mb][CT] = (TH1D*)a4QS[ChComb][0][ED][mb][CT]->Clone();
	  //
	  if(ChComb==0) {
	    b4QSmerged[ED][mb][CT] = (TH1D*)b4QS[0][ED][mb][CT]->Clone();
	    if(CT==0){
	      C3QSBuiltmerged[ED][mb][CT] = (TH1D*)C3QSBuilt[0][ED][mb][CT]->Clone();
	      c3QSBuiltmerged[ED][mb][CT] = (TH1D*)c3QSBuilt[0][ED][mb][CT]->Clone();
	      C4QSBuiltmerged[ED][mb][CT] = (TH1D*)C4QSBuilt[0][ED][mb][CT]->Clone();
	      a4QSBuiltmerged[ED][mb][CT] = (TH1D*)a4QSBuilt[0][ED][mb][CT]->Clone();
	      b4QSBuiltmerged[ED][mb][CT] = (TH1D*)b4QSBuilt[0][ED][mb][CT]->Clone();
	      c4QSBuiltmerged[ED][mb][CT] = (TH1D*)c4QSBuilt[0][ED][mb][CT]->Clone();
	      
	      C3QSBuiltmerged2D[ED][mb][CT] = (TH2D*)C3QSBuilt2D[0][ED][mb][CT]->Clone();
	      C4QSBuiltmerged2D[ED][mb][CT] = (TH2D*)C4QSBuilt2D[0][ED][mb][CT]->Clone();
	      C3QSNegBuiltmerged2D[ED][mb][CT] = (TH2D*)C3QSNegBuilt2D[0][ED][mb][CT]->Clone();
	      C4QSNegBuiltmerged2D[ED][mb][CT] = (TH2D*)C4QSNegBuilt2D[0][ED][mb][CT]->Clone();
	      a4QSBuiltmerged2D[ED][mb][CT] = (TH2D*)a4QSBuilt2D[0][ED][mb][CT]->Clone();
	      b4QSBuiltmerged2D[ED][mb][CT] = (TH2D*)b4QSBuilt2D[0][ED][mb][CT]->Clone();
	      c4QSBuiltmerged2D[ED][mb][CT] = (TH2D*)c4QSBuilt2D[0][ED][mb][CT]->Clone();
	      a4QSNegBuiltmerged2D[ED][mb][CT] = (TH2D*)a4QSNegBuilt2D[0][ED][mb][CT]->Clone();
	      b4QSNegBuiltmerged2D[ED][mb][CT] = (TH2D*)b4QSNegBuilt2D[0][ED][mb][CT]->Clone();
	      c4QSNegBuiltmerged2D[ED][mb][CT] = (TH2D*)c4QSNegBuilt2D[0][ED][mb][CT]->Clone();
	      //
	      c3QSBuiltmerged[ED][mb][CT]->Add(c3QSBuilt[1][ED][mb][CT]); c3QSBuiltmerged[ED][mb][CT]->Scale(0.5);
	      //
	      C4QSBuiltmerged2D[ED][mb][CT]->Add(C4QSBuilt2D[1][ED][mb][CT]);
	      C4QSNegBuiltmerged2D[ED][mb][CT]->Add(C4QSNegBuilt2D[1][ED][mb][CT]);
	      a4QSBuiltmerged2D[ED][mb][CT]->Add(a4QSBuilt2D[1][ED][mb][CT]);
	      a4QSNegBuiltmerged2D[ED][mb][CT]->Add(a4QSNegBuilt2D[1][ED][mb][CT]);
	      b4QSBuiltmerged2D[ED][mb][CT]->Add(b4QSBuilt2D[1][ED][mb][CT]);
	      b4QSNegBuiltmerged2D[ED][mb][CT]->Add(b4QSNegBuilt2D[1][ED][mb][CT]);
	      c4QSBuiltmerged2D[ED][mb][CT]->Add(c4QSBuilt2D[1][ED][mb][CT]);
	      c4QSNegBuiltmerged2D[ED][mb][CT]->Add(c4QSNegBuilt2D[1][ED][mb][CT]);
	    }
	    C4QSBuiltFromFitsmerged3D[ED][mb][CT] = (TH3D*)C4QSBuiltFromFits3D[0][ED][mb][CT]->Clone();
	    a4QSBuiltFromFitsmerged3D[ED][mb][CT] = (TH3D*)a4QSBuiltFromFits3D[0][ED][mb][CT]->Clone();
	    b4QSBuiltFromFitsmerged3D[ED][mb][CT] = (TH3D*)b4QSBuiltFromFits3D[0][ED][mb][CT]->Clone();
	    c4QSBuiltFromFitsmerged3D[ED][mb][CT] = (TH3D*)c4QSBuiltFromFits3D[0][ED][mb][CT]->Clone();
	    C4QSBuiltFromFitsmerged3D[ED][mb][CT]->Add(C4QSBuiltFromFits3D[1][ED][mb][CT]);
	    a4QSBuiltFromFitsmerged3D[ED][mb][CT]->Add(a4QSBuiltFromFits3D[1][ED][mb][CT]);
	    b4QSBuiltFromFitsmerged3D[ED][mb][CT]->Add(b4QSBuiltFromFits3D[1][ED][mb][CT]);
	    c4QSBuiltFromFitsmerged3D[ED][mb][CT]->Add(c4QSBuiltFromFits3D[1][ED][mb][CT]);
	    //
	    TString *DenName=new TString("DenNamePro_"); *DenName += ED; *DenName += mb; *DenName += CT; 
	    TH1D *tempDen2 = (TH1D*)C4QSBuiltFromFitsmerged3D[ED][mb][CT]->ProjectionZ(DenName->Data(), 1,1, 4, 4);
	    TString *FitBuildName=new TString("FitBuildPro_");
	    *FitBuildName += ED; *FitBuildName += mb; *FitBuildName += CT;
	    C4QSBuiltFromFitsmerged1[ED][mb][CT] = (TH1D*)C4QSBuiltFromFitsmerged3D[ED][mb][CT]->ProjectionZ(FitBuildName->Data(),1,1,GbinPlot,GbinPlot);
	    FitBuildName->Append("1");
	    C4QSBuiltFromFitsmerged2[ED][mb][CT] = (TH1D*)C4QSBuiltFromFitsmerged3D[ED][mb][CT]->ProjectionZ(FitBuildName->Data(),2,2,GbinPlot,GbinPlot);
	    FitBuildName->Append("1");
	    a4QSBuiltFromFitsmerged1[ED][mb][CT] = (TH1D*)a4QSBuiltFromFitsmerged3D[ED][mb][CT]->ProjectionZ(FitBuildName->Data(),1,1,GbinPlot,GbinPlot);
	    FitBuildName->Append("1");
	    a4QSBuiltFromFitsmerged2[ED][mb][CT] = (TH1D*)a4QSBuiltFromFitsmerged3D[ED][mb][CT]->ProjectionZ(FitBuildName->Data(),2,2,GbinPlot,GbinPlot);
	    FitBuildName->Append("1");
	    b4QSBuiltFromFitsmerged1[ED][mb][CT] = (TH1D*)b4QSBuiltFromFitsmerged3D[ED][mb][CT]->ProjectionZ(FitBuildName->Data(),1,1,GbinPlot,GbinPlot);
	    FitBuildName->Append("1");
	    b4QSBuiltFromFitsmerged2[ED][mb][CT] = (TH1D*)b4QSBuiltFromFitsmerged3D[ED][mb][CT]->ProjectionZ(FitBuildName->Data(),2,2,GbinPlot,GbinPlot);
	    FitBuildName->Append("1");
	    c4QSBuiltFromFitsmerged1[ED][mb][CT] = (TH1D*)c4QSBuiltFromFitsmerged3D[ED][mb][CT]->ProjectionZ(FitBuildName->Data(),1,1,GbinPlot,GbinPlot);
	    FitBuildName->Append("1");
	    c4QSBuiltFromFitsmerged2[ED][mb][CT] = (TH1D*)c4QSBuiltFromFitsmerged3D[ED][mb][CT]->ProjectionZ(FitBuildName->Data(),2,2,GbinPlot,GbinPlot);
	    //
	    C4QSBuiltFromFitsmerged1[ED][mb][CT]->Add(tempDen2);
	    C4QSBuiltFromFitsmerged1[ED][mb][CT]->Divide(tempDen2);
	    C4QSBuiltFromFitsmerged1[ED][mb][CT]->SetLineColor(4); C4QSBuiltFromFitsmerged1[ED][mb][CT]->SetLineStyle(2);
	    //
	    a4QSBuiltFromFitsmerged1[ED][mb][CT]->Add(tempDen2);
	    a4QSBuiltFromFitsmerged1[ED][mb][CT]->Divide(tempDen2);
	    a4QSBuiltFromFitsmerged1[ED][mb][CT]->SetLineColor(2); a4QSBuiltFromFitsmerged1[ED][mb][CT]->SetLineStyle(2);
	    //
	    b4QSBuiltFromFitsmerged1[ED][mb][CT]->Add(tempDen2);
	    b4QSBuiltFromFitsmerged1[ED][mb][CT]->Divide(tempDen2);
	    b4QSBuiltFromFitsmerged1[ED][mb][CT]->SetLineColor(6); b4QSBuiltFromFitsmerged1[ED][mb][CT]->SetLineStyle(2);
	    //
	    c4QSBuiltFromFitsmerged1[ED][mb][CT]->Add(tempDen2);
	    c4QSBuiltFromFitsmerged1[ED][mb][CT]->Divide(tempDen2);
	    c4QSBuiltFromFitsmerged1[ED][mb][CT]->SetLineColor(1); c4QSBuiltFromFitsmerged1[ED][mb][CT]->SetLineStyle(2);
	    ////
	    C4QSBuiltFromFitsmerged2[ED][mb][CT]->Add(tempDen2);
	    C4QSBuiltFromFitsmerged2[ED][mb][CT]->Divide(tempDen2);
	    C4QSBuiltFromFitsmerged2[ED][mb][CT]->SetLineColor(4);
	    //
	    a4QSBuiltFromFitsmerged2[ED][mb][CT]->Add(tempDen2);
	    a4QSBuiltFromFitsmerged2[ED][mb][CT]->Divide(tempDen2);
	    a4QSBuiltFromFitsmerged2[ED][mb][CT]->SetLineColor(2);
	    //
	    b4QSBuiltFromFitsmerged2[ED][mb][CT]->Add(tempDen2);
	    b4QSBuiltFromFitsmerged2[ED][mb][CT]->Divide(tempDen2);
	    b4QSBuiltFromFitsmerged2[ED][mb][CT]->SetLineColor(6);
	    //
	    c4QSBuiltFromFitsmerged2[ED][mb][CT]->Add(tempDen2);
	    c4QSBuiltFromFitsmerged2[ED][mb][CT]->Divide(tempDen2);
	    c4QSBuiltFromFitsmerged2[ED][mb][CT]->SetLineColor(1);
	  }
	  
	  if(ChComb!=2){
	    for(int bin=1; bin<=C3QSmerged[ChComb][ED][mb][CT]->GetNbinsX(); bin++){
	      double value=0, value_e=0;
	      value = (C3QS[ChComb][0][ED][mb][CT]->GetBinContent(bin) + C3QS[ChComb][1][ED][mb][CT]->GetBinContent(bin)) / 2.;
	      value_e = sqrt(pow(C3QS[ChComb][0][ED][mb][CT]->GetBinError(bin),2) + pow(C3QS[ChComb][1][ED][mb][CT]->GetBinError(bin),2)) / sqrt(2.);
	      C3QSmerged[ChComb][ED][mb][CT]->SetBinContent(bin, value);  C3QSmerged[ChComb][ED][mb][CT]->SetBinError(bin, value_e);
	      //
	      value = (c3QS[ChComb][0][ED][mb][CT]->GetBinContent(bin) + c3QS[ChComb][1][ED][mb][CT]->GetBinContent(bin)) / 2.;
	      value_e = sqrt(pow(c3QS[ChComb][0][ED][mb][CT]->GetBinError(bin),2) + pow(c3QS[ChComb][1][ED][mb][CT]->GetBinError(bin),2)) / sqrt(2.);
	      c3QSmerged[ChComb][ED][mb][CT]->SetBinContent(bin, value);  c3QSmerged[ChComb][ED][mb][CT]->SetBinError(bin, value_e);
	      //
	      if(ChComb==0 && CT==0){
		value = (C3QSBuilt[0][ED][mb][CT]->GetBinContent(bin) + C3QSBuilt[1][ED][mb][CT]->GetBinContent(bin)) / 2.;
		value_e = 0;
		C3QSBuiltmerged[ED][mb][CT]->SetBinContent(bin, value);  C3QSBuiltmerged[ED][mb][CT]->SetBinError(bin, value_e);
		//
		value = (c3QSBuilt[0][ED][mb][CT]->GetBinContent(bin) + c3QSBuilt[1][ED][mb][CT]->GetBinContent(bin)) / 2.;
		value_e = 0;
		c3QSBuiltmerged[ED][mb][CT]->SetBinContent(bin, value);  c3QSBuiltmerged[ED][mb][CT]->SetBinError(bin, value_e);
		//
		for(int binG=1; binG<=C3QSBuiltmerged2D[ED][mb][CT]->GetNbinsX(); binG++){
		  value = (C3QSBuilt2D[0][ED][mb][CT]->GetBinContent(binG, bin) + C3QSBuilt2D[1][ED][mb][CT]->GetBinContent(binG, bin)) / 2.;
		  value_e = 0;
		  C3QSBuiltmerged2D[ED][mb][CT]->SetBinContent(binG, bin, value);  C3QSBuiltmerged2D[ED][mb][CT]->SetBinError(binG, bin, value_e);
		  //
		  value = (C3QSNegBuilt2D[0][ED][mb][CT]->GetBinContent(binG, bin) + C3QSNegBuilt2D[1][ED][mb][CT]->GetBinContent(binG, bin)) / 2.;
		  value_e = 0;
		  C3QSNegBuiltmerged2D[ED][mb][CT]->SetBinContent(binG, bin, value);  C3QSNegBuiltmerged2D[ED][mb][CT]->SetBinError(binG, bin, value_e);
		}
	      }
	    }
	  }
	  
	  //
	  

	  if(ChComb!=2){
	    for(int bin=1; bin<=C4QSmerged[ChComb][ED][mb][CT]->GetNbinsX(); bin++){
	      double value = (C4QS[ChComb][0][ED][mb][CT]->GetBinContent(bin) + C4QS[ChComb][1][ED][mb][CT]->GetBinContent(bin)) / 2.;
	      double value_e = sqrt(pow(C4QS[ChComb][0][ED][mb][CT]->GetBinError(bin),2) + pow(C4QS[ChComb][1][ED][mb][CT]->GetBinError(bin),2)) / sqrt(2.);
	      C4QSmerged[ChComb][ED][mb][CT]->SetBinContent(bin, value);  C4QSmerged[ChComb][ED][mb][CT]->SetBinError(bin, value_e);
	      //
	      value = (c4QS[ChComb][0][ED][mb][CT]->GetBinContent(bin) + c4QS[ChComb][1][ED][mb][CT]->GetBinContent(bin)) / 2.;
	      value_e = sqrt(pow(c4QS[ChComb][0][ED][mb][CT]->GetBinError(bin),2) + pow(c4QS[ChComb][1][ED][mb][CT]->GetBinError(bin),2)) / sqrt(2.);
	      c4QSmerged[ChComb][ED][mb][CT]->SetBinContent(bin, value);  c4QSmerged[ChComb][ED][mb][CT]->SetBinError(bin, value_e);
	      //
	      value = (a4QS[ChComb][0][ED][mb][CT]->GetBinContent(bin) + a4QS[ChComb][1][ED][mb][CT]->GetBinContent(bin)) / 2.;
	      value_e = sqrt(pow(a4QS[ChComb][0][ED][mb][CT]->GetBinError(bin),2) + pow(a4QS[ChComb][1][ED][mb][CT]->GetBinError(bin),2)) / sqrt(2.);
	      a4QSmerged[ChComb][ED][mb][CT]->SetBinContent(bin, value);  a4QSmerged[ChComb][ED][mb][CT]->SetBinError(bin, value_e);
	      //
	      if(ChComb==0){
		value = (b4QS[0][ED][mb][CT]->GetBinContent(bin) + b4QS[1][ED][mb][CT]->GetBinContent(bin)) / 2.;
		value_e = sqrt(pow(b4QS[0][ED][mb][CT]->GetBinError(bin),2) + pow(b4QS[1][ED][mb][CT]->GetBinError(bin),2)) / sqrt(2.);
		b4QSmerged[ED][mb][CT]->SetBinContent(bin, value);  b4QSmerged[ED][mb][CT]->SetBinError(bin, value_e);
		//
		if(CT==0){
		  value = (C4QSBuilt[0][ED][mb][CT]->GetBinContent(bin) + C4QSBuilt[1][ED][mb][CT]->GetBinContent(bin)) / 2.;
		  if(bin<3) value=0;
		  value_e = 0;
		  C4QSBuiltmerged[ED][mb][CT]->SetBinContent(bin, value);  C4QSBuiltmerged[ED][mb][CT]->SetBinError(bin, value_e);
		  //
		  value = (a4QSBuilt[0][ED][mb][CT]->GetBinContent(bin) + a4QSBuilt[1][ED][mb][CT]->GetBinContent(bin)) / 2.;
		  if(bin<3) value=0;
		  value_e = 0;
		  a4QSBuiltmerged[ED][mb][CT]->SetBinContent(bin, value);  a4QSBuiltmerged[ED][mb][CT]->SetBinError(bin, value_e);
		  //
		  value = (b4QSBuilt[0][ED][mb][CT]->GetBinContent(bin) + b4QSBuilt[1][ED][mb][CT]->GetBinContent(bin)) / 2.;
		  if(bin<3) value=0;
		  value_e = 0;
		  b4QSBuiltmerged[ED][mb][CT]->SetBinContent(bin, value);  b4QSBuiltmerged[ED][mb][CT]->SetBinError(bin, value_e);
		  //
		  value = (c4QSBuilt[0][ED][mb][CT]->GetBinContent(bin) + c4QSBuilt[1][ED][mb][CT]->GetBinContent(bin)) / 2.;
		  if(bin<3) value=0;
		  value_e = 0;
		  c4QSBuiltmerged[ED][mb][CT]->SetBinContent(bin, value);  c4QSBuiltmerged[ED][mb][CT]->SetBinError(bin, value_e);
		}
	      }
	    }
	  }
	  for(int bin=1; bin<=20; bin++){
	    if(C4QSmerged[ChComb][ED][mb][CT]->GetBinError(bin) > 0.5*C4QSmerged[ChComb][ED][mb][CT]->GetBinContent(bin)) {
	      C4QSmerged[ChComb][ED][mb][CT]->SetBinContent(bin,100.); C4QSmerged[ChComb][ED][mb][CT]->SetBinError(bin,100.); 
	    }
	    if(a4QSmerged[ChComb][ED][mb][CT]->GetBinError(bin) > 0.5*a4QSmerged[ChComb][ED][mb][CT]->GetBinContent(bin)) {
	      a4QSmerged[ChComb][ED][mb][CT]->SetBinContent(bin,100.); a4QSmerged[ChComb][ED][mb][CT]->SetBinError(bin,100.); 
	    }
	    if(ChComb==0){
	      if(b4QSmerged[ED][mb][CT]->GetBinError(bin) > 0.5*b4QSmerged[ED][mb][CT]->GetBinContent(bin)) {
		b4QSmerged[ED][mb][CT]->SetBinContent(bin,100.); b4QSmerged[ED][mb][CT]->SetBinError(bin,100.); 
	      }
	    }
	    if(c4QSmerged[ChComb][ED][mb][CT]->GetBinError(bin) > 0.3*c4QSmerged[ChComb][ED][mb][CT]->GetBinContent(bin)) {
	      c4QSmerged[ChComb][ED][mb][CT]->SetBinContent(bin,100.); c4QSmerged[ChComb][ED][mb][CT]->SetBinError(bin,100.); 
	    }
	  }
	  
	  
	}// ED
      }// ChComb
    }// mb
  }// CT
  
  // merge r3 histogram centralities
  TH1D *r3merged[2];// ED
  TH1D *r4merged[2];// ED
  for(int ED=0; ED<EDBins; ED++) {
    r3merged[ED]=(TH1D*)r3[0][ED][0]->Clone();
    r4merged[ED]=(TH1D*)r42[0][ED][0]->Clone();
  }
  
  double mergedValue_r3[EDBins][20]={{0}};// ED
  double mergedError_r3[EDBins][20]={{0}};// ED
  double ErrorWeightSum_r3[EDBins][20]={{0}};// ED
  double EnSum_r3[EDBins][20]={{0}};// ED
  //
  double mergedValue_r4[EDBins][20]={{0}};// ED
  double mergedError_r4[EDBins][20]={{0}};// ED
  double ErrorWeightSum_r4[EDBins][20]={{0}};// ED
  double EnSum_r4[EDBins][20]={{0}};// ED
  for(int mb=0; mb<MBins; mb++){
    for(int ch=0; ch<2; ch++) {
      for(int ED=0; ED<EDBins; ED++) {
	for(int bin=1; bin<=20; bin++){
	  if(r3[ch][ED][mb]->GetBinError(bin) == 0) continue;
	  mergedValue_r3[ED][bin] += r3[ch][ED][mb]->GetBinContent(bin) / pow(r3[ch][ED][mb]->GetBinError(bin),2);
	  mergedError_r3[ED][bin] += pow(r3[ch][ED][mb]->GetBinError(bin),2) / pow(r3[ch][ED][mb]->GetBinError(bin),2);
	  ErrorWeightSum_r3[ED][bin] += 1.0 / pow(r3[ch][ED][mb]->GetBinError(bin),2);
	  EnSum_r3[ED][bin]++;
	}// bin
      }// ED
    }// ch
  }// mb
  
  for(int mb=0; mb<MBins; mb++){
    for(int ch=0; ch<2; ch++) {
      for(int ED=0; ED<EDBins; ED++) {
	for(int bin=1; bin<=20; bin++){
	  if(r42[ch][ED][mb]->GetBinError(bin) == 0) continue;
	  mergedValue_r4[ED][bin] += r42[ch][ED][mb]->GetBinContent(bin) / pow(r42[ch][ED][mb]->GetBinError(bin),2);
	  mergedError_r4[ED][bin] += pow(r42[ch][ED][mb]->GetBinError(bin),2) / pow(r42[ch][ED][mb]->GetBinError(bin),2);
	  ErrorWeightSum_r4[ED][bin] += 1.0 / pow(r42[ch][ED][mb]->GetBinError(bin),2);
	  EnSum_r4[ED][bin]++;
	}// bin
      }// ED
    }// ch
  }// mb


  for(int bin=1; bin<=20; bin++){
    for(int ED=0; ED<EDBins; ED++) {
      if(ErrorWeightSum_r3[ED][bin] ==0) continue;
      if(EnSum_r3[ED][bin] == 0) continue;
      r3merged[ED]->SetBinContent(bin, mergedValue_r3[ED][bin] / ErrorWeightSum_r3[ED][bin]);
      r3merged[ED]->SetBinError(bin, sqrt(mergedError_r3[ED][bin] / ErrorWeightSum_r3[ED][bin]));
    }
  }
  
  for(int bin=1; bin<=20; bin++){
    for(int ED=0; ED<EDBins; ED++) {
      if(ErrorWeightSum_r4[ED][bin] ==0) continue;
      if(EnSum_r4[ED][bin] == 0) continue;
      r4merged[ED]->SetBinContent(bin, mergedValue_r4[ED][bin] / ErrorWeightSum_r4[ED][bin]);
      r4merged[ED]->SetBinError(bin, sqrt(mergedError_r4[ED][bin] / ErrorWeightSum_r4[ED][bin]) / sqrt(EnSum_r4[ED][bin]));
    }
  }
  

  gStyle->SetErrorX(0.0001);
  //gStyle->SetErrorX(0.5);

  //////////////////////////////////////////////////////////////////////////
  // 3-pion
  if(ChProdBOI!=2){
    TCanvas *can1 = new TCanvas("can1", "can1",10,0,700,Ysize);// 11,53,700,500
    can1->SetHighLightColor(2);
    gStyle->SetOptFit(0111);
    can1->SetFillColor(0);//10
    can1->SetBorderMode(0);
    can1->SetBorderSize(2);
    can1->SetFrameFillColor(0);
    can1->SetFrameBorderMode(0);
    can1->SetFrameBorderMode(0);
    
    TPad *pad1 = new TPad("pad1","pad1",0.,Ystart,1.,1.);
    TPad *pad1_2 = new TPad("pad1_2","pad1_2",0.,0.,1.,Ystart);
    pad1->Draw();
    if(CollisionType==0 && ChProdBOI==0) pad1_2->Draw();
    pad1->cd();
    gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.03); 
    gPad->SetBottomMargin(Xmarge);
    

    TLegend *legend1 = new TLegend(.39,.4, .89,.75,NULL,"brNDC");//.45 or .4 for x1
    legend1->SetBorderSize(0);
    legend1->SetFillColor(0);
    legend1->SetTextFont(TextFont);
    legend1->SetTextSize(SizeLegend);
    TLegend *legend1_2=(TLegend*)legend1->Clone();
    legend1_2->SetX1(0.8); legend1_2->SetX2(0.9); legend1_2->SetY1(0.4); legend1_2->SetY2(0.65); 
    legend1_2->SetTextSize(2.5*SizeLegend);
    
    C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->SetTitleSize(SizeTitle);
    C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->SetLabelSize(SizeLabel);
    C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetYaxis()->SetTitleSize(SizeTitle);
    C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetYaxis()->SetLabelSize(SizeLabel);
    C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->SetTitleOffset(1.05);
    C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetYaxis()->SetTitleOffset(1.1);
    C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->SetNdivisions(606);
    C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetYaxis()->SetNdivisions(505);
    if(CollisionType!=0) C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->SetRangeUser(0,0.35);
    //
    
       
    //TH1D *C3QS_Syst = (TH1D*)C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Clone();
    TH1D *C3QS_Syst = new TH1D("C3QS_Syst","",C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetNbinsX(),C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->GetBinLowEdge(1)+0.0015, C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->GetBinUpEdge(C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetNbinsX())+0.0015); 
    
    if(CollisionType==0) TH1D *C3QSBuilt_Syst = (TH1D*)C3QSBuiltmerged[EDBin][MBOI][CollisionType]->Clone();
    

    for(int bin=1; bin<=C3QS_Syst->GetNbinsX(); bin++){
      double q3 = C3QS_Syst->GetXaxis()->GetBinCenter(bin);
      C3QS_Syst->SetBinContent(bin, 4.7);
      double syst1 = pow(0.001,2);// cc
      syst1 += pow(0.002 - 0.002*q3/0.1,2);// 11h to 10h
      syst1 += pow(0.01*(1 - q3/0.1),2);// f coefficients, r*<70. was pow(0.9913 - 0.2231*q3 - 1,2)
      syst1 += pow(0.9847 + 0.358*q3 - 2.133*q3*q3 - 1,2);// MRC
      syst1 += pow(0.975 + 0.4189*q3 - 2.055*q3*q3 - 1,2);// Muon, 92%
      syst1 += pow(0.936 + 1.194*q3 - 5.912*q3*q3 - 1,2);// fc2 scale
      //syst1 += pow(0.125*exp(-61.38*q3),2);// K factorization
      //syst1 += pow(0.0858 - 2.808*q3 + 19.72*q3*q3,2);// MC cumulant residue
      syst1 = sqrt(syst1);
      syst1 *= C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetBinContent(bin);
      C3QS_Syst->SetBinError(bin, syst1);
      // Built
      if(CollisionType==0){
	C3QSBuilt_Syst->SetBinContent(bin, 4.7);
	double syst2 = pow(0.002 - 0.002*q3/0.1,2);// 11h to 10h
	syst2 += pow(0.9856 + 0.3285*q3 - 1.897*q3*q3 - 1,2);// MRC
	syst2 += pow(0.9786 + 0.421*q3 - 2.108*q3*q3 - 1,2);// Muon, 92%
	syst2 += pow(0.946 + 0.849*q3 - 3.316*q3*q3 - 1,2);// fc2 scale
	syst2 += pow((0.024-0.276*q3)/2.,2);// Interpolator. was pow(0.0103*exp(-41.68*q3),2)
	if(EDBin==1) syst2 += pow( 0.9973 + 0.462*q3 - 14.2*pow(q3,2) - 2.135*pow(q3,3) + 805.3*pow(q3,4) + 1358*pow(q3,5) - 1,2);// Tij sign
	syst2 = sqrt(syst2);
	syst2 *= C3QSBuiltmerged[EDBin][MBOI][CollisionType]->GetBinContent(bin);
	C3QSBuilt_Syst->SetBinError(bin, syst2);
      }
    }
    double Syst_forChi2_3[15]={0};
    double Syst_forRatio[15]={0};
    for(int bin=1; bin<=15; bin++){
      double q3 = C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->GetBinCenter(bin);
      int SystBin = C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->FindBin(q3);
      //double SystPercent_Diff = pow(0.125*exp(-61.38*q3*sqrt(2.)),2) + pow(0.01*(1 - q3/0.1),2);// K, f coefficients (old way)
      double SystPercent_Diff = pow(0.01*(1 - q3/0.1),2);// f coefficients
      SystPercent_Diff += pow((0.024-0.276*q3)/2.,2);// Interpolator
      SystPercent_Diff += pow(0.9973 + 0.462*q3 - 14.2*pow(q3,2) - 2.135*pow(q3,3) + 805.3*pow(q3,4) + 1358*pow(q3,5) - 1,2);// Tij sign
      SystPercent_Diff += pow( (0.9847 + 0.358*q3 - 2.133*q3*q3) - (0.9856 + 0.3285*q3 - 1.897*q3*q3),2);// MRC
      SystPercent_Diff += pow( (0.975 + 0.4189*q3 - 2.055*q3*q3) - (0.9786 + 0.421*q3 - 2.108*q3*q3),2);// Muon
      SystPercent_Diff += pow( (0.936 + 1.194*q3 - 5.912*q3*q3) - (0.946 + 0.849*q3 - 3.316*q3*q3),2);// fc2 scale
      //SystPercent_Diff += pow( 0.0858 - 2.808*q3 + 19.72*q3*q3,2);// MC cumulant residue
      SystPercent_Diff = sqrt(SystPercent_Diff);
      Syst_forChi2_3[bin-1] = SystPercent_Diff * C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetBinContent(bin);
      if(CollisionType==0){
	// ratio
	double A = C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetBinContent(bin);
	double B = C3QSBuiltmerged[EDBin][MBOI][CollisionType]->GetBinContent(bin);
	if(B<1) continue;
	double ratio = A / B;
	Syst_forRatio[bin-1] = pow( ratio * (0.01*(1 - q3/0.1)) ,2);// f coefficients
	Syst_forRatio[bin-1] += pow( ratio * ((0.024-0.276*q3)/2.) ,2);// Interpolator
	Syst_forRatio[bin-1] += pow( ratio * (0.9973 + 0.462*q3 - 14.2*pow(q3,2) - 2.135*pow(q3,3) + 805.3*pow(q3,4) + 1358*pow(q3,5) - 1) ,2);// Tij sign
	Syst_forRatio[bin-1] += pow( ratio * sqrt(pow(0.9847 + 0.358*q3 - 2.133*q3*q3-1,2) + pow(0.9856 + 0.3285*q3 - 1.897*q3*q3-1,2) - 2*(0.9847 + 0.358*q3 - 2.133*q3*q3-1)*(0.9856 + 0.3285*q3 - 1.897*q3*q3-1) ) ,2);// MRC
	Syst_forRatio[bin-1] += pow( ratio * sqrt(pow(0.975 + 0.4189*q3 - 2.055*q3*q3-1,2) + pow(0.9786 + 0.421*q3 - 2.108*q3*q3-1,2) - 2*(0.975 + 0.4189*q3 - 2.055*q3*q3-1)*(0.9786 + 0.421*q3 - 2.108*q3*q3-1) ) ,2);// Muon
	Syst_forRatio[bin-1] += pow( ratio * sqrt(pow(0.936 + 1.194*q3 - 5.912*q3*q3-1,2) + pow(0.946 + 0.849*q3 - 3.316*q3*q3-1,2) - 2*(0.936 + 1.194*q3 - 5.912*q3*q3-1)*(0.946 + 0.849*q3 - 3.316*q3*q3-1) ) ,2);// fc2 scale
	Syst_forRatio[bin-1] = sqrt(Syst_forRatio[bin-1]);
      }
    }
    C3QS_Syst->SetBinContent(1,100); 
    C3QS_Syst->SetMarkerSize(0); C3QS_Syst->SetFillColor(kBlue-10);
    C3QS_Syst->SetLineColor(kBlue-10); C3QS_Syst->SetLineWidth(5);
    C3QS_Syst->SetMarkerColor(kBlue-10);
    if(CollisionType==0) {
      C3QSBuilt_Syst->SetBinContent(1,100); 
      C3QSBuilt_Syst->SetMarkerSize(0); C3QSBuilt_Syst->SetFillColor(kGray); //C3QSBuilt_Syst->SetFillStyle(3004);
      C3QSBuilt_Syst->SetLineColor(kGray); C3QSBuilt_Syst->SetLineWidth(5);
      C3QSBuilt_Syst->SetMarkerColor(kGray); C3QSBuilt_Syst->GetXaxis()->SetRangeUser(0.01,0.2);
      C3QSBuiltmerged[EDBin][MBOI][CollisionType]->SetLineWidth(1.2);
      C3QSBuiltmerged[EDBin][MBOI][CollisionType]->SetLineStyle(2);
      c3QSBuiltmerged[EDBin][MBOI][CollisionType]->SetLineStyle(2);
    }
    C3QS_Syst->GetXaxis()->SetRangeUser(0.01,0.2);
    //
    if(CollisionType==0) C3QSBuiltmerged[EDBin][MBOI][CollisionType]->GetXaxis()->SetRange(2,15);
    C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->SetMaximum(5.);
    C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->SetBarWidth(0.0001);
    C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Draw();
    
    C3QS_Syst->Draw("same");
    if(CollisionType==0) C3QSBuilt_Syst->Draw("same");// E1
    

    C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Draw("same");
    c3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Draw("same");
    if(ChProdBOI==0 && CollisionType==0){ 
      C3QSBuiltmerged[EDBin][MBOI][CollisionType]->Draw("same");
      c3QSBuiltmerged[EDBin][MBOI][CollisionType]->Draw("same");
      //
      TString *proName=new TString("C3QSbuilt_G"); TString *proNameNeg=new TString("C3QSNegbuilt_G");
      TString *proName_c3=new TString("c3QSbuilt_G"); TString *proNameNeg_c3=new TString("c3QSNegbuilt_G");
      TH1D *C3QSbuilt_G = (TH1D*)C3QSBuiltmerged2D[EDBin][MBOI][CollisionType]->ProjectionY(proName->Data(), GbinPlot, GbinPlot);
      TH1D *C3QSNegbuilt_G = (TH1D*)C3QSNegBuiltmerged2D[EDBin][MBOI][CollisionType]->ProjectionY(proNameNeg->Data(), GbinPlot, GbinPlot);
      //TH1D *c3QSbuilt_G = (TH1D*)c3QSBuiltmerged2D[EDBin][MBOI][CollisionType]->ProjectionY(proName_c3->Data(), GbinPlot, GbinPlot);
      //TH1D *c3QSNegbuilt_G = (TH1D*)c3QSNegBuiltmerged2D[EDBin][MBOI][CollisionType]->ProjectionY(proNameNeg_c3->Data(), GbinPlot, GbinPlot);
      proName->Append("_FullWeightDen"); proNameNeg->Append("_FullWeightDen");
      //proName_c3->Append("_FullWeightDen"); proNameNeg_c3->Append("_FullWeightDen");
      TH1D *tempDen = (TH1D*)C3QSBuiltmerged2D[EDBin][MBOI][CollisionType]->ProjectionY(proName->Data(), 4, 4);
      TH1D *tempDenNeg = (TH1D*)C3QSNegBuiltmerged2D[EDBin][MBOI][CollisionType]->ProjectionY(proNameNeg->Data(), 4, 4);
      //TH1D *tempDen_c3 = (TH1D*)c3QSBuiltmerged2D[EDBin][MBOI][CollisionType]->ProjectionY(proName_c3->Data(), 4, 4);
      //TH1D *tempDenNeg_c3 = (TH1D*)c3QSNegBuiltmerged2D[EDBin][MBOI][CollisionType]->ProjectionY(proNameNeg_c3->Data(), 4, 4);
      // Add Pos with Neg weights
      tempDen->Add(tempDenNeg);
      //tempDen_c3->Add(tempDenNeg_c3);
      C3QSbuilt_G->Add(C3QSNegbuilt_G);
      //c3QSbuilt_G->Add(c3QSNegbuilt_G);
      //
      C3QSbuilt_G->Divide(tempDen);
      //c3QSbuilt_G->Divide(tempDen_c3);
      C3QSbuilt_G->SetLineColor(4); //c3QSbuilt_G->SetLineColor(1);
      C3QSbuilt_G->GetXaxis()->SetRange(2,15); //c3QSbuilt_G->GetXaxis()->SetRange(2,15);
      if(ChProdBOI==0) {
	C3QSbuilt_G->Draw("same");
	//c3QSbuilt_G->Draw("same");
      }
    }

    legend1->AddEntry(C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType],"#font[12]{C}_{3}^{QS}","p");
    legend1->AddEntry(c3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType],"#font[12]{#bf{c}}_{3}^{QS}","p");
    if(ChProdBOI==0 && CollisionType==0) {
      legend1->AddEntry(C3QSBuiltmerged[EDBin][MBOI][CollisionType],"Built from #font[12]{C}_{2} (G=0%)","l");
      legend1->AddEntry(C3QSbuilt_G,BuiltNameC2->Data(),"l");
      legend1->AddEntry(C3QSBuilt_Syst,"Built #font[12]{C}_{3}^{QS} systematic uncertainty","l");
    }  
    legend1->AddEntry(C3QS_Syst,"#font[12]{C}_{3}^{QS} systematic uncertainty","l");
    legend1->Draw("same");
    
    Unity->Draw("same");
    
    ALICEspecif->Draw("same");
    //
    TString *CentAndED = new TString("Centrality 0-5%");
    if(EDBin==0) CentAndED->Append(", 0.16<#font[12]{K}_{T3}<0.3 GeV/#font[12]{c}");
    else CentAndED->Append(", 0.3<#font[12]{K}_{T3}<1.0 GeV/#font[12]{c}");
    TLatex *Cent = new TLatex(0.03,4.8,CentAndED->Data());
    Cent->SetTextFont(TextFont);
    Cent->SetTextSize(SizeSpecif);
    Cent->Draw("same");
    
    
    
    if(ChProdBOI==0 && CollisionType==0){
      pad1_2->cd();
      gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.04);
      gPad->SetTopMargin(0.03); gPad->SetBottomMargin(0.32);

      TH1D *Ratio_3 = (TH1D*)C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Clone();
      Ratio_3->GetXaxis()->SetTitleSize(2.2*SizeTitle);
      Ratio_3->GetXaxis()->SetLabelSize(2.2*SizeLabel);
      Ratio_3->GetYaxis()->SetTitleSize(2.2*SizeTitle);
      Ratio_3->GetYaxis()->SetLabelSize(2.2*SizeLabel);
      Ratio_3->GetXaxis()->SetTitleOffset(1.05);
      Ratio_3->GetYaxis()->SetTitleOffset(0.5);
      Ratio_3->GetXaxis()->SetNdivisions(606);
      Ratio_3->GetYaxis()->SetNdivisions(202);
      Ratio_3->GetYaxis()->SetTitle("#font[12]{C_{3}^{QS}} / built ");
      Ratio_3->SetMinimum(0.85); Ratio_3->SetMaximum(1.05);
      TH1D *Ratio_3_G=(TH1D*)Ratio_3->Clone();
      Ratio_3->SetMarkerStyle(24);
      Ratio_3_G->SetMarkerStyle(20);
      Ratio_3->Divide(C3QSBuiltmerged[EDBin][MBOI][CollisionType]);
      Ratio_3_G->Divide(C3QSbuilt_G);
      TH1D *Ratio_3_Syst = (TH1D*)Ratio_3->Clone();
      for(int bin=1; bin<=15; bin++){
	Ratio_3_Syst->SetBinError(bin, Syst_forRatio[bin-1]);
      }
      Ratio_3_Syst->SetMarkerSize(0); Ratio_3_Syst->SetMarkerColor(kBlue-10); Ratio_3_Syst->SetLineColor(kBlue-10);
      Ratio_3_Syst->SetLineWidth(5);
      Ratio_3_Syst->Draw();
      Ratio_3->Draw("same");

      Ratio_3_G->Draw("same");

      Unity->Draw("same");
      
      legend1_2->AddEntry(Ratio_3,"G=0%","p");
      legend1_2->AddEntry(Ratio_3_G,BuiltNameC2->Data(),"p");
      legend1_2->Draw("same");
      
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
      
      
      
      TH1D *chi2_3 = new TH1D("chi2_3","",100,-0.5,99.5);
      
      chi2_3->SetLineColor(4);
      chi2_3->SetMarkerColor(4);
      chi2_3->GetXaxis()->SetTitle("Coherent fraction (%)"); chi2_3->GetYaxis()->SetTitle("#sqrt{#chi^{2}}");
      chi2_3->GetXaxis()->SetTitleSize(SizeTitle);  chi2_3->GetYaxis()->SetTitleSize(SizeTitle);
      chi2_3->GetXaxis()->SetLabelSize(SizeLabel);  chi2_3->GetYaxis()->SetLabelSize(SizeLabel);
      TH2D *chi2_2D_3[2];// Stat then Syst
      chi2_2D_3[0] = new TH2D("chi2Stat_2D_3","",5,0.5,5.5, 100,-0.5,99.5);
      chi2_2D_3[1] = new TH2D("chi2Syst_2D_3","",5,0.5,5.5, 100,-0.5,99.5);
      
      
      TH1D *tempDen = (TH1D*)C3QSBuiltmerged2D[EDBin][MBOI][CollisionType]->ProjectionY("TPFullWeight3_Den", 4, 4);
      TH1D *tempDenNeg = (TH1D*)C3QSNegBuiltmerged2D[EDBin][MBOI][CollisionType]->ProjectionY("TPNegFullWeight3_Den", 4, 4);
      tempDen->Add(tempDenNeg);// Add Pos and Neg Den
      

      for(int binG=Gbinstart; binG<=Gbinstart+25; binG++){
	TString *proName=new TString("TPFullWeight3_");
	*proName += binG;
	TH1D *tempNum = (TH1D*)C3QSBuiltmerged2D[EDBin][MBOI][CollisionType]->ProjectionY(proName->Data(), binG, binG);
	proName->Append("_Neg");
	TH1D *tempNumNeg = (TH1D*)C3QSNegBuiltmerged2D[EDBin][MBOI][CollisionType]->ProjectionY(proName->Data(), binG, binG);
	//
	// Add Pos and Neg Num
	tempNum->Add(tempNumNeg);
	//
	//tempNum->Add(tempDen);// now done in Plot_FourPion.C
	tempNum->Divide(tempDen);
	
	double value = fabs(C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetBinContent(Q3binChi2) - tempNum->GetBinContent(Q3binChi2));
	double err = pow(C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetBinError(Q3binChi2),2);// stat
	err += pow(Syst_forChi2_3[Q3binChi2-1],2);// syst
	err = sqrt(err);
	if(err <=0) continue;
	double Chi2 = pow(value / err,2);
	//
	
	chi2_3->SetBinContent(1 + 2*(binG-Gbinstart), sqrt(Chi2)); 
	chi2_3->SetBinError(1 + 2*(binG-Gbinstart), 0.001);
	
	//
	Chi2=0;
	double NDF=0;
	for(int binQ3=2; binQ3<=5; binQ3++){
	  if(tempNum->GetBinContent(binQ3) <=0) continue;
	  if(C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetBinError(binQ3) <= 0) continue;
	  double value = fabs(C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetBinContent(binQ3) - tempNum->GetBinContent(binQ3));
	  double err = pow(C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetBinError(binQ3),2);// stat
	  double Chi2 = pow(value / sqrt(err),2);
	  chi2_2D_3[0]->SetBinContent(binQ3, 1+(binG-Gbinstart), sqrt(Chi2));// was value/sqrt(err)
	  //
	  err = pow(Syst_forChi2_3[binQ3-1],2);// syst
	  Chi2 = pow(value / sqrt(err),2);
	  chi2_2D_3[1]->SetBinContent(binQ3, 1+(binG-Gbinstart), sqrt(Chi2));
	  
	}
      }
      chi2_3->SetMarkerStyle(20);
      chi2_3->SetMinimum(0); chi2_3->SetMaximum(5); 
      chi2_3->Draw();
      TString *Q3binName = new TString("0.0");
      *Q3binName += Q3binChi2-1;
      Q3binName->Append(" < #font[12]{Q_{3}} < 0.0");
      *Q3binName += Q3binChi2;
      Q3binName->Append(" GeV/#font[12]{c}");
      legend2->SetHeader(Q3binName->Data());
      legend2->AddEntry(chi2_3,"R_{coh}=0","p");
      legend2->Draw("same");
      
      TString *meanpTName3 = new TString("#LT #font[12]{p}_{T} #GT = 0.");
      *meanpTName3 += Q3_meanpT[EDBin][Q3binChi2-1];
      meanpTName3->Append(" GeV/#font[12]{c}");
      TLatex *Specif_pT3 = new TLatex(0.17,0.9,meanpTName3->Data());
      Specif_pT3->SetNDC();
      Specif_pT3->SetTextFont(TextFont);
      Specif_pT3->SetTextSize(SizeSpecif);
      Specif_pT3->Draw("same");
      
      
      TString *SaveNameChi2_3 = new TString("ChiSq_C3_bin");
      *SaveNameChi2_3 += Q3binChi2;
      SaveNameChi2_3->Append("_K");
      *SaveNameChi2_3 += EDBin;
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
      gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.05);
      gPad->SetTopMargin(0.03); gPad->SetBottomMargin(0.14);
      TLegend *legend2_2 = new TLegend(.15,.15, .35,.35,NULL,"brNDC");//.45 or .4 for x1
      legend2_2->SetBorderSize(0);
      legend2_2->SetFillColor(0);
      legend2_2->SetTextFont(TextFont);
      legend2_2->SetTextSize(SizeLegend);
      
      TH1D *GversusQ3[2];// Stat then Syst
      GversusQ3[0] = new TH1D("GversusQ3_Stat","",5,0,0.05);
      GversusQ3[1] = new TH1D("GversusQ3_Syst","",5,0,0.05);
      
      for(int ErrType=0; ErrType<2; ErrType++){
	for(int binQ3=2; binQ3<=5; binQ3++){
	  double minG = 0;
	  double minG_e1=0, minG_e2=0;
	  double minChi=100;
	  
	  for(int binG=1; binG<=25; binG++){// min
	    if(minChi > chi2_2D_3[ErrType]->GetBinContent(binQ3, binG)) {
	      minChi = chi2_2D_3[ErrType]->GetBinContent(binQ3, binG);
	      minG = 2*(binG-1);
	    }
	  }
	  
	  for(int binG=1; binG<=25; binG++){// error
	    if(minG > 0) {
	      if(fabs(minChi - chi2_2D_3[ErrType]->GetBinContent(binQ3, binG)) < 1.) {
		if(minG>2*(binG-1)) minG_e1 = fabs(minG - 2*(binG-1)); 
		else minG_e2 = fabs(minG - 2*(binG-1));
	      }
	    }else{
	      if(fabs(minChi - chi2_2D_3[ErrType]->GetBinContent(binQ3, binG)) < 1.) {
		minG_e1 = fabs(minG - 2*(binG-1)); 
	      }
	    }
	  }
	  GversusQ3[ErrType]->SetBinContent(binQ3, minG);
	  if(minG_e1 > minG_e2) GversusQ3[ErrType]->SetBinError(binQ3, minG_e1);
	  else GversusQ3[ErrType]->SetBinError(binQ3, minG_e2);
	}
      }// Err Type
      //
      GversusQ3[0]->SetMarkerStyle(20); GversusQ3[0]->SetMarkerColor(4); GversusQ3[0]->SetLineColor(4);
      GversusQ3[1]->SetMarkerSize(0); GversusQ3[1]->SetMarkerColor(kBlue-10); GversusQ3[1]->SetLineColor(kBlue);
      GversusQ3[1]->SetFillStyle(0);
      GversusQ3[0]->SetMinimum(0); GversusQ3[0]->SetMaximum(100); 
      GversusQ3[0]->GetXaxis()->SetTitle("#font[12]{Q_{3}} (GeV/#font[12]{c})"); GversusQ3[0]->GetYaxis()->SetTitle("Coherent fraction (%)");
      GversusQ3[0]->GetXaxis()->SetTitleSize(SizeTitle);  GversusQ3[0]->GetYaxis()->SetTitleSize(SizeTitle);
      GversusQ3[0]->GetXaxis()->SetLabelSize(SizeLabel);  GversusQ3[0]->GetYaxis()->SetLabelSize(SizeLabel);
      GversusQ3[0]->GetXaxis()->SetNdivisions(505); GversusQ3[0]->GetYaxis()->SetNdivisions(505);
      GversusQ3[0]->SetBinContent(1,200); GversusQ3[1]->SetBinContent(1,200);
      GversusQ3[0]->Draw();
      GversusQ3[1]->Draw("E2 same");
      GversusQ3[0]->Draw("same");
      //

      
      //
      ALICEspecif->Draw("same");
      //
      TLatex *Cent_2 = new TLatex(0.022,77,"Centrality 0-5%");
      Cent_2->SetTextFont(TextFont);
      Cent_2->SetTextSize(SizeSpecif);
      if(CollisionType==0) Cent_2->Draw("same");
      TString *KTname_2 = new TString("");
      if(EDBin==0) KTname_2->Append("0.16<#font[12]{K}_{T3}<0.3 GeV/#font[12]{c}");
      else KTname_2->Append("0.3<#font[12]{K}_{T3}<1.0 GeV/#font[12]{c}");
      TLatex *KT_2 = new TLatex(0.022,69,KTname_2->Data());
      KT_2->SetTextFont(TextFont);
      KT_2->SetTextSize(SizeSpecif);
      KT_2->Draw("same");
      //
      legend2_2->AddEntry(GversusQ3[0],"R_{coh}=0","p");
      legend2_2->Draw("same");
      
      
    }// ChProdBOI==0
  }// ChProdBOI!=2
  

  
  
  //////////////////////////////////////////////////////////////////////////
  // 4-pion
  TCanvas *can3 = new TCanvas("can3", "can3",10,700,700,Ysize);// 11,53,700,500
  can3->SetHighLightColor(2);
  gStyle->SetOptFit(0111);
  can3->SetFillColor(0);//10
  can3->SetBorderMode(0);
  can3->SetBorderSize(2);
  can3->SetFrameFillColor(0);
  can3->SetFrameBorderMode(0);
  can3->SetFrameBorderMode(0);
  
  TPad *pad3 = new TPad("pad3","pad3",0.,Ystart,1.,1.);
  TPad *pad3_2 = new TPad("pad3_2","pad3_2",0.,0.,1.,Ystart);
  pad3->Draw();
  if(ChProdBOI==0) pad3_2->Draw();
  pad3->cd();
  gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.03); gPad->SetBottomMargin(Xmarge);

  TLegend *legend3= new TLegend(.6,.35, .82,.78,NULL,"brNDC");//.45 or .4 for x1
  legend3->SetBorderSize(0);
  legend3->SetFillColor(0);
  legend3->SetTextFont(TextFont);
  legend3->SetTextSize(SizeLegend);
  if(ChProdBOI!=0) {legend3->SetY1(0.55); legend3->SetY2(0.75);}
  TLegend *legend3_2=(TLegend*)legend3->Clone();
  legend3_2->SetX1(0.52); legend3_2->SetX2(0.9); legend3_2->SetY1(0.35); legend3_2->SetY2(0.6); 
  legend3_2->SetTextSize(2.5*SizeLegend);

  C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->SetTitleSize(SizeTitle);
  C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->SetLabelSize(SizeLabel);
  C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetYaxis()->SetTitleSize(SizeTitle);
  C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetYaxis()->SetLabelSize(SizeLabel);
  C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->SetTitleOffset(1.05);
  C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetYaxis()->SetTitleOffset(1.1);
  C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->SetNdivisions(606);
  C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetYaxis()->SetNdivisions(505);
 
  //
  TString *proName4=new TString("C4QSbuilt_G"); TString *proNameNeg4=new TString("C4QSNegbuilt_G");
  if(CollisionType==0){
    TH1D *C4QSbuilt_G = (TH1D*)C4QSBuiltmerged2D[EDBin][MBOI][CollisionType]->ProjectionY(proName4->Data(), GbinPlot, GbinPlot);
    TH1D *C4QSNegbuilt_G = (TH1D*)C4QSNegBuiltmerged2D[EDBin][MBOI][CollisionType]->ProjectionY(proNameNeg4->Data(), GbinPlot, GbinPlot);
    TH1D *a4QSbuilt_G = (TH1D*)a4QSBuiltmerged2D[EDBin][MBOI][CollisionType]->ProjectionY("a4_G", GbinPlot, GbinPlot);
    TH1D *a4QSNegbuilt_G = (TH1D*)a4QSNegBuiltmerged2D[EDBin][MBOI][CollisionType]->ProjectionY("a4Neg_G", GbinPlot, GbinPlot);
    TH1D *b4QSbuilt_G = (TH1D*)b4QSBuiltmerged2D[EDBin][MBOI][CollisionType]->ProjectionY("b4_G", GbinPlot, GbinPlot);
    TH1D *b4QSNegbuilt_G = (TH1D*)b4QSNegBuiltmerged2D[EDBin][MBOI][CollisionType]->ProjectionY("b4Neg_G", GbinPlot, GbinPlot);
    TH1D *c4QSbuilt_G = (TH1D*)c4QSBuiltmerged2D[EDBin][MBOI][CollisionType]->ProjectionY("c4_G", GbinPlot, GbinPlot);
    TH1D *c4QSNegbuilt_G = (TH1D*)c4QSNegBuiltmerged2D[EDBin][MBOI][CollisionType]->ProjectionY("c4Neg_G", GbinPlot, GbinPlot);

    proName4->Append("_FullWeightDen"); proNameNeg4->Append("_FullWeightDen");
    TH1D *tempDen4 = (TH1D*)C4QSBuiltmerged2D[EDBin][MBOI][CollisionType]->ProjectionY(proName4->Data(), 4, 4);
    TH1D *tempDenNeg4 = (TH1D*)C4QSNegBuiltmerged2D[EDBin][MBOI][CollisionType]->ProjectionY(proNameNeg4->Data(), 4, 4);
    // Add Pos with Neg weights
    tempDen4->Add(tempDenNeg4);
    C4QSbuilt_G->Add(C4QSNegbuilt_G);
    a4QSbuilt_G->Add(a4QSNegbuilt_G);
    b4QSbuilt_G->Add(b4QSNegbuilt_G);
    c4QSbuilt_G->Add(c4QSNegbuilt_G);
    C4QSbuilt_G->Divide(tempDen4);
    a4QSbuilt_G->Divide(tempDen4);
    b4QSbuilt_G->Divide(tempDen4);
    c4QSbuilt_G->Divide(tempDen4);
    if(ReNormBuiltBaseline){
      int BinL = C4QSbuilt_G->GetXaxis()->FindBin(ReNormL_4);
      int BinH = C4QSbuilt_G->GetXaxis()->FindBin(ReNormH_4);
      C4QSbuilt_G->Scale( C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Integral(BinL,BinH) / C4QSbuilt_G->Integral(BinL,BinH));
      a4QSbuilt_G->Scale( a4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Integral(BinL,BinH) / a4QSbuilt_G->Integral(BinL,BinH));
      b4QSbuilt_G->Scale( b4QSmerged[EDBin][MBOI][CollisionType]->Integral(BinL,BinH) / b4QSbuilt_G->Integral(BinL,BinH));
      c4QSbuilt_G->Scale( c4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Integral(BinL,BinH) / c4QSbuilt_G->Integral(BinL,BinH));
      C4QSbuilt_G->SetLineColor(4);
      a4QSbuilt_G->SetLineColor(2);
      b4QSbuilt_G->SetLineColor(6);
      c4QSbuilt_G->SetLineColor(1);
      a4QSbuilt_G->SetBinContent(2,0);
      b4QSbuilt_G->SetBinContent(2,0);
      c4QSbuilt_G->SetBinContent(2,0);
    }
    for(int bin=1; bin<=C4QSbuilt_G->GetNbinsX(); bin++) {
      C4QSbuilt_G->SetBinError(bin,0);
      a4QSbuilt_G->SetBinError(bin,0);
      b4QSbuilt_G->SetBinError(bin,0);
      c4QSbuilt_G->SetBinError(bin,0);
    }
  }
  if(ReNormBuiltBaseline){
    int BinL = C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->FindBin(ReNormL_4);
    int BinH = C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->FindBin(ReNormH_4);
    C4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->Scale(C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Integral(BinL,BinH) / C4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->Integral(BinL,BinH));
    C4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->Scale(C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Integral(BinL,BinH) / C4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->Integral(BinL,BinH));
    b4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->Scale(b4QSmerged[EDBin][MBOI][CollisionType]->Integral(BinL,BinH) / b4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->Integral(BinL,BinH));
    b4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->Scale(b4QSmerged[EDBin][MBOI][CollisionType]->Integral(BinL,BinH) / b4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->Integral(BinL,BinH));
    c4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->Scale(c4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Integral(BinL,BinH) / c4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->Integral(BinL,BinH));
    c4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->Scale(c4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Integral(BinL,BinH) / c4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->Integral(BinL,BinH));
  }
  for(int bin=1; bin<=C4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->GetNbinsX(); bin++){
    C4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->SetBinError(bin,0);
    C4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->SetBinError(bin,0);
    a4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->SetBinError(bin,0);
    a4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->SetBinError(bin,0);
    b4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->SetBinError(bin,0);
    b4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->SetBinError(bin,0);
    c4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->SetBinError(bin,0);
    c4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->SetBinError(bin,0);
  }
  //
  if(EDBin==1 && ChProdBOI!=0){
    C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->SetBinContent(2,100);
    a4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->SetBinContent(2,100);
    c4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->SetBinContent(2,100);
  }
  if(ChProdBOI==0 && CollisionType!=0){
    C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->SetBinContent(3,100);
    a4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->SetBinContent(3,100);
    b4QSmerged[EDBin][MBOI][CollisionType]->SetBinContent(3,100);
    c4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->SetBinContent(3,100);
    if(CollisionType==2){
      C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->SetBinContent(4,100);
      a4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->SetBinContent(4,100);
      b4QSmerged[EDBin][MBOI][CollisionType]->SetBinContent(4,100);
      c4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->SetBinContent(4,100);
    }
  }
  if(ChProdBOI==0){
    C4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->SetBinContent(2,0);
    a4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->SetBinContent(2,0);
    b4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->SetBinContent(2,0);
    c4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->SetBinContent(2,0);
    C4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->SetBinContent(2,0);
    a4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->SetBinContent(2,0);
    b4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->SetBinContent(2,0);
    c4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->SetBinContent(2,0);
    if(CollisionType!=0){
      C4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->SetBinContent(3,0);
      a4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->SetBinContent(3,0);
      b4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->SetBinContent(3,0);
      c4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->SetBinContent(3,0);
      C4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->SetBinContent(3,0);
      a4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->SetBinContent(3,0);
      b4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->SetBinContent(3,0);
      c4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->SetBinContent(3,0);
    }
    if(CollisionType==2){
      C4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->SetBinContent(4,0);
      a4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->SetBinContent(4,0);
      b4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->SetBinContent(4,0);
      c4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->SetBinContent(4,0);
      C4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->SetBinContent(4,0);
      a4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->SetBinContent(4,0);
      b4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->SetBinContent(4,0);
      c4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->SetBinContent(4,0);
    }
  }
  if(ChProdBOI!=0) C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->SetMinimum(0.8);
  C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->SetMaximum(8.4);
  if(CollisionType==1) C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->SetMaximum(15.0);
  if(CollisionType==2) C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->SetMaximum(19.0);
  if(ChProdBOI==1) C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->SetMaximum(5.0);
  if(ChProdBOI==2) C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->SetMaximum(3.2);
  if(CollisionType==0){
    C4QSBuiltmerged[EDBin][MBOI][CollisionType]->GetXaxis()->SetRange(3,15);
    C4QSbuilt_G->GetXaxis()->SetRange(3,15);
  }
  
  //
  TH1D *C4QS_Syst;
  double shift = 0.002;
  if(CollisionType!=0) shift=0.007;
  if(ChProdBOI==0) C4QS_Syst = new TH1D("C4QS_Syst","",C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetNbinsX(),C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->GetBinLowEdge(1)+shift, C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->GetBinUpEdge(C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetNbinsX())+shift); 
  else C4QS_Syst = (TH1D*)C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Clone();
  TH1D *a4QS_Syst = (TH1D*)a4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Clone();
  TH1D *c4QS_Syst = (TH1D*)c4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Clone();
  TH1D *C4QSBuilt_Syst;
  if(ChProdBOI==0){
    if(CollisionType==0) C4QSBuilt_Syst = (TH1D*)C4QSBuiltmerged[EDBin][MBOI][CollisionType]->Clone();
    else C4QSBuilt_Syst = (TH1D*)C4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->Clone();
  }
  for(int bin=1; bin<=C4QS_Syst->GetNbinsX(); bin++){
    double q4 = C4QS_Syst->GetXaxis()->GetBinCenter(bin);
    if(CollisionType!=0) q4 /= 3.3;// rescale q4 for pp and p-Pb
    if(ChProdBOI==0) {
      if(CollisionType==0) C4QS_Syst->SetBinContent(bin, 7.75);
      else if(CollisionType==1) C4QS_Syst->SetBinContent(bin, 14.);
      else C4QS_Syst->SetBinContent(bin, 17.5);
      double syst1 = pow(0.001,2);// cc
      syst1 += pow(0.004 - 0.004*q4/0.18,2);// 11h to 10h
      syst1 += pow(0.01*(1 - q4/0.18),2);// f coefficients, r*<70. was pow(0.9975 - 0.09*q4 - 1,2)
      syst1 += pow(0.9814 + 0.2471*q4 - 0.8312*q4*q4 - 1,2);// MRC
      syst1 += pow(0.9635 + 0.3475*q4 - 0.9729*q4*q4 - 1,2);// Muon, 92%
      syst1 += pow(0.900 + 1.126*q4 - 3.354*q4*q4 - 1,2);// fc2 scale
      //syst1 += pow(0.125*exp(-61.38*q4),2);// K factorization
      //syst1 +=  pow( 0.028 + 0.254*q4,2);// MC1 cumulant residue;
      syst1 = sqrt(syst1);
      syst1 *= C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetBinContent(bin);
      C4QS_Syst->SetBinError(bin, syst1);
      // Built
      C4QSBuilt_Syst->SetBinContent(bin, C4QS_Syst->GetBinContent(bin));
      double syst2 = pow(0.004 - 0.004*q4/0.18,2);// 11h to 10h
      syst2 += pow(0.9793 + 0.2857*q4 - 0.9888*q4*q4 - 1,2);// MRC
      syst2 += pow(0.9725 + 0.2991*q4 - 0.8058*q4*q4 - 1,2);// Muon, 92%
      syst2 += pow(0.905 + 1.03*q4 - 2.977*q4*q4 - 1,2);// fc2 scale
      syst2 += pow((0.032-0.245*q4)/2.,2);// Interpolator. was pow(0.0379*exp(-42.82*q4),2)
      if(EDBin==0) syst2 += pow( 0.33*(1.005 - 0.437*q4 + 13.28*pow(q4,2) - 174.3*pow(q4,3) + 970.3*pow(q4,4) - 1920*pow(q4,5) - 1),2);// Tij sign
      else syst2 += pow( 0.33*(0.906 + 6.207*q4 - 135*pow(q4,2) + 1133*pow(q4,3) - 3969*pow(q4,4) + 4717*pow(q4,5) - 1),2);// Tij sign
      syst2 = sqrt(syst2);
      if(CollisionType==0) syst2 *= C4QSBuiltmerged[EDBin][MBOI][CollisionType]->GetBinContent(bin);
      else syst2 *= C4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->GetBinContent(bin);
      C4QSBuilt_Syst->SetBinError(bin, syst2);
    }else{// same % systematics for MC1 and MC2
      double syst1 = pow(0.002,2);// cc
      syst1 += pow(0.002 - 0.002*q4/0.18,2);// 11h to 10h
      syst1 += pow(0.01*(1 - q4/0.18),2);// f coefficients, was pow(0.005,2)
      syst1 += pow(0.0417*exp(-34.1*q4),2);// MRC, 10% of full diff
      syst1 += pow(0.9713 + 0.2648*q4 - 0.752*q4*q4 - 1,2);// Muon, 92%
      syst1 += pow(0.908 + 1.118*q4 - 3.612*q4*q4 - 1,2);// fc2 scale
      //if(ChProdBOI==1) syst1 += pow(0.125*exp(-61.38*q4),2);// K factorization
      //else syst1 += pow(-0.008 + 0.102*q4 - 0.329*q4*q4,2);// K factorization
      //
      C4QS_Syst->SetBinError(bin, sqrt(syst1) * C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetBinContent(bin));
      syst1 += pow(0.536*exp(-58.21*q4),2);
      a4QS_Syst->SetBinError(bin, sqrt(syst1) * a4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetBinContent(bin));
      c4QS_Syst->SetBinError(bin, sqrt(syst1) * c4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetBinContent(bin));
    }
  }
  C4QS_Syst->SetMarkerSize(0); C4QS_Syst->SetLineColor(kBlue-10); C4QS_Syst->SetMarkerColor(kBlue-10);
  C4QS_Syst->SetLineWidth(5); a4QS_Syst->SetLineWidth(5); c4QS_Syst->SetLineWidth(5);
  if(ChProdBOI==0) {
    C4QSBuilt_Syst->SetMarkerSize(0); C4QSBuilt_Syst->SetLineColor(kGray); C4QSBuilt_Syst->SetMarkerColor(kGray); //C4QSBuilt_Syst->SetFillStyle(3004);
    C4QSBuilt_Syst->SetLineWidth(5); C4QSBuilt_Syst->SetLineStyle(1);
  }
  c4QS_Syst->SetMarkerSize(0); c4QS_Syst->SetFillColor(kGray); c4QS_Syst->SetMarkerColor(kGray); c4QS_Syst->SetLineColor(kGray);
  a4QS_Syst->SetMarkerSize(0); a4QS_Syst->SetFillColor(kRed-10); a4QS_Syst->SetMarkerColor(kRed-10); a4QS_Syst->SetLineColor(kRed-10);
  
  if(ChProdBOI==0 && CollisionType==0) {
    C4QS_Syst->GetXaxis()->SetRangeUser(0.035,0.2); 
    C4QSBuilt_Syst->GetXaxis()->SetRangeUser(0.035,0.2);
  }
  if(ChProdBOI==0 && CollisionType==1){
    C4QS_Syst->GetXaxis()->SetRangeUser(0.1,0.5);
    C4QSBuilt_Syst->GetXaxis()->SetRangeUser(0.1,0.5);
  }
  if(ChProdBOI==0 && CollisionType==2){
    C4QS_Syst->GetXaxis()->SetRangeUser(0.13,0.5);
    C4QSBuilt_Syst->GetXaxis()->SetRangeUser(0.13,0.5);
  }
  
  double Syst_forChi2_4[15]={0};
  double Syst_forRatio4[15]={0};
  for(int bin=1; bin<=15; bin++){
    double q4 = C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->GetBinCenter(bin);
    if(CollisionType!=0) q4 /= 3.3;// rescale q4 for pp and p-Pb
    int SystBin = C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->FindBin(q4);
    double SystPercent_Diff = 0;
    //SystPercent_Diff += pow(0.125*exp(-61.38*q4),2);// K factorization
    SystPercent_Diff += pow(0.01*(1 - q4/0.18),2);// f coefficients was pow(0.9989 - 0.025*q4 - 1,2)
    SystPercent_Diff += pow((0.032-0.245*q4)/2.,2);// Interpolator. was pow(0.0379*exp(-42.82*q4),2)
    if(EDBin==0) SystPercent_Diff += pow( 0.33*(1.005 - 0.437*q4 + 13.28*pow(q4,2) - 174.3*pow(q4,3) + 970.3*pow(q4,4) - 1920*pow(q4,5) - 1),2);// Tij sign
    else SystPercent_Diff += pow( 0.33*(0.906 + 6.207*q4 - 135*pow(q4,2) + 1133*pow(q4,3) - 3969*pow(q4,4) + 4717*pow(q4,5) - 1),2);// Tij sign
    SystPercent_Diff += pow( (0.9814 + 0.2471*q4 - 0.8312*q4*q4) - (0.9793 + 0.2857*q4 - 0.9888*q4*q4),2);// MRC diff
    SystPercent_Diff += pow( (0.9635 + 0.3475*q4 - 0.9729*q4*q4) - (0.9725 + 0.2991*q4 - 0.8058*q4*q4),2);// Muon diff
    SystPercent_Diff += pow( (0.900 + 1.126*q4 - 3.354*q4*q4) - (0.905 + 1.03*q4 - 2.977*q4*q4),2);// fc2 scale
    //SystPercent_Diff += pow( 0.028 + 0.254*q4,2);// MC1 cumulant residue;
    SystPercent_Diff = sqrt(SystPercent_Diff);
    Syst_forChi2_4[bin-1] = SystPercent_Diff * C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetBinContent(bin);
    // ratio
    if(ChProdBOI==0){
      double A = C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetBinContent(bin);
      double B = 0;
      if(CollisionType==0 && !FitBuild) B = C4QSBuiltmerged[EDBin][MBOI][CollisionType]->GetBinContent(bin);
      else B = C4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->GetBinContent(bin);
      if(B<1) continue;
      double ratio = A / B;
      Syst_forRatio4[bin-1] = pow( ratio * (0.01*(1 - q4/0.18)) ,2);// f coefficients
      Syst_forRatio4[bin-1] += pow( ratio * (0.032-0.245*q4)/2. ,2);// Interpolator
      Syst_forRatio4[bin-1] += pow( ratio * (1.005 - 0.437*q4 + 13.28*pow(q4,2) - 174.3*pow(q4,3) + 970.3*pow(q4,4) - 1920*pow(q4,5) - 1) ,2);// Tij sign
      Syst_forRatio4[bin-1] += pow( ratio * sqrt(pow(0.9814 + 0.2471*q4 - 0.8312*q4*q4-1,2) + pow(0.9793 + 0.2857*q4 - 0.9888*q4*q4-1,2) - 2*(0.9814 + 0.2471*q4 - 0.8312*q4*q4-1)*(0.9793 + 0.2857*q4 - 0.9888*q4*q4-1) ) ,2);// MRC
      Syst_forRatio4[bin-1] += pow( ratio * sqrt(pow(0.9635 + 0.3475*q4 - 0.9729*q4*q4 - 1,2) + pow(0.9725 + 0.2991*q4 - 0.8058*q4*q4 - 1,2) - 2*(0.9635 + 0.3475*q4 - 0.9729*q4*q4 - 1)*(0.9725 + 0.2991*q4 - 0.8058*q4*q4 - 1) ) ,2);// Muon
      Syst_forRatio4[bin-1] += pow( ratio * sqrt(pow(0.900 + 1.126*q4 - 3.354*q4*q4 - 1,2) + pow(0.905 + 1.03*q4 - 2.977*q4*q4 - 1,2) - 2*(0.900 + 1.126*q4 - 3.354*q4*q4 - 1)*(0.905 + 1.03*q4 - 2.977*q4*q4 - 1) ) ,2);// fc2 scale
      Syst_forRatio4[bin-1] = sqrt(Syst_forRatio4[bin-1]);
    }
  }

  //
  
  
  C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->SetBinContent(1,100); C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->SetBinError(1,100);
  C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->SetRangeUser(0,0.16);
  if(CollisionType!=0) C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->SetRangeUser(0,0.5);
  C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Draw();
  
  C4QS_Syst->Draw("same");
  if(ChProdBOI==0) C4QSBuilt_Syst->Draw("same");// E1
  if(ChProdBOI!=0){    
    a4QS_Syst->Draw("same");
    c4QS_Syst->Draw("same");
  }
  
  C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Draw("same");
  a4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Draw("same");
  if(ChProdBOI==0) b4QSmerged[EDBin][MBOI][CollisionType]->Draw("same");
  c4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Draw("same");
  if(ChProdBOI==0) {
    if(FitBuild==kFALSE && CollisionType==0){
      C4QSBuiltmerged[EDBin][MBOI][CollisionType]->SetLineWidth(1.2);
      C4QSBuiltmerged[EDBin][MBOI][CollisionType]->Draw("same");
      a4QSBuiltmerged[EDBin][MBOI][CollisionType]->Draw("same");
      b4QSBuiltmerged[EDBin][MBOI][CollisionType]->Draw("same");
      c4QSBuiltmerged[EDBin][MBOI][CollisionType]->Draw("same");
    }else{
      C4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->Draw("same");
      C4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->Draw("same");
      a4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->Draw("same");
      a4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->Draw("same");
      b4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->Draw("same");
      b4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->Draw("same");
      c4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->Draw("same");
      c4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->Draw("same");
      //
      C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Draw("same");
      a4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Draw("same");
      b4QSmerged[EDBin][MBOI][CollisionType]->Draw("same");
      c4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Draw("same");
    }  
  }
  
  
  if(ChProdBOI==0 && CollisionType==0 && FitBuild==0) {
    C4QSbuilt_G->Draw("same");
    a4QSbuilt_G->Draw("same");
    b4QSbuilt_G->Draw("same");
    c4QSbuilt_G->Draw("same");
  }

  legend3->AddEntry(C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType],"#font[12]{C}_{4}^{QS}","p");
  legend3->AddEntry(a4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType],"#font[12]{#bf{a}}_{4}^{QS}","p");
  if(ChProdBOI==0) legend3->AddEntry(b4QSmerged[EDBin][MBOI][CollisionType],"#font[12]{#bf{b}}_{4}^{QS}","p");
  legend3->AddEntry(c4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType],"#font[12]{#bf{c}}_{4}^{QS}","p");
  if(ChProdBOI==0 && CollisionType==0 && !FitBuild) {
    TH1D *Built1_forLegend = (TH1D*)C4QSBuiltmerged[EDBin][MBOI][CollisionType]->Clone();
    TH1D *Built2_forLegend = (TH1D*)C4QSbuilt_G->Clone();
    Built1_forLegend->SetLineColor(1);
    legend3->AddEntry(Built1_forLegend,"Built from #font[12]{C}_{2} (G=0%)","l");
    Built2_forLegend->SetLineColor(1);
    legend3->AddEntry(Built2_forLegend,BuiltNameC2->Data(),"l");
    //legend3->AddEntry(C4QSbuilt_G,"Built #font[12]{C}_{4}^{QS} (G=30%, R_{coh}=R_{ch})","l");
    //legend3->AddEntry(C4QSBuilt_Syst,"Built #font[12]{C}_{4}^{QS} systematic uncertainty","f");
  }
  if(ChProdBOI==0){
    if(FitBuild || CollisionType!=0){
      legend3->AddEntry(C4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType], BuiltNamec3->Data(),"l");
      legend3->AddEntry(C4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType], BuiltNameC3->Data(),"l");
    }
  }
  //legend3->AddEntry(C4QS_Syst,"#font[12]{C}_{4}^{QS} systematic uncertainty","f");
  
  ALICEspecif->Draw("same");
  if(CollisionType==0) Centspecif->Draw("same");
  KTspecif->Draw("same");

  legend3->Draw("same");
  Unity->Draw("same");

  /*TF1 *Gauss_c4Fit=new TF1("Gauss_c4Fit","[0]*(1+[1]*exp(-pow(x*[2]/0.19733,2)/3.))",0,1);
  Gauss_c4Fit->SetParameter(0,1); Gauss_c4Fit->SetParameter(1,3); Gauss_c4Fit->SetParameter(2,8); 
  Gauss_c4Fit->SetParName(0,"N");  Gauss_c4Fit->SetParName(1,"#lambda_{4}"); Gauss_c4Fit->SetParName(2,"R");
  c4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Fit(Gauss_c4Fit,"IME","",0.03,0.14);
  Gauss_c4Fit->Draw("same");*/

  // hight ED4 reference
  double y_ref[12]={0, 0, 1.00133, 0.980848, 0.988251, 0.994434, 0.999677, 1.00269, 1.00642, 1.00881, 1.01082, 1.01554};
  double y_ref_e[12]={0, 0, 0.054465, 0.00678447, 0.00194947, 0.000799564, 0.00039767, 0.000222628, 0.000135335, 8.75305e-05, 6.31392e-05, 5.53329e-05};
  double y_ref_Syst[12]={0, 0, 0.0167026, 0.00971243, 0.00945824, 0.0121959, 0.014654, 0.0155188, 0.0147369, 0.0128801, 0.0109854, 0.0100858};// without MC1 residue
  if(ChProdBOI==0){
    pad3_2->cd();
    gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.03); gPad->SetBottomMargin(0.32);
    
    TH1D *Ratio_4=(TH1D*)C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Clone();
    Ratio_4->GetXaxis()->SetTitleSize(2.2*SizeTitle);
    Ratio_4->GetXaxis()->SetLabelSize(2.2*SizeLabel);
    Ratio_4->GetYaxis()->SetTitleSize(2.2*SizeTitle);
    Ratio_4->GetYaxis()->SetLabelSize(2.2*SizeLabel);
    Ratio_4->GetXaxis()->SetTitleOffset(1.05);
    Ratio_4->GetYaxis()->SetTitleOffset(0.5);
    Ratio_4->GetXaxis()->SetNdivisions(606);
    Ratio_4->GetYaxis()->SetNdivisions(202);
    Ratio_4->GetYaxis()->SetTitle("#font[12]{C_{4}^{QS}} / built ");
    if(CollisionType==0){Ratio_4->SetMinimum(0.9); Ratio_4->SetMaximum(1.05);}
    else {Ratio_4->SetMinimum(0.8); Ratio_4->SetMaximum(1.18);}
    TH1D *Ratio_4_2=(TH1D*)Ratio_4->Clone();
    TH1D *Ratio_4_G=(TH1D*)Ratio_4->Clone();
    Ratio_4->SetMarkerStyle(24);
    Ratio_4_G->SetMarkerStyle(20);
    if(CollisionType==0 && FitBuild==0){
      Ratio_4->Divide(C4QSBuiltmerged[EDBin][MBOI][CollisionType]);
      Ratio_4_G->Divide(C4QSbuilt_G);
    }else{
      Ratio_4->Divide(C4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]);
      Ratio_4_2->Divide(C4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]);
    }
    TH1D *Ratio_4_Syst=(TH1D*)Ratio_4->Clone();
    Ratio_4_Syst->SetMarkerSize(0); Ratio_4_Syst->SetMarkerColor(kBlue-10); Ratio_4_Syst->SetFillColor(kBlue-10);
    Ratio_4_Syst->SetLineColor(kBlue-10);
    for(int bin=1; bin<=15; bin++){
      Ratio_4_Syst->SetBinError(bin, Syst_forRatio4[bin-1]);
    }
    
    Ratio_4_Syst->SetLineWidth(5);
    Ratio_4_Syst->Draw();
    if(CollisionType==0 && FitBuild==0){
      Ratio_4->Draw("same");
      Ratio_4_G->Draw("same");
      legend3_2->AddEntry(Ratio_4,"Built from #font[12]{C}_{2} (G=0%)","p");
      legend3_2->AddEntry(Ratio_4_G,BuiltNameC2->Data(),"p");
    }else{
     Ratio_4->Draw("same");
     Ratio_4_2->Draw("same");
     legend3_2->AddEntry(Ratio_4,BuiltNamec3->Data(),"p");
     legend3_2->AddEntry(Ratio_4_2,BuiltNameC3->Data(),"p");
     legend3_2->SetX1(0.52); legend3_2->SetY1(0.69); legend3_2->SetY2(0.95);
    }
      
    legend3_2->Draw("same");
  }

  /*if(CollisionType!=0){// FitBuilds from c3/C3
    TH1D *Ratio2_4 = (TH1D*)Ratio_4->Clone();
    TH1D *RatioM_4 = (TH1D*)Ratio_4->Clone();
    Ratio_4->Divide(C4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]);
    Ratio2_4->Divide(C4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]);
    RatioM_4->Divide(C4QSBuiltFromFitsmerged_M[EDBin][MBOI][CollisionType]);
    for(int bin=1; bin<=RatioM_4->GetNbinsX(); bin++){
      if(C4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->GetBinContent(bin) < 1) continue;
      RatioM_4->SetBinError(bin, (Ratio_4->GetBinContent(bin)-Ratio2_4->GetBinContent(bin))/2.);
      Ratio_4->SetBinError(bin, C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetBinError(bin) / C4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->GetBinContent(bin));
      Ratio2_4->SetBinError(bin, C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetBinError(bin) / C4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->GetBinContent(bin));
    }
    Ratio_4->SetMarkerStyle(24);
    RatioM_4->SetMarkerSize(0); RatioM_4->SetFillColor(kBlue-10);
    Ratio_4->SetMaximum(1.1);
    Ratio_4->Draw();
    RatioM_4->Draw("E2 same");
    Ratio2_4->Draw("same");
    Ratio_4->Draw("same");
    }*/
 
  
  
  Unity->Draw("same");
  
  
  if(ChProdBOI==0 && CollisionType==0){// chi2
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
    TLegend *legend4_2=(TLegend*)legend4->Clone();
    legend4_2->SetX1(0.7); legend4_2->SetX2(0.95); legend4_2->SetY1(0.6); legend4_2->SetY2(0.95); 

    
    TH1D *chi2_4 = new TH1D("chi2_4","",100,-0.5,99.5);
    chi2_4->SetLineColor(4);
    chi2_4->SetMarkerColor(4);
    chi2_4->GetXaxis()->SetTitle("Coherent fraction (%)"); chi2_4->GetYaxis()->SetTitle("#sqrt{#chi^{2}}");
    chi2_4->GetXaxis()->SetTitleSize(SizeTitle);  chi2_4->GetYaxis()->SetTitleSize(SizeTitle);
    chi2_4->GetXaxis()->SetLabelSize(SizeLabel);  chi2_4->GetYaxis()->SetLabelSize(SizeLabel);
    chi2_4->GetYaxis()->SetNdivisions(505);
    TH2D *chi2_2D_4[2];// Stat then Syst
    chi2_2D_4[0] = new TH2D("chi2Stat_2D_4","",7,0.5,7.5, 100,-0.5,99.5);
    chi2_2D_4[1] = new TH2D("chi2Syst_2D_4","",7,0.5,7.5, 100,-0.5,99.5);

    TH1D *tempDen = (TH1D*)C4QSBuiltmerged2D[EDBin][MBOI][CollisionType]->ProjectionY("TPFullWeight4_Den", 4, 4);
    TH1D *tempDenNeg = (TH1D*)C4QSNegBuiltmerged2D[EDBin][MBOI][CollisionType]->ProjectionY("TPNegFullWeight4_Den", 4, 4);
    tempDen->Add(tempDenNeg);// Add Pos and Neg Weight

    TH1D *C4QSBuilt_Gdependence1[5];
    TH1D *C4QSBuilt_Gdependence2[5];

    for(int binG=Gbinstart; binG<=Gbinstart+25; binG++){
      TString *proName=new TString("TPFullWeight4_");
      *proName += binG;
      TH1D *tempNum = (TH1D*)C4QSBuiltmerged2D[EDBin][MBOI][CollisionType]->ProjectionY(proName->Data(), binG, binG);
      proName->Append("_Neg");
      TH1D *tempNumNeg = (TH1D*)C4QSNegBuiltmerged2D[EDBin][MBOI][CollisionType]->ProjectionY(proName->Data(), binG, binG);
      //
      // Add Pos and Neg Weights
      tempNum->Add(tempNumNeg);
      //
      //tempNum->Add(tempDen);// now done in Plot_FourPion.C
      tempNum->Divide(tempDen);
      
      if(binG==5) C4QSBuilt_Gdependence1[0] = (TH1D*)tempNum->Clone();
      if(binG==15) C4QSBuilt_Gdependence1[1] = (TH1D*)tempNum->Clone();
      if(binG==25) C4QSBuilt_Gdependence1[2] = (TH1D*)tempNum->Clone();
      if(binG==35) C4QSBuilt_Gdependence1[3] = (TH1D*)tempNum->Clone();
      if(binG==45) C4QSBuilt_Gdependence1[4] = (TH1D*)tempNum->Clone();
      //
      if(binG==55) C4QSBuilt_Gdependence2[0] = (TH1D*)tempNum->Clone();
      if(binG==65) C4QSBuilt_Gdependence2[1] = (TH1D*)tempNum->Clone();
      if(binG==75) C4QSBuilt_Gdependence2[2] = (TH1D*)tempNum->Clone();
      if(binG==85) C4QSBuilt_Gdependence2[3] = (TH1D*)tempNum->Clone();
      if(binG==95) C4QSBuilt_Gdependence2[4] = (TH1D*)tempNum->Clone();

      if(tempNum->GetBinContent(Q4binChi2) <=0) continue;
      double value = fabs(C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetBinContent(Q4binChi2) - tempNum->GetBinContent(Q4binChi2));
      double err = pow(C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetBinError(Q4binChi2),2);
      err += pow(Syst_forChi2_4[Q4binChi2-1],2);
      err = sqrt(err);
      if(err<=0) continue;
      double Chi2 = pow(value / err,2);
      //
      
      
      chi2_4->SetBinContent(1 + 2*(binG-Gbinstart), sqrt(fabs(Chi2))); 
      chi2_4->SetBinError(1 + 2*(binG-Gbinstart), 0.001);
      
      //
      
      for(int binQ4=3; binQ4<=7; binQ4++){
	if(tempNum->GetBinContent(binQ4) <=0) continue;
	if(C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetBinError(binQ4) <= 0) continue;
	double value = fabs(C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetBinContent(binQ4) - tempNum->GetBinContent(binQ4));
	double err = pow(C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetBinError(binQ4),2);
	chi2_2D_4[0]->SetBinContent(binQ4, 1+(binG-Gbinstart), value/sqrt(err));
	//
	err = pow(Syst_forChi2_4[binQ4-1],2);
	chi2_2D_4[1]->SetBinContent(binQ4, 1+(binG-Gbinstart), value/sqrt(err));
      }
      
    }
    chi2_4->SetMarkerStyle(20);
    chi2_4->SetMinimum(0); chi2_4->SetMaximum(10);
    
    TString *Q4binName = new TString("0.0");
    
    if(int((Q4binChi2-1)*1.5*10)%10 == 0) *Q4binName += int((Q4binChi2-1)*1.5);
    else {*Q4binName += int((Q4binChi2-1)*1.5); *Q4binName += 5;}
    Q4binName->Append(" < #font[12]{Q_{4}} < 0.");
    if(Q4binChi2<7) Q4binName->Append("0");
    if(int((Q4binChi2)*1.5*10)%10 == 0) *Q4binName += int((Q4binChi2)*1.5);
    else {*Q4binName += int((Q4binChi2)*1.5); *Q4binName += 5;}
    Q4binName->Append(" GeV/#font[12]{c}");
    legend4->SetHeader(Q4binName->Data());
    chi2_4->Draw();
    legend4->AddEntry(chi2_4,"R_{coh}=0","p");
    legend4->Draw("same");

    TString *meanpTName = new TString("#LT #font[12]{p}_{T} #GT = 0.");
    *meanpTName += Q4_meanpT[EDBin][Q4binChi2-1];
    meanpTName->Append(" GeV/#font[12]{c}");
    TLatex *Specif_pT = new TLatex(0.17,0.9,meanpTName->Data());
    Specif_pT->SetNDC();
    Specif_pT->SetTextFont(TextFont);
    Specif_pT->SetTextSize(SizeSpecif);
    //Specif_pT->Draw("same");
   
    TLatex *ALICEspecif_5 = new TLatex(13,9,"ALICE Preliminary, Pb-Pb #\sqrt{#font[12]{s}_{NN}}=2.76 TeV");// ALICE specifications
    ALICEspecif_5->SetTextFont(TextFont);
    ALICEspecif_5->SetTextSize(SizeSpecif);
    //ALICEspecif_5->Draw("same");
    //
    TString *CentAndED_5 = new TString("Centrality 0-5%");
    if(EDBin==0) CentAndED_5->Append(", 0.16<#font[12]{K}_{T4}<0.3 GeV/#font[12]{c}");
    else CentAndED_5->Append(", 0.3<#font[12]{K}_{T4}<1.0 GeV/#font[12]{c}");
    TLatex *Cent_5 = new TLatex(13,8.2,CentAndED_5->Data());
    Cent_5->SetTextFont(TextFont);
    Cent_5->SetTextSize(SizeSpecif);
    Cent_5->Draw("same");

    TString *SaveNameChi2_4 = new TString("ChiSq_C4_bin");
    *SaveNameChi2_4 += Q4binChi2;
    SaveNameChi2_4->Append("_K");
    *SaveNameChi2_4 += EDBin;
    SaveNameChi2_4->Append("_M");
    *SaveNameChi2_4 += MBOI;
    SaveNameChi2_4->Append(".eps");
    //can4->SaveAs(SaveNameChi2_4->Data());
    
    //
    /*
    C4QSBuilt_Gdependence2[0]->SetLineColor(1); C4QSBuilt_Gdependence2[1]->SetLineColor(2);
    C4QSBuilt_Gdependence2[2]->SetLineColor(4); C4QSBuilt_Gdependence2[3]->SetLineColor(6);
    C4QSBuilt_Gdependence2[4]->SetLineColor(9);
    C4QSBuilt_Gdependence2[0]->GetXaxis()->SetRangeUser(0,.17);
    C4QSBuilt_Gdependence2[0]->GetXaxis()->SetTitleSize(SizeTitle);  C4QSBuilt_Gdependence2[0]->GetYaxis()->SetTitleSize(SizeTitle);
    C4QSBuilt_Gdependence2[0]->GetXaxis()->SetLabelSize(SizeLabel);  C4QSBuilt_Gdependence2[0]->GetYaxis()->SetLabelSize(SizeLabel);
    C4QSBuilt_Gdependence2[0]->GetXaxis()->SetNdivisions(505); C4QSBuilt_Gdependence2[0]->GetYaxis()->SetNdivisions(505);
    C4QSBuilt_Gdependence2[0]->GetXaxis()->SetTitle("#font[12]{Q_{4}} (GeV/#font[12]{c})");
    C4QSBuilt_Gdependence2[0]->GetYaxis()->SetTitle("#font[12]{C_{4}}^{QS}(built)");
    C4QSBuilt_Gdependence2[0]->Draw();
    C4QSBuilt_Gdependence2[1]->Draw("same"); C4QSBuilt_Gdependence2[2]->Draw("same");
    C4QSBuilt_Gdependence2[3]->Draw("same"); C4QSBuilt_Gdependence2[4]->Draw("same");
    legend4_2->SetHeader("R_{coh}=R_{ch}");
    legend4_2->AddEntry(C4QSBuilt_Gdependence2[0],"G=0%","l");
    legend4_2->AddEntry(C4QSBuilt_Gdependence2[1],"G=20%","l");
    legend4_2->AddEntry(C4QSBuilt_Gdependence2[2],"G=40%","l");
    legend4_2->AddEntry(C4QSBuilt_Gdependence2[3],"G=60%","l");
    legend4_2->AddEntry(C4QSBuilt_Gdependence2[4],"G=80%","l");
    legend4_2->Draw("same");
    Unity->Draw("same");
    */
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

    TH1D *GversusQ4[2];// Stat then Syst
    GversusQ4[0] = new TH1D("GversusQ4Stat","",7,0,0.105);
    GversusQ4[1] = new TH1D("GversusQ4Syst","",7,0,0.105);
    
    for(int ErrType=0; ErrType<2; ErrType++){
      for(int binQ4=3; binQ4<=7; binQ4++){
	double minG = 0;
	double minG_e1 = 0, minG_e2=0;
	double minChi=100;
	// Point Source
	for(int binG=1; binG<=25; binG++){// min
	  if(minChi > chi2_2D_4[ErrType]->GetBinContent(binQ4, binG)) {
	    minChi = chi2_2D_4[ErrType]->GetBinContent(binQ4, binG);
	    minG = 2*(binG-1);
	  }
	}
	//cout<<binQ4<<"  "<<minChi<<endl;
	for(int binG=1; binG<=25; binG++){// error
	  if(minG > 0) {
	    if(fabs(minChi - chi2_2D_4[ErrType]->GetBinContent(binQ4, binG)) < 1.) {
	      if(minG>2*(binG-1)) minG_e1 = fabs(minG - 2*(binG-1)); 
	      else minG_e2 = fabs(minG - 2*(binG-1));
	    }
	  }else{
	    if(fabs(minChi - chi2_2D_4[ErrType]->GetBinContent(binQ4, binG)) < 1.) {
	      minG_e1 = fabs(minG - 2*(binG-1)); 
	    }
	  }
	}
	GversusQ4[ErrType]->SetBinContent(binQ4, minG);
	if(minG_e1>minG_e2) GversusQ4[ErrType]->SetBinError(binQ4, minG_e1);
	else GversusQ4[ErrType]->SetBinError(binQ4, minG_e2);
	//
      }
    }// ErrType
    //
    
    GversusQ4[0]->SetMarkerStyle(20); GversusQ4[0]->SetMarkerColor(4); GversusQ4[0]->SetLineColor(4);
    GversusQ4[1]->SetMarkerSize(0); GversusQ4[1]->SetMarkerColor(kBlue-10); GversusQ4[1]->SetLineColor(kBlue);
    GversusQ4[1]->SetFillStyle(0);
        
    GversusQ4[0]->SetMinimum(0); GversusQ4[0]->SetMaximum(100); 
    GversusQ4[0]->GetXaxis()->SetTitle("#font[12]{Q_{4}} (GeV/#font[12]{c})"); GversusQ4[0]->GetYaxis()->SetTitle("Coherent fraction (%)");
    GversusQ4[0]->GetXaxis()->SetTitleSize(SizeTitle);  GversusQ4[0]->GetYaxis()->SetTitleSize(SizeTitle);
    GversusQ4[0]->GetXaxis()->SetLabelSize(SizeLabel);  GversusQ4[0]->GetYaxis()->SetLabelSize(SizeLabel);
    GversusQ4[0]->GetXaxis()->SetNdivisions(606); GversusQ4[0]->GetYaxis()->SetNdivisions(505);
    GversusQ4[0]->SetBinContent(1,200); GversusQ4[1]->SetBinContent(1,200);
    GversusQ4[0]->SetBinContent(2,200); GversusQ4[1]->SetBinContent(2,200);
    GversusQ4[0]->Draw();
    GversusQ4[1]->Draw("E2 same");
    GversusQ4[0]->Draw("same");
    //
    legend5->AddEntry(GversusQ4[0],"R_{coh}=0","p");
    legend5->Draw("same");
    //
    //
    TLatex *System_6 = new TLatex(0.03,82,System->Data());// ALICE specifications
    System_6->SetTextFont(TextFont);
    System_6->SetTextSize(SizeSpecif);
    System_6->Draw("same");
    //
    TLatex *Cent_6 = new TLatex(0.03,74,"Centrality 0-5%");
    Cent_6->SetTextFont(TextFont);
    Cent_6->SetTextSize(SizeSpecif);
    if(CollisionType==0) Cent_6->Draw("same");
    TString *EDname_6 = new TString("");
    if(EDBin==0) EDname_6->Append("0.16<#font[12]{K}_{T4}<0.3 GeV/#font[12]{c}");
    else EDname_6->Append("0.3<#font[12]{K}_{T4}<1.0 GeV/#font[12]{c}");
    TLatex *ED_6 = new TLatex(0.03,66,EDname_6->Data());
    ED_6->SetTextFont(TextFont);
    ED_6->SetTextSize(SizeSpecif);
    ED_6->Draw("same");
    //
    
    
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

