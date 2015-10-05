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

bool SaveFiles=0;
//double G_Cent[10]={30, 29, 27, 27, 29, 29, 29, 25, 29, 22};
//double G_stat_Cent[10]={2, 2, 2, 3, 3, 3, 5, 6, 6, 8};
//double G_syst_Cent[10]={7, 7, 12, 16, 7, 8, 10, 6, 7, 19};
//
//double G_Cent[10]={14, 21, 19, 21, 20, 28, 27, 28, 23, 25};
//double G_stat_Cent[10]={3, 3, 6, 4, 5, 5, 9, 8, 13, 16};
//double G_syst_Cent[10]={14, 8, 10, 18, 10, 6, 8, 5, 5, 18};
const int CollisionType=0;// Pb-Pb(0), p-Pb(1), pp(2)
const int ChProdBOI=0;// 0=SameCharge, 1=MixedCharge1, 2=MixedCharge2
const int EDBin=0;// KT3,4 bin. 0=low KT bin.  1=high KT bin
const int MBOI=0;// Centrality bin: 0-9
const int GValue = 30;// steps of 2
bool ShortSameCharge=1;// Plot short version?
bool FitBuild=0;
bool ReNormBuiltBaseline=1;
//
//
const int Q3binChi2= 3;// 2-5
const int Q4binChi2= 3;// 3-7
const int MBins=10;
const int EDBins=2;
const int RcohIndex=6;// 0(point source), 6(Full Size)
const int Gbinstart= 5 + RcohIndex*25;
const int GbinPlot=int( (GValue) /2. ) + 25*RcohIndex + 5;// +5 (Rcoh=0), +30 (Rcoh=1),....
//
//
int TextFont=42;// 63, or 42
float SizeLabel=0.06;// 20(63 font), 0.08(42 font)
float SizeLegend=0.05;// .08
float SizeTitle=0.06;// 
float SizeSpecif=0.045;// 
float SF1=2/3.*0.95;
float SF2=1/2.*0.95;
float LeftMargin=0.11;


double RightMargin=0.004;// 0.002
//

void Plot_plotsFourPion(){
  
  gStyle->SetOptStat(0);
  gStyle->SetOptDate(0);
  
  if(CollisionType!=0 || ChProdBOI!=0) ShortSameCharge=0;

  int Ysize=600;
  if(ChProdBOI==0) Ysize=800;
  double Ystart=0.3;
  if(ChProdBOI!=0) Ystart=0.005;
  double Xmarge=0.01;
  if(ChProdBOI!=0) Xmarge=0.14;

  TString *Centname=new TString("");
  int CentStart = MBOI*5;
  int CentEnd = (MBOI+1)*5;
  if(MBOI==2) CentEnd=20;
  if(MBOI==4) CentEnd=35;
  if(MBOI==7) CentEnd=50;
  *Centname += CentStart; Centname->Append("-"); *Centname += CentEnd; Centname->Append("%");

  TString *System=new TString("");
  if(CollisionType==0) {
    System->Append("ALICE ");
    System->Append(Centname->Data());
    System->Append(" Pb-Pb #\sqrt{#font[12]{s}_{NN}}=2.76 TeV");
  }else if(CollisionType==1) System->Append("ALICE p-Pb #\sqrt{#font[12]{s}_{NN}}=5.02 TeV");
  else System->Append("ALICE pp #\sqrt{#font[12]{s}}=7 TeV");
  double XstartALICE=0.43;
  if(CollisionType==1) XstartALICE=0.54;
  if(CollisionType==2) XstartALICE=0.65;
  TLatex *ALICEspecif = new TLatex(XstartALICE,.93,System->Data());// ALICE specifications
  ALICEspecif->SetNDC(1);
  ALICEspecif->SetTextFont(TextFont);
  ALICEspecif->SetTextSize(SizeSpecif);
  //
  TString *KT3=new TString("");
  if(EDBin==0) KT3->Append("0.16<#font[12]{K}_{T3}<0.3 GeV/#font[12]{c}");
  else KT3->Append("0.3<#font[12]{K}_{T3}<1.0 GeV/#font[12]{c}");
  TLatex *KT3specif = new TLatex(0.64,0.42,KT3->Data());
  KT3specif->SetNDC(1);
  KT3specif->SetTextFont(TextFont);
  KT3specif->SetTextSize(SizeSpecif);
  //
  TString *KT4=new TString("");
  if(EDBin==0) KT4->Append("0.16<#font[12]{K}_{T4}<0.3 GeV/#font[12]{c}");
  else KT4->Append("0.3<#font[12]{K}_{T4}<1.0 GeV/#font[12]{c}");
  TLatex *KT4specif = new TLatex(0.64,0.42,KT4->Data());
  KT4specif->SetNDC(1);
  KT4specif->SetTextFont(TextFont);
  KT4specif->SetTextSize(SizeSpecif);


  TString *BuiltNameE32 = new TString("#font[12]{E}_{3}(2) (G=");
  TString *BuiltNameE42 = new TString("#font[12]{E}_{4}(2) (G=");
  TString *BuiltNameE43 = new TString("#font[12]{E}_{4}(3) (G=");
  TString *BuiltNamee43 = new TString("#font[12]{e}_{4}(3) (G=");
  *BuiltNameE32 += GValue;
  *BuiltNameE42 += GValue;
  *BuiltNameE43 += GValue;
  *BuiltNamee43 += GValue;
  BuiltNameE32->Append("%)");
  BuiltNameE42->Append("%)");
  BuiltNameE43->Append("%)");
  BuiltNamee43->Append("%)");
  
  TString *ChargeCombSt4 = new TString("");
  if(ChProdBOI==0) ChargeCombSt4->Append("#pi^{-}#pi^{-}#pi^{-}#pi^{-}");
  else if(ChProdBOI==1) ChargeCombSt4->Append("#pi^{-}#pi^{-}#pi^{-}#pi^{+}");
  else ChargeCombSt4->Append("#pi^{-}#pi^{-}#pi^{+}#pi^{+}");
  float xCCT=0.16, yCCT=0.89;
  if(ChProdBOI==0) {xCCT=0.7; yCCT=0.27;}
  TLatex *ChargeCombTitle4 = new TLatex(xCCT,yCCT,ChargeCombSt4->Data());// 0.7,yCCT
  ChargeCombTitle4->SetNDC(1);
  ChargeCombTitle4->SetTextFont(TextFont);
  ChargeCombTitle4->SetTextSize(2*SizeSpecif);


  TF1 *C3ratioFit = new TF1("C3ratioFit","[0] + [1]*exp(-pow(x*[2],2))",0,1);
  C3ratioFit->FixParameter(0,9.98651e-01); C3ratioFit->FixParameter(1,2.48247e-02); C3ratioFit->FixParameter(2,3.13184e+01);
  TF1 *C4ratioFit = new TF1("C4ratioFit","[0] + [1]*exp(-pow(x*[2],2))",0,1);
  C4ratioFit->FixParameter(0,9.98432e-01); C4ratioFit->FixParameter(1,9.43508e-02); C4ratioFit->FixParameter(2,2.89368e+01);


  int Q3_meanpT[2][5]={{0,213,228,233,237},{0,316,336,354,366}};// low then high KT3
  int Q4_meanpT[2][7]={{0,0,221,229,234,238,241},{0,0,325,335,344,351,355}};// low then high KT4

  double ReNormL_3, ReNormL_3_vary;
  double ReNormH_3, ReNormH_3_vary;
  double ReNormL_4, ReNormL_4_vary;
  double ReNormH_4, ReNormH_4_vary;
  if(CollisionType==0){
    ReNormL_3=0.095 + 0.03*(MBOI/9.);// 0.095
    ReNormH_3=0.105 + 0.03*(MBOI/9.);// 0.105
    ReNormL_4=0.125 + 0.04*(MBOI/9.);// 0.125
    ReNormH_4=0.145 + 0.04*(MBOI/9.);// 0.145
    //
    ReNormL_3_vary=ReNormL_3 - 0.02;// 0.075
    ReNormH_3_vary=ReNormH_3 - 0.02;// 0.085
    ReNormL_4_vary=ReNormL_4 + 0.015;// 0.12
    ReNormH_4_vary=ReNormH_4 + 0.015;// 0.13
  }else{
    ReNormL_3=0.3;
    ReNormH_3=0.32;
    ReNormL_3_vary=0.28;
    ReNormH_3_vary=0.3;
    //
    ReNormL_4=0.46;// 0.46
    ReNormH_4=0.49;// 0.49
    //
    ReNormL_4_vary=0.4;// 0.4
    ReNormH_4_vary=0.44;// 0.44
  }
  double C3Builtsyst_ReNorm=1;
  double C4Builtsyst_ReNorm=1;
  double C4BuiltFromFitsyst_ReNorm=1;

  int binStartFinalG=3 ,binEndFinalG=7;
  int BadBinsSC[3] = {0};

  float Q4max=0.16;//Pb-Pb
  if(CollisionType==1) Q4max=0.5;
  if(CollisionType==2) Q4max=0.5;
  
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
	      C3QS[ChComb][ch][ED][mb][CT]->GetYaxis()->SetTitle("Three pion correlation");
	      c3QS[ChComb][ch][ED][mb][CT]->GetYaxis()->SetTitle("#font[12]{#bf{c}}_{3}");
	    }
	    C4QS[ChComb][ch][ED][mb][CT]=(TH1D*)files[ChComb][ch][ED][mb][CT]->Get("C4QS");
	    c4QS[ChComb][ch][ED][mb][CT]=(TH1D*)files[ChComb][ch][ED][mb][CT]->Get("c4QS");
	    a4QS[ChComb][ch][ED][mb][CT]=(TH1D*)files[ChComb][ch][ED][mb][CT]->Get("a4QS");
	    C4QS[ChComb][ch][ED][mb][CT]->SetDirectory(0); c4QS[ChComb][ch][ED][mb][CT]->SetDirectory(0); 
	    a4QS[ChComb][ch][ED][mb][CT]->SetDirectory(0);
	    C4QS[ChComb][ch][ED][mb][CT]->GetXaxis()->SetTitle("#font[12]{Q}_{4} (GeV/#font[12]{c})");
	    c4QS[ChComb][ch][ED][mb][CT]->GetXaxis()->SetTitle("#font[12]{Q}_{4} (GeV/#font[12]{c})");
	    C4QS[ChComb][ch][ED][mb][CT]->GetYaxis()->SetTitle("Four pion correlation");
	    c4QS[ChComb][ch][ED][mb][CT]->GetYaxis()->SetTitle("#font[12]{#bf{c}}_{4}");
	    
	    if(ChComb==0){
	      b4QS[ch][ED][mb][CT]=(TH1D*)files[ChComb][ch][ED][mb][CT]->Get("b4QS");
	      b4QS[ch][ED][mb][CT]->SetDirectory(0);
	      //if(CT==0){
	      //cout<<ch<<"  "<<ED<<" "<<mb<<" "<<CT<<endl;
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
		if(CT==0){
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
  cout<<"Done getting files"<<endl;
  
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
  TH1D *C4QSbuilt_G[10];// Mbin
  TH1D *a4QSbuilt_G[10];// Mbin
  TH1D *b4QSbuilt_G[10];// Mbin
  TH1D *c4QSbuilt_G[10];// Mbin
  TH1D *C4QSNegbuilt_G[10];// Mbin
  TH1D *a4QSNegbuilt_G[10];// Mbin
  TH1D *b4QSNegbuilt_G[10];// Mbin
  TH1D *c4QSNegbuilt_G[10];// Mbin
  

  for(int CT=0; CT<3; CT++){
    for(int mb=0; mb<MBins; mb++){
      if(CT>0 && mb>0) continue;
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
	    //if(CT==0){
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
	      //}
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
	      if(ChComb==0){
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
		//if(CT==0){
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
		  //}
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
	  // additional bad bins
	  int BadBins;
	  if(CT==0){
	    if(ChComb==0){
	      if(ED==0) BadBins=2;
	      else BadBins=2;
	      BadBinsSC[CT] = BadBins+1;
	    }else if(ChComb==1){
	      if(ED==0) BadBins=1;
	      else BadBins=2;
	    }else{
	      if(ED==0) BadBins=1;
	      else BadBins=2;
	    }
	  }else if(CT==1){
	    if(ChComb==0){
	      if(ED==0) BadBins=3;
	      else BadBins=3;
	      BadBinsSC[CT] = BadBins+1;
	    }else if(ChComb==1){
	      if(ED==0) BadBins=2;
	      else BadBins=3;
	    }else{
	      if(ED==0) BadBins=3;
	      else BadBins=3;
	    }
	  }else{
	    if(ChComb==0){
	      if(ED==0) BadBins=4;
	      else BadBins=5;
	      BadBinsSC[CT] = BadBins;
	    }else if(ChComb==1){
	      if(ED==0) BadBins=3;
	      else BadBins=4;
	    }else{
	      if(ED==0) BadBins=4;
	      else BadBins=4;
	    }
	  }

	  for(int BB=1; BB<=BadBins; BB++){
	    C4QSmerged[ChComb][ED][mb][CT]->SetBinContent(BB,100.);
	    a4QSmerged[ChComb][ED][mb][CT]->SetBinContent(BB,100.);
	    b4QSmerged[ED][mb][CT]->SetBinContent(BB,100.);
	    c4QSmerged[ChComb][ED][mb][CT]->SetBinContent(BB,100.);
	    if(ChComb==0){
	    C4QSBuiltmerged[ED][mb][CT]->SetBinContent(BB,0.1);
	    a4QSBuiltmerged[ED][mb][CT]->SetBinContent(BB,0.1);
	    b4QSBuiltmerged[ED][mb][CT]->SetBinContent(BB,0.1);
	    c4QSBuiltmerged[ED][mb][CT]->SetBinContent(BB,0.1);
	    C4QSBuiltFromFitsmerged1[ED][mb][CT]->SetBinContent(BB,0);
	    C4QSBuiltFromFitsmerged2[ED][mb][CT]->SetBinContent(BB,0);
	    a4QSBuiltFromFitsmerged1[ED][mb][CT]->SetBinContent(BB,0);
	    a4QSBuiltFromFitsmerged2[ED][mb][CT]->SetBinContent(BB,0);
	    b4QSBuiltFromFitsmerged1[ED][mb][CT]->SetBinContent(BB,0);
	    b4QSBuiltFromFitsmerged2[ED][mb][CT]->SetBinContent(BB,0);
	    c4QSBuiltFromFitsmerged1[ED][mb][CT]->SetBinContent(BB,0);
	    c4QSBuiltFromFitsmerged2[ED][mb][CT]->SetBinContent(BB,0);
	    }
	  }
	  
	}// ED
      }// ChComb
    }// mb
  }// CT

  
  // make separate G versions for each M bin
  for(int mb=0; mb<MBins; mb++){
    if(CollisionType!=0 && mb>0) continue;
    TString *proName4=new TString("C4QSbuilt_G_M"); TString *proNameNeg4=new TString("C4QSNegbuilt_G_M");
    *proName4 += mb; *proNameNeg4 += mb;
    TString *proName_a4=new TString("a4QSbuilt_G_M"); TString *proNameNeg_a4=new TString("a4QSNegbuilt_G_M");
    *proName_a4 += mb; *proNameNeg_a4 += mb;
    TString *proName_b4=new TString("b4QSbuilt_G_M"); TString *proNameNeg_b4=new TString("b4QSNegbuilt_G_M");
    *proName_b4 += mb; *proNameNeg_b4 += mb;
    TString *proName_c4=new TString("c4QSbuilt_G_M"); TString *proNameNeg_c4=new TString("c4QSNegbuilt_G_M");
    *proName_c4 += mb; *proNameNeg_c4 += mb;
    C4QSbuilt_G[mb] = (TH1D*)C4QSBuiltmerged2D[EDBin][mb][CollisionType]->ProjectionY(proName4->Data(), GbinPlot, GbinPlot);
    C4QSNegbuilt_G[mb] = (TH1D*)C4QSNegBuiltmerged2D[EDBin][mb][CollisionType]->ProjectionY(proNameNeg4->Data(), GbinPlot, GbinPlot);
    a4QSbuilt_G[mb] = (TH1D*)a4QSBuiltmerged2D[EDBin][mb][CollisionType]->ProjectionY(proName_a4->Data(), GbinPlot, GbinPlot);
    a4QSNegbuilt_G[mb] = (TH1D*)a4QSNegBuiltmerged2D[EDBin][mb][CollisionType]->ProjectionY(proNameNeg_a4->Data(), GbinPlot, GbinPlot);
    b4QSbuilt_G[mb] = (TH1D*)b4QSBuiltmerged2D[EDBin][mb][CollisionType]->ProjectionY(proName_b4->Data(), GbinPlot, GbinPlot);
    b4QSNegbuilt_G[mb] = (TH1D*)b4QSNegBuiltmerged2D[EDBin][mb][CollisionType]->ProjectionY(proNameNeg_b4->Data(), GbinPlot, GbinPlot);
    c4QSbuilt_G[mb] = (TH1D*)c4QSBuiltmerged2D[EDBin][mb][CollisionType]->ProjectionY(proName_c4->Data(), GbinPlot, GbinPlot);
    c4QSNegbuilt_G[mb] = (TH1D*)c4QSNegBuiltmerged2D[EDBin][mb][CollisionType]->ProjectionY(proNameNeg_c4->Data(), GbinPlot, GbinPlot);
    proName4->Append("_FullWeightDen"); proNameNeg4->Append("_FullWeightDen");
    TH1D *tempDen4 = (TH1D*)C4QSBuiltmerged2D[EDBin][mb][CollisionType]->ProjectionY(proName4->Data(), 4, 4);
    TH1D *tempDenNeg4 = (TH1D*)C4QSNegBuiltmerged2D[EDBin][mb][CollisionType]->ProjectionY(proNameNeg4->Data(), 4, 4);
    // Add Pos with Neg weights
    tempDen4->Add(tempDenNeg4);
    C4QSbuilt_G[mb]->Add(C4QSNegbuilt_G[mb]);
    a4QSbuilt_G[mb]->Add(a4QSNegbuilt_G[mb]);
    b4QSbuilt_G[mb]->Add(b4QSNegbuilt_G[mb]);
    c4QSbuilt_G[mb]->Add(c4QSNegbuilt_G[mb]);
    C4QSbuilt_G[mb]->Divide(tempDen4);
    a4QSbuilt_G[mb]->Divide(tempDen4);
    b4QSbuilt_G[mb]->Divide(tempDen4);
    c4QSbuilt_G[mb]->Divide(tempDen4);
    C4QSbuilt_G[mb]->SetLineColor(4);
    a4QSbuilt_G[mb]->SetLineColor(2);
    b4QSbuilt_G[mb]->SetLineColor(6);
    c4QSbuilt_G[mb]->SetLineColor(1);
  }

  
  // merge centralities
  if(CollisionType==0){
  for(int mb=0; mb<MBins; mb++){
    for(int bin=1; bin<=20; bin++){
      if(mb==2){
	double value=C4QSmerged[0][EDBin][mb][0]->GetBinContent(bin) + C4QSmerged[0][EDBin][mb+1][0]->GetBinContent(bin);
	value /= 2.;
	double value_e= sqrt(pow(C4QSmerged[0][EDBin][mb][0]->GetBinError(bin),2) + pow(C4QSmerged[0][EDBin][mb+1][0]->GetBinError(bin),2));
	value_e = sqrt(value_e)/2.;
	C4QSmerged[0][EDBin][mb][0]->SetBinContent(bin,value); C4QSmerged[0][EDBin][mb][0]->SetBinError(bin,value_e); 
	//
	value = C4QSBuiltmerged[EDBin][mb][0]->GetBinContent(bin) + C4QSBuiltmerged[EDBin][mb+1][0]->GetBinContent(bin);
	value /= 2.;
	C4QSBuiltmerged[EDBin][mb][0]->SetBinContent(bin,value); C4QSBuiltmerged[EDBin][mb][0]->SetBinError(bin,0); 
	//
	value = C4QSbuilt_G[mb]->GetBinContent(bin) + C4QSbuilt_G[mb+1]->GetBinContent(bin);
	value /= 2.;
	C4QSbuilt_G[mb]->SetBinContent(bin, value); C4QSbuilt_G[mb]->SetBinError(bin,0); 
      }
      if(mb==4 || mb==7){
	double value=C4QSmerged[0][EDBin][mb][0]->GetBinContent(bin);
	value += C4QSmerged[0][EDBin][mb+1][0]->GetBinContent(bin);
	value += C4QSmerged[0][EDBin][mb+2][0]->GetBinContent(bin);
	value /= 3.;
	double value_e= pow(C4QSmerged[0][EDBin][mb][0]->GetBinError(bin),2);
	value_e += pow(C4QSmerged[0][EDBin][mb+1][0]->GetBinError(bin),2);
	value_e += pow(C4QSmerged[0][EDBin][mb+2][0]->GetBinError(bin),2);
	value_e = sqrt(value_e)/3.;
	C4QSmerged[0][EDBin][mb][0]->SetBinContent(bin,value); C4QSmerged[0][EDBin][mb][0]->SetBinError(bin,value_e); 
	//
	value = C4QSBuiltmerged[EDBin][mb][0]->GetBinContent(bin);
	value += C4QSBuiltmerged[EDBin][mb+1][0]->GetBinContent(bin);
	value += C4QSBuiltmerged[EDBin][mb+2][0]->GetBinContent(bin);
	value /= 3.;
	C4QSBuiltmerged[EDBin][mb][0]->SetBinContent(bin,value); C4QSBuiltmerged[EDBin][mb][0]->SetBinError(bin,0); 
	//
	value = C4QSbuilt_G[mb]->GetBinContent(bin);
	value += C4QSbuilt_G[mb+1]->GetBinContent(bin);
	value += C4QSbuilt_G[mb+2]->GetBinContent(bin);
	value /= 3.;
	C4QSbuilt_G[mb]->SetBinContent(bin,value); C4QSbuilt_G[mb]->SetBinError(bin,0); 
      }
    }
  }
  
  
 
  


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
    if(ChProdBOI==0) pad1_2->Draw();
    pad1->cd();
    gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.03); 
    gPad->SetBottomMargin(Xmarge);
    

    TLegend *legend1 = new TLegend(.65,.5, .85,.85,NULL,"brNDC");//.45 or .4 for x1
    legend1->SetBorderSize(0);
    legend1->SetFillColor(0);
    legend1->SetTextFont(TextFont);
    legend1->SetTextSize(SizeLegend);
    TLegend *legend1_2=(TLegend*)legend1->Clone();
    legend1_2->SetX1(0.6); legend1_2->SetX2(0.93); legend1_2->SetY1(0.34); legend1_2->SetY2(0.6); 
    legend1_2->SetTextSize(2.*SizeLegend);
    
    C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->SetTitleSize(SizeTitle);
    C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->SetLabelSize(SizeLabel);
    C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetYaxis()->SetTitleSize(SizeTitle);
    C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetYaxis()->SetLabelSize(SizeLabel);
    C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->SetTitleOffset(1.05);
    C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetYaxis()->SetTitleOffset(0.9);
    C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->SetNdivisions(606);
    C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetYaxis()->SetNdivisions(505);
    if(CollisionType!=0) C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->SetRangeUser(0,0.35);
    //
    
       
    //TH1D *C3QS_Syst = (TH1D*)C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Clone();
    TH1D *C3QS_Syst = new TH1D("C3QS_Syst","",C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetNbinsX(),C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->GetBinLowEdge(1), C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->GetBinUpEdge(C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetNbinsX()));// was 0.0015 offset
    
    TH1D *C3QSBuilt_Syst = (TH1D*)C3QSBuiltmerged[EDBin][MBOI][CollisionType]->Clone();
    
    if(ChProdBOI==0){ 
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
      if(CollisionType!=0) C3QSbuilt_G->GetXaxis()->SetRange(2,35);
      if(ReNormBuiltBaseline){
	int BinL = C3QSbuilt_G->GetXaxis()->FindBin(ReNormL_3);
	int BinH = C3QSbuilt_G->GetXaxis()->FindBin(ReNormH_3);
	C3Builtsyst_ReNorm = 1 / (C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Integral(BinL,BinH) / C3QSbuilt_G->Integral(BinL,BinH));
	C3QSbuilt_G->Scale( C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Integral(BinL,BinH) / C3QSbuilt_G->Integral(BinL,BinH));
      }
      
    }



    for(int bin=1; bin<=C3QS_Syst->GetNbinsX(); bin++){
      double q3 = C3QS_Syst->GetXaxis()->GetBinCenter(bin);
      if(CollisionType!=0) q3 /= 3.1;// rescale q3 for pp and p-Pb
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
      //if(CollisionType==0){
	C3QSBuilt_Syst->SetBinContent(bin, 4.7);
	double syst2 = pow(0.002 - 0.002*q3/0.1,2);// 11h to 10h
	syst2 += pow(0.9856 + 0.3285*q3 - 1.897*q3*q3 - 1,2);// MRC
	syst2 += pow(0.9786 + 0.421*q3 - 2.108*q3*q3 - 1,2);// Muon, 92%
	syst2 += pow(0.946 + 0.849*q3 - 3.316*q3*q3 - 1,2);// fc2 scale
	//syst2 += pow((0.024-0.276*q3)/2.,2);// old way of cubic vs linear Interpolator
	syst2 += pow(0.5*(C3ratioFit->Eval(q3)-1),2);// new Interpolator
	if(EDBin==1) syst2 += pow( 0.9973 + 0.462*q3 - 14.2*pow(q3,2) - 2.135*pow(q3,3) + 805.3*pow(q3,4) + 1358*pow(q3,5) - 1,2);// Tij sign
	syst2 += pow(C3Builtsyst_ReNorm-1,2);// Renormalization 
	syst2 = sqrt(syst2);
	syst2 *= C3QSBuiltmerged[EDBin][MBOI][CollisionType]->GetBinContent(bin);
	C3QSBuilt_Syst->SetBinError(bin, syst2);
	// }
    }
    double Syst_forChi2_3[15]={0};
    double Syst_forRatio[15]={0};
    for(int bin=1; bin<=15; bin++){
      double q3 = C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->GetBinCenter(bin);
      if(CollisionType!=0) q3 /= 3.1;// rescale q3 for pp and p-Pb
      int SystBin = C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->FindBin(q3);
      double SystPercent_Diff = pow(0.01*(1 - q3/0.1),2);// f coefficients
      //SystPercent_Diff += pow(0.0103*exp(-41.68*q3),2);// Old way of cubic vs linear Interpolator
      SystPercent_Diff += pow(0.5*(C3ratioFit->Eval(q3)-1),2);// new Interpolator
      if(EDBin==1) SystPercent_Diff += pow(0.9973 + 0.462*q3 - 14.2*pow(q3,2) - 2.135*pow(q3,3) + 805.3*pow(q3,4) + 1358*pow(q3,5) - 1,2);// Tij sign
      SystPercent_Diff += pow( (0.9847 + 0.358*q3 - 2.133*q3*q3) - (0.9856 + 0.3285*q3 - 1.897*q3*q3),2);// MRC
      SystPercent_Diff += pow( (0.975 + 0.4189*q3 - 2.055*q3*q3) - (0.9786 + 0.421*q3 - 2.108*q3*q3),2);// Muon
      SystPercent_Diff += pow( (0.936 + 1.194*q3 - 5.912*q3*q3) - (0.946 + 0.849*q3 - 3.316*q3*q3),2);// fc2 scale
      SystPercent_Diff += pow( (C3Builtsyst_ReNorm-1),2);// Renormalization
      //SystPercent_Diff += pow( 0.0858 - 2.808*q3 + 19.72*q3*q3,2);// MC cumulant residue
      SystPercent_Diff = sqrt(SystPercent_Diff);
      Syst_forChi2_3[bin-1] = SystPercent_Diff * C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetBinContent(bin);
      //if(CollisionType==0){
	// ratio
	double A = C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetBinContent(bin);
	double B = C3QSBuiltmerged[EDBin][MBOI][CollisionType]->GetBinContent(bin);
	if(B<1) continue;
	double ratio = A / B;
	Syst_forRatio[bin-1] = pow( ratio * (0.01*(1 - q3/0.1)) ,2);// f coefficients
	//Syst_forRatio[bin-1] += pow( ratio * 0.0103*exp(-41.68*q3) ,2);// old Interpolator
	Syst_forRatio[bin-1] += pow( ratio * 0.5*(C3ratioFit->Eval(q3)-1),2);// new Interpolator
	if(EDBin==1) Syst_forRatio[bin-1] += pow( ratio * (0.9973 + 0.462*q3 - 14.2*pow(q3,2) - 2.135*pow(q3,3) + 805.3*pow(q3,4) + 1358*pow(q3,5) - 1) ,2);// Tij sign
	Syst_forRatio[bin-1] += pow( ratio * sqrt(pow(0.9847 + 0.358*q3 - 2.133*q3*q3-1,2) + pow(0.9856 + 0.3285*q3 - 1.897*q3*q3-1,2) - 2*(0.9847 + 0.358*q3 - 2.133*q3*q3-1)*(0.9856 + 0.3285*q3 - 1.897*q3*q3-1) ) ,2);// MRC
	Syst_forRatio[bin-1] += pow( ratio * sqrt(pow(0.975 + 0.4189*q3 - 2.055*q3*q3-1,2) + pow(0.9786 + 0.421*q3 - 2.108*q3*q3-1,2) - 2*(0.975 + 0.4189*q3 - 2.055*q3*q3-1)*(0.9786 + 0.421*q3 - 2.108*q3*q3-1) ) ,2);// Muon
	Syst_forRatio[bin-1] += pow( ratio * sqrt(pow(0.936 + 1.194*q3 - 5.912*q3*q3-1,2) + pow(0.946 + 0.849*q3 - 3.316*q3*q3-1,2) - 2*(0.936 + 1.194*q3 - 5.912*q3*q3-1)*(0.946 + 0.849*q3 - 3.316*q3*q3-1) ) ,2);// fc2 scale
	Syst_forRatio[bin-1] += pow( ratio * (C3Builtsyst_ReNorm-1) ,2);// Renormalization 
	Syst_forRatio[bin-1] = sqrt(Syst_forRatio[bin-1]);
	// }
    }
    C3QS_Syst->SetBinContent(1,100); 
    C3QS_Syst->SetMarkerSize(0); C3QS_Syst->SetFillColor(kBlue-10);
    C3QS_Syst->SetLineColor(kBlue-10); C3QS_Syst->SetMarkerColor(kBlue-10);
    C3QS_Syst->SetLineWidth(57);
    //if(CollisionType==0) {
      C3QSBuilt_Syst->SetBinContent(1,100); 
      C3QSBuilt_Syst->SetMarkerSize(0); C3QSBuilt_Syst->SetFillColor(kGray); //C3QSBuilt_Syst->SetFillStyle(3004);
      C3QSBuilt_Syst->SetLineColor(kGray); C3QSBuilt_Syst->SetLineWidth(5);
      C3QSBuilt_Syst->SetMarkerColor(kGray); C3QSBuilt_Syst->GetXaxis()->SetRangeUser(0.01,0.2);
      C3QSBuiltmerged[EDBin][MBOI][CollisionType]->SetLineWidth(1.2);
      C3QSBuiltmerged[EDBin][MBOI][CollisionType]->SetLineStyle(2);
      c3QSBuiltmerged[EDBin][MBOI][CollisionType]->SetLineStyle(2);
      //}
    C3QS_Syst->GetXaxis()->SetRangeUser(0.01,0.2);
    //
    C3QSBuiltmerged[EDBin][MBOI][CollisionType]->GetXaxis()->SetRange(2,15);
    if(CollisionType!=0) C3QSBuiltmerged[EDBin][MBOI][CollisionType]->GetXaxis()->SetRange(2,35);
    C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->SetMaximum(5.);
    C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->SetBarWidth(0.0001);
    C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Draw();
    
    C3QS_Syst->Draw("same");
    TF1 *SystMeanLine3 = new TF1("SystMeanLine3","[0]",0.01,1); SystMeanLine3->SetParameter(0, C3QS_Syst->GetBinContent(5)); SystMeanLine3->SetLineColor(1);
    SystMeanLine3->Draw("same");
    //if(CollisionType==0) C3QSBuilt_Syst->Draw("same");// E1
    

    C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Draw("same");
    if(!ShortSameCharge) c3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Draw("same");
    if(ChProdBOI==0){ 
      C3QSBuiltmerged[EDBin][MBOI][CollisionType]->Draw("same");
      if(!ShortSameCharge) c3QSBuiltmerged[EDBin][MBOI][CollisionType]->Draw("same");
      if(ShortSameCharge) C3QSbuilt_G->Draw("same");
      //c3QSbuilt_G->Draw("same");
    }

    legend1->AddEntry(C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType],"#font[12]{C}_{3}^{QS}","p");
    if(!ShortSameCharge) legend1->AddEntry(c3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType],"#font[12]{#bf{c}}_{3}^{QS}","p");
    if(ChProdBOI==0) {
      legend1->AddEntry(C3QSBuiltmerged[EDBin][MBOI][CollisionType],"#font[12]{E}_{3}(2) (G=0%)","l");
      if(ShortSameCharge) legend1->AddEntry(C3QSbuilt_G,BuiltNameE32->Data(),"l");
      //legend1->AddEntry(C3QSBuilt_Syst,"#font[12]{E}_{3} systematic","l");
    }  
    //legend1->AddEntry(C3QS_Syst,"#font[12]{C}_{3}^{QS} systematic","l");
    legend1->Draw("same");
    
    Unity->Draw("same");
    
    ALICEspecif->Draw("same");
    KT3specif->Draw("same");
   
    //for(int i=1; i<=10; i++) cout<<c3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetBinContent(i)<<", ";
    //cout<<endl;
    
    if(ChProdBOI==0){
      pad1_2->cd();
      gPad->SetLeftMargin(LeftMargin); gPad->SetRightMargin(0.04);
      gPad->SetTopMargin(0.03); gPad->SetBottomMargin(0.32);

      TH1D *Ratio_3 = (TH1D*)C3QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Clone();
      Ratio_3->GetXaxis()->SetTitleSize(2.2*SizeTitle);
      Ratio_3->GetXaxis()->SetLabelSize(2.2*SizeLabel);
      Ratio_3->GetYaxis()->SetTitleSize(2.2*SizeTitle);
      Ratio_3->GetYaxis()->SetLabelSize(2.2*SizeLabel);
      Ratio_3->GetXaxis()->SetTitleOffset(1.05);
      Ratio_3->GetYaxis()->SetTitleOffset(0.4);
      Ratio_3->GetXaxis()->SetNdivisions(606);
      Ratio_3->GetYaxis()->SetNdivisions(204);
      Ratio_3->GetYaxis()->SetTitle("Ratio");
      Ratio_3->SetMinimum(0.93); Ratio_3->SetMaximum(1.025);
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

      if(ShortSameCharge) Ratio_3_G->Draw("same");

      Unity->Draw("same");
      
      legend1_2->AddEntry(Ratio_3,"#font[12]{C}_{3}^{QS}/#font[12]{E}_{3}(2) (G=0%)","p");
      if(ShortSameCharge) {
	TString *name = new TString("#font[12]{C}_{3}^{QS}/"); name->Append(BuiltNameE32->Data());
	legend1_2->AddEntry(Ratio_3_G,name->Data(),"p");
      }
      legend1_2->Draw("same");
      
      if(SaveFiles && ChProdBOI==0 && !FitBuild && MBOI==0 && CollisionType==0){
	TString *SaveName = new TString("PlotsForPapers/Fig_C3SC_PbPb_");
	//
	if(ShortSameCharge) SaveName->Append("short_");
	else SaveName->Append("long_");
	
	SaveName->Append("K");
	*SaveName += EDBin;
	SaveName->Append("_M");
	*SaveName += MBOI;
	SaveName->Append(".eps");
	//
	can1->SaveAs(SaveName->Data());
      }

      /*for(int i=1; i<=10; i++) cout<<Ratio_3->GetBinContent(i)<<", ";
      cout<<endl;
      for(int i=1; i<=10; i++) cout<<Ratio_3->GetBinError(i)<<", ";
      cout<<endl;
      for(int i=1; i<=10; i++) cout<<Ratio_3_Syst->GetBinError(i)<<", ";
      cout<<endl;*/

      if(CollisionType==0){
      ///////////////////////////////////////////////////
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
      gPad->SetLeftMargin(LeftMargin); gPad->SetRightMargin(0.04);
      gPad->SetTopMargin(0.03); gPad->SetBottomMargin(0.14);
      TLegend *legend2 = new TLegend(.55,.75, .95,.95,NULL,"brNDC");//.45 or .4 for x1
      legend2->SetBorderSize(0);
      legend2->SetFillColor(0);
      legend2->SetTextFont(TextFont);
      legend2->SetTextSize(SizeLegend);
      
      
      
      TH1D *chi2_3 = new TH1D("chi2_3","",100,-0.5,99.5);
      
      chi2_3->SetLineColor(2);
      chi2_3->SetMarkerColor(2);
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
      legend2->AddEntry(chi2_3,"R_{coh}=R_{ch}","p");
      legend2->Draw("same");
      
      TString *meanpTName3 = new TString("#LT #font[12]{p}_{T} #GT = 0.");
      *meanpTName3 += Q3_meanpT[EDBin][Q3binChi2-1];
      meanpTName3->Append(" GeV/#font[12]{c}");
      TLatex *Specif_pT3 = new TLatex(0.17,0.9,meanpTName3->Data());
      Specif_pT3->SetNDC();
      Specif_pT3->SetTextFont(TextFont);
      Specif_pT3->SetTextSize(SizeSpecif);
      //Specif_pT3->Draw("same");
      
      TString *CentAndED_chi3 = new TString(Centname->Data());
      if(EDBin==0) CentAndED_chi3->Append(", 0.16<#font[12]{K}_{T3}<0.3 GeV/#font[12]{c}");
      else CentAndED_chi3->Append(", 0.3<#font[12]{K}_{T3}<1.0 GeV/#font[12]{c}");
      TLatex *Cent_chi3 = new TLatex(45,3,CentAndED_chi3->Data());
      Cent_chi3->SetTextFont(TextFont);
      Cent_chi3->SetTextSize(SizeSpecif);
      Cent_chi3->Draw("same");
      
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
      gPad->SetLeftMargin(LeftMargin); gPad->SetRightMargin(0.05);
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
      TLatex *Cent_2 = new TLatex(0.022,77,Centname->Data());
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
      
      
      }// 
    }// ChProdBOI!=0
  }// ChProdBOI!=2
  if(CollisionType==0 && !FitBuild && MBOI==0 && ChProdBOI==0){
      TCanvas *can2_3 = new TCanvas("can2_3", "can2_3",100,700,700,700);// 11,53,700,500
      can2_3->SetHighLightColor(2);
      gStyle->SetOptFit(0111);
      can2_3->SetFillColor(0);//10
      can2_3->SetBorderMode(0);
      can2_3->SetBorderSize(2);
      can2_3->SetFrameFillColor(0);
      can2_3->SetFrameBorderMode(0);
      can2_3->SetFrameBorderMode(0);
      
      TPad *pad2_2 = new TPad("pad4","pad4",0.,0.,1.,1.);
      pad2_2->Draw();
      
      gPad->SetLeftMargin(LeftMargin); gPad->SetRightMargin(0.04);
      gPad->SetTopMargin(0.03); gPad->SetBottomMargin(0.14);

      // G=0 points
      double RatioPoints_highKT3[10]={0, 0.909203, 1.0032, 1.02495, 1.02467, 1.01893, 1.01324, 1.01007, 1.00914, 1.00985};
      double RatioPoints_highKT3_stat[10]={0, 0.0887926, 0.00432693, 0.00118572, 0.000537269, 0.000304217, 0.000194486, 0.000134348, 9.79633e-05, 7.43569e-05};
      double RatioPoints_highKT3_syst[10]={0, 0.015858, 0.0147842, 0.013576, 0.0144809, 0.0174615, 0.0209702, 0.023193, 0.022295, 0.0167607};
      
      TH1D *DoubleRatio3=(TH1D*)Ratio_3->Clone();
      TH1D *DoubleRatio3_G=(TH1D*)Ratio_3_G->Clone();
      for(int i=1; i<=10; i++) {
	if(RatioPoints_highKT3[i-1] == 0) continue;
	double value = Ratio_3->GetBinContent(i) / RatioPoints_highKT3[i-1];
	double value_stat = sqrt(pow(Ratio_3->GetBinError(i)/pow(RatioPoints_highKT3[i-1],2),2) + pow(RatioPoints_highKT3_stat[i-1]*value/RatioPoints_highKT3[i-1],2));
	double value_syst = value * sqrt(pow(Ratio_3_Syst->GetBinError(i),2) + pow(RatioPoints_highKT3_syst[i-1],2) - 2*Ratio_3_Syst->GetBinError(i)*RatioPoints_highKT3_syst[i-1]);
	DoubleRatio3->SetBinContent(i, value);
	DoubleRatio3->SetBinError(i, sqrt(pow(value_stat,2) + pow(value_syst,2)));
	//
	value = Ratio_3_G->GetBinContent(i) / RatioPoints_highKT3[i-1];
	value_stat = sqrt(pow(Ratio_3_G->GetBinError(i)/pow(RatioPoints_highKT3[i-1],2),2) + pow(RatioPoints_highKT3_stat[i-1]*value/RatioPoints_highKT3[i-1],2));
	value_syst = value * sqrt(pow(Ratio_3_Syst->GetBinError(i),2) + pow(RatioPoints_highKT3_syst[i-1],2) - 2*Ratio_3_Syst->GetBinError(i)*RatioPoints_highKT3_syst[i-1]);
	DoubleRatio3_G->SetBinContent(i, value);
	DoubleRatio3_G->SetBinError(i, sqrt(pow(value_stat,2) + pow(value_syst,2)));
      }
      cout<<"Max significance of 3-pion Double ratio from 1.0 = "<<fabs(DoubleRatio3->GetBinContent(3)-1)/DoubleRatio3->GetBinError(3)<<endl;
      
      DoubleRatio3->GetXaxis()->SetTitleSize(SizeTitle);
      DoubleRatio3->GetXaxis()->SetLabelSize(SizeLabel);
      DoubleRatio3->GetYaxis()->SetTitleSize(SizeTitle);
      DoubleRatio3->GetYaxis()->SetLabelSize(SizeLabel);
      DoubleRatio3->GetYaxis()->SetTitleOffset(1.0);
      DoubleRatio3->SetMinimum(0.9); DoubleRatio3->SetMaximum(1.06);
      DoubleRatio3->GetYaxis()->SetNdivisions(303);
      DoubleRatio3->GetYaxis()->SetTitle("Double ratio");
      DoubleRatio3_G->SetFillColor(kBlue-10); DoubleRatio3_G->SetMarkerColor(kBlue); DoubleRatio3_G->SetLineColor(kBlue-10); DoubleRatio3_G->SetLineWidth(5);
      DoubleRatio3->Draw();
      DoubleRatio3_G->Draw("same");
      DoubleRatio3->Draw("same");
      Unity->Draw("same");

      TF1 *lineFit=new TF1("lineFit","pol0",0,.5);
      //DoubleRatio3_G->Fit(lineFit,"","",0.01,0.07);

      /*legend3_2_2->SetTextSize(SizeLegend);
      legend3_2_2->SetX1(0.6); legend3_2_2->SetY1(0.15); legend3_2_2->SetY2(0.3);
      legend3_2_2->AddEntry(DoubleRatio3,"G=0%","p");
      legend3_2_2->AddEntry(DoubleRatio3_G,"G=30%","p");
      legend3_2_2->Draw("same");*/
  }

  


  
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
  gPad->SetLeftMargin(LeftMargin); gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.03); gPad->SetBottomMargin(Xmarge);

  TLegend *legend3= new TLegend(.65,.47, .85,.86,NULL,"brNDC");//.45 or .4 for x1
  legend3->SetBorderSize(0);
  legend3->SetFillColor(0);
  legend3->SetTextFont(TextFont);
  legend3->SetTextSize(SizeLegend);
  if(ChProdBOI!=0) {legend3->SetY1(0.55); legend3->SetY2(0.75);}
  TLegend *legend3_2=(TLegend*)legend3->Clone();
  legend3_2->SetX1(0.6); legend3_2->SetX2(0.9); legend3_2->SetY1(0.34); legend3_2->SetY2(0.6); 
  if(CollisionType!=0) {legend3_2->SetX1(0.75);}
  if(CollisionType==0 && FitBuild==1) legend3_2->SetX1(0.7);
  legend3_2->SetTextSize(2.*SizeLegend);
  TLegend *legend3_2_2=(TLegend*)legend3_2->Clone();

  C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->SetTitleSize(SizeTitle);
  C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->SetLabelSize(SizeLabel);
  C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetYaxis()->SetTitleSize(SizeTitle);
  C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetYaxis()->SetLabelSize(SizeLabel);
  C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->SetTitleOffset(1.05);
  C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetYaxis()->SetTitleOffset(0.8);
  C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->SetNdivisions(606);
  C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetYaxis()->SetNdivisions(505);
  if(CollisionType==0 && ChProdBOI==2) C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetYaxis()->SetNdivisions(504);
 
  //
 
  if(ReNormBuiltBaseline){
    int BinL = C4QSbuilt_G[MBOI]->GetXaxis()->FindBin(ReNormL_4);
    int BinH = C4QSbuilt_G[MBOI]->GetXaxis()->FindBin(ReNormH_4);
    C4Builtsyst_ReNorm = 1 / (C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Integral(BinL,BinH) / C4QSbuilt_G[MBOI]->Integral(BinL,BinH));
    cout<<"C4 2-pion build ReNorm Syst = "<<C4Builtsyst_ReNorm<<endl;
    //
    C4QSbuilt_G[MBOI]->Scale( C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Integral(BinL,BinH) / C4QSbuilt_G[MBOI]->Integral(BinL,BinH));
    a4QSbuilt_G[MBOI]->Scale( a4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Integral(BinL,BinH) / a4QSbuilt_G[MBOI]->Integral(BinL,BinH));
    b4QSbuilt_G[MBOI]->Scale( b4QSmerged[EDBin][MBOI][CollisionType]->Integral(BinL,BinH) / b4QSbuilt_G[MBOI]->Integral(BinL,BinH));
    c4QSbuilt_G[MBOI]->Scale( c4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Integral(BinL,BinH) / c4QSbuilt_G[MBOI]->Integral(BinL,BinH));
    a4QSbuilt_G[MBOI]->SetBinContent(2,0);
    b4QSbuilt_G[MBOI]->SetBinContent(2,0);
    c4QSbuilt_G[MBOI]->SetBinContent(2,0);
    
    for(int bin=1; bin<=C4QSbuilt_G[MBOI]->GetNbinsX(); bin++) {
      C4QSbuilt_G[MBOI]->SetBinError(bin,0);
      a4QSbuilt_G[MBOI]->SetBinError(bin,0);
      b4QSbuilt_G[MBOI]->SetBinError(bin,0);
      c4QSbuilt_G[MBOI]->SetBinError(bin,0);
    }
  }
  int BinL = C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->FindBin(ReNormL_4);
  int BinH = C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->FindBin(ReNormH_4);
  int BinL_syst = C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->FindBin(ReNormL_4_vary);
  int BinH_syst = C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->FindBin(ReNormH_4_vary);
  C4BuiltFromFitsyst_ReNorm = C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Integral(BinL,BinH) / C4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->Integral(BinL,BinH) / (C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Integral(BinL_syst,BinH_syst) / C4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->Integral(BinL_syst,BinH_syst));
  //cout<<"C4 3-pion build ReNorm = "<<C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Integral(BinL,BinH) / C4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->Integral(BinL,BinH)<<endl;
  //cout<<"C4 3-pion build ReNorm Syst = "<<C4BuiltFromFitsyst_ReNorm<<endl;
  //
  C4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->Scale(C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Integral(BinL,BinH) / C4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->Integral(BinL,BinH));
  C4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->Scale(C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Integral(BinL,BinH) / C4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->Integral(BinL,BinH));
  a4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->Scale(a4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Integral(BinL,BinH) / a4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->Integral(BinL,BinH));
  a4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->Scale(a4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Integral(BinL,BinH) / a4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->Integral(BinL,BinH));
  b4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->Scale(b4QSmerged[EDBin][MBOI][CollisionType]->Integral(BinL,BinH) / b4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->Integral(BinL,BinH));
  b4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->Scale(b4QSmerged[EDBin][MBOI][CollisionType]->Integral(BinL,BinH) / b4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->Integral(BinL,BinH));
  c4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->Scale(c4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Integral(BinL,BinH) / c4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->Integral(BinL,BinH));
  c4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->Scale(c4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Integral(BinL,BinH) / c4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]->Integral(BinL,BinH));
  //cout<<"C4 3-pion build ReNorm Syst = "<<C4BuiltFromFitsyst_ReNorm<<endl;
  //}
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
  if(ChProdBOI!=0) C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->SetMinimum(0.6);
  C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->SetMaximum(8.4);
  if(CollisionType==1) C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->SetMaximum(15.0);
  if(CollisionType==2) C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->SetMaximum(19.);
  if(ChProdBOI==1) {
    if(CollisionType==0) C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->SetMaximum(4.4);
    else C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->SetMaximum(6);
  }
  if(ChProdBOI==2) {
    if(CollisionType==0) C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->SetMaximum(3.);
    else C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->SetMaximum(4);
  }
  if(CollisionType==0){
    C4QSBuiltmerged[EDBin][MBOI][CollisionType]->GetXaxis()->SetRange(3,15);
    C4QSbuilt_G[MBOI]->GetXaxis()->SetRange(3,15);
    C4QSbuilt_G[MBOI]->SetBinContent(2,0); a4QSbuilt_G[MBOI]->SetBinContent(2,0); b4QSbuilt_G[MBOI]->SetBinContent(2,0); c4QSbuilt_G[MBOI]->SetBinContent(2,0);
  }
  
  //
  TH1D *C4QS_Syst;
  double shift = 0.00;// 0.002
  if(CollisionType!=0) shift=0.00;// 0.007
  if(ChProdBOI==0) C4QS_Syst= new TH1D("C4QS_Syst","",C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetNbinsX(),C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->GetBinLowEdge(1)+shift, C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->GetBinUpEdge(C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetNbinsX())+shift); 
  else C4QS_Syst = (TH1D*)C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Clone();
  TH1D *a4QS_Syst = (TH1D*)a4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Clone();
  TH1D *c4QS_Syst = (TH1D*)c4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Clone();
  TH1D *C4QSBuilt_Syst;
  if(ChProdBOI==0){
    C4QSBuilt_Syst = (TH1D*)C4QSBuiltmerged[EDBin][MBOI][CollisionType]->Clone();
  }
  for(int bin=1; bin<=C4QS_Syst->GetNbinsX(); bin++){
    double q4 = C4QS_Syst->GetXaxis()->GetBinCenter(bin);
    if(CollisionType!=0) q4 /= 3.1;// rescale q4 for pp and p-Pb
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
      syst1 = sqrt(syst1);
      syst1 *= C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetBinContent(bin);
      C4QS_Syst->SetBinError(bin, syst1);
      // Built
      C4QSBuilt_Syst->SetBinContent(bin, C4QS_Syst->GetBinContent(bin));
      double syst2 = pow(0.004 - 0.004*q4/0.18,2);// 11h to 10h
      syst2 += pow(0.9793 + 0.2857*q4 - 0.9888*q4*q4 - 1,2);// MRC
      syst2 += pow(0.9725 + 0.2991*q4 - 0.8058*q4*q4 - 1,2);// Muon, 92%
      syst2 += pow(0.905 + 1.03*q4 - 2.977*q4*q4 - 1,2);// fc2 scale
      //if(CollisionType==0) {
	//syst2 += pow((0.032-0.245*q4)/2.,2);// old Interpolator.
	syst2 += pow(0.5*(C4ratioFit->Eval(q4)-1),2);// new Interpolator
	syst2 += pow( (1.005 - 0.437*q4 + 13.28*pow(q4,2) - 174.3*pow(q4,3) + 970.3*pow(q4,4) - 1920*pow(q4,5) - 1),2);// Tij sign
	//}
      if(FitBuild || CollisionType!=0) syst2 += pow(C4BuiltFromFitsyst_ReNorm-1,2);// Renormalization
      else syst2 += pow(C4Builtsyst_ReNorm-1,2);// Renormalization 
      
      syst2 = sqrt(syst2);
      syst2 *= C4QSBuiltmerged[EDBin][MBOI][CollisionType]->GetBinContent(bin);
      C4QSBuilt_Syst->SetBinError(bin, syst2);
      // extra cumulant systematic (low Q4 instabilities)
      double syst_c4=0;
      double unscaled_q4 = C4QS_Syst->GetXaxis()->GetBinCenter(bin);
      if(CollisionType==0){
	syst_c4 = (0.004+9.68*exp(-82.4*unscaled_q4))*c4QS_Syst->GetBinContent(bin);
      }else if(CollisionType==1){
	syst_c4 = (0.04+1.11*exp(-12.6*unscaled_q4))*c4QS_Syst->GetBinContent(bin);
      }else{
	syst_c4 = (0.155+2.58*exp(-18.6*unscaled_q4))*c4QS_Syst->GetBinContent(bin);
      }
      c4QS_Syst->SetBinError(bin, syst_c4);
    
    }else{// same % systematics for MC1 and MC2
      double syst1 = pow(0.002,2);// cc
      syst1 += pow(0.002 - 0.002*q4/0.18,2);// 11h to 10h
      syst1 += pow(0.01*(1 - q4/0.18),2);// f coefficients, was pow(0.005,2)
      syst1 += pow(0.0417*exp(-34.1*q4),2);// MRC, 10% of full diff
      syst1 += pow(0.9713 + 0.2648*q4 - 0.752*q4*q4 - 1,2);// Muon, 92%
      syst1 += pow(0.908 + 1.118*q4 - 3.612*q4*q4 - 1,2);// fc2 scale
      //
      C4QS_Syst->SetBinError(bin, sqrt(syst1) * C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetBinContent(bin));
      syst1 += pow(0.536*exp(-58.21*q4),2);
      a4QS_Syst->SetBinError(bin, sqrt(syst1) * a4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetBinContent(bin));
      c4QS_Syst->SetBinError(bin, sqrt(syst1) * c4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetBinContent(bin));
      // extra cumulant systematic (low Q4 instabilities)
      double syst_c4=pow(c4QS_Syst->GetBinError(bin),2);
      double unscaled_q4 = C4QS_Syst->GetXaxis()->GetBinCenter(bin);
      if(CollisionType==0){
	syst_c4 += pow( (0.002+3.31*exp(-91.1*unscaled_q4))*c4QS_Syst->GetBinContent(bin) ,2);
      }else if(CollisionType==1){
	syst_c4 += pow( (0.015+16.25*exp(-46.9*unscaled_q4))*c4QS_Syst->GetBinContent(bin) ,2);
      }else{
	syst_c4 += pow( (-0.008+3.3*exp(-16.4*unscaled_q4))*c4QS_Syst->GetBinContent(bin) ,2);
      }
      c4QS_Syst->SetBinError(bin, sqrt(syst_c4));
    }
  }
  C4QS_Syst->SetMarkerSize(0); C4QS_Syst->SetLineColor(kBlue-10); C4QS_Syst->SetMarkerColor(kBlue-10); 
  C4QS_Syst->SetLineWidth(5); a4QS_Syst->SetLineWidth(5); c4QS_Syst->SetLineWidth(5);
  if(ChProdBOI==0 && CollisionType==0) C4QS_Syst->SetLineWidth(52); 
  if(ChProdBOI==0 && CollisionType!=0) C4QS_Syst->SetLineWidth(34);
  
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
  
  double Syst_forChi2_4[17]={0};
  double Syst_forRatio4[17]={0};
  for(int bin=1; bin<=17; bin++){
    double q4 = C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->GetBinCenter(bin);
    if(CollisionType!=0) q4 /= 3.1;// rescale q4 for pp and p-Pb
    int SystBin = C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->FindBin(q4);
    double SystPercent_Diff = 0;
    SystPercent_Diff += pow(0.01*(1 - q4/0.18),2);// f coefficients was pow(0.9989 - 0.025*q4 - 1,2)
    //SystPercent_Diff += pow((0.032-0.245*q4)/2.,2);// old Interpolator.
    SystPercent_Diff += pow(0.5*(C4ratioFit->Eval(q4)-1),2);// new Interpolator
    SystPercent_Diff += pow( (1.005 - 0.437*q4 + 13.28*pow(q4,2) - 174.3*pow(q4,3) + 970.3*pow(q4,4) - 1920*pow(q4,5) - 1),2);// Tij sign
    SystPercent_Diff += pow( (0.9814 + 0.2471*q4 - 0.8312*q4*q4) - (0.9793 + 0.2857*q4 - 0.9888*q4*q4),2);// MRC diff
    SystPercent_Diff += pow( (0.9635 + 0.3475*q4 - 0.9729*q4*q4) - (0.9725 + 0.2991*q4 - 0.8058*q4*q4),2);// Muon diff
    SystPercent_Diff += pow( (0.900 + 1.126*q4 - 3.354*q4*q4) - (0.905 + 1.03*q4 - 2.977*q4*q4),2);// fc2 scale
    if(FitBuild || CollisionType!=0) SystPercent_Diff += pow(C4BuiltFromFitsyst_ReNorm-1,2);// Renormalization
    else SystPercent_Diff += pow(C4Builtsyst_ReNorm-1,2);// Renormalization
    
    SystPercent_Diff = sqrt(SystPercent_Diff);
    Syst_forChi2_4[bin-1] = SystPercent_Diff * C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetBinContent(bin);
    // ratio
    if(ChProdBOI==0){
      double A = C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetBinContent(bin);
      double B = 0;
      if(!FitBuild) B = C4QSBuiltmerged[EDBin][MBOI][CollisionType]->GetBinContent(bin);
      else B = C4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]->GetBinContent(bin);
      if(B<1) continue;
      double ratio = A / B;
      Syst_forRatio4[bin-1] = pow( ratio * (0.01*(1 - q4/0.18)) ,2);// f coefficients
      //Syst_forRatio4[bin-1] += pow( ratio * 0.0379*exp(-42.82*q4) ,2);// old Interpolator
      Syst_forRatio4[bin-1] += pow( ratio * 0.5*(C4ratioFit->Eval(q4)-1),2);// new Interpolator
      Syst_forRatio4[bin-1] += pow( ratio * (1.005 - 0.437*q4 + 13.28*pow(q4,2) - 174.3*pow(q4,3) + 970.3*pow(q4,4) - 1920*pow(q4,5) - 1) ,2);// Tij sign
      Syst_forRatio4[bin-1] += pow( ratio * sqrt(pow(0.9814 + 0.2471*q4 - 0.8312*q4*q4 - 1,2) + pow(0.9793 + 0.2857*q4 - 0.9888*q4*q4 - 1,2) - 2*(0.9814 + 0.2471*q4 - 0.8312*q4*q4 - 1)*(0.9793 + 0.2857*q4 - 0.9888*q4*q4 - 1) ) ,2);// MRC
      Syst_forRatio4[bin-1] += pow( ratio * sqrt(pow(0.9635 + 0.3475*q4 - 0.9729*q4*q4 - 1,2) + pow(0.9725 + 0.2991*q4 - 0.8058*q4*q4 - 1,2) - 2*(0.9635 + 0.3475*q4 - 0.9729*q4*q4 - 1)*(0.9725 + 0.2991*q4 - 0.8058*q4*q4 - 1) ) ,2);// Muon
      Syst_forRatio4[bin-1] += pow( ratio * sqrt(pow(0.900 + 1.126*q4 - 3.354*q4*q4 - 1,2) + pow(0.905 + 1.03*q4 - 2.977*q4*q4 - 1,2) - 2*(0.900 + 1.126*q4 - 3.354*q4*q4 - 1)*(0.905 + 1.03*q4 - 2.977*q4*q4 - 1) ) ,2);// fc2 scale
      if(CollisionType==0 && !FitBuild) Syst_forRatio4[bin-1] += pow( ratio * (C4Builtsyst_ReNorm-1) ,2);// Renormalization 
      else Syst_forRatio4[bin-1] += pow( ratio * (C4BuiltFromFitsyst_ReNorm-1) ,2);// Renormalization 
      
      Syst_forRatio4[bin-1] = sqrt(Syst_forRatio4[bin-1]);
    }
  }

  //
  
  
  C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->SetBinContent(1,100); C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->SetBinError(1,100);
  C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetXaxis()->SetRangeUser(0,Q4max);
  C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Draw();
  
  C4QS_Syst->Draw("same"); 
  if(ChProdBOI==0) {
    TF1 *SystMeanLine = new TF1("SystMeanLine","[0]",C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetBinLowEdge(BadBinsSC[CollisionType]),1); SystMeanLine->SetParameter(0, C4QS_Syst->GetBinContent(5)); SystMeanLine->SetLineColor(1);
    SystMeanLine->Draw("same");
  }
  if(!ShortSameCharge) c4QS_Syst->Draw("same");
  //if(ChProdBOI==0) C4QSBuilt_Syst->Draw("same");// E1
  if(ChProdBOI!=0){    
    a4QS_Syst->Draw("same");
  }
  
  C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Draw("same");
  if(!ShortSameCharge) a4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Draw("same");
  if(ChProdBOI==0 && !ShortSameCharge ) b4QSmerged[EDBin][MBOI][CollisionType]->Draw("same");
  if(!ShortSameCharge) c4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Draw("same");
  if(ChProdBOI==0) {
    if(FitBuild==kFALSE){
      C4QSBuiltmerged[EDBin][MBOI][CollisionType]->SetLineWidth(1.2);
      C4QSBuiltmerged[EDBin][MBOI][CollisionType]->Draw("same");
      if(!ShortSameCharge) {
	a4QSBuiltmerged[EDBin][MBOI][CollisionType]->Draw("same");
	b4QSBuiltmerged[EDBin][MBOI][CollisionType]->Draw("same");
	c4QSBuiltmerged[EDBin][MBOI][CollisionType]->Draw("same");
      }
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
  
  
  if(ChProdBOI==0 && FitBuild==0 && ShortSameCharge) {
    C4QSbuilt_G[MBOI]->Draw("same");
    //a4QSbuilt_G[MBOI]->Draw("same");
    //b4QSbuilt_G[MBOI]->Draw("same");
    //c4QSbuilt_G[MBOI]>Draw("same");
  }

  legend3->AddEntry(C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType],"#font[12]{C}_{4}^{QS}","p");
  if(!ShortSameCharge) legend3->AddEntry(a4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType],"#font[12]{#bf{a}}_{4}^{QS}","p");
  if(ChProdBOI==0 && !ShortSameCharge) legend3->AddEntry(b4QSmerged[EDBin][MBOI][CollisionType],"#font[12]{#bf{b}}_{4}^{QS}","p");
  if(!ShortSameCharge) legend3->AddEntry(c4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType],"#font[12]{#bf{c}}_{4}^{QS}","p");
  if(ChProdBOI==0 && !FitBuild) {
    legend3->AddEntry(C4QSBuiltmerged[EDBin][MBOI][CollisionType],"#font[12]{E}_{4}(2) (G=0%)","l");
    if(ShortSameCharge) legend3->AddEntry(C4QSbuilt_G[MBOI], BuiltNameE42->Data(),"l");
    //legend3->AddEntry(C4QSBuilt_Syst,"Built #font[12]{C}_{4}^{QS} systematic","l");
  }
  if(ChProdBOI==0){
    if(FitBuild){
      legend3->AddEntry(C4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType], BuiltNamee43->Data(),"l");
      legend3->AddEntry(C4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType], BuiltNameE43->Data(),"l");
      //legend3->AddEntry(C4QSBuilt_Syst,"Built #font[12]{C}_{4}^{QS} systematic","l");
    }
  }
  //if(ChProdBOI==0) legend3->AddEntry(C4QS_Syst,"#font[12]{C}_{4}^{QS} systematic","l");
  

  
 
  TPad *pad3_3 = new TPad("pad3_3","pad3_3",0.52,0.6,.95,.91);
  if(ChProdBOI!=0) {
    KT4specif->SetX(0.62); KT4specif->SetY(0.52);
    legend3->SetX1(0.42); legend3->SetX2(0.52); legend3->SetY1(0.68); legend3->SetY2(0.92);
    pad3_3->Draw();
    pad3_3->cd();
    gPad->SetLeftMargin(0.22); gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01); gPad->SetBottomMargin(0.14);
    TH1D *c4QS_zoom=(TH1D*)c4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Clone();
    TH1D *c4QS_Syst_zoom=(TH1D*)c4QS_Syst->Clone();
    float ZoomMin=0.84, ZoomMax=1.44;
    if(CollisionType==0) {ZoomMin=0.86, ZoomMax=1.14;}
    c4QS_zoom->SetMinimum(ZoomMin); c4QS_zoom->SetMaximum(ZoomMax);
    c4QS_zoom->GetXaxis()->SetRangeUser(0,Q4max);
    c4QS_zoom->GetXaxis()->SetNdivisions(202); c4QS_zoom->GetYaxis()->SetNdivisions(204);
    if(CollisionType!=0) c4QS_zoom->GetXaxis()->SetNdivisions(204);
    c4QS_zoom->GetXaxis()->SetLabelSize(2.5*SizeLabel); c4QS_zoom->GetXaxis()->SetTitleOffset(100);
    c4QS_zoom->GetYaxis()->SetLabelSize(2.5*SizeLabel); c4QS_zoom->GetYaxis()->SetTitleOffset(100);
    c4QS_zoom->GetYaxis()->SetLabelOffset(0.02);
    c4QS_zoom->Draw();
    c4QS_Syst->Draw("same");
    c4QS_zoom->Draw("same");
    Unity->Draw("same");
  }

  pad3->cd();
  ALICEspecif->Draw("same");
  KT4specif->Draw("same");

  legend3->Draw("same");
  Unity->Draw("same");
  ChargeCombTitle4->Draw("same");

  /*TF1 *Gauss_c4Fit=new TF1("Gauss_c4Fit","[0]*(1+[1]*exp(-pow(x*[2]/0.19733,2)/3.))",0,1);
  Gauss_c4Fit->SetParameter(0,1); Gauss_c4Fit->SetParameter(1,3); Gauss_c4Fit->SetParameter(2,8); 
  Gauss_c4Fit->SetParName(0,"N");  Gauss_c4Fit->SetParName(1,"#lambda_{4}"); Gauss_c4Fit->SetParName(2,"R");
  c4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Fit(Gauss_c4Fit,"IME","",0.03,0.14);
  Gauss_c4Fit->Draw("same");*/

  // high ED4 reference
  double y_ref[12]={0, 0, 1.00133, 0.980848, 0.988251, 0.994434, 0.999677, 1.00269, 1.00642, 1.00881, 1.01082, 1.01554};
  double y_ref_e[12]={0, 0, 0.054465, 0.00678447, 0.00194947, 0.000799564, 0.00039767, 0.000222628, 0.000135335, 8.75305e-05, 6.31392e-05, 5.53329e-05};
  double y_ref_Syst[12]={0, 0, 0.0167026, 0.00971243, 0.00945824, 0.0121959, 0.014654, 0.0155188, 0.0147369, 0.0128801, 0.0109854, 0.0100858};// without MC1 residue
  if(ChProdBOI==0){
    pad3_2->cd();
    gPad->SetLeftMargin(LeftMargin); gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.03); gPad->SetBottomMargin(0.32);
    
    TH1D *Ratio_4=(TH1D*)C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Clone();
    Ratio_4->GetXaxis()->SetTitleSize(2.2*SizeTitle);
    Ratio_4->GetXaxis()->SetLabelSize(2.2*SizeLabel);
    Ratio_4->GetYaxis()->SetTitleSize(2.2*SizeTitle);
    Ratio_4->GetYaxis()->SetLabelSize(2.2*SizeLabel);
    Ratio_4->GetXaxis()->SetTitleOffset(1.05);
    Ratio_4->GetYaxis()->SetTitleOffset(0.4);
    Ratio_4->GetXaxis()->SetNdivisions(606);
    Ratio_4->GetYaxis()->SetNdivisions(204);
    
    Ratio_4->GetYaxis()->SetTitle("Ratio ");
    if(CollisionType==0){
      Ratio_4->SetMinimum(0.9); 
      if(FitBuild) {Ratio_4->SetMinimum(0.86); Ratio_4->SetMaximum(1.12);}
      else Ratio_4->SetMaximum(1.03);
    }else {Ratio_4->SetMinimum(0.7); Ratio_4->SetMaximum(1.28);}
    TH1D *Ratio_4_2=(TH1D*)Ratio_4->Clone();
    TH1D *Ratio_4_G=(TH1D*)Ratio_4->Clone();
    Ratio_4->SetMarkerStyle(24);
    Ratio_4_G->SetMarkerStyle(20);
    if(FitBuild==0){
      Ratio_4->Divide(C4QSBuiltmerged[EDBin][MBOI][CollisionType]);
      Ratio_4_G->Divide(C4QSbuilt_G[MBOI]);
    }else{
      Ratio_4->Divide(C4QSBuiltFromFitsmerged1[EDBin][MBOI][CollisionType]);
      Ratio_4_2->Divide(C4QSBuiltFromFitsmerged2[EDBin][MBOI][CollisionType]);
    }
    TH1D *Ratio_4_Syst=(TH1D*)Ratio_4->Clone();
    Ratio_4_Syst->SetMarkerSize(0); Ratio_4_Syst->SetMarkerColor(kBlue-10); Ratio_4_Syst->SetFillColor(kBlue-10);
    Ratio_4_Syst->SetLineColor(kBlue-10);
    for(int bin=1; bin<=17; bin++){
      Ratio_4_Syst->SetBinError(bin, Syst_forRatio4[bin-1]);
    }
    
    Ratio_4_Syst->SetLineWidth(5);
    Ratio_4_Syst->Draw();
    if(FitBuild==0){
      Ratio_4->Draw("same");
      if(ShortSameCharge) Ratio_4_G->Draw("same");
      legend3_2->AddEntry(Ratio_4,"#font[12]{C}_{4}^{QS}/#font[12]{E}_{4}(2) (G=0%)","p");
      if(ShortSameCharge) {
	TString *name = new TString("#font[12]{C}_{4}^{QS}/"); name->Append(BuiltNameE42->Data());
	legend3_2->AddEntry(Ratio_4_G,name->Data(),"p");
      }
    }else{
     Ratio_4->Draw("same");
     Ratio_4_2->Draw("same");
     legend3_2->AddEntry(Ratio_4,"#font[12]{C}_{4}^{QS}/#font[12]{E}_{4}(3)" ,"p");
     legend3_2->AddEntry(Ratio_4_2,"#font[12]{C}_{4}^{QS}/#font[12]{e}_{4}(3)","p");
    }
      
    legend3_2->Draw("same");
    Unity->Draw("same");

    
    // New Chi2 of Ratio calculation
    double chi2=0, Npoints=0;
    int binStart=3, binEnd=5;
    if(CollisionType==0 && FitBuild) {binStart=4; binEnd=8;}
    if(CollisionType==1) {binStart=4; binEnd=8;}
    if(CollisionType==2) {binStart=5; binEnd=9;}
    for(int bin=binStart; bin<=binEnd; bin++){// 3 to 5 for quoted G value, 3 to 8 for FitBuild
      double TotalErr = sqrt(pow(Ratio_4->GetBinError(bin),2) + pow(Ratio_4_Syst->GetBinError(bin),2));
      if(TotalErr<=0) continue;
      if(FitBuild==0 && CollisionType==0) {
	chi2 += pow( (Ratio_4_G->GetBinContent(bin) - 1.0) / TotalErr ,2);
	Npoints++;
      }else{
	chi2 += pow( (Ratio_4->GetBinContent(bin) - 1.0) / TotalErr ,2);
	chi2 += pow( (Ratio_4_2->GetBinContent(bin) - 1.0) / TotalErr ,2);
	cout<<(Ratio_4->GetBinContent(bin) - 1.0)<<"  "<<(Ratio_4_2->GetBinContent(bin) - 1.0)<<endl;
	Npoints++;
	Npoints++;
      }
    }
    //if(Npoints>0) cout<<"Ratio sqrt(Chi2/NDF) = "<<sqrt(chi2/Npoints)<<endl;
  
    int NbinsDR=11;
    if(CollisionType!=0) NbinsDR=17;
    //
    /*for(int i=1; i<=NbinsDR; i++) cout<<Ratio_4->GetBinContent(i)<<", ";
    cout<<endl;
    for(int i=1; i<=NbinsDR; i++) cout<<Ratio_4->GetBinError(i)<<", ";
    cout<<endl;
    for(int i=1; i<=NbinsDR; i++) cout<<Ratio_4_Syst->GetBinError(i)<<", ";
    cout<<endl;*/
    
    double RatioPoints_highKT4[20]={0};
    double RatioPoints_highKT4_stat[20]={0};
    double RatioPoints_highKT4_syst[20]={0};
    // Draw only the Ratios
    if(!FitBuild && MBOI==0){
      TCanvas *can3_2 = new TCanvas("can3_2", "can3_2",100,700,700,700);// 11,53,700,500
      can3_2->SetHighLightColor(2);
      gStyle->SetOptFit(0111);
      can3_2->SetFillColor(0);//10
      can3_2->SetBorderMode(0);
      can3_2->SetBorderSize(2);
      can3_2->SetFrameFillColor(0);
      can3_2->SetFrameBorderMode(0);
      can3_2->SetFrameBorderMode(0);
      
      TPad *pad4 = new TPad("pad4","pad4",0.,0.,1.,1.);
      pad4->Draw();
      
      gPad->SetLeftMargin(LeftMargin); gPad->SetRightMargin(0.04);
      gPad->SetTopMargin(0.03); gPad->SetBottomMargin(0.14);
      
      // G=0 points
      if(CollisionType==0){
	double RatioPoints_highKT4_CT0[11]={0, 0, 0.945281, 1.01302, 1.01241, 1.00672, 1.00237, 0.9996, 0.999973, 1.00003, 1.00037};
	double RatioPoints_highKT4_stat_CT0[11]={0,0, 0.0444709, 0.00599654, 0.00170795, 0.000693098, 0.000341785, 0.000190292, 0.000115303, 7.44235e-05, 5.36079e-05};
	double RatioPoints_highKT4_syst_CT0[11]={0, 0, 0.0221586, 0.0190284, 0.0181013, 0.017631, 0.0173482, 0.0172557, 0.0172789, 0.0172188, 0.0170602};
	for(int i=1; i<=NbinsDR; i++){
	  RatioPoints_highKT4[i] = RatioPoints_highKT4_CT0[i];
	  RatioPoints_highKT4_stat[i] = RatioPoints_highKT4_stat_CT0[i];
	  RatioPoints_highKT4_syst[i] = RatioPoints_highKT4_syst_CT0[i];
	}
      }else if(CollisionType==1){
	double RatioPoints_highKT4_CT1[17]={0, 0, 0, 1.33908, 1.16369, 1.17244, 1.10655, 1.07037, 1.05464, 1.02833, 1.01907, 1.00715, 1.00157, 0.99995, 0.998861, 0.999809, 1.0002};
	double RatioPoints_highKT4_stat_CT1[17]={0, 0, 61.6273, 0.0995355, 0.0317487, 0.0147504, 0.00761763, 0.00450887, 0.00294661, 0.00204392, 0.00150831, 0.00115896, 0.000923693, 0.000759155, 0.00063769, 0.00054585, 0.000474241};
	double RatioPoints_highKT4_syst_CT1[17]={0, 0, 0, 0.0279524, 0.0164417, 0.01286, 0.0108507, 0.00988162, 0.00928589, 0.00873351, 0.00848218, 0.008343, 0.00832269, 0.00831346, 0.00821794, 0.00803194, 0.00778718};
	for(int i=1; i<=NbinsDR; i++){
	  RatioPoints_highKT4[i] = RatioPoints_highKT4_CT1[i];
	  RatioPoints_highKT4_stat[i] = RatioPoints_highKT4_stat_CT1[i];
	  RatioPoints_highKT4_syst[i] = RatioPoints_highKT4_syst_CT1[i];
	}
      }else{
	double RatioPoints_highKT4_CT2[17]={0, 0, 0, 0, 0, 1.38333, 1.37205, 1.28794, 1.22575, 1.1788, 1.14095, 1.09622, 1.06758, 1.03893, 1.02244, 1.00749, 0.991955};
	double RatioPoints_highKT4_stat_CT2[17]={0, 0, 0, 49.1391, 15.7673, 0.0634993, 0.0341033, 0.0192716, 0.0120421, 0.00817051, 0.00590278, 0.00442957, 0.00348142, 0.00281178, 0.00235022, 0.00200092, 0.00172986};
	double RatioPoints_highKT4_syst_CT2[17]={0, 0, 0, 0, 0, 0.0280514, 0.0269931, 0.0249782, 0.0235274, 0.0224598, 0.0216532, 0.0207852, 0.0202543, 0.0197127, 0.0193613, 0.0189941, 0.0185978};
	for(int i=1; i<=NbinsDR; i++){
	  RatioPoints_highKT4[i] = RatioPoints_highKT4_CT2[i];
	  RatioPoints_highKT4_stat[i] = RatioPoints_highKT4_stat_CT2[i];
	  RatioPoints_highKT4_syst[i] = RatioPoints_highKT4_syst_CT2[i];
	}
      }
      

      TH1D *DoubleRatio=(TH1D*)Ratio_4->Clone();
      TH1D *DoubleRatio_G=(TH1D*)Ratio_4_G->Clone();
      for(int i=1; i<=NbinsDR; i++) {
	if(RatioPoints_highKT4[i-1] == 0) continue;
	double value = Ratio_4->GetBinContent(i) / RatioPoints_highKT4[i-1];
	double value_stat = sqrt(pow(Ratio_4->GetBinError(i)/pow(RatioPoints_highKT4[i-1],2),2) + pow(RatioPoints_highKT4_stat[i-1]*value/RatioPoints_highKT4[i-1],2));
	double value_syst = value * sqrt(pow(Ratio_4_Syst->GetBinError(i),2) + pow(RatioPoints_highKT4_syst[i-1],2) - 2*Ratio_4_Syst->GetBinError(i)*RatioPoints_highKT4_syst[i-1]);
	DoubleRatio->SetBinContent(i, value);
	DoubleRatio->SetBinError(i, sqrt(pow(value_stat,2) + pow(value_syst,2)));
	//
	value = Ratio_4_G->GetBinContent(i) / RatioPoints_highKT4[i-1];
	value_stat = sqrt(pow(Ratio_4_G->GetBinError(i)/pow(RatioPoints_highKT4[i-1],2),2) + pow(RatioPoints_highKT4_stat[i-1]*value/RatioPoints_highKT4[i-1],2));
	value_syst = value * sqrt(pow(Ratio_4_Syst->GetBinError(i),2) + pow(RatioPoints_highKT4_syst[i-1],2) - 2*Ratio_4_Syst->GetBinError(i)*RatioPoints_highKT4_syst[i-1]);
	DoubleRatio_G->SetBinContent(i, value);
	DoubleRatio_G->SetBinError(i, sqrt(pow(value_stat,2) + pow(value_syst,2)));
      }
      cout<<"Max significance of 4-pion Double ratio from 1.0 = "<<fabs(DoubleRatio->GetBinContent(4)-1)/DoubleRatio->GetBinError(4)<<endl;
      

      DoubleRatio->GetXaxis()->SetTitleSize(SizeTitle);
      DoubleRatio->GetXaxis()->SetLabelSize(SizeLabel);
      DoubleRatio->GetYaxis()->SetTitleSize(SizeTitle);
      DoubleRatio->GetYaxis()->SetLabelSize(SizeLabel);
      DoubleRatio->GetYaxis()->SetTitleOffset(0.8);
      DoubleRatio->SetMinimum(0.9); DoubleRatio->SetMaximum(1.06);
      DoubleRatio->GetYaxis()->SetNdivisions(303);
      DoubleRatio->GetYaxis()->SetTitle("Double ratio");
      DoubleRatio_G->SetFillColor(kBlue-10); DoubleRatio_G->SetMarkerColor(kBlue); DoubleRatio_G->SetLineColor(kBlue-10); DoubleRatio_G->SetLineWidth(5);
      DoubleRatio->Draw();
      if(CollisionType==0) DoubleRatio_G->Draw("same");
      DoubleRatio->Draw("same");
      Unity->Draw("same");

      if(CollisionType==0){
	TF1 *lineFit=new TF1("lineFit","pol0",0,.5);
	//DoubleRatio_G->Fit(lineFit,"","",0.03,0.165);
	
	legend3_2_2->SetTextSize(SizeLegend);
	legend3_2_2->SetX1(0.6); legend3_2_2->SetY1(0.15); legend3_2_2->SetY2(0.3);
	legend3_2_2->AddEntry(DoubleRatio,"G=0%","p");
	legend3_2_2->AddEntry(DoubleRatio_G,"G=30%","p");
	legend3_2_2->Draw("same");
      }
      /*
      TH1D *Ratio_4_Syst_2 = (TH1D*)Ratio_4_Syst->Clone();
      Ratio_4_Syst_2->GetXaxis()->SetTitleSize(SizeTitle);
      Ratio_4_Syst_2->GetXaxis()->SetLabelSize(SizeLabel);
      Ratio_4_Syst_2->GetYaxis()->SetTitleSize(SizeTitle);
      Ratio_4_Syst_2->GetYaxis()->SetLabelSize(SizeLabel);
      Ratio_4_Syst_2->GetYaxis()->SetTitleOffset(1.0);
      Ratio_4_Syst_2->SetMinimum(0.94); Ratio_4_Syst_2->SetMaximum(1.06);
      Ratio_4_Syst_2->GetYaxis()->SetNdivisions(303);
      Ratio_4_Syst_2->GetYaxis()->SetTitle("#font[12]{C_{4}^{QS}} / built           ");
      Ratio_4->SetMarkerStyle(20); Ratio_4_2->SetMarkerStyle(21);

      double c3Ratio_points_G0[12]={0, 0, 0.979954, 1.02074, 1.03611, 1.04127, 1.03673, 1.02607, 1.0151, 1.00445, 0.995367, 0.988652};
      double C3Ratio_points_G0[12]={0, 0, 0.985206, 0.974207, 0.962165, 0.964346, 0.976318, 0.98772, 0.995432, 0.998735, 1.00133, 1.00622};
      //for(int i=0; i<12; i++) cout<<Ratio_4->GetBinContent(i+1)<<", ";
      //cout<<endl;
      //for(int i=0; i<12; i++) cout<<Ratio_4_2->GetBinContent(i+1)<<", ";
      //cout<<endl;
      TH1D *Ratio_c3_G0 = (TH1D*)Ratio_4->Clone();
      TH1D *Ratio_C3_G0 = (TH1D*)Ratio_4_2->Clone();
      for(int i=0; i<12; i++) {
	Ratio_c3_G0->SetBinContent(i+1, c3Ratio_points_G0[i]);
	Ratio_C3_G0->SetBinContent(i+1, C3Ratio_points_G0[i]);
      }
      Ratio_c3_G0->SetMarkerStyle(24);
      Ratio_C3_G0->SetMarkerStyle(25);
      
      Ratio_4_Syst_2->SetBinContent(3,0); Ratio_4->SetBinContent(3,0); Ratio_4_2->SetBinContent(3,0);
      Ratio_c3_G0->SetBinContent(3,0); Ratio_C3_G0->SetBinContent(3,0); 

      Ratio_4_Syst_2->Draw();
      Ratio_4->Draw("same");
      Ratio_4_2->Draw("same");
      TLegend *legend3_2_2=(TLegend*)legend3_2->Clone();
      legend3_2_2->SetTextSize(0.7*SizeLegend);
      legend3_2_2->SetX1(0.6); legend3_2_2->SetY1(0.15); legend3_2_2->SetY2(0.3);
      Ratio_c3_G0->Draw("same");
      Ratio_C3_G0->Draw("same");
      legend3_2_2->AddEntry(Ratio_c3_G0,"Built from #font[12]{#bf{c}}_{3} (G=0%)","p");
      legend3_2_2->AddEntry(Ratio_C3_G0,"Built from #font[12]{C}_{3} (G=0%)","p");
      legend3_2_2->Draw("same");
      Unity->Draw("same");
      */


      TLatex *ALICEspecif_3 = (TLatex*)ALICEspecif->Clone();
      //TLatex *KT4specif_3 = (TLatex*)KT4specif->Clone();
      ALICEspecif_3->SetX(0.33);
      //KT4specif_3->SetX(0.58); KT4specif_3->SetY(0.85); 
      ALICEspecif_3->Draw("same");
      //KT4specif_3->Draw("same");
    }
    
  }

 

  if(SaveFiles){
    TString *SaveName = new TString("PlotsForPapers/Fig_C4");
    if(ChProdBOI==0) SaveName->Append("SC_");
    else if(ChProdBOI==1) SaveName->Append("MC1_");
    else SaveName->Append("MC2_");
    //
    if(CollisionType==0) SaveName->Append("PbPb_");
    else if(CollisionType==1) SaveName->Append("pPb_");
    else SaveName->Append("pp_");
    //
    if(ChProdBOI==0){
      if(CollisionType!=0){
	SaveName->Append("FitBuilds_");
      }
      if(CollisionType==0){
	if(FitBuild) SaveName->Append("FitBuilds_");
	else{
	  if(ShortSameCharge) SaveName->Append("short_");
	  else SaveName->Append("long_");
	}
      }
    }
    
    SaveName->Append("K");
    *SaveName += EDBin;
    SaveName->Append("_M");
    *SaveName += MBOI;
    SaveName->Append(".eps");
    //
    can3->SaveAs(SaveName->Data());
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
    gPad->SetLeftMargin(LeftMargin); gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.03); gPad->SetBottomMargin(0.14);

    TLegend *legend4 = new TLegend(.15,.65, .4,.85,NULL,"brNDC");//.45 or .4 for x1
    legend4->SetBorderSize(0);
    legend4->SetFillColor(0);
    legend4->SetTextFont(TextFont);
    legend4->SetTextSize(SizeLegend);
    TLegend *legend4_2=(TLegend*)legend4->Clone();
    legend4_2->SetX1(0.7); legend4_2->SetX2(0.95); legend4_2->SetY1(0.6); legend4_2->SetY2(0.95); 

    
    TH1D *chi2_4 = new TH1D("chi2_4","",100,-0.5,99.5);
    chi2_4->SetLineColor(2);
    chi2_4->SetMarkerColor(2);
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
      tempNum->Add(tempNumNeg);
      tempNum->Divide(tempDen);
      if(MBOI!=0 && MBOI!=1 && MBOI!=2 && MBOI!=4 && MBOI!=7) cout<<"Problem with choice of merged centrality bin!!!!!!!!!!!!!!"<<endl;
      if(MBOI==2){
	TH1D *tempNum2 = (TH1D*)C4QSBuiltmerged2D[EDBin][MBOI+1][CollisionType]->ProjectionY("2ndMB", binG, binG);
	TH1D *tempNumNeg2 = (TH1D*)C4QSNegBuiltmerged2D[EDBin][MBOI+1][CollisionType]->ProjectionY("2ndMB_Neg", binG, binG);
	TH1D *tempDen2 = (TH1D*)C4QSBuiltmerged2D[EDBin][MBOI+1][CollisionType]->ProjectionY("2ndDen", 4, 4);
	TH1D *tempDenNeg2 = (TH1D*)C4QSNegBuiltmerged2D[EDBin][MBOI+1][CollisionType]->ProjectionY("2ndDen_Neg", 4, 4);
	tempNum2->Add(tempNumNeg2); tempDen2->Add(tempDenNeg2); 
	tempNum2->Divide(tempDen2);
	tempNum->Add(tempNum2);
	tempNum->Scale(0.5);
      }
      if(MBOI==4 || MBOI==7){
	TH1D *tempNum2 = (TH1D*)C4QSBuiltmerged2D[EDBin][MBOI+1][CollisionType]->ProjectionY("2ndMB", binG, binG);
	TH1D *tempNumNeg2 = (TH1D*)C4QSNegBuiltmerged2D[EDBin][MBOI+1][CollisionType]->ProjectionY("2ndMB_Neg", binG, binG);
	TH1D *tempDen2 = (TH1D*)C4QSBuiltmerged2D[EDBin][MBOI+1][CollisionType]->ProjectionY("2ndDen", 4, 4);
	TH1D *tempDenNeg2 = (TH1D*)C4QSNegBuiltmerged2D[EDBin][MBOI+1][CollisionType]->ProjectionY("2ndDen_Neg", 4, 4);
	tempNum2->Add(tempNumNeg2); tempDen2->Add(tempDenNeg2); 
	tempNum2->Divide(tempDen2);
	TH1D *tempNum3 = (TH1D*)C4QSBuiltmerged2D[EDBin][MBOI+2][CollisionType]->ProjectionY("3rdMB", binG, binG);
	TH1D *tempNumNeg3 = (TH1D*)C4QSNegBuiltmerged2D[EDBin][MBOI+2][CollisionType]->ProjectionY("3rdMB_Neg", binG, binG);
	TH1D *tempDen3 = (TH1D*)C4QSBuiltmerged2D[EDBin][MBOI+2][CollisionType]->ProjectionY("3rdDen", 4, 4);
	TH1D *tempDenNeg3 = (TH1D*)C4QSNegBuiltmerged2D[EDBin][MBOI+2][CollisionType]->ProjectionY("3rdDen_Neg", 4, 4);
	tempNum3->Add(tempNumNeg3); tempDen3->Add(tempDenNeg3); 
	tempNum3->Divide(tempDen3);
	//
	tempNum->Add(tempNum2);
	tempNum->Add(tempNum3);
	tempNum->Scale(1/3.);
      }
      
      //
      if(ReNormBuiltBaseline){
	int BinL = tempNum->GetXaxis()->FindBin(ReNormL_4);
	int BinH = tempNum->GetXaxis()->FindBin(ReNormH_4);
	if(tempNum->Integral(BinL,BinH)==0) continue;
	tempNum->Scale( C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->Integral(BinL,BinH) / tempNum->Integral(BinL,BinH));
      }
      
      
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
      
      if(tempNum->GetBinContent(Q4binChi2) >0) {
	double value = fabs(C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetBinContent(Q4binChi2) - tempNum->GetBinContent(Q4binChi2));
	double err = pow(C4QSmerged[ChProdBOI][EDBin][MBOI][CollisionType]->GetBinError(Q4binChi2),2);
	err += pow(Syst_forChi2_4[Q4binChi2-1],2);
	err = sqrt(err);
	if(err>0) {
	  double Chi2 = pow(value / err,2);
	  chi2_4->SetBinContent(1 + 2*(binG-Gbinstart), sqrt(fabs(Chi2))); 
	  chi2_4->SetBinError(1 + 2*(binG-Gbinstart), 0.001);
	}
      }
      
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
    //
    chi2_4->Draw();

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
    TString *CentAndED_5 = new TString(Centname->Data());
    if(EDBin==0) CentAndED_5->Append(", 0.16<#font[12]{K}_{T4}<0.3 GeV/#font[12]{c}");
    else CentAndED_5->Append(", 0.3<#font[12]{K}_{T4}<1.0 GeV/#font[12]{c}");
    TLatex *Cent_5 = new TLatex(2,9,CentAndED_5->Data());
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
    gPad->SetLeftMargin(LeftMargin); gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.03); gPad->SetBottomMargin(0.14);

    TLegend *legend5 = new TLegend(.15,.75, .3,.95,NULL,"brNDC");//.45 or .4 for x1
    legend5->SetBorderSize(0);
    legend5->SetFillColor(0);
    legend5->SetTextFont(TextFont);
    legend5->SetTextSize(SizeLegend);

    TH1D *GversusQ4[3];// Stat then Syst then total
    GversusQ4[0] = new TH1D("GversusQ4Stat","",7,0,0.105);// stat
    GversusQ4[1] = new TH1D("GversusQ4Syst","",7,0,0.105);// syst
    GversusQ4[2] = new TH1D("GversusQ4Total","",7,0,0.105);// stat+syst
    for(int ErrType=0; ErrType<2; ErrType++){
      
      for(int binQ4=3; binQ4<=7; binQ4++){
	double minG = 0;
	double minG_e1 = 0, minG_e2=0;
	double minChi=100;
	
	for(int binG=1; binG<=25; binG++){// min
	  if(minChi > chi2_2D_4[ErrType]->GetBinContent(binQ4, binG)) {
	    minChi = chi2_2D_4[ErrType]->GetBinContent(binQ4, binG);
	    minG = 2*(binG-1);
	  }
	}
	double PastDiff=100;
	for(int binG=1; binG<=25; binG++){// error
	  double Diff=fabs(minChi - chi2_2D_4[ErrType]->GetBinContent(binQ4, binG));
	  if(minG > 0) {
	    if(fabs(Diff-1)<PastDiff) {
	      if(minG>2*(binG-1)) minG_e1 = fabs(minG - 2*(binG-1)); 
	      else minG_e2 = fabs(minG - 2*(binG-1));
	      PastDiff=fabs(Diff-1);
	    }
	  }else{
	    if(fabs(Diff-1)<PastDiff) {
	      minG_e1 = fabs(minG - 2*(binG-1)); 
	      PastDiff=fabs(Diff-1);
	    }
	  }
	}
	GversusQ4[ErrType]->SetBinContent(binQ4, minG);
	if(minG==0) GversusQ4[ErrType]->SetBinContent(binQ4, 110);
	if(minG_e1>minG_e2) GversusQ4[ErrType]->SetBinError(binQ4, minG_e1);
	else GversusQ4[ErrType]->SetBinError(binQ4, minG_e2);
	//
      }
    }// ErrType
    //
    for(int binQ4=3; binQ4<=7; binQ4++){
      GversusQ4[2]->SetBinContent(binQ4, GversusQ4[1]->GetBinContent(binQ4));
      GversusQ4[2]->SetBinError(binQ4, sqrt(pow(GversusQ4[0]->GetBinError(binQ4),2)+pow(GversusQ4[1]->GetBinError(binQ4),2)));
    }
    
  
    GversusQ4[0]->SetMarkerStyle(20); GversusQ4[0]->SetMarkerColor(2); GversusQ4[0]->SetLineColor(2);
    GversusQ4[1]->SetMarkerSize(0); GversusQ4[1]->SetMarkerColor(kRed-10); GversusQ4[1]->SetLineColor(kRed-10);
    GversusQ4[1]->SetFillStyle(0); GversusQ4[1]->SetLineWidth(5);
        
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
    //
    //
    TLatex *System_6 = new TLatex(0.03,82,System->Data());// ALICE specifications
    System_6->SetTextFont(TextFont);
    System_6->SetTextSize(SizeSpecif);
    System_6->Draw("same");
    //
    TLatex *Cent_6 = new TLatex(0.03,74,Centname->Data());
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
    double FitStart = GversusQ4[2]->GetXaxis()->GetBinCenter(binStartFinalG);
    double FitEnd = GversusQ4[2]->GetXaxis()->GetBinCenter(binEndFinalG);
    //GversusQ4[0]->Fit("pol0","IMEN","",FitStart, FitEnd);
    //pol0->Draw("same");
    double Avg_G=0, Avg_G_den=0;
    double Avg_G_stat=0, Avg_G_stat_den=0;
    double Avg_G_syst=0, Avg_G_syst_den=0;
    for(int bin=binStartFinalG; bin<=binEndFinalG; bin++){
      if(GversusQ4[0]->GetBinError(bin)==0) continue;
      if(GversusQ4[1]->GetBinError(bin)==0) continue;
      if(GversusQ4[2]->GetBinError(bin)==0) continue;
      if(GversusQ4[0]->GetBinContent(bin)>100) continue;
      double weight = 1 / GversusQ4[2]->GetBinError(bin);
      Avg_G += GversusQ4[2]->GetBinContent(bin) * weight;
      Avg_G_den += weight;
      weight = 1 / GversusQ4[0]->GetBinError(bin);
      Avg_G_stat += GversusQ4[0]->GetBinError(bin) * weight;
      Avg_G_stat_den += weight;
      weight = 1 / GversusQ4[1]->GetBinError(bin);
      Avg_G_syst += GversusQ4[1]->GetBinError(bin) * weight;
      Avg_G_syst_den += weight;
    }
    
    if(Avg_G_den>0) Avg_G = int(Avg_G/Avg_G_den +0.5);
    if(Avg_G_syst_den>0) Avg_G_syst = int(Avg_G_syst/Avg_G_syst_den +0.5);
    if(Avg_G_stat_den>0) Avg_G_stat = int( Avg_G_stat/Avg_G_stat_den +0.5);
    cout<<Avg_G<<" pm "<<Avg_G_stat<<" pm "<<Avg_G_syst<<endl;

    TString *GprintString = new TString("#LT G #GT = ");
    *GprintString += int(Avg_G+0.5);
    GprintString->Append(" #pm ");
    *GprintString += Avg_G_stat;
    GprintString->Append(" #pm ");
    *GprintString += Avg_G_syst;
    TLatex *Gprint = new TLatex(0.03, 90, GprintString->Data());
    Gprint->SetTextFont(TextFont);
    Gprint->SetTextSize(SizeSpecif);
    Gprint->Draw("same");
   
    //TString *Gname=new TString("GversusQ4_K1_M");
    //*Gname += MBOI;
    //Gname->Append(".eps");
    //can5->SaveAs(Gname->Data());

    ///////////////////////////////////////////////////////////////////////////////////////
    // G versus Centrality
    /*
    TCanvas *can6 = new TCanvas("can6", "can6",600,700,700,600);// 11,53,700,500
    can6->SetHighLightColor(2);
    gStyle->SetOptFit(0111);
    can6->SetFillColor(0);//10
    can6->SetBorderMode(0);
    can6->SetBorderSize(2);
    can6->SetFrameFillColor(0);
    can6->SetFrameBorderMode(0);
    can6->SetFrameBorderMode(0);
    
    TPad *pad6 = new TPad("pad6","pad6",0.0,0.0,1.,1.);
    gPad->SetTickx();
    gPad->SetTicky();
    pad6->SetTopMargin(0.0);//0.05
    pad6->SetRightMargin(0.0);//1e-2
    pad6->SetBottomMargin(0.0);//0.12
    pad6->Draw();
    pad6->cd(1);
    gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.03); gPad->SetBottomMargin(0.14);

    // K0
    double G_Cent[5]={28, 28, 28, 16, 24};//30, 29, 27, 27, 29, 29, 29, 25, 29, 22
    double G_stat_Cent[5]={2, 2, 14, 5, 6};//2, 2, 2, 3, 3, 3, 5, 6, 6, 8
    double G_syst_Cent[5]={8, 8, 11, 15, 10};//7, 7, 12, 16, 7, 8, 10, 6, 7, 19
    // K1
    //double G_Cent[5]={0};
    //double G_stat_Cent[5]={0};
    //double G_syst_Cent[5]={0};
    //
    double xpoints[5]={0,1,2,3,4};
    TGraph *Stat_trend=new TGraph(5,xpoints,G_stat_Cent);
    TGraph *Syst_trend=new TGraph(5,xpoints,G_syst_Cent);
    TF1 *stat_fit=new TF1("stat_fit","[0]+[1]*x",0,100);
    TF1 *syst_fit=new TF1("syst_fit","[0]+[1]*x",0,100);
    stat_fit->SetParameter(0,5); stat_fit->SetParLimits(0,1,10);
    syst_fit->SetParameter(0,5); syst_fit->SetParLimits(0,1,10);
    stat_fit->SetParameter(1,.05); stat_fit->SetParLimits(1,0,1.);
    syst_fit->SetParameter(1,.05); syst_fit->SetParLimits(1,0,1.);
    
    Stat_trend->Fit(stat_fit,"IMENQ","",0,5);
    Syst_trend->Fit(syst_fit,"IMENQ","",0,5);

    double CentEdges[6]={0,5,10,20,35,50};
    TH1D *GversusCent_Stat=new TH1D("GversusCent_Stat","",5,CentEdges);
    TH1D *GversusCent_Syst=new TH1D("GversusCent_Syst","",5,CentEdges);
    TH1D *GversusCent_StatSyst=new TH1D("GversusCent_StatSyst","",5,CentEdges);
    for(int i=0; i<5; i++){
      GversusCent_Stat->SetBinContent(i+1, G_Cent[i]);
      GversusCent_Syst->SetBinContent(i+1, G_Cent[i]);
      GversusCent_StatSyst->SetBinContent(i+1, G_Cent[i]);
      //
      GversusCent_Stat->SetBinError(i+1, stat_fit->Eval(i));
      GversusCent_Syst->SetBinError(i+1, syst_fit->Eval(i));
      GversusCent_StatSyst->SetBinError(i+1, sqrt( pow(stat_fit->Eval(i),2) + pow(syst_fit->Eval(i),2)));
    }
    GversusCent_Stat->SetMarkerStyle(20);
    GversusCent_Stat->SetLineColor(2); GversusCent_Stat->SetMarkerColor(2);
    GversusCent_Syst->SetMarkerSize(0); GversusCent_Syst->SetMarkerColor(2); GversusCent_Syst->SetLineColor(kRed-10);
    GversusCent_Syst->SetLineWidth(10);
    GversusCent_Stat->GetXaxis()->SetTitle("Centrality (%)"); GversusCent_Stat->GetYaxis()->SetTitle("Coherent fraction (%)");
    GversusCent_Stat->GetXaxis()->SetTitleSize(SizeTitle);  GversusCent_Stat->GetYaxis()->SetTitleSize(SizeTitle);
    GversusCent_Stat->GetXaxis()->SetLabelSize(SizeLabel);  GversusCent_Stat->GetYaxis()->SetLabelSize(SizeLabel);
    GversusCent_Stat->GetXaxis()->SetNdivisions(606); GversusCent_Stat->GetYaxis()->SetNdivisions(505);
    GversusCent_Stat->SetMinimum(10); GversusCent_Stat->SetMaximum(50);

    GversusCent_Stat->Draw();
    GversusCent_Syst->Draw("same");
    GversusCent_Stat->Draw("same");
    TLatex *ALICEspecif_2=new TLatex(0.5,0.9,"ALICE Pb-Pb #\sqrt{#font[12]{s}_{NN}}=2.76 TeV");
    ALICEspecif_2->SetNDC(1);
    ALICEspecif_2->SetTextFont(TextFont);
    ALICEspecif_2->SetTextSize(SizeSpecif);
    ALICEspecif_2->Draw("same");
    TLatex *KT4specif_2=(TLatex*)KT4specif->Clone();
    KT4specif_2->SetX(0.5); KT4specif_2->SetY(0.8);
    KT4specif_2->Draw("same");
    
    TF1 *GcentFit = new TF1("GcentFit","pol0",0,50);
    GversusCent_Stat->Fit(GcentFit,"IMENQ","",0,50);
    GcentFit->Draw("same");
    */
  }
  
  
  
  
  
  
}

