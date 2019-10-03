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
#define kappa3 0.1
#define kappa4 0.5

using namespace std;


void MuonCorrectionFourPion(){
  
  gStyle->SetOptStat(0);
  gStyle->SetOptDate(0);
  gStyle->SetOptFit(0111);
  
  bool UNITYweights=0;

  bool PbPbcase=0;

  int COI_23=0;// 0=Same-charge, 1=Mixed-charge
  int COI_4=0;// 0=Same-charge, 1=(---+), 2=(--++)
  //
  int TOI_3=1;// Term Of Interest, 3-particle (1-4)
  int TOI_4=11;// Term Of Interest, 4-particle (1-12)
  int Rbin=11;// Radius bin (1-11 fm, 11 bins)

  int ThreeParticleRebin=2;// 2, 3 for a variation
  int FourParticleRebin=3;// 3, 4 for a variation
  if(!PbPbcase) FourParticleRebin=6;

  TFile *_file0;
  //TFile *_file0= new TFile("MyOutput.root","READ");
  //TFile *_file0= new TFile("Results/PDC_12a17a_noTTC.root","READ");
  //TFile *_file0= new TFile("Results/PDC_12a17a_TTCweights.root","READ");
  //TFile *_file0= new TFile("Results/PDC_12a17a_pTSpectrumWeight.root","READ");
  //TFile *_file0= new TFile("Results/PDC_12a17a_10MeVcut.root","READ");
  if(PbPbcase) _file0= new TFile("Results/PDC_12a17a_Qweights.root","READ");
  else _file0= new TFile("Results/PDC_13b2_p1234.root","READ");


  TList *MyList=(TList*)_file0->Get("MyList");
  _file0->Close();

  TH1D *PionCandidates=(TH1D*)MyList->FindObject("fPionCandidates");
  PionCandidates->GetXaxis()->SetTitle("PDG code");
  //PionCandidates->Draw();
  //
  TH1D *MuonParentsPrimary=(TH1D*)MyList->FindObject("fMuonParents");
  MuonParentsPrimary->GetXaxis()->SetTitle("PDG code");
  MuonParentsPrimary->SetFillColor(1);
  //MuonParentsPrimary->Draw();
  //
  TH1D *MuonParentsSecondary=(TH1D*)MyList->FindObject("fSecondaryMuonParents");
  MuonParentsSecondary->GetXaxis()->SetTitle("PDG code");
  MuonParentsSecondary->SetFillColor(1);
  //MuonParentsSecondary->Draw();
  //
  // M0 R10-R6, M6 for R4, M17 for R2
  TH3D *PurityNum_3D = (TH3D*)MyList->FindObject("TwoParticle_Charge1_1_Charge2_1_M_0_ED_0_Term_2_PIDpurityNum");
  TH2D *PurityDen_2D = (TH2D*)MyList->FindObject("TwoParticle_Charge1_1_Charge2_1_M_0_ED_0_Term_2_PIDpurityDen");
  TH1D *PurityNum=PurityNum_3D->ProjectionX("PurityNum",5,6,1,20);
  double PurityNorm=PurityDen_2D->Integral(5,6,1,20);
  PurityNum->Scale(1/PurityNorm);
  char *namesAxis[15]={"e-e","e-mu","e-pi","e-k","e-p","mu-mu","mu-pi","mu-k","mu-p","pi-pi","pi-k","pi-p","k-k","k-p","p-p"};
  for(int i=1; i<=15; i++) PurityNum->GetXaxis()->SetBinLabel(i, namesAxis[i-1]);
  PurityNum->GetXaxis()->SetRange(1,15);
  PurityNum->GetYaxis()->SetTitle("Probability");
  PurityNum->Draw();
  //
  //
  TCanvas *can = new TCanvas("can", "can",800,0,800,800);// 11,53,700,500
  can->SetHighLightColor(2);
  gStyle->SetOptFit(0111);
  can->SetFillColor(10);//10
  can->SetBorderMode(0);
  can->SetBorderSize(2);
  can->SetFrameFillColor(0);
  can->SetFrameBorderMode(0);
  can->SetFrameBorderMode(0);
  can->cd();
  TPad *pad = new TPad("pad","pad",0.,0.,1.,1.);
  gPad->SetTickx();
  gPad->SetTicky();
  pad->SetGridx();
  pad->SetGridy();
  pad->SetTopMargin(0.02);//0.05
  pad->SetRightMargin(0.02);//3e-2
  pad->SetBottomMargin(0.1);//0.12
  pad->SetLeftMargin(0.1);
  pad->Draw();
  pad->cd();
  TLegend *legend = new TLegend(.5,.65, .9,.95,NULL,"brNDC");//.45 or .4 for x1
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.03);
  //

  TH1D *C2muonSmeared;
  TH1D *C2pion;
  TH1D *K2pion;
  TH1D *C3muonSmeared;
  TH1D *C3pion;
  TH1D *C4muonSmeared;
  TH1D *C4pion;
  //
  TH2D *MuonSmeared2[2][2];// ChComb, Term
  TH2D *MuonIdeal2[2][2];
  TH2D *PionPionK2[2][2];
  //
  TH2D *MuonSmeared3[2][4];// ChComb, Term (4 terms only, N1^3 in other dimension)
  TH2D *MuonIdeal3[2][4];
  TH2D *PionPionK3[2][4];
  //
  TH2D *MuonSmeared4[3][12];// ChComb, Term (12 terms only, N1^3 in other dimension)
  TH2D *MuonIdeal4[3][12];
  TH2D *PionPionK4[3][12];

  for(int ch=0; ch<=2; ch++){// same-charge, mixed-charge type 1 then type 2 (4-pion only)
    
    if(ch<=1){
      for(int term=1; term<=2; term++){
	TString *MS2 = new TString("TwoParticle_Charge1_");
	if(ch==0) MS2->Append("1_Charge2_1_M_0_ED_0_Term_");
	else MS2->Append("0_Charge2_1_M_0_ED_0_Term_");
	*MS2 += term;
	TString *MI2 = new TString(MS2->Data());
	TString *K2 = new TString(MS2->Data());
	MS2->Append("_MuonSmeared"); MI2->Append("_MuonIdeal"); K2->Append("_PionPionK2");
	MuonSmeared2[ch][term-1]=(TH2D*)MyList->FindObject(MS2->Data());
	MuonIdeal2[ch][term-1]=(TH2D*)MyList->FindObject(MI2->Data());
	PionPionK2[ch][term-1]=(TH2D*)MyList->FindObject(K2->Data());
	//MuonSmeared2[ch][term-1]->RebinY(2);
	//MuonIdeal2[ch][term-1]->RebinY(2);
	//PionPionK2[ch][term-1]->RebinY(2);
      }
      MuonSmeared2[ch][0]->Divide(MuonSmeared2[ch][1]);
      MuonIdeal2[ch][0]->Divide(MuonIdeal2[ch][1]);
      PionPionK2[ch][0]->Divide(PionPionK2[ch][1]);
      //
      if(ch==COI_23){
	C2muonSmeared = (TH1D*)MuonSmeared2[ch][0]->ProjectionY("C2muonSmeared",Rbin,Rbin);
	C2pion = (TH1D*)MuonIdeal2[ch][0]->ProjectionY("C2muonIdeal",Rbin,Rbin);
	K2pion = (TH1D*)PionPionK2[ch][0]->ProjectionY("K2pion",Rbin,Rbin);
      }
    }
    
    //
    if(ch<=1){
      for(int term=1; term<=4; term++){
	TString *MS3 = new TString("ThreeParticle_Charge1_");
	if(ch==0) MS3->Append("1_Charge2_1_Charge3_1_M_0_ED_0_Term_");
	else MS3->Append("0_Charge2_1_Charge3_1_M_0_ED_0_Term_");
	*MS3 += term;
	TString *MI3 = new TString(MS3->Data());
	TString *K3 = new TString(MS3->Data());
	MS3->Append("_MuonSmeared"); MI3->Append("_MuonIdeal"); K3->Append("_PionPionK3");
	
	TH3D *TempMuonSmeared=(TH3D*)MyList->FindObject(MS3->Data());
	TH3D *TempMuonIdeal=(TH3D*)MyList->FindObject(MI3->Data());
	TH3D *TempPionPionK3=(TH3D*)MyList->FindObject(K3->Data());
	TempMuonSmeared->RebinZ(ThreeParticleRebin);
	TempMuonIdeal->RebinZ(ThreeParticleRebin);
	TempPionPionK3->RebinZ(ThreeParticleRebin);
	//
	TempMuonSmeared->GetXaxis()->SetRange(1,1);
	TempMuonIdeal->GetXaxis()->SetRange(1,1);
	TempPionPionK3->GetXaxis()->SetRange(1,1);
	MuonSmeared3[ch][term-1]=(TH2D*)TempMuonSmeared->Project3D("zy");
	MuonIdeal3[ch][term-1]=(TH2D*)TempMuonIdeal->Project3D("zy");
	PionPionK3[ch][term-1]=(TH2D*)TempPionPionK3->Project3D("zy");
	TString *name1 = new TString("ProNameSmear3"); *name1 += ch; *name1 += term;
	TString *name2 = new TString("ProNameIdeal3"); *name2 += ch; *name2 += term;
	TString *name3 = new TString("ProNamePionK3"); *name3 += ch; *name3 += term;
	MuonSmeared3[ch][term-1]->SetName(name1->Data());
	MuonIdeal3[ch][term-1]->SetName(name2->Data());
	PionPionK3[ch][term-1]->SetName(name3->Data());
	//
	TempMuonSmeared->GetXaxis()->SetRange(2,2);
	TempMuonIdeal->GetXaxis()->SetRange(2,2);
	TempPionPionK3->GetXaxis()->SetRange(2,2);
	TH2D *tempDen1 = (TH2D*)TempMuonSmeared->Project3D("zy");
	TH2D *tempDen2 = (TH2D*)TempMuonIdeal->Project3D("zy");
	TH2D *tempDen3 = (TH2D*)TempPionPionK3->Project3D("zy");
	TString *nameDen1 = new TString("ProNameSmearDen3"); *nameDen1 += ch; *nameDen1 += term;
	TString *nameDen2 = new TString("ProNameIdealDen3"); *nameDen2 += ch; *nameDen2 += term;
	TString *nameDen3 = new TString("ProNamePionDenK3"); *nameDen3 += ch; *nameDen3 += term;
	tempDen1->SetName(nameDen1->Data());
	tempDen2->SetName(nameDen2->Data());
	tempDen3->SetName(nameDen3->Data());
	//
	MuonSmeared3[ch][term-1]->Divide(tempDen1);
	MuonIdeal3[ch][term-1]->Divide(tempDen2);
	PionPionK3[ch][term-1]->Divide(tempDen3);
	if(ch==COI_23 && term==TOI_3){
	  C3muonSmeared = (TH1D*)MuonSmeared3[ch][term-1]->ProjectionY("C3muonSmeared",Rbin,Rbin);
	  C3pion = (TH1D*)MuonIdeal3[ch][term-1]->ProjectionY("C3muonIdeal",Rbin,Rbin);
	}
      }
      
    }
    
    //
    for(int term=1; term<=12; term++){
      TString *MS4 = new TString("FourParticle_Charge1_");
      if(ch==0) MS4->Append("1_Charge2_1_Charge3_1_Charge4_1_M_0_ED_0_Term_");
      else if(ch==1) MS4->Append("0_Charge2_1_Charge3_1_Charge4_1_M_0_ED_0_Term_");
      else MS4->Append("0_Charge2_0_Charge3_1_Charge4_1_M_0_ED_0_Term_");
      *MS4 += term;
      TString *MI4 = new TString(MS4->Data());
      TString *K4 = new TString(MS4->Data());
      MS4->Append("_MuonSmeared"); MI4->Append("_MuonIdeal"); K4->Append("_PionPionK4");
      TH3D *TempMuonSmeared=(TH3D*)MyList->FindObject(MS4->Data());
      TH3D *TempMuonIdeal=(TH3D*)MyList->FindObject(MI4->Data());
      TH3D *TempPionPionK4=(TH3D*)MyList->FindObject(K4->Data());
      TempMuonSmeared->RebinZ(FourParticleRebin);
      TempMuonIdeal->RebinZ(FourParticleRebin);
      TempPionPionK4->RebinZ(FourParticleRebin);
      //
      TempMuonSmeared->GetXaxis()->SetRange(1,1);
      TempMuonIdeal->GetXaxis()->SetRange(1,1);
      TempPionPionK4->GetXaxis()->SetRange(1,1);
      MuonSmeared4[ch][term-1]=(TH2D*)TempMuonSmeared->Project3D("zy");
      MuonIdeal4[ch][term-1]=(TH2D*)TempMuonIdeal->Project3D("zy");
      PionPionK4[ch][term-1]=(TH2D*)TempPionPionK4->Project3D("zy");
      TString *name1 = new TString("ProNameSmear4"); *name1 += ch; *name1 += term;
      TString *name2 = new TString("ProNameIdeal4"); *name2 += ch; *name2 += term;
      TString *name3 = new TString("ProNamePionK4"); *name3 += ch; *name3 += term;
      MuonSmeared4[ch][term-1]->SetName(name1->Data());
      MuonIdeal4[ch][term-1]->SetName(name2->Data());
      PionPionK4[ch][term-1]->SetName(name3->Data());
      //
      TempMuonSmeared->GetXaxis()->SetRange(2,2);
      TempMuonIdeal->GetXaxis()->SetRange(2,2);
      TempPionPionK4->GetXaxis()->SetRange(2,2);
      TH2D *tempDen1 = (TH2D*)TempMuonSmeared->Project3D("zy");
      TH2D *tempDen2 = (TH2D*)TempMuonIdeal->Project3D("zy");
      TH2D *tempDen3 = (TH2D*)TempPionPionK4->Project3D("zy");
      TString *nameDen1 = new TString("ProNameSmearDen3"); *nameDen1 += ch; *nameDen1 += term;
      TString *nameDen2 = new TString("ProNameIdealDen3"); *nameDen2 += ch; *nameDen2 += term;
      TString *nameDen3 = new TString("ProNamePionDenK4"); *nameDen3 += ch; *nameDen3 += term;
      tempDen1->SetName(nameDen1->Data());
      tempDen2->SetName(nameDen2->Data());
      tempDen3->SetName(nameDen3->Data());
      //
      MuonSmeared4[ch][term-1]->Divide(tempDen1);
      MuonIdeal4[ch][term-1]->Divide(tempDen2);
      PionPionK4[ch][term-1]->Divide(tempDen3);
      if(ch==COI_4 && term==TOI_4){
	C4muonSmeared = (TH1D*)MuonSmeared4[ch][term-1]->ProjectionY("C4muonSmeared",Rbin,Rbin);
	C4pion = (TH1D*)MuonIdeal4[ch][term-1]->ProjectionY("C4muonIdeal",Rbin,Rbin);
      }

    }


  }
    
  
  C2pion->SetLineColor(4);
  C2muonSmeared->SetLineColor(2);
  C2pion->GetXaxis()->SetRangeUser(0,0.15);
  C2pion->SetMinimum(0.98); C2pion->SetMaximum(1.35);
  C2pion->GetYaxis()->SetTitleOffset(1.5);
  C2pion->GetXaxis()->SetTitle("q_{inv} (GeV/c)");
  C2pion->GetYaxis()->SetTitle("C_{2}K_{2}");
  //
  C3pion->SetLineColor(4);
  C3muonSmeared->SetLineColor(2);
  C3pion->GetXaxis()->SetRangeUser(0,0.15);
  C3pion->SetMinimum(0.98); C3pion->SetMaximum(2.5);
  C3pion->GetYaxis()->SetTitleOffset(1.5);
  C3pion->GetXaxis()->SetTitle("Q_{3} (GeV/c)");
  C3pion->GetYaxis()->SetTitle("C_{3}K_{3}");
  //
  C4pion->SetLineColor(4);
  C4muonSmeared->SetLineColor(2);
  C4pion->GetXaxis()->SetRangeUser(0,0.15);
  C4pion->SetMinimum(0.98); C4pion->SetMaximum(4.0);
  C4pion->GetYaxis()->SetTitleOffset(1.5);
  C4pion->GetXaxis()->SetTitle("Q_{4} (GeV/c)");
  C4pion->GetYaxis()->SetTitle("C_{4}K_{4}");
  //
  //C2pion->Draw();
  //C2muonSmeared->Draw("same");
  //K2pion->Draw("same");
  //
  //C3pion->Draw();
  //C3muonSmeared->Draw("same");
  //
  C4pion->SetMaximum(1.5);
  C4pion->Draw();
  C4muonSmeared->Draw("same");

  
  // corrections
  float GoodPairFraction=1-0.93*(1-PurityNum->GetBinContent(10));
  //GoodPairFraction = 0.92;// 0.92 and 0.96 as a variation
  cout<<"Pion Pair Purity = "<<PurityNum->GetBinContent(10)<<endl;
  cout<<"Effective Pion Pair Purity = "<<GoodPairFraction<<endl;
  float pionPurity=pow(GoodPairFraction,0.5);
  float muonPurity=1-pionPurity;
  //
  TFile *fout=new TFile("MuonCorrection_temp.root","RECREATE");
  //
  TH2D *C2muonCorrection[2];// ChComb 
  C2muonCorrection[0] = new TH2D("C2muonCorrection_SC","", MuonIdeal2[0][0]->GetNbinsX(),MuonIdeal2[0][0]->GetXaxis()->GetBinLowEdge(1),MuonIdeal2[0][0]->GetXaxis()->GetBinUpEdge(MuonIdeal2[0][0]->GetNbinsX()),  MuonIdeal2[0][0]->GetNbinsY(),MuonIdeal2[0][0]->GetYaxis()->GetBinLowEdge(1),MuonIdeal2[0][0]->GetYaxis()->GetBinUpEdge(MuonIdeal2[0][0]->GetNbinsY()));
  C2muonCorrection[1] = new TH2D("C2muonCorrection_MC","",MuonIdeal2[0][0]->GetNbinsX(),MuonIdeal2[0][0]->GetXaxis()->GetBinLowEdge(1),MuonIdeal2[0][0]->GetXaxis()->GetBinUpEdge(MuonIdeal2[0][0]->GetNbinsX()),  MuonIdeal2[0][0]->GetNbinsY(),MuonIdeal2[0][0]->GetYaxis()->GetBinLowEdge(1),MuonIdeal2[0][0]->GetYaxis()->GetBinUpEdge(MuonIdeal2[0][0]->GetNbinsY()));
  TH2D *WeightmuonCorrection = new TH2D("WeightmuonCorrection","",MuonIdeal2[0][0]->GetNbinsX(),MuonIdeal2[0][0]->GetXaxis()->GetBinLowEdge(1),MuonIdeal2[0][0]->GetXaxis()->GetBinUpEdge(MuonIdeal2[0][0]->GetNbinsX()),  MuonIdeal2[0][0]->GetNbinsY(),MuonIdeal2[0][0]->GetYaxis()->GetBinLowEdge(1),MuonIdeal2[0][0]->GetYaxis()->GetBinUpEdge(MuonIdeal2[0][0]->GetNbinsY()));
  C2muonCorrection[0]->GetXaxis()->SetTitle("Radius (fm)"); C2muonCorrection[0]->GetYaxis()->SetTitle("q_{inv} (GeV/c)"); C2muonCorrection[0]->GetZaxis()->SetTitle("x_{2}"); 
  C2muonCorrection[1]->GetXaxis()->SetTitle("Radius (fm)"); C2muonCorrection[1]->GetYaxis()->SetTitle("q_{inv} (GeV/c)"); C2muonCorrection[1]->GetZaxis()->SetTitle("x_{2}"); 
  WeightmuonCorrection->GetXaxis()->SetTitle("Radius (fm)"); WeightmuonCorrection->GetYaxis()->SetTitle("q_{inv} (GeV/c)"); WeightmuonCorrection->GetZaxis()->SetTitle("x_{2}^{w}");
  //
  TH2D *C3muonCorrection[2][3];// ChComb, terms
  TH2D *C4muonCorrection[3][5];// ChComb, terms
  
  for(int ch=0; ch<3; ch++){
    
    // 2-pion
    if(ch<=1){
      for(int Rb=1; Rb<=MuonIdeal2[ch][0]->GetNbinsX(); Rb++){
	for(int bin=1; bin<=MuonIdeal2[ch][0]->GetNbinsY(); bin++){
	  bool emptybin2=kFALSE;
	  if(PionPionK2[ch][0]->GetBinContent(Rb, bin)==0) {PionPionK2[ch][0]->SetBinContent(Rb, bin, 1.00001);}
	  if(bin > MuonIdeal2[ch][0]->GetNbinsY()) emptybin2=kTRUE;
	  //
	  double value = MuonIdeal2[ch][0]->GetBinContent(Rb, bin)/PionPionK2[ch][0]->GetBinContent(Rb, bin);
	  double den = (GoodPairFraction*MuonIdeal2[ch][0]->GetBinContent(Rb, bin) + (1-GoodPairFraction)*MuonSmeared2[ch][0]->GetBinContent(Rb, bin))/PionPionK2[ch][0]->GetBinContent(Rb, bin);
	  if(den > 0 && !emptybin2) C2muonCorrection[ch]->SetBinContent(Rb, bin,value/den);
	  else C2muonCorrection[ch]->SetBinContent(Rb, bin, 1.);
	  C2muonCorrection[ch]->SetBinError(Rb, bin, 0.);
	  if(ch==0){
	    double valueW = MuonIdeal2[ch][0]->GetBinContent(Rb, bin)/PionPionK2[ch][0]->GetBinContent(Rb, bin) - 1.0;
	    double denW = ((GoodPairFraction*MuonIdeal2[ch][0]->GetBinContent(Rb, bin) + (1-GoodPairFraction)*MuonSmeared2[ch][0]->GetBinContent(Rb, bin))/PionPionK2[ch][0]->GetBinContent(Rb, bin)) - 1.0;
	    //if( valueW/denW > 10) cout<<bin<<"  "<<MuonIdeal2[ch][0]->GetBinContent(Rb, bin)<<"  "<<PionPionK2[ch][0]->GetBinContent(Rb, bin)<<"  "<<MuonSmeared2[ch][0]->GetBinContent(Rb, bin)<<"  "<<valueW<<"  "<<denW<<endl;
	    if(denW != 0 && !emptybin2) {
	      double weight = valueW/denW;
	      if(weight > 1.1) weight = 1.0;
	      if(weight <0) weight = 0;
	      WeightmuonCorrection->SetBinContent(Rb, bin, weight);
	    }	    
	    WeightmuonCorrection->SetBinError(Rb, bin, 0.);
	    //
	    if(PbPbcase && WeightmuonCorrection->GetYaxis()->GetBinCenter(bin) > 0.08) C2muonCorrection[ch]->SetBinContent(Rb, bin, 1.);
	  }
	  if(PbPbcase && C2muonCorrection[ch]->GetYaxis()->GetBinCenter(bin) > 0.08) C2muonCorrection[ch]->SetBinContent(Rb, bin, 1.);
	 
	  if(UNITYweights) {
	    C2muonCorrection[ch]->SetBinContent(Rb, bin, 1.);
	    WeightmuonCorrection->SetBinContent(Rb, bin, 1.);
	  }
	}
      }	
      if(ch==0) WeightmuonCorrection->Write();
    }

    // 3-pion
    if(ch<=1){
      for(int term=1; term<=3; term++){
	if(ch==0 && term>2) continue;
	int realterm=term;
	if(ch==0 && term==2) realterm=3;
	
	if(ch==1 && term==2) realterm=3;
	if(ch==1 && term==3) realterm=4;
	
	TString *name=new TString("C3muonCorrection_");
	if(ch==0) name->Append("SC_term");
	else name->Append("MC_term");
	//
	*name += term;
	C3muonCorrection[ch][term-1] = new TH2D(name->Data(),"",MuonIdeal3[ch][term-1]->GetNbinsX(),MuonIdeal3[ch][term-1]->GetXaxis()->GetBinLowEdge(1),MuonIdeal3[ch][term-1]->GetXaxis()->GetBinUpEdge(MuonIdeal3[ch][term-1]->GetNbinsX()),  MuonIdeal3[ch][term-1]->GetNbinsY(),MuonIdeal3[ch][term-1]->GetYaxis()->GetBinLowEdge(1),MuonIdeal3[ch][term-1]->GetYaxis()->GetBinUpEdge(MuonIdeal3[ch][term-1]->GetNbinsY()));
	C3muonCorrection[ch][term-1]->GetXaxis()->SetTitle("Radius (fm)");
	C3muonCorrection[ch][term-1]->GetYaxis()->SetTitle("Q_{3} (GeV/c)");
	C3muonCorrection[ch][term-1]->GetZaxis()->SetTitle("x_{3} (GeV/c)");
	//
	double f_pions = pow(pionPurity,3);
	double f_muonpions = 3*pow(pionPurity,2)*muonPurity;
	if(ch==1 && term==2 || term==3) {f_pions = pow(pionPurity,2); f_muonpions = 1-f_pions;}
	
	for(int Rb=1; Rb<=MuonIdeal3[ch][realterm-1]->GetNbinsX(); Rb++){
	  for(int bin=1; bin<=MuonIdeal3[ch][realterm-1]->GetNbinsY(); bin++){
	    bool emptybin3=kFALSE;
	    if(PionPionK3[ch][realterm-1]->GetBinContent(Rb, bin)==0) {PionPionK3[ch][realterm-1]->SetBinContent(Rb, bin, 1.00001);}
	    if(bin > MuonIdeal3[ch][realterm-1]->GetNbinsY()) emptybin3=kTRUE;
	    //
	    double value = MuonIdeal3[ch][realterm-1]->GetBinContent(Rb, bin)/PionPionK3[ch][realterm-1]->GetBinContent(Rb, bin);
	    double den = (f_pions*MuonIdeal3[ch][realterm-1]->GetBinContent(Rb, bin) + f_muonpions*MuonSmeared3[ch][realterm-1]->GetBinContent(Rb, bin))/PionPionK3[ch][realterm-1]->GetBinContent(Rb, bin);
	    if(den > 0 && !emptybin3) C3muonCorrection[ch][term-1]->SetBinContent(Rb, bin,value/den);
	    else C3muonCorrection[ch][term-1]->SetBinContent(Rb, bin, 1);
	    C3muonCorrection[ch][term-1]->SetBinError(Rb, bin, 0.);
	    //
	    if(PbPbcase && C3muonCorrection[ch][term-1]->GetYaxis()->GetBinCenter(bin) > 0.13) C3muonCorrection[ch][term-1]->SetBinContent(Rb, bin, 1.);
	   if(UNITYweights) C3muonCorrection[ch][term-1]->SetBinContent(Rb, bin, 1.); 
	  }
	}
	
      }
    }
  
    // 4-pion
    for(int term=1; term<=5; term++){
      int realterm=term;// realterm starts from 1
      if(ch==0 && term>4) continue;// ++++ case
      
      if(ch==0 && term==2) realterm=5;
      if(ch==0 && term==3) realterm=8;
      if(ch==0 && term==4) realterm=12;
      
      if(ch==1 && term==2) realterm=4;
      if(ch==1 && term==3) realterm=5;
      if(ch==1 && term==4) realterm=8;
      if(ch==1 && term==5) realterm=10;
      
      if(ch==2 && term==2) realterm=3;
      if(ch==2 && term==3) realterm=11;
      if(ch==2 && term==4) realterm=8;
      if(ch==2 && term==5) realterm=12;
      //
      TString *name=new TString("C4muonCorrection_");
      if(ch==0) name->Append("SC_term");
      else if(ch==1) name->Append("MC1_term");
      else name->Append("MC2_term");
      *name += term;
      
      C4muonCorrection[ch][term-1] = new TH2D(name->Data(),"",MuonIdeal4[ch][term-1]->GetNbinsX(),MuonIdeal4[ch][term-1]->GetXaxis()->GetBinLowEdge(1),MuonIdeal4[ch][term-1]->GetXaxis()->GetBinUpEdge(MuonIdeal4[ch][term-1]->GetNbinsX()),  MuonIdeal4[ch][term-1]->GetNbinsY(),MuonIdeal4[ch][term-1]->GetYaxis()->GetBinLowEdge(1),MuonIdeal4[ch][term-1]->GetYaxis()->GetBinUpEdge(MuonIdeal4[ch][term-1]->GetNbinsY()));
      C4muonCorrection[ch][term-1]->GetXaxis()->SetTitle("Radius (fm)");
      C4muonCorrection[ch][term-1]->GetYaxis()->SetTitle("Q_{4} (GeV/c)");
      C4muonCorrection[ch][term-1]->GetZaxis()->SetTitle("x_{4} (GeV/c)");
      //
      double f_pions = pow(pionPurity,4);
      double f_muonpions = 4*pow(pionPurity,3)*muonPurity;
      if(ch==0 && term==2) {f_pions = pow(pionPurity,3); f_muonpions = 3*pow(pionPurity,2)*muonPurity;}
      if(ch==0 && (term==3 || term==4) ) {f_pions = pow(pionPurity,2); f_muonpions = 1-f_pions;}
      //
      if(ch==1 && (term==2 || term==3) ) {f_pions = pow(pionPurity,3); f_muonpions = 3*pow(pionPurity,2)*muonPurity;}
      if(ch==1 && (term==4 || term==5) ) {f_pions = pow(pionPurity,2); f_muonpions = 1-f_pions;}
      //
      if(ch==2 && term==2) {f_pions = pow(pionPurity,3); f_muonpions = 3*pow(pionPurity,2)*muonPurity;}
      if(ch==2 && (term==3 || term==4 || term==5) ) {f_pions = pow(pionPurity,2); f_muonpions = 1-f_pions;}

      for(int Rb=1; Rb<=MuonIdeal4[ch][realterm-1]->GetNbinsX(); Rb++){
	for(int bin=1; bin<=MuonIdeal4[ch][realterm-1]->GetNbinsY(); bin++){
	  bool emptybin4=kFALSE;
	  if(PionPionK4[ch][realterm-1]->GetBinContent(Rb, bin)==0) {PionPionK4[ch][realterm-1]->SetBinContent(Rb, bin, 1.00001);}
	  if(bin > MuonIdeal4[ch][realterm-1]->GetNbinsY()) emptybin4=kTRUE;
	  //
	  double value = MuonIdeal4[ch][realterm-1]->GetBinContent(Rb, bin)/PionPionK4[ch][realterm-1]->GetBinContent(Rb, bin);
	  double den = (f_pions*MuonIdeal4[ch][realterm-1]->GetBinContent(Rb, bin) + f_muonpions*MuonSmeared4[ch][realterm-1]->GetBinContent(Rb, bin))/PionPionK4[ch][realterm-1]->GetBinContent(Rb, bin);
	  if(den > 0 && !emptybin4) C4muonCorrection[ch][term-1]->SetBinContent(Rb, bin,value/den);
	  else C4muonCorrection[ch][term-1]->SetBinContent(Rb, bin, 1);
	  C4muonCorrection[ch][term-1]->SetBinError(Rb, bin, 0.);
	  //
	  if(PbPbcase && C4muonCorrection[ch][term-1]->GetYaxis()->GetBinCenter(bin) > 0.19) C4muonCorrection[ch][term-1]->SetBinContent(Rb, bin, 1.);

	  if(UNITYweights) C4muonCorrection[ch][term-1]->SetBinContent(Rb, bin, 1.);
	}
      }
      
    }
    
  }// ch
  
   ///////////////////////////////////////////////////////////////

  for(int ch=0; ch<2; ch++){
    if(C2muonCorrection[ch]) C2muonCorrection[ch]->Write();
  }
  for(int ch=0; ch<2; ch++){
    for(int term=1; term<=3; term++){
      if(C3muonCorrection[ch][term-1]) C3muonCorrection[ch][term-1]->Write();
    }
  }
  for(int ch=0; ch<3; ch++){
    for(int term=1; term<=5; term++){
      if(ch==0 && term>4) continue;// ++++ case
      if(C4muonCorrection[ch][term-1]) C4muonCorrection[ch][term-1]->Write();
    }
  }
  
  fout->Close();
  
  
}
