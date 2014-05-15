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
#define kappa3 0.16
#define kappa4 0.40

using namespace std;


void Plot_muons(){
  
 
  gStyle->SetOptStat(0);
  gStyle->SetOptDate(0);
  gStyle->SetOptFit(0111);
  
  int ChComb=1;// 0=Same-charge, 1=Mixed-charge
  int kTbin_L=5, kTbin_H=6;
  int binKT3=1;//1-2

  //TFile *_file0= new TFile("Results/PDC_12a17a_muons_R6.root","READ");
  //TFile *_file0= new TFile("Results/PDC_13b2_efix_p1_muons_R2.root","READ");
  TFile *_file0= new TFile("Results/PDC_12a17e_muons_R4.root","READ");
  
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
  TH3D *PurityNum_3D = (TH3D*)MyList->FindObject("Explicit2_Charge1_1_Charge2_1_SC_0_M_6_ED_0_Term_1_PIDpurityNum");
  TH2D *PurityDen_2D = (TH2D*)MyList->FindObject("Explicit2_Charge1_1_Charge2_1_SC_0_M_6_ED_0_Term_1_PIDpurityDen");
  TH1D *PurityNum=PurityNum_3D->ProjectionX("PurityNum",kTbin_L,kTbin_H,1,20);
  double PurityNorm=PurityDen_2D->Integral(kTbin_L,kTbin_H,1,20);
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
  TH3D *MuonSmearedNum2_3=(TH3D*)MyList->FindObject("fMuonContamSmearedNum2");
  TH3D *MuonSmearedDen2_3=(TH3D*)MyList->FindObject("fMuonContamSmearedDen2");
  TH3D *PionNum2_3=(TH3D*)MyList->FindObject("fMuonContamIdealNum2");
  TH3D *PionDen2_3=(TH3D*)MyList->FindObject("fMuonContamIdealDen2");
  TH3D *PionPionK2_3=(TH3D*)MyList->FindObject("fPionPionK2");
  //
  TH3D *MuonSmearedNum3_3=(TH3D*)MyList->FindObject("fMuonContamSmearedNum3");
  TH3D *MuonSmearedDen3_3=(TH3D*)MyList->FindObject("fMuonContamSmearedDen3");
  TH3D *PionNum3_3=(TH3D*)MyList->FindObject("fMuonContamIdealNum3");
  TH3D *PionDen3_3=(TH3D*)MyList->FindObject("fMuonContamIdealDen3");
  TH3D *PionPionK3_3=(TH3D*)MyList->FindObject("fPionPionK3");
  TH3D *MuonPionK3_3=(TH3D*)MyList->FindObject("fMuonPionK3");
  
  
  TH1D *MuonSmearedNum2[2];// SC/MC  
  TH1D *MuonSmearedDen2[2];
  TH1D *PionNum2[2];
  TH1D *PionDen2[2];
  TH1D *PionPionK2[2];
  //
  TH1D *MuonSmearedNum3[2];// SC/MC  
  TH1D *MuonSmearedDen3[2];
  TH1D *PionNum3[2];
  TH1D *PionDen3[2];
  TH1D *PionPionK3[2];
  //
  TH1D *C2muonSmeared[2];
  TH1D *C2pion[2];
  TH1D *C3muonSmeared[2];
  TH1D *C3pion[2];
  //
  for(int chtype=0; chtype<2; chtype++){
    TString *names[10];
    for(int i=0; i<10; i++) {names[i]=new TString("name_"); *names[i] += i; *names[i] += chtype;}

    MuonSmearedNum2[chtype]=(TH1D*)MuonSmearedNum2_3->ProjectionZ(names[0]->Data(),chtype+1,chtype+1,kTbin_L,kTbin_H);
    MuonSmearedDen2[chtype]=(TH1D*)MuonSmearedDen2_3->ProjectionZ(names[1]->Data(),chtype+1,chtype+1,kTbin_L,kTbin_H);
    PionNum2[chtype]=(TH1D*)PionNum2_3->ProjectionZ(names[2]->Data(),chtype+1,chtype+1,kTbin_L,kTbin_H);
    PionDen2[chtype]=(TH1D*)PionDen2_3->ProjectionZ(names[3]->Data(),chtype+1,chtype+1,kTbin_L,kTbin_H);
    PionPionK2[chtype]=(TH1D*)PionPionK2_3->ProjectionZ(names[4]->Data(),chtype+1,chtype+1,kTbin_L,kTbin_H);
    PionPionK2[chtype]->Divide(PionDen2[chtype]);
    ////////////////
    MuonSmearedNum3[chtype]=(TH1D*)MuonSmearedNum3_3->ProjectionZ(names[5]->Data(),chtype+1,chtype+1,binKT3,binKT3);
    MuonSmearedDen3[chtype]=(TH1D*)MuonSmearedDen3_3->ProjectionZ(names[6]->Data(),chtype+1,chtype+1,binKT3,binKT3);
    PionNum3[chtype]=(TH1D*)PionNum3_3->ProjectionZ(names[7]->Data(),chtype+1,chtype+1,binKT3,binKT3);
    PionDen3[chtype]=(TH1D*)PionDen3_3->ProjectionZ(names[8]->Data(),chtype+1,chtype+1,binKT3,binKT3);
    PionPionK3[chtype]=(TH1D*)PionPionK3_3->ProjectionZ(names[9]->Data(),chtype+1,chtype+1,binKT3,binKT3);
    PionPionK3[chtype]->Divide(PionDen3[chtype]);
    //
    C2muonSmeared[chtype]=(TH1D*)MuonSmearedNum2[chtype]->Clone();
    C2pion[chtype]=(TH1D*)PionNum2[chtype]->Clone();
    C2muonSmeared[chtype]->Divide(MuonSmearedDen2[chtype]);
    C2pion[chtype]->Divide(PionDen2[chtype]);
    //
    C3muonSmeared[chtype]=(TH1D*)MuonSmearedNum3[chtype]->Clone();
    C3pion[chtype]=(TH1D*)PionNum3[chtype]->Clone();
    C3muonSmeared[chtype]->Divide(MuonSmearedDen3[chtype]);
    C3pion[chtype]->Divide(PionDen3[chtype]);
    //
    //
    C2pion[chtype]->SetLineColor(4);
    C2muonSmeared[chtype]->SetLineColor(2);
    C2pion[chtype]->GetXaxis()->SetRangeUser(0,0.15);
    C2pion[chtype]->SetMinimum(0.98); C2pion[chtype]->SetMaximum(1.35);
    C2pion[chtype]->GetYaxis()->SetTitleOffset(1.5);
    C2pion[chtype]->GetXaxis()->SetTitle("q_{inv} (GeV/c)");
    C2pion[chtype]->GetYaxis()->SetTitle("C_{2}K_{2}");
    //
    C3pion[chtype]->SetLineColor(4);
    C3muonSmeared[chtype]->SetLineColor(2);
    C3pion[chtype]->GetXaxis()->SetRangeUser(0,0.15);
    C3pion[chtype]->SetMinimum(0.90); C3pion[chtype]->SetMaximum(2.0);
    C3pion[chtype]->GetYaxis()->SetTitleOffset(1.5);
    C3pion[chtype]->GetXaxis()->SetTitle("Q_{3} (GeV/c)");
    C3pion[chtype]->GetYaxis()->SetTitle("C_{3}K_{3}");
  }
  //
  //
  
  C2pion[ChComb]->Draw();
  C2muonSmeared[ChComb]->Draw("same");
  legend->AddEntry(C2pion[ChComb],"Input Pion-Pion, C_{2}^{#pi-#pi,QS}K_{2}^{#pi-#pi,QS}","l");
  legend->AddEntry(C2muonSmeared[ChComb],"Pion-Muon residual, C_{2}^{#mu-#pi,QS}K_{2}^{#mu-#pi,QS}","l");
  legend->Draw("same");
  //PionPionK2->Draw("same");
  //
  //
  
  //
  //C3pion[ChComb]->Draw();
  //C3muonSmeared[ChComb]->Draw("same");
  //PionPionK3[0]->Draw("same");
  //legend->AddEntry(C3pion[ChComb],"Input Pion-Pion-Pion, C_{3}^{#pi-#pi-#pi,QS}K_{3}^{#pi-#pi-#pi,QS}","l");
  //legend->AddEntry(C3muonSmeared[ChComb],"Muon-Pion-Pion residual, C_{3}^{#mu-#pi-#pi,QS}K_{3}^{#mu-#pi-#pi,QS}","l");
  //legend->Draw("same");
    
  
  // corrections
  TFile *fout=new TFile("MuonCorrection_temp.root","RECREATE");
  TH1D *C2muonCorrection[2]; 
  C2muonCorrection[0] = new TH1D("C2muonCorrection_SC","",100,0,0.5);
  C2muonCorrection[1] = new TH1D("C2muonCorrection_MC","",100,0,0.5);
  TH1D *WeightmuonCorrection = new TH1D("WeightmuonCorrection","",100,0,0.5);
  TH1D *C3muonCorrection[2]; 
  C3muonCorrection[0] = new TH1D("C3muonCorrection_SC","",50,0,0.5);
  C3muonCorrection[1] = new TH1D("C3muonCorrection_MC","",50,0,0.5);
  //
  C2muonCorrection[0]->GetXaxis()->SetTitle("q_{inv} (GeV/c)"); C2muonCorrection[0]->GetYaxis()->SetTitle("x_{2}"); 
  C2muonCorrection[1]->GetXaxis()->SetTitle("q_{inv} (GeV/c)"); C2muonCorrection[1]->GetYaxis()->SetTitle("x_{2}"); 
  C3muonCorrection[0]->GetXaxis()->SetTitle("Q_{3} (GeV/c)"); C3muonCorrection[0]->GetYaxis()->SetTitle("x_{3}"); 
  C3muonCorrection[1]->GetXaxis()->SetTitle("Q_{3} (GeV/c)"); C3muonCorrection[1]->GetYaxis()->SetTitle("x_{3}"); 
  WeightmuonCorrection->GetXaxis()->SetTitle("q_{inv} (GeV/c)"); WeightmuonCorrection->GetYaxis()->SetTitle("x_{2}^{w}"); 
  
  // 0.944 and 0.959
  float GoodPairFraction=1-0.93*(1-PurityNum->GetBinContent(10));
  cout<<"Pion Pair Purity = "<<PurityNum->GetBinContent(10)<<endl;
  cout<<"Effective Pion Pair Purity = "<<GoodPairFraction<<endl;
  float pionPurity=pow(GoodPairFraction,0.5);
  float muonPurity=1-pionPurity;
  for(int chtype=0; chtype<2; chtype++){
    for(int bin=1; bin<=100; bin++){
      
      bool emptybin2=kFALSE, emptybin3=kFALSE;
      if(PionPionK2[chtype]->GetBinContent(bin)==0) {PionPionK2[chtype]->SetBinContent(bin, 1.00001);}
      if(PionPionK3[chtype]->GetBinContent(bin)==0) {PionPionK3[chtype]->SetBinContent(bin, 1.00001);}
      if(bin > C2pion[chtype]->GetNbinsX()) emptybin2=kTRUE;
      if(bin > C3pion[chtype]->GetNbinsX()) emptybin3=kTRUE;
 
      double value = C2pion[chtype]->GetBinContent(bin)/PionPionK2[chtype]->GetBinContent(bin);
      double den = (GoodPairFraction*C2pion[chtype]->GetBinContent(bin) + (1-GoodPairFraction)*C2muonSmeared[chtype]->GetBinContent(bin))/PionPionK2[chtype]->GetBinContent(bin);
      if(den > 0 && !emptybin2) C2muonCorrection[chtype]->SetBinContent(bin,value/den);
      else C2muonCorrection[chtype]->SetBinContent(bin, 1);
      //
      if(chtype==0){
	value = C2pion[chtype]->GetBinContent(bin)/PionPionK2[chtype]->GetBinContent(bin) - 1.0;
	den = ((GoodPairFraction*C2pion[chtype]->GetBinContent(bin) + (1-GoodPairFraction)*C2muonSmeared[chtype]->GetBinContent(bin))/PionPionK2[chtype]->GetBinContent(bin)) - 1.0;
	if(den > 0 && !emptybin2) WeightmuonCorrection->SetBinContent(bin,value/den);
      }
      //
      value = C3pion[chtype]->GetBinContent(bin)/PionPionK3[chtype]->GetBinContent(bin);
      den = (pow(pionPurity,3)*C3pion[chtype]->GetBinContent(bin) + (3*pow(pionPurity,2)*muonPurity)*C3muonSmeared[chtype]->GetBinContent(bin))/PionPionK3[chtype]->GetBinContent(bin);
      if(den > 0 && !emptybin3) C3muonCorrection[chtype]->SetBinContent(bin,value/den);
      else C3muonCorrection[chtype]->SetBinContent(bin, 1);
    }
    //C2muonCorrection[chtype]->SetBinContent(1, C2muonCorrection[chtype]->GetBinContent(2));
    //C3muonCorrection[chtype]->SetBinContent(1, C3muonCorrection[chtype]->GetBinContent(2));
    C2muonCorrection[chtype]->Write();
    C3muonCorrection[chtype]->Write();
    if(chtype==0) {
      //WeightmuonCorrection->SetBinContent(1, WeightmuonCorrection->GetBinContent(2));
      WeightmuonCorrection->Write();
    }
  }
  
  //
  //C3muonCorrection[0]->SetMinimum(0.99); 
  //C3muonCorrection[0]->SetMaximum(1.05); 
  //C3muonCorrection[0]->GetYaxis()->SetTitleOffset(1.3);
  //C3muonCorrection[0]->Draw();
  //WeightmuonCorrection->GetYaxis()->SetTitleOffset(1.3);
  //WeightmuonCorrection->Draw();
  fout->Close();
  
  
}
