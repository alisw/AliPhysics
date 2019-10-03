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


void MRC_FourPion(){

  gStyle->SetOptStat(0);
  gStyle->SetOptDate(0);
  gStyle->SetOptFit(0111);
  
  double MRCvariation = 1.;// factor to increase/decrease MRC (1.1 = 10% increase)
  bool PbPbcase=0;

  int COI_23=0;// 0=Same-charge, 1=Mixed-charge
  int COI_4=0;// 0=Same-charge, 1=(---+), 2=(--++)
  int TOI_3=1;// Term Of Interest, 3-particle (1-4)
  int TOI_4=1;// Term Of Interest, 4-particle (1-11)
  int Rbin=11;// Radius bin (1-11 fm, 11 bins)

  int ThreeParticleRebin=2;// 2, 3 for variation
  int FourParticleRebin=3;// 3, 4 for variation
  if(!PbPbcase) FourParticleRebin=6;
  
  TFile *_file0;
  //TFile *_file0 = new TFile("MyOutput.root","READ");
  //TFile *_file0= new TFile("Results/PDC_12a17a_noTTC.root","READ");
  //TFile *_file0= new TFile("Results/PDC_12a17a_TTC_lam0p6.root","READ");
  //TFile *_file0= new TFile("Results/PDC_12a17a_lam0p55.root","READ");
  //TFile *_file0= new TFile("Results/PDC_12a17a_TTCweights.root","READ");
  //TFile *_file0= new TFile("Results/PDC_12a17a_pTSpectrumWeight.root","READ");
  //TFile *_file0= new TFile("Results/PDC_12a17a_10MeVcut.root","READ");
  //TFile *_file0= new TFile("Results/PDC_12a17a_8MeVbins.root","READ");
  //TFile *_file0= new TFile("Results/PDC_12a17a_Qweights.root","READ");
  if(PbPbcase) _file0= new TFile("Results/PDC_12a17a_10MeVcut.root","READ");// new default
  else _file0= new TFile("Results/PDC_13b2_p1234.root","READ");

  TList *MyList=(TList*)_file0->Get("MyList");
  _file0->Close();

 
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
  pad->SetLeftMargin(0.12);
  pad->Draw();
  pad->cd();
  TLegend *legend = new TLegend(.5,.65, .9,.95,NULL,"brNDC");//.45 or .4 for x1
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.03);
  //
  
  TH1D *C2Smeared;
  TH1D *C2Ideal;
  TH1D *C3Smeared;
  TH1D *C3Ideal;
  TH1D *C4Smeared;
  TH1D *C4Ideal;
  //
  TH2D *Smeared2[2][2];
  TH2D *Ideal2[2][2];
  //
  TH2D *Smeared3[2][5];
  TH2D *Ideal3[2][5];
  //
  TH2D *Smeared4[3][13];
  TH2D *Ideal4[3][13];
  //
  //
  TFile *fout = new TFile("MomResFile_temp.root","RECREATE");
  TH2D *MRC_2[2][2];
  TH2D *MRC_3[2][5];
  TH2D *MRC_4[3][6];

  for(int ch=0; ch<=2; ch++){// same-charge, mixed-charge type 1 then type 2 (4-pion only)
    
    if(ch<=1){
      for(int term=1; term<=2; term++){
	TString *S2 = new TString("TwoParticle_Charge1_");
	if(ch==0) S2->Append("1_Charge2_1_M_0_ED_0_Term_");
	else S2->Append("0_Charge2_1_M_0_ED_0_Term_");
	*S2 += term;
	TString *I2 = new TString(S2->Data());
	S2->Append("_Smeared"); I2->Append("_Ideal");
	Smeared2[ch][term-1]=(TH2D*)MyList->FindObject(S2->Data());
	Ideal2[ch][term-1]=(TH2D*)MyList->FindObject(I2->Data());
      }
      if(ch==COI_23){
	C2Smeared = (TH1D*)Smeared2[ch][0]->ProjectionY("C2Smeared",Rbin,Rbin);
	TH1D *temp_s_2 = (TH1D*)Smeared2[ch][1]->ProjectionY("temp_s_2",Rbin,Rbin);
	C2Smeared->Divide(temp_s_2);
	C2Ideal = (TH1D*)Ideal2[ch][0]->ProjectionY("C2Ideal",Rbin,Rbin);
	TH1D *temp_i_2 = (TH1D*)Ideal2[ch][1]->ProjectionY("temp_i_2",Rbin,Rbin);
	C2Ideal->Divide(temp_i_2);
      }
      for(int term=1; term<=2; term++){
	TString *name = new TString("MRC_2_");
	if(ch==0) name->Append("SC_term");
	else name->Append("MC_term");
	//
	*name += term;
	MRC_2[ch][term-1] = new TH2D(name->Data(),"", Smeared2[0][0]->GetNbinsX(), Smeared2[0][0]->GetXaxis()->GetBinLowEdge(1), Smeared2[0][0]->GetXaxis()->GetBinUpEdge(Smeared2[0][0]->GetNbinsX()),   Smeared2[0][0]->GetNbinsY(), Smeared2[0][0]->GetYaxis()->GetBinLowEdge(1), Smeared2[0][0]->GetYaxis()->GetBinUpEdge(Smeared2[0][0]->GetNbinsY()));
	MRC_2[ch][term-1]->Add(Ideal2[ch][term-1], MRCvariation);
	MRC_2[ch][term-1]->Add(Smeared2[ch][term-1], 1-MRCvariation);
	MRC_2[ch][term-1]->Divide(Smeared2[ch][term-1]);
	//
	for(int binX=1; binX<=MRC_2[ch][term-1]->GetNbinsX(); binX++){
	  for(int binY=1; binY<=MRC_2[ch][term-1]->GetNbinsY(); binY++){
	    MRC_2[ch][term-1]->SetBinError(binX,binY, 0);
	    if(PbPbcase){
	      if(MRC_2[ch][term-1]->GetYaxis()->GetBinCenter(binY) > 0.08) MRC_2[ch][term-1]->SetBinContent(binX, binY,  MRC_2[ch][term-1]->GetBinContent(binX, MRC_2[ch][term-1]->GetYaxis()->FindBin(0.08)));
	    }
	  }
	}
      }
    }      
    //
   
    if(ch<=1){
      for(int term=1; term<=5; term++){
	TString *S3 = new TString("ThreeParticle_Charge1_");
	if(ch==0) S3->Append("0_Charge2_0_Charge3_0_M_0_ED_0_Term_");
	else S3->Append("0_Charge2_1_Charge3_1_M_0_ED_0_Term_");
	*S3 += term;
	TString *I3 = new TString(S3->Data());
	S3->Append("_Smeared"); I3->Append("_Ideal");
	Smeared3[ch][term-1]=(TH2D*)MyList->FindObject(S3->Data());
	Ideal3[ch][term-1]=(TH2D*)MyList->FindObject(I3->Data());
	//
	Smeared3[ch][term-1]->RebinY(ThreeParticleRebin);
	Ideal3[ch][term-1]->RebinY(ThreeParticleRebin);
      }
      
      for(int term=1; term<=4; term++){
	if(ch==COI_23 && term==TOI_3){
	  C3Smeared = (TH1D*)Smeared3[ch][term-1]->ProjectionY("C3Smeared",Rbin,Rbin);
	  TH1D *temp_s_3 = (TH1D*)Smeared3[ch][4]->ProjectionY("temp_s_3",Rbin,Rbin);
	  C3Smeared->Divide(temp_s_3);
	  C3Ideal = (TH1D*)Ideal3[ch][term-1]->ProjectionY("C3Ideal",Rbin,Rbin);
	  TH1D *temp_i_3 = (TH1D*)Ideal3[ch][4]->ProjectionY("temp_i_3",Rbin,Rbin);
	  C3Ideal->Divide(temp_i_3);
	}
      }
      for(int term=1; term<=4; term++){// 4 unique terms at max
	TString *name = new TString("MRC_3_");
	if(ch==0) name->Append("SC_term");
	else name->Append("MC_term");
	*name += term;
	int realterm = term;
	if(ch==0 && term>3) continue;
	
	if(ch==0 && term==2) realterm=3;
	if(ch==0 && term==3) realterm=5;
	
	if(ch==1 && term==2) realterm=3;
	if(ch==1 && term==3) realterm=4;
	if(ch==1 && term==4) realterm=5;

	MRC_3[ch][term-1] = new TH2D(name->Data(),"", Smeared3[0][0]->GetNbinsX(), Smeared3[0][0]->GetXaxis()->GetBinLowEdge(1), Smeared3[0][0]->GetXaxis()->GetBinUpEdge(Smeared3[0][0]->GetNbinsX()),   Smeared3[0][0]->GetNbinsY(), Smeared3[0][0]->GetYaxis()->GetBinLowEdge(1), Smeared3[0][0]->GetYaxis()->GetBinUpEdge(Smeared3[0][0]->GetNbinsY()));
	MRC_3[ch][term-1]->Add(Ideal3[ch][realterm-1], MRCvariation);
	MRC_3[ch][term-1]->Add(Smeared3[ch][realterm-1], 1-MRCvariation);
	MRC_3[ch][term-1]->Divide(Smeared3[ch][realterm-1]);
	

	for(int binX=1; binX<=MRC_3[ch][term-1]->GetNbinsX(); binX++){
	  for(int binY=1; binY<=MRC_3[ch][term-1]->GetNbinsY(); binY++){
	    MRC_3[ch][term-1]->SetBinError(binX,binY, 0);
	    if(PbPbcase){
	      if(MRC_3[ch][term-1]->GetYaxis()->GetBinCenter(binY) > 0.13) {
		MRC_3[ch][term-1]->SetBinContent(binX, binY,  MRC_3[ch][term-1]->GetBinContent(binX, MRC_3[ch][term-1]->GetYaxis()->FindBin(0.13)));
	      }
	    }
	  }
	}

      }
      
    }
    //
    for(int term=1; term<=13; term++){
      TString *S4 = new TString("FourParticle_Charge1_");
      if(ch==0) S4->Append("1_Charge2_1_Charge3_1_Charge4_1_M_0_ED_0_Term_");
      else if(ch==1) S4->Append("0_Charge2_1_Charge3_1_Charge4_1_M_0_ED_0_Term_");
      else S4->Append("0_Charge2_0_Charge3_1_Charge4_1_M_0_ED_0_Term_");
      *S4 += term;
      TString *I4 = new TString(S4->Data());
      S4->Append("_Smeared"); I4->Append("_Ideal");
      Smeared4[ch][term-1]=(TH2D*)MyList->FindObject(S4->Data());
      Ideal4[ch][term-1]=(TH2D*)MyList->FindObject(I4->Data());
      //
      Smeared4[ch][term-1]->RebinY(FourParticleRebin);
      Ideal4[ch][term-1]->RebinY(FourParticleRebin);
    }
    for(int term=1; term<=12; term++){
      if(ch==COI_4 && term==TOI_4){
	C4Smeared = (TH1D*)Smeared4[ch][term-1]->ProjectionY("C4Smeared",Rbin,Rbin);
	TH1D *temp_s_4 = (TH1D*)Smeared4[ch][12]->ProjectionY("temp_s_4",Rbin,Rbin);
	C4Smeared->Divide(temp_s_4);
	C4Ideal = (TH1D*)Ideal4[ch][term-1]->ProjectionY("C4Ideal",Rbin,Rbin);
	TH1D *temp_i_4 = (TH1D*)Ideal4[ch][12]->ProjectionY("temp_i_4",Rbin,Rbin);
	C4Ideal->Divide(temp_i_4);
      }
    }
    for(int term=1; term<=6; term++){// 6 unique terms at max
      TString *name = new TString("MRC_4_");
      if(ch==0) name->Append("SC_term");
      else if(ch==1) name->Append("MC1_term");
      else name->Append("MC2_term");
      *name += term;
      int realterm = term;// realterm starts from 1
      if(ch==0 && term>5) continue;// ++++ case
            
      if(ch==0 && term==2) realterm=5;
      if(ch==0 && term==3) realterm=8;
      if(ch==0 && term==4) realterm=12;
      if(ch==0 && term==5) realterm=13;

      if(ch==1 && term==2) realterm=4;//4 
      if(ch==1 && term==3) realterm=5;//5 
      if(ch==1 && term==4) realterm=8;//8
      if(ch==1 && term==5) realterm=10;//10
      if(ch==1 && term==6) realterm=13;//13
      
      if(ch==2 && term==2) realterm=3;//3
      if(ch==2 && term==3) realterm=11;//11
      if(ch==2 && term==4) realterm=8;//8
      if(ch==2 && term==5) realterm=12;//12
      if(ch==2 && term==6) realterm=13;//13
      
      MRC_4[ch][term-1] = new TH2D(name->Data(),"", Smeared4[0][0]->GetNbinsX(), Smeared4[0][0]->GetXaxis()->GetBinLowEdge(1), Smeared4[0][0]->GetXaxis()->GetBinUpEdge(Smeared4[0][0]->GetNbinsX()),   Smeared4[0][0]->GetNbinsY(), Smeared4[0][0]->GetYaxis()->GetBinLowEdge(1), Smeared4[0][0]->GetYaxis()->GetBinUpEdge(Smeared4[0][0]->GetNbinsY()));
      MRC_4[ch][term-1]->Add(Ideal4[ch][realterm-1], MRCvariation);
      MRC_4[ch][term-1]->Add(Smeared4[ch][realterm-1], 1-MRCvariation);
      MRC_4[ch][term-1]->Divide(Smeared4[ch][realterm-1]);
      
      for(int binX=1; binX<=MRC_4[ch][term-1]->GetNbinsX(); binX++){
	for(int binY=1; binY<=MRC_4[ch][term-1]->GetNbinsY(); binY++){
	  MRC_4[ch][term-1]->SetBinError(binX,binY, 0);
	  if(PbPbcase){
	    if(MRC_4[ch][term-1]->GetYaxis()->GetBinCenter(binY) > 0.19) MRC_4[ch][term-1]->SetBinContent(binX, binY,  MRC_4[ch][term-1]->GetBinContent(binX, MRC_4[ch][term-1]->GetYaxis()->FindBin(0.19)));
	  }
	}
      }

    }
    //
  }

  /////////////////////////////////////////


  C2Ideal->SetLineColor(4);
  C2Smeared->SetLineColor(2);
  C2Ideal->GetXaxis()->SetRangeUser(0,0.15);
  C2Ideal->SetMinimum(0.98); C2Ideal->SetMaximum(1.35);
  C2Ideal->GetYaxis()->SetTitleOffset(1.5);
  C2Ideal->GetXaxis()->SetTitle("q_{inv} (GeV/c)");
  C2Ideal->GetYaxis()->SetTitle("C_{2}^{Ideal}");
  //
  C3Ideal->SetLineColor(4);
  C3Smeared->SetLineColor(2);
  C3Ideal->GetXaxis()->SetRangeUser(0,0.15);
  C3Ideal->SetMinimum(0.98); C3Ideal->SetMaximum(1.6);
  C3Ideal->GetYaxis()->SetTitleOffset(1.5);
  C3Ideal->GetXaxis()->SetTitle("Q_{3} (GeV/c)");
  C3Ideal->GetYaxis()->SetTitle("C_{3}^{Ideal}");
  //
  C4Ideal->SetLineColor(4);
  C4Smeared->SetLineColor(2);
  C4Ideal->GetXaxis()->SetRangeUser(0,0.15);
  C4Ideal->SetMinimum(0.98); C4Ideal->SetMaximum(4.0);
  C4Ideal->GetYaxis()->SetTitleOffset(1.5);
  C4Ideal->GetXaxis()->SetTitle("Q_{4} (GeV/c)");
  C4Ideal->GetYaxis()->SetTitle("C_{4}^{Ideal}");
  //
  //C2Ideal->Draw();
  //C2Smeared->Draw("same");
  //
  C3Ideal->Draw();
  C3Smeared->Draw("same");
  //
  //C4Ideal->Draw();
  //C4Smeared->Draw("same");
  //
  //
  //TH1D *MRC_3=(TH1D*)C3Ideal->Clone();
  //MRC_3->Divide(C3Smeared);
  //MRC_3->SetMinimum(0.8); MRC_3->SetMaximum(1.25); 
  //MRC_3->Draw();
  //
  //TH1D *MRC_4=(TH1D*)C4Ideal->Clone();
  //MRC_4->Divide(C4Smeared);
  //MRC_4->SetMinimum(0.8); MRC_4->SetMaximum(1.25); 
  //MRC_4->Draw();
  
  //
  //

  for(int ch=0; ch<=1; ch++){
    for(int term=1; term<=2; term++){
      MRC_2[ch][term-1]->Write();
    }
    if(ch==0){
      TH2D *MRC_C2_SC = (TH2D*)MRC_2[ch][0]->Clone();
      MRC_C2_SC->Divide(MRC_2[ch][1]);
      MRC_C2_SC->SetName("MRC_C2_SC");
      MRC_C2_SC->SetTitle("MRC_C2_SC");
      MRC_C2_SC->Write();
    }else{
      TH2D *MRC_C2_MC = (TH2D*)MRC_2[ch][0]->Clone();
      MRC_C2_MC->Divide(MRC_2[ch][1]);
      MRC_C2_MC->SetName("MRC_C2_MC");
      MRC_C2_MC->SetTitle("MRC_C2_MC");
      MRC_C2_MC->Write();
    }
  }
  cout<<"here"<<endl;
  for(int ch=0; ch<=1; ch++){
    for(int term=1; term<=4; term++){
      if(ch==0 && term==1) cout<<MRC_3[ch][term-1]->GetBinContent(10,2)/MRC_3[ch][2]->GetBinContent(10,2)<<endl;
      if(MRC_3[ch][term-1]) MRC_3[ch][term-1]->Write();
    }
  }
  for(int ch=0; ch<=2; ch++){
    for(int term=1; term<=6; term++){
      if(MRC_4[ch][term-1]) MRC_4[ch][term-1]->Write();
    }
  }
  
  //fout->Close();

  
}
