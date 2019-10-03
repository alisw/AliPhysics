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
#include "TPaveStats.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TPad.h"

#define PI 3.1415926
#define BohrR 1963.6885 // Mate's value(1963.6885) ~ 387.5/0.197327(1963.7455)
#define FmToGeV 0.197327 // 0.197327

using namespace std;

void PlotTherm_r3(){

  double epsilon=1.0;// chaoticity

  gStyle->SetOptStat(0);
  //gStyle->SetOptFit(111);
 
  //TFile *infile=new TFile("Therm_r3_b9.root","READ");// For v2 paper (r<100 and r*<80)
  TFile *infile=new TFile("Therm_FSI_b2.root","READ");
  //TFile *infile=new TFile("Therm_dist_temp.root","READ");
 
  //
  TH3D *r3num_init=(TH3D*)infile->Get("r3num");
  TH3D *r3den1_init=(TH3D*)infile->Get("r3den1");
  TH3D *r3den2_init=(TH3D*)infile->Get("r3den2");
  TH3D *r3den3_init=(TH3D*)infile->Get("r3den3");
  TH3D *r3numSq_init=(TH3D*)infile->Get("r3numSq");
  TH3D *r3den1Sq_init=(TH3D*)infile->Get("r3den1Sq");
  TH3D *r3den2Sq_init=(TH3D*)infile->Get("r3den2Sq");
  TH3D *r3den3Sq_init=(TH3D*)infile->Get("r3den3Sq");
  TH3D *r3numEn_init=(TH3D*)infile->Get("r3numEn");
  TH3D *r3den1En_init=(TH3D*)infile->Get("r3den1En");
  TH3D *r3den2En_init=(TH3D*)infile->Get("r3den2En");
  TH3D *r3den3En_init=(TH3D*)infile->Get("r3den3En");
  TH3D *r3numME_init=(TH3D*)infile->Get("r3numME");
  TH3D *r3den1ME_init=(TH3D*)infile->Get("r3den1ME");
  TH3D *r3den2ME_init=(TH3D*)infile->Get("r3den2ME");
  TH3D *r3den3ME_init=(TH3D*)infile->Get("r3den3ME");
  
  TH2D *R2_2D = (TH2D*)infile->Get("Num_Cos_ss");
  //TH2D *R2Sq_2D = (TH2D*)infile->Get("NumSq_Cos_ss");
  TH2D *R2En_2D = (TH2D*)infile->Get("Den_ss");
  TH1D *R2=(TH1D*)R2_2D->ProjectionY();
  //TH1D *R2Sq=(TH1D*)R2Sq_2D->ProjectionY();
  TH1D *R2En=(TH1D*)R2En_2D->ProjectionY();
  R2->Divide(R2En);
  //R2Sq->Divide(R2En);

  TH3D *r3num=(TH3D*)r3num_init->Clone();
  TH3D *r3den1=(TH3D*)r3den1_init->Clone();
  TH3D *r3den2=(TH3D*)r3den2_init->Clone();
  TH3D *r3den3=(TH3D*)r3den3_init->Clone();
  TH3D *r3numSq=(TH3D*)r3numSq_init->Clone();
  TH3D *r3den1Sq=(TH3D*)r3den1Sq_init->Clone();
  TH3D *r3den2Sq=(TH3D*)r3den2Sq_init->Clone();
  TH3D *r3den3Sq=(TH3D*)r3den3Sq_init->Clone();
  TH3D *r3numEn=(TH3D*)r3numEn_init->Clone();
  TH3D *r3den1En=(TH3D*)r3den1En_init->Clone();
  TH3D *r3den2En=(TH3D*)r3den2En_init->Clone();
  TH3D *r3den3En=(TH3D*)r3den3En_init->Clone();
  TH3D *r3numME=(TH3D*)r3numME_init->Clone();
  TH3D *r3den1ME=(TH3D*)r3den1ME_init->Clone();
  TH3D *r3den2ME=(TH3D*)r3den2ME_init->Clone();
  TH3D *r3den3ME=(TH3D*)r3den3ME_init->Clone();


  int FB=40;
  int LB=50;
  /*double SF = r3numEn->Integral(FB,LB,FB,LB,FB,LB)/r3numME->Integral(FB,LB,FB,LB,FB,LB);
  cout<<SF<<endl;
  r3numME->Scale(SF);
  SF = r3den1En->Integral(FB,LB,FB,LB,FB,LB)/r3den1ME->Integral(FB,LB,FB,LB,FB,LB);
  cout<<SF<<endl;
  r3den1ME->Scale(SF);
  SF = r3den2En->Integral(FB,LB,FB,LB,FB,LB)/r3den2ME->Integral(FB,LB,FB,LB,FB,LB);
  cout<<SF<<endl;
  r3den2ME->Scale(SF);
  SF = r3den3En->Integral(FB,LB,FB,LB,FB,LB)/r3den3ME->Integral(FB,LB,FB,LB,FB,LB);
  cout<<SF<<endl;
  r3den3ME->Scale(SF);*/
  //
  /*
  //////////////////////
  // Symmetrize bin contents
  for(int i=1; i<=50; i++){
    for(int j=1; j<=50; j++){
      for(int k=1; k<=50; k++){
	// 3-pion phase
	double num_sum = r3num_init->GetBinContent(i,j,k) + r3num_init->GetBinContent(i,k,j) + r3num_init->GetBinContent(j,i,k);
	num_sum += r3num_init->GetBinContent(j,k,i) + r3num_init->GetBinContent(k,i,j) + r3num_init->GetBinContent(k,j,i);
	double Sq_sum = r3numSq_init->GetBinContent(i,j,k) + r3numSq_init->GetBinContent(i,k,j) + r3numSq_init->GetBinContent(j,i,k);
	Sq_sum += r3numSq_init->GetBinContent(j,k,i) + r3numSq_init->GetBinContent(k,i,j) + r3numSq_init->GetBinContent(k,j,i);
	double den_sum = r3numEn_init->GetBinContent(i,j,k) + r3numEn_init->GetBinContent(i,k,j) + r3numEn_init->GetBinContent(j,i,k);
	den_sum += r3numEn_init->GetBinContent(j,k,i) + r3numEn_init->GetBinContent(k,i,j) + r3numEn_init->GetBinContent(k,j,i);
	if(den_sum>0) {
	  r3num->SetBinContent(i,j,k, num_sum/den_sum); r3num->SetBinContent(i,k,j, num_sum/den_sum); 
	  r3num->SetBinContent(j,i,k, num_sum/den_sum); r3num->SetBinContent(j,k,i, num_sum/den_sum);
	  r3num->SetBinContent(k,i,j, num_sum/den_sum); r3num->SetBinContent(k,j,i, num_sum/den_sum);
	  r3numSq->SetBinContent(i,j,k, Sq_sum/den_sum); r3numSq->SetBinContent(i,k,j, Sq_sum/den_sum); 
	  r3numSq->SetBinContent(j,i,k, Sq_sum/den_sum); r3numSq->SetBinContent(j,k,i, Sq_sum/den_sum);
	  r3numSq->SetBinContent(k,i,j, Sq_sum/den_sum); r3numSq->SetBinContent(k,j,i, Sq_sum/den_sum);
	  r3numEn->SetBinContent(i,j,k, den_sum); r3numEn->SetBinContent(i,k,j, den_sum); 
	  r3numEn->SetBinContent(j,i,k, den_sum); r3numEn->SetBinContent(j,k,i, den_sum);
	  r3numEn->SetBinContent(k,i,j, den_sum); r3numEn->SetBinContent(k,j,i, den_sum);
	}else {
	  r3num->SetBinContent(i,j,k, 0.); r3num->SetBinContent(i,k,j, 0.);
	  r3num->SetBinContent(j,i,k, 0.); r3num->SetBinContent(j,k,i, 0.);
	  r3num->SetBinContent(k,i,j, 0.); r3num->SetBinContent(k,j,i, 0.);
	  r3numEn->SetBinContent(i,j,k, 0.); r3numEn->SetBinContent(i,k,j, 0.); 
	  r3numEn->SetBinContent(j,i,k, 0.); r3numEn->SetBinContent(j,k,i, 0.);
	  r3numEn->SetBinContent(k,i,j, 0.); r3numEn->SetBinContent(k,j,i, 0.);
	}
	// 2-pion phases
	// 12
	num_sum = r3den1_init->GetBinContent(i,j,k) + r3den1_init->GetBinContent(i,k,j);
	Sq_sum = r3den1Sq_init->GetBinContent(i,j,k) + r3den1Sq_init->GetBinContent(i,k,j);
	den_sum = r3numEn_init->GetBinContent(i,j,k) + r3numEn_init->GetBinContent(i,k,j);
	if(den_sum>0) {
	  r3den1->SetBinContent(i,j,k, num_sum/den_sum); r3den1->SetBinContent(i,k,j, num_sum/den_sum); 
	  r3den1Sq->SetBinContent(i,j,k, Sq_sum/den_sum); r3den1Sq->SetBinContent(i,k,j, Sq_sum/den_sum);
	  r3den1En->SetBinContent(i,j,k, den_sum); r3den1En->SetBinContent(i,k,j, den_sum); 
	}else{
	  r3den1->SetBinContent(i,j,k, 0.); r3den1->SetBinContent(i,k,j, 0.); 
	  r3den1Sq->SetBinContent(i,j,k, 0.); r3den1Sq->SetBinContent(i,k,j, 0.);
	  r3den1En->SetBinContent(i,j,k, 0.); r3den1En->SetBinContent(i,k,j, 0.);
	}
	// 23
	num_sum = r3den2_init->GetBinContent(j,i,k) + r3den2_init->GetBinContent(k,i,j);
	Sq_sum = r3den2Sq_init->GetBinContent(j,i,k) + r3den2Sq_init->GetBinContent(k,i,j);
	den_sum = r3numEn_init->GetBinContent(j,i,k) + r3numEn_init->GetBinContent(k,i,j);
	if(den_sum>0) {
	  r3den2->SetBinContent(j,i,k, num_sum/den_sum); r3den2->SetBinContent(k,i,j, num_sum/den_sum); 
	  r3den2Sq->SetBinContent(j,i,k, Sq_sum/den_sum); r3den2Sq->SetBinContent(k,i,j, Sq_sum/den_sum);
	  r3den2En->SetBinContent(j,i,k, den_sum); r3den2En->SetBinContent(k,i,j, den_sum);
	}else{
	  r3den2->SetBinContent(j,i,k, 0.); r3den2->SetBinContent(k,i,j, 0.); 
	  r3den2Sq->SetBinContent(j,i,k, 0.); r3den2Sq->SetBinContent(k,i,j, 0.);
	  r3den2En->SetBinContent(j,i,k, 0.); r3den2En->SetBinContent(k,i,j, 0.);
	}
	// 12
	num_sum = r3den3_init->GetBinContent(j,k,i) + r3den3_init->GetBinContent(k,j,i);
	Sq_sum = r3den3Sq_init->GetBinContent(j,k,i) + r3den3Sq_init->GetBinContent(k,j,i);
	den_sum = r3numEn_init->GetBinContent(j,k,i) + r3numEn_init->GetBinContent(k,j,i);
	if(den_sum>0) {
	  r3den3->SetBinContent(j,k,i, num_sum/den_sum); r3den3->SetBinContent(k,j,i, num_sum/den_sum); 
	  r3den3Sq->SetBinContent(j,k,i, Sq_sum/den_sum); r3den3Sq->SetBinContent(k,j,i, Sq_sum/den_sum);
	  r3den3En->SetBinContent(j,k,i, den_sum); r3den3En->SetBinContent(k,j,i, den_sum);
	}else{
	  r3den3->SetBinContent(j,k,i, 0.); r3den3->SetBinContent(k,j,i, 0.); 
	  r3den3Sq->SetBinContent(j,k,i, 0.); r3den3Sq->SetBinContent(k,j,i, 0.);
	  r3den3En->SetBinContent(j,k,i, 0.); r3den3En->SetBinContent(k,j,i, 0.);
	}
	
      }
    }
  }
  cout<<"Done Symmetrizing"<<endl;
  */

  r3num->Sumw2();// r3numEn->Sumw2();
  r3den1->Sumw2();// r3den1En->Sumw2();
  r3den2->Sumw2();// r3den2En->Sumw2();
  r3den3->Sumw2();// r3den3En->Sumw2();




  r3num->Divide(r3numEn);
  r3den1->Divide(r3den1En);
  r3den2->Divide(r3den2En);
  r3den3->Divide(r3den3En);
  /*r3num->Divide(r3num,r3numEn,1,1,"B");
  r3den1->Divide(r3den1,r3den1En,1,1,"B");
  r3den2->Divide(r3den2,r3den2En,1,1,"B");
  r3den3->Divide(r3den3,r3den3En,1,1,"B");*/
  r3numSq->Divide(r3numEn);
  r3den1Sq->Divide(r3den1En);
  r3den2Sq->Divide(r3den2En);
  r3den3Sq->Divide(r3den3En);
  



  TH1D *r3=new TH1D("r3","",10,0,0.1);
  TH1D *r3den=new TH1D("r3den","",10,0,0.1);
  double r3num_e[20]={0};
  double r3den_e[20]={0};
  double NegativeDenCount[20]={0};
  double TotalDenCount[20]={0};
  double r3den_neg[20]={0};
  double r3den_neg_Sq[20]={0};

  for(int i=1; i<=50; i++){
    for(int j=1; j<=50; j++){
      for(int k=1; k<=50; k++){

	double q3 =sqrt(pow(0.002*(i-0.5),2)+pow(0.002*(j-0.5),2)+pow(0.002*(k-0.5),2));
	if(q3>=0.2) continue;

	if(r3numEn->GetBinContent(i,j,k)==0) continue;
	if(r3den1En->GetBinContent(i,j,k)==0) continue;
	if(r3den2En->GetBinContent(i,j,k)==0) continue;
	if(r3den3En->GetBinContent(i,j,k)==0) continue;

	
	double num = r3num->GetBinContent(i,j,k);
	num *= pow(epsilon,0.5) * (3-2*epsilon)/pow(2-epsilon,1.5);
	double den1 = r3den1->GetBinContent(i,j,k);
	double den2 = r3den2->GetBinContent(i,j,k);
	double den3 = r3den3->GetBinContent(i,j,k);
	double den = den1*den2*den3;
	
	TotalDenCount[int(q3/0.01)]++;
	if(den <= 0) {
	  NegativeDenCount[int(q3/0.01)]++;
	  r3den_neg[int(q3/0.01)] += sqrt(fabs(den));
	  r3den_neg_Sq[int(q3/0.01)] += fabs(den);
	  continue;
	}
	den = sqrt(den);
	double weightWidth = sqrt( (r3den1Sq->GetBinContent(i,j,k) - pow(r3den1->GetBinContent(i,j,k),2))/r3den1En->GetBinContent(i,j,k));
	double den_e = pow(0.5*weightWidth*den2*den3/den,2);
	weightWidth = sqrt( fabs(r3den2Sq->GetBinContent(i,j,k) - pow(r3den2->GetBinContent(i,j,k),2))/r3den2En->GetBinContent(i,j,k));
	den_e += pow(0.5*weightWidth*den1*den3/den,2);
	weightWidth = sqrt( fabs(r3den3Sq->GetBinContent(i,j,k) - pow(r3den3->GetBinContent(i,j,k),2))/r3den3En->GetBinContent(i,j,k));
	den_e += pow(0.5*weightWidth*den1*den2/den,2);
	//
	r3->Fill(q3,num);
	r3den->Fill(q3,den);
	//
	double num_e = (r3numSq->GetBinContent(i,j,k)-pow(r3num->GetBinContent(i,j,k),2))/r3numEn->GetBinContent(i,j,k);
	r3num_e[int(q3/0.01)] += num_e;
	r3den_e[int(q3/0.01)] += den_e;
	

	//if(q3<0.01) cout<<i<<j<<k<<"   "<<r3num_e[0]<<"  "<<r3den_e[0]<<"  "<<r3numSq->GetBinContent(i,j,k)<<"  "<<pow(r3num->GetBinContent(i,j,k),2)<<"  "<<r3numEn->GetBinContent(i,j,k)<<endl;
	//if(q3<0.02) cout<<r3den1->GetBinContent(i,j,k)<<"  "<<r3den1->GetBinContent(i,k,j)<<endl;
      }
    }
  }

  cout<<"Done calculating r3"<<endl;

  TCanvas *can1 = new TCanvas("can1", "can1",10,0,600,600);// 11,53,700,500
  can1->SetHighLightColor(2);
  //can1->Range(-0.7838115,-1.033258,0.7838115,1.033258);
  can1->SetFillColor(10);//10
  can1->SetBorderMode(0);
  can1->SetBorderSize(2);
  can1->SetFrameFillColor(0);
  can1->SetFrameBorderMode(0);
  can1->SetFrameBorderMode(0);
  
  
  TPad *pad1 = new TPad("pad1","pad1",0,0,1.,1.);
  gPad->SetGridx(0);
  gPad->SetGridy(0);
  gPad->SetTickx();
  gPad->SetTicky();
  pad1->SetTopMargin(0.02);
  pad1->SetRightMargin(0.03);
  pad1->SetBottomMargin(0.06);//0.12
  pad1->Draw();
  
  TLegend *legend1 = new TLegend(.15,.25,.65,.45,NULL,"brNDC");
  legend1->SetBorderSize(1);
  legend1->SetTextSize(.04);// small .03; large .036 
  //legend1->SetLineColor(0);
  legend1->SetFillColor(0);

  TF1 *ChaoticLimit=new TF1("ChaoticLimit","2",0,0.2);
  ChaoticLimit->SetLineStyle(3);
  TF1 *Quartic=new TF1("Quartic","[0]*(1-[1]*pow(x,4))",0,0.1);
  Quartic->SetLineColor(4);
  TF1 *Quadratic=new TF1("Quadratic","[0]*(1-[1]*pow(x,2))",0,0.1);
  //Quadratic->SetLineStyle(2);
  Quadratic->SetLineColor(6);

  
  
  for(int i=0; i<10; i++){
    cout<<"Lost Triplet Fraction: "<<NegativeDenCount[i]/TotalDenCount[i]<<endl;
    r3->SetBinError(i+1,sqrt(r3num_e[i]));
    double ExtraDenError = 0;
    if(NegativeDenCount[i]>0) ExtraDenError = sqrt( fabs(r3den_neg_Sq[i]/NegativeDenCount[i] - pow(r3den_neg[i]/NegativeDenCount[i],2))/NegativeDenCount[i]) * NegativeDenCount[i];
    r3den->SetBinContent(i+1, r3den->GetBinContent(i+1) - r3den_neg[i]);
    r3den->SetBinError(i+1,sqrt(r3den_e[i])+r3den_neg[i]);
  }
  //r3->Sumw2(); r3den->Sumw2();
  //r3->Divide(r3, r3den, 1,1,"B");
  r3->Divide(r3den);
  

  // first bin sometimes only has 1 entry (bad error bar)
  if(r3->GetBinError(1)<0.001) r3->SetBinError(1, r3->GetBinError(2));
  
  gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.03); gPad->SetBottomMargin(0.1); gPad->SetTopMargin(0.01);
  r3->SetMarkerStyle(24);
  r3->SetMarkerColor(1);
  r3->SetLineColor(1);
  r3->GetXaxis()->SetLabelFont(63); r3->GetYaxis()->SetLabelFont(63);
  r3->GetXaxis()->SetLabelSize(18); r3->GetYaxis()->SetLabelSize(18);
  r3->GetXaxis()->SetTitleFont(63); r3->GetYaxis()->SetTitleFont(63);
  r3->GetXaxis()->SetTitleSize(32); r3->GetYaxis()->SetTitleSize(32);
  r3->GetXaxis()->SetTitleOffset(0.75); r3->GetYaxis()->SetTitleOffset(1.);
  r3->GetXaxis()->SetTitle("Q_{3} (GeV/c)");
  r3->GetYaxis()->SetTitle("r_{3}");
  r3->SetTitle("");
  r3->SetMinimum(1.1); r3->SetMaximum(2.5);
  r3->GetXaxis()->SetRangeUser(0,0.069);
  r3->Draw();
 
  //r3->Fit(Quartic,"IME","",0,0.09);
  //r3->Fit(Quadratic,"IME","",0,0.09);
  ChaoticLimit->Draw("same");
  //Quartic->Draw("same");
  //Quadratic->Draw("same");
  //legend1->AddEntry(Quartic,"Quartic: #chi^{2}/NDF=1.2","l");
  //legend1->AddEntry(Quadratic,"Quadratic: #chi^{2}/NDF=5.4","l");
  //legend1->Draw("same");
  
  //cout<<"Quartic Chi2/NDF = "<<Quartic->GetChisquare()/Quartic->GetNDF()<<endl;
  //cout<<"Quadratic Chi2/NDF = "<<Quadratic->GetChisquare()/Quadratic->GetNDF()<<endl;
  
  TLatex *Specif = new TLatex(0.005,1.4,"Therminator");// Therminator
  Specif->SetTextSize(.075);
  Specif->Draw("same");
  
}
