//gROOT->Reset();
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <Riostream.h>
#include "RooFit.h"
#include "RooVoigtian.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooComplex.h"
#include "RooMath.h"

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
#include "TLatex.h"
#include "TImage.h"
#include "TGaxis.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"

using namespace std;


double myLevyPtFunc(Double_t *x, Double_t *par);
double myMtExpFunc(Double_t *x, Double_t *par);
double myPtExpFunc(Double_t *x, Double_t *par);
double myBoltzFunc(Double_t *x, Double_t *par);
double myPowerLawPtFunc(Double_t *x, Double_t *par);
double myPowerLawFunc(Double_t *x, Double_t *par);
double myBoltsBWFunc(Double_t *x, Double_t *par);
double IntegrandBG(const double *x, const double *p);

void Plot_XiStarResults(){
  
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(01111);
    

  TFile *filesCutVar[13];
  TH1D *Spectra_CutVar[13];
  TF1 *LevyFits_CutVar[13];
  TH1D *Yields_CutVar[13];
  TH1D *Eff_CutVar[13];
  for(int i=0; i<13; i++){
    TString *name=new TString("CutVariation_");
    *name += i;
    name->Append(".root");
    filesCutVar[i]=new TFile(name->Data());
    Spectra_CutVar[i] = (TH1D*)filesCutVar[i]->Get("h_Spectrum");
    Spectra_CutVar[i]->SetDirectory(0);
    LevyFits_CutVar[i] = (TF1*)filesCutVar[i]->Get("myLevy");
    Yields_CutVar[i] = (TH1D*)filesCutVar[i]->Get("RawYields");
    Eff_CutVar[i] = (TH1D*)filesCutVar[i]->Get("Efficiency");
  }
  //
  TFile *filesFitVar[4];
  TH1D *Spectra_FitVar[4];
  TF1 *LevyFits_FitVar[4];
  TH1D *Yields_FitVar[4];
  for(int i=0; i<4; i++){
    TString *name=new TString("FitVariation_");
    *name += i+1;
    name->Append(".root");
    filesFitVar[i]=new TFile(name->Data());
    Spectra_FitVar[i] = (TH1D*)filesFitVar[i]->Get("h_Spectrum");
    Spectra_FitVar[i]->SetDirectory(0);
    LevyFits_FitVar[i] = (TF1*)filesFitVar[i]->Get("myLevy");
    Yields_FitVar[i] = (TH1D*)filesFitVar[i]->Get("RawYields");
  }
  
  double FVavg=0;
  for(int bin=2; bin<=9; bin++){
    // FV 1 has largest spectrum point variation
    FVavg += fabs(Spectra_FitVar[0]->GetBinContent(bin)-Spectra_FitVar[1]->GetBinContent(bin))/Spectra_FitVar[0]->GetBinContent(bin);
  }
  cout<<"Average point-by-point spectrum variation due to fitting: "<<FVavg/8<<endl;

  double TrackVavg=0;
  for(int bin=2; bin<=9; bin++){
    // CutVar 5 has largest spectrum point variation in the "tracking systematic" category
    TrackVavg += fabs(Spectra_CutVar[0]->GetBinContent(bin)-Spectra_CutVar[5]->GetBinContent(bin))/Spectra_CutVar[0]->GetBinContent(bin);
  }
  cout<<"Average point-by-point spectrum variation due to Tracking: "<<TrackVavg/8<<endl;

  double TopologVavg=0;
  for(int bin=2; bin<=9; bin++){
    // CutVar 8 has largest spectrum point variation in the "topological systematic" category
    TopologVavg += fabs(Spectra_CutVar[0]->GetBinContent(bin)-Spectra_CutVar[8]->GetBinContent(bin))/Spectra_CutVar[0]->GetBinContent(bin);
  }
  cout<<"Average point-by-point spectrum variation due to Topological Cuts: "<<TopologVavg/8<<endl;

  ////////////////////////////////////
  TH1D *FinalSpectrum = (TH1D*)Spectra_CutVar[0]->Clone();
  FinalSpectrum->SetMarkerSize(1.);
  double StatError[8]={0};
  for(int ptbin=1; ptbin<=8; ptbin++) StatError[ptbin-1] = FinalSpectrum->GetBinError(ptbin+1);// add 1 since first bin is a dummy
  

  TH1D *PureSysFromCutVar[12];
  double MeanSysCutVar[12]={0};
  double sigmaSysCutVar[12]={0};
  for(int i=0; i<12; i++) {
    TString *name=new TString("PureSysFromCutVar");
    *name += i+1;
    PureSysFromCutVar[i] = new TH1D(name->Data(),"",8,0.5,8.5);
    PureSysFromCutVar[i]->GetXaxis()->SetTitle("pt bin");
    PureSysFromCutVar[i]->GetYaxis()->SetTitle("|y_{s}-y_{v}|/y_{s}");
  }
  
  //////////////////////////////////////
  // Sys error assignment for pt bins to be used for mean value fits

  int CutVarLow=1, CutVarHigh=13;
  for(int i=CutVarLow; i<CutVarHigh; i++){// 1 to 13
    for(int ptbin=2; ptbin<=9; ptbin++){
      double SysVar = fabs(FinalSpectrum->GetBinContent(ptbin) - Spectra_CutVar[i]->GetBinContent(ptbin))/FinalSpectrum->GetBinContent(ptbin);
      double sigma = sqrt(fabs(pow(FinalSpectrum->GetBinError(ptbin),2) - pow(Spectra_CutVar[i]->GetBinError(ptbin),2)));// stat part
      PureSysFromCutVar[i-1]->SetBinContent(ptbin-1, SysVar);
      PureSysFromCutVar[i-1]->SetBinError(ptbin-1, sigma/FinalSpectrum->GetBinContent(ptbin));
    }
    PureSysFromCutVar[i-1]->Fit("pol0","IMEQ","",0,9);
    MeanSysCutVar[i-1] = pol0->GetParameter(0);
    sigmaSysCutVar[i-1] = pol0->GetParError(0);
    cout<<"Cut Var "<<i<<"  Mean Sys = "<<MeanSysCutVar[i-1]<<"   += "<<pol0->GetParError(0)<<endl;
    //cout<<"Chi2/NDF = "<<pol0->GetChisquare()/pol0->GetNDF()<<endl;
    //TPaveStats *sb2 = (TPaveStats*)(PureSysFromCutVar[i-1]->GetListOfFunctions()->FindObject("stats"));
    //sb2->SetX1NDC(.2);
  }


  double SysError[9]={0};
  // Cut Variations
  //int CutVarN=12;
  for(int i=CutVarLow; i<CutVarHigh; i++){// 1 to 13
    for(int ptbin=1; ptbin<=9; ptbin++){
      if(MeanSysCutVar[i-1] > 0 && MeanSysCutVar[i-1]>sigmaSysCutVar[i-1]) {
	SysError[ptbin-1] += pow(MeanSysCutVar[i-1]*FinalSpectrum->GetBinContent(ptbin),2) - pow(sigmaSysCutVar[i-1]*FinalSpectrum->GetBinContent(ptbin),2);
      }
    }
  }
  // Fit Variations
  for(int i=2; i<3; i++){// now only use i=2 (max deviation).  was 0 to 4
    for(int ptbin=1; ptbin<=9; ptbin++){
      SysError[ptbin-1] += pow(FinalSpectrum->GetBinContent(ptbin) - Spectra_FitVar[i]->GetBinContent(ptbin),2);
    }
  }
  
  // Add Stat errors and set Final Error
  for(int ptbin=1; ptbin<=9; ptbin++){
    double TotError = pow(FinalSpectrum->GetBinError(ptbin),2);// Stat
    TotError += SysError[ptbin-1];// Add on Sys 
    FinalSpectrum->SetBinError(ptbin, sqrt(TotError));
    cout<<"Spectrum ptbin:"<<ptbin-1<<"  "<<FinalSpectrum->GetBinContent(ptbin)<<"  +- "<<FinalSpectrum->GetBinError(ptbin)<<endl;
  }
 
  for(int ptbin=1; ptbin<=9; ptbin++) SysError[ptbin-1] = sqrt(SysError[ptbin-1]);
  
  // print efficiencies
  for(int ptbin=1; ptbin<=8; ptbin++){
    cout<<"Efficiency ptbin:"<<ptbin<<"  "<<Eff_CutVar[0]->GetBinContent(ptbin)<<"  +- "<<Eff_CutVar[0]->GetBinError(ptbin)<<endl;
  }
 
  TCanvas *can = new TCanvas("can","can",13,34,700,500);
  can->Range(-1.25,-0.2625,11.25,2.3625);
  can->SetFillColor(10);
  can->SetBorderMode(0);
  can->SetBorderSize(2);
  //can->SetGridx();
  //can->SetGridy();
  can->SetFrameFillColor(0);
  can->SetFrameBorderMode(0);
  can->SetFrameBorderMode(0);
  
  TLegend *legend = new TLegend(.1,.7,.35,.9,NULL,"brNDC");
  legend->SetBorderSize(1);
  legend->SetTextSize(.04);// small .03; large .036 
  //legend->SetLineColor(0);
  legend->SetFillColor(0);
  TLegend *legend2 = new TLegend(.35,.62,.9,.72,NULL,"brNDC");
  legend2->SetBorderSize(1);
  legend2->SetTextSize(.04);// small .03; large .036 
  //legend2->SetLineColor(0);
  legend2->SetFillColor(0);
  

  double minFitpoint = 0.8; double maxFitpoint = 5.6;
  

  
  TF1 *myLevy=new TF1("myLevy",myLevyPtFunc,0,8.5,3);
  myLevy->SetParName(0,"dN/dy");
  myLevy->SetParName(1,"C");
  myLevy->SetParName(2,"n");
  myLevy->SetParameter(0,.008);
  myLevy->SetParameter(1,.3);
  myLevy->SetParameter(2,15);
  myLevy->SetParLimits(0,.001,.2);
  myLevy->SetParLimits(1,.1,1);
  myLevy->SetParLimits(2,1,500);
  myLevy->SetLineColor(1);
  TF1 *myPtExp = new TF1("myPtExp",myPtExpFunc,0,8.5,2);
  myPtExp->SetParName(0,"dN/dy");
  myPtExp->SetParName(1,"C");
  myPtExp->SetParameter(0,.003);
  myPtExp->SetParameter(1,.3);
  myPtExp->SetParLimits(0,0.0001,0.1);
  myPtExp->SetParLimits(1,.1,1);
  myPtExp->SetLineColor(2);
  TF1 *myMtExp = new TF1("myMtExp",myMtExpFunc,0,8.5,2);
  myMtExp->SetParName(0,"dN/dy");
  myMtExp->SetParName(1,"C");
  myMtExp->SetParameter(0,.003);
  myMtExp->SetParameter(1,.3);
  myMtExp->SetParLimits(0,0.0001,0.1);
  myMtExp->SetParLimits(1,.1,1);
  myMtExp->SetLineColor(4);
  TF1 *myBoltz = new TF1("myBoltz",myBoltzFunc,0,8.5,2);
  myBoltz->SetParName(0,"dN/dy");
  myBoltz->SetParName(1,"C");
  myBoltz->SetParameter(0,.003);
  myBoltz->SetParameter(1,.3);
  myBoltz->SetParLimits(0,0.0001,0.1);
  myBoltz->SetParLimits(1,.1,1);
  myBoltz->SetLineColor(6);
  TF1 *myPowerLawPt=new TF1("myPowerLawPt",myPowerLawPtFunc,0,8.5,3);
  myPowerLawPt->SetParName(0,"dN/dy");
  myPowerLawPt->SetParName(1,"pt0");
  myPowerLawPt->SetParName(2,"n");
  myPowerLawPt->SetParameter(0,.008);
  myPowerLawPt->SetParameter(1,50);
  myPowerLawPt->SetParameter(2,50);
  myPowerLawPt->SetParLimits(0,.001,.01);
  myPowerLawPt->SetParLimits(1,1,500);
  myPowerLawPt->SetParLimits(2,1,500);
  myPowerLawPt->SetLineColor(7);
  TF1 *myPowerLaw=new TF1("myPowerLaw",myPowerLawFunc,0,8.5,3);
  myPowerLaw->SetParName(0,"dN/dy");
  myPowerLaw->SetParName(1,"pt0");
  myPowerLaw->SetParName(2,"n");
  myPowerLaw->SetParameter(0,.008);
  myPowerLaw->SetParameter(1,50);
  myPowerLaw->SetParameter(2,50);
  myPowerLaw->SetParLimits(0,.001,.01);
  myPowerLaw->SetParLimits(1,1,500);
  myPowerLaw->SetParLimits(2,1,500);
  myPowerLaw->SetLineColor(8);
  TF1 *myBGBlastWave=new TF1("myBGBlastWave",myBoltsBWFunc,0,8.5,5);
  myBGBlastWave->SetParName(0,"mass");
  myBGBlastWave->SetParName(1,"beta");
  myBGBlastWave->SetParName(2,"C");
  myBGBlastWave->SetParName(3,"n");
  myBGBlastWave->SetParName(4,"Norm");
  myBGBlastWave->FixParameter(0,1.5318);
  myBGBlastWave->SetParameter(1,.3);
  myBGBlastWave->SetParameter(2,.3);
  myBGBlastWave->SetParameter(3,5);
  myBGBlastWave->SetParameter(4,.3);
  myBGBlastWave->SetParLimits(1,.001,.9);
  myBGBlastWave->SetParLimits(2,.1,200);
  myBGBlastWave->SetParLimits(3,.1,200);
  myBGBlastWave->SetLineColor(9);


  FinalSpectrum->Draw();
  legend->AddEntry(FinalSpectrum,"(#Xi(1530) + cc)/2","p");
  //legend->AddEntry(FinalSpectrum,"(#Xi(1530) + cc)/2: Current Results","p");
  FinalSpectrum->Fit(myLevy,"IME","",minFitpoint,maxFitpoint);
  FinalSpectrum->Fit(myPtExp,"IMEN","",minFitpoint,maxFitpoint);
  FinalSpectrum->Fit(myMtExp,"IMEN","",minFitpoint,maxFitpoint);
  FinalSpectrum->Fit(myBoltz,"IMEN","",minFitpoint,maxFitpoint);
  FinalSpectrum->Fit(myPowerLaw,"IMEN","",minFitpoint,maxFitpoint);
  FinalSpectrum->Fit(myPowerLawPt,"IMEN","",minFitpoint,maxFitpoint);
  //FinalSpectrum->Fit(myBGBlastWave,"IME","",minFitpoint,maxFitpoint);
  myLevy->Draw("same");
  /*myPtExp->Draw("same");
  myMtExp->Draw("same");
  myBoltz->Draw("same");
  myPowerLaw->Draw("same");
  myPowerLawPt->Draw("same");*/
  //myBGBlastWave->Draw("same");
  //
  legend->AddEntry(myLevy,"Levy","l");
  /*legend->AddEntry(myPtExp,"Pt-exp","l");
  legend->AddEntry(myMtExp,"Mt-exp","l");
  legend->AddEntry(myBoltz,"Pt-Boltzmann","l");
  legend->AddEntry(myPowerLaw,"PowerLaw","l");
  legend->AddEntry(myPowerLawPt,"Pt-PowerLaw","l");
  */
  double AvgLowPtExtrap = myLevy->Integral(0,0.8) + myPtExp->Integral(0,0.8) + myPowerLawPt->Integral(0,0.8) + myMtExp->Integral(0,0.8);
  AvgLowPtExtrap /= 4.;
  cout<<"Levy Extrap % = "<<myLevy->Integral(0,0.8)/myLevy->Integral(0,10.)<<"  Avg Extrap % = "<<AvgLowPtExtrap/myLevy->Integral(0,10.)<<endl;
  //double LowPtExtrapUncertaintyPercentage = 0.17;// Levy underestimation of Pythia
  //LowPtExtrapUncertaintyPercentage -= (AvgLowPtExtrap - myLevy->Integral(0,0.8))/myLevy->Integral(0,10.);
  double High_lowptExtrap = (myPtExp->Integral(0,0.8) + myPowerLawPt->Integral(0,0.8))/2.;
  double Low_lowptExtrap = myMtExp->Integral(0,0.8);
  double High_LowPtExtrapUncertaintyPercentage = (High_lowptExtrap - myLevy->Integral(0,0.8))/myLevy->Integral(0,10.);
  double Low_LowPtExtrapUncertaintyPercentage = (myLevy->Integral(0,0.8) - Low_lowptExtrap)/myLevy->Integral(0,10.);
  cout<<"high low-pt uncertainty = "<<High_LowPtExtrapUncertaintyPercentage<<endl;
  cout<<"low low-pt uncertainty = "<<Low_LowPtExtrapUncertaintyPercentage<<endl;

  
  //PureSysFromCutVar->Draw();


  double dNdY = myLevy->Integral(0,10);
  double dNdY_exp = myMtExp->Integral(0,10);
  double dNdY_Stat = LevyFits_CutVar[0]->GetParError(0);
  double dNdY_Covered = myLevy->Integral(0.8,5.6);
  double dNdY_Covered_exp = myMtExp->Integral(0.8,5.6);
  double dNdY_DataPoints = 0, dNdY_DataPoints_e=0;
  for(int ii=0; ii<8; ii++){// to explicitly integrate dNdY without fit (Christina's suggestion)
    dNdY_DataPoints += Spectra_CutVar[0]->GetBinContent(ii+2)*Spectra_CutVar[0]->GetBinWidth(ii+2);
    dNdY_DataPoints_e += pow(Spectra_CutVar[0]->GetBinError(ii+2)*Spectra_CutVar[0]->GetBinWidth(ii+2),2);
  }
  dNdY_DataPoints_e = sqrt(dNdY_DataPoints_e);
  double dNdYSysPlus = 0;
  double dNdYSysMinus = 0;
  double CSys = 0;
  double nSys = 0;
  //////////////////////////////////////
  // Sys errors of Levy parameters
  double sigma_dNdY_cv[12]={0};
  double sigma_C_cv[12]={0};
  double sigma_n_cv[12]={0};
  double SysVar_dNdY_cv[12]={0};
  double SysVar_C_cv[12]={0};
  double SysVar_n_cv[12]={0};
  
  for(int cv=1; cv<13; cv++){
    if(MeanSysCutVar[cv-1] > 0 && MeanSysCutVar[cv-1] > sigmaSysCutVar[cv-1]) {// significant deviation
      sigma_dNdY_cv[cv-1] = sqrt(fabs(pow(LevyFits_CutVar[0]->GetParError(0),2) - pow(LevyFits_CutVar[cv]->GetParError(0),2)));
      sigma_C_cv[cv-1] = sqrt(fabs(pow(LevyFits_CutVar[0]->GetParError(1),2) - pow(LevyFits_CutVar[cv]->GetParError(1),2)));
      sigma_n_cv[cv-1] = sqrt(fabs(pow(LevyFits_CutVar[0]->GetParError(2),2) - pow(LevyFits_CutVar[cv]->GetParError(2),2)));
      SysVar_dNdY_cv[cv-1] = fabs(LevyFits_CutVar[0]->GetParameter(0)-LevyFits_CutVar[cv]->GetParameter(0));
      SysVar_C_cv[cv-1] = fabs(LevyFits_CutVar[0]->GetParameter(1)-LevyFits_CutVar[cv]->GetParameter(1));
      SysVar_n_cv[cv-1] = fabs(LevyFits_CutVar[0]->GetParameter(2)-LevyFits_CutVar[cv]->GetParameter(2));
   
      if(SysVar_dNdY_cv[cv-1] > sigma_dNdY_cv[cv-1]){
	dNdYSysPlus += pow(SysVar_dNdY_cv[cv-1],2)-pow(sigma_dNdY_cv[i-1],2);
	dNdYSysMinus += pow(SysVar_dNdY_cv[cv-1],2)-pow(sigma_dNdY_cv[i-1],2);
      }
      if(SysVar_C_cv[cv-1] > sigma_C_cv[cv-1]){
	CSys += pow(SysVar_C_cv[cv-1],2)-pow(sigma_C_cv[cv-1],2);
      }
      if(SysVar_n_cv[cv-1] > sigma_n_cv[cv-1]){
	nSys += pow(SysVar_n_cv[cv-1],2)-pow(sigma_n_cv[cv-1],2);
      }
   
    }
   
  }
  
 

  double sigma_dNdY_fv[4]={0};
  double sigma_C_fv[4]={0};
  double sigma_n_fv[4]={0};
  double SysVar_dNdY_fv[4]={0};
  double SysVar_C_fv[4]={0};
  double SysVar_n_fv[4]={0};
  for(int fv=0; fv<4; fv++){
    sigma_dNdY_fv[fv] = sqrt(fabs(pow(LevyFits_CutVar[0]->GetParError(0),2) - pow(LevyFits_FitVar[fv]->GetParError(0),2)));
    SysVar_dNdY_fv[fv] = sqrt(pow(LevyFits_CutVar[0]->GetParameter(0)-LevyFits_FitVar[fv]->GetParameter(0),2));
    sigma_C_fv[fv] = sqrt(fabs(pow(LevyFits_CutVar[0]->GetParError(1),2) - pow(LevyFits_FitVar[fv]->GetParError(1),2)));
    SysVar_C_fv[fv] = sqrt(pow(LevyFits_CutVar[0]->GetParameter(1)-LevyFits_FitVar[fv]->GetParameter(1),2));
    sigma_n_fv[fv] = sqrt(fabs(pow(LevyFits_CutVar[0]->GetParError(2),2) - pow(LevyFits_FitVar[fv]->GetParError(2),2)));
    SysVar_n_fv[fv] = sqrt(pow(LevyFits_CutVar[0]->GetParameter(2)-LevyFits_FitVar[fv]->GetParameter(2),2));
  
    if(fv==2){//max case for dN/dy
      //dNdYSysPlus += pow(LevyFits_CutVar[0]->GetParameter(0)-LevyFits_FitVar[fv]->GetParameter(0),2) - pow(LevyFits_CutVar[0]->GetParError(0)-LevyFits_FitVar[fv]->GetParError(0),2);
      //dNdYSysMinus += pow(LevyFits_CutVar[0]->GetParameter(0)-LevyFits_FitVar[fv]->GetParameter(0),2) - pow(LevyFits_CutVar[0]->GetParError(0)-LevyFits_FitVar[fv]->GetParError(0),2);
      dNdYSysPlus += pow(SysVar_dNdY_fv[fv],2) - pow(sigma_dNdY_fv[fv],2);
      dNdYSysMinus += pow(SysVar_dNdY_fv[fv],2) - pow(sigma_dNdY_fv[fv],2);
    }
    //if(fv==2) CSys += pow(LevyFits_CutVar[0]->GetParameter(1)-LevyFits_FitVar[fv]->GetParameter(1),2) - pow(LevyFits_CutVar[0]->GetParError(1)-LevyFits_FitVar[fv]->GetParError(1),2);
    //if(fv==1) nSys += pow(LevyFits_CutVar[0]->GetParameter(2)-LevyFits_FitVar[fv]->GetParameter(2),2) - pow(LevyFits_CutVar[0]->GetParError(2)-LevyFits_FitVar[fv]->GetParError(2),2);
    if(fv==2) CSys += pow(SysVar_C_fv[fv],2) - pow(sigma_C_fv[fv],2);
    if(fv==1) nSys += pow(SysVar_n_fv[fv],2) - pow(sigma_n_fv[fv],2);
  }
  
  
  double dNdYSys_forRatio = dNdYSysPlus;
  dNdYSysPlus += pow(0.062*dNdY,2);// INEL Norm uncertainty
  dNdYSysPlus += pow(0.04*dNdY,2);// Material Budget uncertainty
  dNdYSysPlus += pow(0.01*dNdY,2);// Geant3/Fluka correction uncertainty
  
  dNdYSysPlus += pow(High_LowPtExtrapUncertaintyPercentage * dNdY,2);// extrapolation uncertainty (was LowPtExtrapUncertaintyPercentage * fabs(dNdY-dNdY_Covered))
  dNdYSys_forRatio += pow(High_LowPtExtrapUncertaintyPercentage * dNdY,2);// extrapolation uncertainty 
  dNdYSysPlus = sqrt(dNdYSysPlus);
  dNdYSys_forRatio = sqrt(dNdYSys_forRatio);
  //
  dNdYSysMinus += pow(0.03*dNdY,2);// INEL Norm uncertainty
  dNdYSysMinus += pow(0.04*dNdY,2);// Material Budget uncertainty
  dNdYSysMinus += pow(0.01*dNdY,2);// Geant3/Fluka correction uncertainty
  dNdYSysMinus += pow(Low_LowPtExtrapUncertaintyPercentage * dNdY,2);// extrapolation uncertainty (40% from Pythia study).
  dNdYSysMinus = sqrt(dNdYSysMinus);

  CSys = sqrt(CSys);
  nSys = sqrt(nSys);

  cout<<"dN/dy = "<<dNdY<<"  , Covered range = "<<dNdY_Covered<<"  , DataPoint Sum = "<<dNdY_DataPoints<<"  , stat = "<<dNdY_DataPoints_e<<endl;
  cout<<"dN/dy stat = "<<dNdY_Stat<<"   +sys = "<<dNdYSysPlus<<"  -sys = "<<dNdYSysMinus<<endl;
    
  cout<<"C = "<<myLevy->GetParameter(1)<<"  , stat = "<<LevyFits_CutVar[0]->GetParError(1)<<"  , sys = "<<CSys<<endl;
  cout<<"n = "<<myLevy->GetParameter(2)<<"  , stat = "<<LevyFits_CutVar[0]->GetParError(2)<<"  , sys = "<<nSys<<endl;

  // mean pt
  double ptEdges[9]={0.8, 1.2, 1.6, 2.0, 2.4, 3.2, 4.0, 4.8, 5.6};
  double meanpt_def=0, meanpt_Sys1=0, meanpt_Sys2=0, meanpt_Sys3=0;// Sys1 = Spectrum only, Sys2 = bin centers, Sys3 = low-pt extrapolation
  double yieldSum_def=0, yieldSum_Sys1=0, yieldSum_Sys2=0, yieldSum_Sys3=0;
  double meanpt_statError=0;
  meanpt_def += myLevy->Moment(1,0.,0.8) * myLevy->Integral(0,0.8);// low-pt extrapolation
  yieldSum_def += myLevy->Integral(0,0.8);
  meanpt_Sys1 += meanpt_def; meanpt_Sys2 += meanpt_def; 
  yieldSum_Sys1 += yieldSum_def; yieldSum_Sys2 += yieldSum_def; 
  meanpt_Sys3 += fabs(myLevy->Integral(0,0.8) - High_LowPtExtrapUncertaintyPercentage*dNdY) * myLevy->Moment(1,0.,0.8);// low-pt extrapolation uncertainty
  yieldSum_Sys3 += fabs(myLevy->Integral(0,0.8) - High_LowPtExtrapUncertaintyPercentage*dNdY);// low-pt extrapolation uncertainty
  meanpt_statError += myLevy->Moment(1,0.,0.8) * myLevy->Integral(0,0.8) * myLevy->GetParError(0)/myLevy->GetParameter(0);
  

  double StatErrorPtBins=0;
  for(int bin=1; bin<=8; bin++){
    double bin_width = fabs(ptEdges[bin-1]-ptEdges[bin]);
    meanpt_def += myLevy->Moment(1,ptEdges[bin-1],ptEdges[bin]) * FinalSpectrum->GetBinContent(bin+1)*bin_width;
    yieldSum_def += FinalSpectrum->GetBinContent(bin+1)*bin_width;
    meanpt_Sys1 += myLevy->Moment(1,ptEdges[bin-1],ptEdges[bin]) * myLevy->Integral(ptEdges[bin-1],ptEdges[bin]);
    yieldSum_Sys1 += myLevy->Integral(ptEdges[bin-1],ptEdges[bin]);
    meanpt_Sys2 += (ptEdges[bin]+ptEdges[bin-1])/2. * FinalSpectrum->GetBinContent(bin+1)*bin_width;
    yieldSum_Sys2 += FinalSpectrum->GetBinContent(bin+1)*bin_width;
    meanpt_Sys3 += myLevy->Moment(1,ptEdges[bin-1],ptEdges[bin]) * FinalSpectrum->GetBinContent(bin+1)*bin_width;
    yieldSum_Sys3 += FinalSpectrum->GetBinContent(bin+1)*bin_width;
    //
    StatErrorPtBins += pow(myLevy->Moment(1,ptEdges[bin-1],ptEdges[bin]) * Spectra_CutVar[0]->GetBinError(bin+1)*bin_width,2);
  }
  
  meanpt_def += myLevy->Moment(1,5.6,10.) * myLevy->Integral(5.6,10.);// high-pt extrapolation
  yieldSum_def += myLevy->Integral(5.6,10.);
  meanpt_Sys1 += myLevy->Moment(1,5.6,10.) * myLevy->Integral(5.6,10.);// high-pt extrapolation
  yieldSum_Sys1 += myLevy->Integral(5.6,10.);
  meanpt_Sys2 += myLevy->Moment(1,5.6,10.) * myLevy->Integral(5.6,10.);// high-pt extrapolation
  yieldSum_Sys2 += myLevy->Integral(5.6,10.);
  meanpt_Sys3 += fabs(myLevy->Integral(5.6,10.) - High_LowPtExtrapUncertaintyPercentage*myLevy->Integral(5.6,10.))*myLevy->Moment(1,5.6,10.);// high-pt extrapolation
  yieldSum_Sys3 += fabs(myLevy->Integral(5.6,10.) - High_LowPtExtrapUncertaintyPercentage*myLevy->Integral(5.6,10.));
  //
  meanpt_statError += sqrt(StatErrorPtBins);
  meanpt_statError += myLevy->Moment(1,5.6,10.) * myLevy->Integral(5.6,10.) * myLevy->GetParError(0)/myLevy->GetParameter(0);// high-pt extrapolation
  meanpt_statError /= yieldSum_def;

  meanpt_def /= yieldSum_def;
  meanpt_Sys1 /= yieldSum_Sys1;
  meanpt_Sys2 /= yieldSum_Sys2;
  meanpt_Sys3 /= yieldSum_Sys3;

  double meanpt_Sys = sqrt(pow(meanpt_def-meanpt_Sys1,2)+pow(meanpt_def-meanpt_Sys2,2)+pow(meanpt_def-meanpt_Sys3,2));
  
  cout<<"<pt> = "<<meanpt_def<<"  , stat = "<<meanpt_statError<<"  , sys = "<<meanpt_Sys<<endl;

  /*double meanpt = myLevy->Moment(1,0.,10.);// 0.,30 for full range
  double meanpt_CoveredLevy = myLevy->Moment(1,0.8,5.6);
  double mPt_Stat = 0;
  double mPt_Sys = 0;
  double mPt_StatSys = 0;
  double base_C = myLevy->GetParameter(1);
  double base_n = myLevy->GetParameter(2);
  
  
  // Sys
  myLevy->SetParameter(1, base_C + myLevy->GetParError(1));
  mPt_StatSys += pow(myLevy->Moment(1,0.,10)-meanpt,2);
  //
  //
  myLevy->SetParameter(2, base_n + myLevy->GetParError(2));
  mPt_StatSys += pow(myLevy->Moment(1,0.,10)-meanpt,2);
  mPt_StatSys = sqrt(mPt_StatSys);
  //
  //
  // Stat
  myLevy->SetParameter(1, base_C + LevyFits_CutVar[0]->GetParError(1));
  mPt_Stat += pow(myLevy->Moment(1,0.,10)-meanpt,2);
  //
  //
  myLevy->SetParameter(2, base_n + LevyFits_CutVar[0]->GetParError(2));
  mPt_Stat += pow(myLevy->Moment(1,0.,10)-meanpt,2);
  mPt_Stat = sqrt(mPt_Stat);
  //
  myLevy->SetParameter(1, base_C); myLevy->SetParameter(1, base_n);
  
  cout<<"Levy mean pt = "<<meanpt<<"  , stat = "<<mPt_Stat<<"  , sys = "<<sqrt( pow(mPt_StatSys,2)-pow(mPt_Stat,2) )<<endl;
  cout<<"Levy mean pt in covered range = "<<meanpt_CoveredLevy<<endl;
  */

  double XidNdY = 0.0079;// (Xi^+ + Xi^-)/2. Published
  double XiStarToXi = dNdY/XidNdY;
  double XiStarToXi_Stat = sqrt(pow(dNdY_Stat/XidNdY,2) + pow(dNdY*0.0001/pow(XidNdY,2),2));// 0.0001 is stat error of Xi^{\pm}
  double XiStartoXi_Sys = dNdYSys_forRatio/XidNdY;
  
  
  cout<<"Xi*/Xi = "<<XiStarToXi<<"  , stat = "<<XiStarToXi_Stat<<"  , sys = "<<XiStartoXi_Sys<<endl;
  cout<<endl;
  cout<<endl;
 
  //Old QM2011 values for (Xi* + cc)/2
  double Old_XiStarPoints[8]={0.001259, 0.0008586, 0.0006316, 0.0003604, 1.856e-4, 5.4601e-5, 2.0187e-5, 7.6559e-6};
  double Old_XiStarPoints_e[8]={2.44e-4, 1.224e-4, 7.136e-5, 3.85e-5, 1.59e-5, 7.49e-6, 3.51e-6, 1.822e-6};
  double pt_points[8]={1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.4, 5.2};
  double pt_points_e[8]={.2, .2, .2, .2, .4, .4, .4, .4};
  TGraphErrors *gr_OldXiStar = new TGraphErrors(8, pt_points, Old_XiStarPoints, pt_points_e, Old_XiStarPoints_e);
  gr_OldXiStar->SetMarkerStyle(21);
  gr_OldXiStar->SetMarkerColor(2);
  gr_OldXiStar->SetLineColor(2);
  double NewOldRatio[8]={0};
  for(int ii=0; ii<8; ii++){
    NewOldRatio[ii] = FinalSpectrum->GetBinContent(ii+2)/Old_XiStarPoints[ii];
    //cout<<FinalSpectrum->GetBinError(ii+2)/FinalSpectrum->GetBinContent(ii+2)<<endl;
  }
  TGraph *gr_NewOldRatio = new TGraph(8, pt_points, NewOldRatio);
  gr_NewOldRatio->SetMarkerStyle(20);
  gr_NewOldRatio->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  gr_NewOldRatio->GetYaxis()->SetTitle("New/Old");
  gr_NewOldRatio->SetTitle("New Spectrum Divided by QM11 Spectrum");
  //gr_NewOldRatio->Draw("AP");

  //gr_OldXiStar->Draw("P");
  //legend->AddEntry(gr_OldXiStar,"(#Xi(1530) + cc)/2: QM11 Results","p");
  // To Fit Old points
  /*TH1D *h_OldXiStar = (TH1D*)h_Xispectrum->Clone();
  h_OldXiStar->Add(h_OldXiStar,-1);
  for(int i=0; i<8; i++){ 
    h_OldXiStar->SetBinContent(i+2, Old_XiStarPoints[i]);
    h_OldXiStar->SetBinError(i+2, Old_XiStarPoints_e[i]);
  }
  h_OldXiStar->SetBinContent(1,+100);// first bin has no data
  h_OldXiStar->Fit(myLevy,"IME","",minFitpoint,maxFitpoint);
  */
  //legend->Draw("same");

  double StatPercentage[8]={};
  double SysPercentage[8]={};
  double TotErrPercentage[8]={};
  for(int ii=0; ii<8; ii++){
    StatPercentage[ii] = StatError[ii]/FinalSpectrum->GetBinContent(ii+2);
    SysPercentage[ii] = SysError[ii+1]/FinalSpectrum->GetBinContent(ii+2);
    TotErrPercentage[ii] = sqrt(pow(StatError[ii],2)+pow(SysError[ii+1],2))/FinalSpectrum->GetBinContent(ii+2);
  }
  TGraph *gr_StatPercentage = new TGraph(8, pt_points, StatPercentage);
  gr_StatPercentage->SetMarkerStyle(25);
  gr_StatPercentage->SetMarkerColor(4);
  gr_StatPercentage->SetLineColor(4);
  gr_StatPercentage->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  gr_StatPercentage->GetYaxis()->SetTitle("Uncertainty (fraction)");
  gr_StatPercentage->SetTitle("");
  gr_StatPercentage->SetMinimum(0);
  gr_StatPercentage->SetMaximum(0.2);// 0.35 to compare with Enrico
  //gr_StatPercentage->Draw("APC");
  TGraph *gr_SysPercentage = new TGraph(8, pt_points, SysPercentage);
  gr_SysPercentage->SetMarkerStyle(25);
  gr_SysPercentage->SetMarkerColor(2);
  gr_SysPercentage->SetLineColor(2);
  //gr_SysPercentage->Draw("PC");
  TGraph *gr_TotErrPercentage = new TGraph(8, pt_points, TotErrPercentage);
  gr_TotErrPercentage->SetMarkerStyle(25);
  gr_TotErrPercentage->SetMarkerColor(1);
  gr_TotErrPercentage->SetLineColor(1);
  //gr_TotErrPercentage->Draw("PC");
  //legend->AddEntry(gr_StatPercentage,"Stat. uncertainty","p");
  //legend->AddEntry(gr_SysPercentage,"Sys. uncertainty","p");
  //legend->AddEntry(gr_TotErrPercentage,"#sqrt{Stat^{2}+Sys^{2}}","p");
  legend->Draw("same");
  
  ////////////////////////////////
  // Cut Var
  for(int var=1; var<13; var++){// 0 to 4 for FitVar,  1 to 13 for CutVar
    TString *Varname = new TString();
    if(var==1) Varname->Append("Nclusters TPC");
    if(var==2) Varname->Append("DCA PV proton");
    if(var==3) Varname->Append("DCA PV 1st pion");
    if(var==4) Varname->Append("DCA PV 2nd pion");
    if(var==5) Varname->Append("DCA PV 3rd pion");
    if(var==6) Varname->Append("DCA PV 4th pion");
    if(var==7) Varname->Append("DCA proton-pion");
    if(var==8) Varname->Append("DCA Lambda-pion");
    if(var==9) Varname->Append("Decay length xy Lambda");
    if(var==10) Varname->Append("Decay length xy Xi");
    if(var==11) Varname->Append("Cos PointingAngle Lambda");
    if(var==12) Varname->Append("Cos PointingAngle Xi");

    cout<<Varname->Data()<<" & ";
    cout.precision(2);
    for(int par=0; par<3; par++) {
      cout<<100.*(LevyFits_CutVar[var]->GetParameter(par) - LevyFits_CutVar[0]->GetParameter(par))/LevyFits_CutVar[0]->GetParameter(par)<<"(";
      if(par==0) cout<<100.*sigma_dNdY_cv[var-1]/LevyFits_CutVar[0]->GetParameter(par)<<")"<<" ";
      if(par==1) cout<<100.*sigma_C_cv[var-1]/LevyFits_CutVar[0]->GetParameter(par)<<")"<<" ";
      if(par==2) cout<<100.*sigma_n_cv[var-1]/LevyFits_CutVar[0]->GetParameter(par)<<")"<<" ";
      
      if((par+1)!=3) cout<<"& ";
      else cout<<" \\\\\ \\hline";
    }
    cout<<endl;
  }
  cout<<endl;
 
  //////////////////////////////
  // Fit Var
  for(int var=0; var<4; var++){// 0 to 4 for FitVar,  1 to 13 for CutVar
    TString *Varname = new TString();
    if(var==0) Varname->Append("Pure Voigtian");
    if(var==1) Varname->Append("Pure Voigtian w/o bkg subtraction");
    if(var==2) Varname->Append("bkg norm to the right");
    if(var==3) Varname->Append("1/2 bin-counting range");

    cout<<Varname->Data()<<" & ";
    cout.precision(2);
    for(int par=0; par<3; par++) {
      cout<<100.*(LevyFits_FitVar[var]->GetParameter(par) - LevyFits_CutVar[0]->GetParameter(par))/(LevyFits_CutVar[0]->GetParameter(par))<<"(";
      if(par==0) cout<<100.*sigma_dNdY_fv[var]/LevyFits_CutVar[0]->GetParameter(par)<<")"<<" ";
      if(par==1) cout<<100.*sigma_C_fv[var]/LevyFits_CutVar[0]->GetParameter(par)<<")"<<" ";
      if(par==2) cout<<100.*sigma_n_fv[var]/LevyFits_CutVar[0]->GetParameter(par)<<")"<<" ";

      if((par+1)!=3) cout<<"& ";
      else cout<<" \\\\\ \\hline";
    }
    cout<<endl;
  }
  cout<<endl;

  ///////////////////////////////////////////

  /*
  cout<<"Yields"<<endl;
  // print yields
  for(int var=1; var<=13; var++){
    TString *Varname = new TString();
    if(var==1) Varname->Append("Standard ");
    else {Varname->Append("Cut Var "); *Varname += var-1;}
    
    cout<<Varname->Data()<<" & ";
    cout.precision(4);// 4 or 6
    for(int ptbin=1; ptbin<=8; ptbin++){
      //cout<<Yields_CutVar[var-1]->GetBinContent(ptbin)/Eff_CutVar[var-1]->GetBinContent(ptbin)<<" ";
      cout<<Eff_CutVar[var-1]->GetBinContent(ptbin)<<" ";
      if((ptbin)!=8) cout<<"& ";
      else cout<<" \\\\\ \\hline";
    }
    cout<<endl;
  }
  */
  /* for(int var=1; var<5; var++){
    TString *Varname = new TString("Fit Var ");
    *Varname += var;
    cout<<Varname->Data()<<" & ";
    cout.precision(6);
    for(int ptbin=1; ptbin<=8; ptbin++){
      cout<<Yields_FitVar[var-1]->GetBinContent(ptbin)/Eff_CutVar[0]->GetBinContent(ptbin)<<" ";
      if((ptbin)!=8) cout<<"& ";
      else cout<<" \\\\\ \\hline";
    }
    cout<<endl;
    }*/
  cout<<endl;

  //TString *name=new TString("../AnaNotes/XiStar/XiStar_yield_ErrorsCutVar");
  //*name += CutVarN;
  /*TString *name=new TString("../AnaNotes/XiStar/XiStar_CutVarSys");
  *name += CutVarLow;
  name->Append(".png");
  c1->SaveAs(name->Data());*/
  //TString *name=new TString("../AnaNotes/XiStar/XiStar_CutVarSummary");
  //*name += CutVarLow;
  //name->Append(".png");
  //can->SaveAs(name->Data());
  
}


//________________________________________________________________________
double myLevyPtFunc(Double_t *x, Double_t *par)
{
  Double_t lMass=0;
  lMass = 1.5318; //Xi* mass
  

  Double_t ldNdy  = par[0]; // dN/dy
  Double_t l2pi   = 2*TMath::Pi(); // 2pi (cancels in return statement)
  Double_t lTemp = par[1]; // Temperature
  Double_t lPower = par[2]; // power=n
  
  Double_t lBigCoef = ((lPower-1)*(lPower-2)) / (l2pi*lPower*lTemp*(lPower*lTemp+lMass*(lPower-2)));
  Double_t lInPower = 1 + (TMath::Sqrt(x[0]*x[0]+lMass*lMass)-lMass) / (lPower*lTemp);
  
  return l2pi * ldNdy * x[0] * lBigCoef * TMath::Power(lInPower,(-1)*lPower);
}
//________________________________________________________________________
double myMtExpFunc(Double_t *x, Double_t *par)
{
  Double_t lMass=0;
  lMass = 1.5318; //Xi* mass
  
  return 2*TMath::Pi()*x[0]*par[0]*exp(-TMath::Sqrt(x[0]*x[0]+lMass*lMass)/par[1]);
}
//________________________________________________________________________
double myPtExpFunc(Double_t *x, Double_t *par)
{
  return 2*TMath::Pi()*x[0]*par[0]*exp(-x[0]/par[1]);
}
//________________________________________________________________________
double myBoltzFunc(Double_t *x, Double_t *par)
{
  Double_t lMass=0;
  lMass = 1.5318; //Xi* mass
  Double_t Mt = sqrt(pow(x[0],2)+pow(lMass,2));
  return par[0]*x[0]*Mt*exp(-Mt/par[1]);
  
}
//________________________________________________________________________
double myPowerLawPtFunc(Double_t *x, Double_t *par)
{
  return ( x[0]*par[0]*pow(1+x[0]/par[1], -par[2]) );
}
//________________________________________________________________________
double myPowerLawFunc(Double_t *x, Double_t *par)
{
  return ( par[0]*pow(1+x[0]/par[1], -par[2]) );
}
//________________________________________________________________________
double myBoltsBWFunc(Double_t *x, Double_t *par)
{
  double pT = x[0];
  double mass    = par[0];
  double beta    = par[1];
  double temp    = par[2];
  double n       = par[3];
  double norm    = par[4];

  static TF1 * fIntBG = 0;
  if(!fIntBG)
    fIntBG = new TF1 ("fIntBG", IntegrandBG, 0, 1, 5);

  fIntBG->SetParameters(mass, pT, beta, temp,n);
  double result = fIntBG->Integral(0,1);
  return result*norm;
}
//________________________________________________________________________
double IntegrandBG(const double * x, const double* par){
  double x0 = x[0];
 
  double mass     = par[0];
  double pT       = par[1];
  double beta_max = par[2];
  double temp     = par[3];
  Double_t n      = par[4];

  // Keep beta within reasonable limits
  Double_t beta = beta_max * TMath::Power(x0, n);
  if (beta > 0.9999999999999999) beta = 0.9999999999999999;

  double mT      = TMath::Sqrt(mass*mass+pT*pT);

  double rho0   = TMath::ATanH(beta);  
  double arg00 = pT*TMath::SinH(rho0)/temp;
  if (arg00 > 700.) arg00 = 700.; // avoid FPE
  double arg01 = mT*TMath::CosH(rho0)/temp;
  double f0 = x0*mT*TMath::BesselI0(arg00)*TMath::BesselK1(arg01);

  //  printf("r=%f, pt=%f, beta_max=%f, temp=%f, n=%f, mt=%f, beta=%f, rho=%f, argI0=%f, argK1=%f\n", x0, pT, beta_max, temp, n, mT, beta, rho0, arg00, arg01);

  return f0;
}
