/*************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

///////////////////////////////////////////////////////////////////////////
//                Dielectron SignalFunc                                  //
//                                                                       //
//                                                                       //
/*
Dielectron signal extraction class using functions as input.

A function to describe the signal as well as one to describe the background
has to be deployed by the user. Alternatively on of the default implementaions
can be used.

*/
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <TF1.h>
#include <TH1.h>
#include <TGraph.h>
#include <TMath.h>
#include <TString.h>
#include <TPaveText.h>
#include <RooRealVar.h>
#include <RooCBShape.h>
#include <RooExponential.h>

#include <AliLog.h>

#include "AliDielectronSignalFunc.h"

ClassImp(AliDielectronSignalFunc)

AliDielectronSignalFunc::AliDielectronSignalFunc() :
  AliDielectronSignalBase(),
  fSignal(0x0),
  fBackground(0x0),
  fSigBack(0x0),
  fFitFunc(0x0),
  fVInitParam(0),
  fFitOpt("MNQE"),
  fUseIntegral(kFALSE),
  fFitMin(2.5),
  fFitMax(4),
  fParM(-1),
  fParMres(-1)
{
  //
  // Default Constructor
  //
  
}

//______________________________________________
AliDielectronSignalFunc::AliDielectronSignalFunc(const char* name, const char* title) :
  AliDielectronSignalBase(name, title),
  fSignal(0x0),
  fBackground(0x0),
  fSigBack(0x0),
  fFitFunc(0x0),
  fVInitParam(0),
  fFitOpt("MNQ"),
  fUseIntegral(kFALSE),
  fFitMin(2.5),
  fFitMax(4),
  fParM(-1),
  fParMres(-1)
{
  //
  // Named Constructor
  //
}

//______________________________________________
AliDielectronSignalFunc::~AliDielectronSignalFunc()
{
  //
  // Default Destructor
  //
  if (fSigBack) delete fSigBack;
}


//______________________________________________
void AliDielectronSignalFunc::Process(TH1 * const hist)
{
  //
  // Fit the mass histogram with the function and retrieve the parameters
  //
  
  //reset result arrays
  Reset();

  //sanity check
  if (!fSigBack) {
    AliError("Use 'SetFunctions(TF1*,TF1*)' or 'SetDefaults(Int_t)' to setup the signal functions first'!");
    return;
  }
  
  //initial the fit function to its default parameters
  if (fVInitParam.GetMatrixArray()) fSigBack->SetParameters(fVInitParam.GetMatrixArray());

  //Perform fit and retrieve values
  hist->Fit(fSigBack,fFitOpt.Data(),"",fFitMin,fFitMax);

  //retrieve values and errors
  //TODO: perhpas implement different methods to retrieve the valus
  fSignal->SetParameters(fSigBack->GetParameters());
  fBackground->SetParameters(fSigBack->GetParameters()+fSignal->GetNpar());

  //TODO: proper error estimate?!?
  //      perhaps this is not possible in a general way?
  
  // signal and background histograms
  TH1 *hSignal=(TH1*)hist->Clone("DieleSignalHist");
  TH1 *hBackground=(TH1*)hist->Clone("DieleBackgroundHist");

  fSignal->SetRange(hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
  fBackground->SetRange(hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
  
  hSignal->Add(fBackground,-1);
  hBackground->Add(fSignal,-1);

  
  //get values and errors signal and background
  Double_t signal=0, signal_err=0;
  Double_t background=0, background_err=0;
  Double_t sigPlusBack=0, sigPlusBack_err=0;

  if (!fUseIntegral){
    signal=hSignal->IntegralAndError(hSignal->FindBin(GetIntegralMin()),
                                     hSignal->FindBin(GetIntegralMax()), signal_err);
    background=hBackground->IntegralAndError(hBackground->FindBin(GetIntegralMin()),
                                             hBackground->FindBin(GetIntegralMax()), background_err);
    sigPlusBack=hist->IntegralAndError(hist->FindBin(GetIntegralMin()),
                                       hist->FindBin(GetIntegralMax()), sigPlusBack_err);
  } else {
    signal=fSignal->Integral(GetIntegralMin(),GetIntegralMax())/hist->GetBinWidth(hist->FindBin((GetIntegralMax()-GetIntegralMin())/2.));
    background=fBackground->Integral(GetIntegralMin(),GetIntegralMax())/hist->GetBinWidth(hist->FindBin((GetIntegralMax()-GetIntegralMin())/2.));
    sigPlusBack=fSigBack->Integral(GetIntegralMin(),GetIntegralMax())/hist->GetBinWidth(hist->FindBin((GetIntegralMax()-GetIntegralMin())/2.));
    //TODO: Error calculation
  }
  
  //set values
  SetSignal(signal,signal_err);
  SetBackground(background,background_err);
  SetSignificanceAndSOB();
  
  //cleanup
  delete hSignal;
  delete hBackground;
}

//______________________________________________
void AliDielectronSignalFunc::Process(TObjArray * const arrhist)
{
  //
  // Fit the mass histogram with the function and retrieve the parameters
  //
  
  TH1 *hist=(TH1*)arrhist->At(1);
  Process(hist);
}
//______________________________________________
void AliDielectronSignalFunc::SetFunctions(TF1 * const sig, TF1 * const back, TF1 * const combined,
                                           Int_t parM, Int_t parMres)
{
  //
  // Set the signal, background functions and combined fit function
  // Note: The process method assumes that the first n parameters in the
  //       combined fit function correspond to the n parameters of the signal function
  //       and the n+1 to n+m parameters to the m parameters of the background function!!!
  // if parM and parMres are set this information can be used in the drawing function
  //   to also show in the statistics box the mass and mass resolution

  if (!sig||!back||!combined) {
    AliError("Both, signal and background function need to be set!");
  }
  fSignal=sig;
  fBackground=back;
  fSigBack=combined;

  InitParams();
  fParM=parM;
  fParMres=parMres;
}

//______________________________________________
void AliDielectronSignalFunc::InitParams()
{
  //
  // initilise common fit parameters
  //
  fVInitParam.ResizeTo(fSigBack->GetNpar());
  fVInitParam.SetElements(fSigBack->GetParameters());
}

//______________________________________________
void AliDielectronSignalFunc::SetDefaults(Int_t type)
{
  //
  // Setup some default functions:
  // type = 0: gaus signal + linear background in 2.5 - 4 GeV inv. mass
  // type = 1: gaus signal + exponential background in 2.5 - 4 GeV inv. mass
  // type = 2: half gaussian, half exponential signal function
  // type = 3: Crystal-Ball function
  // type = 4: Crystal-Ball signal + exponential background
  //

  
  if (type==0){
    fSignal=new TF1("DieleSignal","gaus",2.5,4);
    fBackground=new TF1("DieleBackground","pol1",2.5,4);
    fSigBack=new TF1("DieleCombined","gaus+pol1(3)",2.5,4);
    
    fSigBack->SetParameters(1,3.1,.05,2.5,1);
    fSigBack->SetParLimits(0,0,10000000);
    fSigBack->SetParLimits(1,3.05,3.15);
    fSigBack->SetParLimits(2,.02,.1);
    
    SetFunctions(fSignal,fBackground,fSigBack,1,2);

  } else if (type==1){
    fSignal=new TF1("DieleSignal","gaus",2.5,4);
    fBackground=new TF1("DieleBackground","[0]*exp(-(x-[1])/[2])",2.5,4);
    fSigBack=new TF1("DieleCombined","gaus+[3]*exp(-(x-[4])/[5])",2.5,4);
    
    fSigBack->SetParameters(1,3.1,.05,1,2.5,1);
    fSigBack->SetParLimits(0,0,10000000);
    fSigBack->SetParLimits(1,3.05,3.15);
    fSigBack->SetParLimits(2,.02,.1);
    
    SetFunctions(fSignal,fBackground,fSigBack,1,2);

  } else if (type==2){
    // half gaussian, half exponential signal function
    // exponential background
    fSignal = new     TF1("DieleSignal","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",2.5,4);
    fBackground=new TF1("DieleBackground","[0]*exp(-(x-[1])/[2])+[3]",2.5,4);
    fSigBack=new TF1("DieleCombined","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))+[4]*exp(-(x-[5])/[6])+[7]",2.5,4);
    fSigBack->SetParameters(1.,3.1,.05,.1,1,2.5,1,0);
    
    fSigBack->SetParLimits(0,0,10000000);
    fSigBack->SetParLimits(1,3.05,3.15);
    fSigBack->SetParLimits(2,.02,.1);
    fSigBack->FixParameter(6,2.5);
    fSigBack->FixParameter(7,0);
    SetFunctions(fSignal,fBackground,fSigBack,1,2);

  } else if (type==3){
    // Crystal-Ball function
    RooRealVar m("m","Invariant mass",2.5,3.8);
    RooRealVar alpha("alpha","alpha",0.5,0,10);
    RooRealVar n("n","n",1,0,10);
    RooRealVar m0("m0","m0",3.1,1,5);
    RooRealVar sigma("sigma","sigma",0.3,0,3);
    RooCBShape jpsi("jpsi","jpsi",m,m0,sigma,alpha,n);

    //RooAddPdf model("model","jpsi fitting function", RooArgList(jpsi,expo), fsig);
    RooRealVar fsig("fsig","signal fraction",0.5,0.,1.);
    RooAddPdf model("model","jpsi fitting function", RooArgList(jpsi),fsig);
 
    model.Clone("fFitFunc");

  } else if (type==4){
    // Crystal-Ball function
    RooRealVar m("m","Invariant mass",2.5,3.8);
    RooRealVar alpha("alpha","alpha",0.5,0,10);
    RooRealVar n("n","n",1,0,10);
    RooRealVar m0("m0","m0",3.1,1,5);
    RooRealVar sigma("sigma","sigma",0.3,0,3);
    RooCBShape jpsi("jpsi","jpsi",m,m0,sigma,alpha,n);

    // Expoenetial function
    RooRealVar al("al","al",0.5);             
    RooExponential expo("expo","exponential", m, al);

    // add two functions
    RooRealVar fsig("fsig","signal fraction",0.5,0.,1.);
    RooAddPdf model("model","jpsi fitting function", RooArgList(jpsi,expo), fsig);

    model.Clone("fFitFunc");
  }
}


//______________________________________________
void AliDielectronSignalFunc::Draw(const Option_t* option)
{
  //
  // Draw the fitted function
  //

  TString drawOpt(option);
  drawOpt.ToLower();

  Bool_t optStat=drawOpt.Contains("stat");
  
  fSigBack->SetNpx(200);
  fSigBack->SetRange(GetIntegralMin(),GetIntegralMax());
  fBackground->SetNpx(200);
  fBackground->SetRange(GetIntegralMin(),GetIntegralMax());
  
  TGraph *grSig=new TGraph(fSigBack);
  grSig->SetFillColor(kGreen);
  grSig->SetFillStyle(3001);

  TGraph *grBack=new TGraph(fBackground);
  grBack->SetFillColor(kRed);
  grBack->SetFillStyle(3001);

  grSig->SetPoint(0,grBack->GetX()[0],grBack->GetY()[0]);
  grSig->SetPoint(grSig->GetN()-1,grBack->GetX()[grBack->GetN()-1],grBack->GetY()[grBack->GetN()-1]);
  
  grBack->SetPoint(0,grBack->GetX()[0],0.);
  grBack->SetPoint(grBack->GetN()-1,grBack->GetX()[grBack->GetN()-1],0.);
  
  fSigBack->SetRange(fFitMin,fFitMax);
  fBackground->SetRange(fFitMin,fFitMax);
  
  if (!drawOpt.Contains("same")){
    grSig->Draw("af");
  } else {
    grSig->Draw("f");
  }
  grBack->Draw("f");
  fSigBack->Draw("same");
  fBackground->Draw("same");

  if (optStat) {
    TPaveText* pave=DrawStats();
    if (fParM>=0) pave->AddText(Form("Mass: %.3f #pm %.3f", fSigBack->GetParameter(fParM), fSigBack->GetParError(fParM)));
    if (fParMres>=0) pave->AddText(Form("Mass res.: %.3f #pm %.3f", fSigBack->GetParameter(fParMres), fSigBack->GetParError(fParMres)));
  }
}

