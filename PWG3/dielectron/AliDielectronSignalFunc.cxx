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

/* $Id$ */

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
#include <TList.h>
#include <TFitResult.h>
//#include <../hist/hist/src/TF1Helper.h>

#include <AliLog.h>

#include "AliDielectronSignalFunc.h"

ClassImp(AliDielectronSignalFunc)

AliDielectronSignalFunc::AliDielectronSignalFunc() :
AliDielectronSignalBase(),
fFuncSignal(0x0),
fFuncBackground(0x0),
fFuncSigBack(0x0),
fParMass(1),
fParMassWidth(2),
fFitOpt("SMNQE"),
fUseIntegral(kFALSE)
{
  //
  // Default Constructor
  //
  
}

//______________________________________________
AliDielectronSignalFunc::AliDielectronSignalFunc(const char* name, const char* title) :
AliDielectronSignalBase(name, title),
fFuncSignal(0x0),
fFuncBackground(0x0),
fFuncSigBack(0x0),
fParMass(1),
fParMassWidth(2),
fFitOpt("SMNQE"),
fUseIntegral(kFALSE)
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
  if(fFuncSignal) delete fFuncSignal;
  if(fFuncBackground) delete fFuncBackground;
  if(fFuncSigBack) delete fFuncSigBack;
}


//______________________________________________
void AliDielectronSignalFunc::Process(TObjArray * const arrhist)
{
  //
  // Fit the invariant mass histograms and retrieve the signal and background
  //
  switch(fMethod) {
  case kFitted :
    ProcessFit(arrhist);
    break;
    
  case kLikeSign :
    ProcessLS(arrhist);
    break;
    
  case kEventMixing :
    ProcessEM(arrhist);
    break;
    
  default :
    AliError("Background substraction method not known!");
  }
}


//______________________________________________
void AliDielectronSignalFunc::ProcessFit(TObjArray * const arrhist) {
  //
  // Fit the +- invariant mass distribution only
  // Here we assume that the combined fit function is a sum of the signal and background functions
  //    and that the signal function is always the first term of this sum
  //
  
  fHistDataPM = (TH1F*)(arrhist->At(1))->Clone("histPM");  // +-    SE
  fHistDataPM->Sumw2();
  if(fRebin>1)
    fHistDataPM->Rebin(fRebin);
  
  fHistSignal = new TH1F("HistSignal", "Like-Sign substracted signal",
                         fHistDataPM->GetXaxis()->GetNbins(),
                         fHistDataPM->GetXaxis()->GetXmin(), fHistDataPM->GetXaxis()->GetXmax());
  fHistBackground = new TH1F("HistBackground", "Like-sign contribution",
                             fHistDataPM->GetXaxis()->GetNbins(),
                             fHistDataPM->GetXaxis()->GetXmin(), fHistDataPM->GetXaxis()->GetXmax());
  
  // the starting parameters of the fit function and their limits can be tuned
  // by the user in its macro
  fHistDataPM->Fit(fFuncSigBack, fFitOpt.Data(), "", fFitMin, fFitMax);
  TFitResultPtr pmFitPtr = fHistDataPM->Fit(fFuncSigBack, fFitOpt.Data(), "", fFitMin, fFitMax);
  //TFitResult *pmFitResult = pmFitPtr.Get(); // used only with TF1Helper
  fFuncSignal->SetParameters(fFuncSigBack->GetParameters());
  fFuncBackground->SetParameters(fFuncSigBack->GetParameters()+fFuncSignal->GetNpar());
  
  for(Int_t iBin=1; iBin<=fHistDataPM->GetXaxis()->GetNbins(); iBin++) {
    Double_t m = fHistDataPM->GetBinCenter(iBin);
    Double_t pm = fHistDataPM->GetBinContent(iBin);
    Double_t epm = fHistDataPM->GetBinError(iBin);
    Double_t bknd = fFuncBackground->Eval(m);
    Double_t ebknd = 0;
    for(Int_t iPar=fFuncSignal->GetNpar(); iPar<fFuncSigBack->GetNpar(); iPar++) {
/* TF1Helper problem on alien compilation
      for(Int_t jPar=iPar; jPar<fFuncSigBack->GetNpar(); jPar++) {
        TF1 gradientIpar("gradientIpar",
                         ROOT::TF1Helper::TGradientParFunction(iPar-fFuncSignal->GetNpar(),fFuncBackground),0,0,0);
        TF1 gradientJpar("gradientJpar",
                         ROOT::TF1Helper::TGradientParFunction(jPar-fFuncSignal->GetNpar(),fFuncBackground),0,0,0);
        ebknd += pmFitResult->CovMatrix(iPar,jPar)*
          gradientIpar.Eval(m)*gradientJpar.Eval(m)*
          (iPar==jPar ? 1.0 : 2.0);
      }
*/
    }
    Double_t signal = pm-bknd;
    Double_t error = TMath::Sqrt(epm*epm+ebknd);
    fHistSignal->SetBinContent(iBin, signal);
    fHistSignal->SetBinError(iBin, error);
    fHistBackground->SetBinContent(iBin, bknd);
    fHistBackground->SetBinError(iBin, TMath::Sqrt(ebknd));
  }
  
  if(fUseIntegral) {
    // signal
    fValues(0) = fFuncSignal->Integral(fIntMin, fIntMax)/fHistDataPM->GetBinWidth(1);
    fErrors(0) = 0;
    for(Int_t iPar=0; iPar<fFuncSignal->GetNpar(); iPar++) {
/* TF1Helper problem on alien compilation
      for(Int_t jPar=iPar; jPar<fFuncSignal->GetNpar(); jPar++) {
        TF1 gradientIpar("gradientIpar",
                         ROOT::TF1Helper::TGradientParFunction(iPar,fFuncSignal),0,0,0);
        TF1 gradientJpar("gradientJpar",
                         ROOT::TF1Helper::TGradientParFunction(jPar,fFuncSignal),0,0,0);
        fErrors(0) += pmFitResult->CovMatrix(iPar,jPar)*
          gradientIpar.Integral(fIntMin,fIntMax)*gradientJpar.Integral(fIntMin,fIntMax)*
          (iPar==jPar ? 1.0 : 2.0);
      }
*/
    }
    // background
    fValues(1) = fFuncBackground->Integral(fIntMin, fIntMax)/fHistDataPM->GetBinWidth(1);
    fErrors(1) = 0;
    for(Int_t iPar=fFuncSignal->GetNpar(); iPar<fFuncSigBack->GetNpar(); iPar++) {
/* TF1Helper problem on alien compilation
      for(Int_t jPar=iPar; jPar<fFuncSigBack->GetNpar(); jPar++) {
        TF1 gradientIpar("gradientIpar",
                         ROOT::TF1Helper::TGradientParFunction(iPar-fFuncSignal->GetNpar(),fFuncBackground),0,0,0);
        TF1 gradientJpar("gradientJpar",
                         ROOT::TF1Helper::TGradientParFunction(jPar-fFuncSignal->GetNpar(),fFuncBackground),0,0,0);
        fErrors(1) += pmFitResult->CovMatrix(iPar,jPar)*
          gradientIpar.Integral(fIntMin, fIntMax)*gradientJpar.Integral(fIntMin, fIntMax)*
          (iPar==jPar ? 1.0 : 2.0);
      }
*/
    }
  }
  else {
    // signal
    fValues(0) = fHistSignal->IntegralAndError(fHistSignal->FindBin(fIntMin),
                                               fHistSignal->FindBin(fIntMax), fErrors(0));
    // background
    fValues(1) = fHistBackground->IntegralAndError(fHistBackground->FindBin(fIntMin),
      fHistBackground->FindBin(fIntMax),
      fErrors(1));
  }
  // S/B and significance
  SetSignificanceAndSOB();
  fValues(4) = fFuncSigBack->GetParameter(fParMass);
  fErrors(4) = fFuncSigBack->GetParError(fParMass);
  fValues(5) = fFuncSigBack->GetParameter(fParMassWidth);
  fErrors(5) = fFuncSigBack->GetParError(fParMassWidth);
  
  fProcessed = kTRUE;
}

//______________________________________________
void AliDielectronSignalFunc::ProcessLS(TObjArray * const arrhist) {
  //
  // Substract background using the like-sign spectrum
  //
  fHistDataPP = (TH1F*)(arrhist->At(0))->Clone("histPP");  // ++    SE
  fHistDataPM = (TH1F*)(arrhist->At(1))->Clone("histPM");  // +-    SE
  fHistDataMM = (TH1F*)(arrhist->At(2))->Clone("histMM");  // --    SE
  if (fRebin>1) {
    fHistDataPP->Rebin(fRebin);
    fHistDataPM->Rebin(fRebin);
    fHistDataMM->Rebin(fRebin);
  }
  fHistDataPP->Sumw2();
  fHistDataPM->Sumw2();
  fHistDataMM->Sumw2();
  
  fHistSignal = new TH1F("HistSignal", "Like-Sign substracted signal",
                         fHistDataPM->GetXaxis()->GetNbins(),
                         fHistDataPM->GetXaxis()->GetXmin(), fHistDataPM->GetXaxis()->GetXmax());
  fHistBackground = new TH1F("HistBackground", "Like-sign contribution",
                             fHistDataPM->GetXaxis()->GetNbins(),
                             fHistDataPM->GetXaxis()->GetXmin(), fHistDataPM->GetXaxis()->GetXmax());
  
  // fit the +- mass distribution
  fHistDataPM->Fit(fFuncSigBack, fFitOpt.Data(), "", fFitMin, fFitMax);
  fHistDataPM->Fit(fFuncSigBack, fFitOpt.Data(), "", fFitMin, fFitMax);
  // declare the variables where the like-sign fit results will be stored
//   TFitResult *ppFitResult = 0x0;
//   TFitResult *mmFitResult = 0x0;
  // fit the like sign background
  TF1 *funcClonePP = (TF1*)fFuncBackground->Clone("funcClonePP");
  TF1 *funcCloneMM = (TF1*)fFuncBackground->Clone("funcCloneMM");
  fHistDataPP->Fit(funcClonePP, fFitOpt.Data(), "", fFitMin, fFitMax);
  fHistDataPP->Fit(funcClonePP, fFitOpt.Data(), "", fFitMin, fFitMax);
//   TFitResultPtr ppFitPtr = fHistDataPP->Fit(funcClonePP, fFitOpt.Data(), "", fFitMin, fFitMax);
//   ppFitResult = ppFitPtr.Get();
  fHistDataMM->Fit(funcCloneMM, fFitOpt.Data(), "", fFitMin, fFitMax);
  fHistDataMM->Fit(funcCloneMM, fFitOpt.Data(), "", fFitMin, fFitMax);
//   TFitResultPtr mmFitPtr = fHistDataMM->Fit(funcCloneMM, fFitOpt.Data(), "", fFitMin, fFitMax);
//   mmFitResult = mmFitPtr.Get();
  
  for(Int_t iBin=1; iBin<=fHistDataPM->GetXaxis()->GetNbins(); iBin++) {
    Double_t m = fHistDataPM->GetBinCenter(iBin);
    Double_t pm = fHistDataPM->GetBinContent(iBin);
    Double_t pp = funcClonePP->Eval(m);
    Double_t mm = funcCloneMM->Eval(m);
    Double_t epm = fHistDataPM->GetBinError(iBin);
    Double_t epp = 0;
    for(Int_t iPar=0; iPar<funcClonePP->GetNpar(); iPar++) {
/* TF1Helper problem on alien compilation
      for(Int_t jPar=iPar; jPar<funcClonePP->GetNpar(); jPar++) {
        TF1 gradientIpar("gradientIpar",
                         ROOT::TF1Helper::TGradientParFunction(iPar,funcClonePP),0,0,0);
        TF1 gradientJpar("gradientJpar",
                         ROOT::TF1Helper::TGradientParFunction(jPar,funcClonePP),0,0,0);
        epp += ppFitResult->CovMatrix(iPar,jPar)*
          gradientIpar.Eval(m)*gradientJpar.Eval(m)*
          (iPar==jPar ? 1.0 : 2.0);
      }
*/
    }
    Double_t emm = 0;
    for(Int_t iPar=0; iPar<funcCloneMM->GetNpar(); iPar++) {
/* TF1Helper problem on alien compilation
      for(Int_t jPar=iPar; jPar<funcCloneMM->GetNpar(); jPar++) {
        TF1 gradientIpar("gradientIpar",
                         ROOT::TF1Helper::TGradientParFunction(iPar,funcCloneMM),0,0,0);
        TF1 gradientJpar("gradientJpar",
                         ROOT::TF1Helper::TGradientParFunction(jPar,funcCloneMM),0,0,0);
        emm += mmFitResult->CovMatrix(iPar,jPar)*
          gradientIpar.Eval(m)*gradientJpar.Eval(m)*
          (iPar==jPar ? 1.0 : 2.0);
      }
*/
    }
    
    Double_t signal = pm-2.0*TMath::Sqrt(pp*mm);
    Double_t background = 2.0*TMath::Sqrt(pp*mm);
    // error propagation on the signal calculation above
    Double_t esignal = TMath::Sqrt(epm*epm+(mm/pp)*epp+(pp/mm)*emm);
    Double_t ebackground = TMath::Sqrt((mm/pp)*epp+(pp/mm)*emm);
    fHistSignal->SetBinContent(iBin, signal);
    fHistSignal->SetBinError(iBin, esignal);
    fHistBackground->SetBinContent(iBin, background);
    fHistBackground->SetBinError(iBin, ebackground);
  }
  
  // signal
  fValues(0) = fHistSignal->IntegralAndError(fHistSignal->FindBin(fIntMin),
                                             fHistSignal->FindBin(fIntMax), fErrors(0));
  // background
  fValues(1) = fHistBackground->IntegralAndError(fHistBackground->FindBin(fIntMin),
                                                 fHistBackground->FindBin(fIntMax),
                                                 fErrors(1));
  // S/B and significance
  SetSignificanceAndSOB();
  fValues(4) = fFuncSigBack->GetParameter(fParMass);
  fErrors(4) = fFuncSigBack->GetParError(fParMass);
  fValues(5) = fFuncSigBack->GetParameter(fParMassWidth);
  fErrors(5) = fFuncSigBack->GetParError(fParMassWidth);
  
  fProcessed = kTRUE;
}

//______________________________________________
void AliDielectronSignalFunc::ProcessEM(TObjArray * const arrhist) {
  //
  // Substract background with the event mixing technique
  //
  arrhist->GetEntries();   // just to avoid the unused parameter warning
  AliError("Event mixing for background substraction method not implemented!");
}

//______________________________________________
void AliDielectronSignalFunc::SetFunctions(TF1 * const combined, TF1 * const sig, TF1 * const back,
                                           Int_t parM, Int_t parMres)
{
  //
  // Set the signal, background functions and combined fit function
  // Note: The process method assumes that the first n parameters in the
  //       combined fit function correspond to the n parameters of the signal function
  //       and the n+1 to n+m parameters to the m parameters of the background function!!!
  
  if (!sig||!back||!combined) {
    AliError("Both, signal and background function need to be set!");
    return;
  }
  fFuncSignal=sig;
  fFuncBackground=back;
  fFuncSigBack=combined;
  fParMass=parM;
  fParMassWidth=parMres;
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
    fFuncSignal=new TF1("DieleSignal","gaus",2.5,4);
    fFuncBackground=new TF1("DieleBackground","pol1",2.5,4);
    fFuncSigBack=new TF1("DieleCombined","gaus+pol1(3)",2.5,4);
    
    fFuncSigBack->SetParameters(1,3.1,.05,2.5,1);
    fFuncSigBack->SetParLimits(0,0,10000000);
    fFuncSigBack->SetParLimits(1,3.05,3.15);
    fFuncSigBack->SetParLimits(2,.02,.1);
  }
  else if (type==1){
    fFuncSignal=new TF1("DieleSignal","gaus",2.5,4);
    fFuncBackground=new TF1("DieleBackground","[0]*exp(-(x-[1])/[2])",2.5,4);
    fFuncSigBack=new TF1("DieleCombined","gaus+[3]*exp(-(x-[4])/[5])",2.5,4);
    
    fFuncSigBack->SetParameters(1,3.1,.05,1,2.5,1);
    fFuncSigBack->SetParLimits(0,0,10000000);
    fFuncSigBack->SetParLimits(1,3.05,3.15);
    fFuncSigBack->SetParLimits(2,.02,.1);
  }
  else if (type==2){
    // half gaussian, half exponential signal function
    // exponential background
    fFuncSignal = new TF1("DieleSignal","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",2.5,4);
    fFuncBackground = new TF1("DieleBackground","[0]*exp(-(x-[1])/[2])+[3]",2.5,4);
    fFuncSigBack = new TF1("DieleCombined","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))+[4]*exp(-(x-[5])/[6])+[7]",2.5,4);
    fFuncSigBack->SetParameters(1.,3.1,.05,.1,1,2.5,1,0);
    
    fFuncSigBack->SetParLimits(0,0,10000000);
    fFuncSigBack->SetParLimits(1,3.05,3.15);
    fFuncSigBack->SetParLimits(2,.02,.1);
    fFuncSigBack->FixParameter(6,2.5);
    fFuncSigBack->FixParameter(7,0);
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
  
  fFuncSigBack->SetNpx(200);
  fFuncSigBack->SetRange(fIntMin,fIntMax);
  fFuncBackground->SetNpx(200);
  fFuncBackground->SetRange(fIntMin,fIntMax);
  
  TGraph *grSig=new TGraph(fFuncSigBack);
  grSig->SetFillColor(kGreen);
  grSig->SetFillStyle(3001);
  
  TGraph *grBack=new TGraph(fFuncBackground);
  grBack->SetFillColor(kRed);
  grBack->SetFillStyle(3001);
  
  grSig->SetPoint(0,grBack->GetX()[0],grBack->GetY()[0]);
  grSig->SetPoint(grSig->GetN()-1,grBack->GetX()[grBack->GetN()-1],grBack->GetY()[grBack->GetN()-1]);
  
  grBack->SetPoint(0,grBack->GetX()[0],0.);
  grBack->SetPoint(grBack->GetN()-1,grBack->GetX()[grBack->GetN()-1],0.);
  
  fFuncSigBack->SetRange(fFitMin,fFitMax);
  fFuncBackground->SetRange(fFitMin,fFitMax);
  
  if (!drawOpt.Contains("same")){
    if (fHistDataPM){
      fHistDataPM->Draw();
      grSig->Draw("f");
    } else {
      grSig->Draw("af");
    }
  } else {
    grSig->Draw("f");
  }
  if(fMethod==kFitted) grBack->Draw("f");
  fFuncSigBack->Draw("same");
  fFuncSigBack->SetLineWidth(2);
  if(fMethod==kLikeSign) {
    fHistDataPP->SetLineWidth(2);
    fHistDataPP->SetLineColor(6);
    fHistDataPP->Draw("same");
    fHistDataMM->SetLineWidth(2);
    fHistDataMM->SetLineColor(8);
    fHistDataMM->Draw("same");
  }
  
  if(fMethod==kFitted)
    fFuncBackground->Draw("same");
  
  if (optStat) DrawStats();
  
}
