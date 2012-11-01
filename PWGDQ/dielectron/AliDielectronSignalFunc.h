#ifndef ALIDIELECTRONSIGNALFUNC_H
#define ALIDIELECTRONSIGNALFUNC_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#############################################################
//#                                                           # 
//#         Class AliDielectronSignalFunc                     #
//#                                                           #
//#  Authors:                                                 #
//#   Anton     Andronic, GSI / A.Andronic@gsi.de             #
//#   Ionut C.  Arsene,   GSI / I.C.Arsene@gsi.de             #
//#   Julian    Book,     Uni Ffm / Julian.Book@cern.ch       #
//#   Frederick Kramer,   Uni Ffm, / Frederick.Kramer@cern.ch #
//#   Magnus    Mager,    CERN / Magnus.Mager@cern.ch         #
//#   WooJin J. Park,     GSI / W.J.Park@gsi.de               #
//#   Jens      Wiechula, Uni HD / Jens.Wiechula@cern.ch      #
//#                                                           #
//#############################################################

/*
  Class used for extracting the signal from an invariant mass spectrum.
  It implements the AliDielectronSignalBase class and it uses user provided
  functions to fit the unlike-sign spectrum (and the like-sign one).
  
  Example usage:

  AliDielectronSignalFunc *signalProcess = new AliDielectronSignalFunc();
  TObjArray *histoArray = new TObjArray();
  histoArray->Add(signalPP);            // add the spectrum histograms to the array
  histoArray->Add(signalPM);            // the order is important !!!
  histoArray->Add(signalMM);
  // set the extraction method 
  // AliDielectronSignalBase::kFitted       -->  fit only the unlike-sign spectrum and extract everything from that
  // AliDielectronSignalBase::kLikeSign     -->  fit both the unlike- and like-sign spectra
  // AliDielectronSignalBase::kEventMixing  -->  fit both the unlike- and like-sign spectra from event mixing
  signalProcess->SetMethod(AliDielectronSignalBase::kLikeSign); 
  // Initialize the functions to be used and pass them to the signal object
  // External preparation of the functions can(should) be done as this can be a 5 or more parameter fit
  TF1* gaus = new TF1("gaus", "gaus", 0., 4.);
  TF1* expo = new TF1("expo", "[0]*exp([1]*x)", 0., 4.);
  TF1* combined = new TF1("combined", "gaus + [3]*exp([4]*x)", 0.,4.);
  combined->SetParameter(1, 3.1);
  combined->SetParameter(1, 0.1);
  signalPP->Fit(expo, "SME", "", 2.4, 4.0);
  signalPM->Fit(gaus, "SME", "", 3.0, 3.15);
  Double_t pars[5];
  gaus->GetParameters(&pars[0]);
  expo->GetParameters(&pars[3]);
  combined->SetParameters(pars);
  combined->SetParLimits(1, 3.05, 3.15);
  combined->SetParLimits(2, 0.03, 0.1);
  signalProcess->SetFunctions(combined, gaus, expo, 1, 2);

  signalProcess->SetFitRange(2.4,4.0);
  // Use the integral of the fit function to estimate the signal or not
  // The background will always be estimated from the fit
  //  signalProcess->SetUseIntegral(kTRUE);  
  signalProcess->SetFitOption("SME");
  // Give the range where the signal is calculated
  signalProcess->SetIntegralRange(3.0,3.15);
  signalProcess->SetRebin(2);
  signalProcess->Process(histoArray);
  signalProcess->Draw("stat");
  signalProcess->Print();
*/

#include <TVectorT.h>
#include <TString.h>
#include <TH1F.h>

#include "AliDielectronSignalBase.h"

class AliDielectronSignalFunc : public AliDielectronSignalBase {
public:
  AliDielectronSignalFunc();
  AliDielectronSignalFunc(const char*name, const char* title);
  AliDielectronSignalFunc(const AliDielectronSignalFunc &c);
  AliDielectronSignalFunc &operator=(const AliDielectronSignalFunc &c);

  virtual ~AliDielectronSignalFunc();

  virtual void Process(TObjArray * const arrhist);
  void ProcessFit(TObjArray * const arrhist);      // fit the SE +- distribution
  void ProcessLS(TObjArray * const arrhist);       // substract the fitted SE like-sign background
  void ProcessEM(TObjArray * const arrhist);       // substract the fitted SE+ME like-sign background

  void SetMCSignalShape(TH1F* hist) { fgHistSimPM=hist; }
  void SetUseIntegral(Bool_t flag=kTRUE) {fUseIntegral = flag;};
  void SetFunctions(TF1 * const combined, TF1 * const sig=0, TF1 * const back=0, Int_t parM=1, Int_t parMres=2);
  void CombineFunc(TF1 * const peak=0, TF1 * const bgnd=0);
  void SetFitOption(const char* opt) {
    fFitOpt=opt; 
    fFitOpt.ToLower(); 
    if(!fFitOpt.Contains("s")) fFitOpt += "s";
  }
  void SetDefaults(Int_t type);
    
  TF1*  GetSignalFunction()     const { return fFuncSignal;        }
  TF1*  GetBackgroundFunction() const { return fFuncBackground;    }
  TF1*  GetCombinedFunction()   const { return fFuncSigBack;       }
  
  Int_t GetPolDeg()             const { return fPolDeg;            }
  Int_t GetDof()                const { return fDof;            }
  Double_t GetChi2Dof()         const { return fChi2Dof;           }
  
  virtual void Draw(const Option_t* option = "");
  
private:

  // peak functions
  static Double_t PeakFunMC(const Double_t *x, const Double_t *par); // peak function from a mc histo
  static Double_t PeakFunCB(const Double_t *x, const Double_t *par); // crystal ball function
  static Double_t PeakFunGaus(const Double_t *x, const Double_t *par); // gaussian

  static Double_t PeakBgndFun(const Double_t *x, const Double_t *par); // combine any bgrd and any peak function

  static TF1 *fFuncSignal;                // Function for the signal description
  static TF1 *fFuncBackground;            // Function for the background description
  TF1 *fFuncSigBack;               // Combined function signal plus background
  Int_t fParMass;                  // the index of the parameter corresponding to the resonance mass
  Int_t fParMassWidth;             // the index of the parameter corresponding to the resonance mass width
  
  TString fFitOpt;             // fit option used
  Bool_t fUseIntegral;         // use the integral of the fitted functions to extract signal and background

  Int_t    fPolDeg;                // polynomial degree of the background function
  Int_t    fDof;                   // degrees of freedom
  Double_t fChi2Dof;               // chi2/dof of the fitted inv mass spectra
  
  static Int_t    fNparPeak;              // number of parameters for peak function
  static Int_t    fNparBgnd;              // number of parameters for background function

  static TH1F* fgHistSimPM;         // simulated peak shape 
  
  ClassDef(AliDielectronSignalFunc,2)         // Dielectron SignalFunc
};

#endif
