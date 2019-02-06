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

#include <TVectorT.h>
#include <TString.h>
#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>

#include "AliDielectronSignalExt.h"

class AliDielectronSignalFunc : public AliDielectronSignalExt {
public:
  AliDielectronSignalFunc();
  AliDielectronSignalFunc(const char*name, const char* title);
  AliDielectronSignalFunc(const AliDielectronSignalFunc &c);
  AliDielectronSignalFunc &operator=(const AliDielectronSignalFunc &c);

  virtual ~AliDielectronSignalFunc();

  virtual void Process(TObjArray * const arrhist);
  void ProcessCombinatorialPlusFit(TObjArray * const arrhist);      // subtract the combinatorial bg, then fit the SE +- distribution
  void ProcessFit(TObjArray * const arrhist);      // fit the SE +- distribution
  void ProcessFitLS(TObjArray * const arrhist);       // substract the fitted SE like-sign background
  void ProcessFitEM(TObjArray * const arrhist);       // substract the fitted SE+ME like-sign background

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
  
  //  Int_t GetPolDeg()             const { return fPolDeg;            }
  Int_t GetDof()                const { return fDof;            }
  Double_t GetChi2Dof()         const { return fChi2Dof;           }
  
  virtual void Draw(const Option_t* option = "");

  // peak functions
  Double_t PeakFunMC(const Double_t *x, const Double_t *par); // peak function from a mc histo
  Double_t PeakFunCB(const Double_t *x, const Double_t *par); // crystal ball function
  Double_t PeakFunGaus(const Double_t *x, const Double_t *par); // gaussian
  
  
  TH1* GetCombinatorialBackgroundHistogram()      const {return fHistCombinatorialBackground;}

protected:

  static Double_t PeakBgndFun(const Double_t *x, const Double_t *par); // combine any bgrd and any peak function

  static TF1 *fFuncSignal;                // Function for the signal description
  static TF1 *fFuncBackground;            // Function for the background description
  TF1 *fFuncSigBack;               // Combined function signal plus background
  Int_t fParMass;                  // the index of the parameter corresponding to the resonance mass
  Int_t fParMassWidth;             // the index of the parameter corresponding to the resonance mass width
  
  TString fFitOpt;             // fit option used
  Bool_t fUseIntegral;         // use the integral of the fitted functions to extract signal and background

  //  Int_t    fPolDeg;                // polynomial degree of the background function
  Int_t    fDof;                   // degrees of freedom
  Double_t fChi2Dof;               // chi2/dof of the fitted inv mass spectra
  
  static Int_t    fNparPeak;              // number of parameters for peak function
  static Int_t    fNparBgnd;              // number of parameters for background function
  TH1 *fHistCombinatorialBackground;     // histogram of combinatorial bg

  
  ClassDef(AliDielectronSignalFunc, 4)         // class for signal extraction using a combined bgrd+signal fit
};

#endif
