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
#include <RooFit.h>
#include <RooAddPdf.h>

#include "AliDielectronSignalBase.h"

class AliDielectronSignalFunc : public AliDielectronSignalBase {
public:
  AliDielectronSignalFunc();
  AliDielectronSignalFunc(const char*name, const char* title);

  virtual ~AliDielectronSignalFunc();

  void Process(TH1 * const hist);
  virtual void Process(TObjArray * const arrhist);
  
  void SetFunctions(TF1 * const sig, TF1 * const back, TF1 * const combined, Int_t parM=-1, Int_t parMres=-1);
  void InitParams();
  void SetFitOption(const char* opt) {fFitOpt=opt;}
  void SetDefaults(Int_t type);
  void SetFitRange(Double_t min, Double_t max) {fFitMin=min;fFitMax=max;}

  void SetUseIntegral(Bool_t integral=kTRUE) { fUseIntegral=integral; }
  
  const char* GetFitOption()          const { return fFitOpt.Data(); }
  TF1*  GetSignalFunction()     const { return fSignal;        }
  TF1*  GetBackgroundFunction() const { return fBackground;    }
  TF1*  GetCombinedFunction()   const { return fSigBack;       }
  
  Double_t GetFitMin()      const { return fFitMin; }
  Double_t GetFitMax()      const { return fFitMax; }
  
  virtual void Draw(const Option_t* option = "");
  
private:
  TF1 *fSignal;                // Function for the signal description
  TF1 *fBackground;            // Function for the background description
  TF1 *fSigBack;               // Combined function signal plus background
  RooAddPdf fFitFunc;          // RooFit fit function
  TVectorD fVInitParam;        // initial fit parameters
  TString fFitOpt;             // fit option used

  Bool_t fUseIntegral;         // use integral instead of bin counts
  
  Double_t fFitMin;            // fit range min
  Double_t fFitMax;            // fit range max

  Int_t fParM;                 // Paramter which defines the mass
  Int_t fParMres;              // Paramter which defines the resolution of the mass
  
  AliDielectronSignalFunc(const AliDielectronSignalFunc &c);
  AliDielectronSignalFunc &operator=(const AliDielectronSignalFunc &c);

  ClassDef(AliDielectronSignalFunc,1)         // Dielectron SignalFunc
};


#endif
