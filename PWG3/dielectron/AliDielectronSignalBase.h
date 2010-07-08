#ifndef ALIDIELECTRONSIGNALBASE_H
#define ALIDIELECTRONSIGNALBASE_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#############################################################
//#                                                           # 
//#         Class AliDielectronSignalBase                       #
//#         Manage Cuts on the legs of the pair               #
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

#include <TNamed.h>
#include <TVectorT.h>
#include <TMath.h>

class TObjArray;
class TPaveText;
class TH1;

class AliDielectronSignalBase : public TNamed {
public:
  AliDielectronSignalBase();
  AliDielectronSignalBase(const char*name, const char* title);

  virtual ~AliDielectronSignalBase();

  
  void SetIntegralRange(Double_t min, Double_t max) {fIntMin=min;fIntMax=max;}
  
  void SetSignal(Double_t val, Double_t valErr)               {fValues(0)=val; fErrors(0)=valErr;}
  void SetBackground(Double_t val, Double_t valErr)           {fValues(1)=val; fErrors(1)=valErr;}
  void SetSignificance(Double_t val, Double_t valErr)         {fValues(2)=val; fErrors(2)=valErr;}
  void SetSignalOverBackground(Double_t val, Double_t valErr) {fValues(3)=val; fErrors(3)=valErr;}
  void SetMass(Double_t val, Double_t valErr)                 {fValues(4)=val; fErrors(4)=valErr;}
  void SetMassWidth(Double_t val, Double_t valErr)            {fValues(5)=val; fErrors(5)=valErr;}
  
  const TVectorD& GetValues() const {return fValues;}
  const TVectorD& GetErrors() const {return fErrors;}

  Double_t GetIntegralMin() const { return fIntMin; }
  Double_t GetIntegralMax() const { return fIntMax; }
  
  Double_t GetSignal()               {return fValues(0);}
  Double_t GetBackground()           {return fValues(1);}
  Double_t GetSignificance()         {return fValues(2);}
  Double_t GetSignalOverBackground() {return fValues(3);}
  Double_t GetMass()                 {return fValues(4);}
  Double_t GetMassWidth()            {return fValues(5);}
  
  Double_t GetSignalError()               {return fErrors(0);}
  Double_t GetBackgroundError()           {return fErrors(1);}
  Double_t GetSignificanceError()         {return fErrors(2);}
  Double_t GetSignalOverBackgroundError() {return fErrors(3);}
  Double_t GetMassError()                 {return fErrors(4);}
  Double_t GetMassWidthError()            {return fValues(5);}
  
  void GetSignal(Double_t &val, Double_t &valErr)               {val=fValues(0); valErr=fErrors(0);}
  void GetBackground(Double_t &val, Double_t &valErr)           {val=fValues(1); valErr=fErrors(1);}
  void GetSignificance(Double_t &val, Double_t &valErr)         {val=fValues(2); valErr=fErrors(2);}
  void GetSignalOverBackground(Double_t &val, Double_t &valErr) {val=fValues(3); valErr=fErrors(3);}

  /**
  This function needs to be implemented by the signal extraction classes.
  
  The signal extraction is done on the mass spectra.
  The TObjArray should contain the Inv. Mass spectra of the 10 possible combinations
     for single and mixed events defined in AliDielectron.cxx
  In the referece TVectorDs the values and value errors of the result should be stored:
  Parameter - 0: Signal entries
              1: Background entries
              2: Significance ( S/sqr(S+B) )
              3: Signal/Background
              4: Mass
              5: Mass width

  It is enough to calculate the signal and background and then call
            SetSignificanceAndSOB(TVectorD &values, TVectorD &errors)
            Please also read the description there!!!
  */
  virtual void Process(TObjArray * const /*arrhist*/) = 0;

protected:
  
  void SetSignificanceAndSOB();
  TPaveText* DrawStats(Double_t x1=0., Double_t y1=0., Double_t x2=0., Double_t y2=0.);
  void Reset() {fValues.Zero(); fErrors.Zero();}
  
private:
  TVectorD fValues;            // values
  TVectorD fErrors;            // value errors
  
  Double_t fIntMin;            // signal extraction range min
  Double_t fIntMax;            // signal extraction range max
  
  AliDielectronSignalBase(const AliDielectronSignalBase &c);
  AliDielectronSignalBase &operator=(const AliDielectronSignalBase &c);

  ClassDef(AliDielectronSignalBase,2)         // Dielectron SignalBase
};

//
// Inline functions
//

inline void AliDielectronSignalBase::SetSignificanceAndSOB()
{
  //
  // calculate significance and signal over background from values
  // it is expected that:
  // - 'values' and 'errors' have the size 4
  // - Parameters 0 are given: signal and signal error
  // - Parameters 1 are given: background and background error
  // - Optionally parameter 2 can be given: signal+background and its error
  //
  // The calculated significance and S/B and the corresponding errors
  //   are written in parameters 2 and 3, the first two parameters (signal and background)
  //   stay untouched

  Double_t signal=fValues(0),      signal_err=fErrors(0);
  Double_t background=fValues(1),  background_err=fErrors(1);
  Double_t sigPlusBack=fValues(2), sigPlusBack_err=fErrors(2);
  Double_t significance=0,        significance_err=0;
  Double_t sob=0,                 sob_err=0;
  
  if (sigPlusBack<1e-20){
    //calculate signal plus background
    sigPlusBack=signal+background;
    sigPlusBack_err=TMath::Sqrt( signal_err*signal_err + background_err*background_err );
  }

  if (sigPlusBack>1e-30){

    significance=signal/TMath::Sqrt(sigPlusBack);
    
    significance_err=TMath::Sqrt((signal_err/signal)*(signal_err/signal)+
                                 (0.5*sigPlusBack_err/sigPlusBack)*(0.5*sigPlusBack_err/sigPlusBack)
                                )*significance;
    
  }
  
  if (background>1e-30){
    sob=signal/background;
    sob_err=TMath::Sqrt((signal_err/signal)*(signal_err/signal)+
                        (background_err/background)*(background_err/background))*sob;
  }

  fValues(2)=significance;
  fErrors(2)=significance_err;
  
  fValues(3)=sob;
  fErrors(3)=sob_err;
}

#endif
