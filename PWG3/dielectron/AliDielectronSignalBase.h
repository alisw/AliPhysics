#ifndef ALIDIELECTRONSIGNALBASE_H
#define ALIDIELECTRONSIGNALBASE_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

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
class TH1F;

class AliDielectronSignalBase : public TNamed {
public:
  enum EBackgroundMethod {
    kFitted = 0,
    kLikeSign,
    kLikeSignArithm,
    kEventMixing,
    kRotation
  };

  AliDielectronSignalBase();
  AliDielectronSignalBase(const char*name, const char* title);
  
  virtual ~AliDielectronSignalBase();
  
  void SetIntegralRange(Double_t min, Double_t max) {fIntMin=min;fIntMax=max;}
  void SetFitRange(Double_t min, Double_t max) {fFitMin=min; fFitMax=max;}
  void SetRebin(Int_t factor) {fRebin = factor;}
  void SetMethod(EBackgroundMethod method) {fMethod = method;}

  const TVectorD& GetValues() const {return fValues;}
  const TVectorD& GetErrors() const {return fErrors;}

  Double_t GetIntegralMin()          const { return fIntMin; }
  Double_t GetIntegralMax()          const { return fIntMax; }
  Double_t GetSignal()               const { return fValues(0);}
  Double_t GetSignalError()          const { return fErrors(0);}
  Double_t GetBackground()           const { return fValues(1);}
  Double_t GetBackgroundError()      const { return fErrors(1);}
  Double_t GetSignificance()         const { return fValues(2);}
  Double_t GetSignificanceError()    const { return fErrors(2);}
  Double_t GetSB()                   const { return fValues(3);}
  Double_t GetSBError()              const { return fErrors(3);}
  Double_t GetMass()                 const { return fValues(4);}
  Double_t GetMassError()            const { return fErrors(4);}
  Double_t GetMassWidth()            const { return fValues(5);}
  Double_t GetMassWidthError()       const { return fErrors(5);}

  TH1* GetSignalHistogram()      const {return fHistSignal;}
  TH1* GetBackgroundHistogram()  const {return fHistBackground;}
  TH1* GetUnlikeSignHistogram()  const {return fHistDataPM;}

  void SetScaleRawToBackground(Double_t intMin, Double_t intMax) { fScaleMin=intMin; fScaleMax=intMax; }
  Double_t GetScaleFactor() const { return fScaleFactor; }
  
  static Double_t ScaleHistograms(TH1* histRaw, TH1* histBackground, Double_t intMin, Double_t intMax);
  
  virtual void Print(Option_t *option="") const;

  /**
  This function needs to be implemented by the signal extraction classes.
  Here all the work should be done.
  
  The signal extraction is done on the mass spectra.
  The TObjArray should contain the Inv. Mass spectra of the 10 possible combinations
     for single and mixed events defined in AliDielectron.cxx
  */
  virtual void Process(TObjArray * const /*arrhist*/) = 0;
    
protected: 

  TH1 *fHistSignal;                  // histogram of pure signal
  TH1 *fHistBackground;              // histogram of background (fitted=0, like-sign=1, event mixing=2)
  TH1 *fHistDataPM;                  // histogram of selected +- pair candidates
  TH1 *fHistDataPP;                  // histogram of selected ++ pair candidates
  TH1 *fHistDataMM;                  // histogram of selected -- pair candidates

  TVectorD fValues;                   // values
  TVectorD fErrors;                   // value errors

  Double_t fIntMin;                   // signal extraction range min
  Double_t fIntMax;                   // signal extraction range max
  Double_t fFitMin;                   // fit range lowest inv. mass
  Double_t fFitMax;                   // fit range highest inv. mass

  Int_t fRebin;                       // histogram rebin factor
  EBackgroundMethod fMethod;          // method for background substraction
  Double_t fScaleMin;                 // min for scaling of raw and background histogram
  Double_t fScaleMax;                 // max for scaling of raw and background histogram
  Double_t fScaleFactor;              // scale factor of raw to background histogram scaling
  
  Bool_t fProcessed;                  // flag
  
  void SetSignificanceAndSOB();       // calculate the significance and S/B
  TPaveText* DrawStats(Double_t x1=0., Double_t y1=0., Double_t x2=0., Double_t y2=0.);

  AliDielectronSignalBase(const AliDielectronSignalBase &c);
  AliDielectronSignalBase &operator=(const AliDielectronSignalBase &c);

  ClassDef(AliDielectronSignalBase,4) // Dielectron SignalBase
};

inline void AliDielectronSignalBase::SetSignificanceAndSOB()
{
  //
  // Calculate S/B and significance
  //
  // Signal/Background
  fValues(3) = (fValues(1)>0 ? fValues(0)/fValues(1) : 0);
  Float_t epsSig = (fValues(0)>0 ? fErrors(0)/fValues(0) : 0);
  Float_t epsBknd = (fValues(1)>0 ? fErrors(1)/fValues(1) : 0);
  fErrors(3) = fValues(3)*TMath::Sqrt(epsSig*epsSig + epsBknd*epsBknd);
  // Significance
  fValues(2) = ((fValues(0)+fValues(1))>0 ? fValues(0)/TMath::Sqrt(fValues(0)+fValues(1)) : 0);
  Float_t s = fValues(0); Float_t b = fValues(1);
  fErrors(2) = ((s+b)>0 ? TMath::Sqrt((s*(s+2*b)*(s+2*b)+b*s*s)/(4*TMath::Power(s+b,3))) : 0);
}

#endif
