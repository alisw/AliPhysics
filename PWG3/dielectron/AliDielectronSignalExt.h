#ifndef ALIDIELECTRONSIGNALEXT_H
#define ALIDIELECTRONSIGNALEXT_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#############################################################
//#                                                           # 
//#           Class AliDielectronSignalExt                    #
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
#include <TH1.h>

#include "AliDielectronSignalBase.h"

class AliDielectronSignalExt : public AliDielectronSignalBase {

public:
  AliDielectronSignalExt();
  AliDielectronSignalExt(const char*name, const char* title);

  virtual ~AliDielectronSignalExt();

  virtual void Process(TObjArray* const arrhist);
  void ProcessLS(TObjArray* const arrhist);  // like-sign method
  void ProcessEM(TObjArray* const arrhist);  // event mixing method
  
  void SetMethod(Int_t method){ fMethod=method; }
  void SetHistograms(TH1F* const unlike, TH1F* const backg, TH1F* const signal);
  void SetRebin(Int_t rebin) {fRebin = rebin;}
  void SetDrawRange(Double_t min, Double_t max){fDrawMin=min; fDrawMax=max;}
  void SetFitRange(Double_t min, Double_t max) {fFitMin=min;fFitMax=max;}

  Double_t GetFitMin()      const { return fFitMin; }
  Double_t GetFitMax()      const { return fFitMax; }
  void Rebin(Int_t rebin);
  
  virtual void Draw(const Option_t* option = "");


private:
  
  TH1F *fSignPM;              // histogram of unlike sign (plus-minus)
  TH1F *fSignPP;              // histogram of like sign (plus-plus)
  TH1F *fSignMM;              // histogram of like sign (minus-minus)
  TH1F *fBackground;          // histogram of like-sign background
  TH1F *fSignal;              // histogram of subtracted signal
  
  Int_t    fMethod;           // subtraction method. 1(like-sign), 2(event mixing)
  Int_t    fRebin;            // number of histogram rebin iteration
  Int_t    fBins;             // number of bins in X axis
  Double_t fDrawMin;          // minimum X when drawing 
  Double_t fDrawMax;          // maximum X when drawing
  Double_t fFitMin;           // fit range min
  Double_t fFitMax;           // fit range max

  AliDielectronSignalExt(const AliDielectronSignalExt &c);
  AliDielectronSignalExt &operator=(const AliDielectronSignalExt &c);

  ClassDef(AliDielectronSignalExt,1)         // Dielectron SignalFunc
};



#endif
