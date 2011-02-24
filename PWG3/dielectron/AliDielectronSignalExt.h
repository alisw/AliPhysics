#ifndef ALIDIELECTRONSIGNALEXT_H
#define ALIDIELECTRONSIGNALEXT_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

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

/*
  Class used for extracting the signal from an invariant mass spectrum.
  It implements the AliDielectronSignalBase class and uses the like-sign
  substraction method for estimating the signal and background.
  There is no fitting in this class, only bin counting.

  Example usage:
   AliDielectronSignalExt *signalProcess = new AliDielectronSignalExt();
   TObjArray *histoArray = new TObjArray();
   histoArray->Add(signalPP);                  // the order of putting the histograms in the array is important!!
   histoArray->Add(signalPM);
   histoArray->Add(signalMM);
   signalProcess->SetMethod(AliDielectronSignalBase::kLikeSign);  // or kEventMixing
   signalProcess->SetIntegralRange(3.0,3.15);   // J/Psi peak
   signalProcess->SetRebin(2);                  // rebin the histograms
   signalProcess->Process(histoArray);
   signalProcess->Draw("stat");
   signalProcess->Print();

*/

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
  void ProcessRotation(TObjArray* const arrhist);  // event mixing method

  virtual void Draw(const Option_t* option = "");

private:

  AliDielectronSignalExt(const AliDielectronSignalExt &c);
  AliDielectronSignalExt &operator=(const AliDielectronSignalExt &c);

  ClassDef(AliDielectronSignalExt,2)    // Dielectron SignalFunc
};

#endif
