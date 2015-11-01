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
  AliDielectronSignalExt(const char*name, const char* title, bool enummaps);

  virtual ~AliDielectronSignalExt();

  virtual void Process(TObjArray* const arrhist);
  void ProcessLS(TObjArray* const arrhist);  // like-sign method
  void ProcessEM(TObjArray* const arrhist);  // event mixing method
  void ProcessRotation(TObjArray* const arrhist);  // event mixing method

  virtual void Draw(const Option_t* option = "");

private:

  AliDielectronSignalExt(const AliDielectronSignalExt &c);
  AliDielectronSignalExt &operator=(const AliDielectronSignalExt &c);

  ClassDef(AliDielectronSignalExt,2)    // class for signal extraction using LS, ME or ROT
};

#endif
