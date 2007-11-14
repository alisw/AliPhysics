#ifndef ALITPCFITPAD_H
#define ALITPCFITPAD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <iostream>
#include "AliTPCCalPadRegion.h"
#include <TLinearFitter.h>
#include <TIterator.h>

using namespace std;

class TString;

class AliTPCFitPad: public AliTPCCalPadRegion {
public:
   AliTPCFitPad() : AliTPCCalPadRegion() { }
   AliTPCFitPad(const AliTPCFitPad& obj) : AliTPCCalPadRegion(obj), fNdim(obj.fNdim), fFormula(obj.fFormula), fOpt(obj.fOpt) { }
   AliTPCFitPad(Int_t ndim, const char* formula, Option_t* opt = "D");
   AliTPCFitPad& operator=(const AliTPCFitPad& rhs);
   //AliTPCFitPad(const char *name, const char *title) : AliTPCCalPadRegion(name, title) { }
   //AliTPCFitPad(const TString &name, const TString &title) : AliTPCCalPadRegion(name, title) { }
   virtual ~AliTPCFitPad();

   void           Add(AliTPCFitPad* fit);
   TLinearFitter* GetFitter(UInt_t segment, UInt_t padType, Bool_t workaround = kFALSE);
   TLinearFitter* GetFitterSimple(UInt_t segment, UInt_t padType);
   Int_t          Evaluate(Bool_t robust = kFALSE, Double_t frac = -1.);

protected:
   Int_t   fNdim;         // used for generating new TLinearFitter objects
   TString fFormula;      // used for generating new TLinearFitter objects
   TString fOpt;          // used for generating new TLinearFitter objects
   
   ClassDef(AliTPCFitPad, 1)
};


#endif
