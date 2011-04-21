#ifndef ALIRSNVALUE_H
#define ALIRSNVALUE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
////////////////////////////////////////////////////////////////////////////////
//
//  Collection of all values which can be computed within the package
//
////////////////////////////////////////////////////////////////////////////////

#include "TArrayD.h"

#include "AliRsnTarget.h"

class AliRsnValue : public AliRsnTarget {
public:

   AliRsnValue(const char *name = "", AliRsnTarget::ETargetType type = AliRsnTarget::kTargetTypes);
   AliRsnValue(const AliRsnValue& copy);
   AliRsnValue& operator=(const AliRsnValue& copy);
   virtual ~AliRsnValue() { }

   TArrayD&        GetArray()               {return fBinArray;}
   const Double_t* GetArrayValues() const   {return fBinArray.GetArray();}
   Double_t        GetComputedValue() const {return fComputedValue;}

   void            SetUseMCInfo(Bool_t yn = kTRUE) {fUseMCInfo = yn;}
   void            SetBins(Int_t n, Double_t min, Double_t max);
   void            SetBins(Int_t n, Double_t *array);
   void            SetBins(Double_t min, Double_t max, Double_t step);

   virtual Bool_t  Eval(TObject *object);
   virtual void    Print(Option_t *option = "") const;

protected:

   Bool_t   fUseMCInfo;      // flag to choose MC info when choice is possible
   Double_t fComputedValue;  // computed value
   TArrayD  fBinArray;       // array of bins (when used for a histogram axis)

   ClassDef(AliRsnValue, 3)  // AliRsnValue base class
};

#endif
