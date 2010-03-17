//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTDCSARRAY_H
#define ALIHLTDCSARRAY_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTDCSArray.h
    @author Haavard Helstrup
    @date   
    @brief  
*/

#include "TArrayF.h"

class AliHLTDCSArray : public TObject {

 public:
   AliHLTDCSArray(Int_t entries=0);
   virtual ~AliHLTDCSArray();
   void SetTime(UInt_t time) {fTimeStamp=time;}
   UInt_t GetTime() const {return fTimeStamp;}
   void SetValue(Int_t i, Float_t val) {fValues[i]=val;}
   Float_t GetValue(Int_t i) const {return fValues[i];}
   
 protected:  
   UInt_t      fTimeStamp;    // time of DCS reading
   TArrayF     fValues;       // array holding DCS readings
   
 ClassDef(AliHLTDCSArray,1)
};

#endif   

 
