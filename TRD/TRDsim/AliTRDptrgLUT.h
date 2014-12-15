#ifndef AliTRDPTRGLUT_H
#define AliTRDPTRGLUT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// --------------------------------------------------------
// 
// PTRG look up table definition
//
// --------------------------------------------------------

#include "TObject.h"

class AliRunLoader;

class AliTRDptrgLUT : public TObject {
 public:
  AliTRDptrgLUT();
  ~AliTRDptrgLUT();
  Int_t LookUp(UInt_t input); // execute a look up
  Int_t InitTable(Int_t inputWidth, Int_t outputWidth, Int_t *tableData, 
                  Bool_t copy); // load look up table
  Int_t GetInputWidth() { return fInputWidth; } // getter function
 protected:
  Int_t *fLUTData; // lut data storage
  Int_t fInputWidth; // bit width of the input vector
  Int_t fOutputWidth; // bit width of the output vector 
  Int_t fTableEntryCount; // table entry count
  Bool_t fCopiedTable; // indicates whether the look up table was copied
 private:
  AliTRDptrgLUT& operator=(const AliTRDptrgLUT &rhs); // not implemented
  AliTRDptrgLUT(const AliTRDptrgLUT &rhs); // not implemented

  ClassDef(AliTRDptrgLUT, 1);
};

#endif
