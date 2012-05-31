#ifndef ALIACORDERAWDATA_H
#define ALIACORDERAWDATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Converts ACORDE digits into a raw data stream                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"

class AliACORDERawData : public TObject {

 public:

  AliACORDERawData();
  AliACORDERawData(const AliACORDERawData &r); 
  virtual ~AliACORDERawData();
  AliACORDERawData &operator=(const AliACORDERawData &r);      // ass. op.

  void WriteACORDERawData(Bool_t *b,Bool_t multi);
  void SetACORDERawWords(Bool_t *b,Bool_t multi);
  

 private:

  UInt_t fWord9;
  UInt_t fWord10;
  UInt_t fWord11;
  UInt_t fWord12;

  ClassDef(AliACORDERawData,1)             //  ACORDE raw data class

};

typedef AliACORDERawData AliCRTRawData; // for backward compatibility

#endif
