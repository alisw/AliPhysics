#ifndef ALITOFHITDATABUFFER_H
#define ALITOFHITDATABUFFER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////
//                                                           //
//   This class provides the key-reading for TOF raw data.   //
//                                                           //
///////////////////////////////////////////////////////////////

#include "TObject.h"
#include "AliTOFHitData.h"
#include "TClonesArray.h"

class AliTOFHitDataBuffer : 
public TObject
{
  
 public:

  AliTOFHitDataBuffer(); // default constructor 
  AliTOFHitDataBuffer(Int_t size); // overloaded constructor
  ~AliTOFHitDataBuffer();   // default destructor
  AliTOFHitDataBuffer(const AliTOFHitDataBuffer &source) : TObject(source), fBuffer(source.fBuffer) {}; // copy constructor 
  AliTOFHitDataBuffer& operator=(const AliTOFHitDataBuffer & source) {
    if (&source != this) {
      TObject::operator=(source);
      fBuffer = source.fBuffer;
    }
    return *this;
  }; // operator =

  void Reset() {fBuffer.Clear();}; // reset
  Bool_t Add(AliTOFHitData &HitData); // add
  
  TClonesArray *GetBuffer() {return &fBuffer;}; // get buffer
  AliTOFHitData *GetHit(Int_t Hit) const {return (Hit < GetEntries() ? (AliTOFHitData *)fBuffer.At(Hit) : 0x0);}; // get hit
  Int_t GetEntries() const {return fBuffer.GetEntries();}; // get entries
  
 private:

  TClonesArray fBuffer; // buffer

  ClassDef(AliTOFHitDataBuffer, 1);
};

#endif

