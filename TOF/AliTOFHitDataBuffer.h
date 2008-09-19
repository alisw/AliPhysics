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

//class AliTOFHitData;

//data buffer size
#define DATA_BUFFER_SIZE 10000

class AliTOFHitDataBuffer : public TObject{
  
 public:

  AliTOFHitDataBuffer(Int_t BufferSize = DATA_BUFFER_SIZE);
  ~AliTOFHitDataBuffer();  
  AliTOFHitDataBuffer(const AliTOFHitDataBuffer &source);
  AliTOFHitDataBuffer& operator=(const AliTOFHitDataBuffer & source); 
  void           Reset() {fEntries = 0;};
  Bool_t         Add(AliTOFHitData &HitData);
  
  AliTOFHitData *GetBuffer() {return fBuffer;};
  AliTOFHitData *GetHit(Int_t Hit) const {return (Hit < fBufferSize ? &fBuffer[Hit] : 0x0);};
  Int_t          GetEntries() const {return fEntries;};
  
 private:
  Int_t          fBufferSize; // buffer size
  AliTOFHitData *fBuffer;     // buffer
  Int_t          fEntries;    // entries

  ClassDef(AliTOFHitDataBuffer, 1);
};

#endif

