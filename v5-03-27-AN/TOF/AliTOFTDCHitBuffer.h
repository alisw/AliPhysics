#ifndef ALITOFTDCHITBUFFER_H
#define ALITOFTDCHITBUFFER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  author: Roberto Preghenella (R+), preghenella@bo.infn.it
*/

///////////////////////////////////////////////////////////////
//                                                           //
//   This class provides a buffer for TDC hits.              //
//                                                           //
///////////////////////////////////////////////////////////////

#include "TObject.h"
#include "AliTOFTDCHit.h"
#include "TClonesArray.h"

class AliTOFTDCHitBuffer : 
public TObject
{
 public:
  AliTOFTDCHitBuffer(); //default constructor
  AliTOFTDCHitBuffer(const AliTOFTDCHitBuffer &source) : TObject(source), fBuffer(source.fBuffer) {}; //copy constructor
  AliTOFTDCHitBuffer &operator = (const AliTOFTDCHitBuffer &source) {
    if (&source != this) {
      TObject::operator=(source);
      fBuffer = source.fBuffer;
    }
    return *this;
  }; //operator =
  virtual ~AliTOFTDCHitBuffer(); //destructor

  void Reset() {fBuffer.Clear();}; // reset
  void Add(const AliTOFTDCHit &Hit); //add hit
  TClonesArray *GetBuffer() {return &fBuffer;}; //get buffer
  Int_t GetEntries() const {return fBuffer.GetEntries();}; //get entries
  AliTOFTDCHit *GetHit(Int_t Hit) const {return (Hit < GetEntries() ? (AliTOFTDCHit *)fBuffer.At(Hit) : 0x0);}; //get hit

 private:
  
  TClonesArray fBuffer; // buffer
  
  ClassDef(AliTOFTDCHitBuffer, 1);
};

#endif /* ALITOFTDCHITBUFFER_H */
