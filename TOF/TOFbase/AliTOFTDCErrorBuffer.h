#ifndef ALITOFTDCERRORBUFFER_H
#define ALITOFTDCERRORBUFFER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  author: Roberto Preghenella (R+), preghenella@bo.infn.it
*/

///////////////////////////////////////////////////////////////
//                                                           //
//   This class provides a buffer for TDC errors.            //
//                                                           //
///////////////////////////////////////////////////////////////

#include "TObject.h"
#include "AliTOFTDCError.h"
#include "TClonesArray.h"

class AliTOFTDCErrorBuffer : 
public TObject
{
 public:
  AliTOFTDCErrorBuffer(); //default constructor
  AliTOFTDCErrorBuffer(const AliTOFTDCErrorBuffer &source) : TObject(source), fBuffer(source.fBuffer) {}; //copy contructor
  AliTOFTDCErrorBuffer &operator = (const AliTOFTDCErrorBuffer &source) {
    if (&source != this) {
      TObject::operator=(source);
      fBuffer = source.fBuffer;
    }
    return *this;
  }; //operator =
  virtual ~AliTOFTDCErrorBuffer(); //destructor

  void Reset() {fBuffer.Clear();}; // reset
  void Add(const AliTOFTDCError &err); //add hit
  TClonesArray *GetBuffer() {return &fBuffer;}; //get buffer
  Int_t GetEntries() const {return fBuffer.GetEntries();}; //get entries
  AliTOFTDCError *GetError(Int_t ierr) const {return (ierr < GetEntries() ? (AliTOFTDCError *)fBuffer.At(ierr) : 0x0);}; //get error

 private:

  TClonesArray fBuffer; // buffer

  ClassDef(AliTOFTDCErrorBuffer, 1);
};

#endif /* ALITOFTDCERRORBUFFER_H */
