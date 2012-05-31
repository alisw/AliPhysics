#ifndef ALIACORDERAWSTREAM_H
#define ALIACORDERAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliACORDERawStream.h 20210 2007-08-18 08:41:30Z hristov $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Reads ACORDE DDL raw data from raw data stream                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"

class AliRawReader;

class AliACORDERawStream : public TObject {

 public:

  AliACORDERawStream(AliRawReader*);
  AliACORDERawStream(const AliACORDERawStream &r); 
  virtual ~AliACORDERawStream();
  AliACORDERawStream &operator=(const AliACORDERawStream &r);
  //MRC's part 
  Int_t GetNEvents(char* fileName);
  //
  virtual void    Reset();
  virtual Bool_t  Next();

  Int_t           DataSize() const { return fDataSize; }
  UInt_t          GetWord(Int_t index) const;

  enum EACORDERawStreamError {
      kRawDataSizeErr = 1
  };

 private:

  UInt_t          GetNextWord();

  AliRawReader*   fRawReader;    // object for reading the raw data
  Int_t           fPosition;     // current position in the raw-data payload
  UChar_t*        fData;         // pointer to raw data payload

  Int_t           fDataSize;     // data size

  UInt_t          fWord[4];      // data vector

  ClassDef(AliACORDERawStream,0) // class for reading ACORDE DDL raw data

};

typedef AliACORDERawStream AliCRTRawStream; // for backward compatibility

#endif
