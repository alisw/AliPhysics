#ifndef ALIFSTREAM_H
#define ALIFSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
// This is the class which is to be used during the writing of
// simulated raw data (DDL files format).
// It is using the root functionality in order to deal correctly
// with little/big endian issue. By convention the detector raw
// data payload is stored always with little endian (this corresponds
// to the real life situation when the detector data is coming from
// the hardware).
//-------------------------------------------------------------------------

#include <TObject.h>

class AliFstream : public TObject {
public:
  AliFstream();
  AliFstream(const char *fileName);
  virtual ~AliFstream();

  void   Seekp(UInt_t position);
  UInt_t Tellp();
  void   WriteBuffer(const char *buffer, UInt_t size, Bool_t force = kFALSE);

private:

  AliFstream(const AliFstream &source);
  AliFstream &operator =(const AliFstream& source);

  fstream *fFile;       // Output file stream
  UChar_t *fBuffer;     // Pointer to the internal buffer
  UInt_t   fBufferSize; // Internal buffer size
  Bool_t   fSwap;       // Big or little endian

  ClassDef(AliFstream,0)
};

#endif
