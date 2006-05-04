#ifndef ALICTPRAWSTREAM_H
#define ALICTPRAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
///
/// This is a class for reading the CTP (trigger) DDL raw data
/// The format of the raw data is taken form the trigger TDR
///
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliRawReader;

class AliCTPRawStream: public TObject {
  public :
    AliCTPRawStream(AliRawReader* rawReader);
    virtual ~AliCTPRawStream();

    virtual void             Reset();
    virtual Bool_t           Next();

    inline ULong64_t GetClassMask()   const { return fClassMask; }  // Provide the trigger class mask
    inline UChar_t   GetClusterMask() const { return fClusterMask; }// Provide the trigger cluster mask

  protected:
    AliCTPRawStream(const AliCTPRawStream& stream);
    AliCTPRawStream& operator = (const AliCTPRawStream& stream);

  private:

    ULong64_t        fClassMask;   // trigger class mask
    UChar_t          fClusterMask; // trigger cluster mask

    AliRawReader*    fRawReader;   // object for reading the raw data

    enum {kCTPIndex = 17};         //CTP detector index (AliDAQConfig.h)

    ClassDef(AliCTPRawStream, 0)   // class for reading CTP DDL raw data
};

#endif
