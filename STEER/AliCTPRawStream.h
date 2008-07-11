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
#include <TClonesArray.h>

class AliRawReader;
class AliTriggerIR;

class AliCTPRawStream: public TObject {
  public :
    AliCTPRawStream(AliRawReader* rawReader);
    virtual ~AliCTPRawStream();


    virtual void             Reset();
    virtual Bool_t           Next();

    ULong64_t GetClassMask()   const { return fClassMask; }  // Provide the trigger class mask
    UChar_t   GetClusterMask() const { return fClusterMask; }// Provide the trigger cluster mask

    Int_t     GetNIRs() const { return fIRArray.GetEntriesFast(); }
    const TClonesArray &GetAllIRs() const { return fIRArray; }
    const AliTriggerIR *GetIR(Int_t index) const { return (const AliTriggerIR*)fIRArray.UncheckedAt(index); }


  protected:
    AliCTPRawStream(const AliCTPRawStream& stream);
    AliCTPRawStream& operator = (const AliCTPRawStream& stream);

  private:
    TClonesArray     fIRArray;     // array with trigger interaction records

    ULong64_t        fClassMask;   // trigger class mask
    UChar_t          fClusterMask; // trigger cluster mask

    AliRawReader*    fRawReader;   // object for reading the raw data

    ClassDef(AliCTPRawStream, 0)   // class for reading CTP DDL raw data
};

#endif
