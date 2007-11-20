#ifndef ALIRAWREADERMEMORY_H
#define ALIRAWREADERMEMORY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///
/// This is a class for reading raw data memory.
///
///////////////////////////////////////////////////////////////////////////////

#include "AliRawReader.h"
#ifdef __CINT__
class fstream;
#else
#include <Riostream.h>
#endif
#include <TString.h>


class AliRawReaderMemory: public AliRawReader {
  public :
    AliRawReaderMemory();
    AliRawReaderMemory(UChar_t* memory, UInt_t size);
/*     AliRawReaderMemory(const AliRawReaderMemory& rawReader); */
/*     AliRawReaderMemory& operator = (const AliRawReaderMemory& rawReader); */
    virtual ~AliRawReaderMemory();

    virtual void     RequireHeader(Bool_t required);

    virtual UInt_t   GetType() const {return 0;};
    virtual UInt_t   GetRunNumber() const {return 0;};
    virtual const UInt_t* GetEventId() const {return 0;};
    virtual const UInt_t* GetTriggerPattern() const {return 0;};
    virtual const UInt_t* GetDetectorPattern() const {return 0;};
    virtual const UInt_t* GetAttributes() const {return 0;};
    virtual const UInt_t* GetSubEventAttributes() const {return 0;};
    virtual UInt_t   GetLDCId() const {return 0;};
    virtual UInt_t   GetGDCId() const {return 0;};
    virtual UInt_t   GetTimestamp() const {return 0;};

    virtual Int_t    GetEquipmentSize() const {return fBufferSize;};
    virtual Int_t    GetEquipmentType() const {return 0;};
    virtual Int_t    GetEquipmentId() const {return fEquipmentId;};
    virtual const UInt_t* GetEquipmentAttributes() const {return NULL;};
    virtual Int_t    GetEquipmentElementSize() const {return 0;};
    virtual Int_t    GetEquipmentHeaderSize() const {return 0;};

    virtual Bool_t   ReadHeader();
    virtual Bool_t   ReadNextData(UChar_t*& data);
    virtual Bool_t   ReadNext(UChar_t* data, Int_t size);

    virtual Bool_t   Reset();

    virtual Bool_t   NextEvent();
    virtual Bool_t   RewindEvents();

    virtual Bool_t   SetMemory( UChar_t* memory, ULong_t size );

    void             SetEquipmentID(Int_t id) { fEquipmentId = id; }

  protected :

    UChar_t*         fBuffer;      // buffer for payload
    UInt_t           fBufferSize;  // size of fBuffer in bytes
    UInt_t           fPosition;    // Current position in memory
    Int_t            fEquipmentId;  // Equipment id provided by the user

    ClassDef(AliRawReaderMemory, 0) // class for reading raw digits from a memory block

  private:

    AliRawReaderMemory(const AliRawReaderMemory& rawReader);
    AliRawReaderMemory& operator = (const AliRawReaderMemory& rawReader);
};

#endif

