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
#include <vector>

using std::vector;

class AliRawReaderMemory: public AliRawReader {
  public :
    AliRawReaderMemory();
    AliRawReaderMemory(UChar_t* memory, UInt_t size);
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

    virtual Int_t    GetEquipmentSize() const;
    virtual Int_t    GetEquipmentType() const {return 0;};
    virtual Int_t    GetEquipmentId() const;
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

    void             SetEquipmentID(Int_t id);

    Bool_t AddBuffer(UChar_t* memory, ULong_t size, Int_t equipmentId );
    void ClearBuffers();

  protected :

  private:
    class AliRRMBuffer {
    public:
      AliRRMBuffer();
      AliRRMBuffer(UChar_t* pBuffer, UInt_t bufferSize, Int_t equipmentId);
      ~AliRRMBuffer();
      AliRRMBuffer(const AliRRMBuffer& src);
      AliRRMBuffer& operator=(const AliRRMBuffer& src);

      UChar_t* GetBuffer() const {return fBuffer;}
      UInt_t GetBufferSize() const {return fBufferSize;}
      Int_t GetEquipmentId() const {return fEquipmentId;}
      void SetEquipmentId(Int_t equipmentId) {fEquipmentId=equipmentId;}
    private:
      UChar_t*         fBuffer;       //! buffer for payload
      UInt_t           fBufferSize;   //! size of fBuffer in bytes
      Int_t            fEquipmentId;  //! Equipment id
    };

    AliRawReaderMemory(const AliRawReaderMemory& rawReader);
    AliRawReaderMemory& operator = (const AliRawReaderMemory& rawReader);

    Bool_t OpenNextBuffer();

    UInt_t                          fPosition;      //! Current position in current buffer
    vector<AliRRMBuffer>            fBuffers;       //! Current buffer descriptor
    UInt_t                          fCurrent;       //! Current buffer index

    ClassDef(AliRawReaderMemory, 0) // class for reading raw digits from a memory block
};

#endif

