#ifndef ALIRAWREADERDATE_H
#define ALIRAWREADERDATE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliRawReader.h"

struct eventHeaderStruct;


class AliRawReaderDate: public AliRawReader {
  public :
    AliRawReaderDate(void* event);
    AliRawReaderDate(const AliRawReaderDate& rawReader);
    AliRawReaderDate& operator = (const AliRawReaderDate& rawReader);
    virtual ~AliRawReaderDate() {};

    virtual UInt_t   GetType() const;
    virtual UInt_t   GetRunNumber() const;
    virtual const UInt_t* GetEventId() const;
    virtual const UInt_t* GetTriggerPattern() const;
    virtual const UInt_t* GetDetectorPattern() const;
    virtual const UInt_t* GetAttributes() const;
    virtual UInt_t   GetGDCId() const;

    virtual Bool_t   ReadMiniHeader();
    virtual Bool_t   ReadNextData(UChar_t*& data);

    virtual Bool_t   Reset();

    virtual Int_t    CheckData() const;

  protected :
    virtual Bool_t   ReadNext(UChar_t* data, Int_t size);

    eventHeaderStruct* fEvent;      // raw data super event
    eventHeaderStruct* fSubEvent;   // raw data sub event

    UChar_t*         fPosition;     // current position in the raw data
    UChar_t*         fEnd;          // end position of the current subevent

    ClassDef(AliRawReaderDate, 0) // class for reading raw digits from a root file
};

#endif
