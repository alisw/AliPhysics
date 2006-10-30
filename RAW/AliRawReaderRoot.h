#ifndef ALIRAWREADERROOT_H
#define ALIRAWREADERROOT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///
/// This is a class for reading raw data from a root file.
///
///////////////////////////////////////////////////////////////////////////////

#include "AliRawReader.h"

class AliRawEvent;
class AliRawEquipment;
class AliRawData;
class TFile;
class TBranch;


class AliRawReaderRoot: public AliRawReader {
  public :
    AliRawReaderRoot(const char* fileName, Int_t eventNumber = -1);
    AliRawReaderRoot(AliRawEvent* event);
    AliRawReaderRoot(const AliRawReaderRoot& rawReader);
    AliRawReaderRoot& operator = (const AliRawReaderRoot& rawReader);
    virtual ~AliRawReaderRoot();

    virtual const AliRawEventHeaderBase* GetEventHeader() const;

    virtual UInt_t   GetType() const;
    virtual UInt_t   GetRunNumber() const;
    virtual const UInt_t* GetEventId() const;
    virtual const UInt_t* GetTriggerPattern() const;
    virtual const UInt_t* GetDetectorPattern() const;
    virtual const UInt_t* GetAttributes() const;
    virtual const UInt_t* GetSubEventAttributes() const;
    virtual UInt_t   GetLDCId() const;
    virtual UInt_t   GetGDCId() const;

    virtual Int_t    GetEquipmentSize() const;
    virtual Int_t    GetEquipmentType() const;
    virtual Int_t    GetEquipmentId() const;
    virtual const UInt_t* GetEquipmentAttributes() const;
    virtual Int_t    GetEquipmentElementSize() const;
    virtual Int_t    GetEquipmentHeaderSize() const;

    virtual Bool_t   ReadHeader();
    virtual Bool_t   ReadNextData(UChar_t*& data);
    virtual Bool_t   ReadNext(UChar_t* data, Int_t size);

    virtual Bool_t   Reset();

    virtual Bool_t   NextEvent();
    virtual Bool_t   RewindEvents();

    virtual Int_t    CheckData() const;

  protected :
    TFile*           fFile;         // raw data root file
    TBranch*         fBranch;       // branch of raw events
    Int_t            fEventIndex;   // index of the event in the tree
    AliRawEvent*     fEvent;        // (super) event
    Int_t            fSubEventIndex; // index of current sub event
    AliRawEvent*     fSubEvent;     // current sub event
    Int_t            fEquipmentIndex; // index of current equipment
    AliRawEquipment* fEquipment;    // current equipment
    AliRawData*      fRawData;      // current raw data
    UChar_t*         fPosition;     // current position in the raw data
    UChar_t*         fEnd;          // end position of the current subevent

    ClassDef(AliRawReaderRoot, 0) // class for reading raw digits from a root file
};

#endif
