#ifndef ALIRAWREADERROOT_H
#define ALIRAWREADERROOT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>
#include <Riostream.h>
#include "AliRawEvent.h"
#include "AliRawReader.h"


class AliRawReaderRoot: public AliRawReader {
  public :
    AliRawReaderRoot(const char* fileName, Int_t eventNumber);
    AliRawReaderRoot(AliRawEvent* event);
    virtual ~AliRawReaderRoot();

    virtual UInt_t   GetType();
    virtual UInt_t   GetRunNumber();
    virtual const UInt_t* GetEventId();
    virtual const UInt_t* GetTriggerPattern();
    virtual const UInt_t* GetDetectorPattern();
    virtual const UInt_t* GetAttributes();
    virtual UInt_t   GetGDCId();

    virtual Bool_t   ReadMiniHeader();
    virtual Bool_t   ReadNextData(UChar_t*& data);

    virtual Bool_t   Reset();

  protected :
    virtual Bool_t   ReadNext(UChar_t* data, Int_t size);

    TFile*           fFile;         // raw data root file
    AliRawEvent*     fEvent;        // (super) event
    Int_t            fSubEventIndex; // index of current sub event
    AliRawEvent*     fSubEvent;     // current sub event
    AliRawData*      fRawData;      // current raw data
    UChar_t*         fPosition;     // current position in the raw data
    UChar_t*         fEnd;          // end position of the current subevent

    ClassDef(AliRawReaderRoot, 0) // class for reading raw digits from a root file
};

#endif
