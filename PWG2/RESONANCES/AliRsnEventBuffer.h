/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 **************************************************************************/

//
// Class AliRsnEventBuffer
//
// Implements a temporary buffer of many AliRsnEvent objects
// which is useful for event mixing.
//
// author: Martin Vala (Martin.Vala@cern.ch)
// revised by: A. Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#ifndef AliRsnEventBuffer_h
#define AliRsnEventBuffer_h

#include "AliRsnEvent.h"

class AliRsnCutSet;

class AliRsnEventBuffer : public TObject
{
  public:

    AliRsnEventBuffer(Int_t buffsize = 1000, Bool_t deleteBufferWhenReset = kFALSE);
    ~AliRsnEventBuffer();

    void ClearBuffer();
    void ResetIndex();

    void AddEvent(AliRsnEvent *event);
    Int_t IndexOf(AliRsnEvent *event);
    AliRsnEvent *GetEvent(Int_t index) ;
    AliRsnEvent *GetCurrentEvent();
    AliRsnEvent *GetNextEvent();
    AliRsnEvent *GetNextGoodEvent(Int_t &start, AliRsnCutSet *cuts = 0);

    void SetEventsBufferSize(const Int_t& theValue) { fEventsBufferSize = theValue; }
    Int_t GetEventsBufferSize() const { return fEventsBufferSize; }
    Int_t GetEventsBufferIndex() const { return fEventsBufferIndex; }

    void SetDeleteBufferWhenReset(const Bool_t& theValue = kTRUE) { fDeleteBufferWhenReset = theValue; }
    Bool_t GetDeleteBufferWhenReset() const { return fDeleteBufferWhenReset; }

    Int_t NEmptySlots();

  private:

    AliRsnEventBuffer(const AliRsnEventBuffer& buf) :
        TObject(buf), fDeleteBufferWhenReset(0),fEventsBufferSize(0),fEventsBufferIndex(0) {}
    const AliRsnEventBuffer& operator=(const AliRsnEventBuffer& /*buf*/) {return (*this);}

    Bool_t       fDeleteBufferWhenReset;  // flag if buffer should be deleted when reset is done
    Int_t        fEventsBufferSize;       // buffer size
    Int_t        fEventsBufferIndex;      // current buffer index
    AliRsnEvent *fEventsBuffer[10000];    // array of events

    ClassDef(AliRsnEventBuffer, 1)
};

#endif
