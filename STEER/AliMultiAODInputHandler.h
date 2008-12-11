#ifndef ALIMULTIAODINPUTHANDLER_H
#define ALIMULTIAODINPUTHANDLER_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     AOD Input Handler realisation of the AliVEventHandler interface.
//     This class handles multiple events for mixing.
//     Author: Andreas Morsch, CERN
//-------------------------------------------------------------------------

#include "AliInputEventHandler.h"
class AliVEventPool;
class AliAODEvent;

class AliMultiAODInputHandler : public AliInputEventHandler {

 public:
    AliMultiAODInputHandler();
    AliMultiAODInputHandler(Int_t size);
    AliMultiAODInputHandler(const char* name, const char* title, Int_t size);
    virtual ~AliMultiAODInputHandler();
    void   SetBufferSize(Int_t size) {fBufferSize = size;}
    void   SetEventPool(AliVEventPool* pool) {fEventPool = pool;}
    Int_t  GetBufferSize()           const {return fBufferSize;}
    Int_t  GetNBuffered()            const {return fNBuffered;}
    Bool_t IsBufferReady()           const {return (fNBuffered >= (fBufferSize -1));}
    Bool_t IsFreshBuffer()           const {return (fIndex == (fBufferSize - 1));}
    AliVEventPool        *GetEventPool()      const {return fEventPool;}
    virtual AliVEvent    *GetEvent()          const {return 0;}
    virtual AliAODEvent  *GetEvent(Int_t iev) const;
    AliAODEvent          *GetLatestEvent()    const {return fEventBuffer[fIndex];}
    // From the interface
    virtual Bool_t Init(Option_t* /*opt*/)    {return kTRUE;}
    virtual Bool_t Init(TTree* tree, Option_t* /*opt*/);
    virtual Bool_t FinishEvent();
    virtual Bool_t BeginEvent(Long64_t /*entry*/);
    
 private:
    AliMultiAODInputHandler(const AliMultiAODInputHandler& handler);             
    AliMultiAODInputHandler& operator=(const AliMultiAODInputHandler& handler);  
 private:
    Int_t          fBufferSize;   // Size of the buffer
    Int_t          fNBuffered;    // Number of events actually buffered
    Int_t          fIndex;        // Pointer to most recent event
    Int_t          fCurrentBin;   // Current bin from the pool
    TTree*         fTree;         // Pointer to the tree
    AliVEventPool* fEventPool;    // Pointer to the pool
    AliAODEvent**  fEventBuffer;  // The event buffer
    ClassDef(AliMultiAODInputHandler, 1);
};

#endif
