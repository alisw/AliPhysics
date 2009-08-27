#ifndef ALIMULTIEVENTINPUTHANDLER_H
#define ALIMULTIEVENTINPUTHANDLER_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//----------------------------------------------------------------------------
//     Multi VEvent Input Handler realisation of the AliVEventHandler interface.
//     This class handles multiple events for mixing.
//     Author: Andreas Morsch, CERN
//----------------------------------------------------------------------------

#include "AliInputEventHandler.h"
class AliVEventPool;
class AliVEvent;

class AliMultiEventInputHandler : public AliInputEventHandler {

 public:
    AliMultiEventInputHandler();
    AliMultiEventInputHandler(Int_t size, Int_t format = 1);
    AliMultiEventInputHandler(const char* name, const char* title, Int_t size, Int_t format = 1);
    virtual ~AliMultiEventInputHandler();
    void   SetBufferSize(Int_t size) {fBufferSize = size;}
    void   SetEventPool(AliVEventPool* pool) {fEventPool = pool;}
    Int_t  GetBufferSize()           const {return fBufferSize;}
    Int_t  GetNBuffered()            const {return fNBuffered;}
    Bool_t IsBufferReady()           const {return (fNBuffered >= (fBufferSize -1));}
    Bool_t IsFreshBuffer()           const {return (fIndex == (fBufferSize - 1));}
    AliVEventPool           *GetEventPool()      const {return fEventPool;}
    virtual AliVEvent       *GetEvent()          const {return GetLatestEvent();}
    virtual AliVEvent       *GetEvent(Int_t iev) const;
    AliVEvent               *GetLatestEvent()    const {return fEventBuffer[fIndex];}
    // From the interface
    virtual Bool_t Init(Option_t* /*opt*/)    {return kTRUE;}
    virtual Bool_t Init(TTree* tree, Option_t* /*opt*/);
    virtual Bool_t FinishEvent();
    virtual Bool_t BeginEvent(Long64_t /*entry*/);
    
 private:
    AliMultiEventInputHandler(const AliMultiEventInputHandler& handler);             
    AliMultiEventInputHandler& operator=(const AliMultiEventInputHandler& handler);  
 private:
    Int_t          fBufferSize;   // Size of the buffer
    Int_t          fFormat;       // 0: ESD 1: AOD
    Int_t          fNBuffered;    // Number of events actually buffered
    Int_t          fIndex;        // Pointer to most recent event
    Int_t          fCurrentBin;   // Current bin from the pool
    TTree*         fTree;         // Pointer to the tree
    AliVEventPool* fEventPool;    // Pointer to the pool
    AliVEvent**    fEventBuffer;  // The event buffer
    ClassDef(AliMultiEventInputHandler, 1);
};

#endif
