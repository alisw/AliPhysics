#ifndef ALIAODINPUTHANDLER_H
#define ALIAODINPUTHANDLER_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     AOD Input Handler realisation of the AliVEventHandler interface
//     Author: Andreas Morsch, CERN
//-------------------------------------------------------------------------

#include "AliInputEventHandler.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
class TList;
class AliMCEvent;
class TH2F;
class AliMCEvent;


class AliAODInputHandler : public AliInputEventHandler {

 public:
    AliAODInputHandler();
    AliAODInputHandler(const char* name, const char* title);
    virtual ~AliAODInputHandler();
    virtual Bool_t       Init(Option_t* /*opt*/) {return kTRUE;}
    virtual Bool_t       Init(TTree* tree, Option_t* opt);
    AliAODEvent         *GetEvent() const {return fEvent;}
    AliMCEvent          *MCEvent()  const {return fMCEvent;}
    virtual void         AddFriend(char* filename);
    virtual Bool_t       BeginEvent(Long64_t entry);
    virtual Bool_t       Notify() { return AliVEventHandler::Notify();};
    virtual Bool_t       Notify(const char* path);
    virtual Bool_t       FinishEvent();
    Option_t            *GetDataType() const;
    // Get the statistics object (currently TH2). Option can be BIN0.
    virtual TObject     *GetStatistics(Option_t *option="") const;
    // Provisions for event merging
    void                 SetMergeEvents(Bool_t flag) {fMergeEvents = flag;}
    Bool_t               GetMergeEvents() const {return fMergeEvents;}
    AliAODEvent*         GetEventToMerge() {return fAODEventToMerge;}
    void                 SetMergeOffset(Int_t ioff) {fMergeOffset = ioff;}
 private:
    void ConnectFriends();
    AliAODInputHandler(const AliAODInputHandler& handler);             
    AliAODInputHandler& operator=(const AliAODInputHandler& handler);  
 private:
    AliAODEvent    *fEvent;   //! Pointer to the event
    AliMCEvent     *fMCEvent; //! Pointer to the MCEvent
    TList          *fFriends; //  List of friend trees 
// Support for merged events
    Bool_t          fMergeEvents;     // Flag for event merging
    Bool_t          fFriendsConnected;// Friends are connected
    TFile          *fFileToMerge;     //! File for merging
    TTree          *fTreeToMerge;     //! Tree for merging
    AliAODEvent    *fAODEventToMerge; //! Event for merging
    Int_t           fMergeOffset;     //! Event offset for merging
    TH2F*           fHistStatistics[2]; //! how many events are cut away why {all,bin 0}
    ClassDef(AliAODInputHandler, 3);
};

#endif
