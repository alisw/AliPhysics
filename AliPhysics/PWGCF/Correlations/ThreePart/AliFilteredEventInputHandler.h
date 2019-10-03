#ifndef ALIFILTEREDEVENTINPUTHANDLER_H
#define ALIFILTEREDEVENTINPUTHANDLER_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     AliFilteredEvent Input Handler realisation of the AliVEventHandler interface
//     Author: Paul Baetzing, UIO
//-------------------------------------------------------------------------

#include "AliInputEventHandler.h"
#include "AliFilteredEvent.h"
class TList;
class TH2F;


class AliFilteredEventInputHandler : public AliInputEventHandler {
  public:
    AliFilteredEventInputHandler();
    AliFilteredEventInputHandler(const char* name, const char* title);
    virtual ~AliFilteredEventInputHandler();
    virtual Bool_t       Init(Option_t* /*opt*/) {return kTRUE;}
    virtual Bool_t       Init(TTree* tree, Option_t* opt);
    AliFilteredEvent     *GetEvent() const {return fEvent;}
    virtual Bool_t       BeginEvent(Long64_t entry);
    virtual Bool_t       Notify() { return AliVEventHandler::Notify();};
    virtual Bool_t       Notify(const char* path);
    virtual Bool_t       FinishEvent();
    Option_t            *GetDataType() const;
 private:
    AliFilteredEventInputHandler(const AliFilteredEventInputHandler& handler);             
    AliFilteredEventInputHandler& operator=(const AliFilteredEventInputHandler& handler);      
    ClassDef(AliFilteredEventInputHandler, 1);
 private:
    AliFilteredEvent    *fEvent;   //! Pointer to the event
};


#endif