#ifndef ALIDUMMYINPUTHANDLER_H
#define ALIDUMMYINPUTHANDLER_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliDummyHandler.h 58219 2012-08-17 14:18:19Z agheata $ */

//-------------------------------------------------------------------------
//     Dummy implementation of the input handler e.g. for the case of event-loop steered analysis
//     Author: Andrei Gheata, Jan Fiete Grosse-Oetringhaus
//-------------------------------------------------------------------------

#include "AliInputEventHandler.h"

class AliDummyHandler : public AliInputEventHandler {

 public:
    AliDummyHandler();
    AliDummyHandler(const char* name, const char* title);
    virtual ~AliDummyHandler();
    virtual Bool_t       Init(Option_t* opt) { return AliInputEventHandler::Init(opt); }
    virtual Bool_t Init(TTree* tree, Option_t* /*opt*/) { fTree = tree; return kTRUE; }
    Option_t            *GetDataType() const { return "MC"; }

    void                 SetEvent(AliVEvent *e)  { fEvent = e; } 
    AliVEvent           *GetEvent()        const {return fEvent;}

 protected:
    AliVEvent    *fEvent;         // Pointer to the event

 private:
    AliDummyHandler(const AliDummyHandler& handler);             
    AliDummyHandler& operator=(const AliDummyHandler& handler);  

    ClassDef(AliDummyHandler, 3);
};

#endif
