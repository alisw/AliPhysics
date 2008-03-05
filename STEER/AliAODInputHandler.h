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

class AliAODInputHandler : public AliInputEventHandler {

 public:
    AliAODInputHandler();
    AliAODInputHandler(const char* name, const char* title);
    virtual ~AliAODInputHandler();
    virtual Bool_t       Init(Option_t* /*opt*/) {return kTRUE;}
    virtual Bool_t       Init(TTree* tree, Option_t* opt);
    AliAODEvent         *GetEvent() const {return fEvent;}
 private:
    AliAODInputHandler(const AliAODInputHandler& handler);             
    AliAODInputHandler& operator=(const AliAODInputHandler& handler);  
 private:
    AliAODEvent    *fEvent;   //! Pointer to the event 
    ClassDef(AliAODInputHandler, 1);
};

#endif
