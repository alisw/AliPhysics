#ifndef ALIESDINPUTHANDLER_H
#define ALIESDINPUTHANDLER_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     ESD Input Handler realisation of the AliVEventHandler interface
//     Author: Andreas Morsch, CERN
//-------------------------------------------------------------------------

#include "AliInputEventHandler.h"
#include "AliESDEvent.h"

class AliESDInputHandler : public AliInputEventHandler {

 public:
    AliESDInputHandler();
    AliESDInputHandler(const char* name, const char* title);
    virtual ~AliESDInputHandler();
    virtual Bool_t       InitIO(Option_t* opt);
    virtual Bool_t       BeginEvent(Long64_t entry);
    AliESDEvent         *GetEvent() const {return fEvent;}
    private:
    AliESDEvent    *fEvent;   //! Pointer to the event 
    ClassDef(AliESDInputHandler, 1);
};

#endif
