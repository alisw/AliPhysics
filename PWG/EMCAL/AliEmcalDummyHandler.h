#ifndef ALIEMCALDUMMYHANDLER_H
#define ALIEMCALDUMMYHANDLER_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliEmcalDummyHandler.h 58219 2012-08-17 14:18:19Z agheata $ */

//-------------------------------------------------------------------------
//     Dummy implementation of the input handler e.g. for the case of event-loop steered analysis
//     Author: M. Verweij
//-------------------------------------------------------------------------

#include "AliDummyHandler.h"

class AliEmcalDummyHandler : public AliDummyHandler {

 public:
    AliEmcalDummyHandler();
    AliEmcalDummyHandler(const char* name, const char* title);
    virtual ~AliEmcalDummyHandler();

    void                 SetEvent(AliVEvent *e)  { fEvent = e    ; } 
    AliVEvent           *GetEvent()        const { return fEvent ; }


 protected:
    AliVEvent    *fEvent;         //! Pointer to the event

 private:
    AliEmcalDummyHandler(const AliEmcalDummyHandler& handler);             
    AliEmcalDummyHandler& operator=(const AliEmcalDummyHandler& handler);  

    ClassDef(AliEmcalDummyHandler, 1);
};

#endif
