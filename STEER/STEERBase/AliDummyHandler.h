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
    Option_t            *GetDataType() const { return "MC"; }

 private:
    AliDummyHandler(const AliDummyHandler& handler);             
    AliDummyHandler& operator=(const AliDummyHandler& handler);  

    ClassDef(AliDummyHandler, 1);
};

#endif
