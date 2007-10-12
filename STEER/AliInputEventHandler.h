#ifndef ALIINPUTEVENTHANDLER_H
#define ALIINPUTEVENTHANDLER_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     Input Handler realisation of the AliVEventHandler interface
//     Author: Andreas Morsch, CERN
//-------------------------------------------------------------------------

#include "AliVEventHandler.h"

class TChain;
class AliVEvent;

class AliInputEventHandler : public AliVEventHandler {

 public:
    AliInputEventHandler();
    AliInputEventHandler(const char* name, const char* title);
    virtual ~AliInputEventHandler();
    virtual void         SetOutputFileName(char* /*fname*/) {;}
    virtual char        *GetOutputFileName()                {return 0;}
    virtual Bool_t       InitIO(Option_t* /*opt*/)          {return kTRUE;}
    virtual Bool_t       BeginEvent()                       {return kTRUE;}
    virtual Bool_t       Notify(const char */*path*/)       {return kTRUE;}
    virtual Bool_t       FinishEvent()                      {return kTRUE;}        
    virtual Bool_t       Terminate()                        {return kTRUE;}
    virtual Bool_t       TerminateIO()                      {return kTRUE;}
    // Setters
    virtual void         SetInputChain(TChain* chain)       {fChain = chain;}
    // Getters
    virtual AliVEvent *GetEvent() const                     {return fEvent;}
    virtual TChain    *GetChain() const                     {return fChain;}
 protected:
    AliVEvent    *fEvent;   //! Pointer to the event 
    TChain       *fChain;   //! Pointer to the chain
    ClassDef(AliInputEventHandler, 1);
};

#endif
