#ifndef ALIVIRTUALEVENTHANDLER_H
#define ALIVIRTUALEVENTHANDLER_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     Event Handler base class
//     Author: Andreas Morsch, CERN
//-------------------------------------------------------------------------

#include <TNamed.h>


class AliVirtualEventHandler : public TNamed {

 public:
    AliVirtualEventHandler();
    AliVirtualEventHandler(const char* name, const char* title);
    virtual ~AliVirtualEventHandler();
    virtual void         SetOutputFileName(char* fname)  = 0;
    virtual char*        GetOutputFileName()             = 0;
    virtual Bool_t       InitIO(Option_t* opt)           = 0;
    virtual Bool_t       BeginEvent()                    = 0;
    virtual Bool_t       FinishEvent()                   = 0;
    virtual Bool_t       Terminate()                     = 0;
    virtual Bool_t       TerminateIO()                   = 0;
 private :
  ClassDef(AliVirtualEventHandler, 1);
};

#endif
