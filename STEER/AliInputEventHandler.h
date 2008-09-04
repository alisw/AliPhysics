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
#include <TTree.h>


class AliVEvent;
class AliRunTag;

class AliInputEventHandler : public AliVEventHandler {

 public:
    AliInputEventHandler();
    AliInputEventHandler(const char* name, const char* title);
    virtual ~AliInputEventHandler();
    virtual void         SetOutputFileName(const char* /*fname*/) {;}
    virtual const char  *GetOutputFileName()                          {return 0;}
    virtual Bool_t       Init(Option_t* /*opt*/)                      {return kTRUE;}
    virtual Bool_t       Init(TTree* /*tree*/, Option_t* /*opt*/)     {return kTRUE;}
    virtual Bool_t       BeginEvent(Long64_t /*entry*/)               {return kTRUE;}
    virtual Bool_t       Notify() { return AliVEventHandler::Notify(); };
    virtual Bool_t       Notify(const char */*path*/)                 {return kTRUE;}
    virtual Bool_t       FinishEvent()                                {return kTRUE;}        
    virtual Bool_t       Terminate()                                  {return kTRUE;}
    virtual Bool_t       TerminateIO()                                {return kTRUE;}
    // Setters
    virtual void         SetInputTree(TTree* tree)                    {fTree = tree;}
     // Getters
    virtual AliVEvent   *GetEvent()        const                      {return 0;}
    virtual AliRunTag   *GetRunTag()       const                      {return 0;}
    virtual Option_t    *GetAnalysisType() const                      {return 0;}
    virtual TTree       *GetTree( )        const                      {return fTree;}
    virtual Long64_t     GetReadEntry()    const                      {return fTree->GetReadEntry();}
 private:
    AliInputEventHandler(const AliInputEventHandler& handler);             
    AliInputEventHandler& operator=(const AliInputEventHandler& handler);  
 protected:
    TTree        *fTree;    //! Pointer to the tree
    ClassDef(AliInputEventHandler, 1);
};

#endif
