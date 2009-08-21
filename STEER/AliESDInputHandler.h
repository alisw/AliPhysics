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
class TChain;
class TTree;
class AliRunTag;


class AliESDInputHandler : public AliInputEventHandler {

 public:
    AliESDInputHandler();
    AliESDInputHandler(const char* name, const char* title);
    virtual ~AliESDInputHandler();
    virtual Bool_t       Init(Option_t* /*opt*/) {return kTRUE;}
    virtual Bool_t       Init(TTree* tree, Option_t* opt);
    virtual Bool_t       BeginEvent(Long64_t entry);
    virtual Bool_t       Notify() { return AliVEventHandler::Notify(); };
    virtual Bool_t       Notify(const char* path);
    virtual Bool_t       FinishEvent();
    AliESDEvent         *GetEvent()        const {return fEvent;}
    Option_t            *GetAnalysisType() const {return fAnalysisType;}
    Option_t            *GetDataType() const;
    // HLT analysis
    AliESDEvent         *GetHLTEvent()     const {return fHLTEvent;}
    TTree               *GetHLTTree()      const {return fHLTTree;}    
    void                SetReadHLT()             {fUseHLT = kTRUE;}
    // Tag analysis
    void SetReadTags() {fUseTags = kTRUE;}
    AliRunTag           *GetRunTag() const {return fRunTag;}
	    
 private:
    AliESDInputHandler(const AliESDInputHandler& handler);             
    AliESDInputHandler& operator=(const AliESDInputHandler& handler);  
 protected:
    // ESD event
    AliESDEvent    *fEvent;        //! Pointer to the event
    Option_t       *fAnalysisType; //! local, proof, grid
    Int_t           fNEvents;      //! Number of events in the current tree
    // HLT event
    AliESDEvent    *fHLTEvent;     //! Pointer to the HLT Event (if present)
    TTree          *fHLTTree;      //! Pointer to the HLT Event (if present)
    Bool_t          fUseHLT;       //  Flag to access HLT Events
    // ESD Tags (optional)
    Bool_t          fUseTags;    //  Flag to use tags
    TChain         *fChainT;     //! File with event tags
    TTree          *fTreeT;      //! Tree of tags
    AliRunTag      *fRunTag;     //! Pointer to the run tag
    ClassDef(AliESDInputHandler, 4);
};

#endif
