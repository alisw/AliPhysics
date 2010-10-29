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
class AliVCuts;
class AliRunTag;


class AliInputEventHandler : public AliVEventHandler {

 public:
 enum EInputHandlerFlags {
    kUserCallSelectionMask = BIT(14) // Watch out for defining base class bits
 };
    AliInputEventHandler();
    AliInputEventHandler(const char* name, const char* title);
    virtual ~AliInputEventHandler();
    virtual void         SetOutputFileName(const char* /*fname*/) {;}
    virtual const char  *GetOutputFileName()                          {return 0;}
    virtual Bool_t       Init(Option_t* opt) {if(fMixingHandler) fMixingHandler->Init(opt);return kTRUE;}
    virtual Bool_t       Init(TTree* tree, Option_t* opt) {if(fMixingHandler) fMixingHandler->Init(tree,opt);return kTRUE;}
    virtual Bool_t       GetEntry() {if(fMixingHandler) fMixingHandler->GetEntry(); return kTRUE;}
    virtual Bool_t       BeginEvent(Long64_t entry) {if(fMixingHandler) fMixingHandler->BeginEvent(entry);return kTRUE;}

    virtual Bool_t       Notify()      { return AliVEventHandler::Notify();}
    virtual Bool_t       Notify(const char *path) {if(fMixingHandler) fMixingHandler->Notify(path);return kTRUE;}
    virtual Bool_t       FinishEvent() {if(fMixingHandler) fMixingHandler->FinishEvent();return kTRUE;}        
    virtual Bool_t       Terminate()   {if(fMixingHandler) fMixingHandler->Terminate();return kTRUE;}
    virtual Bool_t       TerminateIO() {if(fMixingHandler) fMixingHandler->TerminateIO();return kTRUE;}

    // Setters
    virtual void         SetInputTree(TTree* tree)                    {fTree = tree;}
    virtual void         SetEventSelection(AliVCuts* cuts)            {fEventCuts = cuts;}
    virtual void         SetUserCallSelectionMask(Bool_t flag=kTRUE)  {TObject::SetBit(kUserCallSelectionMask,flag);}
    //
    void SetInactiveBranches(const char* branches) {fBranches   = branches;}
    void SetActiveBranches  (const char* branches) {fBranchesOn = branches;}
     // Getters
    virtual AliVEvent   *GetEvent()        const                      {return 0;}
    virtual AliRunTag   *GetRunTag()       const                      {return 0;}
    virtual Option_t    *GetAnalysisType() const                      {return 0;}
    virtual TTree       *GetTree( )        const                      {return fTree;}
    virtual AliVCuts    *GetEventSelection() const                    {return fEventCuts;}
    virtual Long64_t     GetReadEntry()    const;
    virtual Bool_t       IsUserCallSelectionMask() const              {return TObject::TestBit(kUserCallSelectionMask);}
    virtual Bool_t       NewEvent()
	{Bool_t ne = fNewEvent; fNewEvent = kFALSE; return ne;}
    virtual UInt_t       IsEventSelected() 
        {return fIsSelectedResult;}
    // Mixing
    void SetMixingHandler(AliInputEventHandler* mixing) 
    {fMixingHandler = mixing;}
    AliInputEventHandler* MixingHandler()
    {return fMixingHandler;}

 protected:
    void SwitchOffBranches() const;
    void SwitchOnBranches()  const;
 private:
    AliInputEventHandler(const AliInputEventHandler& handler);             
    AliInputEventHandler& operator=(const AliInputEventHandler& handler);  
 protected:
    TTree          *fTree;         //! Pointer to the tree
    TString         fBranches;     //  List of branches to be switched off (separated by space)
    TString         fBranchesOn;   //  List of branches to be switched on  (separated by space)
    Bool_t          fNewEvent;     //  New event flag 
    AliVCuts*       fEventCuts;    //  Cuts on the event level
    UInt_t          fIsSelectedResult; //  Selection result
    AliInputEventHandler* fMixingHandler; // Optionla plugin for mixing 
    ClassDef(AliInputEventHandler, 5);
};

#endif
