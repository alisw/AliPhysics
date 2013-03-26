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
class AliEventTag;
class AliPIDResponse;

class AliInputEventHandler : public AliVEventHandler {

 public:
 enum EInputHandlerFlags {
    kUserCallSelectionMask = BIT(14), // Watch out for defining base class bits
    kCheckStatistics       = BIT(15)
 };
    AliInputEventHandler();
    AliInputEventHandler(const char* name, const char* title);
    virtual ~AliInputEventHandler();
    virtual void         SetInputFileName(const char* fname);
    virtual const char  *GetInputFileName() const                     {return fInputFileName;}
    virtual void         SetOutputFileName(const char* /*fname*/) {;}
    virtual const char  *GetOutputFileName() const                    {return 0;}
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
    virtual void         SetEventSelection(AliVCuts* cuts)            {if (fEventCuts) Changed(); fEventCuts = cuts;}
    virtual void         SetUserCallSelectionMask(Bool_t flag=kTRUE)  {TObject::SetBit(kUserCallSelectionMask,flag);}
    virtual void         SetCheckStatistics(Bool_t flag=kTRUE)        {Changed(); TObject::SetBit(kCheckStatistics,flag);}
    //
    void SetInactiveBranches(const char* branches) {Changed(); fBranches   = branches;}
    void SetActiveBranches  (const char* branches) {Changed(); fBranchesOn = branches;}
     // Getters
    virtual AliVEvent   *GetEvent()        const                      {return 0;}
    virtual const AliEventTag   *GetEventTag() const                  {return 0;}
    virtual AliRunTag   *GetRunTag()       const                      {return 0;}
    // Get the statistics object (currently TH2). Option can be BIN0.
    virtual TObject     *GetStatistics(Option_t *option="") const;
    virtual Option_t    *GetAnalysisType() const                      {return 0;}
    virtual TTree       *GetTree( )        const                      {return fTree;}
    virtual AliVCuts    *GetEventSelection() const                    {return fEventCuts;}
    virtual Long64_t     GetReadEntry()    const;
    virtual Bool_t       IsCheckStatistics() const                    {return TObject::TestBit(kCheckStatistics);}
    virtual Bool_t       IsUserCallSelectionMask() const              {return TObject::TestBit(kUserCallSelectionMask);}
    virtual Bool_t       NewEvent()
	{Bool_t ne = fNewEvent; fNewEvent = kFALSE; return ne;}
    virtual UInt_t       IsEventSelected() 
        {return fIsSelectedResult;}
    TList       *GetUserInfo() const                         {return fUserInfo;}
    // Mixing
    void SetMixingHandler(AliInputEventHandler* mixing) {Changed(); fMixingHandler = mixing;}
    AliInputEventHandler* MixingHandler()               {return fMixingHandler;}
    // Parent Handler
    void SetParentHandler(AliInputEventHandler* parent) {Changed(); fParentHandler = parent;}
    AliInputEventHandler* ParentHandler()               {return fParentHandler;}

    //PID response
    virtual AliPIDResponse* GetPIDResponse() {return 0x0;}
    virtual void CreatePIDResponse(Bool_t /*isMC*/=kFALSE) {;}
  
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
    TString         fInputFileName; // Name of the input file
    Bool_t          fNewEvent;     //  New event flag 
    AliVCuts*       fEventCuts;    //  Cuts on the event level
    UInt_t          fIsSelectedResult; //  Selection result
    AliInputEventHandler* fMixingHandler; // Optionla plugin for mixing
    AliInputEventHandler* fParentHandler; // optional pointer to parent handlers (used in AliMultiInputEventHandler)
    TList           *fUserInfo;     //! transient user info for current tree
    ClassDef(AliInputEventHandler, 7);
};

#endif
