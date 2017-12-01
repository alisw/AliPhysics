#ifndef ALIVEVENTHANDLER_H
#define ALIVEVENTHANDLER_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     Event Handler base class
//     Author: Andreas Morsch, CERN
//-------------------------------------------------------------------------

#include <TNamed.h>

class TTree;
class TObjArray;
class AliVEvent;
class AliVfriendEvent;
class AliVCuts;
class AliEventTag;
class AliRunTag;
class AliPIDResponse;
class AliMCEvent;

class AliVEventHandler : public TNamed {

 public:
enum EEventHandlerFlags {
   kHandlerLocked       = BIT(14)
};
    AliVEventHandler();
    AliVEventHandler(const char* name, const char* title);
    virtual ~AliVEventHandler();
    // Handled tree
    virtual TTree       *GetTree() const { return NULL; }
    virtual Option_t    *GetDataType() const { return NULL; }
    virtual Bool_t       GetFillAOD() const {return kTRUE;}
    virtual Bool_t       GetFillExtension() const {return kTRUE;}
    virtual void         SetFillAOD(Bool_t) {}
    virtual void         SetFillExtension(Bool_t) {}
    // Input
    virtual void         SetInputFileName(const char*) {}
    virtual const char*  GetInputFileName() const {return 0;}
    // Output
    virtual void         SetOutputFileName(const char* fname)   = 0;
    virtual const char*  GetOutputFileName() const        = 0;
    // Extra outputs as a string separated by commas
    virtual const char*  GetExtraOutputs(Bool_t merge=kFALSE) const;
    // Input
    virtual void         SetInputTree(TTree* tree)        = 0;
    // Steering 
    virtual Bool_t       Init(Option_t* opt)              = 0;
    virtual Bool_t       Init(TTree* tree, Option_t* opt) = 0;
    virtual Bool_t       BeginEvent(Long64_t entry)       = 0;
    virtual Bool_t       GetEntry()                       = 0;
    virtual Bool_t       Notify(const char *path)         = 0;
    virtual Bool_t       FinishEvent()                    = 0;
    virtual Bool_t       Terminate()                      = 0;
    virtual Bool_t       TerminateIO()                    = 0;
    virtual AliVCuts*    GetEventSelection() const        = 0;
    virtual UInt_t       IsEventSelected() {return 0;}
    virtual Long64_t     GetReadEntry() const {return 0;}
    virtual const AliEventTag* GetEventTag() const {return 0;}
    virtual Option_t    *GetAnalysisType() const                      {return 0;}
    virtual AliRunTag   *GetRunTag()       const                      {return 0;}
    virtual AliPIDResponse* GetPIDResponse() {return 0x0;}
    virtual AliMCEvent*  MCEvent() const			      {return 0;}
    virtual void         SetNeedField(Bool_t flag=kTRUE)  {}
    virtual TObject     *GetStatistics(Option_t *option="") const {return NULL;}
    //
    virtual Bool_t       Notify() { return TNamed::Notify(); };
    // Security
    Bool_t               IsLocked() const {return TObject::TestBit(kHandlerLocked);}
    void                 Lock();
    void                 UnLock();
    void                 Changed();
    virtual void         SetCacheSize(Long64_t) {}
    virtual TList        *GetUserInfo() const {return 0x0;};

    // HLT
    virtual Bool_t              InitTaskInputData(AliVEvent* /*event*/, AliVfriendEvent* /*esdFriend*/, TObjArray* /*arrTasks*/) {return kTRUE;};
    virtual AliVEvent*          GetEvent() const {return 0x0;};
    virtual AliVfriendEvent*   GetVfriendEvent() const {return 0x0;};

  ClassDef(AliVEventHandler, 1);
};

#endif
