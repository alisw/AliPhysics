#ifndef ALIRECOINPUTHANDLER_H
#define ALIRECOINPUTHANDLER_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//     Reconstruction-specific input handler
//     Author: Andrei Gheata, CERN
//-------------------------------------------------------------------------

#ifndef ALIESDINPUTHANDLER_H
#include "AliESDInputHandler.h"
#endif

class AliReconstruction;

class AliRecoInputHandler : public AliESDInputHandler {

 public:
    AliRecoInputHandler() {}
    AliRecoInputHandler(const char* name, const char* title);
    virtual ~AliRecoInputHandler() {}
    virtual Bool_t       Notify() { return AliESDInputHandler::Notify(); };
    virtual Bool_t       Notify(const char *) {return kTRUE;}
    virtual Bool_t       Init(Option_t* opt) {return AliESDInputHandler::Init(opt);}
    virtual Bool_t       Init(TTree* tree, Option_t* opt="LOCAL");
    virtual Bool_t       BeginEvent(Long64_t entry);
    virtual Bool_t       FinishEvent() {return kTRUE;}
//    void                 CheckSelectionMask();
//    AliESDEvent         *GetEvent()        const {return fEvent;}
//    Option_t            *GetAnalysisType() const {return fAnalysisType;}
//    Option_t            *GetDataType() const;
    // Tag cut summary analysis
//    Int_t                GetNEventAcceptedInFile();
//    Int_t                GetNEventRejectedInFile();
//    Bool_t               GetCutSummaryForChain(Int_t *aTotal, Int_t *aAccepted, Int_t *aRejected);
//    Int_t                GetNFilesEmpty();
    // HLT  analysis
//    AliESDEvent         *GetHLTEvent()     const {return fHLTEvent;}
//    TTree               *GetHLTTree()      const {return fHLTTree;}    
//    void                 SetReadHLT()            {fUseHLT = kTRUE;}
    // Friends&Co
//    AliESDfriend        *GetESDfriend()    const {return fFriend;}
//    void                 SetReadFriends(Bool_t flag)   {fReadFriends = flag;}
//    void                 SetFriendFileName(const char *fname)  {fFriendFileName = fname;}
    // Tag analysis
//    void                 SetReadTags()           {fUseTags = kTRUE;}
//    AliRunTag           *GetRunTag() const       {return fRunTag;}
//    const AliEventTag   *GetEventTag() const     {return fEventTag;}
    // Get the statistics object (currently TH2). Option can be BIN0.
//    virtual TObject     *GetStatistics(Option_t *option="") const;

    //PID response
//    virtual AliPIDResponse* GetPIDResponse() {return (AliPIDResponse*)fESDpid;}
//    virtual void CreatePIDResponse(Bool_t isMC=kFALSE);
//    AliESDpid           *GetESDpid()       const {return fESDpid;}
//    void                 SetESDpid(AliESDpid* pid)     {fESDpid = pid;}
      
 private:
    AliRecoInputHandler(const AliESDInputHandler& handler);             
    AliRecoInputHandler& operator=(const AliESDInputHandler& handler);  
    // Private setters used by AliReconstruction
    friend class AliReconstruction;
    void                 SetEvent(AliESDEvent *event)          {fEvent = event;}
    void                 SetESDfriend(AliESDfriend *esdfriend) {fFriend = esdfriend;}
    void                 SetHLTEvent(AliESDEvent *hltevent)    {fHLTEvent = hltevent;}
    void                 SetHLTTree(TTree *hlttree)            {fHLTTree = hlttree;}
    
    ClassDef(AliRecoInputHandler, 1);
};

#endif
