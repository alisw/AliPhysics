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
class AliEventTag;
class TMap;
class AliESDfriend;
class AliESDpid;
class AliESDEvent;
class AliPIDResponse;


class AliESDInputHandler : public AliInputEventHandler {

 public:
    AliESDInputHandler();
    AliESDInputHandler(const char* name, const char* title);
    virtual ~AliESDInputHandler();
    virtual Bool_t       Init(Option_t* opt) {return AliInputEventHandler::Init(opt);}
    virtual Bool_t       Init(TTree* tree, Option_t* opt);
    virtual Bool_t       BeginEvent(Long64_t entry);
    virtual Bool_t       Notify() { return AliInputEventHandler::Notify(); };
    virtual Bool_t       Notify(const char* path);
    virtual Bool_t       FinishEvent();
    void                 CheckSelectionMask();
    AliVEvent         *GetEvent()        const {return (AliVEvent*)fEvent;}
    Option_t            *GetAnalysisType() const {return fAnalysisType;}
    Option_t            *GetDataType() const;
    // Tag cut summary analysis
    Int_t                GetNEventAcceptedInFile();
    Int_t                GetNEventRejectedInFile();
    Bool_t               GetCutSummaryForChain(Int_t *aTotal, Int_t *aAccepted, Int_t *aRejected);
    Int_t                GetNFilesEmpty();
    // HLT  analysis
    AliESDEvent         *GetHLTEvent()     const {return fHLTEvent;}
    TTree               *GetHLTTree()      const {return fHLTTree;}    
    void                 SetReadHLT()            {Changed(); fUseHLT = kTRUE;}
    Bool_t               GetReadHLT()      const {return fUseHLT;}
    // Friends&Co
    AliESDfriend        *GetESDfriend()    const {return fFriend;}
    void                 SetReadFriends(Bool_t flag)   {Changed(); fReadFriends = flag;}
    Bool_t               GetReadFriends()  const {return fReadFriends;}
    void                 SetFriendFileName(const char *fname)  {Changed(); fFriendFileName = fname;}
    const char          *GetFriendFileName() const {return fFriendFileName;}
    // Tag analysis
    void                 SetReadTags()           {Changed(); fUseTags = kTRUE;}
    AliRunTag           *GetRunTag() const       {return fRunTag;}
    const AliEventTag   *GetEventTag() const     {return fEventTag;}
    // Get the statistics object (currently TH2). Option can be BIN0.
    virtual TObject     *GetStatistics(Option_t *option="") const;

    //PID response
    virtual AliPIDResponse* GetPIDResponse() {return (AliPIDResponse*)fESDpid;}
    virtual void CreatePIDResponse(Bool_t isMC=kFALSE);
    AliESDpid           *GetESDpid()       const {return fESDpid;}
    void                 SetESDpid(AliESDpid* pid)     {Changed(); fESDpid = pid;}

    //HLT
    virtual AliVfriendEvent*   GetVfriendEvent() const {return fFriend;};
  
 private:
    AliESDInputHandler(const AliESDInputHandler& handler);             
    AliESDInputHandler& operator=(const AliESDInputHandler& handler);  
    void                 ConnectFriends();
 protected:
    // ESD event
    AliESDEvent    *fEvent;         //! Pointer to the event
    AliESDfriend   *fFriend;        //! Pointer to the esd friend
    AliESDpid      *fESDpid;        //! Pointer to PID information
    Option_t       *fAnalysisType;  //! local, proof, grid
    Int_t           fNEvents;       //! Number of events in the current tree
    // HLT event
    AliESDEvent    *fHLTEvent;      //! Pointer to the HLT Event (if present)
    TTree          *fHLTTree;       //! Pointer to the HLT Event (if present)
    Bool_t          fUseHLT;        //  Flag to access HLT Events
    // ESD Tag Cut Summary
    TMap           *fTagCutSumm;    //! Tag cut summary map
    // ESD Tags (optional)
    Bool_t          fUseTags;       //  Flag to use tags
    TChain         *fChainT;        //! File with event tags
    TTree          *fTreeT;         //! Tree of tags
    AliRunTag      *fRunTag;        //! Pointer to the run tag
    const AliEventTag *fEventTag;      //! Current event tag
    // Friends
    Bool_t          fReadFriends;   //  Flag for friends reading 
    TString         fFriendFileName;//  Name of the file containing the frien tree (branch)
    ClassDef(AliESDInputHandler, 6);
};

#endif
