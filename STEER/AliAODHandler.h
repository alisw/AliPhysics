#ifndef ALIAODHANDLER_H
#define ALIAODHANDLER_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     Implementation of the Event Handler Interface for AOD
//     Author: Andreas Morsch, CERN
//-------------------------------------------------------------------------

#include "AliVEventHandler.h"

class AliAODEvent;
class TFile;
class TTree;
class AliMCEventHandler;
class AliAODMCHeader;
class AliGenEventHeader;



class AliAODHandler : public AliVEventHandler {
    
 public:
    AliAODHandler();
    AliAODHandler(const char* name, const char* title);
    virtual ~AliAODHandler();
    virtual void         SetOutputFileName(const char* fname);
    virtual const char*  GetOutputFileName();
    virtual Bool_t       Init(Option_t* option);
    virtual Bool_t       Init(TTree* /*tree*/, Option_t* /*option*/)  {return kTRUE;}
    virtual Bool_t       BeginEvent(Long64_t /*entry*/){return kTRUE;}
    virtual Bool_t       Notify() { return AliVEventHandler::Notify(); };
    virtual Bool_t       Notify(const char * /* path */) {return kTRUE;}
    virtual Bool_t       FinishEvent();
    virtual Bool_t       Terminate();
    virtual Bool_t       TerminateIO();
    //
    virtual void         SetCreateNonStandardAOD()   {fIsStandard = kFALSE;}
    virtual void         SetFillAOD(Bool_t b)      {fFillAOD = b;}
    virtual void         SetNeedsHeaderReplication() {fNeedsHeaderReplication = kTRUE;}
    virtual void         SetNeedsTracksBranchReplication() {fNeedsTracksBranchReplication = kTRUE;}
    virtual void         SetNeedsVerticesBranchReplication() {fNeedsVerticesBranchReplication = kTRUE;}
    virtual void         SetNeedsV0sBranchReplication() {fNeedsV0sBranchReplication = kTRUE;}
    virtual void         SetNeedsTrackletsBranchReplication() {fNeedsTrackletsBranchReplication = kTRUE;}
    virtual void         SetNeedsPMDClustersBranchReplication() {fNeedsPMDClustersBranchReplication = kTRUE;}
    virtual void         SetNeedsJetsBranchReplication() {fNeedsJetsBranchReplication = kTRUE;}
    virtual void         SetNeedsFMDClustersBranchReplication() {fNeedsFMDClustersBranchReplication = kTRUE;}
    virtual void         SetNeedsCaloClustersBranchReplication() {fNeedsCaloClustersBranchReplication = kTRUE;}
    virtual void         SetAODIsReplicated() {fAODIsReplicated = kTRUE;}
    //
    AliAODEvent*         GetAOD()  {return fAODEvent;}
    virtual TTree*       GetTree() const {return fTreeA;}
    void                 CreateTree(Int_t flag);
    void                 FillTree();
    void                 AddAODtoTreeUserInfo();
    void                 AddBranch(const char* cname, void* addobj);
    Bool_t               IsStandard() {return fIsStandard;}
    Bool_t               GetFillAOD(){return fFillAOD;} 
    Bool_t               NeedsHeaderReplication() {return  fNeedsHeaderReplication;}
    Bool_t               NeedsTracksBranchReplication() {return  fNeedsTracksBranchReplication;}
    Bool_t               NeedsVerticesBranchReplication() {return  fNeedsVerticesBranchReplication;}
    Bool_t               NeedsV0sBranchReplication() {return  fNeedsV0sBranchReplication;}
    Bool_t               NeedsTrackletsBranchReplication() {return  fNeedsTrackletsBranchReplication;}
    Bool_t               NeedsPMDClustersBranchReplication() {return  fNeedsPMDClustersBranchReplication;}
    Bool_t               NeedsJetsBranchReplication() {return  fNeedsJetsBranchReplication;}
    Bool_t               NeedsFMDClustersBranchReplication() {return  fNeedsFMDClustersBranchReplication;}
    Bool_t               NeedsCaloClustersBranchReplication() {return  fNeedsCaloClustersBranchReplication;}
    Bool_t               AODIsReplicated() {return fAODIsReplicated;}
    //
    void                 SetInputTree(TTree* /*tree*/) {;}
    void                 SetMCEventHandler(AliMCEventHandler* mcH) {fMCEventH = mcH;} // For internal use
 private:
    void StoreMCParticles();
    void SetMCHeaderInfo(AliAODMCHeader *mcHeader,AliGenEventHeader *genHeader); // Utility function t catch different types of eventheaders
    AliAODHandler(const AliAODHandler&);             // Not implemented
    AliAODHandler& operator=(const AliAODHandler&);  // Not implemented
 private:
    Bool_t                   fIsStandard;                         // Flag for standard aod creation
    Bool_t                   fFillAOD;                          // Flag for filling of the AOD tree at the end (all or nothing)
    Bool_t                   fNeedsHeaderReplication;             // Flag for header replication
    Bool_t                   fNeedsTracksBranchReplication;       // Flag for tracks replication
    Bool_t                   fNeedsVerticesBranchReplication;     // Flag for vertices replication
    Bool_t                   fNeedsV0sBranchReplication;          // Flag for V0s replication
    Bool_t                   fNeedsTrackletsBranchReplication;    // Flag for Tracklets replication
    Bool_t                   fNeedsPMDClustersBranchReplication;  // Flag for PMDClusters replication
    Bool_t                   fNeedsJetsBranchReplication;         // Flag for Jets replication
    Bool_t                   fNeedsFMDClustersBranchReplication;  // Flag for FMDClusters replication
    Bool_t                   fNeedsCaloClustersBranchReplication; // Flag for CaloClusters replication
    Bool_t                   fAODIsReplicated;                    // Flag true if replication as been executed
    AliAODEvent             *fAODEvent;               //! Pointer to the AOD event
    AliMCEventHandler       *fMCEventH;               //! Pointer to mc event handler needed not to depend on the manager
    TTree                   *fTreeA;                  //! tree for AOD persistency
    TFile                   *fFileA;                  //! Output file
    TString                  fFileName;               //  Output file name
    ClassDef(AliAODHandler, 3)
};

#endif
