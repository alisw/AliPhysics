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
class TObjArray;
class AliMCEventHandler;
class AliAODMCHeader;
class AliGenEventHeader;
class TMap;
class AliAODExtension;
class AliAnalysisFilter;

class AliAODHandler : public AliVEventHandler {
    
 public:
    AliAODHandler();
    AliAODHandler(const char* name, const char* title);
    virtual ~AliAODHandler();
    virtual void         SetOutputFileName(const char* fname);
    virtual const char*  GetOutputFileName();
    // Extra outputs as a string separated by commas
    virtual const char*  GetExtraOutputs() const;
    virtual Bool_t       Init(Option_t* option);
    virtual Bool_t       Init(TTree* /*tree*/, Option_t* /*option*/)  {return kTRUE;}
    virtual Bool_t       GetEntry() {return kTRUE;}
	    
    virtual Bool_t       BeginEvent(Long64_t /*entry*/) {fFillAOD=kFALSE; fFillExtension=kFALSE; return kTRUE;}
    virtual Bool_t       Notify() { return AliVEventHandler::Notify(); };
    virtual Bool_t       Notify(const char * /* path */) {return kTRUE;}
    virtual Bool_t       FinishEvent();
    virtual Bool_t       Terminate();
    virtual Bool_t       TerminateIO();
    //
    virtual void         SetCreateNonStandardAOD()               {fIsStandard = kFALSE;}
    virtual void         SetFillAOD(Bool_t b)                    {fFillAOD = b;}
    virtual void         SetFillExtension(Bool_t b)              {fFillExtension = b;}
    virtual void         SetFillAODforRun(Bool_t b)              {fFillAODRun = b;}
    virtual void         SetNeedsHeaderReplication()             {fNeedsHeaderReplication             = kTRUE;}
    virtual void         SetNeedsTracksBranchReplication()       {fNeedsTracksBranchReplication       = kTRUE;}
    virtual void         SetNeedsVerticesBranchReplication()     {fNeedsVerticesBranchReplication     = kTRUE;}
    virtual void         SetNeedsV0sBranchReplication()          {fNeedsV0sBranchReplication          = kTRUE;}
    virtual void         SetNeedsCascadesBranchReplication()     {fNeedsCascadesBranchReplication     = kTRUE;}
    virtual void         SetNeedsTrackletsBranchReplication()    {fNeedsTrackletsBranchReplication    = kTRUE;}
    virtual void         SetNeedsPMDClustersBranchReplication()  {fNeedsPMDClustersBranchReplication  = kTRUE;}
    virtual void         SetNeedsJetsBranchReplication()         {fNeedsJetsBranchReplication         = kTRUE;}
    virtual void         SetNeedsFMDClustersBranchReplication()  {fNeedsFMDClustersBranchReplication  = kTRUE;}
    virtual void         SetNeedsCaloClustersBranchReplication() {fNeedsCaloClustersBranchReplication = kTRUE;}
    virtual void         SetNeedsMCParticlesBranchReplication()  {fNeedsMCParticlesBranchReplication  = kTRUE;}
    virtual void         SetNeedsDimuonsBranchReplication()      {fNeedsDimuonsBranchReplication      = kTRUE;}
    virtual void         SetAODIsReplicated() {fAODIsReplicated = kTRUE;}
    virtual void         SetAODExtensionMode()                   {fIsStandard=kFALSE; fFillAOD=kFALSE; fFillAODRun=kFALSE; fFillExtension=kTRUE;}
    //
    AliAODEvent*         GetAOD()  {return fAODEvent;}
    virtual TTree*       GetTree() const {return fTreeA;}
    TObjArray*           GetExtensions() const {return fExtensions;}
    AliAODExtension*     GetExtension(const char *filename) const;
    TObjArray*           GetFilters() const {return fFilters;}
    AliAODExtension*     GetFilteredAOD(const char *filename) const;
    void                 CreateTree(Int_t flag);
    void                 FillTree();
    void                 AddAODtoTreeUserInfo();
    void                 AddBranch(const char* cname, void* addobj, const char *fname="");
  
    AliAODExtension*     AddExtension(const char *filename, const char *title="");                 
    AliAODExtension*     AddFilteredAOD(const char *filename, const char *filtername);
//    AliAODExtension*     FindExtensionContainingBranch(const char* bname) const;
    Bool_t               IsStandard()                         const {return fIsStandard;}
    Bool_t               GetFillAOD()                         const {return fFillAOD;} 
    Bool_t               NeedsHeaderReplication()             const {return  fNeedsHeaderReplication;}
    Bool_t               NeedsTracksBranchReplication()       const {return  fNeedsTracksBranchReplication;}
    Bool_t               NeedsVerticesBranchReplication()     const {return  fNeedsVerticesBranchReplication;}
    Bool_t               NeedsV0sBranchReplication()          const {return  fNeedsV0sBranchReplication;}
    Bool_t               NeedsCascadesBranchReplication()     const {return  fNeedsCascadesBranchReplication;}
    Bool_t               NeedsTrackletsBranchReplication()    const {return  fNeedsTrackletsBranchReplication;}
    Bool_t               NeedsPMDClustersBranchReplication()  const {return  fNeedsPMDClustersBranchReplication;}
    Bool_t               NeedsJetsBranchReplication()         const {return  fNeedsJetsBranchReplication;}
    Bool_t               NeedsFMDClustersBranchReplication()  const {return  fNeedsFMDClustersBranchReplication;}
    Bool_t               NeedsCaloClustersBranchReplication() const {return  fNeedsCaloClustersBranchReplication;}
    Bool_t               NeedsMCParticlesBranchReplication()  const {return  fNeedsMCParticlesBranchReplication;}
    Bool_t               NeedsDimuonsBranchReplication()      const {return  fNeedsDimuonsBranchReplication;}
    Bool_t               AODIsReplicated()                    const {return  fAODIsReplicated;}
    //
    void                 SetInputTree(TTree* /*tree*/) {;}
    void                 SetMCEventHandler(AliMCEventHandler* mcH) {fMCEventH = mcH;} // For internal use
    void StoreMCParticles(); // Store MC particles, only to be called from AliAnalyisTaskMCParticleFilter

  Bool_t HasExtensions() const;
  
  void Print(Option_t* opt="") const;
  
 private:
  void SetMCHeaderInfo(AliAODMCHeader *mcHeader,AliGenEventHeader *genHeader); // Utility function t catch different types of eventheaders
  AliAODHandler(const AliAODHandler&);             // Not implemented
  AliAODHandler& operator=(const AliAODHandler&);  // Not implemented
  void PrintExtensions(const TObjArray& array) const;

 private:
    Bool_t                   fIsStandard;                         // Flag for standard aod creation
    Bool_t                   fFillAOD;                            // Flag for filling of the AOD tree at the end (all or nothing evt by evt)
    Bool_t                   fFillAODRun;                         // Flag for filling of the AOD tree at the end (run)
    Bool_t                   fFillExtension;                      // Flag for filling or the delta AOD tree at the end
    Bool_t                   fNeedsHeaderReplication;             // Flag for header replication
    Bool_t                   fNeedsTracksBranchReplication;       // Flag for tracks replication
    Bool_t                   fNeedsVerticesBranchReplication;     // Flag for vertices replication
    Bool_t                   fNeedsV0sBranchReplication;          // Flag for V0s replication
    Bool_t                   fNeedsCascadesBranchReplication;     // Flag for Cascade replication
    Bool_t                   fNeedsTrackletsBranchReplication;    // Flag for Tracklets replication
    Bool_t                   fNeedsPMDClustersBranchReplication;  // Flag for PMDClusters replication
    Bool_t                   fNeedsJetsBranchReplication;         // Flag for Jets replication
    Bool_t                   fNeedsFMDClustersBranchReplication;  // Flag for FMDClusters replication
    Bool_t                   fNeedsCaloClustersBranchReplication; // Flag for CaloClusters replication
    Bool_t                   fNeedsMCParticlesBranchReplication;  // Flag for MCParticles replication
    Bool_t                   fNeedsDimuonsBranchReplication;      // Flag for Dimuons replication
    Bool_t                   fAODIsReplicated;                    // Flag true if replication as been executed
    AliAODEvent             *fAODEvent;               //! Pointer to the AOD event
    AliMCEventHandler       *fMCEventH;               //! Pointer to mc event handler needed not to depend on the manager
    TTree                   *fTreeA;                  //! tree for AOD persistency
    TFile                   *fFileA;                  //! Output file
    TString                  fFileName;               //  Output file name
    TObjArray               *fExtensions;             //  List of extensions
    TObjArray               *fFilters;                //  List of filtered AOD's

  ClassDef(AliAODHandler, 7)
};

#endif
