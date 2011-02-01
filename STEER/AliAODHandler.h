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
class AliAODExtension;
class AliGenEventHeader;
class TMap;
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
	    
    virtual Bool_t       BeginEvent(Long64_t /*entry*/) {fFillAOD=kFALSE; return kTRUE;}
    virtual Bool_t       Notify() { return AliVEventHandler::Notify(); };
    virtual Bool_t       Notify(const char * /* path */) {return kTRUE;}
    virtual Bool_t       FinishEvent();
    virtual Bool_t       Terminate();
    virtual Bool_t       TerminateIO();
    //
    virtual void         SetCreateNonStandardAOD()               {fIsStandard = kFALSE;}
    virtual void         SetFillAOD(Bool_t b)                    {fFillAOD = b;}
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
    ClassDef(AliAODHandler, 6)
};

//-------------------------------------------------------------------------
//     Support class for AOD extensions. This is created by the user analysis
//     that requires a separate file for some AOD branches. The name of the 
//     AliAODExtension object is the file name where the AOD branches will be
//     stored.
//     Author: Andrei Gheata, CERN
//-------------------------------------------------------------------------

class AliAODBranchReplicator;

class AliAODExtension : public TNamed {

public:

enum EAliAODExtensionFlags {
   kFilteredAOD      = BIT(14)
};
    
  AliAODExtension();
    AliAODExtension(const char* name, const char* title, Bool_t isfilter=kFALSE);
    virtual ~AliAODExtension();
    void                 AddBranch(const char* cname, void* addobj);
    Bool_t               FinishEvent();
    Int_t                GetNtotal() const         {return fNtotal;}
    Int_t                GetNpassed() const        {return fNpassed;}
    const char*          GetOutputFileName() const {return TNamed::GetName();}
    AliAODEvent*         GetAOD() const            {return fAODEvent;}
    TTree*               GetTree() const           {return fTreeE;}
    Bool_t               Init(Option_t *option);
    Bool_t               IsFilteredAOD() const     {return TObject::TestBit(kFilteredAOD);}
    Bool_t               IsEventSelected() const   {return fSelected;}
    void                 SelectEvent(Bool_t flag=kTRUE)  {fSelected = flag;}
    void                 SetEvent(AliAODEvent *event);
    void                 SetOutputFileName(const char* fname) {TNamed::SetName(fname);}
    Bool_t               TerminateIO();

  void Print(Option_t* opt="") const;
  
  void FilterBranch(const char* branchName, AliAODBranchReplicator* replicator=0x0);

  /* Use DisableReferences if and only if the output AOD contains no TRef or TRefArray,
   otherwise the produced AOD won't be valid.
   */
  void DisableReferences() { fEnableReferences=kFALSE; }
  
  void EnableReferences() { fEnableReferences=kTRUE; }
  
 private:
    AliAODExtension(const AliAODExtension&);             // Not implemented
    AliAODExtension& operator=(const AliAODExtension&);  // Not implemented

 private:
    AliAODEvent             *fAODEvent;               //! Pointer to the AOD event
    TTree                   *fTreeE;                  //! tree for AOD persistency
    TFile                   *fFileE;                  //! Output file
    Int_t                    fNtotal;                 //! Number of processed events
    Int_t                    fNpassed;                //! Number of events that passed the filter
    Bool_t                   fSelected;               //! Select current event for filtered AOD's. Made false at event start.

    TMap*                    fRepFiMap; // which branch(es) to filter out / and or replicate
    TList*                   fRepFiList; // list of unique filter/replicator
  
    Bool_t                   fEnableReferences; // whether or not to enable the TRefTable branch
    TList*                   fObjectList; //! internal list of which objects to keep 
  
  ClassDef(AliAODExtension, 2) // Support for extra AOD branches in a separate AOD file
};
#endif
