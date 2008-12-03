#ifndef ALIANALYSISMANAGER_H
#define ALIANALYSISMANAGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Author: Andrei Gheata, 31/05/2006

//==============================================================================
//   AliAnalysysManager - Manager analysis class. Allows creation of several
// analysis tasks and data containers storing their input/output. Allows 
// connecting/chaining tasks via shared data containers. Serializes the current
// event for all tasks depending only on initial input data.
//==============================================================================

#ifndef ROOT_TNamed
#include <TNamed.h>
#endif

class TClass;
class TTree;
class TFile;
class AliAnalysisSelector;
class AliAnalysisDataContainer;
class AliAnalysisTask;
class AliVEventHandler;
class AliVEventPool;
class AliAnalysisGrid;


class AliAnalysisManager : public TNamed {

public:

enum EAliAnalysisContType {
   kExchangeContainer  = 0,
   kInputContainer   = 1,
   kOutputContainer  = 2
};   

enum EAliAnalysisExecMode {
   kLocalAnalysis    = 0,
   kProofAnalysis    = 1,
   kGridAnalysis     = 2,
   kMixingAnalysis   = 3
};

enum EAliAnalysisFlags {
   kEventLoop        = BIT(14),
   kDisableBranches  = BIT(15),
   kUseDataSet       = BIT(16)
};   

   AliAnalysisManager(const char *name = "mgr", const char *title="");
   virtual            ~AliAnalysisManager();

   AliAnalysisManager(const AliAnalysisManager& other);
   AliAnalysisManager& operator=(const AliAnalysisManager& other);
   
   // Management
   void                StartAnalysis(const char *type="local", TTree *tree=0, Long64_t nentries=1234567890, Long64_t firstentry=0);
   void                StartAnalysis(const char *type, const char *dataset, Long64_t nentries=1234567890, Long64_t firstentry=0);

   virtual Bool_t      Init(TTree *tree);   
   virtual Bool_t      Notify();
   virtual void        SlaveBegin(TTree *tree);
   virtual Bool_t      ProcessCut(Long64_t entry) {return Process(entry);}
   virtual Bool_t      Process(Long64_t entry);
   virtual Int_t       GetEntry(Long64_t entry, Int_t getall = 0);
   TFile              *OpenProofFile(const char *name, const char *option);
   void                PackOutput(TList *target);
   void                UnpackOutput(TList *source);
   virtual void        Terminate();

   // Getters/Setters
   static AliAnalysisManager *GetAnalysisManager() {return fgAnalysisManager;}
   TObjArray          *GetContainers() const {return fContainers;}
   UInt_t              GetDebugLevel() const {return fDebug;}
   TObjArray          *GetInputs() const     {return fInputs;}
   TObjArray          *GetOutputs() const    {return fOutputs;}
   TObjArray          *GetTasks() const      {return fTasks;}
   TObjArray          *GetTopTasks() const   {return fTopTasks;}
   TTree              *GetTree() const       {return fTree;}
   TObjArray          *GetZombieTasks() const {return fZombies;}
   Long64_t            GetCurrentEntry() const {return fCurrentEntry;}
   EAliAnalysisExecMode 
                       GetAnalysisType() const {return fMode;}
   Bool_t              IsUsingDataSet() const  {return TObject::TestBit(kUseDataSet);}

   void                SetAnalysisType(EAliAnalysisExecMode mode) {fMode = mode;}
   void                SetCurrentEntry(Long64_t entry) {fCurrentEntry = entry;}
   void                SetDebugLevel(UInt_t level) {fDebug = level;}
   void                SetSpecialOutputLocation(const char *location) {fSpecialOutputLocation = location;}
   void                SetDisableBranches(Bool_t disable=kTRUE) {TObject::SetBit(kDisableBranches,disable);}
   void                SetCollectSysInfoEach(Int_t nevents=0) {fNSysInfo = nevents;}
   void                SetInputEventHandler(AliVEventHandler*  handler)  {fInputEventHandler   = handler;}
   void                SetOutputEventHandler(AliVEventHandler*  handler) {fOutputEventHandler  = handler;}
   void                SetMCtruthEventHandler(AliVEventHandler* handler) {fMCtruthEventHandler = handler;}
   void                SetGridHandler(AliAnalysisGrid *handler) {fGridHandler = handler;}
   void                SetEventPool(AliVEventPool* epool) {fEventPool = epool;}
   void                SetNSysInfo(Long64_t nevents) {fNSysInfo = nevents;}
   void                SetSelector(AliAnalysisSelector *sel) {fSelector = sel;}
   AliVEventHandler*   GetInputEventHandler()   {return fInputEventHandler;}
   AliVEventHandler*   GetOutputEventHandler()  {return fOutputEventHandler;}
   AliVEventHandler*   GetMCtruthEventHandler() {return fMCtruthEventHandler;}
   AliAnalysisGrid*    GetGridHandler()         {return fGridHandler;}
   AliVEventPool*      GetEventPool()           {return fEventPool;}

   // Container handling
   AliAnalysisDataContainer *CreateContainer(const char *name, TClass *datatype, 
                       EAliAnalysisContType type     = kExchangeContainer, 
                       const char          *filename = NULL);
   
   // Including tasks and getting them
   void                 AddTask(AliAnalysisTask *task);
   AliAnalysisTask     *GetTask(const char *name) const;
   
   // Connecting data containers to task inputs/outputs
   Bool_t               ConnectInput(AliAnalysisTask *task, Int_t islot,
                                     AliAnalysisDataContainer *cont);
   Bool_t               ConnectOutput(AliAnalysisTask *task, Int_t islot,
                                     AliAnalysisDataContainer *cont);
   // Garbage collection
   void                 CleanContainers();
   
   // Analysis initialization and execution, status
   Bool_t               InitAnalysis();
   Bool_t               IsInitialized() const {return fInitOK;}
   Bool_t               IsEventLoop() const {return TObject::TestBit(kEventLoop);}
   void                 ResetAnalysis();
   void                 ExecAnalysis(Option_t *option="");
   void                 FinishAnalysis();
   void                 PrintStatus(Option_t *option="all") const;

protected:
   void                 ImportWrappers(TList *source);
   void                 SetEventLoop(Bool_t flag=kTRUE) {TObject::SetBit(kEventLoop,flag);}

private:
   TTree                  *fTree;                //! Input tree in case of TSelector model
   AliVEventHandler       *fInputEventHandler;   //  Optional common input  event handler
   AliVEventHandler       *fOutputEventHandler;  //  Optional common output event handler
   AliVEventHandler       *fMCtruthEventHandler; //  Optional common MC Truth event handler
   AliVEventPool          *fEventPool;           //  Event pool for mixing analysis
   Long64_t                fCurrentEntry;        //! Current processed entry in the tree
   Long64_t                fNSysInfo;            // Event frequency for collecting system information
   EAliAnalysisExecMode    fMode;                // Execution mode
   Bool_t                  fInitOK;              // Initialisation done
   UInt_t                  fDebug;               // Debug level
   TString                 fSpecialOutputLocation; // URL/path where the special outputs will be copied
   TObjArray              *fTasks;               // List of analysis tasks
   TObjArray              *fTopTasks;            // List of top tasks
   TObjArray              *fZombies;             // List of zombie tasks
   TObjArray              *fContainers;          // List of all containers
   TObjArray              *fInputs;              // List of containers with input data
   TObjArray              *fOutputs;             // List of containers with results
   AliAnalysisSelector    *fSelector;            //! Current selector
   AliAnalysisGrid        *fGridHandler;         //! Grid handler plugin

   static AliAnalysisManager *fgAnalysisManager; //! static pointer to object instance
   ClassDef(AliAnalysisManager,3)  // Analysis manager class
};   
#endif
