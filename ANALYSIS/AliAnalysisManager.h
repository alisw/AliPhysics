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
#ifndef ROOT_THashTable
#include <THashTable.h>
#endif
#ifndef ROOT_Riostream
#include <Riostream.h>
#endif

class TClass;
class TTree;
class TFile;
class TStopwatch;
class AliAnalysisSelector;
class AliAnalysisDataContainer;
class AliAnalysisTask;
class AliVEventHandler;
class AliVEventPool;
class AliAnalysisGrid;
class AliAnalysisStatistics;


class AliAnalysisManager : public TNamed {

public:

enum EAliAnalysisContType {
   kExchangeContainer  = 0,   // use to exchange data between tasks
   kInputContainer     = 1,   // use for the task private input data
   kOutputContainer    = 2,   // use for the task private output data
   kParamContainer     = 3    // use for storing task private parameters/cuts per analysis session
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
   kUseDataSet       = BIT(16),
   kSaveCanvases     = BIT(17),
   kExternalLoop     = BIT(18),
   kSkipTerminate    = BIT(19),
   kUseProgressBar   = BIT(20),
   kTrueNotify       = BIT(21)
};   

   AliAnalysisManager(const char *name = "mgr", const char *title="");
   virtual            ~AliAnalysisManager();

   AliAnalysisManager(const AliAnalysisManager& other);
   AliAnalysisManager& operator=(const AliAnalysisManager& other);
   
   // Event loop control
   virtual Int_t       GetEntry(Long64_t entry, Int_t getall = 0);
   virtual Bool_t      Init(TTree *tree);
   virtual Bool_t      Notify();
   virtual Bool_t      ProcessCut(Long64_t entry) {return Process(entry);}
   virtual Bool_t      Process(Long64_t entry);
   TFile              *OpenProofFile(AliAnalysisDataContainer *cont, const char *option, const char *extaod="");
   static TFile       *OpenFile(AliAnalysisDataContainer *cont, const char *option, Bool_t ignoreProof=kFALSE);
   void                PackOutput(TList *target);
   void                RegisterExtraFile(const char *fname);
   Long64_t            StartAnalysis(const char *type, TTree * const tree, Long64_t nentries=1234567890, Long64_t firstentry=0);
   Long64_t            StartAnalysis(const char *type, const char *dataset, Long64_t nentries=1234567890, Long64_t firstentry=0);
   Long64_t            StartAnalysis(const char *type, Long64_t nentries=1234567890, Long64_t firstentry=0);
   virtual void        SlaveBegin(TTree *tree);
   virtual void        Terminate();
   void                UnpackOutput(TList *source);

   // Getters/Setters
   static AliAnalysisManager *GetAnalysisManager() {return fgAnalysisManager;}
   EAliAnalysisExecMode 
                       GetAnalysisType() const    {return fMode;}
   void                GetAnalysisTypeString(TString &type) const;                    
   Bool_t              GetAutoBranchLoading() const {return fAutoBranchHandling;} 
   static const char  *GetCommonFileName()        {return fgCommonFileName.Data();}
   AliAnalysisDataContainer *
                       GetCommonInputContainer() const  {return fCommonInput;}
   AliAnalysisDataContainer *
                       GetCommonOutputContainer() const {return fCommonOutput;}
   TObjArray          *GetContainers() const      {return fContainers;}
   Long64_t            GetCurrentEntry() const    {return fCurrentEntry;}
   UInt_t              GetDebugLevel() const      {return fDebug;}
   TString             GetExtraFiles() const      {return fExtraFiles;}
   AliVEventPool*      GetEventPool()  const      {return fEventPool;}
   Bool_t              GetFileFromWrapper(const char *filename, const TList *source);
   static Int_t        GetRunFromAlienPath(const char *path);
   AliAnalysisGrid*    GetGridHandler()           {return fGridHandler;}
   TObjArray          *GetInputs() const          {return fInputs;}
   AliVEventHandler*   GetInputEventHandler() const   {return fInputEventHandler;}
   AliVEventHandler*   GetMCtruthEventHandler() const {return fMCtruthEventHandler;}
   AliVEventHandler*   GetOutputEventHandler() const  {return fOutputEventHandler;}
   TObjArray          *GetOutputs() const         {return fOutputs;}
   TObjArray          *GetParamOutputs() const    {return fParamCont;}
   Int_t               GetRunFromPath() const     {return fRunFromPath;}
   const char         *GetRequestedBranches() const {return fRequestedBranches.Data();}
   TObjArray          *GetTasks() const           {return fTasks;}
   TObjArray          *GetTopTasks() const        {return fTopTasks;}
   TTree              *GetTree() const            {return fTree;}
   TObjArray          *GetZombieTasks() const     {return fZombies;}
   static Bool_t       IsPipe(std::ostream &out);
   Bool_t              IsProofMode() const        {return (fMode==kProofAnalysis)?kTRUE:kFALSE;}
   Bool_t              IsRemote() const           {return fIsRemote;}
   Bool_t              IsUsingDataSet() const     {return TObject::TestBit(kUseDataSet);}
   void                LoadBranch(const char *n)  { if(fAutoBranchHandling) return; DoLoadBranch(n); }
   void                SetAnalysisType(EAliAnalysisExecMode mode) {fMode = mode;}
   void                SetAutoBranchLoading(Bool_t b) { fAutoBranchHandling = b; }
   void                SetCurrentEntry(Long64_t entry)            {fCurrentEntry = entry;}
   void                SetCollectSysInfoEach(Int_t nevents=0)     {fNSysInfo = nevents;}
   static void         SetCommonFileName(const char *name)        {fgCommonFileName = name;}
   void                SetDebugLevel(UInt_t level);
   void                SetDisableBranches(Bool_t disable=kTRUE)   {TObject::SetBit(kDisableBranches,disable);}
   void                SetExternalLoop(Bool_t flag)               {TObject::SetBit(kExternalLoop,flag);}
   void                SetEventPool(AliVEventPool* const epool)   {fEventPool = epool;}
   void                SetGridHandler(AliAnalysisGrid * const handler) {fGridHandler = handler;}
   void                SetInputEventHandler(AliVEventHandler* const handler);
   void                SetMCtruthEventHandler(AliVEventHandler* const handler) {fMCtruthEventHandler = handler;}
   void                SetNSysInfo(Long64_t nevents)              {fNSysInfo = nevents;}
   void                SetOutputEventHandler(AliVEventHandler* const handler);
   void                SetRunFromPath(Int_t run)                  {fRunFromPath = run;}
   void                SetSelector(AliAnalysisSelector * const sel)      {fSelector = sel;}
   void                SetSaveCanvases(Bool_t flag=kTRUE)         {TObject::SetBit(kSaveCanvases,flag);}
   void                SetSkipTerminate(Bool_t flag)              {TObject::SetBit(kSkipTerminate,flag);}
   void                SetUseProgressBar(Bool_t flag, Int_t freq=1);
   void                SetSpecialOutputLocation(const char *loc)  {fSpecialOutputLocation = loc;}

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
   void                 AddBranches(const char *branches);
   void                 AddStatisticsTask(UInt_t offlineMask=0);
   void                 CheckBranches(Bool_t load=kFALSE);
   Bool_t               CheckTasks() const;
   void                 CountEvent(Int_t ninput, Int_t nprocessed, Int_t nfailed, Int_t naccepted);
   Bool_t               InitAnalysis();
   Bool_t               IsInitialized() const {return fInitOK;}
   Bool_t               IsExternalLoop() const {return TObject::TestBit(kExternalLoop);}
   Bool_t               IsEventLoop() const {return TObject::TestBit(kEventLoop);}
   Bool_t               IsSkipTerminate() const {return TObject::TestBit(kSkipTerminate);}
   void                 ResetAnalysis();
   void                 ExecAnalysis(Option_t *option="");
   void                 PrintStatus(Option_t *option="all") const;
   void                 ProfileTask(const char *name, const char *option="VM") const;
   void                 ProfileTask(Int_t itop, const char *option="VM") const;
   static void          ProgressBar(const char *opname, Long64_t current, Long64_t size, TStopwatch * const watch=0, Bool_t last=kFALSE, Bool_t refresh=kFALSE);
   void                 AddStatisticsMsg(const char *line);
   const char          *GetStatisticsMsg() const {return fStatisticsMsg.Data();}
   const AliAnalysisStatistics *GetStatistics() const {return fStatistics;}
   void                 SetStatistics(AliAnalysisStatistics *stat) {fStatistics = stat;}
   void                 WriteStatisticsMsg(Int_t nevents);
   Int_t                GetNcalls() const {return fNcalls;}
   Bool_t               ValidateOutputFiles() const;
   
   static const char*   GetOADBPath();

protected:
   void                 ImportWrappers(TList *source);
   void                 SetEventLoop(Bool_t flag=kTRUE) {TObject::SetBit(kEventLoop,flag);}
   void                 DoLoadBranch(const char *name);

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
   Bool_t                  fIsRemote;            //! Flag is set for remote analysis
   UInt_t                  fDebug;               // Debug level
   TString                 fSpecialOutputLocation; // URL/path where the special outputs will be copied
   TObjArray              *fTasks;               // List of analysis tasks
   TObjArray              *fTopTasks;            // List of top tasks
   TObjArray              *fZombies;             // List of zombie tasks
   TObjArray              *fContainers;          // List of all containers
   TObjArray              *fInputs;              // List of containers with input data
   TObjArray              *fOutputs;             // List of containers with results
   TObjArray              *fParamCont;           // List of containers with results
   AliAnalysisDataContainer *fCommonInput;       // Common input container
   AliAnalysisDataContainer *fCommonOutput;      // Common output container
   AliAnalysisSelector    *fSelector;            //! Current selector
   AliAnalysisGrid        *fGridHandler;         //! Grid handler plugin
   TString                 fExtraFiles;          // List of extra files to be merged
   Bool_t                  fAutoBranchHandling;  // def=kTRUE, turn off if you use LoadBranch
   THashTable              fTable;               // keep branch ptrs in case of manual branch loading
   Int_t                   fRunFromPath;         // Run number retrieved from path to input data
   Int_t                   fNcalls;              // Total number of calls (events) of ExecAnalysis
   Long64_t                fMaxEntries;          // Maximum number of entries
   static Int_t            fPBUpdateFreq;        // Progress bar update freq.
   TString                 fStatisticsMsg;       // Statistics user message
   TString                 fRequestedBranches;   // Requested branch names
   AliAnalysisStatistics  *fStatistics;          // Statistics info about input events
   static TString          fgCommonFileName;     //! Common output file name (not streamed)
   static AliAnalysisManager *fgAnalysisManager; //! static pointer to object instance
   ClassDef(AliAnalysisManager,13)  // Analysis manager class
};   
#endif
