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
class AliAnalysisDataContainer;
class AliAnalysisTask;

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
   kGridAnalysis     = 2
};

   AliAnalysisManager();
   AliAnalysisManager(const char *name, const char *title="");
   virtual            ~AliAnalysisManager();

   AliAnalysisManager(const AliAnalysisManager& other);
   AliAnalysisManager& operator=(const AliAnalysisManager& other);
   
   void                StartAnalysis(const char *type="local", TTree *tree=0);

   virtual void        Init(TTree *tree);   
   virtual void        Begin(TTree *tree);
   virtual Bool_t      Notify();
   virtual void        SlaveBegin(TTree *tree);
   virtual Bool_t      ProcessCut(Long64_t entry) {return Process(entry);}
   virtual Bool_t      Process(Long64_t entry);
   virtual Int_t       GetEntry(Long64_t entry, Int_t getall = 0);
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

   void                SetAnalysisType(EAliAnalysisExecMode mode) {fMode = mode;}
   void                SetCurrentEntry(Long64_t entry) {fCurrentEntry = entry;}
   void                SetDebugLevel(UInt_t level) {fDebug = level;}

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
   void                 ResetAnalysis();
   void                 ExecAnalysis(Option_t *option="");
   void                 FinishAnalysis();
   void                 PrintStatus(Option_t *option="all") const;

protected:
   void                 ReplaceOutputContainers(TList *source);

private:
   TTree               *fTree;         //! Input tree in case of TSelector model
   Long64_t             fCurrentEntry; //! Current processed entry in the tree
   EAliAnalysisExecMode fMode;         // Execution mode
   Bool_t               fInitOK;       // Initialisation done
   UInt_t               fDebug;        // Debug level
   TObjArray           *fTasks;        // List of analysis tasks
   TObjArray           *fTopTasks;     // List of top tasks
   TObjArray           *fZombies;      // List of zombie tasks
   TObjArray           *fContainers;   // List of all containers
   TObjArray           *fInputs;       // List of containers with input data
   TObjArray           *fOutputs;      // List of containers with results

   static AliAnalysisManager *fgAnalysisManager; //! static pointer to object instance
   ClassDef(AliAnalysisManager,1)  // Analysis manager class
};   
#endif
