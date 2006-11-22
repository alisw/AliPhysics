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

#ifndef ROOT_TSelector
#include "TSelector.h"
#endif

class TClass;
class AliAnalysisDataContainer;
class AliAnalysisTask;

class AliAnalysisManager : public TSelector {

public:

enum EAliAnalysisContType {
   kNormalContainer  = 0,
   kInputContainer   = 1,
   kOutputContainer  = 2
};   
   AliAnalysisManager();
   virtual            ~AliAnalysisManager();
   
   // Selector-specific methods
   virtual void        Init(TTree *tree);   
   virtual void        Begin(TTree *tree);
   virtual Bool_t      Notify();
   virtual void        SlaveBegin(TTree *tree);
   virtual Bool_t      ProcessCut(Long64_t entry) {return Process(entry);}
   virtual Bool_t      Process(Long64_t entry);
   virtual Int_t       GetEntry(Long64_t entry, Int_t getall = 0);
   virtual void        SlaveTerminate();
   virtual void        Terminate();

   // Getters
   TTree              *GetTree() const       {return fTree;}
   TObjArray          *GetContainers() const {return fContainers;}
   TObjArray          *GetInputs() const     {return fInputs;}
   TObjArray          *GetOutputs() const    {return fOutputs;}
   TObjArray          *GetTasks() const      {return fTasks;}
   TObjArray          *GetTopTasks() const   {return fTopTasks;}
   TObjArray          *GetZombieTasks() const {return fZombies;}
//   AliAnalysisInfo    *GetStatus() const     {return fStatus;}  

   // Container handling
   AliAnalysisDataContainer *CreateContainer(const char *name, TClass *datatype, 
                                EAliAnalysisContType type=kNormalContainer);
   
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
   TTree               *fTree;         // Input tree in case of TSelector model
   Bool_t               fInitOK;       // Initialisation done
   TObjArray           *fContainers;   //-> List of all containers
   TObjArray           *fInputs;       //-> List of containers with input data
   TObjArray           *fOutputs;      //-> List of containers with results
   TObjArray           *fTasks;        //-> List of analysis tasks
   TObjArray           *fTopTasks;     //-> List of top tasks
   TObjArray           *fZombies;      //-> List of zombie tasks
//   AliAnalysisInfo     *fStatus;       // Analysis info object

private:
   AliAnalysisManager(const AliAnalysisManager& other);
   AliAnalysisManager& operator=(const AliAnalysisManager& other);


   ClassDef(AliAnalysisManager,1)  // Analysis manager class
};   
#endif
