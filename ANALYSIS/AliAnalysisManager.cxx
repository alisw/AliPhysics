/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */
// Author: Andrei Gheata, 31/05/2006

//==============================================================================
//   AliAnalysysManager - Manager analysis class. Allows creation of several
// analysis tasks and data containers storing their input/output. Allows 
// connecting/chaining tasks via shared data containers. Serializes the current
// event for all tasks depending only on initial input data.
//==============================================================================
//
//==============================================================================

#include "Riostream.h"

#include "TClass.h"
#include "TFile.h"
#include "TTree.h"
//#include "AliLog.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisDataSlot.h"

ClassImp(AliAnalysisManager)

//______________________________________________________________________________
AliAnalysisManager::AliAnalysisManager() : TSelector(),
                    fTree(NULL),
                    fInitOK(kFALSE),
                    fContainers(NULL),
                    fInputs(NULL),
                    fOutputs(NULL),
                    fTasks(NULL),
                    fTopTasks(NULL),
                    fZombies(NULL)
{
// Default constructor.
   if (TClass::IsCallingNew() != TClass::kDummyNew) {
      fContainers = new TObjArray();
      fInputs     = new TObjArray();
      fOutputs    = new TObjArray();  
      fTasks      = new TObjArray();
      fTopTasks   = new TObjArray();
      fZombies    = new TObjArray();
//      fStatus     = new AliAnalysisInfo(this);
   }
}

#ifdef NEVER
//______________________________________________________________________________
AliAnalysisManager::AliAnalysisManager(const AliAnalysisManager& other)
                   :TSelector(other),
                    fTree(NULL),
                    fInitOK(kFALSE),
                    fContainers(NULL),
                    fInputs(NULL),
                    fOutputs(NULL),
                    fTasks(NULL),
                    fTopTasks(NULL),
                    fZombies(NULL)
{
// Copy constructor.
   fInitOK     = other.fInitOK;
   fContainers = new TObjArray(*other.fContainers);
   fInputs     = new TObjArray(*other.fInputs);
   fOutputs    = new TObjArray(*other.fOutputs);
   fTasks      = new TObjArray(*other.fTasks);
   fTopTasks   = new TObjArray(*other.fTopTasks);
   fZombies    = new TObjArray(*other.fZombies);
//   fStatus     = new AliAnalysisInfo(this);
}
   
//______________________________________________________________________________
AliAnalysisManager& AliAnalysisManager::operator=(const AliAnalysisManager& other)
{
// Assignment
   if (&other != this) {
      TSelector::operator=(other);
      fTree       = other.fTree;
      fInitOK     = other.fInitOK;
      fContainers = new TObjArray(*other.fContainers);
      fInputs     = new TObjArray(*other.fInputs);
      fOutputs    = new TObjArray(*other.fOutputs);
      fTasks      = new TObjArray(*other.fTasks);
      fTopTasks   = new TObjArray(*other.fTopTasks);
      fZombies    = new TObjArray(*other.fZombies);
//      fStatus     = new AliAnalysisInfo(this);
   }
   return *this;
}
#endif

//______________________________________________________________________________
AliAnalysisManager::~AliAnalysisManager()
{
// Destructor.
   if (fContainers) {fContainers->Delete(); delete fContainers;}
   if (fInputs) delete fInputs;
   if (fOutputs) delete fOutputs;
   if (fTasks) {fTasks->Delete(); delete fTasks;}
   if (fTopTasks) delete fTopTasks;
   if (fZombies) delete fZombies;
}
//______________________________________________________________________________
Int_t AliAnalysisManager::GetEntry(Long64_t entry, Int_t getall)
{
// Read one entry of the tree or a whole branch.
   return fTree ? fTree->GetTree()->GetEntry(entry, getall) : 0;
}
   
//______________________________________________________________________________
void AliAnalysisManager::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses of the tree
  // will be set. It is normaly not necessary to make changes to the
  // generated code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running with PROOF.
   printf("AliAnalysisManager::Init(%s)\n", tree->GetName());
   if (!fInitOK) {
     cout<<"You have to call InitAnalysis first"<<endl;
     //AliError("You have to call InitAnalysis first");
      return;
   }   
   if (!tree) return;
   fTree = tree;
   AliAnalysisDataContainer *top = (AliAnalysisDataContainer*)fInputs->At(0);
   top->SetData(tree);
}

//______________________________________________________________________________
void AliAnalysisManager::Begin(TTree *tree)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
   printf("AliAnalysisManager::Begin(%s)\n", tree->GetName());
   Init(tree);
}

//______________________________________________________________________________
void AliAnalysisManager::SlaveBegin(TTree *tree)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).
   printf("AliAnalysisManager::SlaveBegin(%s)\n", tree->GetName());
   Init(tree);
}

//______________________________________________________________________________
Bool_t AliAnalysisManager::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.
   if (fTree) {
      TFile *curfile = fTree->GetCurrentFile();
      if (curfile) printf("AliAnalysisManager::Notify() file: %s\n", curfile->GetName());
   }
   return kTRUE;
}    

//______________________________________________________________________________
Bool_t AliAnalysisManager::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either TTree::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing
  // keyed objects with PROOF, the object is already loaded and is available
  // via the fObject pointer.
  //
  // This function should contain the "body" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.

  // WARNING when a selector is used with a TChain, you must use
  //  the pointer to the current TTree to call GetEntry(entry).
  //  The entry is always the local entry number in the current tree.
  //  Assuming that fChain is the pointer to the TChain being processed,
  //  use fChain->GetTree()->GetEntry(entry).
  
//   printf("AliAnalysisManager::Process(%lld)\n", entry);
   GetEntry(entry);
   ExecAnalysis();
   return kTRUE;
}

//______________________________________________________________________________
void AliAnalysisManager::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

   printf("AliAnalysisManager::SlaveTerminate()\n");
   if (!fOutput)
   {
     cout<<"ERROR: Output list not initialized."<<endl;
     //AliError("ERROR: Output list not initialized.");
     return;
   }
   TIter next(fOutputs);
   AliAnalysisDataContainer *output;
   while ((output=(AliAnalysisDataContainer *)next())) {
      output->SetDataOwned(kFALSE);
      fOutput->Add(output->GetData());
   }   
}

//______________________________________________________________________________
void AliAnalysisManager::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
   printf("AliAnalysisManager::Terminate()\n");
   AliAnalysisDataContainer *output;
   AliAnalysisTask *task;
   TIter next(fOutputs);
   while ((output=(AliAnalysisDataContainer *)next())) output->WriteData();
   TIter next1(fTasks);
   // Call Terminate() for tasks
   while ((task=(AliAnalysisTask*)next1())) task->Terminate();
}

//______________________________________________________________________________
void AliAnalysisManager::AddTask(AliAnalysisTask *task)
{
// Adds a user task to the global list of tasks.
   task->SetActive(kFALSE);
   fTasks->Add(task);
}  

//______________________________________________________________________________
AliAnalysisTask *AliAnalysisManager::GetTask(const char *name) const
{
// Retreive task by name.
   if (!fTasks) return NULL;
   return (AliAnalysisTask*)fTasks->FindObject(name);
}

//______________________________________________________________________________
AliAnalysisDataContainer *AliAnalysisManager::CreateContainer(const char *name, 
                                TClass *datatype, EAliAnalysisContType type)
{
// Create a data container of a certain type. Types can be:
//   kNormalContainer  = 0, used to exchange date between tasks
//   kInputContainer   = 1, used to store input data
//   kOutputContainer  = 2, used for posting results
   AliAnalysisDataContainer *cont = new AliAnalysisDataContainer(name, datatype);
   fContainers->Add(cont);
   switch (type) {
      case kInputContainer:
         fInputs->Add(cont);
         break;
      case kOutputContainer:
         fOutputs->Add(cont);
         break;
      case kNormalContainer:
         break;   
   }
   return cont;
}
         
//______________________________________________________________________________
Bool_t AliAnalysisManager::ConnectInput(AliAnalysisTask *task, Int_t islot,
                                        AliAnalysisDataContainer *cont)
{
// Connect input of an existing task to a data container.
   if (!fTasks->FindObject(task)) {
      AddTask(task);
      cout<<"Task "<<task->GetName()<<" not registered. Now owned by analysis manager"<<endl;
      //AliInfo(Form("Task %s not registered. Now owned by analysis manager", task->GetName()));
   } 
   Bool_t connected = task->ConnectInput(islot, cont);
   return connected;
}   

//______________________________________________________________________________
Bool_t AliAnalysisManager::ConnectOutput(AliAnalysisTask *task, Int_t islot,
                                        AliAnalysisDataContainer *cont)
{
// Connect output of an existing task to a data container.
   if (!fTasks->FindObject(task)) {
      AddTask(task);
      cout<<"Task "<<task->GetName()<<"not registered. Now owned by analysis manager"<<endl;
      //AliInfo(Form("Task %s not registered. Now owned by analysis manager", task->GetName()));
   } 
   Bool_t connected = task->ConnectOutput(islot, cont);
   return connected;
}   
                               
//______________________________________________________________________________
void AliAnalysisManager::CleanContainers()
{
// Clean data from all containers that have already finished all client tasks.
   TIter next(fContainers);
   AliAnalysisDataContainer *cont;
   while ((cont=(AliAnalysisDataContainer *)next())) {
      if (cont->IsOwnedData() && 
          cont->IsDataReady() && 
          cont->ClientsExecuted()) cont->DeleteData();
   }
}

//______________________________________________________________________________
Bool_t AliAnalysisManager::InitAnalysis()
{
// Initialization of analysis chain of tasks. Should be called after all tasks
// and data containers are properly connected
   // Check for input/output containers
   fInitOK = kFALSE;
//   if (!fInputs->GetEntriesFast()) {
//      AliError("No input container defined. At least one container should store input data");
//      return kFALSE;
//   }   
//   if (!fOutputs->GetEntriesFast()) {
//      AliError("No output container defined. At least one container should store output data");
//      return kFALSE;
//   }   
   // Check for top tasks (depending only on input data containers)
   if (!fTasks->First()) {
     cout<<"Analysis has no tasks !"<<endl;
     //AliError("Analysis have no tasks !");
      return kFALSE;
   }   
   TIter next(fTasks);
   AliAnalysisTask *task;
   AliAnalysisDataContainer *cont;
   Int_t ntop = 0;
   Int_t nzombies = 0;
   Bool_t iszombie = kFALSE;
   Bool_t istop = kTRUE;
   Int_t i;
   while ((task=(AliAnalysisTask*)next())) {
      istop = kTRUE;
      iszombie = kFALSE;
      Int_t ninputs = task->GetNinputs();
//      if (!ninputs) {
//         task->SetZombie();
//         fZombies->Add(task);
//         nzombies++;
//         AliWarning(Form("Task %s has no input slots defined ! Declared zombie...",task->GetName()));
//         continue;
//      }
      for (i=0; i<ninputs; i++) {
         cont = task->GetInputSlot(i)->GetContainer();
         if (!cont) {
            if (!iszombie) {
               task->SetZombie();
               fZombies->Add(task);
               nzombies++;
               iszombie = kTRUE;
            }   
            cout<<"Input slot "<<i<<" of task "<<task->GetName()<<" has no container connected ! Declared zombie..."<<endl;
	    //AliWarning(Form("Input slot %i of task %s has no container connected ! Declared zombie...",i,task->GetName()));
         }
         if (iszombie) continue;
         // Check if cont is an input container
         if (istop && !fInputs->FindObject(cont)) istop=kFALSE;
         // Connect to parent task
      }
      if (istop) {
         ntop++;
         fTopTasks->Add(task);
      }
   }
   if (!ntop) {
     cout<<"No top task defined. At least one task should be connected only to input containers"<<endl;
     //AliError("No top task defined. At least one task should be connected only to input containers");
      return kFALSE;
   }                        
   // Check now if there are orphan tasks
   for (i=0; i<ntop; i++) {
      task = (AliAnalysisTask*)fTopTasks->At(i);
      task->SetUsed();
   }
   Int_t norphans = 0;
   next.Reset();
   while ((task=(AliAnalysisTask*)next())) {
      if (!task->IsUsed()) {
         norphans++;
         cout<<"Task "<<task->GetName()<<" is orphan"<<endl;
	 //AliWarning(Form("Task %s is orphan",task->GetName()));
      }   
   }          
   // Check the task hierarchy (no parent task should depend on data provided
   // by a daughter task)
   for (i=0; i<ntop; i++) {
      task = (AliAnalysisTask*)fTopTasks->At(i);
      if (task->CheckCircularDeps()) {
	cout<<"Found illegal circular dependencies between following tasks:"<<endl;
	//AliError("Found illegal circular dependencies between following tasks:");
         PrintStatus("dep");
         return kFALSE;
      }   
   }
   fInitOK = kTRUE;
   return kTRUE;
}   

//______________________________________________________________________________
void AliAnalysisManager::PrintStatus(Option_t *option) const
{
// Print task hierarchy.
   TIter next(fTopTasks);
   AliAnalysisTask *task;
   while ((task=(AliAnalysisTask*)next()))
      task->PrintTask(option);
}

//______________________________________________________________________________
void AliAnalysisManager::ResetAnalysis()
{
// Reset all execution flags and clean containers.
   CleanContainers();
}

//______________________________________________________________________________
void AliAnalysisManager::ExecAnalysis(Option_t *option)
{
// Execute analysis.
   if (!fInitOK) {
     cout<<"Analysis manager was not initialized !"<<endl;
     //AliError("Analysis manager was not initialized !");
      return;
   }   
   AliAnalysisTask *task;
   // Check if the top tree is active.
   if (fTree) {
      TIter next(fTasks);
   // De-activate all tasks
      while ((task=(AliAnalysisTask*)next())) task->SetActive(kFALSE);
      AliAnalysisDataContainer *cont = (AliAnalysisDataContainer*)fInputs->At(0);
      if (!cont) {
	cout<<"Cannot execute analysis in TSelector mode without at least one top container"<<endl;
	//AliError("Cannot execute analysis in TSelector mode without at least one top container");
         return;
      }   
      cont->SetData(fTree); // This will notify all consumers
      TIter next1(cont->GetConsumers());
      while ((task=(AliAnalysisTask*)next1())) {
//         task->SetActive(kTRUE);
         task->ExecuteTask(option);
      }
      return;
   }   
   // The event loop is not controlled by TSelector   
   TIter next2(fTopTasks);
   while ((task=(AliAnalysisTask*)next2())) {
      task->SetActive(kTRUE);
      printf("executing %s\n", task->GetName());
      task->ExecuteTask(option);
   }   
}

//______________________________________________________________________________
void AliAnalysisManager::FinishAnalysis()
{
// Finish analysis.
}
