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

#include <Riostream.h>

#include <TClass.h>
#include <TFile.h>
#include <TMethodCall.h>
#include <TChain.h>
#include <TSystem.h>
#include <TROOT.h>

#include "AliAnalysisTask.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisDataSlot.h"
#include "AliVEventHandler.h"
#include "AliAnalysisManager.h"

ClassImp(AliAnalysisManager)

AliAnalysisManager *AliAnalysisManager::fgAnalysisManager = NULL;

//______________________________________________________________________________
AliAnalysisManager::AliAnalysisManager(const char *name, const char *title)
                   :TNamed(name,title),
                    fTree(NULL),
		    fInputEventHandler(NULL),
		    fOutputEventHandler(NULL),
		    fMCtruthEventHandler(NULL),
                    fCurrentEntry(-1),
                    fMode(kLocalAnalysis),
                    fInitOK(kFALSE),
                    fDebug(0),
                    fTasks(NULL),
                    fTopTasks(NULL),
                    fZombies(NULL),
                    fContainers(NULL),
                    fInputs(NULL),
                    fOutputs(NULL)
{
// Default constructor.
   fgAnalysisManager = this;
   fTasks      = new TObjArray();
   fTopTasks   = new TObjArray();
   fZombies    = new TObjArray();
   fContainers = new TObjArray();
   fInputs     = new TObjArray();
   fOutputs    = new TObjArray();
   SetEventLoop(kTRUE);
}

//______________________________________________________________________________
AliAnalysisManager::AliAnalysisManager(const AliAnalysisManager& other)
                   :TNamed(other),
                    fTree(NULL),
		    fInputEventHandler(NULL),
		    fOutputEventHandler(NULL),
		    fMCtruthEventHandler(NULL),
                    fCurrentEntry(-1),
                    fMode(other.fMode),
                    fInitOK(other.fInitOK),
                    fDebug(other.fDebug),
                    fTasks(NULL),
                    fTopTasks(NULL),
                    fZombies(NULL),
                    fContainers(NULL),
                    fInputs(NULL),
                    fOutputs(NULL)
{
// Copy constructor.
   fTasks      = new TObjArray(*other.fTasks);
   fTopTasks   = new TObjArray(*other.fTopTasks);
   fZombies    = new TObjArray(*other.fZombies);
   fContainers = new TObjArray(*other.fContainers);
   fInputs     = new TObjArray(*other.fInputs);
   fOutputs    = new TObjArray(*other.fOutputs);
   fgAnalysisManager = this;
}
   
//______________________________________________________________________________
AliAnalysisManager& AliAnalysisManager::operator=(const AliAnalysisManager& other)
{
// Assignment
   if (&other != this) {
      TNamed::operator=(other);
      fInputEventHandler   = other.fInputEventHandler;
      fOutputEventHandler  = other.fOutputEventHandler;
      fMCtruthEventHandler = other.fMCtruthEventHandler;
      fTree       = NULL;
      fCurrentEntry = -1;
      fMode       = other.fMode;
      fInitOK     = other.fInitOK;
      fDebug      = other.fDebug;
      fTasks      = new TObjArray(*other.fTasks);
      fTopTasks   = new TObjArray(*other.fTopTasks);
      fZombies    = new TObjArray(*other.fZombies);
      fContainers = new TObjArray(*other.fContainers);
      fInputs     = new TObjArray(*other.fInputs);
      fOutputs    = new TObjArray(*other.fOutputs);
      fgAnalysisManager = this;
   }
   return *this;
}

//______________________________________________________________________________
AliAnalysisManager::~AliAnalysisManager()
{
// Destructor.
   if (fTasks) {fTasks->Delete(); delete fTasks;}
   if (fTopTasks) delete fTopTasks;
   if (fZombies) delete fZombies;
   if (fContainers) {fContainers->Delete(); delete fContainers;}
   if (fInputs) delete fInputs;
   if (fOutputs) delete fOutputs;
   if (fgAnalysisManager==this) fgAnalysisManager = NULL;
}

//______________________________________________________________________________
Int_t AliAnalysisManager::GetEntry(Long64_t entry, Int_t getall)
{
// Read one entry of the tree or a whole branch.
   if (fDebug > 1) {
      cout << "== AliAnalysisManager::GetEntry()" << endl;
   }   
   fCurrentEntry = entry;
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
   if (!tree) return;
   if (fDebug > 1) {
      printf("->AliAnalysisManager::Init(%s)\n", tree->GetName());
   }

   if (fInputEventHandler) {
       fInputEventHandler->SetInputTree(tree);
       fInputEventHandler->InitIO("proof");
   }

   if (!fInitOK) InitAnalysis();
   if (!fInitOK) return;
   fTree = tree;
   AliAnalysisDataContainer *top = (AliAnalysisDataContainer*)fInputs->At(0);
   if (!top) {
      cout<<"Error: No top input container !" <<endl;
      return;
   }
   top->SetData(tree);
   if (fDebug > 1) {
      printf("<-AliAnalysisManager::Init(%s)\n", tree->GetName());
   }
}

//______________________________________________________________________________
void AliAnalysisManager::Begin(TTree *tree)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
   if (fDebug > 1) {
      cout << "AliAnalysisManager::Begin()" << endl;
   }   
   Init(tree);
}

//______________________________________________________________________________
void AliAnalysisManager::SlaveBegin(TTree *tree)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).
   if (fDebug > 1) {
      cout << "->AliAnalysisManager::SlaveBegin()" << endl;
   }
   // Call InitIO of EventHandler
   if (fOutputEventHandler) {
       if (fMode == kProofAnalysis) {
	   fOutputEventHandler->InitIO("proof");
       } else {
	   fOutputEventHandler->InitIO("local");
       }
   }
   if (fInputEventHandler) {
       if (fMode == kProofAnalysis) {
	   fInputEventHandler->SetInputTree(tree);
	   fInputEventHandler->InitIO("proof");
       } else {
	   fInputEventHandler->SetInputTree(tree);
	   fInputEventHandler->InitIO("local");
       }
   }

   if (fMCtruthEventHandler) {
       if (fMode == kProofAnalysis) {
	   fMCtruthEventHandler->InitIO("proof");
       } else {
	   fMCtruthEventHandler->InitIO("local");
       }
   }
   
   //
   TIter next(fTasks);
   AliAnalysisTask *task;
   // Call CreateOutputObjects for all tasks
   while ((task=(AliAnalysisTask*)next())) {
      TDirectory *curdir = gDirectory;
      task->CreateOutputObjects();
      if (curdir) curdir->cd();
   }   
   if (fMode == kLocalAnalysis) 
       Init(tree);   
   if (fDebug > 1) {
      cout << "<-AliAnalysisManager::SlaveBegin()" << endl;
   }
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
	if (curfile && fDebug>1) printf("AliAnalysisManager::Notify() file: %s\n", curfile->GetName());
	TIter next(fTasks);
	AliAnalysisTask *task;
	// Call Notify for all tasks
	while ((task=(AliAnalysisTask*)next())) 
	    task->Notify();
	
	// Call Notify of the event handlers
	if (fInputEventHandler) {
	    fInputEventHandler->Notify(curfile->GetName());
	}

	if (fOutputEventHandler) {
	    fOutputEventHandler->Notify(curfile->GetName());
	}

	if (fMCtruthEventHandler) {
	    fMCtruthEventHandler->Notify(curfile->GetName());
	}

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
   if (fDebug > 1) {
      cout << "->AliAnalysisManager::Process()" << endl;
   }
   if (fInputEventHandler)   fInputEventHandler  ->BeginEvent(entry);
   if (fOutputEventHandler)  fOutputEventHandler ->BeginEvent(entry);
   if (fMCtruthEventHandler) fMCtruthEventHandler->BeginEvent(entry);
   
   GetEntry(entry);
   ExecAnalysis();
   if (fDebug > 1) {
      cout << "<-AliAnalysisManager::Process()" << endl;
   }
   return kTRUE;
}

//______________________________________________________________________________
void AliAnalysisManager::PackOutput(TList *target)
{
  // Pack all output data containers in the output list. Called at SlaveTerminate
  // stage in PROOF case for each slave.
   if (fDebug > 1) {
      cout << "->AliAnalysisManager::PackOutput()" << endl;
   }   
   if (!target) {
      Error("PackOutput", "No target. Aborting.");
      return;
   }
   if (fInputEventHandler)   fInputEventHandler  ->Terminate();
   if (fOutputEventHandler)  fOutputEventHandler ->Terminate();
   if (fMCtruthEventHandler) fMCtruthEventHandler->Terminate();
   
   if (fMode == kProofAnalysis) {
      TIter next(fOutputs);
      AliAnalysisDataContainer *output;
      while ((output=(AliAnalysisDataContainer*)next())) {
         if (output->GetData()) {
            if (output->GetProducer()->IsPostEventLoop()) continue;
            AliAnalysisDataWrapper *wrap = output->ExportData();
            // Output wrappers must delete data after merging (AG 13/11/07)
            wrap->SetDeleteData(kTRUE);
            if (fDebug > 1) printf("   Packing container %s...\n", output->GetName());
            target->Add(wrap);
         }   
      }
   } 
   if (fDebug > 1) {
      printf("<-AliAnalysisManager::PackOutput: output list contains %d containers\n", target->GetSize());
   }
}

//______________________________________________________________________________
void AliAnalysisManager::ImportWrappers(TList *source)
{
// Import data in output containers from wrappers coming in source.
   if (fDebug > 1) {
      cout << "->AliAnalysisManager::ImportWrappers()" << endl;
   }   
   TIter next(fOutputs);
   AliAnalysisDataContainer *cont;
   AliAnalysisDataWrapper   *wrap;
   Int_t icont = 0;
   while ((cont=(AliAnalysisDataContainer*)next())) {
      if (cont->GetProducer()->IsPostEventLoop()) continue;
      wrap = (AliAnalysisDataWrapper*)source->FindObject(cont->GetName());
      if (!wrap && fDebug>1) {
         printf("(WW) ImportWrappers: container %s not found in analysis output !\n", cont->GetName());
         continue;
      }
      icont++;
      if (fDebug > 1) printf("   Importing data for container %s\n", wrap->GetName());
      if (cont->GetFileName()) printf("    -> %s\n", cont->GetFileName());
      cont->ImportData(wrap);
   }         
   if (fDebug > 1) {
      cout << "<-AliAnalysisManager::ImportWrappers(): "<< icont << " containers imported" << endl;
   }   
}

//______________________________________________________________________________
void AliAnalysisManager::UnpackOutput(TList *source)
{
  // Called by AliAnalysisSelector::Terminate. Output containers should
  // be in source in the same order as in fOutputs.
   if (fDebug > 1) {
      cout << "->AliAnalysisManager::UnpackOutput()" << endl;
   }   
   if (!source) {
      Error("UnpackOutput", "No target. Aborting.");
      return;
   }
   if (fDebug > 1) {
      printf("   Source list contains %d containers\n", source->GetSize());
   }   

   if (fMode == kProofAnalysis) ImportWrappers(source);

   TIter next(fOutputs);
   AliAnalysisDataContainer *output;
   while ((output=(AliAnalysisDataContainer*)next())) {
      if (!output->GetData()) continue;
      // Check if there are client tasks that run post event loop
      if (output->HasConsumers()) {
         // Disable event loop semaphore
         output->SetPostEventLoop(kTRUE);
         TObjArray *list = output->GetConsumers();
         Int_t ncons = list->GetEntriesFast();
         for (Int_t i=0; i<ncons; i++) {
            AliAnalysisTask *task = (AliAnalysisTask*)list->At(i);
            task->CheckNotify(kTRUE);
            // If task is active, execute it
            if (task->IsPostEventLoop() && task->IsActive()) {
               if (fDebug > 1) {
                  cout << "== Executing post event loop task " << task->GetName() << endl;
               }                  
               task->ExecuteTask();
            }   
         }
      }   
      // Check if the output need to be written to a file.
      const char *filename = output->GetFileName();
      if (!(strcmp(filename, "default"))) {
	  if (fOutputEventHandler) filename = fOutputEventHandler->GetOutputFileName();
      }
      
      if (!filename || !strlen(filename)) continue;
      TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
      if (file) file->cd();
      else      file = new TFile(filename, "RECREATE");
      if (file->IsZombie()) continue;
      // Reparent data to this file
      TMethodCall callEnv;
      if (output->GetData()->IsA())
         callEnv.InitWithPrototype(output->GetData()->IsA(), "SetDirectory", "TDirectory*");
      if (callEnv.IsValid()) {
         callEnv.SetParam((Long_t) file);
         callEnv.Execute(output->GetData());
      }
      output->GetData()->Write();
   }
   if (fDebug > 1) {
      cout << "<-AliAnalysisManager::UnpackOutput()" << endl;
   }   
}

//______________________________________________________________________________
void AliAnalysisManager::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically.
   if (fDebug > 1) {
      cout << "->AliAnalysisManager::Terminate()" << endl;
   }   
   AliAnalysisTask *task;
   TIter next(fTasks);
   // Call Terminate() for tasks
   while ((task=(AliAnalysisTask*)next())) task->Terminate();
   if (fDebug > 1) {
      cout << "<-AliAnalysisManager::Terminate()" << endl;
   }   
   //
   if (fInputEventHandler)   fInputEventHandler  ->TerminateIO();
   if (fOutputEventHandler)  fOutputEventHandler ->TerminateIO();
   if (fMCtruthEventHandler) fMCtruthEventHandler->TerminateIO();
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
                                TClass *datatype, EAliAnalysisContType type, const char *filename)
{
// Create a data container of a certain type. Types can be:
//   kExchangeContainer  = 0, used to exchange date between tasks
//   kInputContainer   = 1, used to store input data
//   kOutputContainer  = 2, used for posting results
   if (fContainers->FindObject(name)) {
      Error("CreateContainer","A container named %s already defined !\n",name);
      return NULL;
   }   
   AliAnalysisDataContainer *cont = new AliAnalysisDataContainer(name, datatype);
   fContainers->Add(cont);
   switch (type) {
      case kInputContainer:
         fInputs->Add(cont);
         break;
      case kOutputContainer:
         fOutputs->Add(cont);
         if (filename && strlen(filename)) cont->SetFileName(filename);
         break;
      case kExchangeContainer:
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
      Warning("ConnectInput", "Task %s not registered. Now owned by analysis manager", task->GetName());
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
      Warning("ConnectOutput", "Task %s not registered. Now owned by analysis manager", task->GetName());
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
   // Check for top tasks (depending only on input data containers)
   if (!fTasks->First()) {
      Error("InitAnalysis", "Analysis has no tasks !");
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
      for (i=0; i<ninputs; i++) {
         cont = task->GetInputSlot(i)->GetContainer();
         if (!cont) {
            if (!iszombie) {
               task->SetZombie();
               fZombies->Add(task);
               nzombies++;
               iszombie = kTRUE;
            }   
            Error("InitAnalysis", "Input slot %d of task %s has no container connected ! Declared zombie...", 
                  i, task->GetName()); 
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
      Error("InitAnalysis", "No top task defined. At least one task should be connected only to input containers");
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
         Warning("InitAnalysis", "Task %s is orphan", task->GetName());
      }   
   }          
   // Check the task hierarchy (no parent task should depend on data provided
   // by a daughter task)
   for (i=0; i<ntop; i++) {
      task = (AliAnalysisTask*)fTopTasks->At(i);
      if (task->CheckCircularDeps()) {
         Error("InitAnalysis", "Found illegal circular dependencies between following tasks:");
         PrintStatus("dep");
         return kFALSE;
      }   
   }
   // Check that all containers feeding post-event loop tasks are in the outputs list
   TIter nextcont(fContainers); // loop over all containers
   while ((cont=(AliAnalysisDataContainer*)nextcont())) {
      if (!cont->IsPostEventLoop() && !fOutputs->FindObject(cont)) {
         if (cont->HasConsumers()) {
         // Check if one of the consumers is post event loop
            TIter nextconsumer(cont->GetConsumers());
            while ((task=(AliAnalysisTask*)nextconsumer())) {
               if (task->IsPostEventLoop()) {
                  fOutputs->Add(cont);
                  break;
               }
            }
         }
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
void AliAnalysisManager::StartAnalysis(const char *type, TTree *tree)
{
// Start analysis for this manager. Analysis task can be: LOCAL, PROOF or GRID.
   if (!fInitOK) {
      Error("StartAnalysis","Analysis manager was not initialized !");
      return;
   }
   if (fDebug>1) {
      cout << "StartAnalysis: " << GetName() << endl;   
   }   
   TString anaType = type;
   anaType.ToLower();
   fMode = kLocalAnalysis;
   if (tree) {
      if (anaType.Contains("proof"))     fMode = kProofAnalysis;
      else if (anaType.Contains("grid")) fMode = kGridAnalysis;
   }
   if (fMode == kGridAnalysis) {
      Warning("StartAnalysis", "GRID analysis mode not implemented. Running local.");
      fMode = kLocalAnalysis;
   }
   char line[128];
   SetEventLoop(kFALSE);
   // Disable all branches if requested and set event loop mode
   if (tree) {
      if (TestBit(kDisableBranches)) {
         printf("Disabling all branches...\n");
//         tree->SetBranchStatus("*",0); // not yet working
      }   
      SetEventLoop(kTRUE);
   }   

   TChain *chain = dynamic_cast<TChain*>(tree);

   // Initialize locally all tasks
   TIter next(fTasks);
   AliAnalysisTask *task;
   while ((task=(AliAnalysisTask*)next())) {
      task->LocalInit();
   }
   
   switch (fMode) {
      case kLocalAnalysis:
         if (!tree) {
            TIter next(fTasks);
            AliAnalysisTask *task;
            // Call CreateOutputObjects for all tasks
            while ((task=(AliAnalysisTask*)next())) {
               TDirectory *curdir = gDirectory;
               task->CreateOutputObjects();
               if (curdir) curdir->cd();
            }   
            ExecAnalysis();
            Terminate();
            return;
         } 
         // Run tree-based analysis via AliAnalysisSelector  
//         gROOT->ProcessLine(".L $ALICE_ROOT/ANALYSIS/AliAnalysisSelector.cxx+");
         cout << "===== RUNNING LOCAL ANALYSIS " << GetName() << " ON TREE " << tree->GetName() << endl;
         sprintf(line, "AliAnalysisSelector *selector = new AliAnalysisSelector((AliAnalysisManager*)0x%lx);",(ULong_t)this);
         gROOT->ProcessLine(line);
         sprintf(line, "((TTree*)0x%lx)->Process(selector);",(ULong_t)tree);
         gROOT->ProcessLine(line);
         break;
      case kProofAnalysis:
         if (!gROOT->GetListOfProofs() || !gROOT->GetListOfProofs()->GetEntries()) {
            printf("StartAnalysis: no PROOF!!!\n");
            return;
         }   
         sprintf(line, "gProof->AddInput((TObject*)0x%lx);", (ULong_t)this);
         gROOT->ProcessLine(line);
         if (chain) {
            chain->SetProof();
            cout << "===== RUNNING PROOF ANALYSIS " << GetName() << " ON CHAIN " << chain->GetName() << endl;
            chain->Process("AliAnalysisSelector");
         } else {
            printf("StartAnalysis: no chain\n");
            return;
         }      
         break;
      case kGridAnalysis:
         Warning("StartAnalysis", "GRID analysis mode not implemented. Running local.");
   }   
}   

//______________________________________________________________________________
void AliAnalysisManager::ExecAnalysis(Option_t *option)
{
// Execute analysis.
   if (!fInitOK) {
     Error("ExecAnalysis", "Analysis manager was not initialized !");
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
	      Error("ExecAnalysis","Cannot execute analysis in TSelector mode without at least one top container");
         return;
      }   
      cont->SetData(fTree); // This will notify all consumers
      Long64_t entry = fTree->GetTree()->GetReadEntry();
      
//
//    Call BeginEvent() for optional input/output and MC services 
      if (fInputEventHandler)   fInputEventHandler  ->BeginEvent(entry);
      if (fOutputEventHandler)  fOutputEventHandler ->BeginEvent(entry);
      if (fMCtruthEventHandler) fMCtruthEventHandler->BeginEvent(entry);
//
//    Execute the tasks
      TIter next1(cont->GetConsumers());
      while ((task=(AliAnalysisTask*)next1())) {
         if (fDebug >1) {
            cout << "    Executing task " << task->GetName() << endl;
         }   
	 
         task->ExecuteTask(option);
      }
//
//    Call FinishEvent() for optional output and MC services 
      if (fInputEventHandler)   fInputEventHandler  ->FinishEvent();
      if (fOutputEventHandler)  fOutputEventHandler ->FinishEvent();
      if (fMCtruthEventHandler) fMCtruthEventHandler->FinishEvent();
//
      return;
   }   
   // The event loop is not controlled by TSelector   
//
//  Call BeginEvent() for optional input/output and MC services 
   if (fInputEventHandler)   fInputEventHandler  ->BeginEvent(-1);
   if (fOutputEventHandler)  fOutputEventHandler ->BeginEvent(-1);
   if (fMCtruthEventHandler) fMCtruthEventHandler->BeginEvent(-1);
   TIter next2(fTopTasks);
   while ((task=(AliAnalysisTask*)next2())) {
      task->SetActive(kTRUE);
      if (fDebug > 1) {
         cout << "    Executing task " << task->GetName() << endl;
      }   
      task->ExecuteTask(option);
   }   
//
// Call FinishEvent() for optional output and MC services 
   if (fInputEventHandler)   fInputEventHandler  ->FinishEvent();
   if (fOutputEventHandler)  fOutputEventHandler ->FinishEvent();
   if (fMCtruthEventHandler) fMCtruthEventHandler->FinishEvent();
}

//______________________________________________________________________________
void AliAnalysisManager::FinishAnalysis()
{
// Finish analysis.
}
