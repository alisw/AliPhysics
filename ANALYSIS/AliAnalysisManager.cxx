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
#include <TKey.h>
#include <TMethodCall.h>
#include <TChain.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TCanvas.h>

#include "AliAnalysisSelector.h"
#include "AliAnalysisGrid.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisDataSlot.h"
#include "AliVEventHandler.h"
#include "AliVEventPool.h"
#include "AliSysInfo.h"
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
                    fEventPool(NULL),
                    fCurrentEntry(-1),
                    fNSysInfo(0),
                    fMode(kLocalAnalysis),
                    fInitOK(kFALSE),
                    fDebug(0),
                    fSpecialOutputLocation(""), 
                    fTasks(NULL),
                    fTopTasks(NULL),
                    fZombies(NULL),
                    fContainers(NULL),
                    fInputs(NULL),
                    fOutputs(NULL),
                    fCommonInput(NULL),
                    fCommonOutput(NULL),
                    fSelector(NULL),
                    fGridHandler(NULL),
                    fExtraFiles("")
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
		    fEventPool(NULL),
                    fCurrentEntry(-1),
                    fNSysInfo(0),
                    fMode(other.fMode),
                    fInitOK(other.fInitOK),
                    fDebug(other.fDebug),
                    fSpecialOutputLocation(""), 
                    fTasks(NULL),
                    fTopTasks(NULL),
                    fZombies(NULL),
                    fContainers(NULL),
                    fInputs(NULL),
                    fOutputs(NULL),
                    fCommonInput(NULL),
                    fCommonOutput(NULL),
                    fSelector(NULL),
                    fGridHandler(NULL),
                    fExtraFiles()
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
      fEventPool           = other.fEventPool;
      fTree       = NULL;
      fCurrentEntry = -1;
      fNSysInfo   = other.fNSysInfo;
      fMode       = other.fMode;
      fInitOK     = other.fInitOK;
      fDebug      = other.fDebug;
      fTasks      = new TObjArray(*other.fTasks);
      fTopTasks   = new TObjArray(*other.fTopTasks);
      fZombies    = new TObjArray(*other.fZombies);
      fContainers = new TObjArray(*other.fContainers);
      fInputs     = new TObjArray(*other.fInputs);
      fOutputs    = new TObjArray(*other.fOutputs);
      fCommonInput = NULL;
      fCommonOutput = NULL;
      fSelector   = NULL;
      fGridHandler = NULL;
      fExtraFiles = other.fExtraFiles;
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
   if (fGridHandler) delete fGridHandler;
   if (fgAnalysisManager==this) fgAnalysisManager = NULL;
}

//______________________________________________________________________________
Int_t AliAnalysisManager::GetEntry(Long64_t entry, Int_t getall)
{
// Read one entry of the tree or a whole branch.
   if (fDebug > 0) printf("== AliAnalysisManager::GetEntry(%lld)\n", entry);
   fCurrentEntry = entry;
   return fTree ? fTree->GetTree()->GetEntry(entry, getall) : 0;
}
   
//______________________________________________________________________________
Bool_t AliAnalysisManager::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses of the tree
  // will be set. It is normaly not necessary to make changes to the
  // generated code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running with PROOF.
   Bool_t init = kFALSE;
   if (!tree) return kFALSE; // Should not happen - protected in selector caller
   if (fDebug > 0) {
      printf("->AliAnalysisManager::Init(%s)\n", tree->GetName());
   }
   // Call InitTree of EventHandler
   if (fOutputEventHandler) {
      if (fMode == kProofAnalysis) {
         init = fOutputEventHandler->Init(0x0, "proof");
      } else {
         init = fOutputEventHandler->Init(0x0, "local");
      }
      if (!init) {
         Error("Init", "Output event handler failed to initialize");
         return kFALSE;
      }         
   }
   
   if (fInputEventHandler) {
      if (fMode == kProofAnalysis) {
         init = fInputEventHandler->Init(tree, "proof");
      } else {
         init = fInputEventHandler->Init(tree, "local");
      }
      if (!init) {
         Error("Init", "Input event handler failed to initialize tree"); 
         return kFALSE;
      }         
   } else {
      // If no input event handler we need to get the tree once
      // for the chain
      if(!tree->GetTree()) {
         Long64_t readEntry = tree->LoadTree(0);
         if (readEntry == -2) {
            Error("Init", "Input tree has no entry. Aborting");
            return kFALSE;
         }
      }   
   }

   if (fMCtruthEventHandler) {
      if (fMode == kProofAnalysis) {
         init = fMCtruthEventHandler->Init(0x0, "proof");
      } else {
         init = fMCtruthEventHandler->Init(0x0, "local");
      }
      if (!init) {
         Error("Init", "MC event handler failed to initialize"); 
         return kFALSE;
      }         
   }

   if (!fInitOK) InitAnalysis();
   if (!fInitOK) return kFALSE;
   fTree = tree;
   AliAnalysisDataContainer *top = fCommonInput;
   if (!top) top = (AliAnalysisDataContainer*)fInputs->At(0);
   if (!top) {
      Error("Init","No top input container !");
      return kFALSE;
   }
   top->SetData(tree);
   if (fDebug > 0) {
      printf("<-AliAnalysisManager::Init(%s)\n", tree->GetName());
   }
   return kTRUE;
}

//______________________________________________________________________________
void AliAnalysisManager::SlaveBegin(TTree *tree)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).
   if (fDebug > 0) printf("->AliAnalysisManager::SlaveBegin()\n");
   static Bool_t isCalled = kFALSE;
   Bool_t init = kFALSE;
   Bool_t initOK = kTRUE;
   TString msg;
   TDirectory *curdir = gDirectory;
   // Call SlaveBegin only once in case of mixing
   if (isCalled && fMode==kMixingAnalysis) return;
   // Call Init of EventHandler
   if (fOutputEventHandler) {
      if (fMode == kProofAnalysis) {
         // Merging AOD's in PROOF via TProofOutputFile
         if (fDebug > 1) printf("   Initializing AOD output file %s...\n", fOutputEventHandler->GetOutputFileName());
         init = fOutputEventHandler->Init("proof");
         if (!init) msg = "Failed to initialize output handler on worker";
      } else {
         init = fOutputEventHandler->Init("local");
         if (!init) msg = "Failed to initialize output handler";
      }
      initOK &= init;
      if (!fSelector) Error("SlaveBegin", "Selector not set");
      else if (!init) {fSelector->Abort(msg); fSelector->SetStatus(-1);}
   }

   if (fInputEventHandler) {
      fInputEventHandler->SetInputTree(tree);
      if (fMode == kProofAnalysis) {
         init = fInputEventHandler->Init("proof");
         if (!init) msg = "Failed to initialize input handler on worker";
      } else {
         init = fInputEventHandler->Init("local");
         if (!init) msg = "Failed to initialize input handler";
      }
      initOK &= init;
      if (!fSelector) Error("SlaveBegin", "Selector not set");      
      else if (!init) {fSelector->Abort(msg); fSelector->SetStatus(-1);}
   }

   if (fMCtruthEventHandler) {
      if (fMode == kProofAnalysis) {
         init = fMCtruthEventHandler->Init("proof");
         if (!init) msg = "Failed to initialize MC handler on worker";
      } else {
         init = fMCtruthEventHandler->Init("local");
         if (!init) msg = "Failed to initialize MC handler";
      }
      initOK &= init;
      if (!fSelector) Error("SlaveBegin", "Selector not set");      
      else if (!init) {fSelector->Abort(msg); fSelector->SetStatus(-1);}
   }
   if (curdir) curdir->cd();
   isCalled = kTRUE;
   if (!initOK) return;   
   TIter next(fTasks);
   AliAnalysisTask *task;
   // Call CreateOutputObjects for all tasks
   while ((task=(AliAnalysisTask*)next())) {
      curdir = gDirectory;
      task->CreateOutputObjects();
      if (curdir) curdir->cd();
   }
   if (fDebug > 0) printf("<-AliAnalysisManager::SlaveBegin()\n");
}

//______________________________________________________________________________
Bool_t AliAnalysisManager::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.
   if (!fTree) return kFALSE;

   TFile *curfile = fTree->GetCurrentFile();
   if (!curfile) {
      Error("Notify","No current file");
      return kFALSE;
   }   
   
   if (fDebug > 0) printf("->AliAnalysisManager::Notify() file: %s\n", curfile->GetName());
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
   if (fDebug > 0) printf("<-AliAnalysisManager::Notify()\n");
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
   if (fDebug > 0) printf("->AliAnalysisManager::Process(%lld)\n", entry);

   if (fInputEventHandler)   fInputEventHandler  ->BeginEvent(entry);
   if (fOutputEventHandler)  fOutputEventHandler ->BeginEvent(entry);
   if (fMCtruthEventHandler) fMCtruthEventHandler->BeginEvent(entry);
   
   GetEntry(entry);
   ExecAnalysis();
   if (fDebug > 0) printf("<-AliAnalysisManager::Process()\n");
   return kTRUE;
}

//______________________________________________________________________________
void AliAnalysisManager::PackOutput(TList *target)
{
  // Pack all output data containers in the output list. Called at SlaveTerminate
  // stage in PROOF case for each slave.
   if (fDebug > 0) printf("->AliAnalysisManager::PackOutput()\n");
   if (!target) {
      Error("PackOutput", "No target. Aborting.");
      return;
   }
   if (fInputEventHandler)   fInputEventHandler  ->Terminate();
   if (fOutputEventHandler)  fOutputEventHandler ->Terminate();
   if (fMCtruthEventHandler) fMCtruthEventHandler->Terminate();

   // Call FinishTaskOutput() for each event loop task (not called for 
   // post-event loop tasks - use Terminate() fo those)
   TIter nexttask(fTasks);
   AliAnalysisTask *task;
   while ((task=(AliAnalysisTask*)nexttask())) {
      if (!task->IsPostEventLoop()) {
         if (fDebug > 0) printf("->FinishTaskOutput: task %s\n", task->GetName());
         task->FinishTaskOutput();
         if (fDebug > 0) printf("<-FinishTaskOutput: task %s\n", task->GetName());
      }
   }      
   
   if (fMode == kProofAnalysis) {
      TIter next(fOutputs);
      AliAnalysisDataContainer *output;
      Bool_t isManagedByHandler = kFALSE;
      while ((output=(AliAnalysisDataContainer*)next())) {
         // Do not consider outputs of post event loop tasks
         isManagedByHandler = kFALSE;
         if (output->GetProducer()->IsPostEventLoop()) continue;
         const char *filename = output->GetFileName();
         if (!(strcmp(filename, "default")) && fOutputEventHandler) {
            isManagedByHandler = kTRUE;
            filename = fOutputEventHandler->GetOutputFileName();
         }
         // Check if data was posted to this container. If not, issue an error.
         if (!output->GetData() && !isManagedByHandler) {
            Error("PackOutput", "No data for output container %s. Forgot to PostData ?\n", output->GetName());
            continue;
         }   
         if (!output->IsSpecialOutput()) {
            // Normal outputs
            if (strlen(filename) && !isManagedByHandler) {
               // File resident outputs
               TFile *file = output->GetFile();
               // Backup current folder
               TDirectory *opwd = gDirectory;
               // Create file if not existing and register to container.
               if (file) file->cd();
               else      file = new TFile(filename, "RECREATE"); 
               if (file->IsZombie()) {
                  Fatal("PackOutput", "Could not recreate file %s\n", filename);
                  return;
               }   
               output->SetFile(file);
               // Clear file list to release object ownership to user.
               file->Clear();
               // Save data to file, then close.
               if (output->GetData()->InheritsFrom(TCollection::Class())) {
                  // If data is a collection, we set the name of the collection 
                  // as the one of the container and we save as a single key.
                  TCollection *coll = (TCollection*)output->GetData();
                  coll->SetName(output->GetName());
                  coll->Write(output->GetName(), TObject::kSingleKey);
               } else {
                  if (output->GetData()->InheritsFrom(TTree::Class())) {
                     TTree *tree = (TTree*)output->GetData();
                     tree->SetDirectory(file);
                     tree->AutoSave();
                  } else {
                     output->GetData()->Write();
                  }   
               }      
               if (fDebug > 1) printf("PackOutput %s: memory merge, file resident output\n", output->GetName());
               if (fDebug > 2) {
                  printf("   file %s listing content:\n", filename);
                  file->ls();
               }   
               file->Close();
               // Restore current directory
               if (opwd) opwd->cd();
            } else {
               // Memory-resident outputs   
               if (fDebug > 1) printf("PackOutput %s: memory merge memory resident output\n", filename);
            }   
            AliAnalysisDataWrapper *wrap = 0;
            if (isManagedByHandler) {
               wrap = new AliAnalysisDataWrapper(fOutputEventHandler->GetTree());
               wrap->SetName(output->GetName());
            }   
            else                    wrap =output->ExportData();
            // Output wrappers must NOT delete data after merging - the user owns them
            wrap->SetDeleteData(kFALSE);
            target->Add(wrap);
         } else {
         // Special outputs
            TDirectory *opwd = gDirectory;
            TFile *file = output->GetFile();
            if (fDebug > 1 && file) printf("PackOutput %s: file merge, special output\n", output->GetName());
            if (isManagedByHandler) {
               // Terminate IO for files managed by the output handler
               if (file) file->Write();
               if (file && fDebug > 2) {
                  printf("   handled file %s listing content:\n", file->GetName());
                  file->ls();
               }   
               fOutputEventHandler->TerminateIO();
               continue;
            }   
            
            if (!file) {
               AliAnalysisTask *producer = output->GetProducer();
               Error("PackOutput", 
                     "File %s for special container %s was NOT opened in %s::CreateOutputObjects !!!",
                     output->GetFileName(), output->GetName(), producer->ClassName());
               continue;
            }   
            file->cd();
            // Release object ownership to users after writing data to file
            if (output->GetData()->InheritsFrom(TCollection::Class())) {
               // If data is a collection, we set the name of the collection 
               // as the one of the container and we save as a single key.
               TCollection *coll = (TCollection*)output->GetData();
               coll->SetName(output->GetName());
               coll->Write(output->GetName(), TObject::kSingleKey);
            } else {
               if (output->GetData()->InheritsFrom(TTree::Class())) {
                  TTree *tree = (TTree*)output->GetData();
                  tree->SetDirectory(file);
                  tree->AutoSave();
               } else {
                  output->GetData()->Write();
               }   
            }      
            file->Clear();
            if (fDebug > 2) {
               printf("   file %s listing content:\n", output->GetFileName());
               file->ls();
            }
            TString outFilename = file->GetName();
            file->Close();
            // Restore current directory
            if (opwd) opwd->cd();
            // Check if a special output location was provided or the output files have to be merged
            if (strlen(fSpecialOutputLocation.Data())) {
               TString remote = fSpecialOutputLocation;
               remote += "/";
               Int_t gid = gROOT->ProcessLine("gProofServ->GetGroupId();");
               remote += Form("%s_%d_", gSystem->HostName(), gid);
               remote += output->GetFileName();
               TFile::Cp ( outFilename.Data(), remote.Data() );
            } else {
            // No special location specified-> use TProofOutputFile as merging utility
            // The file at this output slot must be opened in CreateOutputObjects
               if (fDebug > 1) printf("   File %s to be merged...\n", output->GetFileName());
            }
         }      
      }
   } 
   if (fDebug > 0) printf("<-AliAnalysisManager::PackOutput: output list contains %d containers\n", target->GetSize());
}

//______________________________________________________________________________
void AliAnalysisManager::ImportWrappers(TList *source)
{
// Import data in output containers from wrappers coming in source.
   if (fDebug > 0) printf("->AliAnalysisManager::ImportWrappers()\n");
   TIter next(fOutputs);
   AliAnalysisDataContainer *cont;
   AliAnalysisDataWrapper   *wrap;
   Int_t icont = 0;
   Bool_t inGrid = (fMode == kGridAnalysis)?kTRUE:kFALSE;
   while ((cont=(AliAnalysisDataContainer*)next())) {
      wrap = 0;
      if (cont->GetProducer()->IsPostEventLoop() && !inGrid) continue;
      const char *filename = cont->GetFileName();
      Bool_t isManagedByHandler = kFALSE;
      if (!(strcmp(filename, "default")) && fOutputEventHandler) {
         isManagedByHandler = kTRUE;
         filename = fOutputEventHandler->GetOutputFileName();
      }
      if (cont->IsSpecialOutput() || inGrid) {
         if (strlen(fSpecialOutputLocation.Data()) && !isManagedByHandler) continue;
         // Copy merged file from PROOF scratch space. 
         // In case of grid the files are already in the current directory.
         if (!inGrid) {
            if (isManagedByHandler && fExtraFiles.Length()) {
               // Copy extra registered dAOD files.
               TObjArray *arr = fExtraFiles.Tokenize(" ");
               TObjString *os;
               TIter nextfilename(arr);
               while ((os=(TObjString*)nextfilename())) GetFileFromWrapper(os->GetString(), source);
               delete arr;
            }
            if (!GetFileFromWrapper(filename, source)) continue;
         }   
         // Normally we should connect data from the copied file to the
         // corresponding output container, but it is not obvious how to do this
         // automatically if several objects in file...
         TFile *f = TFile::Open(filename, "READ");
         if (!f) {
            Error("ImportWrappers", "Cannot open file %s in read-only mode", filename);
            continue;
         }   
         TObject *obj = 0;
         // Try to fetch first a list object having the container name.
         obj = f->Get(cont->GetName());
         if (!obj) {
         // Fetch first object from file having the container type.
            TIter nextkey(f->GetListOfKeys());
            TKey *key;
            while ((key=(TKey*)nextkey())) {
               obj = f->Get(key->GetName());
               if (obj && obj->IsA()->InheritsFrom(cont->GetType())) break;
            }                     
         }
         if (!obj) {
            Error("ImportWrappers", "Could not find object for container %s in file %s", cont->GetName(), filename);
            continue;
         }  
         wrap = new AliAnalysisDataWrapper(obj);
         wrap->SetDeleteData(kFALSE);
      }   
      if (!wrap) wrap = (AliAnalysisDataWrapper*)source->FindObject(cont->GetName());
      if (!wrap) {
         Error("ImportWrappers","Container %s not found in analysis output !", cont->GetName());
         continue;
      }
      icont++;
      if (fDebug > 1) {
         printf("   Importing data for container %s", cont->GetName());
         if (strlen(filename)) printf("    -> file %s\n", filename);
         else printf("\n");
      }   
      cont->ImportData(wrap);
   }         
   if (fDebug > 0) printf("<-AliAnalysisManager::ImportWrappers(): %d containers imported\n", icont);
}

//______________________________________________________________________________
void AliAnalysisManager::UnpackOutput(TList *source)
{
  // Called by AliAnalysisSelector::Terminate only on the client.
   if (fDebug > 0) printf("->AliAnalysisManager::UnpackOutput()\n");
   if (!source) {
      Error("UnpackOutput", "No target. Aborting.");
      return;
   }
   if (fDebug > 1) printf("   Source list contains %d containers\n", source->GetSize());

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
               if (fDebug > 0) printf("== Executing post event loop task %s\n", task->GetName());
               task->ExecuteTask();
            }   
         }
      }   
   }
   if (fDebug > 0) printf("<-AliAnalysisManager::UnpackOutput()\n");
}

//______________________________________________________________________________
void AliAnalysisManager::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically.
   if (fDebug > 0) printf("->AliAnalysisManager::Terminate()\n");
   AliAnalysisTask *task;
   TIter next(fTasks);
   // Call Terminate() for tasks
   while ((task=(AliAnalysisTask*)next())) task->Terminate();
   //
   TIter next1(fOutputs);
   AliAnalysisDataContainer *output;
   while ((output=(AliAnalysisDataContainer*)next1())) {
      // Special outputs or grid files have the files already closed and written.
      if (fMode == kGridAnalysis) continue;
      if (output->IsSpecialOutput()&&(fMode == kProofAnalysis)) continue;
      const char *filename = output->GetFileName();
      if (!(strcmp(filename, "default"))) {
         if (fOutputEventHandler) filename = fOutputEventHandler->GetOutputFileName();
         TFile *aodfile = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
         if (aodfile) {
            if (fDebug > 1) printf("Writing output handler file: %s\n", filename);
            aodfile->Write();
            continue;
         }   
      }      
      if (!strlen(filename)) continue;
      if (!output->GetData()) continue;
      TFile *file = output->GetFile();
      TDirectory *opwd = gDirectory;
      file = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
      if (!file) file = new TFile(filename, "RECREATE");
      if (file->IsZombie()) continue;
      output->SetFile(file);
      file->cd();
      if (fDebug > 1) printf("   writing output data %s to file %s\n", output->GetData()->GetName(), file->GetName());
      if (output->GetData()->InheritsFrom(TCollection::Class())) {
      // If data is a collection, we set the name of the collection 
      // as the one of the container and we save as a single key.
         TCollection *coll = (TCollection*)output->GetData();
         coll->SetName(output->GetName());
         coll->Write(output->GetName(), TObject::kSingleKey);
      } else {
         if (output->GetData()->InheritsFrom(TTree::Class())) {
            TTree *tree = (TTree*)output->GetData();
            tree->SetDirectory(file);
            tree->AutoSave();
         } else {
            output->GetData()->Write();
         }   
      }      
      if (opwd) opwd->cd();
   }   
   next1.Reset();
   while ((output=(AliAnalysisDataContainer*)next1())) {
      // Close all files at output
      TDirectory *opwd = gDirectory;
      if (output->GetFile()) output->GetFile()->Close();
      if (opwd) opwd->cd();
   }   

   if (fInputEventHandler)   fInputEventHandler  ->TerminateIO();
   if (fOutputEventHandler)  fOutputEventHandler ->TerminateIO();
   if (fMCtruthEventHandler) fMCtruthEventHandler->TerminateIO();

   Bool_t getsysInfo = ((fNSysInfo>0) && (fMode==kLocalAnalysis))?kTRUE:kFALSE;
   if (getsysInfo) {
      TDirectory *cdir = gDirectory;
      TFile f("syswatch.root", "RECREATE");
      if (!f.IsZombie()) {
         TTree *tree = AliSysInfo::MakeTree("syswatch.log");
         tree->SetMarkerStyle(kCircle);
         tree->SetMarkerColor(kBlue);
         tree->SetMarkerSize(0.5);
         if (!gROOT->IsBatch()) {
            tree->SetAlias("event", "id0");
            tree->SetAlias("memUSED", "mi.fMemUsed");
            new TCanvas("SysInfo","SysInfo",10,10,800,600);
            tree->Draw("memUSED:event","","", 1234567890, 0);
         }   
         tree->Write();
         f.Close();
         delete tree;
      }
      if (cdir) cdir->cd();
   }      
   if (fDebug > 0) printf("<-AliAnalysisManager::Terminate()\n");
}

//______________________________________________________________________________
void AliAnalysisManager::AddTask(AliAnalysisTask *task)
{
// Adds a user task to the global list of tasks.
   if (fTasks->FindObject(task)) {
      Warning("AddTask", "Task %s: the same object already added to the analysis manager. Not adding.", task->GetName());
      return;
   }   
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
         if (filename && strlen(filename)) {
            cont->SetFileName(filename);
            cont->SetDataOwned(kFALSE);  // data owned by the file
         }   
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
   if (!task) {
      Error("ConnectInput", "Task pointer is NULL");
      return kFALSE;
   }   
   if (!fTasks->FindObject(task)) {
      AddTask(task);
      Info("ConnectInput", "Task %s was not registered. Now owned by analysis manager", task->GetName());
   } 
   Bool_t connected = task->ConnectInput(islot, cont);
   return connected;
}   

//______________________________________________________________________________
Bool_t AliAnalysisManager::ConnectOutput(AliAnalysisTask *task, Int_t islot,
                                        AliAnalysisDataContainer *cont)
{
// Connect output of an existing task to a data container.
   if (!task) {
      Error("ConnectOutput", "Task pointer is NULL");
      return kFALSE;
   }   
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
   // Check if all special output containers have a file name provided
   TIter nextout(fOutputs);
   while ((cont=(AliAnalysisDataContainer*)nextout())) {
      if (cont->IsSpecialOutput() && !strlen(cont->GetFileName())) {
         Error("InitAnalysis", "Wrong container %s : a file name MUST be provided for special outputs", cont->GetName());
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
   if (!fInitOK) {
      Info("PrintStatus", "Analysis manager %s not initialized : call InitAnalysis() first", GetName());
      return;
   }   
   Bool_t getsysInfo = ((fNSysInfo>0) && (fMode==kLocalAnalysis))?kTRUE:kFALSE;
   if (getsysInfo)
      Info("PrintStatus", "System information will be collected each %lld events", fNSysInfo);
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
void AliAnalysisManager::StartAnalysis(const char *type, TTree *tree, Long64_t nentries, Long64_t firstentry)
{
// Start analysis for this manager. Analysis task can be: LOCAL, PROOF, GRID or
// MIX. Process nentries starting from firstentry
   if (!fInitOK) {
      Error("StartAnalysis","Analysis manager was not initialized !");
      return;
   }
   if (fDebug > 0) printf("StartAnalysis %s\n",GetName());
   TString anaType = type;
   anaType.ToLower();
   fMode = kLocalAnalysis;
   Bool_t runlocalinit = kTRUE;
   if (anaType.Contains("file")) runlocalinit = kFALSE;
   if (anaType.Contains("proof"))     fMode = kProofAnalysis;
   else if (anaType.Contains("grid")) fMode = kGridAnalysis;
   else if (anaType.Contains("mix"))  fMode = kMixingAnalysis;

   if (fMode == kGridAnalysis) {
      if (!fGridHandler) {
         Error("StartAnalysis", "Cannot start grid analysis without a grid handler.");
         Info("===", "Add an AliAnalysisAlien object as plugin for this manager and configure it.");
         return;
      }
      // Write analysis manager in the analysis file
      cout << "===== RUNNING GRID ANALYSIS: " << GetName() << endl;
      // run local task configuration
      TIter nextTask(fTasks);
      AliAnalysisTask *task;
      while ((task=(AliAnalysisTask*)nextTask())) {
         task->LocalInit();
      }
      fGridHandler->StartAnalysis(nentries, firstentry);

      // Terminate grid analysis
      if (fSelector && fSelector->GetStatus() == -1) return;
      if (fGridHandler->GetRunMode() == AliAnalysisGrid::kOffline) return;
      cout << "===== MERGING OUTPUTS REGISTERED BY YOUR ANALYSIS JOB: " << GetName() << endl;
      if (!fGridHandler->MergeOutputs()) {
         // Return if outputs could not be merged or if it alien handler
         // was configured for offline mode or local testing.
         return;
      }
      ImportWrappers(NULL);
      Terminate();
      return;
   }
   char line[256];
   SetEventLoop(kFALSE);
   // Enable event loop mode if a tree was provided
   if (tree || fMode==kMixingAnalysis) SetEventLoop(kTRUE);

   TChain *chain = 0;
   TString ttype = "TTree";
   if (tree && tree->IsA() == TChain::Class()) {
      chain = (TChain*)tree;
      if (!chain || !chain->GetListOfFiles()->First()) {
         Error("StartAnalysis", "Cannot process null or empty chain...");
         return;
      }   
      ttype = "TChain";
   }   

   // Initialize locally all tasks (happens for all modes)
   TIter next(fTasks);
   AliAnalysisTask *task;
   if (runlocalinit) {
      while ((task=(AliAnalysisTask*)next())) {
         task->LocalInit();
      }
   }   
   
   switch (fMode) {
      case kLocalAnalysis:
         if (!tree) {
            TIter nextT(fTasks);
            // Call CreateOutputObjects for all tasks
            while ((task=(AliAnalysisTask*)nextT())) {
               TDirectory *curdir = gDirectory;
               task->CreateOutputObjects();
               if (curdir) curdir->cd();
            }   
            ExecAnalysis();
            Terminate();
            return;
         } 
         // Run tree-based analysis via AliAnalysisSelector  
         cout << "===== RUNNING LOCAL ANALYSIS " << GetName() << " ON TREE " << tree->GetName() << endl;
         fSelector = new AliAnalysisSelector(this);
         tree->Process(fSelector, "", nentries, firstentry);
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
            chain->Process("AliAnalysisSelector", "", nentries, firstentry);
         } else {
            printf("StartAnalysis: no chain\n");
            return;
         }      
         break;
      case kGridAnalysis:
         Warning("StartAnalysis", "GRID analysis mode not implemented. Running local.");
         break;
      case kMixingAnalysis:   
         // Run event mixing analysis
         if (!fEventPool) {
            Error("StartAnalysis", "Cannot run event mixing without event pool");
            return;
         }
         cout << "===== RUNNING EVENT MIXING ANALYSIS " << GetName() << endl;
         fSelector = new AliAnalysisSelector(this);
         while ((chain=fEventPool->GetNextChain())) {
            next.Reset();
            // Call NotifyBinChange for all tasks
            while ((task=(AliAnalysisTask*)next()))
               if (!task->IsPostEventLoop()) task->NotifyBinChange();
            chain->Process(fSelector);
         }
         PackOutput(fSelector->GetOutputList());
         Terminate();
   }   
}   

//______________________________________________________________________________
void AliAnalysisManager::StartAnalysis(const char *type, const char *dataset, Long64_t nentries, Long64_t firstentry)
{
// Start analysis for this manager on a given dataset. Analysis task can be: 
// LOCAL, PROOF or GRID. Process nentries starting from firstentry.
   if (!fInitOK) {
      Error("StartAnalysis","Analysis manager was not initialized !");
      return;
   }
   if (fDebug > 0) printf("StartAnalysis %s\n",GetName());
   TString anaType = type;
   anaType.ToLower();
   if (!anaType.Contains("proof")) {
      Error("StartAnalysis", "Cannot process datasets in %s mode. Try PROOF.", type);
      return;
   }   
   fMode = kProofAnalysis;
   char line[256];
   SetEventLoop(kTRUE);
   // Set the dataset flag
   TObject::SetBit(kUseDataSet);
   fTree = 0;

   // Initialize locally all tasks
   TIter next(fTasks);
   AliAnalysisTask *task;
   while ((task=(AliAnalysisTask*)next())) {
      task->LocalInit();
   }
   
   if (!gROOT->GetListOfProofs() || !gROOT->GetListOfProofs()->GetEntries()) {
      printf("StartAnalysis: no PROOF!!!\n");
      return;
   }   
   sprintf(line, "gProof->AddInput((TObject*)0x%lx);", (ULong_t)this);
   gROOT->ProcessLine(line);
   sprintf(line, "gProof->GetDataSet(\"%s\");", dataset);
   if (!gROOT->ProcessLine(line)) {
      Error("StartAnalysis", "Dataset %s not found", dataset);
      return;
   }   
   sprintf(line, "gProof->Process(\"%s\", \"AliAnalysisSelector\", \"\", %lld, %lld);",
           dataset, nentries, firstentry);
   cout << "===== RUNNING PROOF ANALYSIS " << GetName() << " ON DATASET " << dataset << endl;
   gROOT->ProcessLine(line);
}   

//______________________________________________________________________________
TFile *AliAnalysisManager::OpenProofFile(const char *filename, const char *option)
{
// Opens a special output file used in PROOF.
   char line[256];
   if (fMode!=kProofAnalysis || !fSelector) {
      Error("OpenProofFile","Cannot open PROOF file %s",filename);
      return NULL;
   }   
   sprintf(line, "TProofOutputFile *pf = new TProofOutputFile(\"%s\");", filename);
   if (fDebug > 1) printf("=== %s\n", line);
   gROOT->ProcessLine(line);
   sprintf(line, "pf->OpenFile(\"%s\");", option);
   gROOT->ProcessLine(line);
   if (fDebug > 1) {
      gROOT->ProcessLine("pf->Print()");
      printf(" == proof file name: %s\n", gFile->GetName());
   }   
   sprintf(line, "((TList*)0x%lx)->Add(pf);",(ULong_t)fSelector->GetOutputList());
   if (fDebug > 1) printf("=== %s\n", line);
   gROOT->ProcessLine(line);
   return gFile;
}   

//______________________________________________________________________________
void AliAnalysisManager::ExecAnalysis(Option_t *option)
{
// Execute analysis.
   static Long64_t ncalls = 0;
   Bool_t getsysInfo = ((fNSysInfo>0) && (fMode==kLocalAnalysis))?kTRUE:kFALSE;
   if (getsysInfo && ncalls==0) AliSysInfo::AddStamp("Start", (Int_t)ncalls);
   ncalls++;
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
      AliAnalysisDataContainer *cont = fCommonInput;
      if (!cont) cont = (AliAnalysisDataContainer*)fInputs->At(0);
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
//      TIter next1(cont->GetConsumers());
      TIter next1(fTopTasks);
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
      // Gather system information if requested
      if (getsysInfo && ((ncalls%fNSysInfo)==0)) 
         AliSysInfo::AddStamp(Form("Event#%lld",ncalls),(Int_t)ncalls);
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

//______________________________________________________________________________
void AliAnalysisManager::SetInputEventHandler(AliVEventHandler*  handler)
{
// Set the input event handler and create a container for it.
   fInputEventHandler   = handler;
   fCommonInput = CreateContainer("cAUTO_INPUT", TChain::Class(), AliAnalysisManager::kInputContainer);
   Warning("SetInputEventHandler", " An automatic input container for the input chain was created.\nPlease use: mgr->GetCommonInputContainer() to access it.");
}

//______________________________________________________________________________
void AliAnalysisManager::SetOutputEventHandler(AliVEventHandler*  handler)
{
// Set the input event handler and create a container for it.
   fOutputEventHandler   = handler;
   fCommonOutput = CreateContainer("cAUTO_OUTPUT", TTree::Class(), AliAnalysisManager::kOutputContainer, "default");
   fCommonOutput->SetSpecialOutput();
   Warning("SetOutputEventHandler", " An automatic output container for the output tree was created.\nPlease use: mgr->GetCommonOutputContainer() to access it.");
}

//______________________________________________________________________________
void AliAnalysisManager::RegisterExtraFile(const char *fname)
{
// This method is used externally to register output files which are not
// connected to any output container, so that the manager can properly register,
// retrieve or merge them when running in distributed mode. The file names are
// separated by blancs. The method has to be called in MyAnalysisTask::LocalInit().
   if (fExtraFiles.Length()) fExtraFiles += " ";
   fExtraFiles += fname;
}

//______________________________________________________________________________
Bool_t AliAnalysisManager::GetFileFromWrapper(const char *filename, TList *source)
{
// Copy a file from the location specified ina the wrapper with the same name from the source list.
   char full_path[512];
   char ch_url[512];
   TObject *pof =  source->FindObject(filename);
   if (!pof || !pof->InheritsFrom("TProofOutputFile")) {
      Error("GetFileFromWrapper", "TProofOutputFile object not found in output list for file %s", filename);
      return kFALSE;
   }
   gROOT->ProcessLine(Form("sprintf((char*)0x%lx, \"%%s\", ((TProofOutputFile*)0x%lx)->GetOutputFileName();)", full_path, pof));
   gROOT->ProcessLine(Form("sprintf((char*)0x%lx, \"%%s\", gProof->GetUrl();)", ch_url));
   TString clientUrl(ch_url);
   TString full_path_str(full_path);
   if (clientUrl.Contains("localhost")){
      TObjArray* array = full_path_str.Tokenize ( "//" );
      TObjString *strobj = ( TObjString *)array->At(1);
      TObjArray* arrayPort = strobj->GetString().Tokenize ( ":" );
      TObjString *strobjPort = ( TObjString *) arrayPort->At(1);
      full_path_str.ReplaceAll(strobj->GetString().Data(),"localhost:PORT");
      full_path_str.ReplaceAll(":PORT",Form(":%s",strobjPort->GetString().Data()));
      if (fDebug > 1) Info("GetFileFromWrapper","Using tunnel from %s to %s",full_path_str.Data(),filename);
      delete arrayPort;
      delete array;
   }
   if (fDebug > 1) 
      Info("GetFileFromWrapper","Copying file %s from PROOF scratch space", full_path_str.Data());
   Bool_t gotit = TFile::Cp(full_path_str.Data(), filename); 
   if (!gotit)
      Error("GetFileFromWrapper", "Could not get file %s from proof scratch space", filename);
   return gotit;
}
