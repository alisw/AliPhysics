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
//   AliAnalysisManager - Manager analysis class. Allows creation of several
// analysis tasks and data containers storing their input/output. Allows
// connecting/chaining tasks via shared data containers. Serializes the current
// event for all tasks depending only on initial input data.
//==============================================================================
//
//==============================================================================

#include "AliAnalysisManager.h"

#include <cerrno>
#include <Riostream.h>
#include <TError.h>
#include <TMap.h>
#include <TClass.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1.h>
#include <TMethodCall.h>
#include <TChain.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TStopwatch.h>

#include "AliLog.h"
#include "AliAnalysisSelector.h"
#include "AliAnalysisGrid.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisDataSlot.h"
#include "AliVEventHandler.h"
#include "AliVEventPool.h"
#include "AliSysInfo.h"
#include "AliAnalysisStatistics.h"

ClassImp(AliAnalysisManager)

AliAnalysisManager *AliAnalysisManager::fgAnalysisManager = NULL;
TString AliAnalysisManager::fgCommonFileName = "";
Int_t AliAnalysisManager::fPBUpdateFreq = 1;

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
                    fIsRemote(kFALSE),
                    fDebug(0),
                    fSpecialOutputLocation(""), 
                    fTasks(NULL),
                    fTopTasks(NULL),
                    fZombies(NULL),
                    fContainers(NULL),
                    fInputs(NULL),
                    fOutputs(NULL),
                    fParamCont(NULL),
                    fCommonInput(NULL),
                    fCommonOutput(NULL),
                    fSelector(NULL),
                    fGridHandler(NULL),
                    fExtraFiles(""),
                    fAutoBranchHandling(kTRUE), 
                    fTable(),
                    fRunFromPath(0),
                    fNcalls(0),
                    fMaxEntries(0),
                    fStatisticsMsg(),
                    fRequestedBranches(),
                    fStatistics(0),
                    fGlobals(0)
{
// Default constructor.
   fgAnalysisManager = this;
   fgCommonFileName  = "AnalysisResults.root";
   if (TClass::IsCallingNew() != TClass::kDummyNew) {
     fTasks      = new TObjArray();
     fTopTasks   = new TObjArray();
     fZombies    = new TObjArray();
     fContainers = new TObjArray();
     fInputs     = new TObjArray();
     fOutputs    = new TObjArray();
     fParamCont  = new TObjArray();
     fGlobals    = new TMap();
   }  
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
                    fIsRemote(other.fIsRemote),
                    fDebug(other.fDebug),
                    fSpecialOutputLocation(""), 
                    fTasks(NULL),
                    fTopTasks(NULL),
                    fZombies(NULL),
                    fContainers(NULL),
                    fInputs(NULL),
                    fOutputs(NULL),
                    fParamCont(NULL),
                    fCommonInput(NULL),
                    fCommonOutput(NULL),
                    fSelector(NULL),
                    fGridHandler(NULL),
                    fExtraFiles(),
                    fAutoBranchHandling(other.fAutoBranchHandling), 
                    fTable(),
                    fRunFromPath(0),
                    fNcalls(other.fNcalls),
                    fMaxEntries(other.fMaxEntries),
                    fStatisticsMsg(other.fStatisticsMsg),
                    fRequestedBranches(other.fRequestedBranches),
                    fStatistics(other.fStatistics),
                    fGlobals(other.fGlobals)
{
// Copy constructor.
   fTasks      = new TObjArray(*other.fTasks);
   fTopTasks   = new TObjArray(*other.fTopTasks);
   fZombies    = new TObjArray(*other.fZombies);
   fContainers = new TObjArray(*other.fContainers);
   fInputs     = new TObjArray(*other.fInputs);
   fOutputs    = new TObjArray(*other.fOutputs);
   fParamCont  = new TObjArray(*other.fParamCont);
   fgCommonFileName  = "AnalysisResults.root";
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
      fIsRemote   = other.fIsRemote;
      fDebug      = other.fDebug;
      fTasks      = new TObjArray(*other.fTasks);
      fTopTasks   = new TObjArray(*other.fTopTasks);
      fZombies    = new TObjArray(*other.fZombies);
      fContainers = new TObjArray(*other.fContainers);
      fInputs     = new TObjArray(*other.fInputs);
      fOutputs    = new TObjArray(*other.fOutputs);
      fParamCont  = new TObjArray(*other.fParamCont);
      fCommonInput = NULL;
      fCommonOutput = NULL;
      fSelector   = NULL;
      fGridHandler = NULL;
      fExtraFiles = other.fExtraFiles;
      fgCommonFileName = "AnalysisResults.root";
      fgAnalysisManager = this;
      fAutoBranchHandling = other.fAutoBranchHandling;
      fTable.Clear("nodelete");
      fRunFromPath = other.fRunFromPath;
      fNcalls     = other. fNcalls;
      fMaxEntries = other.fMaxEntries;
      fStatisticsMsg = other.fStatisticsMsg;
      fRequestedBranches = other.fRequestedBranches;
      fStatistics = other.fStatistics;
      fGlobals = new TMap();
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
   if (fParamCont) delete fParamCont;
   if (fGridHandler) delete fGridHandler;
   if (fInputEventHandler) delete fInputEventHandler;
   if (fOutputEventHandler) delete fOutputEventHandler;
   if (fMCtruthEventHandler) delete fMCtruthEventHandler;
   if (fEventPool) delete fEventPool;
   if (fgAnalysisManager==this) fgAnalysisManager = NULL;
   if (fGlobals) {fGlobals->DeleteAll(); delete fGlobals;}
}

//______________________________________________________________________________
Int_t AliAnalysisManager::GetEntry(Long64_t entry, Int_t getall)
{
// Read one entry of the tree or a whole branch.
   fCurrentEntry = entry;
   if (!fAutoBranchHandling)
     return 123456789;
   return fTree ? fTree->GetTree()->GetEntry(entry, getall) : -1;
}

//______________________________________________________________________________
Int_t AliAnalysisManager::GetRunFromAlienPath(const char *path)
{
// Attempt to extract run number from input data path. Works only for paths to
// alice data in alien.
//    sim:  /alice/sim/<production>/run_no/...
//    data: /alice/data/year/period/000run_no/... (ESD or AOD)
   TString type = "unknown";
   TString s(path);
   if (s.Contains("/alice/data")) type = "real";
   else if (s.Contains("/alice/sim")) type = "simulated";
   TString srun;
   Int_t ind1, ind2;
   ind1 = s.Index("/00");
   if (ind1>0) {
      ind2 = s.Index("/",ind1+1);
      if (ind2-ind1>8) srun = s(ind1+1, ind2-ind1-1);
   }   
   if (srun.IsNull()) {
      ind1 = s.Index("/LHC");
      if (ind1>0) {
         ind1 = s.Index("/",ind1+1);
         if (ind1>0) {
            ind2 = s.Index("/",ind1+1);
            if (ind2>0) srun = s(ind1+1, ind2-ind1-1);
         }
      }
   }         
   Int_t run = srun.Atoi();
   if (run>0) printf("=== GetRunFromAlienPath: run %d of %s data ===\n", run, type.Data());
   return run;
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
   if (fDebug > 1) {
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
            Error("Init", "Input tree has no entry. Exiting");
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
   fTable.Rehash(100);
   AliAnalysisDataContainer *top = fCommonInput;
   if (!top) top = (AliAnalysisDataContainer*)fInputs->At(0);
   if (!top) {
      Error("Init","No top input container !");
      return kFALSE;
   }
   top->SetData(tree);
   CheckBranches(kFALSE);
   if (fDebug > 1) {
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
   if (fDebug > 1) printf("->AliAnalysisManager::SlaveBegin()\n");
   if (!CheckTasks()) Fatal("SlaveBegin", "Not all needed libraries were loaded");
   static Bool_t isCalled = kFALSE;
   Bool_t init = kFALSE;
   Bool_t initOK = kTRUE;
   TString msg;
   TDirectory *curdir = gDirectory;
   // Call SlaveBegin only once in case of mixing
   if (isCalled && fMode==kMixingAnalysis) return;
   gROOT->cd();
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
   gROOT->cd();
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
   gROOT->cd();
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
   Bool_t getsysInfo = ((fNSysInfo>0) && (fMode==kLocalAnalysis))?kTRUE:kFALSE;
   Bool_t dirStatus = TH1::AddDirectoryStatus();
   Int_t itask = 0;
   while ((task=(AliAnalysisTask*)next())) {
      gROOT->cd();
      // Start with memory as current dir and make sure by default histograms do not get attached to files.
      TH1::AddDirectory(kFALSE);
      task->CreateOutputObjects();
      if (!task->CheckPostData()) {
         Error("SlaveBegin","####### IMPORTANT! ####### \n\n\n\
                Task %s (%s) did not call PostData() for all its outputs in (User)CreateOutputObjects()\n\n\
                ####### FIX YOUR CODE, THIS WILL PRODUCE A FATAL ERROR IN FUTURE! ##########", task->GetName(), task->ClassName());
      }
      if (getsysInfo) AliSysInfo::AddStamp(Form("%s_CREATEOUTOBJ",task->ClassName()), 0, itask, 0);
      itask++;
   }
   TH1::AddDirectory(dirStatus);
   if (curdir) curdir->cd();
   if (fDebug > 1) printf("<-AliAnalysisManager::SlaveBegin()\n");
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
   if (!TObject::TestBit(AliAnalysisManager::kTrueNotify)) return kFALSE;

   fTable.Clear("nodelete"); // clearing the hash table may not be needed -> C.L.
   if (fMode == kProofAnalysis) fIsRemote = kTRUE;

   TFile *curfile = fTree->GetCurrentFile();
   if (!curfile) {
      Error("Notify","No current file");
      return kFALSE;
   }   
   
   if (fDebug > 1) printf("->AliAnalysisManager::Notify() file: %s\n", curfile->GetName());
   Int_t run = AliAnalysisManager::GetRunFromAlienPath(curfile->GetName());
   if (run && (run != fRunFromPath)) {
      fRunFromPath = run;
      if (fDebug > 1) printf("   ### run found from path: %d\n", run);
   }
   TIter next(fTasks);
   AliAnalysisTask *task;
	
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

   // Call Notify for all tasks
   while ((task=(AliAnalysisTask*)next())) 
      task->Notify();

   if (fDebug > 1) printf("<-AliAnalysisManager::Notify()\n");
   return kTRUE;
}    

//______________________________________________________________________________
Bool_t AliAnalysisManager::Process(Long64_t)
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

   // This method is obsolete. ExecAnalysis is called instead.
   return kTRUE;
}

//______________________________________________________________________________
void AliAnalysisManager::PackOutput(TList *target)
{
  // Pack all output data containers in the output list. Called at SlaveTerminate
  // stage in PROOF case for each slave.
   if (fDebug > 1) printf("->AliAnalysisManager::PackOutput()\n");
   if (!target) {
      Error("PackOutput", "No target. Exiting.");
      return;
   }
   TDirectory *cdir = gDirectory;
   gROOT->cd();
   if (fInputEventHandler)   fInputEventHandler  ->Terminate();
   if (fOutputEventHandler)  fOutputEventHandler ->Terminate();
   if (fMCtruthEventHandler) fMCtruthEventHandler->Terminate();
   gROOT->cd();

   // Call FinishTaskOutput() for each event loop task (not called for 
   // post-event loop tasks - use Terminate() fo those)
   TIter nexttask(fTasks);
   AliAnalysisTask *task;
   while ((task=(AliAnalysisTask*)nexttask())) {
      if (!task->IsPostEventLoop()) {
         if (fDebug > 1) printf("->FinishTaskOutput: task %s\n", task->GetName());
         task->FinishTaskOutput();
         gROOT->cd();
         if (fDebug > 1) printf("<-FinishTaskOutput: task %s\n", task->GetName());
      }
   }
   // Write statistics message on the workers.
   if (fStatistics) WriteStatisticsMsg(fNcalls);
   
   if (fMode == kProofAnalysis) {
      TIter next(fOutputs);
      AliAnalysisDataContainer *output;
      Bool_t isManagedByHandler = kFALSE;
      TList filestmp;
      filestmp.SetOwner();
      while ((output=(AliAnalysisDataContainer*)next())) {
         // Do not consider outputs of post event loop tasks
         isManagedByHandler = kFALSE;
         if (output->GetProducer() && output->GetProducer()->IsPostEventLoop()) continue;
         const char *filename = output->GetFileName();
         if (!(strcmp(filename, "default")) && fOutputEventHandler) {
            isManagedByHandler = kTRUE;
            printf("#### Handler output. Extra: %s\n", fExtraFiles.Data());
            filename = fOutputEventHandler->GetOutputFileName();
         }
         // Check if data was posted to this container. If not, issue an error.
         if (!output->GetData() && !isManagedByHandler) {
            Error("PackOutput", "No data for output container %s. Forgot to PostData ?", output->GetName());
            continue;
         }   
         if (!output->IsSpecialOutput()) {
            // Normal outputs
            if (strlen(filename) && !isManagedByHandler) {
               // Backup current folder
               TDirectory *opwd = gDirectory;
               // File resident outputs. 
               // Check first if the file exists.
               TString openoption = "RECREATE";
               Bool_t firsttime = kTRUE;
               if (filestmp.FindObject(output->GetFileName())) {
                  firsttime = kFALSE;
               } else {   
                  filestmp.Add(new TNamed(output->GetFileName(),""));
               }   
               if (!gSystem->AccessPathName(output->GetFileName()) && !firsttime) openoption = "UPDATE";
//               TFile *file = AliAnalysisManager::OpenFile(output, openoption, kTRUE);
               // Save data to file, then close.
               if (output->GetData()->InheritsFrom(TCollection::Class())) {
                  // If data is a collection, we set the name of the collection 
                  // as the one of the container and we save as a single key.
                  TCollection *coll = (TCollection*)output->GetData();
                  coll->SetName(output->GetName());
//                  coll->Write(output->GetName(), TObject::kSingleKey);
               } else {
                  if (output->GetData()->InheritsFrom(TTree::Class())) {
                     TFile *file = AliAnalysisManager::OpenFile(output, openoption, kTRUE);
                     // Save data to file, then close.
                     TTree *tree = (TTree*)output->GetData();
                     // Check if tree is in memory
                     if (tree->GetDirectory()==gROOT) tree->SetDirectory(gDirectory);
                     tree->AutoSave();
                     file->Close();
                  } else {
//                     output->GetData()->Write();
                  }   
               }      
               if (fDebug > 1) printf("PackOutput %s: memory merge, file resident output\n", output->GetName());
//               if (fDebug > 2) {
//                  printf("   file %s listing content:\n", filename);
//                  file->ls();
//               }   
               // Clear file list to release object ownership to user.
//               file->Clear();
//               file->Close();
               output->SetFile(NULL);
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
         // Special outputs. The file must be opened and connected to the container.
            TDirectory *opwd = gDirectory;
            TFile *file = output->GetFile();
            if (!file) {
               AliAnalysisTask *producer = output->GetProducer();
               Fatal("PackOutput", 
                     "File %s for special container %s was NOT opened in %s::CreateOutputObjects !!!",
                     output->GetFileName(), output->GetName(), producer->ClassName());
               continue;
            }   
            TString outFilename = file->GetName();
            if (fDebug > 1) printf("PackOutput %s: special output\n", output->GetName());
            if (isManagedByHandler) {
               // Terminate IO for files managed by the output handler
               // file->Write() moved to AOD handler (A.G. 11.01.10)
//               if (file) file->Write();
               if (file && fDebug > 2) {
                  printf("   handled file %s listing content:\n", file->GetName());
                  file->ls();
               }   
               fOutputEventHandler->TerminateIO();
            } else {               
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
               if (fDebug > 2) {
                  printf("   file %s listing content:\n", output->GetFileName());
                  file->ls();
               }
               // Clear file list to release object ownership to user.
//               file->Clear();
               file->Close();
               output->SetFile(NULL);
            }
            // Restore current directory
            if (opwd) opwd->cd();
            // Check if a special output location was provided or the output files have to be merged
            if (strlen(fSpecialOutputLocation.Data())) {
               TString remote = fSpecialOutputLocation;
               remote += "/";
               Int_t gid = gROOT->ProcessLine("gProofServ->GetGroupId();");
               if (remote.BeginsWith("alien:")) {
                  gROOT->ProcessLine("TGrid::Connect(\"alien:\", gProofServ->GetUser());");
                  remote += outFilename;
                  remote.ReplaceAll(".root", Form("_%d.root", gid));
               } else {   
                  remote += Form("%s_%d_", gSystem->HostName(), gid);
                  remote += outFilename;
               }   
               if (fDebug > 1) 
                  Info("PackOutput", "Output file for container %s to be copied \n   at: %s. No merging.",
                       output->GetName(), remote.Data());
               TFile::Cp ( outFilename.Data(), remote.Data() );
               // Copy extra outputs
               if (fExtraFiles.Length() && isManagedByHandler) {
                  TObjArray *arr = fExtraFiles.Tokenize(" ");
                  TObjString *os;
                  TIter nextfilename(arr);
                  while ((os=(TObjString*)nextfilename())) {
                     outFilename = os->GetString();
                     remote = fSpecialOutputLocation;
                     remote += "/";
                     if (remote.BeginsWith("alien://")) {
                        remote += outFilename;
                        remote.ReplaceAll(".root", Form("_%d.root", gid));
                     } else {   
                        remote += Form("%s_%d_", gSystem->HostName(), gid);
                        remote += outFilename;
                     }   
                     if (fDebug > 1) 
                        Info("PackOutput", "Extra AOD file %s to be copied \n   at: %s. No merging.",
                             outFilename.Data(), remote.Data());
                     TFile::Cp ( outFilename.Data(), remote.Data() );
                  }   
                  delete arr;
               }   
            } else {
            // No special location specified-> use TProofOutputFile as merging utility
            // The file at this output slot must be opened in CreateOutputObjects
               if (fDebug > 1) printf("   File for container %s to be merged via file merger...\n", output->GetName());
            }
         }      
      }
   } 
   cdir->cd();
   if (fDebug > 1) printf("<-AliAnalysisManager::PackOutput: output list contains %d containers\n", target->GetSize());
}

//______________________________________________________________________________
void AliAnalysisManager::ImportWrappers(TList *source)
{
// Import data in output containers from wrappers coming in source.
   if (fDebug > 1) printf("->AliAnalysisManager::ImportWrappers()\n");
   TIter next(fOutputs);
   AliAnalysisDataContainer *cont;
   AliAnalysisDataWrapper   *wrap;
   Int_t icont = 0;
   Bool_t inGrid = (fMode == kGridAnalysis)?kTRUE:kFALSE;
   TDirectory *cdir = gDirectory;
   while ((cont=(AliAnalysisDataContainer*)next())) {
      wrap = 0;
      if (cont->GetProducer() && cont->GetProducer()->IsPostEventLoop() && !inGrid) continue;
      if (cont->IsRegisterDataset()) continue;
      const char *filename = cont->GetFileName();
      Bool_t isManagedByHandler = kFALSE;
      if (!(strcmp(filename, "default")) && fOutputEventHandler) {
         isManagedByHandler = kTRUE;
         filename = fOutputEventHandler->GetOutputFileName();
      }
      if (cont->IsSpecialOutput() || inGrid) {
         if (strlen(fSpecialOutputLocation.Data())) continue;
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
         TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
         if (!f) f = TFile::Open(filename, "READ");
         if (!f) {
            Error("ImportWrappers", "Cannot open file %s in read-only mode", filename);
            continue;
         }   
         f->cd();
         TObject *obj = 0;
         // Cd to the directory pointed by the container
         TString folder = cont->GetFolderName();
         if (!folder.IsNull()) f->cd(folder);
         // Try to fetch first an object having the container name.
         obj = gDirectory->Get(cont->GetName());
         if (!obj) {
            Warning("ImportWrappers", "Could not import object of type:%s for container %s in file %s:%s.\n Object will not be available in Terminate(). Try if possible to name the output object as the container (%s) or to embed it in a TList", 
                    cont->GetType()->GetName(), cont->GetName(), filename, cont->GetFolderName(), cont->GetName());
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
         printf("   Importing data for container %s\n", cont->GetName());
         if (strlen(filename)) printf("    -> file %s\n", filename);
         else printf("\n");
      }   
      cont->ImportData(wrap);
   }
   if (cdir) cdir->cd();
   if (fDebug > 1) printf("<-AliAnalysisManager::ImportWrappers(): %d containers imported\n", icont);
}

//______________________________________________________________________________
void AliAnalysisManager::UnpackOutput(TList *source)
{
  // Called by AliAnalysisSelector::Terminate only on the client.
   if (fDebug > 1) printf("->AliAnalysisManager::UnpackOutput()\n");
   if (!source) {
      Error("UnpackOutput", "No target. Exiting.");
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
               if (fDebug > 1) printf("== Executing post event loop task %s\n", task->GetName());
               task->ExecuteTask();
            }   
         }
      }   
   }
   if (fDebug > 1) printf("<-AliAnalysisManager::UnpackOutput()\n");
}

//______________________________________________________________________________
void AliAnalysisManager::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically.
   if (fDebug > 1) printf("->AliAnalysisManager::Terminate()\n");
   TDirectory *cdir = gDirectory;
   gROOT->cd();
   AliAnalysisTask *task;
   AliAnalysisDataContainer *output;
   TIter next(fTasks);
   TStopwatch timer;
   Bool_t getsysInfo = ((fNSysInfo>0) && (fMode==kLocalAnalysis))?kTRUE:kFALSE;
   // Call Terminate() for tasks
   Int_t itask = 0;
   while (!IsSkipTerminate() && (task=(AliAnalysisTask*)next())) {
      // Save all the canvases produced by the Terminate
      TString pictname = Form("%s_%s", task->GetName(), task->ClassName());
      task->Terminate();
      gROOT->cd();
      if (getsysInfo) 
         AliSysInfo::AddStamp(Form("%s_TERMINATE",task->ClassName()),0, itask, 2);
      itask++;   
      if (TObject::TestBit(kSaveCanvases)) {
         if (!gROOT->IsBatch()) {
            if (fDebug>1) printf("Waiting 5 sec for %s::Terminate() to finish drawing ...\n", task->ClassName());
            timer.Start();
            while (timer.CpuTime()<5) {
               timer.Continue();
               gSystem->ProcessEvents();
            }
         }
         Int_t iend = gROOT->GetListOfCanvases()->GetEntries();
         if (iend==0) continue;
         TCanvas *canvas;
         for (Int_t ipict=0; ipict<iend; ipict++) {
            canvas = (TCanvas*)gROOT->GetListOfCanvases()->At(ipict);
            if (!canvas) continue;         
            canvas->SaveAs(Form("%s_%02d.gif", pictname.Data(),ipict));
         } 
         gROOT->GetListOfCanvases()->Delete(); 
      }
   }   
   //
   if (fInputEventHandler)   fInputEventHandler  ->TerminateIO();
   if (fOutputEventHandler)  fOutputEventHandler ->TerminateIO();
   if (fMCtruthEventHandler) fMCtruthEventHandler->TerminateIO();
   gROOT->cd();
   TObjArray *allOutputs = new TObjArray();
   Int_t icont;
   for (icont=0; icont<fOutputs->GetEntriesFast(); icont++) allOutputs->Add(fOutputs->At(icont));
   if (!IsSkipTerminate())
      for (icont=0; icont<fParamCont->GetEntriesFast(); icont++) allOutputs->Add(fParamCont->At(icont));
   TIter next1(allOutputs);
   TString handlerFile = "";
   TString extraOutputs = "";
   if (fOutputEventHandler) {
      handlerFile = fOutputEventHandler->GetOutputFileName();
      extraOutputs = fOutputEventHandler->GetExtraOutputs();
   }
   icont = 0;
   TList filestmp;
   while ((output=(AliAnalysisDataContainer*)next1())) {
      // Special outputs or grid files have the files already closed and written.
      icont++;
      if (fMode == kGridAnalysis && icont<=fOutputs->GetEntriesFast()) continue;
      if (fMode == kProofAnalysis) {
        if (output->IsSpecialOutput() || output->IsRegisterDataset()) continue;
      }  
      const char *filename = output->GetFileName();
      TString openoption = "RECREATE";
      if (!(strcmp(filename, "default"))) continue;
      if (!strlen(filename)) continue;
      if (!output->GetData()) continue;
      TDirectory *opwd = gDirectory;
      TFile *file = output->GetFile();
      if (!file) file = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
      if (!file) {
	      //if (handlerFile == filename && !gSystem->AccessPathName(filename)) openoption = "UPDATE";
         Bool_t firsttime = kTRUE;
         if (filestmp.FindObject(filename) || extraOutputs.Contains(filename)) {
            firsttime = kFALSE;
         } else {   
            filestmp.Add(new TNamed(filename,""));
         }   
         if (!gSystem->AccessPathName(filename) && !firsttime) openoption = "UPDATE";
	      if (fDebug>1) printf("Opening file: %s  option=%s\n",filename, openoption.Data());
         file = new TFile(filename, openoption);
      } else {
         if (fDebug>1) printf("File <%s> already opened with option: <%s> \n", filename, file->GetOption());
         openoption = file->GetOption();
         if (openoption == "READ") {
            if (fDebug>1) printf("...reopening in UPDATE mode\n");
            file->ReOpen("UPDATE");            
         }
      }   
      if (file->IsZombie()) {
         Error("Terminate", "Cannot open output file %s", filename);
         continue;
      }   
      output->SetFile(file);
      file->cd();
      // Check for a folder request
      TString dir = output->GetFolderName();
      if (!dir.IsNull()) {
         if (!file->GetDirectory(dir)) file->mkdir(dir);
         file->cd(dir);
      }  
      if (fDebug > 1) printf("...writing container %s to file %s:%s\n", output->GetName(), file->GetName(), output->GetFolderName());
      if (output->GetData()->InheritsFrom(TCollection::Class())) {
      // If data is a collection, we set the name of the collection 
      // as the one of the container and we save as a single key.
         TCollection *coll = (TCollection*)output->GetData();
         coll->SetName(output->GetName());
         coll->Write(output->GetName(), TObject::kSingleKey);
      } else {
         if (output->GetData()->InheritsFrom(TTree::Class())) {
            TTree *tree = (TTree*)output->GetData();
            tree->SetDirectory(gDirectory);
            tree->AutoSave();
         } else {
            output->GetData()->Write();
         }   
      }      
      if (opwd) opwd->cd();
   }
   gROOT->cd();
   next1.Reset();
   TString copiedFiles;
   while ((output=(AliAnalysisDataContainer*)next1())) {
      // Close all files at output
      TDirectory *opwd = gDirectory;
      if (output->GetFile()) {
         // Clear file list to release object ownership to user.
//         output->GetFile()->Clear();
         output->GetFile()->Close();
         // Copy merged outputs in alien if requested
         if (fSpecialOutputLocation.BeginsWith("alien://")) {
            if (copiedFiles.Contains(output->GetFile()->GetName())) {
               if (opwd) opwd->cd();
               output->SetFile(NULL);
               continue;
            } 
            Info("Terminate", "Copy file %s to %s", output->GetFile()->GetName(),fSpecialOutputLocation.Data()); 
            gROOT->ProcessLine("if (!gGrid) TGrid::Connect(\"alien:\");");
            TFile::Cp(output->GetFile()->GetName(), 
                      Form("%s/%s", fSpecialOutputLocation.Data(), output->GetFile()->GetName()));
            copiedFiles += output->GetFile()->GetName();
         }             
         output->SetFile(NULL);
      }   
      if (opwd) opwd->cd();
   }   
   delete allOutputs;
   //Write statistics information on the client
   if (fStatistics) WriteStatisticsMsg(fNcalls);
   if (getsysInfo) {
      TDirectory *crtdir = gDirectory;
      TFile f("syswatch.root", "RECREATE");
      TH1 *hist;
      TString cut;
      if (!f.IsZombie()) {
         TTree *tree = AliSysInfo::MakeTree("syswatch.log");
         tree->SetName("syswatch");
         tree->SetMarkerStyle(kCircle);
         tree->SetMarkerColor(kBlue);
         tree->SetMarkerSize(0.5);
         if (!gROOT->IsBatch()) {
            tree->SetAlias("event", "id0");
            tree->SetAlias("task",  "id1");
            tree->SetAlias("stage", "id2");
            // Already defined aliases
            // tree->SetAlias("deltaT","stampSec-stampOldSec");
            // tree->SetAlias("T","stampSec-first");
            // tree->SetAlias("deltaVM","(pI.fMemVirtual-pIOld.fMemVirtual)");
            // tree->SetAlias("VM","pI.fMemVirtual");
            TCanvas *canvas = new TCanvas("SysInfo","SysInfo",10,10,1200,800);
            Int_t npads = 1 /*COO plot for all tasks*/ +
                          fTopTasks->GetEntries() /*Exec plot per task*/ +
                          1 /*Terminate plot for all tasks*/ +
                          1; /*vm plot*/
                          
            Int_t iopt = (Int_t)TMath::Sqrt((Double_t)npads);
            if (npads<iopt*(iopt+1))
               canvas->Divide(iopt, iopt+1, 0.01, 0.01);
            else
               canvas->Divide(iopt+1, iopt+1, 0.01, 0.01);
            Int_t ipad = 1;
            // draw the plot of deltaVM for Exec for each task
            for (itask=0; itask<fTopTasks->GetEntriesFast(); itask++) {
               task = (AliAnalysisTask*)fTopTasks->At(itask);
               canvas->cd(ipad++);
               cut = Form("task==%d && stage==1", itask);
               tree->Draw("deltaVM:event",cut,"", 1234567890, 0);
               hist = (TH1*)gPad->GetListOfPrimitives()->FindObject("htemp");            
               if (hist) {
                  hist->SetTitle(Form("%s: Exec dVM[MB]/event", task->GetName()));
                  hist->GetYaxis()->SetTitle("deltaVM [MB]");
               }   
            }
            // Draw the plot of deltaVM for CreateOutputObjects for all tasks
            canvas->cd(ipad++);
            tree->SetMarkerStyle(kFullTriangleUp);
            tree->SetMarkerColor(kRed);
            tree->SetMarkerSize(0.8);
            cut = "task>=0 && task<1000 && stage==0";
            tree->Draw("deltaVM:sname",cut,"", 1234567890, 0);
            hist = (TH1*)gPad->GetListOfPrimitives()->FindObject("htemp");            
            if (hist) {
               hist->SetTitle("Memory in CreateOutputObjects()");
               hist->GetYaxis()->SetTitle("deltaVM [MB]");
               hist->GetXaxis()->SetTitle("task");
            }   
            // draw the plot of deltaVM for Terminate for all tasks
            canvas->cd(ipad++);
            tree->SetMarkerStyle(kOpenSquare);
            tree->SetMarkerColor(kMagenta);
            cut = "task>=0 && task<1000 && stage==2";
            tree->Draw("deltaVM:sname",cut,"", 1234567890, 0);
            hist = (TH1*)gPad->GetListOfPrimitives()->FindObject("htemp");
            if (hist) {
               hist->SetTitle("Memory in Terminate()");
               hist->GetYaxis()->SetTitle("deltaVM [MB]");
               hist->GetXaxis()->SetTitle("task");
            }   
            // Full VM profile
            canvas->cd(ipad++);
            tree->SetMarkerStyle(kFullCircle);
            tree->SetMarkerColor(kGreen);
            cut = Form("task==%d && stage==1",fTopTasks->GetEntriesFast()-1);            
            tree->Draw("VM:event",cut,"", 1234567890, 0);
            hist = (TH1*)gPad->GetListOfPrimitives()->FindObject("htemp");
            if (hist) {
               hist->SetTitle("Virtual memory");
               hist->GetYaxis()->SetTitle("VM [MB]");
            }
            canvas->Modified();   
         }   
         tree->SetMarkerStyle(kCircle);
         tree->SetMarkerColor(kBlue);
         tree->SetMarkerSize(0.5);
         tree->Write();
         f.Close();
         delete tree;
      }
      if (crtdir) crtdir->cd();
   }
   // Validate the output files
   if (ValidateOutputFiles() && fIsRemote && fMode!=kProofAnalysis) {
      ofstream out;
      out.open("outputs_valid", ios::out);
      out.close();
   }
   cdir->cd();      
   if (fDebug > 1) printf("<-AliAnalysisManager::Terminate()\n");
}
//______________________________________________________________________________
void AliAnalysisManager::ProfileTask(Int_t itop, const char *option) const
{
// Profiles the task having the itop index in the list of top (first level) tasks.
   AliAnalysisTask *task = (AliAnalysisTask*)fTopTasks->At(itop);
   if (!task) {
      Error("ProfileTask", "There are only %d top tasks in the manager", fTopTasks->GetEntries());
      return;
   }
   ProfileTask(task->GetName(), option);
}      

//______________________________________________________________________________
void AliAnalysisManager::ProfileTask(const char *name, const char */*option*/) const
{
// Profile a managed task after the execution of the analysis in case NSysInfo
// was used.
   if (gSystem->AccessPathName("syswatch.root")) {
      Error("ProfileTask", "No file syswatch.root found in the current directory");
      return;
   }
   if (gROOT->IsBatch()) return;
   AliAnalysisTask *task = (AliAnalysisTask*)fTopTasks->FindObject(name);
   if (!task) {
      Error("ProfileTask", "No top task named %s known by the manager.", name);
      return;
   }
   Int_t itop = fTopTasks->IndexOf(task);
   Int_t itask = fTasks->IndexOf(task);
   // Create canvas with 2 pads: first draw COO + Terminate, second Exec
   TDirectory *cdir = gDirectory;
   TFile f("syswatch.root");
   TTree *tree = (TTree*)f.Get("syswatch");
   if (!tree) {
      Error("ProfileTask", "No tree named <syswatch> found in file syswatch.root");
      return;
   }   
   if (fDebug > 1) printf("=== Profiling task %s (class %s)\n", name, task->ClassName());
   TCanvas *canvas = new TCanvas(Form("profile_%d",itop),Form("Profile of task %s (class %s)",name,task->ClassName()),10,10,800,600);
   canvas->Divide(2, 2, 0.01, 0.01);
   Int_t ipad = 1;
   TString cut;
   TH1 *hist;
   // VM profile for COO and Terminate methods
   canvas->cd(ipad++);
   cut = Form("task==%d && (stage==0 || stage==2)",itask);
   tree->Draw("deltaVM:sname",cut,"", 1234567890, 0);
   hist = (TH1*)gPad->GetListOfPrimitives()->FindObject("htemp");
   if (hist) {
      hist->SetTitle("Alocated VM[MB] for COO and Terminate");
      hist->GetYaxis()->SetTitle("deltaVM [MB]");
      hist->GetXaxis()->SetTitle("method");
   }   
   // CPU profile per event
   canvas->cd(ipad++);
   cut = Form("task==%d && stage==1",itop);
   tree->Draw("deltaT:event",cut,"", 1234567890, 0);
   hist = (TH1*)gPad->GetListOfPrimitives()->FindObject("htemp");
   if (hist) {
      hist->SetTitle("Execution time per event");
      hist->GetYaxis()->SetTitle("CPU/event [s]");
   }   
   // VM profile for Exec
   canvas->cd(ipad++);
   cut = Form("task==%d && stage==1",itop);
   tree->Draw("deltaVM:event",cut,"", 1234567890, 0);
   hist = (TH1*)gPad->GetListOfPrimitives()->FindObject("htemp");
   if (hist) {
      hist->SetTitle("Alocated VM[MB] per event");
      hist->GetYaxis()->SetTitle("deltaVM [MB]");
   }   
   canvas->Modified();
   delete tree;
   f.Close();
   if (cdir) cdir->cd();
}     

//______________________________________________________________________________
void AliAnalysisManager::AddTask(AliAnalysisTask *task)
{
// Adds a user task to the global list of tasks.
   if (fInitOK) {
      Error("AddTask", "Cannot add task %s since InitAnalysis was already called", task->GetName());
      return;
   }   
      
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
//   kExchangeContainer  = 0, used to exchange data between tasks
//   kInputContainer   = 1, used to store input data
//   kOutputContainer  = 2, used for writing result to a file
// filename: composed by file#folder (e.g. results.root#INCLUSIVE) - will write
// the output object to a folder inside the output file
   if (fContainers->FindObject(name)) {
      Error("CreateContainer","A container named %s already defined !",name);
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
      case kParamContainer:
         fParamCont->Add(cont);
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
   // Reset flag and remove valid_outputs file if exists
   if (fInitOK) return kTRUE;
   if (!gSystem->AccessPathName("outputs_valid"))
      gSystem->Unlink("outputs_valid");
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
   // Initialize requested branch list if needed
   if (!fAutoBranchHandling) {
      next.Reset();
      while ((task=(AliAnalysisTask*)next())) {
         if (!task->HasBranches()) {
            Error("InitAnalysis", "Manual branch loading requested but task %s of type %s does not define branches.\nUse: fBranchNames = \"ESD:br1,br2,...,brN AOD:bra1,bra2,...,braM\"",
                  task->GetName(), task->ClassName());
            return kFALSE;
         }
         if (!fInputEventHandler || !strlen(fInputEventHandler->GetDataType())) {
            Error("InitAnalysis", "Manual branch loading requested but no input handler defined or handler does not define data type.");
            return kFALSE;
         }
         TString taskbranches;
         task->GetBranches(fInputEventHandler->GetDataType(), taskbranches);
         if (taskbranches.IsNull()) {
            Error("InitAnalysis", "Manual branch loading requested but task %s of type %s does not define branches of type %s:",
                  task->GetName(), task->ClassName(), fInputEventHandler->GetDataType());
            return kFALSE;      
         }
         AddBranches(taskbranches);
      }         
   }
   fInitOK = kTRUE;
   return kTRUE;
}   

//______________________________________________________________________________
void AliAnalysisManager::AddBranches(const char *branches)
{
// Add branches to the existing fRequestedBranches.
   TString br(branches);
   TObjArray *arr = br.Tokenize(",");
   TIter next(arr);
   TObject *obj;
   while ((obj=next())) {
      if (!fRequestedBranches.Contains(obj->GetName())) {
         if (!fRequestedBranches.IsNull()) fRequestedBranches += ",";
         fRequestedBranches += obj->GetName();
      }
   }
  delete arr;
}   

//______________________________________________________________________________
void AliAnalysisManager::CheckBranches(Bool_t load)
{
// The method checks the input branches to be loaded during the analysis.
   if (fAutoBranchHandling || fRequestedBranches.IsNull() || !fTree) return;   
   TObjArray *arr = fRequestedBranches.Tokenize(",");
   TIter next(arr);
   TObject *obj;
   while ((obj=next())) {
      TBranch *br = dynamic_cast<TBranch*>(fTable.FindObject(obj->GetName()));
      if (!br) {
         br = fTree->GetBranch(obj->GetName());
         if (!br) {
            Error("CheckBranches", "Could not find branch %s",obj->GetName());
            continue;
         }
      }   
      fTable.Add(br);
      if (load && br->GetReadEntry()!=GetCurrentEntry()) br->GetEntry(GetCurrentEntry());
   }
  delete arr;
}

//______________________________________________________________________________
Bool_t AliAnalysisManager::CheckTasks() const
{
// Check consistency of tasks.
   Int_t ntasks = fTasks->GetEntries();
   if (!ntasks) {
      Error("CheckTasks", "No tasks connected to the manager. This may be due to forgetting to compile the task or to load their library.");
      return kFALSE;
   }
   // Get the pointer to AliAnalysisTaskSE::Class()
   TClass *badptr = (TClass*)gROOT->ProcessLine("AliAnalysisTaskSE::Class()");
   // Loop all tasks to check if their corresponding library was loaded
   TIter next(fTasks);
   TObject *obj;
   while ((obj=next())) {
      if (obj->IsA() == badptr) {
         Error("CheckTasks", "##################\n \
         Class for task %s NOT loaded. You probably forgot to load the library for this task (or compile it dynamically).\n###########################\n",obj->GetName());
         return kFALSE;
      }
   }
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
  
   if (!fAutoBranchHandling && !fRequestedBranches.IsNull()) 
      printf("Requested input branches:\n%s\n", fRequestedBranches.Data());
  
  TString sopt(option);
  sopt.ToUpper();
  
  if (sopt.Contains("ALL"))
  {
    if ( fOutputEventHandler )
    {
      cout << TString('_',78) << endl;
      cout << "OutputEventHandler:" << endl;
      fOutputEventHandler->Print("   ");
    }
  }
}

//______________________________________________________________________________
void AliAnalysisManager::ResetAnalysis()
{
// Reset all execution flags and clean containers.
   CleanContainers();
}

//______________________________________________________________________________
Long64_t AliAnalysisManager::StartAnalysis(const char *type, Long64_t nentries, Long64_t firstentry)
{
// Start analysis having a grid handler.
   if (!fGridHandler) {
      Error("StartAnalysis", "Cannot start analysis providing just the analysis type without a grid handler.");
      Info("===", "Add an AliAnalysisAlien object as plugin for this manager and configure it.");
      return -1;
   }
   TTree *tree = NULL;
   return StartAnalysis(type, tree, nentries, firstentry);
}

//______________________________________________________________________________
Long64_t AliAnalysisManager::StartAnalysis(const char *type, TTree * const tree, Long64_t nentries, Long64_t firstentry)
{
// Start analysis for this manager. Analysis task can be: LOCAL, PROOF, GRID or
// MIX. Process nentries starting from firstentry
   Long64_t retv = 0;
   // Backup current directory and make sure gDirectory points to gROOT
   TDirectory *cdir = gDirectory;
   gROOT->cd();
   if (!fInitOK) {
      Error("StartAnalysis","Analysis manager was not initialized !");
      cdir->cd();
      return -1;
   }
   if (!CheckTasks()) Fatal("StartAnalysis", "Not all needed libraries were loaded");
   if (fDebug > 1) {
      printf("StartAnalysis %s\n",GetName());
      AliLog::SetGlobalLogLevel(AliLog::kInfo);
   }   
   fMaxEntries = nentries;
   fIsRemote = kFALSE;
   TString anaType = type;
   anaType.ToLower();
   fMode = kLocalAnalysis;
   Bool_t runlocalinit = kTRUE;
   if (anaType.Contains("file")) {
      runlocalinit = kFALSE;
      fIsRemote = kTRUE;
   }   
   if (anaType.Contains("proof"))     fMode = kProofAnalysis;
   else if (anaType.Contains("grid")) fMode = kGridAnalysis;
   else if (anaType.Contains("mix"))  fMode = kMixingAnalysis;

   if (fMode == kGridAnalysis) {
      fIsRemote = kTRUE;
      if (!anaType.Contains("terminate")) {
         if (!fGridHandler) {
            Error("StartAnalysis", "Cannot start grid analysis without a grid handler.");
            Info("===", "Add an AliAnalysisAlien object as plugin for this manager and configure it.");
            cdir->cd();
            return -1;
         }
         // Write analysis manager in the analysis file
         cout << "===== RUNNING GRID ANALYSIS: " << GetName() << endl;
         // run local task configuration
         TIter nextTask(fTasks);
         AliAnalysisTask *task;
         while ((task=(AliAnalysisTask*)nextTask())) {
            task->LocalInit();
            gROOT->cd();
         }
         if (!fGridHandler->StartAnalysis(nentries, firstentry)) {
            Info("StartAnalysis", "Grid analysis was stopped and cannot be terminated");
            cdir->cd();
            return -1;
         }   

         // Terminate grid analysis
         if (fSelector && fSelector->GetStatus() == -1) {cdir->cd(); return -1;}
         if (fGridHandler->GetRunMode() == AliAnalysisGrid::kOffline) {cdir->cd(); return 0;}
         cout << "===== MERGING OUTPUTS REGISTERED BY YOUR ANALYSIS JOB: " << GetName() << endl;
         if (!fGridHandler->MergeOutputs()) {
            // Return if outputs could not be merged or if it alien handler
            // was configured for offline mode or local testing.
            cdir->cd();
            return 0;
         }
      }   
      cout << "===== TERMINATING GRID ANALYSIS JOB: " << GetName() << endl;
      ImportWrappers(NULL);
      Terminate();
      cdir->cd();
      return 0;
   }
   TString line;
   SetEventLoop(kFALSE);
   // Enable event loop mode if a tree was provided
   if (tree || fGridHandler || fMode==kMixingAnalysis) SetEventLoop(kTRUE);

   TChain *chain = 0;
   TString ttype = "TTree";
   if (tree && tree->IsA() == TChain::Class()) {
      chain = (TChain*)tree;
      if (!chain || !chain->GetListOfFiles()->First()) {
         Error("StartAnalysis", "Cannot process null or empty chain...");
         cdir->cd();
         return -1;
      }   
      ttype = "TChain";
   }   

   Bool_t getsysInfo = ((fNSysInfo>0) && (fMode==kLocalAnalysis))?kTRUE:kFALSE;
   if (getsysInfo) AliSysInfo::AddStamp("Start", 0);
   // Initialize locally all tasks (happens for all modes)
   TIter next(fTasks);
   AliAnalysisTask *task;
   if (runlocalinit) {
      while ((task=(AliAnalysisTask*)next())) {
         task->LocalInit();
         gROOT->cd();
      }
      if (getsysInfo) AliSysInfo::AddStamp("LocalInit_all", 0);
   }   
   
   switch (fMode) {
      case kLocalAnalysis:
         if (!tree && !fGridHandler) {
            TIter nextT(fTasks);
            // Call CreateOutputObjects for all tasks
            Int_t itask = 0;
            Bool_t dirStatus = TH1::AddDirectoryStatus();
            while ((task=(AliAnalysisTask*)nextT())) {
               TH1::AddDirectory(kFALSE);
               task->CreateOutputObjects();
               if (!task->CheckPostData()) {
                  Error("SlaveBegin","####### IMPORTANT! ####### \n\n\n\
                        Task %s (%s) did not call PostData() for all its outputs in (User)CreateOutputObjects()\n\n\
                        ########### FIX YOUR CODE, THIS WILL PRODUCE A FATAL ERROR IN FUTURE! ###########", task->GetName(), task->ClassName());
               }
               if (getsysInfo) AliSysInfo::AddStamp(Form("%s_CREATEOUTOBJ",task->ClassName()), 0, itask, 0);
               gROOT->cd();
               itask++;
            }   
            TH1::AddDirectory(dirStatus);
            if (IsExternalLoop()) {
               Info("StartAnalysis", "Initialization done. Event loop is controlled externally.\
                     \nSetData for top container, call ExecAnalysis in a loop and then Terminate manually");
               return 0;
            }
            ExecAnalysis();
            Terminate();
            return 0;
         } 
         fSelector = new AliAnalysisSelector(this);
         // Check if a plugin handler is used
         if (fGridHandler) {
            // Get the chain from the plugin
            TString dataType = "esdTree";
            if (fInputEventHandler) {
               dataType = fInputEventHandler->GetDataType();
               dataType.ToLower();
               dataType += "Tree";
            }   
            chain = fGridHandler->GetChainForTestMode(dataType);
            if (!chain) {
               Error("StartAnalysis", "No chain for test mode. Aborting.");
               return -1;
            }
            cout << "===== RUNNING LOCAL ANALYSIS" << GetName() << " ON CHAIN " << chain->GetName() << endl;
            retv = chain->Process(fSelector, "", nentries, firstentry);
            break;
         }
         // Run tree-based analysis via AliAnalysisSelector  
         cout << "===== RUNNING LOCAL ANALYSIS " << GetName() << " ON TREE " << tree->GetName() << endl;
         retv = tree->Process(fSelector, "", nentries, firstentry);
         break;
      case kProofAnalysis:
         fIsRemote = kTRUE;
         // Check if the plugin is used
         if (fGridHandler) {
            return StartAnalysis(type, fGridHandler->GetProofDataSet(), nentries, firstentry);
         }
         if (!gROOT->GetListOfProofs() || !gROOT->GetListOfProofs()->GetEntries()) {
            Error("StartAnalysis", "No PROOF!!! Exiting.");
            cdir->cd();
            return -1;
         }   
         line = Form("gProof->AddInput((TObject*)%p);", this);
         gROOT->ProcessLine(line);
         if (chain) {
            chain->SetProof();
            cout << "===== RUNNING PROOF ANALYSIS " << GetName() << " ON CHAIN " << chain->GetName() << endl;
            retv = chain->Process("AliAnalysisSelector", "", nentries, firstentry);
         } else {
            Error("StartAnalysis", "No chain!!! Exiting.");
            cdir->cd();
            return -1;
         }      
         break;
      case kGridAnalysis:
         fIsRemote = kTRUE;
         if (!anaType.Contains("terminate")) {
            if (!fGridHandler) {
               Error("StartAnalysis", "Cannot start grid analysis without a grid handler.");
               Info("===", "Add an AliAnalysisAlien object as plugin for this manager and configure it.");
               cdir->cd();
               return -1;
            }
            // Write analysis manager in the analysis file
            cout << "===== RUNNING GRID ANALYSIS: " << GetName() << endl;
            // Start the analysis via the handler
            if (!fGridHandler->StartAnalysis(nentries, firstentry)) {
               Info("StartAnalysis", "Grid analysis was stopped and cannot be terminated");
               cdir->cd();
               return -1;
            }   

            // Terminate grid analysis
            if (fSelector && fSelector->GetStatus() == -1) {cdir->cd(); return -1;}
            if (fGridHandler->GetRunMode() == AliAnalysisGrid::kOffline) {cdir->cd(); return 0;}
            cout << "===== MERGING OUTPUTS REGISTERED BY YOUR ANALYSIS JOB: " << GetName() << endl;
            if (!fGridHandler->MergeOutputs()) {
               // Return if outputs could not be merged or if it alien handler
               // was configured for offline mode or local testing.
               cdir->cd();
               return 0;
            }
         }   
         cout << "===== TERMINATING GRID ANALYSIS JOB: " << GetName() << endl;
         ImportWrappers(NULL);
         Terminate();
         cdir->cd();
         return 0;
      case kMixingAnalysis:   
         // Run event mixing analysis
         if (!fEventPool) {
            Error("StartAnalysis", "Cannot run event mixing without event pool");
            cdir->cd();
            return -1;
         }
         cout << "===== RUNNING EVENT MIXING ANALYSIS " << GetName() << endl;
         fSelector = new AliAnalysisSelector(this);
         while ((chain=fEventPool->GetNextChain())) {
            next.Reset();
            // Call NotifyBinChange for all tasks
            while ((task=(AliAnalysisTask*)next()))
               if (!task->IsPostEventLoop()) task->NotifyBinChange();
            retv = chain->Process(fSelector);
            if (retv < 0) {
               Error("StartAnalysis", "Mixing analysis failed");
               cdir->cd();
               return retv;
            }   
         }
         PackOutput(fSelector->GetOutputList());
         Terminate();
   }
   cdir->cd();
   return retv;
}   

//______________________________________________________________________________
Long64_t AliAnalysisManager::StartAnalysis(const char *type, const char *dataset, Long64_t nentries, Long64_t firstentry)
{
// Start analysis for this manager on a given dataset. Analysis task can be: 
// LOCAL, PROOF or GRID. Process nentries starting from firstentry.
   if (!fInitOK) {
      Error("StartAnalysis","Analysis manager was not initialized !");
      return -1;
   }
   fIsRemote = kTRUE;
   if (fDebug > 1) printf("StartAnalysis %s\n",GetName());
   TString anaType = type;
   anaType.ToLower();
   if (!anaType.Contains("proof")) {
      Error("StartAnalysis", "Cannot process datasets in %s mode. Try PROOF.", type);
      return -1;
   }   
   fMode = kProofAnalysis;
   TString line;
   SetEventLoop(kTRUE);
   // Set the dataset flag
   TObject::SetBit(kUseDataSet);
   fTree = 0;
   if (fGridHandler) {
      // Start proof analysis using the grid handler
      if (!fGridHandler->StartAnalysis(nentries, firstentry)) {
         Error("StartAnalysis", "The grid plugin could not start PROOF analysis");
         return -1;
      }
      // Check if the plugin is in test mode
      if (fGridHandler->GetRunMode() == AliAnalysisGrid::kTest) {
         dataset = "test_collection";
      } else {
         dataset = fGridHandler->GetProofDataSet();
      }
   }   

   if (!gROOT->GetListOfProofs() || !gROOT->GetListOfProofs()->GetEntries()) {
      Error("StartAnalysis", "No PROOF!!! Exiting.");
      return -1;
   }   

   // Initialize locally all tasks
   TIter next(fTasks);
   AliAnalysisTask *task;
   while ((task=(AliAnalysisTask*)next())) {
      task->LocalInit();
   }
   
   line = Form("gProof->AddInput((TObject*)%p);", this);
   gROOT->ProcessLine(line);
   Long_t retv;
   line = Form("gProof->Process(\"%s\", \"AliAnalysisSelector\", \"\", %lld, %lld);",
               dataset, nentries, firstentry);
   cout << "===== RUNNING PROOF ANALYSIS " << GetName() << " ON DATASET " << dataset << endl;
   retv = (Long_t)gROOT->ProcessLine(line);
   return retv;
}   

//______________________________________________________________________________
TFile *AliAnalysisManager::OpenFile(AliAnalysisDataContainer *cont, const char *option, Bool_t ignoreProof)
{
// Opens according the option the file specified by cont->GetFileName() and changes
// current directory to cont->GetFolderName(). If the file was already opened, it
// checks if the option UPDATE was preserved. File open via TProofOutputFile can
// be optionally ignored.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  TString filename = cont->GetFileName();
  TFile *f = NULL;
  if (filename.IsNull()) {
    ::Error("AliAnalysisManager::OpenFile", "No file name specified for container %s", cont->GetName());
    return NULL;
  }
  if (mgr->GetAnalysisType()==AliAnalysisManager::kProofAnalysis && cont->IsSpecialOutput()
      && !ignoreProof)
    f = mgr->OpenProofFile(cont,option);
  else {
    // Check first if the file is already opened
    f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
    if (f) {
      // Check if option "UPDATE" was preserved 
      TString opt(option);
      opt.ToUpper();
      if ((opt=="UPDATE") && (opt!=f->GetOption())) 
        ::Info("AliAnalysisManager::OpenFile", "File %s already opened in %s mode!", cont->GetFileName(), f->GetOption());
    } else {
      f = TFile::Open(filename, option);
    }    
  }   
  if (f && !f->IsZombie() && !f->TestBit(TFile::kRecovered)) {
    cont->SetFile(f);
    // Cd to file
    f->cd();
    // Check for a folder request
    TString dir = cont->GetFolderName(); 
    if (!dir.IsNull()) {
      if (!f->GetDirectory(dir)) f->mkdir(dir);
      f->cd(dir);
    }
    return f;
  }
  ::Fatal("AliAnalysisManager::OpenFile", "File %s could not be opened", filename.Data());
  cont->SetFile(NULL);
  return NULL;
}    
 
//______________________________________________________________________________
TFile *AliAnalysisManager::OpenProofFile(AliAnalysisDataContainer *cont, const char *option, const char *extaod)
{
// Opens a special output file used in PROOF.
  TString line;
  TString filename = cont->GetFileName();
  if (cont == fCommonOutput) {
     if (fOutputEventHandler) {
        if (strlen(extaod)) filename = extaod;
        filename = fOutputEventHandler->GetOutputFileName();
     }   
     else Fatal("OpenProofFile","No output container. Exiting.");
  }   
  TFile *f = NULL;
  if (fMode!=kProofAnalysis || !fSelector) {
    Fatal("OpenProofFile","Cannot open PROOF file %s: no PROOF or selector",filename.Data());
    return NULL;
  } 
  if (fSpecialOutputLocation.Length()) {
    f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
    if (f) {
      // Check if option "UPDATE" was preserved 
      TString opt(option);
      opt.ToUpper();
      if ((opt=="UPDATE") && (opt!=f->GetOption()))
        ::Info("OpenProofFile", "File %s already opened in %s mode!", cont->GetFileName(), f->GetOption());
    } else {
      f = new TFile(filename, option);
    }
    if (f && !f->IsZombie() && !f->TestBit(TFile::kRecovered)) {
      cont->SetFile(f);
      // Cd to file
      f->cd();
      // Check for a folder request
      TString dir = cont->GetFolderName(); 
      if (dir.Length()) {
        if (!f->GetDirectory(dir)) f->mkdir(dir);
        f->cd(dir);
      }      
      return f;
    }
    Fatal("OpenProofFile", "File %s could not be opened", cont->GetFileName());
    cont->SetFile(NULL);
    return NULL;       
  }
  // Check if there is already a proof output file in the output list
  TObject *pof = fSelector->GetOutputList()->FindObject(filename);
  if (pof) {
    // Get the actual file
    line = Form("((TProofOutputFile*)%p)->GetFileName();", pof);
    filename = (const char*)gROOT->ProcessLine(line);
    if (fDebug>1) {
      printf("File: %s already booked via TProofOutputFile\n", filename.Data());
    }  
    f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
    if (!f) {
       Fatal("OpenProofFile", "Proof output file found but no file opened for %s", filename.Data());
       return NULL;
    }   
    // Check if option "UPDATE" was preserved 
    TString opt(option);
    opt.ToUpper();
    if ((opt=="UPDATE") && (opt!=f->GetOption())) 
      Fatal("OpenProofFile", "File %s already opened, but not in UPDATE mode!", cont->GetFileName());
  } else {
    if (cont->IsRegisterDataset()) {
      TString dsetName = filename;
      dsetName.ReplaceAll(".root", cont->GetTitle());
      dsetName.ReplaceAll(":","_");
      if (fDebug>1) printf("Booking dataset: %s\n", dsetName.Data());
      line = Form("TProofOutputFile *pf = new TProofOutputFile(\"%s\", \"DROV\", \"%s\");", filename.Data(), dsetName.Data());
    } else {
      if (fDebug>1) printf("Booking TProofOutputFile: %s to be merged\n", filename.Data());
      line = Form("TProofOutputFile *pf = new TProofOutputFile(\"%s\");", filename.Data());
    }
    if (fDebug > 1) printf("=== %s\n", line.Data());
    gROOT->ProcessLine(line);
    line = Form("pf->OpenFile(\"%s\");", option);
    gROOT->ProcessLine(line);
    f = gFile;
    if (fDebug > 1) {
      gROOT->ProcessLine("pf->Print()");
      printf(" == proof file name: %s", f->GetName());
    }   
    // Add to proof output list
    line = Form("((TList*)%p)->Add(pf);",fSelector->GetOutputList());
    if (fDebug > 1) printf("=== %s\n", line.Data());
    gROOT->ProcessLine(line);
  }
  if (f && !f->IsZombie() && !f->TestBit(TFile::kRecovered)) {
    cont->SetFile(f);
    // Cd to file
    f->cd();
    // Check for a folder request
    TString dir = cont->GetFolderName(); 
    if (!dir.IsNull()) {
      if (!f->GetDirectory(dir)) f->mkdir(dir);
      f->cd(dir);
    }
    return f;
  }
  Fatal("OpenProofFile", "File %s could not be opened", cont->GetFileName());
  cont->SetFile(NULL);  
  return NULL;
}   

//______________________________________________________________________________
void AliAnalysisManager::ExecAnalysis(Option_t *option)
{
// Execute analysis.
   static Long64_t nentries = 0;
   static TTree *lastTree = 0;
   static TStopwatch *timer = new TStopwatch();
   // Only the first call to Process will trigger a true Notify. Other Notify
   // coming before is ignored.
   if (!TObject::TestBit(AliAnalysisManager::kTrueNotify)) {
      TObject::SetBit(AliAnalysisManager::kTrueNotify);
      Notify();
   }   
   if (fDebug > 0) printf("MGR: Processing event #%d\n", fNcalls);
   else {
      if (fTree && (fTree != lastTree)) {
         nentries += fTree->GetEntries();
         lastTree = fTree;
      }   
      if (!fNcalls) timer->Start();
      if (!fIsRemote && TObject::TestBit(kUseProgressBar)) ProgressBar("Processing event", fNcalls, TMath::Min(fMaxEntries,nentries), timer, kFALSE);
   }
   gROOT->cd();
   TDirectory *cdir = gDirectory;
   Bool_t getsysInfo = ((fNSysInfo>0) && (fMode==kLocalAnalysis))?kTRUE:kFALSE;
   if (getsysInfo && ((fNcalls%fNSysInfo)==0)) AliSysInfo::AddStamp("Exec_start", (Int_t)fNcalls);
   if (!fInitOK) {
      Error("ExecAnalysis", "Analysis manager was not initialized !");
      cdir->cd();
      return;
   }
   fNcalls++;
   AliAnalysisTask *task;
   // Check if the top tree is active.
   if (fTree) {
      if (getsysInfo && ((fNcalls%fNSysInfo)==0)) 
         AliSysInfo::AddStamp("Handlers_BeginEventGroup",fNcalls, 1002, 0);
      TIter next(fTasks);
   // De-activate all tasks
      while ((task=(AliAnalysisTask*)next())) task->SetActive(kFALSE);
      AliAnalysisDataContainer *cont = fCommonInput;
      if (!cont) cont = (AliAnalysisDataContainer*)fInputs->At(0);
      if (!cont) {
	      Error("ExecAnalysis","Cannot execute analysis in TSelector mode without at least one top container");
         cdir->cd();
         return;
      }   
      cont->SetData(fTree); // This will notify all consumers
      Long64_t entry = fTree->GetTree()->GetReadEntry();      
//
//    Call BeginEvent() for optional input/output and MC services 
      if (fInputEventHandler)   fInputEventHandler  ->BeginEvent(entry);
      if (fOutputEventHandler)  fOutputEventHandler ->BeginEvent(entry);
      if (fMCtruthEventHandler) fMCtruthEventHandler->BeginEvent(entry);
      gROOT->cd();
      if (getsysInfo && ((fNcalls%fNSysInfo)==0)) 
         AliSysInfo::AddStamp("Handlers_BeginEvent",fNcalls, 1000, 0);
//
//    Execute the tasks
//      TIter next1(cont->GetConsumers());
      TIter next1(fTopTasks);
      Int_t itask = 0;
      while ((task=(AliAnalysisTask*)next1())) {
         if (fDebug >1) {
            cout << "    Executing task " << task->GetName() << endl;
         }   	 
         task->ExecuteTask(option);
         gROOT->cd();
         if (getsysInfo && ((fNcalls%fNSysInfo)==0)) 
            AliSysInfo::AddStamp(task->ClassName(), fNcalls, itask, 1);
         itask++;   
      }
//
//    Call FinishEvent() for optional output and MC services 
      if (fInputEventHandler)   fInputEventHandler  ->FinishEvent();
      if (fOutputEventHandler)  fOutputEventHandler ->FinishEvent();
      if (fMCtruthEventHandler) fMCtruthEventHandler->FinishEvent();
      // Gather system information if requested
      if (getsysInfo && ((fNcalls%fNSysInfo)==0)) 
         AliSysInfo::AddStamp("Handlers_FinishEvent",fNcalls, 1001, 1);
      cdir->cd();   
      return;
   }   
   // The event loop is not controlled by TSelector   
//
//  Call BeginEvent() for optional input/output and MC services 
   if (fInputEventHandler)   fInputEventHandler  ->BeginEvent(-1);
   if (fOutputEventHandler)  fOutputEventHandler ->BeginEvent(-1);
   if (fMCtruthEventHandler) fMCtruthEventHandler->BeginEvent(-1);
   gROOT->cd();
   if (getsysInfo && ((fNcalls%fNSysInfo)==0)) 
      AliSysInfo::AddStamp("Handlers_BeginEvent",fNcalls, 1000, 0);
   TIter next2(fTopTasks);
   while ((task=(AliAnalysisTask*)next2())) {
      task->SetActive(kTRUE);
      if (fDebug > 1) {
         cout << "    Executing task " << task->GetName() << endl;
      }   
      task->ExecuteTask(option);
      gROOT->cd();
   }   
//
// Call FinishEvent() for optional output and MC services 
   if (fInputEventHandler)   fInputEventHandler  ->FinishEvent();
   if (fOutputEventHandler)  fOutputEventHandler ->FinishEvent();
   if (fMCtruthEventHandler) fMCtruthEventHandler->FinishEvent();
   if (getsysInfo && ((fNcalls%fNSysInfo)==0)) 
      AliSysInfo::AddStamp("Handlers_FinishEvent",fNcalls, 1000, 1);
   cdir->cd();   
}

//______________________________________________________________________________
Bool_t AliAnalysisManager::IsPipe(std::ostream &out)
{
// Check if the stdout is connected to a pipe (C.Holm)
  Bool_t ispipe = kFALSE;
  out.seekp(0, std::ios_base::cur);
  if (out.fail()) {
    out.clear();
    if (errno == ESPIPE) ispipe = kTRUE;
  }
  return ispipe;
}
   
//______________________________________________________________________________
void AliAnalysisManager::SetInputEventHandler(AliVEventHandler* const handler)
{
// Set the input event handler and create a container for it.
   fInputEventHandler   = handler;
   if (!fCommonInput) fCommonInput = CreateContainer("cAUTO_INPUT", TChain::Class(), AliAnalysisManager::kInputContainer);
}

//______________________________________________________________________________
void AliAnalysisManager::SetOutputEventHandler(AliVEventHandler* const handler)
{
// Set the input event handler and create a container for it.
   fOutputEventHandler   = handler;
   if (!fCommonOutput) fCommonOutput = CreateContainer("cAUTO_OUTPUT", TTree::Class(), AliAnalysisManager::kOutputContainer, "default");
   fCommonOutput->SetSpecialOutput();
}

//______________________________________________________________________________
void AliAnalysisManager::SetDebugLevel(UInt_t level)
{
// Set verbosity of the analysis manager. If the progress bar is used, the call is ignored
   if (TObject::TestBit(kUseProgressBar)) {
      Info("SetDebugLevel","Ignored. Disable the progress bar first.");
      return;
   }
   fDebug = level;
}
   
//______________________________________________________________________________
void AliAnalysisManager::SetUseProgressBar(Bool_t flag, Int_t freq)
{
// Enable a text mode progress bar. Resets debug level to 0.
   Info("SetUseProgressBar", "Progress bar enabled, updated every %d events.\n  ### NOTE: Debug level reset to 0 ###", freq);
   TObject::SetBit(kUseProgressBar,flag);
   fPBUpdateFreq = freq;
   fDebug = 0;
}   

//______________________________________________________________________________
void AliAnalysisManager::RegisterExtraFile(const char *fname)
{
// This method is used externally to register output files which are not
// connected to any output container, so that the manager can properly register,
// retrieve or merge them when running in distributed mode. The file names are
// separated by blancs. The method has to be called in MyAnalysisTask::LocalInit().
   if (fExtraFiles.Contains(fname)) return;
   if (fExtraFiles.Length()) fExtraFiles += " ";
   fExtraFiles += fname;
}

//______________________________________________________________________________
Bool_t AliAnalysisManager::GetFileFromWrapper(const char *filename, const TList *source)
{
// Copy a file from the location specified ina the wrapper with the same name from the source list.
   char fullPath[512];
   char chUrl[512];
   char tmp[1024];
   TObject *pof =  source->FindObject(filename);
   if (!pof || !pof->InheritsFrom("TProofOutputFile")) {
      Error("GetFileFromWrapper", "TProofOutputFile object not found in output list for file %s", filename);
      return kFALSE;
   }
   gROOT->ProcessLine(Form("sprintf((char*)%p, \"%%s\", ((TProofOutputFile*)%p)->GetOutputFileName());", fullPath, pof));
   gROOT->ProcessLine(Form("sprintf((char*)%p, \"%%s\", gProof->GetUrl());",chUrl));
   TString clientUrl(chUrl);
   TString fullPath_str(fullPath);
   if (clientUrl.Contains("localhost")){
      TObjArray* array = fullPath_str.Tokenize ( "//" );
      TObjString *strobj = ( TObjString *)array->At(1);
      TObjArray* arrayPort = strobj->GetString().Tokenize ( ":" );
      TObjString *strobjPort = ( TObjString *) arrayPort->At(1);
      fullPath_str.ReplaceAll(strobj->GetString().Data(),"localhost:PORT");
      fullPath_str.ReplaceAll(":PORT",Form(":%s",strobjPort->GetString().Data()));
      if (fDebug > 1) Info("GetFileFromWrapper","Using tunnel from %s to %s",fullPath_str.Data(),filename);
      delete arrayPort;
      delete array;
   }
   else if (clientUrl.Contains("__lite__")) { 
     // Special case for ProofLite environement - get file info and copy. 
     gROOT->ProcessLine(Form("sprintf((char*)%p,\"%%s\",((TProofOutputFile*)%p)->GetDir());", tmp, pof));
     fullPath_str = Form("%s/%s", tmp, fullPath);
   }
   if (fDebug > 1) 
     Info("GetFileFromWrapper","Copying file %s from PROOF scratch space to %s", fullPath_str.Data(),filename);
   Bool_t gotit = TFile::Cp(fullPath_str.Data(), filename); 
   if (!gotit)
      Error("GetFileFromWrapper", "Could not get file %s from proof scratch space", filename);
   return gotit;
}

//______________________________________________________________________________
void AliAnalysisManager::GetAnalysisTypeString(TString &type) const
{
// Fill analysis type in the provided string.
   switch (fMode) {
      case kLocalAnalysis:
         type = "local";
         return;
      case kProofAnalysis:
         type = "proof";
         return;
      case kGridAnalysis:
         type = "grid";
         return;
      case kMixingAnalysis:
         type = "mix";
   }
}

//______________________________________________________________________________
Bool_t AliAnalysisManager::ValidateOutputFiles() const
{
// Validate all output files.
   TIter next(fOutputs);
   AliAnalysisDataContainer *output;
   TDirectory *cdir = gDirectory;
   TString openedFiles;
   while ((output=(AliAnalysisDataContainer*)next())) {
      if (output->IsRegisterDataset()) continue;
      TString filename = output->GetFileName();
      if (filename == "default") {
         if (!fOutputEventHandler) continue;
         filename = fOutputEventHandler->GetOutputFileName();
         // Main AOD may not be there
         if (gSystem->AccessPathName(filename)) continue;
      }
      // Check if the file is closed
      if (openedFiles.Contains(filename)) continue;;
      TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
      if (file) {
         Warning("ValidateOutputs", "File %s was not closed. Closing.", filename.Data());
         // Clear file list to release object ownership to user.
//         file->Clear();
         file->Close();
      }
      file = TFile::Open(filename);
      if (!file || file->IsZombie() || file->TestBit(TFile::kRecovered)) {
         Error("ValidateOutputs", "Output file <%s> was not created or invalid", filename.Data());
         cdir->cd();
         return kFALSE;
      }
      file->Close();
      openedFiles += filename;
      openedFiles += " ";
   }
   cdir->cd();
   return kTRUE;
}   

//______________________________________________________________________________
void AliAnalysisManager::ProgressBar(const char *opname, Long64_t current, Long64_t size, TStopwatch * const watch, Bool_t last, Bool_t refresh)
{
// Implements a nice text mode progress bar.
   static Long64_t icount = 0;
   static TString oname;
   static TString nname;
   static Long64_t ocurrent = 0;
   static Long64_t osize = 0;
   static Int_t oseconds = 0;
   static TStopwatch *owatch = 0;
   static Bool_t oneoftwo = kFALSE;
   static Int_t nrefresh = 0;
   static Int_t nchecks = 0;
   static char lastChar = 0;
   const char symbol[4] = {'-','\\','|','/'}; 
   
   if (!lastChar) lastChar = (IsPipe(std::cerr))?'\r':'\n';
   if (!refresh) {
      nrefresh = 0;
      if (!size) return;
      owatch = watch;
      oname = opname;
      ocurrent = TMath::Abs(current);
      osize = TMath::Abs(size);
      if (ocurrent > osize) ocurrent=osize;
   } else {
      nrefresh++;
      if (!osize) return;
   }     
   if ((current % fPBUpdateFreq) != 0) return;
   icount++;
   char progress[11] = "          ";
   Int_t ichar = icount%4;
   Double_t time = 0.;
   Int_t hours = 0;
   Int_t minutes = 0;
   Int_t seconds = 0;
   if (owatch && !last) {
      owatch->Stop();
      time = owatch->RealTime();
      seconds   = int(time) % 60;
      minutes   = (int(time) / 60) % 60;
      hours     = (int(time) / 60 / 60);
      if (refresh)  {
         if (oseconds==seconds) {
            owatch->Continue();
            return;
         }
         oneoftwo = !oneoftwo;   
      }
      oseconds = seconds;   
   }
   if (refresh && oneoftwo) {
      nname = oname;
      if (nchecks <= 0) nchecks = nrefresh+1;
      Int_t pctdone = (Int_t)(100.*nrefresh/nchecks);
      oname = Form("     == %d%% ==", pctdone);
   }         
   Double_t percent = 100.0*ocurrent/osize;
   Int_t nchar = Int_t(percent/10);
   if (nchar>10) nchar=10;
   Int_t i;
   for (i=0; i<nchar; i++)  progress[i] = '=';
   progress[nchar] = symbol[ichar];
   for (i=nchar+1; i<10; i++) progress[i] = ' ';
   progress[10] = '\0';
   oname += "                    ";
   oname.Remove(20);
   if(size<10000) fprintf(stderr, "%s [%10s] %4lld ", oname.Data(), progress, ocurrent);
   else if(size<100000) fprintf(stderr, "%s [%10s] %5lld ",oname.Data(), progress, ocurrent);
   else fprintf(stderr, "%s [%10s] %7lld ",oname.Data(), progress, ocurrent);
   if (time>0.) {
     Int_t full   = Int_t(ocurrent > 0 ? 
			  time * (float(osize)/ocurrent) + .5 : 
			  99*3600+59*60+59); 
     Int_t remain = Int_t(full - time);
     Int_t rsec   = remain % 60;
     Int_t rmin   = (remain / 60) % 60;
     Int_t rhour  = (remain / 60 / 60);
     fprintf(stderr, "[%6.2f %%]   TIME %.2d:%.2d:%.2d  ETA %.2d:%.2d:%.2d%c",
	     percent, hours, minutes, seconds, rhour, rmin, rsec, lastChar);
   }
   else fprintf(stderr, "[%6.2f %%]%c", percent, lastChar);
   if (refresh && oneoftwo) oname = nname;
   if (owatch) owatch->Continue();
   if (last) {
      icount = 0;
      owatch = 0;
      ocurrent = 0;
      osize = 0;
      oseconds = 0;
      oneoftwo = kFALSE;
      nrefresh = 0;
      fprintf(stderr, "\n");
   }   
}

//______________________________________________________________________________
void AliAnalysisManager::DoLoadBranch(const char *name) 
{
  // Get tree and load branch if needed.
  static Long64_t crtEntry = -100;

  if (fAutoBranchHandling || !fTree)
    return;

  TBranch *br = dynamic_cast<TBranch*>(fTable.FindObject(name));
  if (!br) {
    br = fTree->GetBranch(name);
    if (!br) {
      Error("DoLoadBranch", "Could not find branch %s",name);
      return;
    }
    fTable.Add(br);
  }
  if (br->GetReadEntry()==fCurrentEntry) return;
  Int_t ret = br->GetEntry(GetCurrentEntry());
  if (ret<0) {
    Error("DoLoadBranch", "Could not load entry %lld from branch %s",GetCurrentEntry(), name);
    if (crtEntry != fCurrentEntry) {
      CountEvent(1,0,1,0);
      crtEntry = fCurrentEntry;
    }  
  } else {
    if (crtEntry != fCurrentEntry) {
      CountEvent(1,1,0,0);
      crtEntry = fCurrentEntry;
    }
  }
}

//______________________________________________________________________________
void AliAnalysisManager::AddStatisticsTask(UInt_t offlineMask)
{
// Add the statistics task to the manager.
  if (fStatistics) {
     Info("AddStatisticsTask", "Already added");
     return;
  }
  TString line = Form("AliAnalysisTaskStat::AddToManager(%u);", offlineMask);
  gROOT->ProcessLine(line);
}  

//______________________________________________________________________________
void AliAnalysisManager::CountEvent(Int_t ninput, Int_t nprocessed, Int_t nfailed, Int_t naccepted)
{
// Bookkeep current event;
   if (!fStatistics) return;
   fStatistics->AddInput(ninput);
   fStatistics->AddProcessed(nprocessed);
   fStatistics->AddFailed(nfailed);
   fStatistics->AddAccepted(naccepted);
}   

//______________________________________________________________________________
void AliAnalysisManager::AddStatisticsMsg(const char *line)
{
// Add a line in the statistics message. If available, the statistics message is written
// at the end of the SlaveTerminate phase on workers AND at the end of Terminate
// on the client.
   if (!strlen(line)) return;
   if (!fStatisticsMsg.IsNull()) fStatisticsMsg += "\n";
   fStatisticsMsg += line;
}

//______________________________________________________________________________
void AliAnalysisManager::WriteStatisticsMsg(Int_t)
{
// If fStatistics is present, write the file in the format ninput_nprocessed_nfailed_naccepted.stat
   static Bool_t done = kFALSE;
   if (done) return;
   done = kTRUE;
   if (!fStatistics) return;
   ofstream out;
   AddStatisticsMsg(Form("Number of input events:        %lld",fStatistics->GetNinput()));
   AddStatisticsMsg(Form("Number of processed events:    %lld",fStatistics->GetNprocessed()));      
   AddStatisticsMsg(Form("Number of failed events (I/O): %lld",fStatistics->GetNfailed()));
   AddStatisticsMsg(Form("Number of accepted events for mask %s: %lld", AliAnalysisStatistics::GetMaskAsString(fStatistics->GetOfflineMask()), fStatistics->GetNaccepted()));
   out.open(Form("%lld_%lld_%lld_%lld.stat",fStatistics->GetNinput(),
                 fStatistics->GetNprocessed(),fStatistics->GetNfailed(),
                 fStatistics->GetNaccepted()), ios::out);      
   out << fStatisticsMsg << endl;
   out.close();
}

//______________________________________________________________________________
const char* AliAnalysisManager::GetOADBPath()
{
// returns the path of the OADB
// this static function just depends on environment variables

   static TString oadbPath;

   if (gSystem->Getenv("OADB_PATH"))
      oadbPath = gSystem->Getenv("OADB_PATH");
   else if (gSystem->Getenv("ALICE_ROOT"))
      oadbPath.Form("%s/OADB", gSystem->Getenv("ALICE_ROOT"));
   else
      ::Fatal("AliAnalysisManager::GetOADBPath", "Cannot figure out AODB path. Define ALICE_ROOT or OADB_PATH!");
      
   return oadbPath;
}

//______________________________________________________________________________
void AliAnalysisManager::SetGlobalStr(const char *key, const char *value)
{
// Define a custom string variable mapped to a global unique name. The variable
// can be then retrieved by a given analysis macro via GetGlobalStr(key).
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AliAnalysisManager::SetGlobalStr", "No analysis manager defined");
      return;
   }   
   Bool_t valid = kFALSE;
   TString existing = AliAnalysisManager::GetGlobalStr(key, valid);
   if (valid) {
      ::Error("AliAnalysisManager::SetGlobalStr", "Global %s = %s already defined.", key, existing.Data());
      return;
   }
   mgr->GetGlobals()->Add(new TObjString(key), new TObjString(value));
}

//______________________________________________________________________________
const char *AliAnalysisManager::GetGlobalStr(const char *key, Bool_t &valid)
{
// Static method to retrieve a global variable defined via SetGlobalStr.
   valid = kFALSE;
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) return 0;
   TObject *value = mgr->GetGlobals()->GetValue(key);
   if (!value) return 0;
   valid = kTRUE;
   return value->GetName();
}

//______________________________________________________________________________
void AliAnalysisManager::SetGlobalInt(const char *key, Int_t value)
{
// Define a custom integer variable mapped to a global unique name. The variable
// can be then retrieved by a given analysis macro via GetGlobalInt(key).
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AliAnalysisManager::SetGlobalStr", "No analysis manager defined");
      return;
   }   
   Bool_t valid = kFALSE;
   Int_t existing = AliAnalysisManager::GetGlobalInt(key, valid);
   if (valid) {
      ::Error("AliAnalysisManager::SetGlobalInt", "Global %s = %i already defined.", key, existing);
      return;
   }
   mgr->GetGlobals()->Add(new TObjString(key), new TObjString(TString::Format("%i",value)));
}

//______________________________________________________________________________
Int_t AliAnalysisManager::GetGlobalInt(const char *key, Bool_t &valid)
{
// Static method to retrieve a global variable defined via SetGlobalInt.
   valid = kFALSE;
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) return 0;
   TObject *value = mgr->GetGlobals()->GetValue(key);
   if (!value) return 0;
   valid = kTRUE;
   TString s = value->GetName();
   return s.Atoi();
}

//______________________________________________________________________________
void AliAnalysisManager::SetGlobalDbl(const char *key, Double_t value)
{
// Define a custom double precision variable mapped to a global unique name. The variable
// can be then retrieved by a given analysis macro via GetGlobalInt(key).
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AliAnalysisManager::SetGlobalStr", "No analysis manager defined");
      return;
   }   
   Bool_t valid = kFALSE;
   Double_t existing = AliAnalysisManager::GetGlobalDbl(key, valid);
   if (valid) {
      ::Error("AliAnalysisManager::SetGlobalInt", "Global %s = %g already defined.", key, existing);
      return;
   }
   mgr->GetGlobals()->Add(new TObjString(key), new TObjString(TString::Format("%f.16",value)));
}

//______________________________________________________________________________
Double_t AliAnalysisManager::GetGlobalDbl(const char *key, Bool_t &valid)
{
// Static method to retrieve a global variable defined via SetGlobalDbl.
   valid = kFALSE;
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) return 0;
   TObject *value = mgr->GetGlobals()->GetValue(key);
   if (!value) return 0;
   valid = kTRUE;
   TString s = value->GetName();
   return s.Atof();
}
