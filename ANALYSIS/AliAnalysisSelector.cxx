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
//   AliAnalysisSelector - A transparent selector to be created by 
// AliAnalysisManager to handle analysis.
//==============================================================================

#include <Riostream.h>
#include <TProcessID.h>

#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisSelector.h"

ClassImp(AliAnalysisSelector)

//______________________________________________________________________________
AliAnalysisSelector::AliAnalysisSelector()
                    :TSelector(), 
                     fInitialized(kFALSE), 
                     fAnalysis(NULL)
{
// Dummy ctor.
   fAnalysis = AliAnalysisManager::GetAnalysisManager();
   if (fAnalysis) fAnalysis->SetSelector(this);
}   

//______________________________________________________________________________
AliAnalysisSelector::AliAnalysisSelector(AliAnalysisManager *mgr)
                    :TSelector(),
                     fInitialized(kFALSE),
                     fAnalysis(mgr)
{
// Constructor. Called by AliAnalysisManager which registers itself on the
// selector running on the master.
   mgr->SetSelector(this);
}

//______________________________________________________________________________
AliAnalysisSelector::~AliAnalysisSelector()
{
// Dtor. The analysis manager object is sent in the input list and duplicated
// on the workers - it needs to be deleted (?)
//   if (fAnalysis) delete fAnalysis;
}

//______________________________________________________________________________
void AliAnalysisSelector::Init(TTree *tree)
{
// Called after Begin/SlaveBegin, assumes that fAnalysis is already initialized.
// Is Init called on workers in case of PROOF.
   if (!fAnalysis) {
      Error("Init", "Analysis manager NULL !");
      Abort("Cannot initialize without analysis manager. Aborting.");
      SetStatus(-1);
      return;
   }
   if (fAnalysis->GetDebugLevel()>1) {
      cout << "->AliAnalysisSelector->Init()" << endl;
   }   
   if (!tree) {
      Error("Init", "Input tree is NULL !");
      Abort("Cannot initialize without tree. Aborting.");
      SetStatus(-1);
      return;
   }
   fInitialized = fAnalysis->Init(tree);
   if (!fInitialized) {
      Error("Init", "Some error occured during analysis manager initialization. Aborting.");
      Abort("Error during AliAnalysisManager::Init()");
      SetStatus(-1);
      return;
   }   
   if (fAnalysis->GetDebugLevel()>1) {
      cout << "<-AliAnalysisSelector->Init()" << endl;
   }   
}

//______________________________________________________________________________
void AliAnalysisSelector::Begin(TTree *)
{
// Assembly the input list.
   RestoreAnalysisManager();
   if (fAnalysis && fAnalysis->GetDebugLevel()>1) {
      cout << "->AliAnalysisSelector->Begin: Analysis manager restored" << endl;
   }
}

//______________________________________________________________________________
void AliAnalysisSelector::SlaveBegin(TTree *tree)
{
// Called on each worker. We "unpack" analysis manager here and call InitAnalysis.
   RestoreAnalysisManager();
   if (fAnalysis) {
      if (fAnalysis->GetDebugLevel()>1) {
         cout << "->AliAnalysisSelector->SlaveBegin() after Restore" << endl;
      }   
      fAnalysis->SlaveBegin(tree);   
      if (fAnalysis->GetDebugLevel()>1) {
         cout << "<-AliAnalysisSelector->SlaveBegin()" << endl;
      }   
   }   
}

//______________________________________________________________________________
Bool_t AliAnalysisSelector::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.
   if (fAnalysis) return fAnalysis->Notify();
   return kFALSE;
}

//______________________________________________________________________________
Bool_t AliAnalysisSelector::Process(Long64_t entry)
{
// Event loop.
   static Int_t count = 0;
   count++;
   if (fAnalysis->GetDebugLevel() > 1) {
      cout << "->AliAnalysisSelector::Process()" << endl;
   }
   static Bool_t init=kTRUE;
   static Int_t nobjCount = 0;
   if(init) {
     nobjCount = TProcessID::GetObjectCount();
     init=kFALSE;
   }
   TProcessID::SetObjectCount(nobjCount);
   Int_t returnCode = fAnalysis->GetEntry(entry);
   if (returnCode <= 0) {
      cout << "Error retrieving event:" << entry << " Skipping ..." << endl;
      fAnalysis->CountEvent(1,0,1,0);
      // Try to skip file
      Abort("Bad stream to file. Trying next image.", kAbortFile);
      return kFALSE;
   } else {
      fAnalysis->ExecAnalysis();
      if (returnCode<100000000) fAnalysis->CountEvent(1,1,0,0);
   }   
   if (fAnalysis->GetDebugLevel() > 1) {
      cout << "<-AliAnalysisSelector::Process()" << endl;
   }   
   return kTRUE;
}   

//______________________________________________________________________________
void AliAnalysisSelector::RestoreAnalysisManager()
{
// Restores analysis manager from the input list.
   if (!fAnalysis) {
      TIter next(fInput);
      TObject *obj;
      while ((obj=next())) {
         if (obj->IsA() == AliAnalysisManager::Class()) {
            fAnalysis = (AliAnalysisManager*)obj;
            fAnalysis->SetSelector(this);
            if (fAnalysis->GetDebugLevel()>1) {
               cout << "->AliAnalysisSelector->RestoreAnalysisManager: Analysis manager restored" << endl;
            }   
            break;
         }
      }
      if (!fAnalysis) {
         Error("SlaveBegin", "Analysis manager not found in the input list");
         return;
      }   
   }
}

//______________________________________________________________________________
void AliAnalysisSelector::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
   if (fStatus == -1) return;  // TSelector won't abort...
   if (fAnalysis->GetAnalysisType() == AliAnalysisManager::kMixingAnalysis) return;
   if (fAnalysis->GetDebugLevel() > 1) {
      cout << "->AliAnalysisSelector::SlaveTerminate()" << endl;
   }   
   fAnalysis->PackOutput(fOutput);
   if (fAnalysis->GetDebugLevel() > 1) {
      cout << "<-AliAnalysisSelector::SlaveTerminate()" << endl;
   }   
}  

//______________________________________________________________________________
void AliAnalysisSelector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
   if (fStatus == -1) return;  // TSelector won't abort...
   if (!fAnalysis) {
      Error("Terminate","AliAnalysisSelector::Terminate: No analysis manager!!!");
      return;
   }   
   // No Terminate() in case of event mixing
   if (fAnalysis->GetAnalysisType() == AliAnalysisManager::kMixingAnalysis) return;
   if (fAnalysis->GetDebugLevel() > 1) {
      cout << "->AliAnalysisSelector::Terminate()" << endl;
   }   
   fAnalysis->UnpackOutput(fOutput);
   fAnalysis->Terminate();   
   if (fAnalysis->GetDebugLevel() > 1) {
      cout << "<-AliAnalysisSelector::Terminate()" << endl;
   }   
}
