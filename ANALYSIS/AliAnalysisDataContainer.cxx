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
//   AliAnalysysDataContainer - Container of data of arbitrary type deriving
//      from TObject used for analysis. A container must be connected to the 
//      output data slot of a single analysis task (producer) , but also as 
//      input slot for possibly several other tasks (consumers). The connected 
//      slots must enforce the same data type as the container (or a derived type).
//      A container becomes the owner of the contained data once this was produced.
//
// Containers should be defined by the analysis module using:
//
//   AliAnalysisModule::AddContainer(const char *name, TClass *type);
//
// A container should be connected to a producer:

//   AliAnalysisModule::ConnectOutput(AliAnalysisTask *task,
//                                    AliAnalysisDataContainer *cont)
// and to its consumers:
//
//   AliAnalysisModule::ConnectInput(AliAnalysisTask *task, Int_t islot,
//                                   AliAnalysisDataContainer *cont)
//
// The container will create an implicit connection between the producer task 
// and all consumers, which will become sub-tasks of the producer.
//
//==============================================================================

#include "TClass.h"
#include "TTree.h"
#include "TFile.h"
#include "AliLog.h"

#include "AliAnalysisDataContainer.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisTask.h"

ClassImp(AliAnalysisDataContainer)

//______________________________________________________________________________
AliAnalysisDataContainer::AliAnalysisDataContainer() : TNamed(),
                          fDataReady(kFALSE),
                          fOwnedData(kFALSE),
                          fFile(NULL),
                          fData(NULL),
                          fType(NULL),
                          fProducer(NULL),
                          fConsumers(NULL)
{
// Default ctor.
}

//______________________________________________________________________________
AliAnalysisDataContainer::AliAnalysisDataContainer(const char *name, TClass *type)
                         :TNamed(name,""),
                          fDataReady(kFALSE),
                          fOwnedData(kTRUE),
                          fFile(NULL),
                          fData(NULL),
                          fType(type),
                          fProducer(NULL),
                          fConsumers(NULL)
{
// Normal constructor.
}

//______________________________________________________________________________
AliAnalysisDataContainer::AliAnalysisDataContainer(const AliAnalysisDataContainer &cont)
                         :TNamed(cont),
                          fDataReady(cont.fDataReady),
                          fOwnedData(kFALSE),
                          fFile(cont.fFile),
                          fData(cont.fData),
                          fType(cont.fType),
                          fProducer(cont.fProducer),
                          fConsumers(NULL)
{
// Copy ctor.
   if (cont.fConsumers) {
      fConsumers = new TObjArray(2);
      Int_t ncons = cont.fConsumers->GetEntriesFast();
      for (Int_t i=0; i<ncons; i++) fConsumers->Add(cont.fConsumers->At(i));
   }   
}

//______________________________________________________________________________
AliAnalysisDataContainer::~AliAnalysisDataContainer()
{
// Destructor. Deletes data ! (What happens if data is a container ???)
   if (fData && fOwnedData) delete fData;
   if (fConsumers) delete fConsumers;
   if (fFile) {
      fFile->Close();
      delete fFile;
   }   
}

//______________________________________________________________________________
AliAnalysisDataContainer &AliAnalysisDataContainer::operator=(const AliAnalysisDataContainer &cont)
{
// Assignment.
   if (&cont != this) {
      TNamed::operator=(cont);
      fDataReady = cont.fDataReady;
      fOwnedData = kFALSE;  // !!! Data owned by cont.
      fFile = cont.fFile;
      fData = cont.fData;
      fType = cont.fType;
      fProducer = cont.fProducer;
      if (cont.fConsumers) {
         fConsumers = new TObjArray(2);
         Int_t ncons = cont.fConsumers->GetEntriesFast();
         for (Int_t i=0; i<ncons; i++) fConsumers->Add(cont.fConsumers->At(i));
      }   
   }   
   return *this;
}      

//______________________________________________________________________________
Bool_t AliAnalysisDataContainer::SetData(TObject *data, Option_t *)
{
// Set the data as READY only if it was published by the producer.
   // If there is no producer declared, this is a top level container.
   AliAnalysisTask *task;
   Bool_t init = kFALSE;
   Int_t i, nc;
   if (!fProducer) {
      if (data != fData) init = kTRUE;
      fData = data;
      fDataReady = kTRUE;
      if (fConsumers) {
         nc = fConsumers->GetEntriesFast();
         for (i=0; i<nc; i++) {
            task = (AliAnalysisTask*)fConsumers->At(i);
            task->CheckNotify(init);
         }
      }      
      return kTRUE;
   } 
   // Check if it is the producer who published the data     
   if (fProducer->GetPublishedData()==data) {
      fData = data;
      fDataReady = kTRUE;
      if (fConsumers) {
         nc = fConsumers->GetEntriesFast();
         for (i=0; i<nc; i++) {
            task = (AliAnalysisTask*)fConsumers->At(i);
            task->CheckNotify();
         }
      }      
      return kTRUE;   
   } else {
      AliWarning(Form("Data for container %s can be published only by producer task %s", 
                 GetName(), fProducer->GetName()));   
      return kFALSE;           
   }              
}

//______________________________________________________________________________
void AliAnalysisDataContainer::OpenFile(const char *name, Option_t *option)
{
// Data will be written to this file at the end of processing.
// Option represent the way the file is accessed: NEW, APPEND, ...
   if (fFile) {
      fFile->Close();
      delete fFile;
   }
   fFile =  new TFile(name, option);
   if (fFile->IsZombie()) {
      AliError(Form("Cannot open file %s with option %s",name,option));
      fFile = 0;
   }   
}   

//______________________________________________________________________________
void AliAnalysisDataContainer::WriteData()
{
// Write data to the file.
   if (fFile) {
      TDirectory *cursav = gDirectory;
      fFile->cd();
      fData->Write();
//      fFile->Write();
      if (cursav) cursav->cd();
   }   
}

//______________________________________________________________________________
void AliAnalysisDataContainer::GetEntry(Long64_t ientry)
{
// If data is ready and derives from TTree or from TBranch, this will get the
// requested entry in memory if not already loaded.
   if (!fDataReady) return;
   Bool_t istree = fType->InheritsFrom(TTree::Class());
   if (istree) {
      TTree *tree = (TTree*)fData;
      if (tree->GetReadEntry() != ientry) tree->GetEntry(ientry);
      return;
   }   
   Bool_t isbranch = fType->InheritsFrom(TBranch::Class());
   if (isbranch) {
      TBranch *branch = (TBranch*)fData;
      if (branch->GetReadEntry() != ientry) branch->GetEntry(ientry);
      return;
   }   
}   

//______________________________________________________________________________
void AliAnalysisDataContainer::SetProducer(AliAnalysisTask *prod, Int_t islot)
{
// Set the producer of data. The slot number is required for data type checking.
   if (fProducer) {
      AliWarning(Form("Data container %s already has a producer: %s",
                 GetName(),fProducer->GetName()));
   } 
   if (fDataReady) {
      AliError(Form("%s container contains data - cannot change producer!", GetName()));
      return;
   }   
   AliAnalysisDataSlot *slot = prod->GetOutputSlot(islot);
   if (!slot) {
      AliError(Form("Producer task %s does not have an output #%i", prod->GetName(),islot));
      return;
   }   
   if (!slot->GetType()->InheritsFrom(fType)) {
      AliError(Form("Data type %s for output slot %i of task %s does not match container type %s", 
                     slot->GetType()->GetName(),islot,prod->GetName(),fType->GetName()));
      return;
   }   
   
   fProducer = prod;
   // Add all consumers as daughter tasks
   TIter next(fConsumers);
   AliAnalysisTask *cons;
   while ((cons=(AliAnalysisTask*)next())) {
      if (!prod->GetListOfTasks()->FindObject(cons)) prod->Add(cons);
   }   
}   

//______________________________________________________________________________
void AliAnalysisDataContainer::AddConsumer(AliAnalysisTask *consumer, Int_t islot)
{
// Add a consumer for contained data;
   AliAnalysisDataSlot *slot = consumer->GetInputSlot(islot);
   if (!slot) {
      AliError(Form("Consumer task %s does not have an input #%i", consumer->GetName(),islot));
      return;
   }   
   if (!slot->GetType()->InheritsFrom(fType)) {
      AliError(Form("Data type %s for input slot %i of task %s does not match container type %s", 
                     slot->GetType()->GetName(),islot,consumer->GetName(),fType->GetName()));
      return;
   }   

   if (!fConsumers) fConsumers = new TObjArray(2);   
   fConsumers->Add(consumer);
   // Add the consumer task to the list of task of the producer
   if (fProducer && !fProducer->GetListOfTasks()->FindObject(consumer)) 
      fProducer->Add(consumer);
}      

//______________________________________________________________________________
Bool_t AliAnalysisDataContainer::ClientsExecuted() const
{
// Check if all client tasks have executed.
   TIter next(fConsumers);
   AliAnalysisTask *task;
   while ((task=(AliAnalysisTask*)next())) {
      if (!task->HasExecuted()) return kFALSE;
   }
   return kTRUE;
}   

//______________________________________________________________________________
void AliAnalysisDataContainer::DeleteData()
{
// Delete data if not needed anymore.
   if (!fDataReady || !ClientsExecuted()) {
      AliWarning(Form("Data not ready or not all clients of container %s executed. Data not deleted.", GetName()));
      return;
   }
   if (!fOwnedData) {
      AliWarning(Form("Data not owned by container %s. Not deleted.", GetName()));
      return;
   }
   delete fData;
   fData = 0;
   fDataReady = kFALSE;
}   
      
//______________________________________________________________________________
void AliAnalysisDataContainer::PrintContainer(Option_t *option, Int_t indent) const
{
// Print info about this container.
   TString ind;
   for (Int_t i=0; i<indent; i++) ind += " ";
   TString opt(option);
   opt.ToLower();
   Bool_t dep = (opt.Contains("dep"))?kTRUE:kFALSE;
   if (!dep) {
      printf("%s\n", Form("%sContainer: %s  type: %s", ind.Data(), GetName(), fType->GetName()));
      if (fProducer) 
         printf("%s\n", Form("%s = Data producer: task %s",ind.Data(),fProducer->GetName()));
      else
         printf("%s\n", Form("%s= No data producer"));
      printf("%s", Form("%s = Consumer tasks: "));
      if (!fConsumers || !fConsumers->GetEntriesFast()) printf("-none-\n");
      else printf("\n");
   }   
   TIter next(fConsumers);
   AliAnalysisTask *task;
   while ((task=(AliAnalysisTask*)next())) task->PrintTask(option, indent+3);
}   
