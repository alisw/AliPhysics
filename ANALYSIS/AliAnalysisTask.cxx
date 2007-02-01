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
//   AliAnalysysTask - Class representing a basic analysis task. Any
// user-defined task should derive from it and implement the Exec() virtual
// method.
//==============================================================================
//
// A specific user analysis task have to derive from this class. The list of
// specific input and output slots have to be defined in the derived class ctor:
//
//   UserTask::UserTask(name, title)
//   {
//      DefineInput(0, TTree::Class());
//      DefineInput(1, TH1::Class());
//      ...
//      DefineOutput(0, TTree::Class());
//      DefineOutput(1, MyObject::Class());
//      ...
//   }
//
// An existing data contaner (AliAnalysisDataContainer) can be connected to the
// input/output slots of an analysis task. Containers should not be defined and
// connected by the derived analysis task, but from the level of AliAnalysisManager:
//
//   AliAnalysisManager::ConnectInput(AliAnalysisTask *task, Int_t islot,
//                                   AliAnalysisDataContainer *cont)
//   AliAnalysisManager::ConnectOutput(AliAnalysisTask *task, Int_t islot,
//                                    AliAnalysisDataContainer *cont)
// To connect a slot to a data container, the data types declared by both must
// match.
//
// The method Init will be called once per session at the moment when the data is
// available at all input slots. 
// The method Init() has to be overloaded by the derived class in order to:
// 1. Define objects that should be created only once per session (e.g. output
//    histograms)
// 2. Set the branch address or connect to a branch address in case the input
//    slots are connected to trees.
//
// Example:
// MyAnalysisTask::Init(Option_t *)
// {
//    if (!fHist1) fHist1 = new TH1F("h1", ....);
//    if (!fHist2) fHist2 = new TH1F("h2", ....);
//    fESD = GetBranchAddress(islot=0, "ESD");
//    if (!fESD) SetBranchAddress(islot=0, "ESD", &fESD);
// }
//
// The method Terminate() will be called by the framework once at the end of
// data processing. Overload this if needed.
//
//==============================================================================

#include <Riostream.h>
#include <TClass.h>
#include <TDirectory.h>

//#include "AliLog.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"

ClassImp(AliAnalysisTask)

//______________________________________________________________________________
AliAnalysisTask::AliAnalysisTask()
                :fReady(kFALSE),
                 fInitialized(kFALSE),
                 fNinputs(0),
                 fNoutputs(0),
                 fOutputReady(NULL),
                 fPublishedData(NULL),
                 fInputs(NULL),
                 fOutputs(NULL)
{
// Default constructor.
}

//______________________________________________________________________________
AliAnalysisTask::AliAnalysisTask(const char *name, const char *title)
                :TTask(name,title),
                 fReady(kFALSE),
                 fInitialized(kFALSE),
                 fNinputs(0),
                 fNoutputs(0),
                 fOutputReady(NULL),
                 fPublishedData(NULL),
                 fInputs(NULL),
                 fOutputs(NULL)                 
{
// Constructor.
   fInputs      = new TObjArray(2);
   fOutputs     = new TObjArray(2);
}

//______________________________________________________________________________
AliAnalysisTask::AliAnalysisTask(const AliAnalysisTask &task)
                :TTask(task),
                 fReady(task.fReady),
                 fInitialized(task.fInitialized),
                 fNinputs(task.fNinputs),
                 fNoutputs(task.fNoutputs),                 
                 fOutputReady(NULL),
                 fPublishedData(NULL),
                 fInputs(NULL),
                 fOutputs(NULL)                 
{
// Copy ctor.
   fInputs      = new TObjArray((fNinputs)?fNinputs:2);
   fOutputs     = new TObjArray((fNoutputs)?fNoutputs:2);
   fPublishedData = 0;
   Int_t i;
   for (i=0; i<fNinputs; i++) fInputs->AddAt(task.GetInputSlot(i),i);
   fOutputReady = new Bool_t[(fNoutputs)?fNoutputs:2];
   for (i=0; i<fNoutputs; i++) {
      fOutputReady[i] = IsOutputReady(i);
      fOutputs->AddAt(task.GetOutputSlot(i),i);
   }   
}

//______________________________________________________________________________
AliAnalysisTask::~AliAnalysisTask()
{
// Dtor.
   if (fInputs)  {fInputs->Delete(); delete fInputs;}
   if (fOutputs) {fOutputs->Delete(); delete fOutputs;}
}   
  
//______________________________________________________________________________
AliAnalysisTask& AliAnalysisTask::operator=(const AliAnalysisTask& task)
{
// Assignment
   if (&task == this) return *this;
   TTask::operator=(task);
   fReady       = task.IsReady();
   fInitialized = task.IsInitialized();
   fNinputs     = task.GetNinputs();
   fNoutputs    = task.GetNoutputs();
   fInputs      = new TObjArray((fNinputs)?fNinputs:2);
   fOutputs     = new TObjArray((fNoutputs)?fNoutputs:2);
   fPublishedData = 0;
   Int_t i;
   for (i=0; i<fNinputs; i++) fInputs->AddAt(new AliAnalysisDataSlot(*task.GetInputSlot(i)),i);
   fOutputReady = new Bool_t[(fNoutputs)?fNoutputs:2];
   for (i=0; i<fNoutputs; i++) {
      fOutputReady[i] = IsOutputReady(i);
      fOutputs->AddAt(new AliAnalysisDataSlot(*task.GetOutputSlot(i)),i);
   }         
   return *this;
}

//______________________________________________________________________________
Bool_t AliAnalysisTask::AreSlotsConnected()
{
// Check if all input/output slots are connected. If this is the case fReady=true
   fReady = kFALSE;
   if (!fNinputs || !fNoutputs) return kFALSE;
   Int_t i;
   AliAnalysisDataSlot *slot;
   for (i=0; i<fNinputs; i++) {
      slot = (AliAnalysisDataSlot*)fInputs->At(i);
      if (!slot) {
	cout<<"Input slot "<<i<<" of task "<<GetName()<<" not defined !"<<endl;
	//AliError(Form("Input slot %i of task %s not defined !",i,GetName()));
         return kFALSE;
      }   
      if (!slot->IsConnected()) return kFALSE;
   }   
   for (i=0; i<fNoutputs; i++) {
      slot = (AliAnalysisDataSlot*)fOutputs->At(i);
      if (!slot) {
	cout<<"Output slot "<<i<<" of task "<<GetName()<<" not defined !"<<endl;
	//AliError(Form("Output slot %i of task %s not defined !",i,GetName()));
         return kFALSE;
      }   
      if (!slot->IsConnected()) return kFALSE;
   } 
   fReady = kTRUE;  
   return kTRUE;
}

//______________________________________________________________________________
void AliAnalysisTask::CheckNotify(Bool_t init)
{
// Check if data is available from all inputs. Change the status of the task
// accordingly. This method is called automatically for all tasks connected
// to a container where the data was published.
   if (init) fInitialized = kFALSE;
   for (Int_t islot=0; islot<fNinputs; islot++) {
      if (!GetInputData(islot)) {
         SetActive(kFALSE);
         return;
      }   
   }   
   SetActive(kTRUE);
   if (fInitialized) return;
   TDirectory *cursav = gDirectory;
   Init();
   if (cursav) cursav->cd();
   fInitialized = kTRUE;
}

//______________________________________________________________________________
Bool_t AliAnalysisTask::ConnectInput(Int_t islot, AliAnalysisDataContainer *cont)
{
// Connect an input slot to a data container.
   AliAnalysisDataSlot *input = GetInputSlot(islot);
   if (!input) {
     cout<<"Input slot "<<islot<<" not defined for analysis task "<<GetName()<<endl;
     //AliError(Form("Input slot %i not defined for analysis task %s", islot, GetName()));
      return kFALSE;
   }
   // Check type matching          
   if (!input->GetType()->InheritsFrom(cont->GetType())) {
     cout<<"Data type "<<input->GetType()->GetName()<<" for input "<<islot<<" of task "<<GetName()<<" not matching container "<<cont->GetName()<<" of type "<<cont->GetType()->GetName()<<endl;
     //AliError(Form("Data type %s for input %i of task %s not matching container %s of type %s",input->GetType()->GetName(), islot, GetName(), cont->GetName(), cont->GetType()->GetName()));
      return kFALSE;
   }  
   // Connect the slot to the container as input          
   if (!input->ConnectContainer(cont)) return kFALSE;
   // Add this to the list of container consumers
   cont->AddConsumer(this, islot);
   AreSlotsConnected();
   return kTRUE;
}   

//______________________________________________________________________________
Bool_t AliAnalysisTask::ConnectOutput(Int_t islot, AliAnalysisDataContainer *cont)
{
// Connect an output slot to a data container.
   AliAnalysisDataSlot *output = GetOutputSlot(islot);
   if (!output) {
     cout<<"Output slot "<<islot<<" not defined for analysis task "<<GetName()<<endl;
     //AliError(Form("Output slot %i not defined for analysis task %s", islot, GetName()));
      return kFALSE;
   }
   // Check type matching          
   if (!output->GetType()->InheritsFrom(cont->GetType())) {
     cout<<"Data type "<<output->GetType()->GetName()<<" for output "<<islot<<" of task "<<GetName()<<" not matching container "<<cont->GetName()<<" of type "<<cont->GetType()->GetName()<<endl;
     //AliError(Form("Data type %s for output %i of task %s not matching container %s of type %s",output->GetType()->GetName(), islot, GetName(), cont->GetName(), cont->GetType()->GetName()));
      return kFALSE;
   }            
   // Connect the slot to the container as output         
   if (!output->ConnectContainer(cont)) return kFALSE;
   // Declare this as the data producer
   cont->SetProducer(this, islot);
   AreSlotsConnected();
   return kTRUE;
}   

//______________________________________________________________________________
void AliAnalysisTask::DefineInput(Int_t islot, TClass *type)
{
// Define an input slot and its type.
   AliAnalysisDataSlot *input = new AliAnalysisDataSlot(type, this);
   if (fNinputs<islot+1) fNinputs = islot+1;
   fInputs->AddAtAndExpand(input, islot);
}

//______________________________________________________________________________
void AliAnalysisTask::DefineOutput(Int_t islot, TClass *type)
{
// Define an output slot and its type.
   if (islot<0) {
     cout<<"Cannot define negative output slot number for task "<<GetName()<<endl;
     //AliError(Form("Cannot define negative output slot number for task %s", GetName()));
      return;
   }   
   AliAnalysisDataSlot *output = new AliAnalysisDataSlot(type, this);
   if (fNoutputs<islot+1) {
      fNoutputs = islot+1;
      if (fOutputReady) delete [] fOutputReady;
      fOutputReady = new Bool_t[fNoutputs];
      memset(fOutputReady, 0, fNoutputs*sizeof(Bool_t));
   } 
   fOutputs->AddAtAndExpand(output, islot);
}

//______________________________________________________________________________
TClass *AliAnalysisTask::GetInputType(Int_t islot) const
{
// Retreive type of a given input slot.
   AliAnalysisDataSlot *input = GetInputSlot(islot);
   if (!input) {
     cout<<"Input slot "<<islot<<" not defined for analysis task "<<GetName()<<endl;
     //AliError(Form("Input slot %i not defined for analysis task %s", islot, GetName()));
      return NULL;
   }
   return (input->GetType());
}

//______________________________________________________________________________
TClass *AliAnalysisTask::GetOutputType(Int_t islot) const
{
// Retreive type of a given output slot.
   AliAnalysisDataSlot *output = GetOutputSlot(islot);
   if (!output) {
     cout<<"Output slot "<<islot<<" not defined for analysis task "<<GetName()<<endl;
     //AliError(Form("Output slot %i not defined for analysis task %s", islot, GetName()));
      return NULL;
   }
   return (output->GetType());
}

//______________________________________________________________________________
TObject *AliAnalysisTask::GetInputData(Int_t islot) const
{
// Retreive input data for a slot if ready. Normally called by Exec() and
// the object has to be statically cast to the appropriate type.
   AliAnalysisDataSlot *input = GetInputSlot(islot);
   if (!input) {
     cout<<"Input slot "<<islot<<" not defined for analysis task "<<GetName()<<endl;
     //AliError(Form("Input slot %i not defined for analysis task %s", islot, GetName()));
      return NULL;
   }
   return (input->GetData()); 
}

//______________________________________________________________________________
char *AliAnalysisTask::GetBranchAddress(Int_t islot, const char *branch) const
{
// Check if a branch with a given name from the specified input is connected
// to some address. Call this in Init() before trying to call SetBranchAddress()
// since the adress may be set by other task.
   return (char *)GetInputSlot(islot)->GetBranchAddress(branch);
}

//______________________________________________________________________________
Bool_t AliAnalysisTask::SetBranchAddress(Int_t islot, const char *branch, void *address) const
{
// Connect an object address to a branch of the specified input.
   return GetInputSlot(islot)->SetBranchAddress(branch, address);
}   

//______________________________________________________________________________
void AliAnalysisTask::Init(Option_t *)
{
// Branch address initialization.
}

//______________________________________________________________________________
void AliAnalysisTask::Terminate(Option_t *)
{
// Method called by the framework at the end of data processing.
}

//______________________________________________________________________________
void AliAnalysisTask::OpenFile(Int_t iout, const char *name, Option_t *option) const
{
// Set data at output iout to be written in the specified file.
   GetOutputSlot(iout)->GetContainer()->OpenFile(name, option);
}   

//______________________________________________________________________________
Bool_t AliAnalysisTask::PostData(Int_t iout, TObject *data, Option_t *option)
{
// Post output data for a given ouput slot in the corresponding data container.
// Published data becomes owned by the data container.
// If option is specified, the container connected to the output slot must have
// an associated file name defined. The option represents the method to open the file.
   fPublishedData = 0;
   AliAnalysisDataSlot *output = GetOutputSlot(iout);
   if (!output) {
     cout<<"Output slot "<<iout<<" not defined for analysis task "<<GetName()<<endl;
     //AliError(Form("Output slot %i not defined for analysis task %s", iout, GetName()));
      return kFALSE;
   }
   if (!output->IsConnected()) {
     cout<<"Output slot "<<iout<<" of analysis task "<<GetName()<<" not connected to any data container"<<endl;
     //AliError(Form("Output slot %i of analysis task %s not connected to any data container", iout, GetName()));
      return kFALSE;
   }
   if (!fOutputReady) {
      fOutputReady = new Bool_t[fNoutputs];
      memset(fOutputReady, 0, fNoutputs*sizeof(Bool_t));
   }   
   fOutputReady[iout] = kTRUE;
   fPublishedData = data;
   return (output->GetContainer()->SetData(data, option));
}

//______________________________________________________________________________
void AliAnalysisTask::SetUsed(Bool_t flag)
{
// Set 'used' flag recursively to task and all daughter tasks.
   if (TestBit(kTaskUsed)==flag) return;
   TObject::SetBit(kTaskUsed,flag);
   Int_t nd = fTasks->GetSize();
   AliAnalysisTask *task;
   for (Int_t i=0; i<nd; i++) {
      task = (AliAnalysisTask*)fTasks->At(i);
      task->SetUsed(flag);
   }
}   

//______________________________________________________________________________
Bool_t AliAnalysisTask::CheckCircularDeps()
{
// Check for illegal circular dependencies, e.g. a daughter task should not have
// a hierarchical parent as subtask.
   if (IsChecked()) return kTRUE;
   SetChecked();
   TList *tasks = GetListOfTasks();
   Int_t ntasks = tasks->GetSize();
   AliAnalysisTask *task;
   for (Int_t i=0; i<ntasks; i++) {
      task = (AliAnalysisTask*)tasks->At(i);
      if (task->CheckCircularDeps()) return kTRUE;
   }
   SetChecked(kFALSE);
   return kFALSE;
}   
   
//______________________________________________________________________________
void AliAnalysisTask::PrintTask(Option_t *option, Int_t indent) const
{
// Print task info.
   AliAnalysisTask *thistask = (AliAnalysisTask*)this;
   TString opt(option);
   opt.ToLower();
   Bool_t dep = (opt.Contains("dep"))?kTRUE:kFALSE;
   TString ind;
   Int_t islot;
   AliAnalysisDataContainer *cont;
   for (Int_t i=0; i<indent; i++) ind += " ";
   if (!dep || (dep && IsChecked())) {
      printf("%s\n", Form("%stask: %s  ACTIVE=%i", ind.Data(), GetName(),IsActive()));
      if (dep) thistask->SetChecked(kFALSE);
      else {
         for (islot=0; islot<fNinputs; islot++) {
            printf("%s", Form("%s   INPUT #%i: %s <- ",ind.Data(),islot, GetInputType(islot)->GetName()));
            cont = GetInputSlot(islot)->GetContainer();
            if (cont) printf(" [%s]\n", cont->GetName());
            else printf(" [NO CONTAINER]\n");
         }
         for (islot=0; islot<fNoutputs; islot++) {
            printf("%s", Form("%s   OUTPUT #%i: %s -> ",ind.Data(),islot, GetOutputType(islot)->GetName()));
            cont = GetOutputSlot(islot)->GetContainer();
            if (cont) printf(" [%s]\n", cont->GetName());
            else printf(" [NO CONTAINER]\n");
         }            
      }
   }
   PrintContainers(option, indent+3);
}      

//______________________________________________________________________________
void AliAnalysisTask::PrintContainers(Option_t *option, Int_t indent) const
{
// Print containers info.
   AliAnalysisDataContainer *cont;
   TString ind;
   for (Int_t i=0; i<indent; i++) ind += " ";
   Int_t islot;
   for (islot=0; islot<fNoutputs; islot++) {
      cont = GetOutputSlot(islot)->GetContainer();
      cont->PrintContainer(option, indent);
   }   
}
