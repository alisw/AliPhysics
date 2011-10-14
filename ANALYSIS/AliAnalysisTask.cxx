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
// The method ConnectInputData() has to be overloaded by the derived class in order to
// set the branch address or connect to a branch address in case the input
// slots are connected to trees.
// Example:
// MyAnalysisTask::ConnectInputData(Option_t *)
// {
//  // One should first check if the branch address was taken by some other task
//    char ** address = (char **)GetBranchAddress(0, "ESD");
//    if (address) {
//      fESD = (AliESD*)(*address);
//    } else {
//      fESD = new AliESD();
//      SetBranchAddress(0, "ESD", &fESD);
//    }
// }
// 
// The method LocalInit() may be implemented to call locally (on the client)
// all initialization methods of the class. It is not mandatory and was created
// in order to minimize the complexity and readability of the analysis macro.
// DO NOT create in this method the histigrams or task output objects that will
// go in the task output containers. Use CreateOutputObjects for that.
//
// The method CreateOutputObjects() has to be implemented an will contain the
// objects that should be created only once per session (e.g. output
// histograms)
//
// void MyAnalysisTask::CreateOutputObjects()
//{
  // create histograms 
//  fhPt = new TH1F("fhPt","This is the Pt distribution",15,0.1,3.1);
//  fhPt->SetStats(kTRUE);
//  fhPt->GetXaxis()->SetTitle("P_{T} [GeV]");
//  fhPt->GetYaxis()->SetTitle("#frac{dN}{dP_{T}}");
//  fhPt->GetXaxis()->SetTitleColor(1);
//  fhPt->SetMarkerStyle(kFullCircle);
// }
//
// The method Terminate() will be called by the framework once at the end of
// data processing. Overload this if needed. DO NOT ASSUME that the pointers
// to histograms defined in  CreateOutputObjects() are valid, since this is
// not true in case of PROOF. Restore the pointer values like:
//
//void MyAnalysisTask::Terminate(Option_t *) 
//{
//  fhPt = (TH1F*)GetOutputData(0);
// ...
//}

//
//==============================================================================

#include <Riostream.h>
#include <TFile.h>
#include <TClass.h>
#include <TCollection.h>
#include <TTree.h>
#include <TROOT.h>

#include "AliAnalysisTask.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"

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
                 fOutputs(NULL),
                 fBranchNames()
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
                 fOutputs(NULL),
                 fBranchNames()                 
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
                 fOutputs(NULL),
                 fBranchNames(task.fBranchNames)
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
   if (fTasks) fTasks->Clear();
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
   fBranchNames = task.fBranchNames;
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
	      Error("AreSlotsConnected", "Input slot %d of task %s not defined !",i,GetName());
         return kFALSE;
      }   
      if (!slot->IsConnected()) return kFALSE;
   }   
   for (i=0; i<fNoutputs; i++) {
      slot = (AliAnalysisDataSlot*)fOutputs->At(i);
      if (!slot) {
         Error("AreSlotsConnected", "Output slot %d of task %s not defined !",i,GetName());
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
   Bool_t single_shot = IsPostEventLoop();
   AliAnalysisDataContainer *cinput;
   for (Int_t islot=0; islot<fNinputs; islot++) {
      cinput = GetInputSlot(islot)->GetContainer();
      if (!cinput->GetData() || (single_shot && !cinput->IsPostEventLoop())) {
         SetActive(kFALSE);
         return;
      }   
   }   
   SetActive(kTRUE);
   if (fInitialized) return;
   TDirectory *cursav = gDirectory;
   ConnectInputData();
   if (cursav) cursav->cd();
   fInitialized = kTRUE;
}

//______________________________________________________________________________
Bool_t AliAnalysisTask::CheckPostData() const
{
// Checks if data was posted to all outputs defined by the task. If task does
// not have output slots this returns always kTRUE.
   Bool_t dataPosted = kTRUE;
   AliAnalysisDataContainer *coutput;
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   for (Int_t islot=0; islot<fNoutputs; islot++) {
      coutput = GetOutputSlot(islot)->GetContainer();
      if (!mgr->GetOutputs()->FindObject(coutput)) continue;
      if (!coutput->GetData()) {
         Error("CheckPostData", "Data not posted for slot #%d of task %s (%s)", 
               islot, GetName(), ClassName());
         dataPosted = kFALSE;
      }   
   }
   CheckOwnership();
   return dataPosted;
}

//______________________________________________________________________________
Bool_t AliAnalysisTask::CheckOwnership() const
{
// Check ownership of containers posted on output slots (1 level only)
   TObject *outdata;
   for (Int_t islot=0; islot<fNoutputs; islot++) {
      outdata = GetOutputData(islot);
      if (outdata && outdata->InheritsFrom(TCollection::Class())) {
         TCollection *coll = (TCollection*)outdata;
         if (!coll->IsOwner()) {
            Error("CheckOwnership","####### IMPORTANT! ####### \n\n\n\
                Task %s (%s) posts a container that is not owner at output #%d. This may apply for other embedded containers. \n\n\
                ####### FIX YOUR CODE, THIS WILL PRODUCE A FATAL ERROR IN FUTURE! ##########", GetName(), ClassName(), islot);
             return kFALSE;   
         }
      }
   }
   return kTRUE;
}
      
//______________________________________________________________________________
Bool_t AliAnalysisTask::ConnectInput(Int_t islot, AliAnalysisDataContainer *cont)
{
// Connect an input slot to a data container.
   AliAnalysisDataSlot *input = GetInputSlot(islot);
   if (!input) {
      Error("ConnectInput","Input slot %i not defined for analysis task %s", islot, GetName());
      return kFALSE;
   }
   // Check type matching          
   if (!input->GetType()->InheritsFrom(cont->GetType())) {
      Error("ConnectInput","Data type %s for input %i of task %s not matching container %s of type %s",input->GetType()->GetName(), islot, GetName(), cont->GetName(), cont->GetType()->GetName());
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
      Error("ConnectOutput","Output slot %i not defined for analysis task %s", islot, GetName());
      return kFALSE;
   }
   // Check type matching          
   if (!output->GetType()->InheritsFrom(cont->GetType())) {
      Error("ConnectOutput","Data type %s for output %i of task %s not matching container %s of type %s",output->GetType()->GetName(), islot, GetName(), cont->GetName(), cont->GetType()->GetName());
      return kFALSE;
   }            
   // Connect the slot to the container as output         
   if (!output->ConnectContainer(cont)) return kFALSE;
   // Set event loop type the same as for the task
   cont->SetPostEventLoop(IsPostEventLoop());
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
      Error("GetInputType","Input slot %d not defined for analysis task %s", islot, GetName());
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
      Error("GetOutputType","Output slot %d not defined for analysis task %s", islot, GetName());
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
      Error("GetInputData","Input slot %d not defined for analysis task %s", islot, GetName());
      return NULL;
   }
   return (input->GetData()); 
}

//______________________________________________________________________________
TObject *AliAnalysisTask::GetOutputData(Int_t islot) const
{
// Retreive output data for a slot. Normally called in UserTask::Terminate to
// get a valid pointer to data even in case of Proof.
   AliAnalysisDataSlot *output = GetOutputSlot(islot);
   if (!output) {
      Error("GetOutputData","Input slot %d not defined for analysis task %s", islot, GetName());
      return NULL;
   }
   return (output->GetData()); 
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
void AliAnalysisTask::EnableBranch(Int_t islot, const char *bname) const
{
// Call this in ConnectInputData() to enable only the branches needed by this 
// task. "*" will enable everything.
   AliAnalysisDataSlot *input = GetInputSlot(islot);
   if (!input || !input->GetType()->InheritsFrom(TTree::Class())) {
      Error("EnableBranch", "Wrong slot type #%d for task %s: not TTree-derived type", islot, GetName());
      return;
   }   
   TTree *tree = (TTree*)input->GetData();
   if (!strcmp(bname, "*")) {
      tree->SetBranchStatus("*",1);
      return;
   }
   AliAnalysisDataSlot::EnableBranch(bname, tree);
}

//______________________________________________________________________________
void AliAnalysisTask::FinishTaskOutput()
{
// Optional method that is called in SlaveTerminate phase. 
// Used for calling aditional methods just after the last event was processed ON
// THE WORKING NODE. The call is made also in local case.
// Do NOT delete output objects here since they will have to be sent for 
// merging in PROOF mode - use class destructor for cleanup.
}
      
//______________________________________________________________________________
void AliAnalysisTask::ConnectInputData(Option_t *)
{
// Overload and connect your branches here.
}

//______________________________________________________________________________
void AliAnalysisTask::LocalInit()
{
// The method LocalInit() may be implemented to call locally (on the client)
// all initialization methods of the class. It is not mandatory and was created
// in order to minimize the complexity and readability of the analysis macro.
// DO NOT create in this method the histigrams or task output objects that will
// go in the task output containers. Use CreateOutputObjects for that.
}

//______________________________________________________________________________
void AliAnalysisTask::CreateOutputObjects()
{
// Called once per task either in PROOF or local mode. Overload to put some 
// task initialization and/or create your output objects here.
}

//______________________________________________________________________________
TFile *AliAnalysisTask::OpenFile(Int_t iout, Option_t *option) const
{
// This method has to be called INSIDE the user redefined CreateOutputObjects
// method, before creating each object corresponding to the output containers
// that are to be written to a file. This need to be done in general for the big output
// objects that may not fit memory during processing. 
// - 'option' is the file opening option.
//=========================================================================
// NOTE !: The method call will be ignored in PROOF mode, in which case the 
// results have to be streamed back to the client and written just before Terminate()
//=========================================================================
//
// Example:
// void MyAnaTask::CreateOutputObjects() {
//    OpenFile(0);   // Will open the file for the object to be written at output #0
//    fAOD = new TTree("AOD for D0toKPi");
//    OpenFile(1);
// now some histos that should go in the file of the second output container
//    fHist1 = new TH1F("my quality check hist1",...);
//    fHist2 = new TH2F("my quality check hist2",...);
// }
   
   if (iout<0 || iout>=fNoutputs) {
      Error("OpenFile", "No output slot for task %s with index %d", GetName(), iout);
      return NULL;
   }   
   // Method delegated to the analysis manager (A.G. 02/11/09)
   AliAnalysisDataContainer *cont = GetOutputSlot(iout)->GetContainer();
   return AliAnalysisManager::OpenFile(cont, option);
}

//______________________________________________________________________________
Bool_t AliAnalysisTask::Notify()
{
// Overload this IF you need to treat input file change.
   return kTRUE;
}

//______________________________________________________________________________
Bool_t AliAnalysisTask::NotifyBinChange()
{
// Overload this IF you need to treat bin change in event mixing.
   return kTRUE;
}

//______________________________________________________________________________
void AliAnalysisTask::Terminate(Option_t *)
{
// Method called by the framework at the end of data processing.
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
      Error("PostData","Output slot %i not defined for analysis task %s", iout, GetName());
      return kFALSE;
   }
   if (!output->IsConnected()) {
      Error("PostData","Output slot %i of analysis task %s not connected to any data container", iout, GetName());
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
   TString opt(option);
   opt.ToLower();
   Bool_t dep = (opt.Contains("dep"))?kTRUE:kFALSE;
   TString ind;
   Int_t islot;
   AliAnalysisDataContainer *cont;
   for (Int_t i=0; i<indent; i++) ind += " ";
   if (!dep || (dep && IsChecked())) {
      printf("______________________________________________________________________________\n");
      printf("%s\n", Form("%stask: %s  ACTIVE=%i POST_LOOP=%i", ind.Data(), GetName(),IsActive(),IsPostEventLoop()));
      if (dep) const_cast<AliAnalysisTask*>(this)->SetChecked(kFALSE);
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
   if (!fBranchNames.IsNull()) printf("Requested branches:   %s\n", fBranchNames.Data());
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
      if (cont) cont->PrintContainer(option, indent);
   }   
}

//______________________________________________________________________________
void AliAnalysisTask::SetPostEventLoop(Bool_t flag)
{
// Set the task execution mode - run after event loop or not. All output
// containers of this task will get the same type.
   TObject::SetBit(kTaskPostEventLoop,flag);
   AliAnalysisDataContainer *cont;
   Int_t islot;
   for (islot=0; islot<fNoutputs; islot++) {
      cont = GetOutputSlot(islot)->GetContainer();
      if (cont) cont->SetPostEventLoop(flag);
   }   
}
   
//______________________________________________________________________________
void AliAnalysisTask::GetBranches(const char *type, TString &result) const
{
// Get the list of branches for a given type (ESD, AOD). The list of branches
// requested by a task has to ve declared in the form:
//   SetBranches("ESD:branch1,branch2,...,branchN AOD:branch1,branch2,...,branchM")
   result = "";
   if (fBranchNames.IsNull()) return;
   Int_t index1 = fBranchNames.Index(type);
   if (index1<0) return;
   index1 += 1+strlen(type);
   Int_t index2 = fBranchNames.Index(" ", index1);
   if (index2<0) index2 = fBranchNames.Length();
   result = fBranchNames(index1, index2-index1);
}
