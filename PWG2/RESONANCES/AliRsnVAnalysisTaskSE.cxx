//
// Class AliRsnVAnalysisTaskSE
//
// Virtual Class derivated from AliAnalysisTaskSE which will be base class
// for all RSN SE tasks
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include <Riostream.h>

#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliAODEvent.h"
#include "AliRsnEvent.h"
#include "AliRsnTarget.h"

#include "AliRsnVAnalysisTaskSE.h"

ClassImp(AliRsnVAnalysisTaskSE)

//_____________________________________________________________________________
AliRsnVAnalysisTaskSE::AliRsnVAnalysisTaskSE
(const char *name, Bool_t mcOnly) :
   AliAnalysisTaskSE(name),
   fLogType(AliLog::kInfo),
   fLogClassesString(""),
   fESDEvent(0x0),
   fMCEvent(0x0),
   fAODEventIn(0x0),
   fAODEventOut(0x0),
   fMCOnly(mcOnly),
   fRsnEvent(),
   fInfoList(0x0),
   fTaskInfo(name)
{
//
// Default constructor.
// Define the output slot for histograms.
//

   DefineOutput(1, TList::Class());
   DefineOutput(2, TList::Class());
}

//_____________________________________________________________________________
AliRsnVAnalysisTaskSE::AliRsnVAnalysisTaskSE(const AliRsnVAnalysisTaskSE& copy) :
   AliAnalysisTaskSE(copy),
   fLogType(copy.fLogType),
   fLogClassesString(copy.fLogClassesString),
   fESDEvent(copy.fESDEvent),
   fMCEvent(copy.fMCEvent),
   fAODEventIn(copy.fAODEventIn),
   fAODEventOut(copy.fAODEventOut),
   fMCOnly(copy.fMCOnly),
   fRsnEvent(),
   fInfoList(0x0),
   fTaskInfo(copy.fTaskInfo)
{
//
// Copy constructor.
// Defined for coding conventions compliance but never used.
//
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskSE::LocalInit()
{
//
// Local initialization.
// Defines the debug message level and calls the mother class LocalInit().
//

   AliAnalysisTaskSE::LocalInit();
   SetDebugForAllClasses();
}

//_____________________________________________________________________________
Bool_t AliRsnVAnalysisTaskSE::UserNotify()
{
//
// Calls the mother class Notify()
//

   return AliAnalysisTaskSE::UserNotify();
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskSE::ConnectInputData(Option_t *opt)
{
//
// Connect input data, which consist in initializing properly
// the pointer to the input event, which is dynamically casted
// to all available types, and this allows to know its type.
// Calls also the mother class omonyme method.
//

   AliAnalysisTaskSE::ConnectInputData(opt);

   // get AliESDEvent and, if successful
   // retrieve the corresponding MC if exists
   fESDEvent = dynamic_cast<AliESDEvent *>(fInputEvent);
   if (fESDEvent) {
      fMCEvent = (AliMCEvent*) MCEvent();
      AliInfo(Form("Input event is of type ESD   (%p)", fESDEvent));
      if (fMCEvent) AliInfo(Form("Input has an associated MC (%p)", fMCEvent));
   }

   // get AliAODEvent from input and, if successful
   // it will contain both the reconstructed and MC informations
   fAODEventIn = dynamic_cast<AliAODEvent *>(fInputEvent);
   if (fAODEventIn) {
      AliInfo(Form("Input event if of type native AOD (%p)", fAODEventIn));
   }

   // get AliAODEvent from output of previous task
   fAODEventOut = dynamic_cast<AliAODEvent *>(AODEvent());
   if (fAODEventOut) {
      AliInfo(Form("Input event if of type produced AOD from previous step (%p)", fAODEventOut));
   }
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskSE::UserCreateOutputObjects()
{
//
// Creates and links to task all output objects.
// Does explicitly the initialization for the event info class,
// and then calls the customized function which must be overloaded
// in the applications of this base class.
//

   SetDebugForAllClasses();

   // set event info outputs
   fInfoList = new TList();
   fInfoList->SetOwner();
   fTaskInfo.GenerateInfoList(fInfoList);

   // create customized outputs
   RsnUserCreateOutputObjects();

   PostData(1, fInfoList);
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskSE::UserExec(Option_t* opt)
{
//
// Prepares for execution, setting the correct pointers of the
// RSN package event interface, which will point to the not NULL
// objects obtained from dynamic-casts called in ConnectInputData().
//

   if (fMCOnly && fMCEvent) {
      fRsnEvent.SetRef(fMCEvent);
      fRsnEvent.SetRefMC(fMCEvent);
   } else if (fESDEvent) {
      fRsnEvent.SetRef(fESDEvent);
      fRsnEvent.SetRefMC(fMCEvent);
   } else if (fAODEventOut) {
      fRsnEvent.SetRef(fAODEventOut);
      fRsnEvent.SetRefMC(fAODEventOut);
   } else if (fAODEventIn) {
      fRsnEvent.SetRef(fAODEventIn);
      fRsnEvent.SetRefMC(fAODEventIn);
   } else {
      AliError("Unknown input event format. Skipping");
      return;
   }

   // call event preprocessing...
   Bool_t preCheck = EventProcess();
   // ...then fill the information object and print informations...
   fTaskInfo.FillInfo(&fRsnEvent);
   fTaskInfo.PrintInfo(fTaskInfo.GetNumerOfEventsProcessed());
   // ...and return if event did not pass selections
   if (!preCheck) {
      AliDebug(AliLog::kDebug, "Event preprocessing has failed. Skipping event");
      return;
   }


   // call customized implementation for execution
   RsnUserExec(opt);

   // post outputs for the info object
   // (eventually others are done in the derived classes)
   PostData(1, fInfoList);
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskSE::Terminate(Option_t* opt)
{
//
// Termination routines.
// Stores all histograms (after checking they exist)
// and includes to the TList all task informations.
//

   AliAnalysisTask::Terminate();

   TList* list  = dynamic_cast<TList*>(GetOutputData(1));
   if (!list) {
      AliError(Form("At end of analysis, fOutList is %p", list));
      return;
   }

   RsnTerminate(opt);

   TH1I *hEventInfo = (TH1I*) list->FindObject(fTaskInfo.GetEventHistogramName());
   if (!hEventInfo) {
      AliError(Form("hEventInfo is %p", hEventInfo));
      return;
   }
   AliInfo(Form("=== %s ==================", GetName()));
   AliInfo(Form("Number Of Events Processed : %10lld", (Long64_t)hEventInfo->Integral()));
   AliInfo(Form("Number Of Events Accepted  : %10lld", (Long64_t)hEventInfo->GetBinContent(2)));
   AliInfo(Form("Number Of Events Skipped   : %10lld", (Long64_t)hEventInfo->GetBinContent(1)));
   AliInfo(Form("=== end %s ==============", GetName()));

   AliDebug(AliLog::kDebug + 2, "->");
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskSE::RsnUserCreateOutputObjects()
{
//
// Define here all instructions to create output objects.
// This method will be called inside the "UserCreateOutputObjects"
// in the used task.
//

   AliWarning("Implement this in derived classes");
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskSE::RsnUserExec(Option_t*)
{
//
//
//

   AliWarning("Implement this in derived classes");
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskSE::RsnTerminate(Option_t*)
{
//
// Overload this to add additional termination operations
//

   AliWarning("Implement this in derived classes");
}

//_____________________________________________________________________________
Bool_t AliRsnVAnalysisTaskSE::EventProcess()
{
//
// Performs some pre-processing of current event,
// which is useful for all the operations which
// need to be done only once for each event.
//

   // in this case, return always a success
   return kTRUE;
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskSE::SetDebugForAllClasses()
{
//
// Set debug level for all classes for which it is required
//

   TObjArray  *array = fLogClassesString.Tokenize(":");
   TObjString *objStr;
   TString     str;
   Int_t       i, n = array->GetEntriesFast();

   for (i = 0; i < n; i++) {
      objStr = (TObjString*)array->At(i);
      str    = objStr->GetString();
      AliLog::SetClassDebugLevel(str.Data(), fLogType);
      AliInfo(Form("Setting Debug to %s", str.Data()));
   }
}

