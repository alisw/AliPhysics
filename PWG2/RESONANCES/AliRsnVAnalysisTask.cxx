//
// Class AliRsnVAnalysisTask
//
// Virtual Class derivated from AliAnalysisTaskSE which will be base class
// for all RSN SE tasks
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include <TH1.h>

#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMultiInputEventHandler.h"
#include "AliMixInputEventHandler.h"


#include "AliRsnEvent.h"
#include "AliRsnTarget.h"

#include "AliRsnVAnalysisTask.h"

ClassImp(AliRsnVAnalysisTask)

//_____________________________________________________________________________
AliRsnVAnalysisTask::AliRsnVAnalysisTask
(const char *name, Bool_t mcOnly) :
   AliAnalysisTaskSE(name),
   fLogType(AliLog::kInfo),
   fLogClassesString(""),
   fIsMixing(kFALSE),
   fMCOnly(mcOnly),
   fInfoList(0x0),
   fTaskInfo(name),
   fMixedEH(0),
   fUseMixingRange(kTRUE)
{
//
// Default constructor.
// Define the output slot for histograms.
//

   DefineOutput(1, TList::Class());
   DefineOutput(2, TList::Class());
   
   for (Int_t i = 0; i < 2; i++) {
      fESDEvent[i] = 0;
      fMCEvent[i] = 0;
      fAODEventIn[i] = 0;
      fAODEventOut[i] = 0;
   }
}

//_____________________________________________________________________________
AliRsnVAnalysisTask::AliRsnVAnalysisTask(const AliRsnVAnalysisTask& copy) :
   AliAnalysisTaskSE(copy),
   fLogType(copy.fLogType),
   fLogClassesString(copy.fLogClassesString),
   fIsMixing(copy.fIsMixing),
   fMCOnly(copy.fMCOnly),
   fInfoList(0x0),
   fTaskInfo(copy.fTaskInfo),
   fMixedEH(copy.fMixedEH),
   fUseMixingRange(copy.fUseMixingRange)
{
//
// Copy constructor.
// Defined for coding conventions compliance but never used.
//
   AliDebug(AliLog::kDebug + 2, "<-");
   
   for (Int_t i = 0; i < 2; i++) {
      fESDEvent[i] = copy.fESDEvent[i];
      fMCEvent[i] = copy.fMCEvent[i];
      fAODEventIn[i] = copy.fAODEventIn[i];
      fAODEventOut[i] = copy.fAODEventOut[i];
      fRsnEvent[i] = copy.fRsnEvent[i];
   }
   
   AliDebug(AliLog::kDebug + 2, "->");
}

//_____________________________________________________________________________
void AliRsnVAnalysisTask::LocalInit()
{
//
// Local initialization.
// Defines the debug message level and calls the mother class LocalInit().
//

   AliAnalysisTaskSE::LocalInit();
   SetDebugForAllClasses();
}

//_____________________________________________________________________________
Bool_t AliRsnVAnalysisTask::UserNotify()
{
//
// Calls the mother class Notify()
//

   return AliAnalysisTaskSE::UserNotify();
}

//_____________________________________________________________________________
void AliRsnVAnalysisTask::ConnectInputData(Option_t *opt)
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
   fESDEvent[0] = dynamic_cast<AliESDEvent *>(fInputEvent);
   if (fESDEvent[0]) {
      fMCEvent[0] = (AliMCEvent*) MCEvent();
      AliInfo(Form("Input event is of type ESD   (%p)", fESDEvent[0]));
      if (fMCEvent[0]) AliInfo(Form("Input has an associated MC (%p)", fMCEvent[0]));
   }

   // get AliAODEvent from input and, if successful
   // it will contain both the reconstructed and MC informations
   fAODEventIn[0] = dynamic_cast<AliAODEvent *>(fInputEvent);
   if (fAODEventIn[0]) {
      AliInfo(Form("Input event if of type native AOD (%p)", fAODEventIn[0]));
   }

   // get AliAODEvent from output of previous task
   fAODEventOut[0] = dynamic_cast<AliAODEvent *>(AODEvent());
   if (fAODEventOut[0]) {
      AliInfo(Form("Input event if of type produced AOD from previous step (%p)", fAODEventOut[0]));
   }
}

//_____________________________________________________________________________
void AliRsnVAnalysisTask::UserCreateOutputObjects()
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
void AliRsnVAnalysisTask::UserExec(Option_t* opt)
{
//
// Prepares for execution, setting the correct pointers of the
// RSN package event interface, which will point to the not NULL
// objects obtained from dynamic-casts called in ConnectInputData().
//

   if (!IsMixing()) {

      if (fMCOnly && fMCEvent) {
         fRsnEvent[0].SetRef(fMCEvent[0]);
         fRsnEvent[0].SetRefMC(fMCEvent[0]);
      } else if (fESDEvent[0]) {
         fRsnEvent[0].SetRef(fESDEvent[0]);
         fRsnEvent[0].SetRefMC(fMCEvent[0]);
      } else if (fAODEventOut[0]) {
         fRsnEvent[0].SetRef(fAODEventOut[0]);
         fRsnEvent[0].SetRefMC(fAODEventOut[0]);
      } else if (fAODEventIn[0]) {
         fRsnEvent[0].SetRef(fAODEventIn[0]);
         fRsnEvent[0].SetRefMC(fAODEventIn[0]);
      } else {
         AliError("Unknown input event format. Skipping");
         return;
      }

      // call event preprocessing...
      Bool_t preCheck = RsnEventProcess();
      // ...then fill the information object and print informations...
      fTaskInfo.FillInfo(&fRsnEvent[0]);
      fTaskInfo.PrintInfo(fTaskInfo.GetNumerOfEventsProcessed());
      // ...and return if event did not pass selections
      if (!preCheck) {
         AliDebug(AliLog::kDebug, "Event preprocessing has failed. Skipping event");
         return;
      }
      
      // call customized implementation for execution
      RsnUserExec(opt);
   }
   
   // post outputs for the info object
   // (eventually others are done in the derived classes)
   PostData(1, fInfoList);
}

void AliRsnVAnalysisTask::UserExecMix(Option_t* option)
{
   AliDebug(AliLog::kDebug + 2, "<-");


   if (!IsMixing()) return;

   SetupMixingEvents();

   if (!fMixedEH) return;

   if (fMCOnly && fMCEvent[0]) {
      fRsnEvent[0].SetRef(fMCEvent[0]);
      fRsnEvent[0].SetRefMC(fMCEvent[0]);
      fRsnEvent[1].SetRef(fMCEvent[1]);
      fRsnEvent[1].SetRefMC(fMCEvent[1]);
   } else if (fESDEvent[0]) {
      fRsnEvent[0].SetRef(fESDEvent[0]);
      fRsnEvent[0].SetRefMC(fMCEvent[0]);
      fRsnEvent[1].SetRef(fESDEvent[1]);
      fRsnEvent[1].SetRefMC(fMCEvent[1]);
   } else if (fAODEventOut) {
      fRsnEvent[0].SetRef(fAODEventOut[0]);
      fRsnEvent[0].SetRefMC(fAODEventOut[0]);
      fRsnEvent[1].SetRef(fAODEventOut[1]);
      fRsnEvent[1].SetRefMC(fAODEventOut[1]);
   } else if (fAODEventIn) {
      fRsnEvent[0].SetRef(fAODEventIn[0]);
      fRsnEvent[0].SetRefMC(fAODEventIn[0]);
      fRsnEvent[1].SetRef(fAODEventIn[1]);
      fRsnEvent[1].SetRefMC(fAODEventIn[1]);
   } else {
      AliError("NO ESD or AOD object!!! Skipping ...");
      return;
   }

   // call event preprocessing...
   Bool_t preCheck = RsnEventProcess();
   // ...then fill the information object and print informations...
   fTaskInfo.FillInfo(&fRsnEvent[0]);
   fTaskInfo.PrintInfo(fTaskInfo.GetNumerOfEventsProcessed());
   // ...and return if event did not pass selections
   if (!preCheck) {
      AliDebug(AliLog::kDebug, "Event preprocessing has failed. Skipping event");
      return;
   }

   RsnUserExecMix(option);
   AliDebug(AliLog::kDebug + 2, "->");
}


//_____________________________________________________________________________
void AliRsnVAnalysisTask::Terminate(Option_t* opt)
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
void AliRsnVAnalysisTask::RsnUserCreateOutputObjects()
{
//
// Define here all instructions to create output objects.
// This method will be called inside the "UserCreateOutputObjects"
// in the used task.
//

   AliWarning("Implement this in derived classes");
}

//_____________________________________________________________________________
void AliRsnVAnalysisTask::RsnUserExec(Option_t*)
{
//
//
//

   AliWarning("Implement this in derived classes");
}

void AliRsnVAnalysisTask::RsnUserExecMix(Option_t*)
{
   //
   //
   //

   AliWarning("Implement this in derived classes");
}

//_____________________________________________________________________________
void AliRsnVAnalysisTask::RsnTerminate(Option_t*)
{
//
// Overload this to add additional termination operations
//

   AliWarning("Implement this in derived classes");
}

//_____________________________________________________________________________
Bool_t AliRsnVAnalysisTask::RsnEventProcess()
{
//
// Performs some pre-processing of current event,
// which is useful for all the operations which
// need to be done only once for each event.
//

   // if not using mixing cuts return kTRUE
   if (!IsUsingMixingRange()) return kTRUE;

   // cut if event was in range of mixing cuts
   AliVEventHandler *inh = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
   if (inh->InheritsFrom(AliMultiInputEventHandler::Class())) {
      AliMultiInputEventHandler *inEvHMain = static_cast<AliMultiInputEventHandler *>(inh);
      if (inEvHMain->GetFirstMultiInputHandler()->InheritsFrom(AliMixInputEventHandler::Class())) {
         fMixedEH = static_cast<AliMixInputEventHandler *>(inEvHMain->GetFirstMultiInputHandler());
         if (fMixedEH) {
            if (fMixedEH->CurrentBinIndex() < 0) return kFALSE;
         }
      }
   }

   // in this case, return always a success
   return kTRUE;
}

//_____________________________________________________________________________
void AliRsnVAnalysisTask::SetupMixingEvents()
{
//
// Setup the pointers to the mixed event.
// This requires to retrieve them from the mixed event handler
//
   
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   AliMultiInputEventHandler *inEvHMain = dynamic_cast<AliMultiInputEventHandler *>(mgr->GetInputEventHandler());
   if (inEvHMain) {
      fMixedEH = dynamic_cast<AliMixInputEventHandler *>(inEvHMain->GetFirstMultiInputHandler());
      if (fMixedEH) {
         AliMultiInputEventHandler *inEvHMainTmpMix = dynamic_cast<AliMultiInputEventHandler *>(fMixedEH->InputEventHandler(0));
         if (!inEvHMainTmpMix) return;
         AliInputEventHandler *inEvHMixTmp = dynamic_cast<AliInputEventHandler *>(inEvHMainTmpMix->GetFirstInputEventHandler());
         AliMCEventHandler *inEvHMCMixTmp = dynamic_cast<AliMCEventHandler *>(inEvHMainTmpMix->GetFirstMCEventHandler());
         if (!inEvHMixTmp) return;
         fESDEvent[1] = dynamic_cast<AliESDEvent *>(inEvHMixTmp->GetEvent());
         if (fESDEvent[1]) AliDebug(AliLog::kDebug, Form("Input is ESD (%p) MIXED", fESDEvent[1]));
         // getting AliAODEvent from input
         fAODEventIn[1] = dynamic_cast<AliAODEvent *>(inEvHMixTmp->GetEvent());
         if (fAODEventIn[1]) AliDebug(AliLog::kDebug, Form("Input is AOD (%p) MIXED", fAODEventIn[1]));
         // getting AliAODEvent if it is output from previous task (not supported)
         fAODEventOut[1] = 0;

         if (inEvHMCMixTmp) {
            fMCEvent[1] = inEvHMCMixTmp->MCEvent();
            if (fMCEvent[1]) AliDebug(AliLog::kDebug, Form("Input is ESDMC (%p) MIXED", fMCEvent[1]));
         }

      }

   }
}

//_____________________________________________________________________________
void AliRsnVAnalysisTask::SetDebugForAllClasses()
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

