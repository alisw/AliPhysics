//
// Class AliRsnVAnalysisTaskSE
//
// Virtual Class derivated from AliAnalysisTaskSE which will be base class
// for all RSN SE tasks
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include "AliRsnVAnalysisTaskSE.h"

#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliAODEvent.h"

ClassImp(AliRsnVAnalysisTaskSE)

//_____________________________________________________________________________
AliRsnVAnalysisTaskSE::AliRsnVAnalysisTaskSE
(const char *name, Int_t numOfOutputs, Bool_t mcOnly) :
    AliAnalysisTaskSE(name),
    fLogType(AliLog::kInfo),
    fLogClassesString(""),
    fESDEvent(0x0),
    fMCEvent(0x0),
    fAODEventIn(0x0),
    fAODEventOut(0x0),
    fMCOnly(mcOnly),
    fRsnEvent(),
    fRsnPIDIndex(),
    fNumberOfOutputs(numOfOutputs),
    fTaskInfo(name)
{
//
// Default constructor.
// Define the output slot for histograms.
//

  AliDebug(AliLog::kDebug+2,"<-");

  if (fNumberOfOutputs<0) fNumberOfOutputs = 0;
  if (fNumberOfOutputs>kMaxNumberOfOutputs) {
    AliWarning(Form("We support only %d outputs. If you need more ask for it.",kMaxNumberOfOutputs));
    AliWarning(Form("For now we are setting it to %d.",kMaxNumberOfOutputs));
    fNumberOfOutputs = kMaxNumberOfOutputs;
  }

  DefineOutput(1, TList::Class());

  AliDebug(AliLog::kDebug+2,"->");
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
    fRsnPIDIndex(),
    fNumberOfOutputs(copy.fNumberOfOutputs),
    fTaskInfo(copy.fTaskInfo)
{
//
// Copy constructor.
// Defined for coding conventions compliance but never used.
//

  AliDebug(AliLog::kDebug+2, "<-");
  AliDebug(AliLog::kDebug+2, "->");
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskSE::LocalInit()
{
//
// Local initialization.
// Defines the debug message level and calls the mother class LocalInit().
//

  SetDebugForAllClasses();

  AliDebug(AliLog::kDebug+2, "<-");
  AliAnalysisTaskSE::LocalInit();
  AliDebug(AliLog::kDebug+2, "->");
}

//_____________________________________________________________________________
Bool_t AliRsnVAnalysisTaskSE::Notify()
{
//
// Calls the mother class Notify()
//

  AliDebug(AliLog::kDebug+2,"<-");
  AliDebug(AliLog::kDebug+2,"->");

  return AliAnalysisTaskSE::Notify();
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskSE::ConnectInputData(Option_t *opt)
{
//
// Connect input data.
// Links the data member pointers to any possible AliVEvenb input
// to the appropriate object belonging to the mother class,
// for a fast retrieval of informations from it through the
// data interface classes provided in this package.
// Makes use of dynamic_cast, in order to know the kind of input
// just checking if the casted pointers are NULL or not.
//

  AliDebug(AliLog::kDebug+2,"<-");
  AliAnalysisTaskSE::ConnectInputData(opt);

  // getting AliESDEvent
  fESDEvent = dynamic_cast<AliESDEvent *>(fInputEvent);

  if (fESDEvent) {
    AliInfo(Form("Input is ESD (%p)", fESDEvent));

    // getting AliMCEvent
    fMCEvent = (AliMCEvent*) MCEvent();
    if (fMCEvent) AliInfo(Form("Input is MC (%p)", fMCEvent));
  }

  // getting AliAODEvent from input
  fAODEventIn = dynamic_cast<AliAODEvent *>(fInputEvent);
  if (fAODEventIn) AliInfo(Form("Input is AOD INPUT (%p)",fAODEventIn));

  // getting AliAODEvent if it is output from previous task
  fAODEventOut = dynamic_cast<AliAODEvent *>(AODEvent());
  if (fAODEventOut) AliInfo(Form("Input is AOD OUTPUT (%p)",fAODEventOut));

  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskSE::RsnUserCreateOutputObjects()
{
//
// Define here all instructions to create output objects.
// This method will be called inside the "UserCreateOutputObjects"
// in the used task.
//

  AliDebug(AliLog::kDebug+2,"<-");
  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskSE::UserCreateOutputObjects()
{
//
// Creates and links to task all output objects.
// They are all stored inside a unique TList which will be saved
// in output slot #1.
//

  SetDebugForAllClasses();

  AliDebug(AliLog::kDebug+2, "<-");

  fOutList[0] = new TList();
  fOutList[0]->SetOwner();
  fTaskInfo.GenerateInfoList(fOutList[0]);
  RsnUserCreateOutputObjects();

  AliDebug(AliLog::kDebug+2,"<-");
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskSE::UserExec(Option_t* opt)
{
//
//
//

  AliDebug(AliLog::kDebug+2,"<-");

  // sets properly the RSN package event interface:
  // if an ESD event is available, it has priority,
  // otherwise the AOD event is used;
  // if the MC information is available, it is linked
  if (fMCOnly && fMCEvent)
    fRsnEvent.SetRef(fMCEvent, fMCEvent);
  else if (fESDEvent)
    fRsnEvent.SetRef(fESDEvent, fMCEvent);
  else if (fAODEventOut)
    fRsnEvent.SetRef(fAODEventOut);
  else if (fAODEventIn)
    fRsnEvent.SetRef(fAODEventIn);
  else {
    AliError("NO ESD or AOD object!!! Skipping ...");
    return;
  }

  // sort tracks w.r. to PID...
  fRsnPIDIndex.FillFromEvent(&fRsnEvent);

  RsnUserExec(opt);

  FillInfo();

  fTaskInfo.PrintInfo(fTaskInfo.GetNumerOfEventsProcessed());

  PostData(1, fOutList[0]);

  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskSE::RsnUserExec(Option_t*)
{
//
//
//

  if (fESDEvent) {
    AliDebug(AliLog::kDebug+1, Form("fESDEvent is %p", fESDEvent));
    AliDebug(AliLog::kDebug,   Form("ESD tracks %d", fESDEvent->GetNumberOfTracks()));
  }
  if (fMCEvent) {
    AliDebug(AliLog::kDebug+1, Form("fMCEvent is %p", fMCEvent));
    AliDebug(AliLog::kDebug,   Form("MC tracks %d", fMCEvent->GetNumberOfTracks()));
  }
  if (fAODEventIn) {
    AliDebug(AliLog::kDebug+1, Form("fAODEventIn is %p", fAODEventIn));
    AliDebug(AliLog::kDebug,   Form("AOD (in) tracks %d", fAODEventIn->GetNumberOfTracks()));
  }

  if (fAODEventOut) {
    AliDebug(AliLog::kDebug+1, Form("fAODEventOut if %p", fAODEventOut));
    AliDebug(AliLog::kDebug,   Form("AOD (out) tracks %d", fAODEventOut->GetNumberOfTracks()));
  }
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskSE::Terminate(Option_t* opt)
{
//
// Termination routines.
// Stores all histograms (after checking they exist)
// and includes to the TList all task informations.
//

  AliDebug(AliLog::kDebug+2,"<-");
  AliAnalysisTask::Terminate();

  TList* list  = dynamic_cast<TList*>(GetOutputData(1));
  if (!list) {
    AliError(Form("At end of analysis, fOutList is %p", list));
    return;
  }

  RsnTerminate(opt);

  TH1I *hEventInfo = (TH1I*) list->FindObject(fTaskInfo.GetEventHistogramName());
  if (!hEventInfo) {
    AliError(Form("hEventInfo is %p",hEventInfo));
    return;
  }

  AliInfo(Form("=== %s ==================",GetName()));
  AliInfo(Form("Number Of Events Processed : %10l",(Long64_t)hEventInfo->Integral()));
  AliInfo(Form("Number Of Events Accepted  : %10l",(Long64_t)hEventInfo->GetBinContent(2)));
  AliInfo(Form("Number Of Events Skipped   : %10l",(Long64_t)hEventInfo->GetBinContent(1)));
  AliInfo(Form("=== end %s ==============",GetName()));

  AliDebug(AliLog::kDebug+2, "->");
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskSE::RsnTerminate(Option_t*)
{
//
// Overload this to add additional termination operations
//

  AliDebug(AliLog::kDebug+2, "<-");
  AliDebug(AliLog::kDebug+2, "->");
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskSE::FillInfo()
{
//
// Fill information object with statistics of analysis
//

  AliDebug(AliLog::kDebug+2, "<-");

  fTaskInfo.FillInfo();

  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskSE::SetLogType(AliLog::EType_t type, TString allClasses)
{
//
// Set Log level for this and other classes (list of their names)
//

  AliDebug(AliLog::kDebug+2,"<-");
  fLogType = type;
  fLogClassesString = allClasses;
  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskSE::SetDebugForAllClasses()
{
//
// Set debug level for all classes for which it is required
//

  AliDebug(AliLog::kDebug+2, "<-");
  TObjArray* array = fLogClassesString.Tokenize(":");
  TObjString *str;
  TString strr;
  for (Int_t i=0;i< array->GetEntriesFast();i++) {
    str = (TObjString *) array->At(i);
    strr = str->GetString();
    AliLog::SetClassDebugLevel(strr.Data(), fLogType);
    AliInfo(Form("Setting Debug to %s",strr.Data()));
  }
  AliDebug(AliLog::kDebug+2,"->");
}
