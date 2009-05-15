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

ClassImp(AliRsnVAnalysisTaskSE)

//_____________________________________________________________________________
AliRsnVAnalysisTaskSE::AliRsnVAnalysisTaskSE(const char *name) :
  AliAnalysisTaskSE(name),
  fLogType(AliLog::kInfo),
  fLogClassesString(""),
  fESDEvent(0x0),
  fMCEvent(0x0),
  fAODEventIn(0x0),
  fAODEventOut(0x0),
  fOutList(0x0),
  fTaskInfo(name)
{
//
// Default constructor.
// Define the output slot for histograms.
//

  AliDebug(AliLog::kDebug+2,"<-");

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
  fOutList(copy.fOutList),
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

  AliLog::SetClassDebugLevel(GetName(), fLogType);
  SetDebugForOtherClasses();

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
  fESDEvent = dynamic_cast<AliESDEvent *> (fInputEvent);

  if (fESDEvent) {
    AliInfo(Form("Input is ESD (%p)", fESDEvent));

    // getting AliMCEvent
    fMCEvent = (AliMCEvent*) MCEvent();
    if (fMCEvent) AliInfo(Form("Input is MC (%p)", fMCEvent));
  }

  // getting AliAODEvent from input
  fAODEventIn = dynamic_cast<AliAODEvent *> (fInputEvent);
  if (fAODEventIn) AliInfo(Form("Input is AOD INPUT (%p)",fAODEventIn));

  // getting AliAODEvent if it is output from previous task
  fAODEventOut = dynamic_cast<AliAODEvent *> (AODEvent());
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

  AliLog::SetClassDebugLevel(GetName(), fLogType);

  SetDebugForOtherClasses();

  AliDebug(AliLog::kDebug+2, "<-");

  fOutList = new TList();
  fOutList->SetOwner();

  fOutList->Add(fTaskInfo.GenerateInfoList());

  RsnUserCreateOutputObjects();

  AliDebug(AliLog::kDebug+2,"<-");
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskSE::UserExec(Option_t* opt)
{

  AliDebug(AliLog::kDebug+2,"<-");

  RsnUserExec(opt);

  FillInfo();

  fTaskInfo.PrintInfo(fEntry);

  PostData(1, fOutList);

  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskSE::RsnUserExec(Option_t* )
{

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

  fOutList = dynamic_cast<TList*>(GetOutputData(1));
  if (!fOutList) {
    AliError(Form("At end of analysis, fOutList is %p", fOutList));
    return;
  }

  RsnTerminate(opt);

  TList* lEventInfo = (TList*) fOutList->FindObject(fTaskInfo.GetName());

  TH1I *hEventInfo = (TH1I*) lEventInfo->FindObject(fTaskInfo.GetEventHistogramName());
  if (!hEventInfo) {
    AliError(Form("hEventInfo is %p",hEventInfo));
    return;
  }

  AliInfo(Form("=== %s ==================",GetName()));
  AliInfo(Form("Number Of Events Processed : %10d",(Long64_t)hEventInfo->Integral()));
  AliInfo(Form("Number Of Events Accepted  : %10d",(Long64_t)hEventInfo->GetBinContent(2)));
  AliInfo(Form("Number Of Events Skipped   : %10d",(Long64_t)hEventInfo->GetBinContent(1)));
  AliInfo(Form("=== end %s ==============",GetName()));

  AliDebug(AliLog::kDebug+2, "->");
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskSE::RsnTerminate(Option_t* )
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
  if (fAODEventOut) {
    fTaskInfo.SetNumberOfTracks(fAODEventOut->GetNumberOfTracks());
  }
  else if (fESDEvent) {
    fTaskInfo.SetNumberOfTracks(fESDEvent->GetNumberOfTracks());
  }
  else if (fAODEventIn) {
    fTaskInfo.SetNumberOfTracks(fAODEventIn->GetNumberOfTracks());
  }

  fTaskInfo.FillInfo();
  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskSE::SetLogType(AliLog::EType_t type, TString otherClasses)
{
//
// Set Log level for this and other classes (list of their names)
//

  AliDebug(AliLog::kDebug+2,"<-");
  fLogType = type;
  fLogClassesString = otherClasses;
  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskSE::SetDebugForOtherClasses()
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
  }
  AliDebug(AliLog::kDebug+2,"->");
}
