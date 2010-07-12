//
// Class AliRsnVAnalysisTaskME
//
// Virtual Class derivated from AliAnalysisTaskME which will be base class
// for all RSN ME tasks
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include <TH1.h>
#include <AliLog.h>
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliMultiEventInputHandler.h"
#include "AliRsnVAnalysisTaskME.h"

ClassImp(AliRsnVAnalysisTaskME)

//_____________________________________________________________________________
AliRsnVAnalysisTaskME::AliRsnVAnalysisTaskME(const char *name, Int_t numOfOutputs) :
    AliAnalysisTaskME(name),
    fLogType(AliLog::kInfo),
    fLogClassesString(""),
    fESDEvent(0x0),
    fMCEvent(0x0),
    fAODEvent(0x0),
    fNumberOfOutputs(numOfOutputs),
    fTaskInfo(name)
{
//
// Default constructor
//
  AliDebug(AliLog::kDebug+2,"<-");

  if (fNumberOfOutputs < 0) fNumberOfOutputs = 0;
  if (fNumberOfOutputs > kMaxNumberOfOutputs) {
    AliWarning(Form("We support only %d outputs. If you need more ask for it.", kMaxNumberOfOutputs));
    AliWarning(Form("For now we are setting it to %d.", kMaxNumberOfOutputs));
    fNumberOfOutputs = kMaxNumberOfOutputs;
  }

  DefineOutput(1, TList::Class());

  AliDebug(AliLog::kDebug+2,"->");
}

AliRsnVAnalysisTaskME::AliRsnVAnalysisTaskME(const AliRsnVAnalysisTaskME& copy) : AliAnalysisTaskME(copy),
    fLogType(copy.fLogType),
    fLogClassesString(copy.fLogClassesString),
    fESDEvent(copy.fESDEvent),
    fMCEvent(copy.fMCEvent),
    fAODEvent(copy.fAODEvent),
    fNumberOfOutputs(copy.fNumberOfOutputs),
    fTaskInfo(copy.fTaskInfo)
{
  AliDebug(AliLog::kDebug+2,"<-");
  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskME::LocalInit()
{
//
// LocalInit()
//

  SetDebugForAllClasses();
  AliDebug(AliLog::kDebug+2,"<-");
  AliAnalysisTaskME::LocalInit();
  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
Bool_t AliRsnVAnalysisTaskME::Notify()
{
//
// Notify()
//

  AliDebug(AliLog::kDebug+2,"<-");
  if (!AliAnalysisTaskME::Notify()) return kFALSE;
  AliDebug(AliLog::kDebug+2,"->");

  return kTRUE;
}


//_____________________________________________________________________________
void AliRsnVAnalysisTaskME::ConnectInputData(Option_t *opt)
{
//
// ConnectInputData
//

  SetDebugForAllClasses();

  AliDebug(AliLog::kDebug+2,"<-");

  AliAnalysisTaskME::ConnectInputData(opt);

  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskME::RsnUserCreateOutputObjects()
{
//
// Rsn User Create Output Objects
//

  AliDebug(AliLog::kDebug+2,"<-");

  AliDebug(AliLog::kDebug+2,"->");
}


//_____________________________________________________________________________
void AliRsnVAnalysisTaskME::UserCreateOutputObjects()
{
//
// User Create Output Objects
//

  SetDebugForAllClasses();

  AliDebug(AliLog::kDebug+2,"<-");

  fOutList[0] = new TList();
  fOutList[0]->SetOwner();
  fTaskInfo.GenerateInfoList(fOutList[0]);
  RsnUserCreateOutputObjects();

  AliDebug(AliLog::kDebug+2,"<-");
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskME::UserExec(Option_t* opt)
{
//
// User Exec
//

  AliDebug(AliLog::kDebug+2,"<-");

  RsnUserExec(opt);

  FillInfo();

  fTaskInfo.PrintInfo(fEntry);

  PostData(1, fOutList[0]);

  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskME::RsnUserExec(Option_t*)
{
//
// Rsn User Exec
//

  AliDebug(AliLog::kDebug+2,"<-");

  if (!CheckAndPrintEvents()) return;

  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
Bool_t AliRsnVAnalysisTaskME::CheckAndPrintEvents()
{
//
// Check for supported events
// return false in the case of unkown format
// or number of events is less or equal 1
//

  AliInfo(Form("Current Entry %l", Entry()));
  Int_t nEvents = fInputHandler->GetBufferSize();
  if (nEvents <= 1) return kFALSE;
  fESDEvent = dynamic_cast<AliESDEvent*>(GetEvent(0));
  fAODEvent = dynamic_cast<AliAODEvent*>(GetEvent(0));

  if (fESDEvent) {
    AliESDEvent **allESDEvents = new AliESDEvent*[nEvents];
    for (Int_t i = 0; i < nEvents; i++) {
      allESDEvents[i] = dynamic_cast<AliESDEvent*>(GetEvent(i));
      if (!allESDEvents[i]) {
        AliWarning(Form("Null ESD event in index %d", i));
        continue;
      }
      AliDebug(AliLog::kDebug, Form("ESD event %d has %d tracks", i, allESDEvents[i]->GetNumberOfTracks()));
    }
  } else if (fAODEvent) {
    AliAODEvent **allAODEvents = new AliAODEvent*[nEvents];
    for (Int_t i = 0; i < nEvents; i++) {
      allAODEvents[i] = dynamic_cast<AliAODEvent*>(GetEvent(i));
      if (!allAODEvents[i]) {
        AliWarning(Form("Null AOD event in index %d", i));
        continue;
      }
      AliDebug(AliLog::kDebug, Form("AOD event %d has %d tracks", i, allAODEvents[i]->GetNumberOfTracks()));
    }
  } else {
    AliWarning("Unknown input format");
    return kFALSE;
  }

  return kTRUE;
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskME::Terminate(Option_t* opt)
{
//
// Terminate
//

  AliDebug(AliLog::kDebug+2,"<-");
  AliAnalysisTask::Terminate();

  fOutList[0] = dynamic_cast<TList*>(GetOutputData(1));
  if (!fOutList[0]) {
    AliError(Form("At end of analysis, fOutList is %p", fOutList));
    return;
  }

  RsnTerminate(opt);

  TH1I *hEventInfo = (TH1I*) fOutList[0]->FindObject(fTaskInfo.GetEventHistogramName());
  if (!hEventInfo) {
    AliError(Form("hEventInfo is %p", hEventInfo));
    return;
  }

  AliInfo(Form("=== %s ==================", GetName()));
  AliInfo(Form("Number Of Events Processed : %10l", (Long64_t)hEventInfo->Integral()));
  AliInfo(Form("Number Of Events Accepted  : %10l", (Long64_t)hEventInfo->GetBinContent(2)));
  AliInfo(Form("Number Of Events Skipped   : %10l", (Long64_t)hEventInfo->GetBinContent(1)));
  AliInfo(Form("=== end %s ==============", GetName()));

  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskME::RsnTerminate(Option_t*)
{
//
// RsnTerminate
//

  AliDebug(AliLog::kDebug+2,"<-");
  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskME::FillInfo()
{
//
// Fills Info
//

  fTaskInfo.FillInfo();
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskME::SetLogType(AliLog::EType_t type, TString allClasses)
{
//
// Sets Log Type
//

  AliDebug(AliLog::kDebug+2,"<-");
  fLogType = type;
  fLogClassesString = allClasses;
  AliDebug(AliLog::kDebug+2,"->");
}
//_____________________________________________________________________________
void AliRsnVAnalysisTaskME::SetDebugForAllClasses()
{
//
// Set Debug For All Classes
//

  AliDebug(AliLog::kDebug+2,"<-");
  TObjArray* array = fLogClassesString.Tokenize(":");
  TObjString *str;
  TString strr;
  for (Int_t i = 0;i < array->GetEntriesFast();i++) {
    str = (TObjString *) array->At(i);
    strr = str->GetString();
    AliLog::SetClassDebugLevel(strr.Data(), fLogType);
    AliInfo(Form("Setting Debug to %s", strr.Data()));
  }
  AliDebug(AliLog::kDebug+2,"->");
}
