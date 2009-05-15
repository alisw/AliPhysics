//
// Class AliRsnVAnalysisTaskME
//
// Virtual Class derivated from AliAnalysisTaskME which will be base class
// for all RSN ME tasks
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

// #include "AliAnalysisManager.h"
// #include "AliAODHandler.h"

#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliAODEvent.h"

#include "AliRsnVAnalysisTaskME.h"

ClassImp(AliRsnVAnalysisTaskME)

//_____________________________________________________________________________
AliRsnVAnalysisTaskME::AliRsnVAnalysisTaskME(const char *name) :
    AliAnalysisTaskME(name),
    fLogType(AliLog::kInfo),
    fLogClassesString(""),
    fESDEvent(0x0),
    fMCEvent(0x0),
    fAODEvent(0x0),
    fOutList(0x0),
    fTaskInfo(name)
{
//
// Default constructor
//
  AliDebug(AliLog::kDebug+2,"<-");
  DefineOutput(1, TList::Class());
  AliDebug(AliLog::kDebug+2,"->");
}

AliRsnVAnalysisTaskME::AliRsnVAnalysisTaskME(const AliRsnVAnalysisTaskME& copy) : AliAnalysisTaskME(copy),
    fLogType(copy.fLogType),
    fLogClassesString(copy.fLogClassesString),
    fESDEvent(copy.fESDEvent),
    fMCEvent(copy.fMCEvent),
    fAODEvent(copy.fAODEvent),
    fOutList(copy.fOutList),
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
  AliLog::SetClassDebugLevel(GetName(), fLogType);
  SetDebugForOtherClasses();
  
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
  return AliAnalysisTaskME::Notify();
  AliDebug(AliLog::kDebug+2,"->");
}


//_____________________________________________________________________________
void AliRsnVAnalysisTaskME::ConnectInputData(Option_t *opt)
{
  AliDebug(AliLog::kDebug+2,"<-");
  AliAnalysisTaskME::ConnectInputData(opt);

//     fESDEvent = dynamic_cast<AliESDEvent *> (fInputEvent);
//     if (fESDEvent)
//         AliInfo(Form("Input is ESD (%p)",fESDEvent));
//
//     fMCEvent = (AliMCEvent*) MCEvent();
//     if (fMCEvent)
//         AliInfo(Form("Input is MC (%p)",fMCEvent));
//
//     fAODEvent = dynamic_cast<AliAODEvent *> (fInputEvent);
//     if (fAODEvent)
//         AliInfo(Form("Input is AOD INPUT (%p)",fAODEvent));
//
//     if (!fAODEvent) {
//         fAODEvent = dynamic_cast<AliAODEvent *> (AODEvent());
//         if (fAODEvent)
//             AliInfo(Form("Input is AOD OUTPUT (%p)",fAODEvent));
//     }



  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskME::RsnUserCreateOutputObjects()
{
  AliDebug(AliLog::kDebug+2,"<-");

  AliDebug(AliLog::kDebug+2,"->");
}


//_____________________________________________________________________________
void AliRsnVAnalysisTaskME::UserCreateOutputObjects()
{
  AliLog::SetClassDebugLevel(GetName(), fLogType);
  
  SetDebugForOtherClasses();
  
  AliDebug(AliLog::kDebug+2,"<-");

  fOutList = new TList();
  fOutList->SetOwner();

  fOutList->Add(fTaskInfo.GenerateInfoList());

  RsnUserCreateOutputObjects();

  AliDebug(AliLog::kDebug+2,"<-");
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskME::UserExec(Option_t* opt)
{

  AliDebug(AliLog::kDebug+2,"<-");

  RsnUserExec(opt);

  FillInfo();

  fTaskInfo.PrintInfo(fEntry);

  PostData(1, fOutList);

  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskME::RsnUserExec(Option_t* )
{
  AliDebug(AliLog::kDebug+2,"<-");
//     if (fESDEvent) {
//         AliDebug(AliLog::kDebug+1,Form("fESDEvent if %p",fESDEvent));
//         AliDebug(AliLog::kDebug,Form("ESD tracks %d",fESDEvent->GetNumberOfTracks()));
//     }
//     if (fMCEvent) {
//         AliDebug(AliLog::kDebug+1,Form("fMCEvent if %p",fMCEvent));
//         AliDebug(AliLog::kDebug,Form("MC tracks %d",fMCEvent->GetNumberOfTracks()));
//     }
//     if (fAODEvent) {
//         AliDebug(AliLog::kDebug+1,Form("fAODEvent if %p",fAODEvent));
//         AliDebug(AliLog::kDebug,Form("AOD tracks %d",fAODEvent->GetNumberOfTracks()));
//
//
//     }
  AliAODEvent* aod1 = (AliAODEvent*)GetEvent(0);
  AliAODEvent* aod2 = (AliAODEvent*)GetEvent(1);
  AliDebug(AliLog::kDebug,Form("AOD tracks %d",aod1->GetNumberOfTracks()));
  AliDebug(AliLog::kDebug,Form("AOD tracks %d",aod2->GetNumberOfTracks()));
  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskME::Terminate(Option_t* opt)
{
  AliDebug(AliLog::kDebug+2,"<-");
  AliAnalysisTask::Terminate();

  fOutList = dynamic_cast<TList*>(GetOutputData(1));
  if (!fOutList) {
    AliError(Form("At end of analysis, fOutList is %p",fOutList));
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

  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskME::RsnTerminate(Option_t* )
{
  AliDebug(AliLog::kDebug+2,"<-");
  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskME::FillInfo()
{

  if (fESDEvent) {
    fTaskInfo.SetNumberOfTracks(fESDEvent->GetNumberOfTracks());
  } else if (fAODEvent) {
    fTaskInfo.SetNumberOfTracks(fAODEvent->GetNumberOfTracks());
  }

  fTaskInfo.FillInfo();
}

//_____________________________________________________________________________
void AliRsnVAnalysisTaskME::SetLogType(AliLog::EType_t type,TString otherClasses)
{
  AliDebug(AliLog::kDebug+2,"<-");
  fLogType = type;
  fLogClassesString = otherClasses;
  AliDebug(AliLog::kDebug+2,"->");
}
//_____________________________________________________________________________
void AliRsnVAnalysisTaskME::SetDebugForOtherClasses()
{
  AliDebug(AliLog::kDebug+2,"<-");
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
