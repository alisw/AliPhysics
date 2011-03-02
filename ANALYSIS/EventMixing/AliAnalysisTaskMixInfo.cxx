//
// Class AliAnalysisTaskMixInfo
//
// AliAnalysisTaskMixInfo is task
// for mixing info
//
// Mixing info can be enabled by setting one of following lines in UserCreateOutput() in your task
//       // prints mixing info
//       AliLog::SetClassDebugLevel("AliAnalysisTaskMixInfo", AliLog::kDebug);
//       // prints mixing info + event info for both (main and mixed) events
//       AliLog::SetClassDebugLevel("AliAnalysisTaskMixInfo", AliLog::kDebug+1);
//
// authors:
//          Martin Vala (martin.vala@cern.ch)
//

#include <TList.h>
#include <TObjString.h>

#include "AliAnalysisManager.h"

#include "AliMixInputEventHandler.h"
#include "AliAnalysisTaskMixInfo.h"
#include "AliMixInfo.h"
#include "AliMixEventPool.h"
#include "AliMixEventCutObj.h"


ClassImp(AliAnalysisTaskMixInfo)

//________________________________________________________________________
AliAnalysisTaskMixInfo::AliAnalysisTaskMixInfo(const char *name)
   : AliAnalysisTaskSE(name),
     fInputEHMain(0),
     fInputEHMix(0),
     fOutputList(0),
     fMixInfo(0),
     fCurrentEntryTmp(-1),
     fLogType(AliLog::kInfo),
     fLogClassesString()
{
   //
   // Constructor
   //
   DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskMixInfo::~AliAnalysisTaskMixInfo()
{
   //
   // Destructor
   //
   if (fOutputList && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) delete fOutputList;
}

//________________________________________________________________________
void AliAnalysisTaskMixInfo::UserCreateOutputObjects()
{
   // Create histograms
   // Called once

   SetDebugForAllClasses();
   fOutputList = new TList();
   fOutputList->SetOwner(kTRUE);

   // sets all Inuput Handler pointers
   InitInputHandlers();

   // inits mix info
   InitMixInfo();

   if (fInputEHMix) {
      AliMixEventPool *evPool = fInputEHMix->GetEventPool();
      if (evPool) {
         evPool->SetBufferSize(fInputEHMix->BufferSize());
         evPool->SetMixNumber(fInputEHMix->MixNumber());
         fMixInfo->SetEventPool(evPool);
      }
   }
   if (fMixInfo) fOutputList->Add(fMixInfo);

   // Post output data.
   PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskMixInfo::UserExec(Option_t *)
{
   // Main loop
   // Called for each event
   if (fMixInfo && fInputEHMix) {
      if (fInputEHMix->BufferSize() > 1) {
         if (fInputEHMix->NumberMixedTimes() >= fInputEHMix->BufferSize())
            fMixInfo->FillHistogram(AliMixInfo::kMainEvents, fInputEHMix->CurrentBinIndex());
      } else {
         if ((!fInputEHMix->IsMixingIfNotEnoughEvents())) {
            if (fInputEHMix->NumberMixed() == fInputEHMix->MixNumber())
               // add main entry only when there was enough mixed events mixed
               fMixInfo->FillHistogram(AliMixInfo::kMainEvents, fInputEHMix->CurrentBinIndex());
         } else {
            fMixInfo->FillHistogram(AliMixInfo::kMainEvents, fInputEHMix->CurrentBinIndex());
         }
      }
      AliDebug(AliLog::kDebug, Form("Main %lld %d [%lld,%lld] %d", fInputEHMix->CurrentEntry(), fInputEHMix->CurrentBinIndex(), fInputEHMix->CurrentEntryMain(), fInputEHMix->CurrentEntryMix(), fInputEHMix->NumberMixed()));
   }
   // Post output data.
   PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskMixInfo::UserExecMix(Option_t *)
{
   // UserExecMix

   if (!fInputEHMix) return;

   // fills bin index (even when it is -1, so we know out of range combinations)
   if (fMixInfo) fMixInfo->FillHistogram(AliMixInfo::kMixedEvents, fInputEHMix->CurrentBinIndex());

   // just test
   if (fInputEHMix->CurrentEntryMix() < 0) {
      AliError("Mix entry is -1 and it should not happen !!!!!");
      return ;
   }
   AliDebug(AliLog::kDebug, Form("Mixing %lld %d [%lld,%lld] %d", fInputEHMix->CurrentEntry(), fInputEHMix->CurrentBinIndex(), fInputEHMix->CurrentEntryMain(), fInputEHMix->CurrentEntryMix(), fInputEHMix->NumberMixed()));
   if (AliLog::GetDebugLevel("", IsA()->GetName()) > AliLog::kDebug) PrintEventInfo();
   // Post output data.
   PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskMixInfo::FinishTaskOutput()
{
   // FinishTaskOutput
   if (fMixInfo) fMixInfo->Print();
}


//________________________________________________________________________
void AliAnalysisTaskMixInfo::Terminate(Option_t *)
{
   // Draw result to the screen
   // Called once at the end of the query
   fOutputList = dynamic_cast<TList *>(GetOutputData(1));
   if (!fOutputList) {
      AliError("fOutputList not available");
      return;
   }
   fOutputList->Print();
   fMixInfo = (AliMixInfo *) fOutputList->FindObject("mixInfo");
   if (fMixInfo) {
      fMixInfo->Draw("HIST");
      AliMixEventPool *evPool = (AliMixEventPool *) fMixInfo->GetEventPool("mixEventPool");
      if (evPool) evPool->Print();
   }
}

//_____________________________________________________________________________
void AliAnalysisTaskMixInfo::SetLogType(AliLog::EType_t type, TString allClasses)
{
   //
   // Set Log level for this and other classes (list of their names)
   //
   AliDebug(AliLog::kDebug + 10, "<-");
   fLogType = type;
   fLogClassesString = allClasses;
   SetDebugForAllClasses();
   AliDebug(AliLog::kDebug + 10, "->");
}

//_____________________________________________________________________________
void AliAnalysisTaskMixInfo::SetDebugForAllClasses()
{
   //
   // Set debug level for all classes for which it is required
   //
   AliDebug(AliLog::kDebug + 10, "<-");
   TObjArray *array = fLogClassesString.Tokenize(":");
   TObjString *str;
   TString strr;
   for (Int_t i = 0; i < array->GetEntriesFast(); i++) {
      str = (TObjString *) array->At(i);
      strr = str->GetString();
      AliLog::SetClassDebugLevel(strr.Data(), fLogType);
      AliDebug(AliLog::kDebug + 5, Form("Setting Debug level %d to %s ...", (Int_t)fLogType - AliLog::kDebug, strr.Data()));
   }
   AliDebug(AliLog::kDebug + 10, "->");
}

//_____________________________________________________________________________
void AliAnalysisTaskMixInfo::InitMixInfo()
{
   //
   // Init mixing info
   //
   if (fInputEHMix) {
      fMixInfo = new AliMixInfo("mixInfo", "Mix title");
      AliMixEventPool *evPool = fInputEHMix->GetEventPool();
      if (!evPool) {
         //             TList *list = new TList;
         if (fMixInfo) fMixInfo->CreateHistogram(AliMixInfo::kMainEvents, 1, 1, 2);
         if (fMixInfo) fMixInfo->CreateHistogram(AliMixInfo::kMixedEvents, 1, 1, 2);
      } else {
         if (evPool->NeedInit()) evPool->Init();
         Int_t num = evPool->GetListOfEntryLists()->GetEntriesFast();
         if (fMixInfo) fMixInfo->CreateHistogram(AliMixInfo::kMainEvents, num, 1, num + 1);
         if (fMixInfo) fMixInfo->CreateHistogram(AliMixInfo::kMixedEvents, num, 1, num + 1);
      }
   }
}

//_____________________________________________________________________________
void AliAnalysisTaskMixInfo::PrintEventInfo()
{
   //
   // Prints event info for all event mxing cuts
   //
   if (fInputEHMix) {
      TObjArrayIter next(fInputEHMix->GetEventPool()->GetListOfEventCuts());
      AliMixEventCutObj *cut;
      AliInputEventHandler *ihMain = fInputEHMain->GetFirstInputEventHandler();
      AliMultiInputEventHandler *ihMultiMix = fInputEHMix->GetFirstMultiInputHandler();
      AliInputEventHandler *ihMix = 0;
      if (ihMultiMix) ihMix = ihMultiMix->GetFirstInputEventHandler();
      if (!ihMix) return;
      while ((cut = (AliMixEventCutObj *) next())) {
         cut->PrintValues(ihMain->GetEvent(), ihMix->GetEvent());
      }
   }
}

//_____________________________________________________________________________
void AliAnalysisTaskMixInfo::InitInputHandlers()
{
   //
   // Sets needed input handlers
   //
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   fInputEHMain = dynamic_cast<AliMultiInputEventHandler *>(mgr->GetInputEventHandler());
   if (fInputEHMain) {
      fInputEHMix = dynamic_cast<AliMixInputEventHandler *>(fInputEHMain->GetFirstMultiInputHandler());
   }
}
