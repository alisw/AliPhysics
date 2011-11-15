//
// Class AliMixEventInputHandler
//
// Mixing input handler prepare N events before UserExec
// TODO example
// author:
//        Martin Vala (martin.vala@cern.ch)
//

#include <TFile.h>
#include <TChain.h>
#include <TEntryList.h>
#include <TChainElement.h>
#include <TSystem.h>

#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

#include "AliMixEventPool.h"
#include "AliMixInputEventHandler.h"
#include "AliMixInputHandlerInfo.h"

#include "AliAnalysisTaskSE.h"

ClassImp(AliMixInputEventHandler)

AliMixInputEventHandler::AliMixInputEventHandler(const Int_t size, const Int_t mixNum): AliMultiInputEventHandler(size),
   fMixTrees(),
   fTreeMap(size > 0 ? size : 1),
   fMixIntupHandlerInfoTmp(0),
   fEntryCounter(0),
   fEventPool(0),
   fNumberMixed(0),
   fMixNumber(mixNum),
   fUseDefautProcess(kFALSE),
   fDoMixExtra(kTRUE),
   fDoMixIfNotEnoughEvents(kTRUE),
   fCurrentEntry(0),
   fCurrentEntryMain(0),
   fCurrentEntryMix(0),
   fCurrentBinIndex(-1),
   fOfflineTriggerMask(0)
{
   //
   // Default constructor.
   //
   AliDebug(AliLog::kDebug + 10, "<-");
   fMixTrees.SetOwner(kTRUE);
   SetMixNumber(mixNum);
   AliDebug(AliLog::kDebug + 10, "->");
}

//_____________________________________________________________________________
void AliMixInputEventHandler::SetInputHandlerForMixing(const AliInputEventHandler *const inHandler)
{
   //
   // Create N (fBufferSize) copies of input handler
   //
   AliDebug(AliLog::kDebug + 5, Form("<-"));
   fInputHandlers.Clear();
   AliDebug(AliLog::kDebug + 5, Form("Creating %d input event handlers ...", fBufferSize));
   for (Int_t i = 0; i < fBufferSize; i++) {
      AliDebug(AliLog::kDebug + 5, Form("Adding %d ...", i));
      fInputHandlers.Add((AliInputEventHandler *) inHandler->Clone());
   }

   AliDebug(AliLog::kDebug + 5, Form("->"));
}

//_____________________________________________________________________________
Bool_t AliMixInputEventHandler::Init(TTree *tree, Option_t *opt)
{
   //
   // Init(const char*path) is called for all mix input handlers.
   // Create event pool if needed
   //
   AliDebug(AliLog::kDebug + 5, Form("<- %p %s", (void *)tree, opt));
   fAnalysisType = opt;
   if (!tree) {
      AliDebug(AliLog::kDebug + 5, Form("->"));
      return kFALSE;
   }

   if (!fDoMixIfNotEnoughEvents) {
      fDoMixExtra = kFALSE;
      AliWarning("fDoMixIfNotEnoughEvents=kTRUE -> setting fDoMixExtra=kFALSE");
   }

   // clears array of input handlers
   fMixTrees.Clear();
   // create AliMixInputHandlerInfo
   if (!fMixIntupHandlerInfoTmp) {
      // loads first file TChain (tree)
      tree->LoadTree(0);
      fMixIntupHandlerInfoTmp = new AliMixInputHandlerInfo(tree->GetName());
   }

   AliInputEventHandler *ih = 0;
   for (Int_t i = 0; i < fInputHandlers.GetEntries(); i++) {
      ih = (AliInputEventHandler *) fInputHandlers.At(i);
      ih->SetParentHandler(this);
//       ih->Init(tree,opt);
   }

   AliDebug(AliLog::kDebug + 5, Form("->"));
   return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliMixInputEventHandler::Notify()
{
   //
   // Notify() is called for all mix input handlers
   //
   AliDebug(AliLog::kDebug + 5, Form("<-"));
   if (fUseDefautProcess) {
      AliDebug(AliLog::kDebug, Form("-> SKIPPED"));
      return AliMultiInputEventHandler::Notify();
   }
   AliDebug(AliLog::kDebug + 5, Form("->"));
   return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliMixInputEventHandler::Notify(const char *path)
{
   //
   // Notify(const char*path) is called for all mix input handlers
   //
   AliDebug(AliLog::kDebug + 5, Form("<- %s", path));

   Bool_t doPrepareEntry=kTRUE;
   TString anType = fAnalysisType;

   // in case of local doPrepareEntry only first time
   if (anType.CompareTo("proof")) doPrepareEntry = (fMixIntupHandlerInfoTmp->GetChain()->GetEntries()<=0);

   // adds current file
   fMixIntupHandlerInfoTmp->AddTreeToChain(path);
   Int_t lastIndex = fMixIntupHandlerInfoTmp->GetChain()->GetListOfFiles()->GetEntries();
   TChainElement *che = (TChainElement *)fMixIntupHandlerInfoTmp->GetChain()->GetListOfFiles()->At(lastIndex - 1);
   AliMixInputHandlerInfo *mixIHI = 0;
   for (Int_t i = 0; i < fInputHandlers.GetEntries(); i++) {
      AliDebug(AliLog::kDebug + 5, Form("fInputHandlers[%d]", i));
      mixIHI = new AliMixInputHandlerInfo(fMixIntupHandlerInfoTmp->GetName(), fMixIntupHandlerInfoTmp->GetTitle());
//       mixIHI->PrepareEntry(che, -1, (AliInputEventHandler *)InputEventHandler(i), fAnalysisType);
      if (doPrepareEntry) mixIHI->PrepareEntry(che, -1, (AliInputEventHandler *)InputEventHandler(i), fAnalysisType);
      AliDebug(AliLog::kDebug + 5, Form("chain[%d]->GetEntries() = %lld", i, mixIHI->GetChain()->GetEntries()));
      fMixTrees.Add(mixIHI);
   }
   AliDebug(AliLog::kDebug + 5, Form("fEntryCounter=%lld", fEntryCounter));
   if (fEventPool && fEventPool->NeedInit())
      fEventPool->Init();
   if (fUseDefautProcess) {
      AliDebug(AliLog::kDebug, Form("-> SKIPPED"));
      return AliMultiInputEventHandler::Notify(path);
   }
   AliDebug(AliLog::kDebug + 5, Form("->"));
   return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliMixInputEventHandler::BeginEvent(Long64_t entry)
{
   //
   // BeginEvent(Long64_t entry) is called for all mix input handlers
   //
   AliDebug(AliLog::kDebug + 5, Form("-> %lld", entry));
   if (fUseDefautProcess) {
      AliDebug(AliLog::kDebug, Form("-> SKIPPED"));
      AliMultiInputEventHandler::BeginEvent(entry);/* return GetEntry();*/
   }
   AliDebug(AliLog::kDebug + 5, Form("->"));
   return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliMixInputEventHandler::GetEntry()
{
   //
   // All mixed events are set
   //
   AliDebug(AliLog::kDebug + 5, Form("<-"));

   if (!fEventPool) {
      MixStd();
   }
   // if buffer size is higher then 1
   else if (fBufferSize > 1) {
      MixBuffer();
   }
   // if mix number is higher then 0 and buffer size is 1
   else if (fMixNumber > 0) {
      MixEventsMoreTimesWithOneEvent();
   } else {
      AliWarning("Not supported Mixing !!!");
   }

   AliDebug(AliLog::kDebug + 5, Form("->"));
   return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliMixInputEventHandler::MixStd()
{
   //
   // Mix std - No event pool
   //
   AliDebug(AliLog::kDebug + 5, Form("<-"));
   AliDebug(AliLog::kDebug + 1, "Mix method");
   // get correct handler
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   AliMultiInputEventHandler *mh = dynamic_cast<AliMultiInputEventHandler *>(mgr->GetInputEventHandler());
   AliInputEventHandler *inEvHMain = 0;
   if (mh) inEvHMain = dynamic_cast<AliInputEventHandler *>(mh->GetFirstInputEventHandler());
   else inEvHMain = dynamic_cast<AliInputEventHandler *>(mgr->GetInputEventHandler());
   if (!inEvHMain) return kFALSE;

   // check for PhysSelection
   if (!IsEventCurrentSelected()) return kFALSE;

   // return in case of 0 entry in full chain
   if (!fEntryCounter) {
      AliDebug(AliLog::kDebug + 3, Form("-> fEntryCounter == 0"));
      // runs UserExecMix for all tasks, if needed
      UserExecMixAllTasks(fEntryCounter, 1, fEntryCounter, -1, 0);
      return kTRUE;
   }
   // pre mix evetns
   Int_t mixNum = fMixNumber;
   if (fDoMixExtra) {
      if (fEntryCounter <= 2 * fMixNumber) mixNum = 2 * fMixNumber + 2;
   }
   // start of
   AliDebug(AliLog::kDebug + 3, Form("++++++++++++++ BEGIN SETUP EVENT %lld +++++++++++++++++++", fEntryCounter));
   // reset mix number
   fNumberMixed = 0;
   AliMixInputHandlerInfo *mihi = 0;
   Long64_t entryMix = 0, entryMixReal = 0;
   Int_t counter = 0;
   for (counter = 0; counter < mixNum; counter++) {
      entryMix = fEntryCounter - 1 - counter ;
      AliDebug(AliLog::kDebug + 5, Form("Handler[%d] entryMix %lld ", counter, entryMix));
      if (entryMix < 0) break;
      entryMixReal = entryMix;
      mihi = (AliMixInputHandlerInfo *) fMixTrees.At(0);
      TChainElement *te = fMixIntupHandlerInfoTmp->GetEntryInTree(entryMix);
      if (!te) {
         AliError("te is null. this is error. tell to developer (#1)");
      } else {
         mihi->PrepareEntry(te, entryMix, (AliInputEventHandler *)InputEventHandler(0), fAnalysisType);
         // runs UserExecMix for all tasks
         fNumberMixed++;
         UserExecMixAllTasks(fEntryCounter, 1, fEntryCounter, entryMixReal, fNumberMixed);
         InputEventHandler(0)->FinishEvent();
      }
   }
   AliDebug(AliLog::kDebug + 3, Form("fEntryCounter=%lld fMixEventNumber=%d", fEntryCounter, fNumberMixed));
   AliDebug(AliLog::kDebug + 3, Form("++++++++++++++ END SETUP EVENT %lld +++++++++++++++++++", fEntryCounter));
   AliDebug(AliLog::kDebug + 5, Form("->"));
   AliDebug(AliLog::kDebug + 5, Form("->"));
   return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliMixInputEventHandler::MixBuffer()
{
   //
   // Mix in event buffer
   //
   AliDebug(AliLog::kDebug + 5, Form("<-"));
   AliDebug(AliLog::kDebug + 1, "Mix method");
   // get correct handler
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   AliMultiInputEventHandler *mh = dynamic_cast<AliMultiInputEventHandler *>(mgr->GetInputEventHandler());
   AliInputEventHandler *inEvHMain = 0;
   if (mh) inEvHMain = dynamic_cast<AliInputEventHandler *>(mh->GetFirstInputEventHandler());
   else inEvHMain = dynamic_cast<AliInputEventHandler *>(mgr->GetInputEventHandler());
   if (!inEvHMain) return kFALSE;

   // check for PhysSelection
   if (!IsEventCurrentSelected()) return kFALSE;

   // find out zero chain entries
   Long64_t zeroChainEntries = fMixIntupHandlerInfoTmp->GetChain()->GetEntries() - inEvHMain->GetTree()->GetTree()->GetEntries();
   // fill entry
   Long64_t currentMainEntry = inEvHMain->GetTree()->GetTree()->GetReadEntry() + zeroChainEntries;
   // fills entry
   if (fEventPool && inEvHMain) fEventPool->AddEntry(currentMainEntry, inEvHMain->GetEvent());
   // start of
   AliDebug(AliLog::kDebug + 3, Form("++++++++++++++ BEGIN SETUP EVENT %lld +++++++++++++++++++", fEntryCounter));
   // reset mix number
   fNumberMixed = 0;
   Long64_t elNum = 0;
   TEntryList *el = 0;
   Int_t idEntryList = -1;
   if (fEventPool) el = fEventPool->FindEntryList(inEvHMain->GetEvent(), idEntryList);
   // return in case of 0 entry in full chain
   if (!fEntryCounter) {
      AliDebug(AliLog::kDebug + 3, Form("-> fEntryCounter == 0"));
      // runs UserExecMix for all tasks, if needed
      if (el) UserExecMixAllTasks(fEntryCounter, idEntryList, currentMainEntry, -1, 0);
      else UserExecMixAllTasks(fEntryCounter, -1, currentMainEntry, -1, 0);
      return kTRUE;
   }
   if (!el) {
      AliDebug(AliLog::kDebug + 3, Form("++++++++++++++ END SETUP EVENT %lld SKIPPED (el null) +++++++++++++++++++", fEntryCounter));
      UserExecMixAllTasks(fEntryCounter, -1, fEntryCounter, -1, 0);
      return kTRUE;
   } else {
      elNum = el->GetN();
      if (elNum < fBufferSize + 1) {
         UserExecMixAllTasks(fEntryCounter, idEntryList, currentMainEntry, -1, 0);
         AliDebug(AliLog::kDebug + 3, Form("++++++++++++++ END SETUP EVENT %lld SKIPPED (%lld) LESS THEN BUFFER +++++++++++++++++++", fEntryCounter, elNum));
         return kTRUE;
      }
   }
   AliMixInputHandlerInfo *mihi = 0;
   Long64_t entryMix = 0, entryMixReal = 0;
   Int_t counter = 0;
   AliInputEventHandler *eh = 0;
   TObjArrayIter next(&fInputHandlers);
   while ((eh = dynamic_cast<AliInputEventHandler *>(next()))) {
      if (fEventPool && fEventPool->GetListOfEventCuts()->GetEntries() > 0) {
         entryMix = -1;
         if (elNum >= fBufferSize) {
            Long64_t entryInEntryList =  elNum - 2 - counter;
            if (entryInEntryList < 0) break;
            entryMix = el->GetEntry(entryInEntryList);
         }
      }
      AliDebug(AliLog::kDebug + 5, Form("Handler[%d] entryMix %lld ", counter, entryMix));
      if (entryMix < 0) {
         UserExecMixAllTasks(fEntryCounter, -1, currentMainEntry, -1, 0);
         break;
      }
      entryMixReal = entryMix;
      mihi = (AliMixInputHandlerInfo *) fMixTrees.At(counter);
      TChainElement *te = fMixIntupHandlerInfoTmp->GetEntryInTree(entryMix);
      if (!te) {
         AliError("te is null. this is error. tell to developer (#1)");
      } else {
         AliDebug(AliLog::kDebug + 3, Form("Preparing InputEventHandler(%d)", counter));
         mihi->PrepareEntry(te, entryMix, (AliInputEventHandler *)InputEventHandler(counter), fAnalysisType);
         // runs UserExecMix for all tasks
         UserExecMixAllTasks(fEntryCounter, idEntryList, fEntryCounter, entryMixReal, counter);
         fNumberMixed++;
      }
      counter++;
   }
   AliDebug(AliLog::kDebug + 3, Form("fEntryCounter=%lld fMixEventNumber=%d", fEntryCounter, fNumberMixed));
   AliDebug(AliLog::kDebug + 3, Form("++++++++++++++ END SETUP EVENT %lld +++++++++++++++++++", fEntryCounter));
   AliDebug(AliLog::kDebug + 5, Form("->"));
   AliDebug(AliLog::kDebug + 5, Form("->"));
   return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliMixInputEventHandler::MixEventsMoreTimesWithOneEvent()
{
   //
   // Mix in history with one event in buffer
   //
   AliDebug(AliLog::kDebug + 5, "<-");
   AliDebug(AliLog::kDebug + 1, "Mix method");
   // get correct handler
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   AliMultiInputEventHandler *mh = dynamic_cast<AliMultiInputEventHandler *>(mgr->GetInputEventHandler());
   AliInputEventHandler *inEvHMain = 0;
   if (mh) inEvHMain = dynamic_cast<AliInputEventHandler *>(mh->GetFirstInputEventHandler());
   else inEvHMain = dynamic_cast<AliInputEventHandler *>(mgr->GetInputEventHandler());
   if (!inEvHMain) return kFALSE;

   // check for PhysSelection
   if (!IsEventCurrentSelected()) return kFALSE;

   // find out zero chain entries
   Long64_t zeroChainEntries = fMixIntupHandlerInfoTmp->GetChain()->GetEntries() - inEvHMain->GetTree()->GetTree()->GetEntries();
   // fill entry
   Long64_t currentMainEntry = inEvHMain->GetTree()->GetTree()->GetReadEntry() + zeroChainEntries;
   if (fEventPool && inEvHMain) fEventPool->AddEntry(currentMainEntry, inEvHMain->GetEvent());
   // start of
   AliDebug(AliLog::kDebug + 3, Form("++++++++++++++ BEGIN SETUP EVENT %lld +++++++++++++++++++", fEntryCounter));
   // reset mix number
   fNumberMixed = 0;
   Long64_t elNum = 0;
   Int_t idEntryList = -1;
   TEntryList *el = 0;
   if (fEventPool) el = fEventPool->FindEntryList(inEvHMain->GetEvent(), idEntryList);
   // return in case of 0 entry in full chain
   if (!fEntryCounter) {
      // runs UserExecMix for all tasks, if needed
      if (el && fDoMixIfNotEnoughEvents) {
         UserExecMixAllTasks(fEntryCounter, idEntryList, currentMainEntry, -1, 0);
      } else {
         idEntryList = -1;
         UserExecMixAllTasks(fEntryCounter, idEntryList, currentMainEntry, -1, 0);
      }
      AliDebug(AliLog::kDebug + 3, Form("++++++++++++++ END SETUP EVENT %lld SKIPPED (fEntryCounter=0, idEntryList=%d) +++++++++++++++++++", fEntryCounter, idEntryList));
      return kTRUE;
   }
   if (!el) {
      if (fEventPool) {
         AliDebug(AliLog::kDebug + 3, Form("++++++++++++++ END SETUP EVENT %lld SKIPPED (el null, idEntryList=%d) +++++++++++++++++++", fEntryCounter, idEntryList));
         UserExecMixAllTasks(fEntryCounter, -1, currentMainEntry, -1, 0);
         return kTRUE;
      }
   } else {
      elNum = el->GetN();
      if (elNum < fBufferSize + 1) {
         if (fDoMixIfNotEnoughEvents) {
            // include main event in to counter in this case (so idEntryList>0)
            UserExecMixAllTasks(fEntryCounter, idEntryList, currentMainEntry, -1, 0);
         }  else {
            // dont include it in main event counter (idEntryList = -1)
            idEntryList = -1;
            UserExecMixAllTasks(fEntryCounter, idEntryList, currentMainEntry, -1, 0);
         }
         AliDebug(AliLog::kDebug + 3, Form("++++++++++++++ END SETUP EVENT %lld SKIPPED [FIRST ENTRY in el] (elnum=%lld, idEntryList=%d) +++++++++++++++++++", fEntryCounter, elNum, idEntryList));
         return kTRUE;
      }
      if (!fDoMixIfNotEnoughEvents) {
         if (elNum <= fMixNumber + 1) {
            UserExecMixAllTasks(fEntryCounter, idEntryList, currentMainEntry, -1, 0);
            AliDebug(AliLog::kDebug + 3, Form("++++++++++++++ END SETUP EVENT %lld SKIPPED (%lld) NOT ENOUGH EVENTS TO MIX => NEED=%d +++++++++++++++++++", fEntryCounter, elNum, fMixNumber + 1));
            return kTRUE;
         }
      }
   }
   // pre mix evetns
   Int_t mixNum = fMixNumber;
   if (fDoMixExtra) {
      if (elNum <= 2 * fMixNumber + 1) mixNum = elNum + 1;
   }
   AliMixInputHandlerInfo *mihi = 0;
   Long64_t entryMix = 0, entryMixReal = 0;
   Int_t counter = 0;
   mihi = (AliMixInputHandlerInfo *) fMixTrees.At(0);
   // fills num for main events
   for (counter = 0; counter < mixNum; counter++) {
      Long64_t entryInEntryList =  elNum - 2 - counter;
      AliDebug(AliLog::kDebug + 3, Form("entryInEntryList=%lld", entryInEntryList));
      if (entryInEntryList < 0) break;
      entryMix = el->GetEntry(entryInEntryList);
      AliDebug(AliLog::kDebug + 3, Form("entryMix=%lld", entryMix));
      if (entryMix < 0) break;
      entryMixReal = entryMix;
      TChainElement *te = fMixIntupHandlerInfoTmp->GetEntryInTree(entryMix);
      if (!te) {
         AliError("te is null. this is error. tell to developer (#2)");
      } else {
         mihi->PrepareEntry(te, entryMix, (AliInputEventHandler *)InputEventHandler(0), fAnalysisType);
         // runs UserExecMix for all tasks
         fNumberMixed++;
         UserExecMixAllTasks(fEntryCounter, idEntryList, currentMainEntry, entryMixReal, fNumberMixed);
         InputEventHandler(0)->FinishEvent();
      }
   }
   AliDebug(AliLog::kDebug + 3, Form("fEntryCounter=%lld fMixEventNumber=%d", fEntryCounter, fNumberMixed));
   AliDebug(AliLog::kDebug + 3, Form("++++++++++++++ END SETUP EVENT %lld +++++++++++++++++++", fEntryCounter));
   AliDebug(AliLog::kDebug + 5, Form("->"));
   AliDebug(AliLog::kDebug + 5, Form("->"));
   return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliMixInputEventHandler::MixEventsMoreTimesWithBuffer()
{
   //
   // Mix more events in buffer with mixing with history
   //
   AliWarning("Not implemented");
   return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliMixInputEventHandler::FinishEvent()
{
   //
   // FinishEvent() is called for all mix input handlers
   //
   AliDebug(AliLog::kDebug + 5, Form("<-"));
   AliMultiInputEventHandler::FinishEvent();
   fEntryCounter++;
   AliDebug(AliLog::kDebug + 5, Form("->"));
   return kTRUE;
}

//_____________________________________________________________________________
void AliMixInputEventHandler::AddInputEventHandler(AliVEventHandler *)
{
   //
   // AddInputEventHandler will not be used
   //
   AliWarning("Function AddInputEventHandler is disabled for AliMixEventInputHandler !!!");
   AliWarning("Use AliMixEventInputHandler::SetInputHandlerForMixing instead. Exiting ...");
}

//_____________________________________________________________________________
void AliMixInputEventHandler::UserExecMixAllTasks(Long64_t entryCounter, Int_t idEntryList, Long64_t entryMainReal, Long64_t entryMixReal, Int_t numMixed)
{
   //
   // Execute all task and sets mixing parameters
   //
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   AliAnalysisTaskSE *mixTask = 0;
   TObjArrayIter next(mgr->GetTasks());
   while ((mixTask = dynamic_cast<AliAnalysisTaskSE *>(next()))) {
         AliDebug(AliLog::kDebug, Form("%s %lld %d [%lld,%lld] %d", mixTask->GetName(), entryCounter, numMixed, entryMainReal, entryMixReal, idEntryList));
         fCurrentEntry = entryCounter;
         fCurrentEntryMain = entryMainReal;
         fCurrentEntryMix = entryMixReal;
         fCurrentBinIndex = idEntryList;
         if (entryMixReal >= 0) mixTask->UserExecMix("");
   }
}

//_____________________________________________________________________________
void AliMixInputEventHandler::SetMixNumber(const Int_t mixNum)
{
   //
   // Sets mix number
   //
   if (fMixNumber > 1 && fBufferSize > 1) {
      AliWarning("Sleeping 10 sec to show Warning Message ...");
      AliWarning("=========================================================================================");
      AliWarning(Form("BufferSize(%d) higher > 1 and fMixNumber(%d) > 1, which is not supported", fBufferSize, mixNum));
      AliWarning("");
      AliWarning("\tBufferSize will be set to 1");
      AliWarning("");
      AliWarning("Hints:");
      AliWarning("");
      AliWarning("\t1.If you want to use buffer do:");
      AliWarning(Form("\t\tAliMixInputEventHandler *mixH = new AliMixInputEventHandler(%d,1)", fBufferSize));
      AliWarning("");
      AliWarning("\t2.If you want to use mix more time with buffer size 1, then do:");
      AliWarning(Form("\t\tAliMixInputEventHandler *mixH = new AliMixInputEventHandler(1,%d)", mixNum));
      AliWarning("");
      AliWarning("=========================================================================================");
      gSystem->Sleep(10000);
      fBufferSize = 1;
   }
   fMixNumber = mixNum;
}

Bool_t AliMixInputEventHandler::IsEventCurrentSelected()
{
   AliDebug(AliLog::kDebug + 5, Form("<-"));
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   AliMultiInputEventHandler *mh = dynamic_cast<AliMultiInputEventHandler *>(mgr->GetInputEventHandler());
   Bool_t isSelected = kTRUE;
   if (mh) {
      if (fOfflineTriggerMask && mh->GetEventSelection()) {
         isSelected = fOfflineTriggerMask & mh->IsEventSelected();
      }
   }
   AliDebug(AliLog::kDebug + 1, Form("isSelected=%d", isSelected));
   AliDebug(AliLog::kDebug + 5, Form("-> %d", isSelected));
   return isSelected;
}
