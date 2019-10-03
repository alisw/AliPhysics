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
#include "AliLog.h"
#include "AliMixEventPool.h"
#include "AliMixEventInputHandler.h"
#include "AliAnalysisManager.h"

#include "AliMixInputHandlerInfo.h"
#include <TChainElement.h>

ClassImp(AliMixEventInputHandler)

//_____________________________________________________________________________
AliMixEventInputHandler::AliMixEventInputHandler(const Int_t size) :
   AliInputEventHandler(),
   fBufferSize(size),
   fInputHandlers(),
   fMixTrees(),
   fTreeMap(size),
   fMixIntupHandlerInfoTmp(0),
   fEntryCounter(0),
   fEventPool(0),
   fMixEventNumber(0) {
//
// Default constructor.
//
   AliDebug(AliLog::kDebug + 10, "<-");
   AliDebug(AliLog::kDebug + 10, "->");
}

//_____________________________________________________________________________
AliInputEventHandler *AliMixEventInputHandler::InputEventHandler(const Int_t index) {
   //
   // Returns input handler
   //
   AliDebug(AliLog::kDebug, Form("<-"));
   if ((index >= 0) && (index < fBufferSize)) {
      AliDebug(AliLog::kDebug, Form("->"));
      return (AliInputEventHandler *) fInputHandlers.At(index);
   }
   AliDebug(AliLog::kDebug, Form("->"));
   return 0;
}
//_____________________________________________________________________________
void AliMixEventInputHandler::SetInputHandlerForMixing(const AliInputEventHandler *const inHandler) {
   //
   // Create N (fBufferSize) copies of input handler
   //
   AliDebug(AliLog::kDebug, Form("<-"));
   AliDebug(AliLog::kDebug, Form("Creating %d input event handlers ...", fBufferSize));
   for (Int_t i = 0; i < fBufferSize; i++) {
      AliDebug(AliLog::kDebug + 5, Form("Adding %d ...", i));
      fInputHandlers.Add((AliInputEventHandler *) inHandler->Clone());
   }
   AliDebug(AliLog::kDebug, Form("->"));
}
//_____________________________________________________________________________
Bool_t AliMixEventInputHandler::Init(Option_t *opt) {
   //
   // Init() is called for all mix input handlers.
   //
   AliDebug(AliLog::kDebug, Form("<- \"%s\"", opt));
   AliInputEventHandler *eh = 0;
   TObjArrayIter next(&fInputHandlers);
   while ((eh = (AliInputEventHandler *) next())) {
      eh->Init(opt);
   }
   AliDebug(AliLog::kDebug, Form("->"));
   return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliMixEventInputHandler::Init(TTree *tree, Option_t *) {
   //
   // Init(const char*path) is called for all mix input handlers.
   // Create event pool if needed
   //
   AliDebug(AliLog::kDebug, Form("<- %p", tree));
   if (!tree) {
      AliDebug(AliLog::kDebug, Form("->"));
      return kFALSE;
   }
   AliDebug(AliLog::kDebug, Form("%s", tree->GetCurrentFile()->GetName()));

   // clears array of input handlers
   fMixTrees.SetOwner(kTRUE);
   fMixTrees.Clear();

   // create AliMixInputHandlerInfo
   if (!fMixIntupHandlerInfoTmp) fMixIntupHandlerInfoTmp = new AliMixInputHandlerInfo(tree->GetName());

   // adds current file
   fMixIntupHandlerInfoTmp->AddTreeToChain(tree);
   Int_t lastIndex = fMixIntupHandlerInfoTmp->GetChain()->GetListOfFiles()->GetEntries();
   TChainElement *che = (TChainElement *)fMixIntupHandlerInfoTmp->GetChain()->GetListOfFiles()->At(lastIndex - 1);
   AliMixInputHandlerInfo *mixIHI = 0;
   for (Int_t i = 0; i < fInputHandlers.GetEntries(); i++) {

      AliDebug(AliLog::kDebug, Form("fInputHandlers[%d]", i));
      mixIHI = new AliMixInputHandlerInfo(fMixIntupHandlerInfoTmp->GetName(), fMixIntupHandlerInfoTmp->GetTitle());
      mixIHI->PrepareEntry(che, -1, InputEventHandler(i));
      AliDebug(AliLog::kDebug, Form("chain[%d]->GetEntries() = %lld", i, mixIHI->GetChain()->GetEntries()));
      fMixTrees.Add(mixIHI);
   }
   AliDebug(AliLog::kDebug, Form("fEntryCounter=%lld", fEntryCounter));

   if (fEventPool->NeedInit())
      fEventPool->Init();

   AliDebug(AliLog::kDebug, Form("->"));
   return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliMixEventInputHandler::Notify() {
   //
   // Notify() is called for all mix input handlers
   //
   AliDebug(AliLog::kDebug, Form("<-"));
   AliInputEventHandler *eh = 0;
   TObjArrayIter next(&fInputHandlers);
   while ((eh = (AliInputEventHandler *) next())) {
      eh->Notify();
   }
   AliDebug(AliLog::kDebug, Form("->"));
   return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliMixEventInputHandler::Notify(const char *path) {
   //
   // Notify(const char*path) is called for all mix input handlers
   //
   AliDebug(AliLog::kDebug, Form("<- %s", path));
   AliInputEventHandler *eh = 0;
   TObjArrayIter next(&fInputHandlers);
   while ((eh = (AliInputEventHandler *) next())) {
      eh->Notify(path);
   }
   AliDebug(AliLog::kDebug, Form("->"));
   return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliMixEventInputHandler::BeginEvent(Long64_t entry) {
   //
   // BeginEvent(Long64_t entry) is called for all mix input handlers
   //

   AliDebug(AliLog::kDebug, Form("-> %lld", entry));
   AliInputEventHandler *eh = 0;
   TObjArrayIter next(&fInputHandlers);
   while ((eh = (AliInputEventHandler *) next())) {
      eh->BeginEvent(entry);
   }
   AliDebug(AliLog::kDebug, Form("->"));
   return kTRUE;
}


//_____________________________________________________________________________
Bool_t AliMixEventInputHandler::GetEntry() {
   //
   // Sets correct events to every mix events
   //
   AliDebug(AliLog::kDebug, Form("<-"));

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   AliInputEventHandler *inEvHMain = dynamic_cast<AliInputEventHandler *>(mgr->GetInputEventHandler());

   Long64_t zeroChainEntries = fMixIntupHandlerInfoTmp->GetChain()->GetEntries() - inEvHMain->GetTree()->GetEntries();

   // if fEntryCounter is 0 just add entry
   if (!fEntryCounter) {
      if (inEvHMain) {

         fEventPool->AddEntry(inEvHMain->GetTree()->GetReadEntry() + zeroChainEntries, inEvHMain->GetEvent());
      }
      return kTRUE;
   }

   AliDebug(AliLog::kDebug, Form("++++++++++++++ BEGIN SETUP EVENT %lld +++++++++++++++++++", fEntryCounter));

   fMixEventNumber = 0;
   TEntryList *el = fEventPool->FindEntryList(inEvHMain->GetEvent());
   Long64_t elNum = 0;
   if (el)
      elNum = el->GetN();
   else
      AliDebug(AliLog::kDebug, "el is null");

   AliInputEventHandler *eh = 0;
   AliMixInputHandlerInfo *mihi = 0;
   TObjArrayIter next(&fInputHandlers);
   Int_t counter = 0;
   Long64_t entryMix = 0;
   while ((eh = (AliInputEventHandler *) next())) {
      if (fEventPool->GetListOfEventCuts()->GetEntries() > 0) {
         entryMix = -1;
         if (el && el->GetN() >= fBufferSize) {
            Long64_t entryInEntryList =  elNum - 1 - counter;
            if (entryInEntryList < 0) break;
            entryMix = el->GetEntry(entryInEntryList);
         }
      } else {
         entryMix = fEntryCounter - 1 - counter ;
      }

      AliDebug(AliLog::kDebug, Form("Handler[%d] entryMix %lld ", counter, entryMix));
      if (entryMix < 0) break;

      mihi = (AliMixInputHandlerInfo *) fMixTrees.At(counter);
      TChainElement *te = fMixIntupHandlerInfoTmp->GetEntryInTree(entryMix);
      mihi->PrepareEntry(te, entryMix, InputEventHandler(counter));
      fMixEventNumber++;
      counter++;
   }

   if (inEvHMain) {
      fEventPool->AddEntry(inEvHMain->GetTree()->GetReadEntry() + zeroChainEntries, inEvHMain->GetEvent());
   }

   AliDebug(AliLog::kDebug, Form("fEntryCounter=%lld fMixEventNumber=%d", fEntryCounter, fMixEventNumber));
   AliDebug(AliLog::kDebug, Form("++++++++++++++ END SETUP EVENT %lld +++++++++++++++++++", fEntryCounter));
   AliDebug(AliLog::kDebug, Form("->"));
   return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliMixEventInputHandler::FinishEvent() {
   //
   // FinishEvent() is called for all mix input handlers
   //
   AliDebug(AliLog::kDebug, Form("<-"));
   AliInputEventHandler *eh = 0;
   TObjArrayIter next(&fInputHandlers);
   while ((eh = (AliInputEventHandler *) next())) {
      eh->FinishEvent();
   }
   fEntryCounter++;
   AliDebug(AliLog::kDebug, Form("->"));
   return kTRUE;
}
