//
// Class AliMixEventInputHandler
//
// Mixing input handler prepare N events before UserExec
// TODO example
// author:
//        Martin Vala (martin.vala@cern.ch)
//

#ifndef ALIMIXINPUTEVENTHANDLER_H
#define ALIMIXINPUTEVENTHANDLER_H

#include <TObjArray.h>
#include <TEntryList.h>
#include <TArrayI.h>

#include <AliVEvent.h>

#include "AliMultiInputEventHandler.h"

class TChain;
class TChainElement;
class AliMixEventPool;
class AliMixInputHandlerInfo;
class AliInputEventHandler;
class AliMixInputEventHandler : public AliMultiInputEventHandler {

public:
   AliMixInputEventHandler(const Int_t size = 1, const Int_t mixNum = 1);
   virtual ~AliMixInputEventHandler();
   // From the interface
   virtual Bool_t  Init(Option_t *opt) { return AliMultiInputEventHandler::Init(opt); }
   virtual Bool_t  Init(TTree *tree, Option_t *opt);
   virtual Bool_t  Notify();
   virtual Bool_t  Notify(const char *path);
   virtual Bool_t  BeginEvent(Long64_t entry);
   virtual Bool_t  GetEntry();
   virtual Bool_t  FinishEvent();

   // removing default impementation
   virtual void            AddInputEventHandler(AliVEventHandler */*inHandler*/);

   void                    SetInputHandlerForMixing(const AliInputEventHandler *const inHandler);
   void                    SetEventPool(AliMixEventPool *const evPool) { fEventPool = evPool; }

   AliMixEventPool        *GetEventPool() const { return fEventPool; }
   Int_t                   BufferSize() const { return fBufferSize; }
   Int_t                   NumberMixedTimes() const { return fNumberMixed; }
   Int_t                   MixNumber() const { return fMixNumber; }
   Long64_t                EntryAll() const { return fEntryCounter; }
   void                    UseDefaultProcess(Bool_t b = kTRUE) { fUseDefautProcess = b; }
   void                    DoMixExtra(Bool_t b = kTRUE) { fDoMixExtra = b; }
   void                    DoMixIfNotEnoughEvents(Bool_t b = kTRUE) { fDoMixIfNotEnoughEvents = b; }
   void                    SetMixNumber(const Int_t mixNum);

   void                    SetCurrentBinIndex(Int_t const index) { fCurrentBinIndex = index; }
   void                    SetCurrentEntry(Long64_t const entry) { fCurrentEntry = entry ; }
   void                    SetCurrentEntryMain(Long64_t const entry) { fCurrentEntryMain = entry ; }
   void                    SetCurrentEntryMix(Long64_t const entry) { fCurrentEntryMix = entry ; }
   void                    SetNumberMixed(Int_t const index) { fNumberMixed = index; }

   Int_t                   CurrentBinIndex() const { return fCurrentBinIndex; }
   Long64_t                CurrentEntry() const { return fCurrentEntry; }
   Long64_t                CurrentEntryMain() const { return fCurrentEntryMain; }
   Long64_t                CurrentEntryMix() const { return fCurrentEntryMix; }
   Int_t                   NumberMixed() const { return fNumberMixed; }

   void                    SelectCollisionCandidates(UInt_t offlineTriggerMask = AliVEvent::kMB) {fOfflineTriggerMask = offlineTriggerMask;}
   Bool_t                  IsEventCurrentSelected();
   Bool_t                  IsMixingIfNotEnoughEvents() { return fDoMixIfNotEnoughEvents;}

   void                    DoMixEventGetEntryAuto(Bool_t doAuto=kTRUE) { fDoMixEventGetEntryAuto = doAuto; }

   Bool_t                  GetEntryMainEvent();
   Bool_t                  GetEntryMixedEvent(Int_t idHandler=0);
protected:

   TObjArray               fMixTrees;              // buffer of input handlers
   TArrayI                 fTreeMap;               // tree map
   AliMixInputHandlerInfo *fMixIntupHandlerInfoTmp;//! mix input handler info full chain
   Long64_t                fEntryCounter;          // entry counter
   AliMixEventPool        *fEventPool;             // event pool
   Int_t                   fNumberMixed;           // number of mixed events with current event
   Int_t                   fMixNumber;             // user's mix number request

private:

   Bool_t                  fUseDefautProcess;      // use default process
   Bool_t                  fDoMixExtra;            // mix extra events to get enough combinations
   Bool_t                  fDoMixIfNotEnoughEvents;// mix events if they don't have enough events to mix
   Bool_t                  fDoMixEventGetEntryAuto;// flag for preparing mixed events automatically (default on)

   // mixing info
   Long64_t fCurrentEntry;       //! current entry number (adds 1 for every event processed on each worker)
   Long64_t fCurrentEntryMain;   //! current entry in chain of processed files
   Long64_t fCurrentEntryMix;    //! current mixed entry in chain of processed files
   Int_t    fCurrentBinIndex;    //! current bin index
   ULong64_t fOfflineTriggerMask;   //  Task processes collision candidates only

   TEntryList fCurrentMixEntry;    //! array of mix entries currently used (user should touch)
   Long64_t fCurrentEntryMainTree; //! current entry in current tree (main event)

   virtual Bool_t          MixStd();
   virtual Bool_t          MixBuffer();
   virtual Bool_t          MixEventsMoreTimesWithOneEvent();
   virtual Bool_t          MixEventsMoreTimesWithBuffer();

   void                    UserExecMixAllTasks(Long64_t entryCounter, Int_t idEntryList, Long64_t entryMainReal, Long64_t entryMixReal, Int_t numMixed);

   AliMixInputEventHandler(const AliMixInputEventHandler &handler);
   AliMixInputEventHandler &operator=(const AliMixInputEventHandler &handler);

   ClassDef(AliMixInputEventHandler, 5)
};

#endif
