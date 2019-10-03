//
// Class AliMixEventInputHandler
//
// Mixing input handler prepare N events before UserExec
// TODO example
// author:
//        Martin Vala (martin.vala@cern.ch)
//

#ifndef ALIMIXEVENTINPUTHANDLER_H
#define ALIMIXEVENTINPUTHANDLER_H

#include <TObjArray.h>

#include "AliInputEventHandler.h"
#include <TArrayI.h>

class TChain;
class AliMixEventPool;
class AliMixInputHandlerInfo;
class AliMixEventInputHandler : public AliInputEventHandler {

public:
   AliMixEventInputHandler(const Int_t size=1);
   virtual ~AliMixEventInputHandler() {;};

   // From the interface
   virtual Bool_t  Init(Option_t *opt /*opt*/);
   virtual Bool_t  Init(TTree *tree, Option_t* /*opt*/);
   virtual Bool_t  BeginEvent(Long64_t entry /*entry*/);
   virtual Bool_t  GetEntry();
   virtual Bool_t  FinishEvent();
   virtual Bool_t  Notify();
   virtual Bool_t  Notify(const char *path);


   void SetInputHandlerForMixing(const AliInputEventHandler *const inHandler);
   void SetEventPool(AliMixEventPool *const evPool) {fEventPool = evPool;}

   AliInputEventHandler *InputEventHandler(const Int_t index);
   AliMixEventPool *GetEventPool() const { return fEventPool;}
   Int_t           BufferSize() const {return fBufferSize;}
   Int_t           MixedEventNumber() const {return fMixEventNumber;}
   Long64_t        EntryAll() const { return fEntryCounter;}
protected:

   Int_t             fBufferSize;          // Size of the buffer
   TObjArray         fInputHandlers;       // buffer of input handlers
   TObjArray         fMixTrees;            // buffer of input handlers
   TArrayI           fTreeMap;             // tree map
   AliMixInputHandlerInfo *fMixIntupHandlerInfoTmp;    //! mix input handler info full chain
   Long64_t          fEntryCounter;        // entry counter
   AliMixEventPool  *fEventPool;           // event pool

   Int_t             fMixEventNumber;      // number mix

private:

   AliMixEventInputHandler(const AliMixEventInputHandler &handler);
   AliMixEventInputHandler &operator=(const AliMixEventInputHandler &handler);

   ClassDef(AliMixEventInputHandler, 1)
};

#endif
