//
// Class AliRsnInputHandler
//
// AliRsnInputHandler
// TODO example
// author:
//        Martin Vala (martin.vala@cern.ch)
//

#ifndef ALIRSNINPUTHANDLER_H
#define ALIRSNINPUTHANDLER_H

#include "AliInputEventHandler.h"
//#include "AliRsnPIDManager.h"
#include "AliRsnDaughterSelector.h"

class AliRsnEvent;
class AliRsnCutSet;
class AliRsnInputHandler : public AliInputEventHandler {

public:
   AliRsnInputHandler(const char *name = "rsnInputHandler");
   AliRsnInputHandler(const AliRsnInputHandler& handler);
   AliRsnInputHandler &operator=(const AliRsnInputHandler &handler);
   virtual ~AliRsnInputHandler();

   // From the interface
   virtual Bool_t  Init(Option_t *opt);
   virtual Bool_t  Init(TTree *tree, Option_t *opt);
   virtual Bool_t  BeginEvent(Long64_t entry);
   virtual Bool_t  FinishEvent();
   virtual Bool_t  Notify();
   virtual Bool_t  Notify(const char *path);
   virtual Bool_t  GetEntry();

   void SetEventCuts(AliRsnCutSet *cuts) {fRsnEventCuts = cuts;}
   AliRsnEvent *GetRsnEvent() { return fRsnEvent; }
   //AliRsnPIDManager *GetPIDManager() { return &fRsnPIDManager; }
   AliRsnDaughterSelector *GetSelector() { return &fRsnSelector; }

private:

   AliRsnEvent *fRsnEvent;               //! pointer to current event
   AliRsnDaughterSelector fRsnSelector;  //  pointer to selector for building lists
   AliRsnCutSet *fRsnEventCuts;          //  unique set of cuts to be used for whole analysis

   ClassDef(AliRsnInputHandler, 1)
};

#endif
