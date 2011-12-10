//
// Class AliRsnInputHandler
//
// AliRsnInputHandler
// TODO example
// author:
//        Martin Vala (martin.vala@cern.ch)
//

#include <Riostream.h>
#include "AliLog.h"

#include "AliRsnEvent.h"
#include "AliMultiInputEventHandler.h"
#include "AliMixInputEventHandler.h"
#include "AliMCEventHandler.h"

#include "AliRsnCutSet.h"
#include "AliRsnInputHandler.h"
ClassImp(AliRsnInputHandler)

//_____________________________________________________________________________
AliRsnInputHandler::AliRsnInputHandler(const char *name) :
   AliInputEventHandler(name, name),
   fRsnEvent(0),
   fRsnSelector(),
   fRsnEventCuts(0)
{
//
// Default constructor.
//
   AliDebug(AliLog::kDebug + 10, "<-");
   AliDebug(AliLog::kDebug + 10, "->");
}

//_____________________________________________________________________________
AliRsnInputHandler::AliRsnInputHandler(const AliRsnInputHandler &copy) :
   AliInputEventHandler(),
   fRsnEvent(0),
   fRsnSelector(),
   fRsnEventCuts(copy.fRsnEventCuts)
{
//
// Default constructor.
//
   AliDebug(AliLog::kDebug + 10, "<-");
   AliDebug(AliLog::kDebug + 10, "->");
}

//_____________________________________________________________________________
AliRsnInputHandler &AliRsnInputHandler::operator=(const AliRsnInputHandler &copy)
{
//
// Default constructor.
//
   if (this == &copy)
      return *this;
   fRsnEventCuts = copy.fRsnEventCuts;
   return *this;
}


//_____________________________________________________________________________
AliRsnInputHandler::~AliRsnInputHandler()
{
//
// Destructor
//
   AliDebug(AliLog::kDebug + 10, "<-");
   delete fRsnEvent;
   AliDebug(AliLog::kDebug + 10, "->");
}

//_____________________________________________________________________________
Bool_t AliRsnInputHandler::Init(Option_t *opt)
{
//
// Init() is called for all mix input handlers.
//
   AliDebug(AliLog::kDebug + 5, Form("<- opt=%s", opt));

   AliDebug(AliLog::kDebug + 5, Form("->"));
   return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliRsnInputHandler::Init(TTree *tree, Option_t *opt)
{
//
// Init(const char*path) is called for all mix input handlers.
// Create event pool if needed
//
   AliDebug(AliLog::kDebug + 5, Form("<- %p %s opt=%s", (void *) tree, tree->GetName(), opt));
   AliDebug(AliLog::kDebug + 5, Form("->"));
   return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliRsnInputHandler::Notify()
{
//
// Notify() is called for all mix input handlers
//
   AliDebug(AliLog::kDebug + 5, Form("<-"));
   AliDebug(AliLog::kDebug + 5, Form("->"));
   return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliRsnInputHandler::Notify(const char *path)
{
//
// Notify(const char*path) is called for all mix input handlers
//
   AliDebug(AliLog::kDebug + 5, Form("<- %s", path));
   AliDebug(AliLog::kDebug + 5, "->");
   return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliRsnInputHandler::BeginEvent(Long64_t entry)
{
//
// BeginEvent(Long64_t entry) is called for all mix input handlers
//
   AliDebug(AliLog::kDebug + 5, Form("<- %lld", entry));

   if (fParentHandler) {
      TString tmp = "";
      AliInputEventHandler *ih = 0;
      AliMultiInputEventHandler *multiIH = dynamic_cast<AliMultiInputEventHandler *>(fParentHandler);
      if (multiIH) {
         ih = multiIH->GetFirstInputEventHandler();
         if (ih) {
            if (!fRsnEvent) fRsnEvent = new AliRsnEvent();
            fRsnEvent->SetRef(ih->GetEvent());
            fRsnEvent->SetPIDResponse(ih->GetPIDResponse());
            if (fRsnEvent->GetRefESD()) {
               AliMCEventHandler *mcH =  multiIH->GetFirstMCEventHandler();
               if (mcH) fRsnEvent->SetRefMC(mcH->MCEvent());
            } else if (fRsnEvent->GetRefAOD()) {
               AliAODEvent *aod = fRsnEvent->GetRefAOD();
               TClonesArray *listAOD = (TClonesArray *)(aod->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
               if (listAOD) fRsnEvent->SetRefMC(fRsnEvent->GetRefAOD());
            }
            if (fParentHandler->ParentHandler()) tmp = "MIX";
            // applying pid cuts
            //fRsnPIDManager.Reset();

            //fRsnPIDManager.ApplyCuts(fRsnEvent);
            fRsnSelector.Reset();

            // reject event if needed
            if (fRsnEventCuts) if (!fRsnEventCuts->IsSelected(fRsnEvent)) return kTRUE;
            fRsnSelector.ScanEvent(fRsnEvent);
         }
      }
   }
   AliDebug(AliLog::kDebug + 5, "->");
   return kTRUE;
}

Bool_t AliRsnInputHandler::GetEntry()
{
   AliDebug(AliLog::kDebug + 5, "<-");
   AliDebug(AliLog::kDebug + 5, "->");
   return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliRsnInputHandler::FinishEvent()
{
   //
   // FinishEvent() is called for all mix input handlers
   //
   AliDebug(AliLog::kDebug + 5, Form("<-"));
   AliDebug(AliLog::kDebug + 5, Form("->"));
   return kTRUE;
}
