//
// Computator for single daughters.
// Implements a simple loop on tracks from one of the entry lists
// filled by the task AliRsnInputHandler, adding a check on their
// definition specified in the daughter def.
//

#include <TEntryList.h>

#include "AliLog.h"

#include "AliRsnEvent.h"
#include "AliRsnDaughterDef.h"
#include "AliRsnDaughterSelector.h"

#include "AliRsnLoopDaughter.h"

ClassImp(AliRsnLoopDaughter)

//_____________________________________________________________________________
AliRsnLoopDaughter::AliRsnLoopDaughter(const char *name, Int_t listID, AliRsnDaughterDef *def) :
   AliRsnLoop(name),
   fOnlyTrue(kFALSE),
   fListID(listID),
   fDef(def),
   fDaughter()
{
//
// Default constructor
//
}

//_____________________________________________________________________________
AliRsnLoopDaughter::AliRsnLoopDaughter(const AliRsnLoopDaughter& copy) :
   AliRsnLoop(copy),
   fOnlyTrue(copy.fOnlyTrue),
   fListID(copy.fListID),
   fDef(copy.fDef),
   fDaughter(copy.fDaughter)
{
//
// Copy constructor
//
}

//_____________________________________________________________________________
AliRsnLoopDaughter& AliRsnLoopDaughter::operator=(const AliRsnLoopDaughter& copy)
{
//
// Assignment operator
//

   AliRsnLoop::operator=(copy);
   fOnlyTrue = copy.fOnlyTrue;
   fListID = copy.fListID;
   fDaughter = copy.fDaughter;
   fDef = copy.fDef;

   return (*this);
}

//_____________________________________________________________________________
AliRsnLoopDaughter::~AliRsnLoopDaughter()
{
//
// Destructor
//
}

//_____________________________________________________________________________
void AliRsnLoopDaughter::Print(Option_t* /*option*/) const
{
//
// Prints info about pair
//
}

//_____________________________________________________________________________
Bool_t AliRsnLoopDaughter::Init(const char *prefix, TList *list)
{
//
// Initialization function.
// Loops on all functions and eventual the ntuple, to initialize output objects.
//

   return AliRsnLoop::Init(Form("%s_%s", prefix, GetName()), list);
}

//_____________________________________________________________________________
Int_t AliRsnLoopDaughter::DoLoop
(AliRsnEvent *evMain, AliRsnDaughterSelector *selMain, AliRsnEvent *, AliRsnDaughterSelector *)
{
//
// Loop function.
// Computes what is needed from passed events.
// Returns the number of pairs successfully processed.
//

   if (!OkEvent(evMain)) return 0;

   Int_t i, il, nadd = 0, nlist = 0;
   TEntryList *list[2] = {0, 0};
   
   if (fDef->IsChargeDefined()) {
      list[0] = selMain->GetSelected(fListID, fDef->GetChargeC());
      list[1] = 0x0;
      nlist = 1;
   } else {
      list[0] = selMain->GetSelected(fListID, '+');
      if (list[0]) {
         list[1] = selMain->GetSelected(fListID, '-');
         nlist = 2;
      } else {
         list[0] = selMain->GetSelected(fListID, '0');
         list[1] = 0x0;
         nlist = 1;
      }
   }
   
   TObjArrayIter next(&fOutputs);
   AliRsnListOutput *out = 0x0;
   
   for (il = 0; il < nlist; il++) {
      if (!list[il]) {
         AliError(Form("List #%d is null", il));
         continue;
      }
      for (i = 0; i < list[il]->GetN(); i++) {
         evMain->SetDaughter(fDaughter, (Int_t)list[il]->GetEntry(i));
         // check matching
         if (fOnlyTrue && !fDef->MatchesPID(&fDaughter)) continue;
         if (!fDef->MatchesCharge(&fDaughter)) continue;
         if (!fDef->MatchesRefType(&fDaughter)) continue;
         // fill outputs
         nadd++;
         next.Reset();
         while ( (out = (AliRsnListOutput*)next()) ) {
            out->Fill(&fDaughter);
         }
      }
   }
   
   return nadd;
}
