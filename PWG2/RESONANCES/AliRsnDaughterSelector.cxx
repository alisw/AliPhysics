#include <Riostream.h>
#include <TEntryList.h>

#include "AliLog.h"

#include "AliRsnCutSet.h"
#include "AliRsnDaughterDef.h"

#include "AliRsnDaughterSelector.h"

ClassImp(AliRsnDaughterSelector)

//__________________________________________________________________________________________________
AliRsnDaughterSelector::AliRsnDaughterSelector(const char *name, const char *title) : 
   TNamed(name, title),
   fCutSetsN("AliRsnCutSet", 0),
   fCutSetsC("AliRsnCutSet", 0),
   fEntryListsN("TEntryList", 0),
   fEntryListsP("TEntryList", 0),
   fEntryListsM("TEntryList", 0)
{
//
// Default constructor.
//

   AliDebug(AliLog::kDebug + 10, "<-");
   AliDebug(AliLog::kDebug + 10, "->");
}

//__________________________________________________________________________________________________
AliRsnDaughterSelector::AliRsnDaughterSelector(const AliRsnDaughterSelector& copy) : 
   TNamed(copy),
   fCutSetsN(copy.fCutSetsN),
   fCutSetsC(copy.fCutSetsC),
   fEntryListsN(copy.fEntryListsN),
   fEntryListsP(copy.fEntryListsP),
   fEntryListsM(copy.fEntryListsM)
{
//
// Copy constructor.
//

   AliDebug(AliLog::kDebug + 10, "<-");
   AliDebug(AliLog::kDebug + 10, "->");
}

//__________________________________________________________________________________________________
AliRsnDaughterSelector& AliRsnDaughterSelector::operator=(const AliRsnDaughterSelector& copy)
{
//
// Copy constructor.
//

   AliDebug(AliLog::kDebug + 10, "<-");

   TNamed::operator=(copy);
   
   fCutSetsN = copy.fCutSetsN;
   fCutSetsC = copy.fCutSetsC;
   fEntryListsN = copy.fEntryListsN;
   fEntryListsP = copy.fEntryListsP;
   fEntryListsM = copy.fEntryListsM;

   AliDebug(AliLog::kDebug + 10, "->");
   
   return (*this);
}

//__________________________________________________________________________________________________
AliRsnDaughterSelector::~AliRsnDaughterSelector()
{
//
// Destructor
//

   AliDebug(AliLog::kDebug + 10, "<-");
   
   fCutSetsN.Delete();
   fCutSetsC.Delete();
   fEntryListsN.Delete();
   fEntryListsP.Delete();
   fEntryListsM.Delete();
   
   AliDebug(AliLog::kDebug + 10, "->");
}

//__________________________________________________________________________________________________
void AliRsnDaughterSelector::Print(Option_t* option) const
{
//
// Override TObject::Print()
//

   TNamed::Print(option);
   
   Int_t i, nSets;
   AliRsnCutSet *set = 0x0;
   TEntryList *list = 0x0;
   
   // neutral
   nSets = fCutSetsN.GetEntries();
   for (i = 0; i < nSets; i++) {
      set = (AliRsnCutSet*)fCutSetsN[i];
      list = (TEntryList*)fEntryListsN[i];
      AliInfo(Form("Neutral entry list for cut set '%s' has %d entries", set->GetName(), (Int_t)list->GetN()));
   }
   
   // charged
   nSets = fCutSetsC.GetEntries();
   for (i = 0; i < nSets; i++) {
      set = (AliRsnCutSet*)fCutSetsC[i];
      list = (TEntryList*)fEntryListsP[i];
      AliInfo(Form("Positive entry list for cut set '%s' has %d entries", set->GetName(), (Int_t)list->GetN()));
      list = (TEntryList*)fEntryListsM[i];
      AliInfo(Form("Negative entry list for cut set '%s' has %d entries", set->GetName(), (Int_t)list->GetN()));
   }
}

//__________________________________________________________________________________________________
void AliRsnDaughterSelector::Init()
{
//
// Initialize the arrays of entry lists to the same size
// of the corresponding arrays of cut sets.
// If they are not empty, they are cleared.
//
   
   Int_t i, nSets;
   
   // neutral
   nSets = fCutSetsN.GetEntries();
   if (!fEntryListsN.IsEmpty()) fEntryListsN.Delete();
   for (i = 0; i < nSets; i++) {
      AliRsnCutSet *set = (AliRsnCutSet*)fCutSetsN[i];
      new (fEntryListsN[i]) TEntryList;
      AliInfo(Form("Adding 1 entry list for neutrals --> cut set '%s' [scheme = '%s']", set->GetName(), set->GetCutScheme().Data()));
   }
   
   // charged
   nSets = fCutSetsC.GetEntries();
   if (!fEntryListsP.IsEmpty()) fEntryListsP.Delete();
   if (!fEntryListsM.IsEmpty()) fEntryListsM.Delete();
   for (i = 0; i < nSets; i++) {
      AliRsnCutSet *set = (AliRsnCutSet*)fCutSetsC[i];
      new (fEntryListsP[i]) TEntryList;
      new (fEntryListsM[i]) TEntryList;
      AliInfo(Form("Adding 2 entry lists for charged --> cut set '%s' [scheme = '%s']", set->GetName(), set->GetCutScheme().Data()));
   }
}

//__________________________________________________________________________________________________
void AliRsnDaughterSelector::Reset()
{
   TEntryList *el;
   Int_t i, nSets;
   
   // N
   nSets = fCutSetsN.GetEntries();
   for (i = 0; i < nSets; i++) {
      el = (TEntryList*)fEntryListsN.At(i);
      el->Reset();
   }
   
   // charged
   nSets = fCutSetsC.GetEntries();
   for (i = 0; i < nSets; i++) {
      el = (TEntryList*)fEntryListsP.At(i);
      el->Reset();
      el = (TEntryList*)fEntryListsM.At(i);
      el->Reset();
   }
}

//__________________________________________________________________________________________________
Int_t AliRsnDaughterSelector::Add(AliRsnCutSet *cuts, Bool_t charged)
{
//
// Add a new selection slot defined by a set of cuts and daughter definition
//

   Int_t n = 0;
   
   if (!charged) {
      n = fCutSetsN.GetEntries();
      new (fCutSetsN[n]) AliRsnCutSet(*cuts);
   } else {
      n = fCutSetsC.GetEntries();
      new (fCutSetsC[n]) AliRsnCutSet(*cuts);
   }
   
   return n;
}

//__________________________________________________________________________________________________
Int_t AliRsnDaughterSelector::GetID(const char *name, Bool_t charged)
{
//
// Add a new selection slot defined by a set of cuts and daughter definition
//

   AliRsnCutSet *cuts;
   
   if (!charged) {
      cuts = (AliRsnCutSet*)fCutSetsN.FindObject(name);
      if (cuts) return fCutSetsN.IndexOf(cuts);
   } else {
      cuts = (AliRsnCutSet*)fCutSetsC.FindObject(name);
      if (cuts) return fCutSetsC.IndexOf(cuts);
   }
   
   return -1;
}

//__________________________________________________________________________________________________
TEntryList* AliRsnDaughterSelector::GetSelected(Int_t i, Char_t charge)
{
//
// Retrieve a given entry list (needs charge specified as a char)
//

   if (charge == '+')
      return (TEntryList*)fEntryListsP.At(i);
   else if (charge == '-')
      return (TEntryList*)fEntryListsM.At(i);
   else
      return (TEntryList*)fEntryListsN.At(i);
}

//__________________________________________________________________________________________________
TEntryList* AliRsnDaughterSelector::GetSelected(Int_t i, Short_t charge)
{
//
// Retrieve a given entry list passing charge as short
//

   if (charge > 0)
      return (TEntryList*)fEntryListsP.At(i);
   else if (charge < 0)
      return (TEntryList*)fEntryListsM.At(i);
   else
      return (TEntryList*)fEntryListsN.At(i);
}

//__________________________________________________________________________________________________
void AliRsnDaughterSelector::ScanEvent(AliRsnEvent* ev)
{
//
// Loop over event and fill all entry lists
//

   Int_t id, is;
   Int_t nSel, nTot = ev->GetAbsoluteSum();
   AliRsnDaughter check;
   TClonesArray *cutsArray = 0x0, *entryArray = 0x0;
   
   for (id = 0; id < nTot; id++) {
      ev->SetDaughterAbs(check, id);
      if (!check.IsOK()) continue;
      // set pointers according to charge
      switch (check.ChargeS()) {
         case 1:
            cutsArray = &fCutSetsC;
            entryArray = &fEntryListsP;
            break;
         case -1:
            cutsArray = &fCutSetsC;
            entryArray = &fEntryListsM;
            break;
         default:
            cutsArray = &fCutSetsN;
            entryArray = &fEntryListsN;
            break;
      }
      // check with all cuts in that charge
      nSel = cutsArray->GetEntriesFast();
      for (is = 0; is < nSel; is++) {
         AliRsnCutSet *cuts = (AliRsnCutSet*)cutsArray->At(is);
         if (cuts->IsSelected(&check)) {
            TEntryList *el = (TEntryList*)entryArray->At(is);
            el->Enter(id);
         } 
      }
   }
   
   //Print();
}
