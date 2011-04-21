//
// *** Class AliRsnLoopPair ***
//
// "Core" method for defining the work on a pari of particles.
// For one analysis, one must setup one of this for each pair he wants to analyze,
// adding to it all analysis which he desires to do.
// Here he defines the cuts, and the particle types and charges, and can add
// functions which do different operations on the same pair, and some binning
// with respect to some kinematic variables (eta, momentum)
//
// authors: A. Pulvirenti (email: alberto.pulvirenti@ct.infn.it)
//          M. Vala (email: martin.vala@cern.ch)
//

#include <Riostream.h>
#include <TList.h>
#include <TEntryList.h>

#include "AliLog.h"

#include "AliRsnMother.h"
#include "AliRsnCutSet.h"
#include "AliRsnDaughterSelector.h"

#include "AliRsnLoopPair.h"

ClassImp(AliRsnLoopPair)

//_____________________________________________________________________________
AliRsnLoopPair::AliRsnLoopPair(const char *name, AliRsnPairDef *def, Bool_t isMixed) :
   AliRsnLoop(name, isMixed),
   fOnlyTrue(kFALSE),
   fCheckDecay(kFALSE),
   fPairDef(def),
   fPairCuts(0x0),
   fMother()
{
//
// Default constructor
//

   fListID[0] = -1;
   fListID[1] = -1;
}

//_____________________________________________________________________________
AliRsnLoopPair::AliRsnLoopPair(const AliRsnLoopPair& copy) :
   AliRsnLoop(copy),
   fOnlyTrue(copy.fOnlyTrue),
   fCheckDecay(copy.fCheckDecay),
   fPairDef(copy.fPairDef),
   fPairCuts(copy.fPairCuts),
   fMother(copy.fMother)
{
//
// Copy constructor
//

   fListID[0] = copy.fListID[0];
   fListID[1] = copy.fListID[1];
}

//_____________________________________________________________________________
AliRsnLoopPair& AliRsnLoopPair::operator=(const AliRsnLoopPair& copy)
{
   AliRsnLoop::operator=(copy);
   
   fOnlyTrue = copy.fOnlyTrue;
   fCheckDecay = copy.fCheckDecay;
   fPairDef = copy.fPairDef;
   fPairCuts = copy.fPairCuts;
   fMother = copy.fMother;
   fListID[0] = copy.fListID[0];
   fListID[1] = copy.fListID[1];

   return (*this);
}

//_____________________________________________________________________________
AliRsnLoopPair::~AliRsnLoopPair()
{
//
// Destructor
//
}

//_____________________________________________________________________________
void AliRsnLoopPair::Print(Option_t* /*option*/) const
{
//
// Prints info about pair
//

   AliRsnLoop::Print();
}

//_____________________________________________________________________________
Bool_t AliRsnLoopPair::Init(const char *prefix, TList *list)
{
//
// Initialization function.
// Loops on all functions and eventual the ntuple, to initialize output objects.
//

   // assign data members relations
   fMother.SetDaughter(0, &fDaughter[0]);
   fMother.SetDaughter(1, &fDaughter[1]);
   AliInfo(Form("[%s] Initialization", GetName()));
   
   TString name(prefix);
   name += '_';
   name += GetName();
   if (IsMixed()) name.Prepend("mix_");
   
   return AliRsnLoop::Init(name.Data(), list);
}

//_____________________________________________________________________________
Int_t AliRsnLoopPair::DoLoop
(AliRsnEvent *evMain, AliRsnDaughterSelector *selMain, AliRsnEvent *evMix, AliRsnDaughterSelector *selMix)
{
//
// Loop function.
// Computes what is needed from passed events.
// Returns the number of pairs successfully processed.
//

   if (fIsMixed) {
      AliDebugClass(1, Form("[%s]: event-mixing loop", GetName()));
      if (!evMix || !selMix) {
         AliError(Form("[%s] NULL mixed event when mixing is required: cannot process", GetName()));
         return 0;
      }
   } else {
      AliDebugClass(1, Form("[%s]: single-event loop", GetName()));
      evMix = evMain;
      selMix = selMain;
   }
   fMother.SetRefEvent(evMain);
   
   // check cuts
   if (!OkEvent(evMain)) {
      AliDebugClass(1, Form("[%s]: main event not accepted", GetName()));
      return 0;
   }
   if (!OkEvent(evMix)) {
      AliDebugClass(1, Form("[%s]: mixed event not accepted", GetName()));
      return 0;
   }
   
   Int_t i0, i1, start, npairs = 0;
   
   TEntryList *list0 = selMain->GetSelected(fListID[0], fPairDef->GetDef1().GetChargeC());
   TEntryList *list1 = selMix ->GetSelected(fListID[1], fPairDef->GetDef2().GetChargeC());
   if (!list0 || !list1) {
      AliError("Can't process NULL lists");
      return 0;
   }
   AliDebugClass(2, Form("[%s]: list counts: %d, %d", GetName(), (Int_t)list0->GetN(), (Int_t)list1->GetN()));
   if (!list0->GetN() || !list1->GetN()) {
      AliDebugClass(2, Form("[%s]: at least one list is empty", GetName()));
      return 0;
   }
   
   TObjArrayIter next(&fOutputs);
   AliRsnListOutput *out = 0x0;
   
   for (i0 = 0; i0 < list0->GetN(); i0++) {
      evMain->SetDaughter(fDaughter[0], (Int_t)list0->GetEntry(i0));
      fDaughter[0].FillP(fPairDef->GetDef1().GetMass());
      start = 0;
      if (!fIsMixed && list0 == list1) start = i0 + 1;
      for (i1 = start; i1 < list1->GetN(); i1++) {
         AliDebugClass(4, Form("Checking entries pair: %d (%d) with %d (%d)", (Int_t)i0, (Int_t)list0->GetEntry(i0), (Int_t)i1, (Int_t)list1->GetEntry(i1)));
         evMix->SetDaughter(fDaughter[1], (Int_t)list1->GetEntry(i1));
         fDaughter[1].FillP(fPairDef->GetDef2().GetMass());
         fMother.Sum(0) = fDaughter[0].Prec() + fDaughter[1].Prec();
         fMother.Sum(1) = fDaughter[0].Psim() + fDaughter[1].Psim();
         fMother.Ref(0).SetXYZM(fMother.Sum(0).X(), fMother.Sum(0).Y(), fMother.Sum(0).Z(), fPairDef->GetMotherMass());
         fMother.Ref(1).SetXYZM(fMother.Sum(1).X(), fMother.Sum(1).Y(), fMother.Sum(1).Z(), fPairDef->GetMotherMass());
         // check cuts
         if (fPairCuts) {
            if (!fPairCuts->IsSelected(&fMother)) {
               AliDebugClass(2, Form("[%s]: candidate mother didn't pass the cuts", GetName()));
               continue;
            }
         }
         // check mother
         if (fOnlyTrue) {
            if (!IsTrueMother()) {
               AliDebugClass(2, Form("[%s]: candidate mother is not true", GetName()));
               continue;
            }
         }
         // fill outputs
         next.Reset();
         while ( (out = (AliRsnListOutput*)next()) ) {
            if (out->Fill(&fMother)) npairs++;
            else AliDebugClass(3, Form("[%s]: failed computation", GetName()));
         }
      }
   }

   return npairs;
}

//_____________________________________________________________________________
Bool_t AliRsnLoopPair::IsTrueMother()
{
//
// Checks to see if the mother comes from a true resonance.
// It is triggered by the 'SetOnlyTrue()' function
//
      
   // check #1:
   // daughters have same mother with the right PDG code
   Int_t commonPDG = fMother.CommonMother();
   if (commonPDG != fPairDef->GetMotherPDG()) return kFALSE;
   AliDebugClass(4, "Found a true mother");
   
   // check #2:
   // checks if daughter have the right particle type
   // (activated by fCheckDecay)
   if (fCheckDecay) {
      AliRsnDaughterDef &def1 = fPairDef->GetDef1();
      AliRsnDaughterDef &def2 = fPairDef->GetDef2();
      if (!def1.MatchesPID(&fDaughter[0])) return kFALSE;
      if (!def2.MatchesPID(&fDaughter[1])) return kFALSE;
   }
   AliDebugClass(4, "Decay products match");
   
   return kTRUE;
}
