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
      evMain->SetDaughterAbs(fDaughter[0], (Int_t)list0->GetEntry(i0));
      start = 0;
      if (!fDaughter[0].GetRef()) {
         AliDebugClass(3, Form("[%s]: daughte #1 has NULL ref", GetName()));
         continue;
      }
      if (!fDaughter[0].IsOK()) {
         AliDebugClass(3, Form("[%s]: daughte #1 is BAD", GetName()));
         continue;
      }
      if (!fIsMixed && list0 == list1) start = i0 + 1;
      for (i1 = start; i1 < list1->GetN(); i1++) {
         evMix->SetDaughterAbs(fDaughter[1], (Int_t)list1->GetEntry(i1));
         if (!fDaughter[1].GetRef()) {
            AliDebugClass(3, Form("[%s]: daughte #2 has NULL ref", GetName()));
            continue;
         }
         if (!fDaughter[1].IsOK()) {
            AliDebugClass(3, Form("[%s]: daughte #2 is BAD", GetName()));
            continue;
         }
         // check mother
         if (!MotherOK()) {
            AliDebugClass(2, Form("[%s]: candidate mother didn't pass the cuts", GetName()));
            continue;
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
Bool_t AliRsnLoopPair::MotherOK()
{
//
// Checks that first argument matches definitions for first daughter
// and the same for second argument, where the order is defined by
// the AliRsnPairDef data member.
// If the matching is successful, the AliRsnMother data member is 
// initialized using the mass hypotheses defined here and the momenta
// in the passed daughters.
// The third argument is necessary to choose which one of the possible two
// events owning the two daughter will be used as reference.
//
   
   // check matching and exit if one of them fails
   // if true pair is required, this is taken into account:
   // if both true pairs and correct decay tree is required,
   // then we must be sure that also the true PID of daughters matches,
   // instead if correct decay tree is not required this additional check is not done
   fPairDef->GetDef1().SetOnlyTrue(fOnlyTrue && fCheckDecay);
   fPairDef->GetDef2().SetOnlyTrue(fOnlyTrue && fCheckDecay);
   if (!fPairDef->GetDef1().MatchesDaughter(&fDaughter[0])) return kFALSE;
   if (!fPairDef->GetDef2().MatchesDaughter(&fDaughter[1])) return kFALSE;
   
   // if matching is successful
   // compute 4-momenta of daughters and mother
   fMother.ComputeSum(fPairDef->GetDef1().GetMass(), fPairDef->GetDef2().GetMass());
   
   // if required a true pair, check this here and eventually return a fail message
   // this is done using the method AliRsnMother::CommonMother with 2 arguments
   // passed by reference, where the real GEANT label of the particle is stored
   // and one can check if these tracks are both really secondaries (ID >= 0)
   if (fOnlyTrue) {
      Int_t m0, m1, common;
      common = fMother.CommonMother(m0, m1);
      if (m0 < 0 || m1 < 0) return kFALSE;
      if (common != fPairDef->GetMotherPDG()) return kFALSE;
   }
   
   // point to first event as reference
   // and checks the pair cuts,
   // (done first because it is more likely 
   // that it is not passed and execution is faster)
   if (fPairCuts)
      return fPairCuts->IsSelected(&fMother);
   else
      return kTRUE;
}
