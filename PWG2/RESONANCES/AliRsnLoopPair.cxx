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

#include "AliRsnEvent.h"
#include "AliRsnPairDef.h"
#include "AliRsnMother.h"
#include "AliRsnCutSet.h"
#include "AliRsnDaughterSelector.h"

#include "AliRsnLoopPair.h"

ClassImp(AliRsnLoopPair)

//_____________________________________________________________________________
AliRsnLoopPair::AliRsnLoopPair(const char *name, AliRsnPairDef *def, Bool_t isMixed) :
   AliRsnLoop(name, isMixed),
   fTrueMC(kFALSE),
   fOnlyTrue(kFALSE),
   fUseMCRef(kFALSE),
   fCheckDecay(kFALSE),
   fRangeY(1E20),
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
AliRsnLoopPair::AliRsnLoopPair(const AliRsnLoopPair &copy) :
   AliRsnLoop(copy),
   fTrueMC(copy.fTrueMC),
   fOnlyTrue(copy.fOnlyTrue),
   fUseMCRef(copy.fUseMCRef),
   fCheckDecay(copy.fCheckDecay),
   fRangeY(copy.fRangeY),
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
AliRsnLoopPair &AliRsnLoopPair::operator=(const AliRsnLoopPair &copy)
{
//
// Assignment operator
//

   AliRsnLoop::operator=(copy);

   fTrueMC = copy.fTrueMC;
   fOnlyTrue = copy.fOnlyTrue;
   fUseMCRef = copy.fUseMCRef;
   fCheckDecay = copy.fCheckDecay;
   fRangeY = copy.fRangeY;
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
//    if (IsMixed()) name.Prepend("mix_");
   if (IsMixed()) name.Append("_mix");

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
      AliDebugClass(3, Form("[%s]: event-mixing loop", GetName()));
      if (!evMix || !selMix) {
         AliError(Form("[%s] NULL mixed event when mixing is required: cannot process", GetName()));
         return 0;
      }
   } else {
      AliDebugClass(3, Form("[%s]: single-event loop", GetName()));
      evMix = evMain;
      selMix = selMain;
   }
   fMother.SetRefEvent(evMain);

   // check cuts
   if (!OkEvent(evMain)) {
      AliDebugClass(3, Form("[%s]: main event not accepted", GetName()));
      return 0;
   }
   if (!OkEvent(evMix)) {
      AliDebugClass(3, Form("[%s]: mixed event not accepted", GetName()));
      return 0;
   }

   // if it is required to loop over True MC, do this here and skip the rest of the method
   if (fTrueMC) return LoopTrueMC(evMain);

   Int_t i0, i1, start, npairs = 0;

   TEntryList *list0 = selMain->GetSelected(fListID[0], fPairDef->GetDef1().GetChargeC());
   TEntryList *list1 = selMix ->GetSelected(fListID[1], fPairDef->GetDef2().GetChargeC());
   if (!list0 || !list1) {
      AliError("Can't process NULL lists");
      return 0;
   }
   AliDebugClass(3, Form("[%s]: list counts: %lld, %lld", GetName(), list0->GetN(), list1->GetN()));
   if (!list0->GetN() || !list1->GetN()) {
      AliDebugClass(3, Form("[%s]: at least one list is empty", GetName()));
      return 0;
   }

   TObjArrayIter next(&fOutputs);
   AliRsnListOutput *out = 0x0;
   Long64_t iEntry1,iEntry2;
   for (i0 = 0; i0 < list0->GetN(); i0++) {
      iEntry1 = list0->GetEntry(i0);
      evMain->SetDaughter(fDaughter[0], iEntry1,fUseMCRef);
      fDaughter[0].FillP(fPairDef->GetDef1().GetMass());
      start = 0;
      if (!fIsMixed && list0 == list1) start = i0 + 1;
      for (i1 = start; i1 < list1->GetN(); i1++) {
         iEntry2 = list1->GetEntry(i1);
	 if (iEntry1 == iEntry2) continue;
         AliDebugClass(4, Form("Checking entries pair: %d (%lld) with %d (%lld)", i0, iEntry1, i1, iEntry2));
         evMix->SetDaughter(fDaughter[1], iEntry2,fUseMCRef);
         fDaughter[1].FillP(fPairDef->GetDef2().GetMass());
         fMother.Sum(0) = fDaughter[0].Prec() + fDaughter[1].Prec();
         fMother.Sum(1) = fDaughter[0].Psim() + fDaughter[1].Psim();
         fMother.Ref(0).SetXYZM(fMother.Sum(0).X(), fMother.Sum(0).Y(), fMother.Sum(0).Z(), fPairDef->GetMotherMass());
         fMother.Ref(1).SetXYZM(fMother.Sum(1).X(), fMother.Sum(1).Y(), fMother.Sum(1).Z(), fPairDef->GetMotherMass());
         // check rapidity range
         if (TMath::Abs(fMother.Rapidity(0)) > fRangeY) {
            AliDebugClass(2, Form("[%s]: Outside rapidity range", GetName()));
            continue;
         }
         // check mother
         if (fOnlyTrue) {
            if (!IsTrueMother()) {
               AliDebugClass(2, Form("[%s]: candidate mother is not true", GetName()));
               continue;
            }
         }
         // check cuts
         if (fPairCuts) {
            if (!fPairCuts->IsSelected(&fMother)) {
               AliDebugClass(2, Form("[%s]: candidate mother didn't pass the cuts", GetName()));
               continue;
            }
         }
         // fill outputs
         next.Reset();
         while ( (out = (AliRsnListOutput *)next()) ) {
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
   AliDebugClass(1, "Found a true mother");

   // check #2:
   // checks if daughter have the right particle type
   // (activated by fCheckDecay)
   if (fCheckDecay) {
      AliRsnDaughterDef &def1 = fPairDef->GetDef1();
      AliRsnDaughterDef &def2 = fPairDef->GetDef2();
      if (!def1.MatchesPID(&fDaughter[0])) return kFALSE;
      if (!def2.MatchesPID(&fDaughter[1])) return kFALSE;
   }
   AliDebugClass(1, "Decay products match");

   return kTRUE;
}

//__________________________________________________________________________________________________
Bool_t AliRsnLoopPair::AssignMotherAndDaughters(AliRsnEvent *rsnEvent, Int_t ipart)
{
//
// Calls the appropriate assignment method
//

   // setup pointers
   fMother.SetDaughter(0, &fDaughter[0]);
   fMother.SetDaughter(1, &fDaughter[1]);
   fMother.SetRefEvent(rsnEvent);
   fDaughter[0].SetOwnerEvent(rsnEvent);
   fDaughter[1].SetOwnerEvent(rsnEvent);

   if (rsnEvent->IsESD())
      return AssignMotherAndDaughtersESD(rsnEvent, ipart);
   else if (rsnEvent->IsAOD())
      return AssignMotherAndDaughtersAOD(rsnEvent, ipart);
   else {
      AliError("Unrecognized input event");
      return kFALSE;
   }
}

//__________________________________________________________________________________________________
Bool_t AliRsnLoopPair::AssignMotherAndDaughtersESD(AliRsnEvent *rsnEvent, Int_t ipart)
{
//
// Gets a particle in the MC event and try to assign it to the mother.
// If it has two daughters, retrieve them and assign also them.
// NOTE: assignment is done only for MC, since reconstructed match is assigned in the same way
//       for ESD and AOD, if available
// ---
// Implementation for ESD inputs
//

   AliMCEvent    *mc      = rsnEvent->GetRefMCESD();
   AliStack      *stack   = mc->Stack();
   AliMCParticle *mother  = (AliMCParticle *)mc->GetTrack(ipart);
   TParticle     *motherP = mother->Particle();
   Int_t          ntracks = stack->GetNtrack();

   // check PDG code and exit if it is wrong
   if (TMath::Abs(motherP->GetPdgCode()) != fPairDef->GetMotherPDG()) return kFALSE;

   // check number of daughters and exit if it is not 2
   if (motherP->GetNDaughters() < 2) return kFALSE;
   //if (!stack->IsPhysicalPrimary(ipart)) return kFALSE;

   /*
   // check distance from primary vertex
   TLorentzVector vprod;
   motherP->ProductionVertex(vprod);
   if (!CheckDistanceFromPV(vprod.X(), vprod.Y(), vprod.Z())) {
      AliDebugClass(1, "Distant production vertex");
      return kFALSE;
   }
   */

   // get the daughters and check their PDG code and charge:
   // if they match one of the pair daughter definitions,
   // assign them as MC reference of the 'fDaughter' objects
   fDaughter[0].Reset();
   fDaughter[1].Reset();
   Int_t index[2] = {motherP->GetDaughter(0), motherP->GetDaughter(1)};
   Int_t i, pdg;
   Short_t charge;
   AliMCParticle *daughter = 0x0;
   for (i = 0; i < 2; i++) {
      // check index for stack
      if (index[i] < 0 || index[i] > ntracks) {
         AliError(Form("Index %d overflow: value = %d, max = %d", i, index[i], ntracks));
         return kFALSE;
      }
      // get daughter and its PDG and charge
      daughter = (AliMCParticle *)mc->GetTrack(index[i]);
      pdg      = TMath::Abs(daughter->Particle()->GetPdgCode());
      charge   = (Short_t)(daughter->Particle()->GetPDG()->Charge() / 3);
      // check if it matches one definition
      if (fPairDef->GetDef1().MatchesPDG(pdg) && fPairDef->GetDef1().MatchesChargeS(charge)) {
         fDaughter[0].SetGood();
         fDaughter[0].SetRefMC(daughter);
         fDaughter[0].SetLabel(index[i]);
      } else if (fPairDef->GetDef2().MatchesPDG(pdg) && fPairDef->GetDef2().MatchesChargeS(charge)) {
         fDaughter[1].SetGood();
         fDaughter[1].SetRefMC(daughter);
         fDaughter[1].SetLabel(index[i]);
      }
   }

   // return success if both daughters were assigned
   if (fDaughter[0].IsOK() && fDaughter[1].IsOK()) {
      return kTRUE;
   } else {
      fDaughter[0].Reset();
      fDaughter[1].Reset();
      return kFALSE;
   }
}

//__________________________________________________________________________________________________
Bool_t AliRsnLoopPair::AssignMotherAndDaughtersAOD(AliRsnEvent *rsnEvent, Int_t ipart)
{
//
// Gets a particle in the MC event and try to assign it to the mother.
// If it has two daughters, retrieve them and assign also them.
// NOTE: assignment is done only for MC, since reconstructed match is assigned in the same way
//       for ESD and AOD, if available
// ---
// Implementation for AOD inputs
//

   AliAODEvent      *aod     = rsnEvent->GetRefAOD();
   TClonesArray     *listAOD = (TClonesArray *)(aod->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
   AliAODMCParticle *mother  = (AliAODMCParticle *)listAOD->At(ipart);

   // check PDG code and exit if it is wrong
   if (TMath::Abs(mother->GetPdgCode()) != fPairDef->GetMotherPDG()) return kFALSE;

   // check number of daughters and exit if it is not 2
   if (mother->GetNDaughters() < 2) return kFALSE;
   if (!mother->IsPrimary()) return kFALSE;

   /*
   // check distance from primary vertex
   Double_t vprod[3] = {(Double_t)mother->Xv(), (Double_t)mother->Yv(), (Double_t)mother->Zv()};
   Double_t dv = DistanceFromPV(vprod[0], vprod[1], vprod[2]);
   if (dv > fMaxDistPV) {
      AliDebugClass(1, "Distant production vertex");
      return kFALSE;
   }
   */

   // get the daughters and check their PDG code and charge:
   // if they match one of the pair daughter definitions,
   // assign them as MC reference of the 'fDaughter' objects
   fDaughter[0].Reset();
   fDaughter[1].Reset();
   Int_t index[2] = {(Int_t)mother->GetDaughter(0), (Int_t)mother->GetDaughter(1)};
   Int_t ntracks = listAOD->GetEntriesFast();
   Int_t i, pdg;
   Short_t charge;
   AliAODMCParticle *daughter = 0x0;
   for (i = 0; i < 2; i++) {
      // check index for stack
      if (index[i] < 0 || index[i] > ntracks) {
         AliError(Form("Index %d overflow: value = %d, max = %d", i, index[i], ntracks));
         return kFALSE;
      }
      // get daughter and its PDG and charge
      daughter = (AliAODMCParticle *)listAOD->At(index[i]);
      pdg      = TMath::Abs(daughter->GetPdgCode());
      charge   = (Short_t)(daughter->Charge() / 3);
      // check if it matches one definition
      if (fPairDef->GetDef1().MatchesPDG(pdg) && fPairDef->GetDef1().MatchesChargeS(charge)) {
         fDaughter[0].SetGood();
         fDaughter[0].SetRefMC(daughter);
         fDaughter[0].SetLabel(index[i]);
      } else if (fPairDef->GetDef2().MatchesPDG(pdg) && fPairDef->GetDef2().MatchesChargeS(charge)) {
         fDaughter[1].SetGood();
         fDaughter[1].SetRefMC(daughter);
         fDaughter[1].SetLabel(index[i]);
      }
   }

   // return success if both daughters were assigned
   if (fDaughter[0].IsOK() && fDaughter[1].IsOK()) {
      return kTRUE;
   } else {
      fDaughter[0].Reset();
      fDaughter[1].Reset();
      return kFALSE;
   }
}

//_____________________________________________________________________________
Int_t AliRsnLoopPair::LoopTrueMC(AliRsnEvent *rsn)
{
//
// Loop on event and fill containers
//

   // check presence of MC reference
   if (!rsn->GetRefMC()) {
      AliError("Need a MC to compute efficiency");
      return 0;
   }

   // check event type:
   // must be ESD or AOD, and then use a bool to know in the rest
   if (!rsn->IsESD() && !rsn->IsAOD()) {
      AliError("Need to process ESD or AOD input");
      return 0;
   }

   // retrieve the MC primary vertex position
   // and do some additional coherence checks
   //Double_t fVertex[3] = {0.0, 0.0, 0.0};
   Int_t npart = 0;
   if (rsn->IsESD()) {
      //TArrayF mcVertex(3);
      //AliGenEventHeader *genEH = rsn->GetRefMCESD()->GenEventHeader();
      //genEH->PrimaryVertex(mcVertex);
      //for (i = 0; i < 3; i++) fVertex[i] = (Double_t)mcVertex[i];
      npart = rsn->GetRefMCESD()->GetNumberOfTracks();
   } else {
      //for (i = 0; i < 3; i++) fVertex[i] = 0.0;
      AliAODEvent *aod = rsn->GetRefMCAOD();
      TClonesArray *listAOD = (TClonesArray *)(aod->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
      if (listAOD) npart = listAOD->GetEntries();
      //AliAODMCHeader *mcH = static_cast<AliAODMCHeader*>(aod->FindListObject(AliAODMCHeader::StdBranchName()));
      //if (mcH) mcH->GetVertex(fVertex);
   }

   // check number of particles
   if (!npart) {
      AliInfo("Empty event");
      return 0;
   }

   // utility variables
   Int_t ipart, count = 0;
   AliRsnDaughter check;

   TObjArrayIter next(&fOutputs);
   AliRsnListOutput *out = 0x0;

   // loop over particles
   for (ipart = 0; ipart < npart; ipart++) {
      // check i-th particle
      if (!AssignMotherAndDaughters(rsn, ipart)) continue;
      // if assignment was successful, for first step, use MC info
      fDaughter[0].SetRef(fDaughter[0].GetRefMC());
      fDaughter[1].SetRef(fDaughter[1].GetRefMC());
      fMother.SetRefEvent(rsn);
      fMother.ComputeSum(fPairDef->GetDef1().GetMass(), fPairDef->GetDef2().GetMass(), fPairDef->GetMotherMass());
      // check rapidity range
      if (TMath::Abs(fMother.Rapidity(0)) > fRangeY) {
         AliDebugClass(2, Form("[%s]: Outside rapidity range", GetName()));
         continue;
      }
      // fill outputs
      next.Reset();
      while ( (out = (AliRsnListOutput *)next()) ) {
         if (out->Fill(&fMother)) count++;
         else AliDebugClass(3, Form("[%s]: failed computation", GetName()));
      }
   }

   return count;
}

