//
// Class AliRsnLoopEffPair
//
// Inherits from basic AliRsnLoopEff for efficiency,
// and computed efficiencies for single-tracks
//
// author: Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include <Riostream.h>

#include "AliStack.h"
#include "AliGenEventHeader.h"
#include "AliAODMCHeader.h"

#include "AliRsnEvent.h"
#include "AliRsnCutManager.h"
#include "AliRsnDaughterDef.h"
#include "AliRsnPairDef.h"
#include "AliRsnLoopEffPair.h"

ClassImp(AliRsnLoopEffPair)

//_____________________________________________________________________________
AliRsnLoopEffPair::AliRsnLoopEffPair(const char *name, AliRsnPairDef *def) :
   AliRsnLoopEff(name, 1),
   fDef(def),
   fMother()
{
//
// Default constructor.
// Do not repeat 'DefineOutput' since it is done in base class and we don't add new ones.
//
}

//_____________________________________________________________________________
AliRsnLoopEffPair::AliRsnLoopEffPair(const AliRsnLoopEffPair& copy) :
   AliRsnLoopEff(copy),
   fDef(copy.fDef),
   fMother(copy.fMother)
{
//
// Copy constrtuctor.
//
}

//_____________________________________________________________________________
AliRsnLoopEffPair& AliRsnLoopEffPair::operator=(const AliRsnLoopEffPair& copy)
{
//
// Assignment operator.
// Owned data members are meaningless for this operator.
//

   AliRsnLoopEff::operator=(copy);
   fDef = copy.fDef;

   return (*this);
}

//__________________________________________________________________________________________________
Bool_t AliRsnLoopEffPair::AssignMotherAndDaughters(AliRsnEvent *rsnEvent, Int_t ipart)
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
Bool_t AliRsnLoopEffPair::AssignMotherAndDaughtersESD(AliRsnEvent *rsnEvent, Int_t ipart)
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
   AliMCParticle *mother  = (AliMCParticle*)mc->GetTrack(ipart);
   TParticle     *motherP = mother->Particle();
   Int_t          ntracks = stack->GetNtrack();
   
   // check PDG code and exit if it is wrong
   if (TMath::Abs(motherP->GetPdgCode()) != fDef->GetMotherPDG()) return kFALSE;
   
   // check number of daughters and exit if it is not 2
   if (motherP->GetNDaughters() < 2) return kFALSE;
   
   // check distance from primary vertex
   TLorentzVector vprod;
   motherP->ProductionVertex(vprod);
   if (!CheckDistanceFromPV(vprod.X(), vprod.Y(), vprod.Z())) {
      AliDebugClass(1, "Distant production vertex");
      return kFALSE;
   }
   
   // get the daughters and check their PDG code and charge:
   // if they match one of the pair daughter definitions, 
   // assign them as MC reference of the 'fDaughter' objects
   fDaughter[0].Reset();
   fDaughter[1].Reset();
   Int_t index[2] = {motherP->GetFirstDaughter(), motherP->GetLastDaughter()};
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
      daughter = (AliMCParticle*)mc->GetTrack(index[i]);
      pdg      = TMath::Abs(daughter->Particle()->GetPdgCode());
      charge   = (Short_t)(daughter->Particle()->GetPDG()->Charge() / 3);
      // check if it matches one definition
      if (fDef->GetDef1().MatchesPDG(pdg) && fDef->GetDef1().MatchesChargeS(charge)) {
         fDaughter[0].SetGood();
         fDaughter[0].SetRefMC(daughter);
         fDaughter[0].SetLabel(index[i]);
      } else if (fDef->GetDef2().MatchesPDG(pdg) && fDef->GetDef2().MatchesChargeS(charge)) {
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
Bool_t AliRsnLoopEffPair::AssignMotherAndDaughtersAOD(AliRsnEvent *rsnEvent, Int_t ipart)
{
   AliAODEvent      *aod     = rsnEvent->GetRefAOD();
   TClonesArray     *listAOD = (TClonesArray*)(aod->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
   AliAODMCParticle *mother  = (AliAODMCParticle*)listAOD->At(ipart);
   Int_t             ntracks = listAOD->GetEntries();
   
   // check PDG code and exit if it is wrong
   if (TMath::Abs(mother->GetPdgCode()) != fDef->GetMotherPDG()) return kFALSE;
   
   // check number of daughters and exit if it is not 2
   if (mother->GetNDaughters() < 2) return kFALSE;
   
   // check distance from primary vertex
   Double_t vprod[3] = {(Double_t)mother->Xv(), (Double_t)mother->Yv(), (Double_t)mother->Zv()};
   Double_t dv = DistanceFromPV(vprod[0], vprod[1], vprod[2]);
   if (dv > fMaxDistPV) {
      AliDebugClass(1, "Distant production vertex");
      return kFALSE;
   }
   
   // get the daughters and check their PDG code and charge:
   // if they match one of the pair daughter definitions, 
   // assign them as MC reference of the 'fDaughter' objects
   fDaughter[0].Reset();
   fDaughter[1].Reset();
   Int_t index[2] = {(Int_t)mother->GetDaughter(0), (Int_t)mother->GetDaughter(1)};
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
      daughter = (AliAODMCParticle*)listAOD->At(index[i]);
      pdg      = TMath::Abs(daughter->GetPdgCode());
      charge   = (Short_t)(daughter->Charge() / 3);
      // check if it matches one definition
      if (fDef->GetDef1().MatchesPDG(pdg) && fDef->GetDef1().MatchesChargeS(charge)) {
         fDaughter[0].SetGood();
         fDaughter[0].SetRefMC(daughter);
         fDaughter[0].SetLabel(index[i]);
      } else if (fDef->GetDef2().MatchesPDG(pdg) && fDef->GetDef2().MatchesChargeS(charge)) {
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
Int_t AliRsnLoopEffPair::DoLoop(AliRsnEvent *rsn, AliRsnDaughterSelector*, AliRsnEvent*, AliRsnDaughterSelector*)
{
//
// Loop on event and fill containers
//

   // check event cuts
   if (!OkEvent(rsn)) return 0;
   
   // retrieve output
   fOutput = (AliRsnListOutput*)fOutputs[0];
   
   // check presence of MC reference
   if (!rsn->GetRefMC()) {
      AliError("Need a MC to compute efficiency");
      return 0;
   }
   
   // check presence of event
   if (!rsn->GetRef()) {
      AliError("Need an event to compute efficiency");
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
   Int_t i, npart = 0;
   if (rsn->IsESD()) {
      TArrayF mcVertex(3);
      AliGenEventHeader *genEH = rsn->GetRefMCESD()->GenEventHeader();
      genEH->PrimaryVertex(mcVertex);
      for (i = 0; i < 3; i++) fVertex[i] = (Double_t)mcVertex[i];
      npart = rsn->GetRefMCESD()->GetNumberOfTracks();
   } else {
      for (i = 0; i < 3; i++) fVertex[i] = 0.0;
      AliAODEvent *aod = rsn->GetRefMCAOD();
      TClonesArray *listAOD = (TClonesArray*)(aod->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
      if (listAOD) npart = listAOD->GetEntries();
      AliAODMCHeader *mcH = static_cast<AliAODMCHeader*>(aod->FindListObject(AliAODMCHeader::StdBranchName()));
      if (mcH) mcH->GetVertex(fVertex);
   }
   
   // check number of particles
   if (!npart) {
      AliInfo("Empty event");
      return 0;
   }
   
   // utility variables
   Int_t ipart, istep, count = 0, nsteps = fSteps.GetEntries();
   Int_t ntracks = rsn->GetAbsoluteSum();
   AliRsnDaughter check;
   
   // loop over particles
   for (ipart = 0; ipart < npart; ipart++) {
      // check i-th particle
      if (!AssignMotherAndDaughters(rsn, ipart)) continue;
      // if assignment was successful, for first step, use MC info
      fDaughter[0].SetRef(fDaughter[0].GetRefMC());
      fDaughter[1].SetRef(fDaughter[1].GetRefMC());
      fMother.SetRefEvent(rsn);
      fMother.ComputeSum(fDef->GetDef1().GetMass(), fDef->GetDef2().GetMass(), fDef->GetMotherMass());
      fOutput->Fill(&fMother, 0);
      count++;
      // for each further step, try to find two tracks which pass the related cuts
      for (istep = 0; istep < nsteps; istep++) {
         AliRsnCutManager *cuts = (AliRsnCutManager*)fSteps[istep];
         fDaughter[0].SetBad();
         fDaughter[1].SetBad();
         for (i = 0; i < ntracks; i++) {
            rsn->SetDaughter(check, i);
            if (!cuts->PassCommonDaughterCuts(&check)) continue;
            if (TMath::Abs(check.GetLabel()) == fDaughter[0].GetLabel()) {
               if (!cuts->PassDaughter1Cuts(&check)) continue;
               fDaughter[0].SetRef(check.GetRef());
               fDaughter[0].SetGood();
            } else if (TMath::Abs(check.GetLabel()) == fDaughter[1].GetLabel()) {
               if (!cuts->PassDaughter2Cuts(&check)) continue;
               fDaughter[1].SetRef(check.GetRef());
               fDaughter[1].SetGood();
            }
            if (fDaughter[0].IsOK() && fDaughter[1].IsOK()) {
               fMother.ComputeSum(fDef->GetDef1().GetMass(), fDef->GetDef2().GetMass(), fDef->GetMotherMass());
               if (cuts->PassMotherCuts(&fMother)) {
                  fOutput->Fill(&fMother, 1 + istep);
                  break;
               }
            }
         }
      }
   }
   
   return count;
}
