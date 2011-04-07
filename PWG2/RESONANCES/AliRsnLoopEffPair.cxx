//
// Class AliRsnLoopEffPair
//
// Inherits from basic AliRsnLoopEff for efficiency,
// and computed efficiencies for single-tracks
//
// author: Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include "AliStack.h"

#include "AliRsnEvent.h"
#include "AliRsnCutManager.h"
#include "AliRsnDaughterDef.h"
#include "AliRsnPairDef.h"
#include "AliRsnLoopEffPair.h"

ClassImp(AliRsnLoopEffPair)

//_____________________________________________________________________________
AliRsnLoopEffPair::AliRsnLoopEffPair(const char *name, AliRsnPairDef *def) :
   AliRsnLoopEff(name, 2),
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

//_____________________________________________________________________________
Int_t AliRsnLoopEffPair::DoLoop(AliRsnEvent *rsn, AliRsnDaughterSelector*, AliRsnEvent*, AliRsnDaughterSelector*)
{
//
// Loop on event and fill containers
//

   // check event cuts
   if (!OkEvent(rsn)) return 0;
   
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
   
   // additional coherence checks
   AliVEvent    *mcEvent = 0x0;
   AliStack     *listESD = 0x0;
   TClonesArray *listAOD = 0x0;
   Int_t         npart   = 0;
   if (rsn->IsESD()) {
      mcEvent = rsn->GetRefMCESD();
      listESD = ((AliMCEvent*)mcEvent)->Stack();
      if (!listESD) {
         AliError("Stack is not present");
         return 0;
      }
      npart = listESD->GetNprimary();
   } else {
      mcEvent = rsn->GetRefMCAOD();
      listAOD = (TClonesArray*)((AliAODEvent*)mcEvent)->GetList()->FindObject(AliAODMCParticle::StdBranchName());
      if (!listAOD) {
         AliError("Stack is not present");
         return 0;
      }
      npart = listAOD->GetEntries();
   }
   
   // check number of particles
   if (!npart) {
      AliInfo("Empty event");
      return 0;
   }

   // by default, assign daughters to mother
   // in the correct order (ref. pairDef)
   fMother.SetDaughter(0, &fDaughter[0]);
   fMother.SetDaughter(1, &fDaughter[1]);
   fMother.SetRefEvent(rsn);
   fDaughter[0].SetOwnerEvent(rsn);
   fDaughter[1].SetOwnerEvent(rsn);
   
   // utility variables
   Bool_t            stop;
   Int_t             ipart, pdg, ndaughters, dindex[2], id, label, charge, imatch, index, istep, count = 0, nsteps = NStepsArray();
   AliVParticle     *mother = 0x0, *daughter = 0x0;
   AliMCParticle    *motherESD = 0x0, *daughterESD = 0x0;
   AliAODMCParticle *motherAOD = 0x0, *daughterAOD = 0x0;
   
   // loop over particles
   for (ipart = 0; ipart < npart; ipart++) {
      
      // get next particle and take some quantities
      // in different way from ESD and AOD MC
      if (rsn->IsESD()) {
         mother = mcEvent->GetTrack(ipart);
         motherESD = static_cast<AliMCParticle*>(mother);
         pdg = TMath::Abs(motherESD->Particle()->GetPdgCode());
         ndaughters = motherESD->Particle()->GetNDaughters();
      } else {
         mother = (AliVParticle*)listAOD->At(ipart);
         motherAOD = static_cast<AliAODMCParticle*>(mother);
         pdg = TMath::Abs(motherAOD->GetPdgCode());
         ndaughters = motherAOD->GetNDaughters();
      }
      
      // skip particles with wrong PDG code,
      if (pdg != fDef->GetMotherPDG()) continue;
      
      // fill first step
      fMother.SumMC().SetXYZM(mother->Px(), mother->Py(), mother->Pz(), fDef->GetMotherMass());
      GetOutput()->Fill(&fMother, 0);
      count++;
      
      // reject particles with more/less than 2 daughters
      if (ndaughters != 2) continue;
      
      // retrieve daughters
      // if one of them matches the definition for daughte #1 or #2 in pairDef,
      // assign it to the corresponding AliRsnDaughter member, and set it to 'good'
      // then, if both daughters are not 'good', skip the pair
      fDaughter[0].Reset();
      fDaughter[1].Reset();
      if (rsn->IsESD()) {
         dindex[0] = motherESD->Particle()->GetFirstDaughter();
         dindex[1] = motherESD->Particle()->GetLastDaughter();
      } else {
         dindex[0] = motherAOD->GetDaughter(0);
         dindex[1] = motherAOD->GetDaughter(1);
      }
      
      // try to assign a daughter to each fDaughter[] data member
      stop = kFALSE;
      for (id = 0; id < 2; id++) {
         // avoid overflows
         if (dindex[id] < 0 || dindex[id] > npart) {
            AliWarning(Form("Found a stack overflow in dindex[%d]: value = %d, max = %d", id, dindex[id], npart));
            stop = kTRUE;
            break;
         }
         // retrieve daughter and copy some info
         if (rsn->IsESD()) {
            daughter = mcEvent->GetTrack(dindex[id]);
            daughterESD = static_cast<AliMCParticle*>(daughter);
            pdg = TMath::Abs(daughterESD->Particle()->GetPdgCode());
            charge = (Short_t)(daughterESD->Particle()->GetPDG()->Charge() / 3);
         } else {
            daughter = (AliAODMCParticle*)listAOD->At(dindex[id]);
            daughterAOD = static_cast<AliAODMCParticle*>(daughter);
            pdg = TMath::Abs(daughterAOD->GetPdgCode());
            charge = (Short_t)daughterAOD->Charge();
         }
         label = TMath::Abs(daughter->GetLabel());
         // check if can be assigned
         imatch = -1;
         if ( fDef->GetDef1().MatchesPDG(pdg) && fDef->GetDef1().MatchesChargeS(charge) )
            imatch = 0;
         else if ( fDef->GetDef2().MatchesPDG(pdg) && fDef->GetDef2().MatchesChargeS(charge) )
            imatch = 1;
         if (imatch < 0 || imatch > 1) continue;
         // assign
         label = daughter->GetLabel();
         fDaughter[imatch].SetRefMC(daughter);
         fDaughter[imatch].SetGood();
         index = FindTrack(label, mcEvent);
         if (index > 0) {
            fDaughter[imatch].SetRef(rsn->GetRef()->GetTrack(index));
         }
      }
      
      // if both daughters were assigned, this means that the decay is correct
      // and then we fill second step
      if (fDaughter[0].IsOK() && fDaughter[1].IsOK()) {
         GetOutput()->Fill(&fMother, 1);
         for (istep = 0; istep < nsteps; istep++) {
            AliRsnCutManager *cuts = (AliRsnCutManager*)fSteps[istep];
            if (!cuts->IsSelected(&fMother)) break;
            GetOutput()->Fill(&fMother, 2 + istep);
         }
      }
   }
   
   return count;
}
