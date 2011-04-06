//
// Class AliRsnLoopEffDaughter
//
// Inherits from basic AliRsnLoopEff for efficiency,
// and computed efficiencies for single-tracks
//
// author: Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include "AliStack.h"

#include "AliRsnDaughterDef.h"
#include "AliRsnCutSet.h"

#include "AliRsnLoopEffDaughter.h"

ClassImp(AliRsnLoopEffDaughter)

//_____________________________________________________________________________
AliRsnLoopEffDaughter::AliRsnLoopEffDaughter(const char *name, AliRsnDaughterDef *def) :
   AliRsnLoopEff(name),
   fDef(def)
{
//
// Default constructor.
// Do not repeat 'DefineOutput' since it is done in base class and we don't add new ones.
//
}

//_____________________________________________________________________________
AliRsnLoopEffDaughter::AliRsnLoopEffDaughter(const AliRsnLoopEffDaughter& copy) :
   AliRsnLoopEff(copy),
   fDef(copy.fDef)
{
//
// Copy constrtuctor.
//
}

//_____________________________________________________________________________
AliRsnLoopEffDaughter& AliRsnLoopEffDaughter::operator=(const AliRsnLoopEffDaughter& copy)
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
Bool_t AliRsnLoopEffDaughter::OkStepMC(TObject *target, Int_t istep)
{
//
// Check step with MC
//

   AliRsnCutSet *cuts = (AliRsnCutSet*)fStepsMC[istep];
   return cuts->IsSelected(target);
}

//_____________________________________________________________________________
Bool_t AliRsnLoopEffDaughter::OkStepRec(TObject *target, Int_t istep)
{
//
// Check step with MC
//

   AliRsnCutSet *cuts = (AliRsnCutSet*)fStepsRec[istep];
   return cuts->IsSelected(target);
}

//_____________________________________________________________________________
Int_t AliRsnLoopEffDaughter::ProcessEventESD(AliRsnEvent *rsn)
{
//
// Process current event with the definitions of the specified step in MC list
// and store results in the container slot defined in second argument.
// It is associated with the AliCFContainer with the name of the pair.
//

   AliESDEvent *esd   = rsn->GetRefESD();
   AliMCEvent  *mc    = rsn->GetRefMCESD();
   TArrayI      indexes;
   Int_t        imax, istep, icheck, itrack, ipart;
   Int_t        nsteps = NSteps();

   static AliRsnDaughter daughter;
   daughter.SetOwnerEvent(rsn);

   // loop on the MC list of particles
   for (ipart = 0; ipart < stack->GetNprimary(); ipart++) {

      // MC particle
      daughter.SetRefMC((AliMCParticle*)mc->GetTrack(ipart));

      // search for reconstructed track
      // if no tracks are found with that label, rec ref is set to zero
      // if more than one tracks are found we use the one which passes
      // most cut steps
      indexes = FindTracks(ipart, esd);
      if (indexes.GetSize() < 1)
         daughter.SetRef(0x0);
      else if (indexes.GetSize() == 1)
         daughter.SetRef(esd->GetTrack(indexes[0]));
      else {
         imax = istep = itrack = 0;
         for (icheck = 0; icheck < indexes.GetSize(); icheck++) {
            daughter.SetRef(esd->GetTrack(indexes[icheck]));
            daughter.SetMass(def->GetMass());
            istep = NGoodSteps(&daughter);
            if (istep > imax) itrack = icheck;
         }
         daughter.SetRef(esd->GetTrack(indexes[itrack]));
      }

      // compute 4-momenta
      daughter.SetMass(def->GetMass());

      // fill MC container
      for (istep = 0; istep < nsteps; istep++) {
         if (!OkStep(&daughter, istep)) break;
         GetOutput()->Fill(&daughter, istep);
      }
   }
}

//_____________________________________________________________________________
Int_t AliRsnLoopEffDaughter::ProcessEventAOD(AliRsnEvent *rsn)
{
//
// Process current event with the definitions of the specified step in MC list
// and store results in the container slot defined in second argument.
// It is associated with the AliCFContainer with the name of the pair.
//

   AliAODEvent  *aod     = fRsnEvent[0].GetRefAOD();
   TClonesArray *mcArray = (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
   if (!mcArray) return;
   TArrayI       indexes;
   Int_t         imax, istep, icheck, itrack, ipart;

   static AliRsnDaughter daughter;
   daughter->SetOwnerEvent(rsn);

   // loop on the MC list of particles
   TObjArrayIter next(mcArray);
   AliAODMCParticle *particle;
   while ((particle = (AliAODMCParticle*)next())) {

      // MC particle
      daughter.SetRefMC(particle);

      // search for reconstructed track
      // if no tracks are found with that label, rec ref is set to zero
      // if more than one tracks are found we use the one which passes
      // most cut steps
      ipart = particle->GetLabel();
      indexes = FindTracks(ipart, aod);
      if (indexes.GetSize() < 1)
         daughter.SetRef(0x0);
      else if (indexes.GetSize() == 1)
         daughter.SetRef(aod->GetTrack(indexes[0]));
      else {
         imax = istep = itrack = 0;
         for (icheck = 0; icheck < indexes.GetSize(); icheck++) {
            daughter.SetRef(aod->GetTrack(indexes[icheck]));
            daughter.SetMass(def->GetMass());
            istep = NGoodSteps();
            if (istep > imax) {
               itrack = icheck;
               imax = istep;
            }
         }
         daughter.SetRef(aod->GetTrack(indexes[itrack]));
      }

      // compute 4-momenta
      daughter.SetMass(def->GetMass());

      // fill MC container
      for (istep = 0; istep < nsteps; istep++) {
         if (!OkStep(&daughter, istep)) break;
         GetOutput()->Fill(&daughter, istep);
      }
   }
}
