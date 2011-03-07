//
// Class AliRsnAnalysisTaskEffMonitor
//
// Inherits from basic AliRsnAnalysisTaskEff for efficiency,
// and computed efficiencies for single-tracks
//
// author: Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"

#include "AliRsnDaughterDef.h"
#include "AliRsnAnalysisTaskEffMonitor.h"

ClassImp(AliRsnAnalysisTaskEffMonitor)

//_____________________________________________________________________________
AliRsnAnalysisTaskEffMonitor::AliRsnAnalysisTaskEffMonitor(const char *name) :
   AliRsnAnalysisTaskEff(name),
   fDaughter()
{
//
// Default constructor.
// Do not repeat 'DefineOutput' since it is done in base class and we don't add new ones.
//
}

//_____________________________________________________________________________
AliRsnAnalysisTaskEffMonitor::AliRsnAnalysisTaskEffMonitor(const AliRsnAnalysisTaskEffMonitor& copy) :
   AliRsnAnalysisTaskEff(copy),
   fDaughter()
{
//
// Copy constrtuctor.
//
}

//_____________________________________________________________________________
AliRsnAnalysisTaskEffMonitor& AliRsnAnalysisTaskEffMonitor::operator=(const AliRsnAnalysisTaskEffMonitor& copy)
{
//
// Assignment operator.
// Owned data members are meaningless for this operator.
//
   
   AliRsnAnalysisTaskEff::operator=(copy);
   return (*this);
}

//_____________________________________________________________________________
Int_t AliRsnAnalysisTaskEffMonitor::NGoodSteps()
{
//
// Checks how many 'reconstruction' steps are passed by current daughter
//

   Int_t istep, count = 0;
   Int_t nSteps = fStepsRec.GetEntries();
   
   for (istep = 0; istep < nSteps; istep++) {
      AliRsnCutSet *cutSet = (AliRsnCutSet*)fStepsRec[istep];
      if (!cutSet->IsSelected(&fDaughter)) break;
      count++;
   }
   
   return count;
}

//_____________________________________________________________________________
void AliRsnAnalysisTaskEffMonitor::ProcessEventESD()
{
//
// Process current event with the definitions of the specified step in MC list
// and store results in the container slot defined in second argument.
// It is associated with the AliCFContainer with the name of the pair.
//

   AliESDEvent   *esd   = fRsnEvent[0].GetRefESD();
   AliMCEvent    *mc    = fRsnEvent[0].GetRefMCESD();
   AliStack      *stack = mc->Stack();
   TArrayI        indexes;
   Int_t          imax, istep, icheck, itrack, ipart;
   
   // loop on definitions
   AliRsnDaughterDef *def = 0x0;
   TObjArrayIter nextDef(&fDefs);
   while ( (def = (AliRsnDaughterDef*)nextDef()) ) {

      // loop on the MC list of particles
      for (ipart = 0; ipart < stack->GetNprimary(); ipart++) {
         
         // MC particle
         fDaughter.SetRefMC((AliMCParticle*)mc->GetTrack(ipart));
         
         // search for reconstructed track
         // if no tracks are found with that label, rec ref is set to zero
         // if more than one tracks are found we use the one which passes
         // most cut steps
         indexes = FindTracks(ipart, esd);
         if (indexes.GetSize() < 1)
            fDaughter.SetRef(0x0);
         else if (indexes.GetSize() == 1)
            fDaughter.SetRef(esd->GetTrack(indexes[0]));
         else {
            imax = istep = itrack = 0;
            for (icheck = 0; icheck < indexes.GetSize(); icheck++) {
               fDaughter.SetRef(esd->GetTrack(indexes[icheck]));
               fDaughter.SetMass(def->GetMass());
               istep = NGoodSteps();
               if (istep > imax) itrack = icheck;
            }
            fDaughter.SetRef(esd->GetTrack(indexes[itrack]));
         }
         
         // compute 4-momenta
         fDaughter.SetMass(def->GetMass());

         // fill MC container
         FillContainer(kTRUE, def);
         if (fDaughter.GetRef()) FillContainer(kFALSE, def);
      }
   }
}

//_____________________________________________________________________________
void AliRsnAnalysisTaskEffMonitor::ProcessEventAOD()
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
   
   // loop on definitions
   AliRsnDaughterDef *def = 0x0;
   TObjArrayIter nextDef(&fDefs);
   while ( (def = (AliRsnDaughterDef*)nextDef()) ) {

      // loop on the MC list of particles
      TObjArrayIter next(mcArray);
      AliAODMCParticle *particle;
      while ((particle = (AliAODMCParticle*)next())) {
         
         // MC particle
         fDaughter.SetRefMC(particle);
         
         // search for reconstructed track
         // if no tracks are found with that label, rec ref is set to zero
         // if more than one tracks are found we use the one which passes
         // most cut steps
         ipart = particle->GetLabel();
         indexes = FindTracks(ipart, aod);
         if (indexes.GetSize() < 1)
            fDaughter.SetRef(0x0);
         else if (indexes.GetSize() == 1)
            fDaughter.SetRef(aod->GetTrack(indexes[0]));
         else {
            imax = istep = itrack = 0;
            for (icheck = 0; icheck < indexes.GetSize(); icheck++) {
               fDaughter.SetRef(aod->GetTrack(indexes[icheck]));
               fDaughter.SetMass(def->GetMass());
               istep = NGoodSteps();
               if (istep > imax) {
                  itrack = icheck;
                  imax = istep;
               }
            }
            fDaughter.SetRef(aod->GetTrack(indexes[itrack]));
         }
         
         // compute 4-momenta
         fDaughter.SetMass(def->GetMass());

         // fill MC container
         FillContainer(kTRUE, def);
         if (fDaughter.GetRef()) FillContainer(kFALSE, def);
      }
   }
}

//_____________________________________________________________________________
void AliRsnAnalysisTaskEffMonitor::FillContainer(Bool_t mcList, TObject *defObj)
{
//
// Fill the container corresponding to current definition.
//

   // cast the object into the def
   if (!defObj->InheritsFrom(AliRsnDaughterDef::Class())) {
      AliError("Def object does not inherit from AliRsnDaughterDef!");
      return;
   }
   AliRsnDaughterDef *def = static_cast<AliRsnDaughterDef*>(defObj);

   // retrieve container
   AliCFContainer *cont = (AliCFContainer*)fOutList->FindObject(def->GetName());
   if (!cont) return;

   TObjArray &stepList =  (mcList ? fStepsMC : fStepsRec);
   Int_t      firstStep = (mcList ? 0 : ((Int_t)fStepsMC.GetEntries()));
   Int_t      iaxis, nAxes  = fAxes.GetEntries();
   Int_t      istep, nSteps = stepList.GetEntries();
   Bool_t     computeOK;
   
   // compute values for all axes
   for (iaxis = 0; iaxis < nAxes; iaxis++) {
      AliRsnValue *fcnAxis = (AliRsnValue*)fAxes.At(iaxis);
      fVar[iaxis] = -1E10;
      switch (fcnAxis->GetTargetType()) {
         case AliRsnTarget::kDaughter:
            fcnAxis->SetSupportObject(def);
            computeOK = fcnAxis->Eval(&fDaughter, mcList);
            break;
         case AliRsnTarget::kEvent:
            computeOK = fcnAxis->Eval(&fRsnEvent[0]);
            break;
         default:
            AliError(Form("Allowed targets are daughters and events; cannot use axis '%s' which has target '%s'", fcnAxis->GetName(), fcnAxis->GetTargetTypeName()));
            computeOK = kFALSE;
      }
      if (computeOK) fVar[iaxis] = ((Float_t)fcnAxis->GetComputedValue());
   }

   // fill all successful steps
   for (istep = 0; istep < nSteps; istep++) {
      AliRsnCutSet *cutSet = (AliRsnCutSet*)stepList[istep];
      
      if (!cutSet->IsSelected(&fDaughter)) break;
      
      AliDebug(AliLog::kDebug + 2, Form("DEF: %s --> filling step %d", def->GetName(), istep));
      cont->Fill(fVar.GetArray(), istep + firstStep);
   }
}
