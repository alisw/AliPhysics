//
// Class AliRsnAnalysisTaskEffPair
//
// Inherits from basic AliRsnAnalysisTaskEff for efficiency,
// and computed efficiencies for pairs
//
// author: Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"

#include "AliRsnPairDef.h"
#include "AliRsnCutManager.h"
#include "AliRsnAnalysisTaskEffPair.h"

ClassImp(AliRsnAnalysisTaskEffPair)

//_____________________________________________________________________________
AliRsnAnalysisTaskEffPair::AliRsnAnalysisTaskEffPair(const char *name) :
   AliRsnAnalysisTaskEff(name),
   fMother()
{
//
// Default constructor.
// Do not repeat 'DefineOutput' since it is done in base class and we don't add new ones.
//
}

//_____________________________________________________________________________
AliRsnAnalysisTaskEffPair::AliRsnAnalysisTaskEffPair(const AliRsnAnalysisTaskEffPair& copy) :
   AliRsnAnalysisTaskEff(copy),
   fMother()
{
//
// Copy constrtuctor.
//
}

//_____________________________________________________________________________
AliRsnAnalysisTaskEffPair& AliRsnAnalysisTaskEffPair::operator=(const AliRsnAnalysisTaskEffPair& copy)
{
//
// Assignment operator.
// Owned data members are meaningless for this operator.
//
   
   AliRsnAnalysisTaskEff::operator=(copy);
   return (*this);
}

//_____________________________________________________________________________
Int_t AliRsnAnalysisTaskEffPair::NGoodSteps()
{
//
// Checks how many 'reconstruction' steps are passed by current daughter
//

   Int_t istep, count = 0;
   Int_t nSteps = fStepsRec.GetEntries();
   
   for (istep = 0; istep < nSteps; istep++) {
      AliRsnCutManager *cutMgr = (AliRsnCutManager*)fStepsRec[istep];
      
      if (!cutMgr->PassCommonDaughterCuts(&fDaughter[0])) break;
      if (!cutMgr->PassCommonDaughterCuts(&fDaughter[1])) break;
      if (!cutMgr->PassDaughter1Cuts(&fDaughter[0])) break;
      if (!cutMgr->PassDaughter2Cuts(&fDaughter[1])) break;
      if (!cutMgr->PassMotherCuts(&fMother)) break;
      
      count++;
   }
   
   return count;
}

//_____________________________________________________________________________
void AliRsnAnalysisTaskEffPair::ProcessEventESD()
{
//
// Process current event with the definitions of the specified step in MC list
// and store results in the container slot defined in second argument.
// It is associated with the AliCFContainer with the name of the pair.
//

   AliESDEvent   *esd   = fRsnEvent[0].GetRefESD();
   AliMCEvent    *mc    = fRsnEvent[0].GetRefMCESD();
   AliStack      *stack = mc->Stack();
   TArrayI        indexes[2];
   Int_t          i, j, istep, imax, icheck, itrack[2], ipart;
   Int_t          pdg, label[2];
   Short_t        charge, pairDefMatch[2];
   TParticle     *part = 0x0;
   AliMCParticle *mother = 0x0;
   
   // set pointers
   fMother.SetDaughter(0, &fDaughter[0]);
   fMother.SetDaughter(1, &fDaughter[1]);
   
   // loop on definitions
   AliRsnPairDef *def = 0x0;
   TObjArrayIter nextDef(&fDefs);
   while ( (def = (AliRsnPairDef*)nextDef()) ) {
      
      // loop on the MC list of particles
      for (ipart = 0; ipart < stack->GetNprimary(); ipart++) {
         
         // reset daughters
         fDaughter[0].Reset();
         fDaughter[1].Reset();
         
         // taks MC particle
         mother = (AliMCParticle*)mc->GetTrack(ipart);
         
         // check that it is a binary decay and the PDG code matches
         if (mother->Particle()->GetNDaughters() != 2) continue;
         if (mother->Particle()->GetPdgCode() != def->GetMotherPDG()) continue;
         
         // store the labels of the two daughters
         // and check that they are in the stack
         label[0] = mother->Particle()->GetFirstDaughter();
         label[1] = mother->Particle()->GetLastDaughter();
         if (label[0] < 0 || label[0] > stack->GetNtrack()) continue;
         if (label[1] < 0 || label[1] > stack->GetNtrack()) continue;
         
         // for each daughter, check what slot in the pair definition it matches
         // or if it does not match any of them
         for (i = 0; i < 2; i++) {
            pairDefMatch[i] = -1;
            part = stack->Particle(label[i]);
            if (part) {
               pdg    = TMath::Abs(part->GetPdgCode());
               charge = (Short_t)(part->GetPDG()->Charge() / 3);
               if (def->GetDef1()->MatchesPDG(pdg) && def->GetDef1()->MatchesCharge(charge)) pairDefMatch[i] = 0;
               if (def->GetDef2()->MatchesPDG(pdg) && def->GetDef2()->MatchesCharge(charge)) pairDefMatch[i] = 1;
            }
         }
         
         // if the two label match the two definitions for the pair
         // and if they are in the wrong order, swap them,
         // otherwise, if they don't match the definition in any order, skip
         if (pairDefMatch[0] == 1 && pairDefMatch[1] == 0) {
            icheck = label[0];
            label[0] = label[1];
            label[1] = icheck;
         }
         else if (pairDefMatch[0] < 0 || pairDefMatch[1] < 0) continue;
         
         // from now on, we are sure that label[0] refers to the particle
         // that matches definitions of first daughter, and label[1] to
         // the particle that matches definitions of second daughter
         fDaughter[0].SetRefMC(mc->GetTrack(label[0]));
         fDaughter[1].SetRefMC(mc->GetTrack(label[1]));
         
         // assign masses and fill the MC steps,
         // where reconstruction is not taken into account
         fMother.ComputeSum(def->GetMass1(), def->GetMass2());
         FillContainer(kTRUE, def);
         
         // search for all reconstructed tracks which have these labels
         for (i = 0; i < 2; i++) indexes[i] = FindTracks(label[i], esd);
         
         // if not both tracks have been reconstructed, stop here
         if (indexes[0].GetSize() < 1 || indexes[1].GetSize() < 1) continue;
         
         // if both daughters were reconstructed
         // search for the best combination of indexes for this pair
         imax = itrack[0] = itrack[1] = 0;
         for (i = 0; i < indexes[0].GetSize(); i++) {
            for (j = 0; j < indexes[1].GetSize(); j++) {
               fDaughter[0].SetRef(esd->GetTrack(indexes[0][i]));
               fDaughter[1].SetRef(esd->GetTrack(indexes[1][j]));
               fMother.ComputeSum(def->GetMass1(), def->GetMass2());
               istep = NGoodSteps();
               if (istep > imax) {
                  itrack[0] = indexes[0][i];
                  itrack[1] = indexes[1][j];
                  imax = istep;
               }
            }
         }
         
         // then assign definitely the best combination and fill rec container
         fDaughter[0].SetRef(esd->GetTrack(itrack[0]));
         fDaughter[1].SetRef(esd->GetTrack(itrack[1]));
         fMother.ComputeSum(def->GetMass1(), def->GetMass2());
         FillContainer(kFALSE, def);
      }
   }
}

//_____________________________________________________________________________
void AliRsnAnalysisTaskEffPair::ProcessEventAOD()
{
//
// Process current event with the definitions of the specified step in MC list
// and store results in the container slot defined in second argument.
// It is associated with the AliCFContainer with the name of the pair.
//

   AliAODEvent  *aod     = fRsnEvent[0].GetRefAOD();
   TClonesArray *mcArray = (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
   if (!mcArray) return;
   
   TArrayI           indexes[2];
   Int_t             i, j, pdg, imax, istep, icheck, itrack[2], label[2];
   AliAODMCParticle *part = 0x0, *mother = 0x0;
   Short_t           charge = 0, pairDefMatch[2];
   
   // set pointers
   fMother.SetDaughter(0, &fDaughter[0]);
   fMother.SetDaughter(1, &fDaughter[1]);
   
   // loop on definitions
   AliRsnPairDef *def = 0x0;
   TObjArrayIter nextDef(&fDefs);
   while ( (def = (AliRsnPairDef*)nextDef()) ) {
      
      // loop on the MC list of particles
      TObjArrayIter next(mcArray);
      while ((mother = (AliAODMCParticle*)next())) {
         
         // check that it is a binary decay and the PDG code matches
         if (mother->GetNDaughters() != 2) continue;
         if (mother->GetPdgCode() != def->GetMotherPDG()) continue;
         
         // store the labels of the two daughters
         label[0] = mother->GetDaughter(0);
         label[1] = mother->GetDaughter(1);
         
         // for each daughter, check what slot in the pair definition it matches
         // or if it does not match any of them
         for (i = 0; i < 2; i++) {
            pairDefMatch[i] = -1;
            part = (AliAODMCParticle*)mcArray->At(label[i]);
            if (part) {
               pdg    = TMath::Abs(part->GetPdgCode());
               charge = (Short_t)part->Charge();
               if (def->GetDef1()->MatchesPDG(pdg) && def->GetDef1()->MatchesCharge(charge)) pairDefMatch[i] = 0;
               if (def->GetDef2()->MatchesPDG(pdg) && def->GetDef2()->MatchesCharge(charge)) pairDefMatch[i] = 1;
            }
         }
         
         // if the two label match the two definitions for the pair
         // and if they are in the wrong order, swap them,
         // otherwise, if they don't match the definition in any order, skip
         if (pairDefMatch[0] == 1 && pairDefMatch[1] == 0) {
            icheck = label[0];
            label[0] = label[1];
            label[1] = icheck;
         }
         else if (pairDefMatch[0] < 0 || pairDefMatch[1] < 0) continue;
         
         // from now on, we are sure that label[0] refers to the particle
         // that matches definitions of first daughter, and label[1] to
         // the particle that matches definitions of second daughter
         fDaughter[0].SetRefMC((AliAODMCParticle*)mcArray->At(label[0]));
         fDaughter[1].SetRefMC((AliAODMCParticle*)mcArray->At(label[1]));
         
         // assign masses and fill the MC steps,
         // where reconstruction is not taken into account
         fMother.ComputeSum(def->GetMass1(), def->GetMass2());
         FillContainer(kTRUE, def);
         
         // search for all reconstructed tracks which have these labels
         for (i = 0; i < 2; i++) indexes[i] = FindTracks(label[i], aod);
         
         // if not both tracks have been reconstructed, stop here
         if (indexes[0].GetSize() < 1 || indexes[1].GetSize() < 1) continue;
         
         // if both daughters were reconstructed
         // search for the best combination of indexes for this pair
         imax = itrack[0] = itrack[1] = 0;
         for (i = 0; i < indexes[0].GetSize(); i++) {
            for (j = 0; j < indexes[1].GetSize(); j++) {
               fDaughter[0].SetRef(aod->GetTrack(indexes[0][i]));
               fDaughter[1].SetRef(aod->GetTrack(indexes[1][j]));
               fMother.ComputeSum(def->GetMass1(), def->GetMass2());
               istep = NGoodSteps();
               if (istep > imax) {
                  itrack[0] = indexes[0][i];
                  itrack[1] = indexes[1][j];
               }
            }
         }
         
         // then assign definitely the best combination and fill rec container
         fDaughter[0].SetRef(aod->GetTrack(itrack[0]));
         fDaughter[1].SetRef(aod->GetTrack(itrack[1]));
         fMother.ComputeSum(def->GetMass1(), def->GetMass2());
         FillContainer(kFALSE, def);
      }
   }
}

//_____________________________________________________________________________
void AliRsnAnalysisTaskEffPair::FillContainer(Bool_t mcList, TObject *defObj)
{
//
// Fill the container corresponding to current definition.
//
   
   // cast the object into the def
   if (!defObj->InheritsFrom(AliRsnPairDef::Class())) {
      AliError("Def object does not inherit from AliRsnPairDef!");
      return;
   }
   AliRsnPairDef *def = static_cast<AliRsnPairDef*>(defObj);
   
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
         case AliRsnTarget::kMother:
            fcnAxis->SetSupportObject(def);
            computeOK = fcnAxis->Eval(&fMother, mcList);
            break;
         case AliRsnTarget::kEvent:
            computeOK = fcnAxis->Eval(&fRsnEvent[0]);
            break;
         default:
            AliError(Form("Allowed targets are mothers and events; cannot use axis '%s' which has target '%s'", fcnAxis->GetName(), fcnAxis->GetTargetTypeName()));
            computeOK = kFALSE;
      }
      if (computeOK) fVar[iaxis] = ((Float_t)fcnAxis->GetComputedValue());
   }

   // fill all successful steps
   for (istep = 0; istep < nSteps; istep++) {
      AliRsnCutManager *cutMgr = (AliRsnCutManager*)stepList[istep];
      
      if (!cutMgr->PassCommonDaughterCuts(&fDaughter[0])) break;
      if (!cutMgr->PassCommonDaughterCuts(&fDaughter[1])) break;
      if (!cutMgr->PassDaughter1Cuts(&fDaughter[0])) break;
      if (!cutMgr->PassDaughter2Cuts(&fDaughter[1])) break;
      if (!cutMgr->PassMotherCuts(&fMother)) break;
      
      AliDebug(AliLog::kDebug + 2, Form("DEF: %s --> filling step %d", def->GetName(), istep));
      cont->Fill(fVar.GetArray(), istep + firstStep);
   }
}
