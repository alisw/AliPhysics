//
// Class AliRsnAnalysisEffSE
//
// TODO
// TODO
//
// authors: Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//          Martin Vala (martin.vala@cern.ch)
//
#include <Riostream.h>
#include "AliStack.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"

#include "AliCFContainer.h"
#include "AliRsnCutManager.h"
#include "AliRsnValue.h"
#include "AliRsnAnalysisEffSE.h"
#include "AliRsnPairDef.h"
#include "AliRsnCutSet.h"

ClassImp(AliRsnAnalysisEffSE)

//_____________________________________________________________________________
AliRsnAnalysisEffSE::AliRsnAnalysisEffSE(const char *name) :
   AliRsnVAnalysisTaskSE(name),
   fUseITSSA(kTRUE),
   fUseGlobal(kTRUE),
   fStepListMC(0),
   fStepListESD(0),
   fAxisList("AliRsnValue", 0),
   fPairDefList(0),
   fContainerList(0x0),
   fOutList(0x0),
   fVar(0),
   fPair(),
   fEventCuts("eventCuts", AliRsnCut::kEvent)
{
//
// Default constructor.
//

   AliDebug(AliLog::kDebug + 2, "<-");

   DefineOutput(2, TList::Class());

   AliDebug(AliLog::kDebug + 2, "->");
}

//_____________________________________________________________________________
AliRsnAnalysisEffSE::AliRsnAnalysisEffSE(const AliRsnAnalysisEffSE& copy) :
   AliRsnVAnalysisTaskSE(copy),
   fUseITSSA(copy.fUseITSSA),
   fUseGlobal(copy.fUseGlobal),
   fStepListMC(copy.fStepListMC),
   fStepListESD(copy.fStepListESD),
   fAxisList(copy.fAxisList),
   fPairDefList(copy.fPairDefList),
   fContainerList(copy.fContainerList),
   fOutList(0x0),
   fVar(0),
   fPair(),
   fEventCuts(copy.fEventCuts)
{
//
// Copy constrtuctor
//
}

//_____________________________________________________________________________
void AliRsnAnalysisEffSE::AddStepMC(AliRsnCutManager *mgr)
{
//
// Add a step on montecarlo
//

   fStepListMC.AddLast(mgr);
}

//_____________________________________________________________________________
void AliRsnAnalysisEffSE::AddStepESD(AliRsnCutManager *mgr)
{
//
// Add a step on ESD
//

   fStepListESD.AddLast(mgr);
}

//_____________________________________________________________________________
void AliRsnAnalysisEffSE::AddAxis(AliRsnValue *axis)
{
//
// Add a new axis
//

   Int_t n = fAxisList.GetEntries();
   new(fAxisList[n]) AliRsnValue(*axis);
}

//_____________________________________________________________________________
void AliRsnAnalysisEffSE::RsnUserCreateOutputObjects()
{
//
// Creation of output objects.
// These are created through the utility methods in the analysis manager,
// which produces a list of histograms for each specified set of pairs.
// Each of these lists is added to the main list of this task.
//

   AliDebug(AliLog::kDebug + 2, "<-");

   // get number of steps and axes
   Int_t iaxis  = 0;
   Int_t nAxes  = fAxisList.GetEntries();
   Int_t nSteps = (Int_t)fStepListMC.GetEntries() + (Int_t)fStepListESD.GetEntries();

   if (!nSteps) {
      AliError("No steps defined");
      return;
   }
   if (!nAxes) {
      AliError("No axes defined");
      return;
   }

   // initialize variable list
   fVar.Set(nAxes);

   // retrieve number of bins for each axis
   Int_t *nBins = new Int_t[nAxes];
   for (iaxis = 0; iaxis < nAxes; iaxis++) {
      AliRsnValue *fcnAxis = (AliRsnValue*)fAxisList.At(iaxis);
      nBins[iaxis] = fcnAxis->GetArray().GetSize() - 1;
   }

   // create ouput list of containers
   fContainerList = new TList;
   fContainerList->SetOwner();
   fContainerList->SetName(Form("%s_containers", GetName()));

   // initialize output list
   OpenFile(2);
   fOutList = new TList();
   fOutList->SetOwner();

   // create the containers
   Int_t i = 0, nDef = (Int_t)fPairDefList.GetEntries();
   for (i = 0; i < nDef; i++) {
      AliRsnPairDef *def = (AliRsnPairDef*)fPairDefList[i];
      AliCFContainer *cont = new AliCFContainer(Form("%s", def->GetPairName()), "", nSteps, nAxes, nBins);
      // set the bin limits for each axis
      for (iaxis = 0; iaxis < nAxes; iaxis++) {
         AliRsnValue *fcnAxis = (AliRsnValue*)fAxisList.At(iaxis);
         cont->SetBinLimits(iaxis, fcnAxis->GetArray().GetArray());
      }
      // add the container to output list
      fContainerList->Add(cont);
   }

   fOutList->Add(fContainerList);
   fOutList->Print();

   PostData(2, fOutList);

   // clear heap
   delete [] nBins;

   AliDebug(AliLog::kDebug + 2, "->");
}

//_____________________________________________________________________________
void AliRsnAnalysisEffSE::RsnUserExec(Option_t*)
{
//
// Execution of the analysis task.
// Recovers the input event and processes it with all included pair objects.
// In this case, we NEED to have ESD and MC, otherwise we cannod do anything.
//

   // process first MC steps and then ESD steps
   AliRsnPairDef *pairDef = 0;
   TObjArrayIter iter(&fPairDefList);
   while ((pairDef = (AliRsnPairDef*)iter.Next())) {
      if (fRsnEvent.IsESD()) ProcessEventESD(pairDef);
      if (fRsnEvent.IsAOD()) ProcessEventAOD(pairDef);
   }

   // Post the data
   PostData(2, fOutList);
}

//_____________________________________________________________________________
void AliRsnAnalysisEffSE::ProcessEventESD(AliRsnPairDef *pairDef)
{
//
// Process current event with the definitions of the specified step in MC list
// and store results in the container slot defined in second argument.
// It is associated with the AliCFContainer with the name of the pair.
//

   AliStack      *stack = fRsnEvent.GetRefMCESD()->Stack();
   AliESDEvent   *esd   = fRsnEvent.GetRefESD();
   AliMCEvent    *mc    = fRsnEvent.GetRefMCESD();
   AliMCParticle *mother;

   if (!pairDef) return;
   AliCFContainer *cont = (AliCFContainer*)fContainerList->FindObject(pairDef->GetPairName());
   if (!cont) return;

   // get informations from pairDef
   Int_t   pdgM = 0;
   Int_t   pdgD[2] = {0, 0};
   Short_t chargeD[2] = {0, 0};
   pdgM       = pairDef->GetMotherPDG();
   pdgD   [0] = AliPID::ParticleCode(pairDef->GetPID(0));
   pdgD   [1] = AliPID::ParticleCode(pairDef->GetPID(1));
   chargeD[0] = pairDef->GetChargeS(0);
   chargeD[1] = pairDef->GetChargeS(1);

   // other utility variables
   Int_t      first, j, ipart;
   Int_t      label[2] = { -1, -1};
   Short_t    charge[2] = {0, 0};
   Short_t    pairDefMatch[2] = { -1, -1};
   Int_t      esdIndex[2] = { -1, -1};
   TParticle *part[2] = {0, 0};

   // in this case, we first take the resonance from MC
   // and then we find its daughters and compute cuts on them
   for (ipart = 0; ipart < stack->GetNprimary(); ipart++) {
      // take a track from the MC event
      mother = (AliMCParticle*) fMCEvent->GetTrack(ipart);

      // check that it is a binary decay and the PDG code matches
      if (mother->Particle()->GetNDaughters() != 2) continue;
      if (mother->Particle()->GetPdgCode() != pdgM) continue;

      // store the labels of the two daughters
      label[0] = mother->Particle()->GetFirstDaughter();
      label[1] = mother->Particle()->GetLastDaughter();

      // retrieve the particles and other stuff
      // check if they match the order in the pairDef
      for (j = 0; j < 2; j++) {
         if (label[j] < 0) continue;
         part[j]   = stack->Particle(label[j]);
         pdgD[j]   = TMath::Abs(part[j]->GetPdgCode());
         charge[j] = (Short_t)(part[j]->GetPDG()->Charge() / 3);
         if (pdgD[j] == AliPID::ParticleCode(pairDef->GetPID(0)) && charge[j] == pairDef->GetChargeS(0))
            pairDefMatch[j] = 0;
         else if (pdgD[j] == AliPID::ParticleCode(pairDef->GetPID(1)) && charge[j] == pairDef->GetChargeS(1))
            pairDefMatch[j] = 1;
         else
            pairDefMatch[j] = -1;

         // find corresponding ESD particle: first try rejecting fakes,
         // and in case of failure, try accepting fakes
         esdIndex[j] = FindESDtrack(label[j], esd, kTRUE);
         //TArrayI idx = FindESDtracks(label[j], esd);
         //for (Int_t kk = 0; kk < idx.GetSize(); kk++) cout << "DAUGHTER " << j << " --> FOUND INDEX: " << idx[kk] << endl;
         if (esdIndex[j] < 0) esdIndex[j] = FindESDtrack(label[j], esd, kFALSE);
         //cout << "DAUGHTER " << j << " SINGLE FOUND INDEX = " << esdIndex[j] << endl;
      }

      // since each candidate good resonance is taken once, we must check
      // that it matches the pair definition in any order, and reject in case
      // in none of them the pair is OK
      // anyway, we must associate the correct daughter to the correct data member
      if (pairDefMatch[0] == 0 && pairDefMatch[1] == 1) {
         // 1st track --> 1st member of PairDef
         fDaughter[0].SetRef(mc->GetTrack(label[0]));
         fDaughter[0].SetRefMC((AliMCParticle*)mc->GetTrack(label[0]));
         fDaughter[0].SetGood();
         // 2nd track --> 2nd member of PairDef
         fDaughter[1].SetRef(mc->GetTrack(label[1]));
         fDaughter[1].SetRefMC((AliMCParticle*)mc->GetTrack(label[1]));
         fDaughter[1].SetGood();
      } else if ((pairDefMatch[0] == 1 && pairDefMatch[1] == 0)) {
         // 1st track --> 2nd member of PairDef
         fDaughter[0].SetRef(mc->GetTrack(label[1]));
         fDaughter[0].SetRefMC((AliMCParticle*)mc->GetTrack(label[1]));
         fDaughter[0].SetGood();
         // 2nd track --> 1st member of PairDef
         fDaughter[1].SetRef(mc->GetTrack(label[0]));
         fDaughter[1].SetRefMC((AliMCParticle*)mc->GetTrack(label[0]));
         fDaughter[1].SetGood();
      } else {
         fDaughter[0].SetBad();
         fDaughter[1].SetBad();
      }

      // reject the pair if the matching was unsuccessful
      if (!fDaughter[0].IsOK() || !fDaughter[1].IsOK()) continue;

      // first, we set the internal AliRsnMother object to
      // the MC particles and then operate the selections on MC
      fPair.SetDaughters(&fDaughter[0], pairDef->GetMass(0), &fDaughter[1], pairDef->GetMass(1));
      FillContainer(cont, &fStepListMC, pairDef, 0);

      // then, if both particles found a good match in the ESD
      // reassociate the ESD tracks to the pair and fill ESD containers
      if (esdIndex[0] < 0 || esdIndex[1] < 0) continue;
      if (pairDefMatch[0] == 0 && pairDefMatch[1] == 1) {
         // 1st track --> 1st member of PairDef
         fDaughter[0].SetRef(esd->GetTrack(esdIndex[0]));
         // 2nd track --> 2nd member of PairDef
         fDaughter[1].SetRef(esd->GetTrack(esdIndex[1]));
         //cout << "****** MATCHING SCHEME 1" << endl;
      } else if ((pairDefMatch[0] == 1 && pairDefMatch[1] == 0)) {
         // 1st track --> 2nd member of PairDef
         fDaughter[0].SetRef(esd->GetTrack(esdIndex[1]));
         // 2nd track --> 1st member of PairDef
         fDaughter[1].SetRef(esd->GetTrack(esdIndex[0]));
         //cout << "****** MATCHING SCHEME 2" << endl;
      }
      //cout << "****** IDs = " << fDaughter[0].GetID() << ' ' << fDaughter[1].GetID() << endl;
      // here we must remember how many steps were already filled
      first = (Int_t)fStepListMC.GetEntries();
      FillContainer(cont, &fStepListESD, pairDef, first);
   }
}

//_____________________________________________________________________________
void AliRsnAnalysisEffSE::ProcessEventAOD(AliRsnPairDef *pairDef)
{
//
// Process current event with the definitions of the specified step in MC list
// and store results in the container slot defined in second argument.
// It is associated with the AliCFContainer with the name of the pair.
//

   AliAODEvent      *aod   = fRsnEvent.GetRefAOD();
   AliAODMCParticle *mother;
   TClonesArray     *mcArray = (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
   if (!mcArray) return;

   if (!pairDef) return;
   AliCFContainer *cont = (AliCFContainer*)fContainerList->FindObject(pairDef->GetPairName());
   if (!cont) return;

   // get informations from pairDef
   Int_t   pdgM = 0;
   Int_t   pdgD[2] = {0, 0};
   Short_t chargeD[2] = {0, 0};
   pdgM       = pairDef->GetMotherPDG();
   pdgD   [0] = AliPID::ParticleCode(pairDef->GetPID(0));
   pdgD   [1] = AliPID::ParticleCode(pairDef->GetPID(1));
   chargeD[0] = pairDef->GetChargeS(0);
   chargeD[1] = pairDef->GetChargeS(1);

   // other utility variables
   Int_t             first, j;
   Int_t             label [2] = { -1, -1};
   Short_t           charge[2] = {0, 0};
   Short_t           pairDefMatch[2] = { -1, -1};
   Int_t             aodIndex[2] = { -1, -1};
   AliAODMCParticle *part[2] = {0, 0};

   // loop on MC particles
   TObjArrayIter next(mcArray);
   while ((mother = (AliAODMCParticle*)next())) {
      // check that it is a binary decay and the PDG code matches
      if (mother->GetNDaughters() != 2) continue;
      if (mother->GetPdgCode() != pdgM) continue;

      // store the labels of the two daughters
      label[0] = mother->GetDaughter(0);
      label[1] = mother->GetDaughter(1);

      // retrieve the particles and other stuff
      // check if they match the order in the pairDef
      for (j = 0; j < 2; j++) {
         if (label[j] < 0) continue;
         part[j]   = (AliAODMCParticle*)mcArray->At(label[j]);
         pdgD[j]   = TMath::Abs(part[j]->GetPdgCode());
         charge[j] = (Short_t)(part[j]->Charge());
         if (pdgD[j] == AliPID::ParticleCode(pairDef->GetPID(0)) && charge[j] == pairDef->GetChargeS(0))
            pairDefMatch[j] = 0;
         else if (pdgD[j] == AliPID::ParticleCode(pairDef->GetPID(1)) && charge[j] == pairDef->GetChargeS(1))
            pairDefMatch[j] = 1;
         else
            pairDefMatch[j] = -1;

         // find corresponding ESD particle: first try rejecting fakes,
         // and in case of failure, try accepting fakes
         aodIndex[j] = FindAODtrack(label[j], aod, kTRUE);
         if (aodIndex[j] < 0) aodIndex[j] = FindAODtrack(label[j], aod, kFALSE);
      }

      // since each candidate good resonance is taken once, we must check
      // that it matches the pair definition in any order, and reject in case
      // in none of them the pair is OK
      // anyway, we must associate the correct daughter to the correct data member
      if (pairDefMatch[0] == 0 && pairDefMatch[1] == 1) {
         // 1st track --> 1st member of PairDef
         fDaughter[0].SetRef((AliAODMCParticle*)mcArray->At(label[0]));
         fDaughter[0].SetRefMC((AliAODMCParticle*)mcArray->At(label[0]));
         fDaughter[0].SetGood();
         // 2nd track --> 2nd member of PairDef
         fDaughter[1].SetRef((AliAODMCParticle*)mcArray->At(label[1]));
         fDaughter[1].SetRefMC((AliAODMCParticle*)mcArray->At(label[1]));
         fDaughter[1].SetGood();
      } else if ((pairDefMatch[0] == 1 && pairDefMatch[1] == 0)) {
         // 1st track --> 2nd member of PairDef
         fDaughter[0].SetRef((AliAODMCParticle*)mcArray->At(label[1]));
         fDaughter[0].SetRefMC((AliAODMCParticle*)mcArray->At(label[1]));
         fDaughter[0].SetGood();
         // 2nd track --> 1st member of PairDef
         fDaughter[1].SetRef((AliAODMCParticle*)mcArray->At(label[0]));
         fDaughter[1].SetRefMC((AliAODMCParticle*)mcArray->At(label[0]));
         fDaughter[1].SetGood();
      } else {
         fDaughter[0].SetBad();
         fDaughter[1].SetBad();
      }

      // reject the pair if the matching was unsuccessful
      if (!fDaughter[0].IsOK() || !fDaughter[1].IsOK()) continue;

      // first, we set the internal AliRsnMother object to
      // the MC particles and then operate the selections on MC
      fPair.SetDaughters(&fDaughter[0], pairDef->GetMass(0), &fDaughter[1], pairDef->GetMass(1));
      FillContainer(cont, &fStepListMC, pairDef, 0);

      // then, if both particles found a good match in the AOD
      // reassociate the AOD tracks to the pair and fill AOD containers
      if (aodIndex[0] < 0 || aodIndex[1] < 0) continue;
      if (pairDefMatch[0] == 0 && pairDefMatch[1] == 1) {
         // 1st track --> 1st member of PairDef
         fDaughter[0].SetRef(aod->GetTrack(aodIndex[0]));
         // 2nd track --> 2nd member of PairDef
         fDaughter[1].SetRef(aod->GetTrack(aodIndex[1]));
         //cout << "****** MATCHING SCHEME 1" << endl;
      } else if ((pairDefMatch[0] == 1 && pairDefMatch[1] == 0)) {
         // 1st track --> 2nd member of PairDef
         fDaughter[0].SetRef(aod->GetTrack(aodIndex[1]));
         // 2nd track --> 1st member of PairDef
         fDaughter[1].SetRef(aod->GetTrack(aodIndex[0]));
         //cout << "****** MATCHING SCHEME 2" << endl;
      }
      //cout << "****** IDs = " << fDaughter[0].GetID() << ' ' << fDaughter[1].GetID() << endl;
      // here we must remember how many steps were already filled
      first = (Int_t)fStepListMC.GetEntries();
      FillContainer(cont, &fStepListESD, pairDef, first);
   }
}

//_____________________________________________________________________________
void AliRsnAnalysisEffSE::FillContainer(AliCFContainer *cont, const TObjArray *stepList, AliRsnPairDef *pd, Int_t firstOutStep)
{
//
// Fill the containers
//

   Int_t  iaxis, nAxes  = fAxisList.GetEntries();
   Int_t  istep, nSteps = stepList->GetEntries();
   Bool_t computeOK;

   // set daughters to pair
   fPair.SetDaughters(&fDaughter[0], pd->GetMass(0), &fDaughter[1], pd->GetMass(1));

   // compute values for all axes
   for (iaxis = 0; iaxis < nAxes; iaxis++) {
      AliRsnValue *fcnAxis = (AliRsnValue*)fAxisList.At(iaxis);
      fVar[iaxis] = -1E10;
      switch (fcnAxis->GetTargetType()) {
         case AliRsnTarget::kMother:
            fcnAxis->SetSupportObject(pd);
            computeOK = fcnAxis->Eval(&fPair);
            break;
         case AliRsnTarget::kEvent:
            computeOK = fcnAxis->Eval(&fRsnEvent);
            break;
         default:
            AliError(Form("Allowed targets are mothers and events; cannot use axis '%s' which has target '%s'", fcnAxis->GetName(), fcnAxis->GetTargetTypeName()));
            computeOK = kFALSE;
      }
      if (computeOK) fVar[iaxis] = ((Float_t)fcnAxis->GetComputedValue());
   }

   // fill all steps
   for (istep = 0; istep < nSteps; istep++) {
      AliRsnCutManager *cutMgr = (AliRsnCutManager*)stepList->At(istep);
      AliRsnTarget::SwitchToFirst();
      if (!cutMgr->PassCommonDaughterCuts(&fDaughter[0])) break;
      if (!cutMgr->PassCommonDaughterCuts(&fDaughter[1])) break;
      if (!cutMgr->PassDaughter1Cuts(&fDaughter[0])) break;
      if (!cutMgr->PassDaughter2Cuts(&fDaughter[1])) break;
      if (!cutMgr->PassMotherCuts(&fPair)) break;
      //cout << "**************************************** FILLING STEP " << istep << endl;
      cont->Fill(fVar.GetArray(), istep + firstOutStep);
   }
}

//_____________________________________________________________________________
void AliRsnAnalysisEffSE::RsnTerminate(Option_t*)
{
//
// Termination.
// Could be added some monitor histograms here.
//

   AliDebug(AliLog::kDebug + 2, "<-");
   AliDebug(AliLog::kDebug + 2, "->");
}

//_____________________________________________________________________________
void AliRsnAnalysisEffSE::AddPairDef(AliRsnPairDef* pairDef)
{
//
//  Adds pair definition
//
   fPairDefList.AddLast(pairDef);
}

//_____________________________________________________________________________
Int_t AliRsnAnalysisEffSE::FindESDtrack(Int_t label, AliESDEvent *esd, Bool_t rejectFakes)
{
//
// Finds in the ESD a track whose label corresponds to that in argument.
// When global tracks are enabled, tries first to find a global track
// satisfying that requirement.
// If no global tracks are found, if ITS-SA are enable, tries to search among them
// otherwise return a negative number.
// If global tracks are disabled, search only among ITS SA
//

   Int_t   i = 0;
   Int_t   ntracks = esd->GetNumberOfTracks();
   ULong_t status;
   Bool_t  isTPC;
   Bool_t  isITSSA;

   // loop for global tracks
   if (fUseGlobal) {
      for (i = 0; i < ntracks; i++) {
         AliESDtrack *track = esd->GetTrack(i);
         status  = (ULong_t)track->GetStatus();
         isTPC   = ((status & AliESDtrack::kTPCin)  != 0);
         if (!isTPC) continue;

         // check that label match
         if (TMath::Abs(track->GetLabel()) != label) continue;

         // if required, reject fakes
         if (rejectFakes && track->GetLabel() < 0) continue;

         // if all checks are passed and we are searching among global
         // this means that thie track is a global one with the right label
         // then, the return value is set to this, and returned
         return i;
      }
   }

   // loop for ITS-SA tracks (this happens only if no global tracks were found
   // or searching among globals is disabled)
   if (fUseITSSA) {
      for (i = 0; i < ntracks; i++) {
         AliESDtrack *track = esd->GetTrack(i);
         status  = (ULong_t)track->GetStatus();
         isITSSA = ((status & AliESDtrack::kTPCin)  == 0 && (status & AliESDtrack::kITSrefit) != 0 && (status & AliESDtrack::kITSpureSA) == 0 && (status & AliESDtrack::kITSpid) != 0);
         if (!isITSSA) continue;

         // check that label match
         if (TMath::Abs(track->GetLabel()) != label) continue;

         // if required, reject fakes
         if (rejectFakes && track->GetLabel() < 0) continue;

         // if all checks are passed and we are searching among global
         // this means that thie track is a global one with the right label
         // then, the return value is set to this, and returned
         return i;
      }
   }

   // if we reach this point, no match were found
   return -1;
}

//_____________________________________________________________________________
TArrayI AliRsnAnalysisEffSE::FindESDtracks(Int_t label, AliESDEvent *esd)
{
//
// Finds in the ESD a track whose label corresponds to that in argument.
// When global tracks are enabled, tries first to find a global track
// satisfying that requirement.
// If no global tracks are found, if ITS-SA are enable, tries to search among them
// otherwise return a negative number.
// If global tracks are disabled, search only among ITS SA
//

   Int_t   i = 0;
   Int_t   ntracks = esd->GetNumberOfTracks();
   ULong_t status;
   Bool_t  isTPC;
   Bool_t  isITSSA;
   TArrayI array(100);
   Int_t   nfound = 0;

   // loop for global tracks
   if (fUseGlobal) {
      for (i = 0; i < ntracks; i++) {
         AliESDtrack *track = esd->GetTrack(i);
         status  = (ULong_t)track->GetStatus();
         isTPC   = ((status & AliESDtrack::kTPCin)  != 0);
         if (!isTPC) continue;

         // check that label match
         if (TMath::Abs(track->GetLabel()) != label) continue;

         array[nfound++] = i;
      }
   }

   // loop for ITS-SA tracks (this happens only if no global tracks were found
   // or searching among globals is disabled)
   if (fUseITSSA) {
      for (i = 0; i < ntracks; i++) {
         AliESDtrack *track = esd->GetTrack(i);
         status  = (ULong_t)track->GetStatus();
         isITSSA = ((status & AliESDtrack::kTPCin)  == 0 && (status & AliESDtrack::kITSrefit) != 0 && (status & AliESDtrack::kITSpureSA) == 0 && (status & AliESDtrack::kITSpid) != 0);
         if (!isITSSA) continue;

         // check that label match
         if (TMath::Abs(track->GetLabel()) != label) continue;

         array[nfound++] = i;
      }
   }

   array.Set(nfound);
   return array;
}

//_____________________________________________________________________________
Int_t AliRsnAnalysisEffSE::FindAODtrack(Int_t label, AliAODEvent *aod, Bool_t rejectFakes)
{
//
// Finds in the ESD a track whose label corresponds to that in argument.
// When global tracks are enabled, tries first to find a global track
// satisfying that requirement.
// If no global tracks are found, if ITS-SA are enable, tries to search among them
// otherwise return a negative number.
// If global tracks are disabled, search only among ITS SA
//

   Int_t   i = 0;
   Int_t   ntracks = aod->GetNumberOfTracks();
   ULong_t status;
   Bool_t  isTPC;
   Bool_t  isITSSA;

   // loop for global tracks
   if (fUseGlobal) {
      for (i = 0; i < ntracks; i++) {
         AliAODTrack *track = aod->GetTrack(i);
         status  = (ULong_t)track->GetStatus();
         isTPC   = ((status & AliESDtrack::kTPCin)  != 0);
         if (!isTPC) continue;

         // check that label match
         if (TMath::Abs(track->GetLabel()) != label) continue;

         // if required, reject fakes
         if (rejectFakes && track->GetLabel() < 0) continue;

         // if all checks are passed and we are searching among global
         // this means that thie track is a global one with the right label
         // then, the return value is set to this, and returned
         return i;
      }
   }

   // loop for ITS-SA tracks (this happens only if no global tracks were found
   // or searching among globals is disabled)
   if (fUseITSSA) {
      for (i = 0; i < ntracks; i++) {
         AliAODTrack *track = aod->GetTrack(i);
         status  = (ULong_t)track->GetStatus();
         isITSSA = ((status & AliESDtrack::kTPCin)  == 0 && (status & AliESDtrack::kITSrefit) != 0 && (status & AliESDtrack::kITSpureSA) == 0 && (status & AliESDtrack::kITSpid) != 0);
         if (!isITSSA) continue;

         // check that label match
         if (TMath::Abs(track->GetLabel()) != label) continue;

         // if required, reject fakes
         if (rejectFakes && track->GetLabel() < 0) continue;

         // if all checks are passed and we are searching among global
         // this means that thie track is a global one with the right label
         // then, the return value is set to this, and returned
         return i;
      }
   }

   // if we reach this point, no match were found
   return -1;
}

//_____________________________________________________________________________
TArrayI AliRsnAnalysisEffSE::FindAODtracks(Int_t label, AliAODEvent *aod)
{
//
// Finds in the ESD a track whose label corresponds to that in argument.
// When global tracks are enabled, tries first to find a global track
// satisfying that requirement.
// If no global tracks are found, if ITS-SA are enable, tries to search among them
// otherwise return a negative number.
// If global tracks are disabled, search only among ITS SA
//

   Int_t   i = 0;
   Int_t   ntracks = aod->GetNumberOfTracks();
   ULong_t status;
   Bool_t  isTPC;
   Bool_t  isITSSA;
   TArrayI array(100);
   Int_t   nfound = 0;

   // loop for global tracks
   if (fUseGlobal) {
      for (i = 0; i < ntracks; i++) {
         AliAODTrack *track = aod->GetTrack(i);
         status  = (ULong_t)track->GetStatus();
         isTPC   = ((status & AliESDtrack::kTPCin)  != 0);
         if (!isTPC) continue;

         // check that label match
         if (TMath::Abs(track->GetLabel()) != label) continue;

         array[nfound++] = i;
      }
   }

   // loop for ITS-SA tracks (this happens only if no global tracks were found
   // or searching among globals is disabled)
   if (fUseITSSA) {
      for (i = 0; i < ntracks; i++) {
         AliAODTrack *track = aod->GetTrack(i);
         status  = (ULong_t)track->GetStatus();
         isITSSA = ((status & AliESDtrack::kTPCin)  == 0 && (status & AliESDtrack::kITSrefit) != 0 && (status & AliESDtrack::kITSpureSA) == 0 && (status & AliESDtrack::kITSpid) != 0);
         if (!isITSSA) continue;

         // check that label match
         if (TMath::Abs(track->GetLabel()) != label) continue;

         array[nfound++] = i;
      }
   }

   array.Set(nfound);
   return array;
}

//______________________________________________________________________________
Bool_t AliRsnAnalysisEffSE::EventProcess()
{
//
// Customized event pre-processing.
// First checks if the current event passes all cuts,
// and if it does, updates the informations and then
// call the operations which are already defined in the
// omonyme function in mother class
//

   // initially, an event is expected to be bad
   fTaskInfo.SetEventUsed(kFALSE);

   // check the event cuts and update the info data accordingly
   // events not passing the cuts must be rejected
   if (!fEventCuts.IsSelected(&fRsnEvent)) {
      fTaskInfo.SetEventUsed(kFALSE);
      return kFALSE;
   }

   // if we reach this point, cuts were passed;
   // then additional operations can be done

   // find leading particle (without any PID/momentum restriction)
   fRsnEvent.SelectLeadingParticle(0);

   // final return value is positive
   // but call the mother class method which updates info object
   fTaskInfo.SetEventUsed(kTRUE);
   return AliRsnVAnalysisTaskSE::EventProcess();
}
