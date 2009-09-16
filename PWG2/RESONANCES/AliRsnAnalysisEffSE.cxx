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
#include "AliMCEvent.h"

#include "AliCFContainer.h"
#include "AliRsnCutMgr.h"
#include "AliRsnFunctionAxis.h"
#include "AliRsnAnalysisEffSE.h"
#include "AliRsnPairDef.h"
#include "AliRsnCutSet.h"

ClassImp(AliRsnAnalysisEffSE)

//_____________________________________________________________________________
AliRsnAnalysisEffSE::AliRsnAnalysisEffSE(const char *name) :
  AliRsnVAnalysisTaskSE(name),
  fEventCuts(0x0),
  fStepListMC(0),
  fStepListESD(0),
  fAxisList(0),
  fPairDefList(0),
  fContainerList(0x0),
  fVar(0),
  fPair()
{
//
// Default constructor.
//

  AliDebug(AliLog::kDebug+2,"<-");

  DefineOutput(2, TList::Class());

  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
AliRsnAnalysisEffSE::AliRsnAnalysisEffSE(const AliRsnAnalysisEffSE& copy) :
  AliRsnVAnalysisTaskSE(copy),
  fEventCuts(copy.fEventCuts),
  fStepListMC(copy.fStepListMC),
  fStepListESD(copy.fStepListESD),
  fAxisList(copy.fAxisList),
  fPairDefList(copy.fPairDefList),
  fContainerList(copy.fContainerList),
  fVar(0),
  fPair()
{
//
// Copy constrtuctor
//
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

  AliDebug(AliLog::kDebug+2,"<-");

  // get number of steps and axes
  Int_t iaxis;
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
  Int_t   *nBins     = new Int_t[nAxes];
  TArrayD *binLimits = new TArrayD[nAxes];
  for (iaxis = 0; iaxis < nAxes; iaxis++) {
    AliRsnFunctionAxis *fcnAxis = (AliRsnFunctionAxis*)fAxisList.At(iaxis);
    binLimits[iaxis] = fcnAxis->GetArray();
    nBins[iaxis] = binLimits[iaxis].GetSize() - 1;
  }

  // create ouput list of containers
  fContainerList = new TList;
  fContainerList->SetOwner();
  fContainerList->SetName(Form("%s_containers", GetName()));

  // initialize output list
  OpenFile(2);
  fOutList[1] = new TList();
  fOutList[1]->SetOwner();

  // create the containers
  Int_t i, nDef = (Int_t)fPairDefList.GetEntries();
  for (i = 0; i < nDef; i++) {
    AliRsnPairDef *def = (AliRsnPairDef*)fPairDefList[i];
    AliCFContainer *cont = new AliCFContainer(Form("%s", def->GetPairName().Data()), "", nSteps, nAxes, nBins);
    // set the bin limits for each axis
    for (iaxis = 0; iaxis < nAxes; iaxis++) cont->SetBinLimits(iaxis, binLimits[iaxis].GetArray());
    // add the container to output list
    fContainerList->Add(cont);
  }

  fOutList[1]->Add(fContainerList);

  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnAnalysisEffSE::RsnUserExec(Option_t*)
{
//
// Execution of the analysis task.
// Recovers the input event and processes it with all included pair objects.
// In this case, we NEED to have ESD and MC, otherwise we cannod do anything.
//ULTIMO UNDO
/*
  AliRsnDaughter trk;
  for (Int_t i = 0; i <= AliPID::kSPECIES; i++)
  {
    cout << AliPID::ParticleName((AliPID::EParticleType)i) << ": " << endl;
    for (Int_t m = 0; m < AliRsnDaughter::kMethods; m++)
    {
      cout << "-- method: " << AliRsnDaughter::MethodName((AliRsnDaughter::EPIDMethod)m) << endl;
      Char_t   sign[2] = {'+', '-'};
      for (Int_t s = 0; s < 2; s++)
      {
        TArrayI *a = fRsnPIDIndex.GetTracksArray((AliRsnDaughter::EPIDMethod)m, sign[s], (AliPID::EParticleType)i);
        Int_t n = a->GetSize();
        for (Int_t j = 0; j < n; j++)
        {
          Int_t k = a->At(j);
          cout << "-- -- track " << Form("%4d ", k) << ": ";
          fRsnEvent.SetDaughter(trk, k);
          cout << "charge = " << (trk.IsPos() ? "+ " : (trk.IsNeg() ? "- " : "0 "));
          cout << "truePID = " << Form("%10s ", AliPID::ParticleName(trk.PerfectPID()));
          cout << "realPID = " << Form("%10s ", AliPID::ParticleName(trk.RealisticPID()));
          cout << endl;
          cout << "-- -- -- weights (computed): ";
          for (Int_t q = 0; q < AliPID::kSPECIES; q++)
            cout << Form("%15.12f", trk.ComputedWeights()[q]) << ' ';
          cout << endl;
          cout << "-- -- -- weights (original): ";
          for (Int_t q = 0; q < AliPID::kSPECIES; q++)
            cout << Form("%15.12f", trk.GetRef()->PID()[q]) << ' ';
          cout << endl;
        }
      }
    }
  } return;
  */

  AliDebug(AliLog::kDebug+2,"<-");
  if (!fESDEvent || !fMCEvent) {
    AliError("This task can process only ESD + MC events");
    return;
  }
  fRsnEvent.SetRef(fESDEvent, fMCEvent);

  // if general event cuts are added to the task (recommended)
  // they are checked here on the RSN event interface and,
  // if the event does not pass them, it is skipped and ProcessInfo
  // is updated accordingly
  if (fEventCuts) {
    if (!fEventCuts->IsSelected(AliRsnCut::kEvent, &fRsnEvent)) {
      fTaskInfo.SetEventUsed(kFALSE);
      return;
    }
  }

  // if cuts are passed or not cuts were defined,
  // update the task info before processing the event
  fTaskInfo.SetEventUsed(kTRUE);

  // process first MC steps and then ESD steps
  AliRsnPairDef *pairDef = 0;
  TObjArrayIter iter(&fPairDefList);
  while ( (pairDef = (AliRsnPairDef*)iter.Next()) )
  {
    ProcessEventMC(pairDef);
    ProcessEventESD(pairDef);
  }

  // Post the data
  PostData(2, fOutList[1]);

  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnAnalysisEffSE::ProcessEventMC(AliRsnPairDef *pairDef)
{
//
// Process current event with the definitions of the specified step in MC list
// and store results in the container slot defined in second argument
//

  AliStack      *stack = fMCEvent->Stack();
  AliMCParticle *mother, *daughter;

  if (!pairDef) return;
  AliCFContainer *cont = (AliCFContainer*)fContainerList->FindObject(pairDef->GetPairName().Data());
  if (!cont) return;

  // other utility variables
  Int_t i[2], j, ipart;

  // in this case, we first take the resonance from MC
  // and then we find its daughters and compute cuts on them
  for (ipart = 0; ipart < stack->GetNprimary(); ipart++) {
    mother = (AliMCParticle*) fMCEvent->GetTrack(ipart);
    if (mother->Particle()->GetNDaughters() != 2) continue;

    i[0] = mother->Particle()->GetFirstDaughter();
    i[1] = mother->Particle()->GetLastDaughter();

    for (j = 0; j < 2; j++) {
      daughter = (AliMCParticle*) fMCEvent->GetTrack(i[j]);
      fDaughter[j].SetRef(daughter);
      fDaughter[j].SetParticle(daughter->Particle());
      fDaughter[j].FindMotherPDG(stack);
    }

    if (fDaughter[0].ChargeC() != pairDef->GetCharge(0)) continue;
    if (fDaughter[1].ChargeC() != pairDef->GetCharge(1)) continue;
    if (fDaughter[0].PerfectPID() != pairDef->GetType(0)) continue;
    if (fDaughter[1].PerfectPID() != pairDef->GetType(1)) continue;

    fPair.SetPair(&fDaughter[0], &fDaughter[1]);

    // create pair
    FillContainer(cont, &fStepListMC, pairDef, 0);
  }
}

//_____________________________________________________________________________
void AliRsnAnalysisEffSE::ProcessEventESD(AliRsnPairDef *pairDef)
{
//
// Process current event with the definitions of the specified step in ESD list
// and store results in the container slot defined in second argument
//

  Int_t i0, i1, first = (Int_t)fStepListMC.GetEntries();

  if (!pairDef) return;
  AliCFContainer *cont = (AliCFContainer*)fContainerList->FindObject(pairDef->GetPairName().Data());
  if (!cont) return;

  // get arrays of all charged tracks
  TArrayI *a0 = fRsnPIDIndex.GetTracksArray(AliRsnDaughter::kPerfect, pairDef->GetCharge(0), pairDef->GetType(0));
  TArrayI *a1 = fRsnPIDIndex.GetTracksArray(AliRsnDaughter::kPerfect, pairDef->GetCharge(1), pairDef->GetType(1));

  // external loop on tracks
  for (i0 = 0; i0 < a0->GetSize(); i0++) {
    // connect interface
    fRsnEvent.SetDaughter(fDaughter[0], a0->At(i0));
    if (!fDaughter[0].IsOK()) continue;
    fDaughter[0].SetRequiredPID(pairDef->GetType(0));

    // internal loop on tracks
    for (i1 = 0; i1 < a1->GetSize(); i1++) {
      // connect interface
      fRsnEvent.SetDaughter(fDaughter[1], a1->At(i1));
      if (!fDaughter[1].IsOK()) continue;
      fDaughter[1].SetRequiredPID(pairDef->GetType(1));
      // build pair
      fPair.SetPair(&fDaughter[0], &fDaughter[1]);
      if (TMath::Abs(fPair.CommonMother()) != pairDef->GetMotherPDG()) continue;
      // fill containers
      FillContainer(cont, &fStepListESD, pairDef, first);
    }
  }
}

//_____________________________________________________________________________
void AliRsnAnalysisEffSE::FillContainer(AliCFContainer *cont, const TObjArray *stepList, AliRsnPairDef *pd, Int_t firstOutStep)
{
//
// Fill the containers
//

  Int_t iaxis, nAxes  = fAxisList.GetEntries();
  Int_t istep, nSteps = stepList->GetEntries();

  // compute values for all axes
  for (iaxis = 0; iaxis < nAxes; iaxis++) {
    AliRsnFunctionAxis *fcnAxis = (AliRsnFunctionAxis*)fAxisList.At(iaxis);
    switch (fcnAxis->GetAxisObject()) {
    case AliRsnFunctionAxis::kPair:
      fVar[iaxis] = (Double_t)fcnAxis->Eval(&fPair, pd);
      break;
    case AliRsnFunctionAxis::kEvent:
      fVar[iaxis] = (Double_t)fcnAxis->Eval(&fRsnEvent);
      break;
    default:
      fVar[iaxis] = 0.0;
    }
  }

  // fill all steps
  for (istep = 0; istep < nSteps; istep++) {
    AliRsnCutMgr *cutMgr = (AliRsnCutMgr*)stepList->At(istep);
    if (!cutMgr->IsSelected(AliRsnCut::kParticle, &fDaughter[0])) break;
    if (!cutMgr->IsSelected(AliRsnCut::kParticle, &fDaughter[1])) break;
    if (!cutMgr->IsSelected(AliRsnCut::kPair,     &fPair)) break;
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

  AliDebug(AliLog::kDebug+2,"<-");
  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnAnalysisEffSE::AddPairDef(AliRsnPairDef* pairDef)
{
//
//  Adds pair definition
//
  fPairDefList.AddLast(pairDef);
}
