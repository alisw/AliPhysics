//
// Class AliRsnAnalysisTrackEffSE
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
#include "AliRsnCutSet.h"
#include "AliRsnFunctionAxis.h"
#include "AliRsnAnalysisTrackEffSE.h"
#include "AliRsnCutSet.h"

ClassImp(AliRsnAnalysisTrackEffSE)

//_____________________________________________________________________________
AliRsnAnalysisTrackEffSE::AliRsnAnalysisTrackEffSE(const char *name) :
  AliRsnVAnalysisTaskSE(name),
  fEventCuts(0x0),
  fStepListMC(0),
  fStepListESD(0),
  fAxisList(0),
  fVar(0),
  fDaughter()
{
//
// Default constructor.
//

  AliDebug(AliLog::kDebug+2,"<-");

  DefineOutput(2, TList::Class());

  Int_t i;
  for (i = 0; i <= AliPID::kSPECIES; i++) fContainer[i] = 0x0;

  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
AliRsnAnalysisTrackEffSE::AliRsnAnalysisTrackEffSE(const AliRsnAnalysisTrackEffSE& copy) :
  AliRsnVAnalysisTaskSE(copy),
  fEventCuts(copy.fEventCuts),
  fStepListMC(copy.fStepListMC),
  fStepListESD(copy.fStepListESD),
  fAxisList(copy.fAxisList),
  fVar(0),
  fDaughter()
{
//
// Copy constrtuctor
//
}

//_____________________________________________________________________________
void AliRsnAnalysisTrackEffSE::RsnUserCreateOutputObjects()
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

  // initialize output list
  OpenFile(2);
  fOutList[1] = new TList();
  fOutList[1]->SetOwner();

  // create the containers
  Int_t i;
  for (i = 0; i <= AliPID::kSPECIES; i++)
  {
    if (i < AliPID::kSPECIES) fContainer[i] = new AliCFContainer(Form("%s", AliPID::ParticleName((AliPID::EParticleType)i)), "", nSteps, nAxes, nBins);
    else fContainer[i] = new AliCFContainer("all", "", nSteps, nAxes, nBins);
    // set the bin limits for each axis
    for (iaxis = 0; iaxis < nAxes; iaxis++) fContainer[i]->SetBinLimits(iaxis, binLimits[iaxis].GetArray());
    // add the container to output list
    fOutList[1]->Add(fContainer[i]);
  }

  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnAnalysisTrackEffSE::RsnUserExec(Option_t*)
{
//
// Execution of the analysis task.
// Recovers the input event and processes it with all included pair objects.
// In this case, we NEED to have ESD and MC, otherwise we cannod do anything.
//

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
  ProcessEventMC();
  ProcessEventESD();

  // Post the data
  PostData(2, fOutList[1]);

  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnAnalysisTrackEffSE::ProcessEventMC()
{
//
// Process current event with the definitions of the specified step in MC list
// and store results in the container slot defined in second argument
//

  AliStack      *stack = fMCEvent->Stack();
  AliMCParticle *part;

  // other utility variables
  Int_t ipart;

  // in this case, we first take the resonance from MC
  // and then we find its daughters and compute cuts on them
  for (ipart = 0; ipart < fMCEvent->GetNumberOfTracks(); ipart++)
  {
    part = (AliMCParticle*) fMCEvent->GetTrack(ipart);

    fDaughter.SetRef(part);
    fDaughter.SetParticle(part->Particle());
    fDaughter.FindMotherPDG(stack);
    
    FillContainer(&fStepListMC, 0);
  }
}

//_____________________________________________________________________________
void AliRsnAnalysisTrackEffSE::ProcessEventESD()
{
//
// Process current event with the definitions of the specified step in ESD list
// and store results in the container slot defined in second argument
//

  Int_t ic, i, first = (Int_t)fStepListMC.GetEntries();

  // get arrays of all charged tracks
  TArrayI *a[2];
  a[0] = fRsnPIDIndex.GetTracksArray(AliRsnDaughter::kNoPID, '+', AliPID::kUnknown);
  a[1] = fRsnPIDIndex.GetTracksArray(AliRsnDaughter::kNoPID, '-', AliPID::kUnknown);

  // external loop on tracks
  for (ic = 0; ic < 2; ic++)
  {
    for (i = 0; i < a[ic]->GetSize(); i++)
    {
      // connect interface
      fRsnEvent.SetDaughter(fDaughter, a[ic]->At(i));
      if (!fDaughter.IsOK()) continue;
      fDaughter.SetRequiredPID(fDaughter.PerfectPID());

      FillContainer(&fStepListESD, first);
    }
  }
}

//_____________________________________________________________________________
void AliRsnAnalysisTrackEffSE::FillContainer
(const TObjArray *stepList, Int_t firstOutStep)
{
//
// Fill the containers
//

  Int_t ipid;
  Int_t iaxis, nAxes  = fAxisList.GetEntries();
  Int_t istep, nSteps = stepList->GetEntries();

  // compute values for all axes
  for (iaxis = 0; iaxis < nAxes; iaxis++) {
    AliRsnFunctionAxis *fcnAxis = (AliRsnFunctionAxis*)fAxisList.At(iaxis);
    switch (fcnAxis->GetAxisObject())
    {
      case AliRsnFunctionAxis::kParticle:
        fVar[iaxis] = (Double_t)fcnAxis->Eval(&fDaughter);
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
    AliRsnCutSet *cutSet = (AliRsnCutSet*)stepList->At(istep);
    if (!cutSet->IsSelected(AliRsnCut::kParticle, &fDaughter)) break;
    if (stepList == &fStepListESD && !PassedAllCutsMC()) break;

    ipid = (Int_t)fDaughter.PerfectPID();
    if (ipid == (Int_t)AliPID::kUnknown) ipid = (Int_t)AliPID::kSPECIES;
    fContainer[ipid]->Fill(fVar.GetArray(), istep + firstOutStep);
    fContainer[AliPID::kSPECIES]->Fill(fVar.GetArray(), istep + firstOutStep);
  }
}

//_____________________________________________________________________________
Bool_t AliRsnAnalysisTrackEffSE::PassedAllCutsMC()
{
//
// Check if this daughter passes all cuts MC
//

  Int_t istep, nSteps = fStepListMC.GetEntries();

  for (istep = 0; istep < nSteps; istep++) {
    AliRsnCutSet *cutSet = (AliRsnCutSet*)fStepListMC.At(istep);
    if (!cutSet->IsSelected(AliRsnCut::kParticle, &fDaughter)) return kFALSE;
  }

  return kTRUE;
}

//_____________________________________________________________________________
void AliRsnAnalysisTrackEffSE::RsnTerminate(Option_t*)
{
//
// Termination.
// Could be added some monitor histograms here.
//

  AliDebug(AliLog::kDebug+2,"<-");
  AliDebug(AliLog::kDebug+2,"->");
}

