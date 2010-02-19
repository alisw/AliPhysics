//
// Class AliRsnAnalysisME
//
//
// Virtual Class derivated from AliRsnVAnalysisTaskME which will be base class
// for all RSN SE tasks
//
// authors: Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//          Martin Vala (martin.vala@cern.ch)
//

#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliAODEvent.h"
#include "AliMultiEventInputHandler.h"
#include "AliRsnAnalysisME.h"

ClassImp(AliRsnAnalysisME)

//_____________________________________________________________________________
AliRsnAnalysisME::AliRsnAnalysisME(const char *name, Int_t numOfOutputs) :
    AliRsnVAnalysisTaskME(name, numOfOutputs),
    fRsnAnalysisManager(),
    fPIDIndex(0),
    fPIDIndexMix(0),
    fEvent(),
    fEventMix()
{
//
// Default constructor
//
  AliDebug(AliLog::kDebug+2, "<-");
  Int_t i = 0;
  for (i = 0;i < fNumberOfOutputs;i++) {
    DefineOutput(i + 2, TList::Class());
  }
  AliDebug(AliLog::kDebug+2,"->");
}

AliRsnAnalysisME::AliRsnAnalysisME(const AliRsnAnalysisME& copy) : AliRsnVAnalysisTaskME(copy),
    fRsnAnalysisManager(copy.fRsnAnalysisManager),
    fPIDIndex(copy.fPIDIndex),
    fPIDIndexMix(copy.fPIDIndexMix),
    fEvent(copy.fEvent),
    fEventMix(copy.fEvent)
{
  AliDebug(AliLog::kDebug+2, "<-");
  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnAnalysisME::RsnUserCreateOutputObjects()
{
//
// Creation of output objects.
// These are created through the utility methods in the analysis manager,
// which produces a list of histograms for each specified set of pairs.
// Each of these lists is added to the main list of this task.
//

  AliDebug(AliLog::kDebug+2, "<-");

  Int_t i;
  for (i = 1; i < kMaxNumberOfOutputs + 1; i++) {
    if (i <= fNumberOfOutputs + 1) OpenFile(i/* + 1*/);
    fOutList[i] = new TList();
    fOutList[i]->SetOwner();
  }

  for (i = 0; i < fNumberOfOutputs; i++) {
    fRsnAnalysisManager[i].InitAllPairMgrs(fOutList[i+1]);
  }

  AliDebug(AliLog::kDebug+2,"->");
}

void AliRsnAnalysisME::RsnUserExec(Option_t*)
{
//
// Rsn User Exec
//

  fTaskInfo.SetEventUsed(kFALSE);

  AliDebug(AliLog::kDebug+2, "<-");
  if (!CheckAndPrintEvents()) return;

  DoMixing(GetEvent(0));


  // if cuts are passed or not cuts were defined,
  // update the task info...
  fTaskInfo.SetEventUsed(kTRUE);

  // the virtual class has already sorted tracks in the PID index
  // so we need here just to call the execution of analysis
  for (Int_t i = 0; i < fNumberOfOutputs; i++) {
//     fRsnAnalysisManager[i].ProcessAllPairMgrs(&fRsnPIDIndex, &fRsnEvent);
    PostData(i + 2, fOutList[i+1]);
  }

  AliDebug(AliLog::kDebug+2,"->");
}


//_____________________________________________________________________________
void AliRsnAnalysisME::DoMixing(AliVEvent* ev)
{
//
// Do Mixing
//

  Int_t nEvents = fInputHandler->GetBufferSize();
  fESDEvent = dynamic_cast<AliESDEvent*>(ev);
  fAODEvent = dynamic_cast<AliAODEvent*>(ev);

  if (fESDEvent) {
    AliESDEvent **esdEvent = new AliESDEvent*[nEvents];
    for (Int_t i = 0; i < nEvents; i++) {
      esdEvent[i] = dynamic_cast<AliESDEvent*>(GetEvent(i));
      if (!esdEvent[i]) {
        AliWarning(Form("Null ESD event in index %d", i));
        continue;
      }
      if (i > 0)
        DoESDMixing(esdEvent[0], esdEvent[i]);
    }
  } else if (fAODEvent) {
    AliAODEvent **aodEvent = new AliAODEvent*[nEvents];
    for (Int_t i = 0; i < nEvents; i++) {
      aodEvent[i] = dynamic_cast<AliAODEvent*>(GetEvent(i));
      if (!aodEvent[i]) {
        AliWarning(Form("Null AOD event in index %d", i));
        continue;
      }
      if (i > 0)
        DoAODMixing(aodEvent[0], aodEvent[i]);
    }
  }

}


//_____________________________________________________________________________
void AliRsnAnalysisME::DoAODMixing(AliAODEvent* aod1, AliAODEvent* aod2)
{
//
// mixing of two aod events
//

  // assign events
  fEvent.SetRef(aod1);
  fEventMix.SetRef(aod2);
  if (fEvent.GetMultiplicity() < 2) return;
  if (fEventMix.GetMultiplicity() < 2) return;

  // sort tracks w.r. to PID
  fPIDIndex.ResetAll(fEvent.GetMultiplicity());
  fEvent.SetPriorProbability(fPrior);
  fPIDIndex.FillFromEvent(&fEvent);
  fPIDIndex.SetCorrectIndexSize();

  fPIDIndexMix.ResetAll(fEventMix.GetMultiplicity());
  fEventMix.SetPriorProbability(fPrior);
  fPIDIndexMix.FillFromEvent(&fEventMix);
  fPIDIndexMix.SetCorrectIndexSize();

  for (Int_t i = 0; i < fNumberOfOutputs; i++) {
    fRsnAnalysisManager[i].ProcessAllPairMgrs(&fPIDIndex, &fEvent, &fPIDIndexMix, &fEventMix);
    PostData(i + 2, fOutList[i+1]);
  }
  AliDebug(AliLog::kDebug, Form("AOD tracks %d", aod1->GetNumberOfTracks()));
  AliDebug(AliLog::kDebug, Form("AOD tracks %d", aod2->GetNumberOfTracks()));

}


//_____________________________________________________________________________
void AliRsnAnalysisME::DoESDMixing(AliESDEvent* esd1, AliESDEvent* esd2)
{
//
// mixing of two esd events
//

  AliWarning(Form("ESD mixing not supported yet !!! (%p,%p)", esd1, esd2));
  return;
}



//_____________________________________________________________________________
void AliRsnAnalysisME::RsnTerminate(Option_t*)
{
//
// Rsn Terminate
//

  AliDebug(AliLog::kDebug+2, "<-");
  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
AliRsnAnalysisManager* AliRsnAnalysisME::GetAnalysisManager(Int_t index, TString name)
{
//
// Recovery the analysis manager
//

  if (!name.IsNull()) {
    SetAnalysisManagerName(name.Data(), index);
  }

  return &fRsnAnalysisManager[index];
}


//_____________________________________________________________________________
void AliRsnAnalysisME::SetPriorProbability(AliPID::EParticleType type, Double_t p)
{
//
// Sets the prior probability for Realistic PID, for a
// given particle species.
//

  if (type >= 0 && type < (Int_t)AliPID::kSPECIES) {
    fPrior[type] = p;
  }

}

//_____________________________________________________________________________
void AliRsnAnalysisME::DumpPriors()
{
  //
  // Print all prior probabilities
  //

  Int_t i;
  for (i = 0; i < AliPID::kSPECIES; i++) {
    AliInfo(Form("Prior probability for %10s = %3.5f", AliPID::ParticleName((AliPID::EParticleType)i), fPrior[i]));
  }
}

//_____________________________________________________________________________
void AliRsnAnalysisME::GetPriorProbability(Double_t *out) const
{
//
// Gets all prior probabilities to out
//

  Int_t i;
  for (i = 0;i < AliPID::kSPECIES;i++) {
    out[i] = fPrior[i];
  }
}
