//
// Class AliRsnAnalysisME
//
// TODO
// TODO
//
// authors: Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//          Martin Vala (martin.vala@cern.ch)
//

#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliAODEvent.h"

#include "AliRsnAnalysisME.h"

ClassImp(AliRsnAnalysisME)

//_____________________________________________________________________________
AliRsnAnalysisME::AliRsnAnalysisME(const char *name) :
    AliRsnVAnalysisTaskME(name),
    fRsnAnalysisManager(),
    fPIDIndex(0),
    fPIDIndexMix(0),
    fEvent(),
    fEventMix()
{
//
// Default constructor
//
  AliDebug(AliLog::kDebug+2,"<-");
  AliDebug(AliLog::kDebug+2,"->");
}

AliRsnAnalysisME::AliRsnAnalysisME(const AliRsnAnalysisME& copy) : AliRsnVAnalysisTaskME(copy),
    fRsnAnalysisManager(copy.fRsnAnalysisManager),
    fPIDIndex(copy.fPIDIndex),
    fPIDIndexMix(copy.fPIDIndexMix),
    fEvent(copy.fEvent),
    fEventMix(copy.fEvent)
{
  AliDebug(AliLog::kDebug+2,"<-");
  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnAnalysisME::RsnUserCreateOutputObjects()
{
//
// TODO
//
  AliLog::SetClassDebugLevel(GetName(), fLogType);
  AliDebug(AliLog::kDebug+2,"<-");
//     AliRsnVAnalysisTaskME::UserCreateOutputObjects();

//       fRsnEvent = new AliRsnEvent();
//       fRsnEvent->SetName("rsnEvents");
//       fRsnEvent->Init();
//       AddAODBranch("AliRsnEvent", &fRsnEvent);

//     fOutList = new TList();
//
//     fOutList->Add(fTaskInfo.GenerateInfoList());

//   fOutList->Add();
  fRsnAnalysisManager.InitAllPairMgrs(fOutList);
//     fRsnAnalysisManager.Print();

  AliDebug(AliLog::kDebug+2,"->");
}

void AliRsnAnalysisME::RsnUserExec(Option_t*)
{
//
// TODO
//

  AliDebug(AliLog::kDebug+2,"<-");
  if (fESDEvent) {
    AliDebug(AliLog::kDebug+1,Form("fESDEvent if %p",fESDEvent));
    AliDebug(AliLog::kDebug,Form("ESD tracks %d",fESDEvent->GetNumberOfTracks()));
  }
  if (fMCEvent) {
    AliDebug(AliLog::kDebug+1,Form("fMCEvent if %p",fMCEvent));
    AliDebug(AliLog::kDebug,Form("MC tracks %d",fMCEvent->GetNumberOfTracks()));
  }
  if (fAODEvent) {
    AliDebug(AliLog::kDebug+1,Form("fAODEvent if %p",fAODEvent));
    AliDebug(AliLog::kDebug,Form("AOD tracks %d",fAODEvent->GetNumberOfTracks()));
  }

  AliAODEvent* aod1 = (AliAODEvent*)GetEvent(0);
  AliAODEvent* aod2 = (AliAODEvent*)GetEvent(1);

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

  fRsnAnalysisManager.ProcessAllPairMgrs(&fPIDIndex, &fEvent, &fPIDIndexMix, &fEventMix);

  AliDebug(AliLog::kDebug,Form("AOD tracks %d",aod1->GetNumberOfTracks()));
  AliDebug(AliLog::kDebug,Form("AOD tracks %d",aod2->GetNumberOfTracks()));
  AliDebug(AliLog::kDebug+2,"->");
}


//_____________________________________________________________________________
void AliRsnAnalysisME::RsnTerminate(Option_t*)
{
//
// TODO
//

  AliDebug(AliLog::kDebug+2,"<-");
  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
AliRsnAnalysisManager* AliRsnAnalysisME::GetAnalysisManager(TString name)
{
//
// Gets all prior probabilities to out
//

  if (!name.IsNull()) {
    SetAnalysisManagerName(name.Data());
  }
  return &fRsnAnalysisManager;
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
  for (i=0;i<AliPID::kSPECIES;i++) {
    out[i] = fPrior[i];
  }
}
