//
// Class AliRsnAnalysisME
//
// TODO
//
// authors: Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//          Martin Vala (martin.vala@cern.ch)
//

#include "AliRsnAnalysisSE.h"

ClassImp(AliRsnAnalysisSE)

//_____________________________________________________________________________
AliRsnAnalysisSE::AliRsnAnalysisSE(const char *name) :
    AliRsnVAnalysisTaskSE(name),
    fRsnAnalysisManager(),
    fPIDIndex(0),
    fEvent(),
    fESDCuts(0)
{
//
// Default constructor
//
  AliDebug(AliLog::kDebug+2,"<-");
  AliDebug(AliLog::kDebug+2,"->");
}

AliRsnAnalysisSE::AliRsnAnalysisSE(const AliRsnAnalysisSE& copy) : AliRsnVAnalysisTaskSE(copy),
    fRsnAnalysisManager(copy.fRsnAnalysisManager),
    fPIDIndex(copy.fPIDIndex),
    fEvent(copy.fEvent),
    fESDCuts(copy.fESDCuts)
{
  AliDebug(AliLog::kDebug+2,"<-");
  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnAnalysisSE::RsnUserCreateOutputObjects()
{
  AliDebug(AliLog::kDebug+2,"<-");

  fOutList->Add(fRsnAnalysisManager.InitAllPairMgrs());

  AliDebug(AliLog::kDebug+2,"->");
}

void AliRsnAnalysisSE::RsnUserExec(Option_t* )
{
  AliDebug(AliLog::kDebug+2,"<-");
  if (fESDEvent) {
    AliDebug(AliLog::kDebug+1,Form("fESDEvent if %p",fESDEvent));
    AliDebug(AliLog::kDebug,Form("ESD tracks %d",fESDEvent->GetNumberOfTracks()));
  }
  if (fMCEvent) {
    AliDebug(AliLog::kDebug+1,Form("fMCEvent if %p",fMCEvent));
    AliDebug(AliLog::kDebug,Form("MC tracks %d",fMCEvent->GetNumberOfTracks()));
  }
  if (fAODEventIn) {
    AliDebug(AliLog::kDebug+1,Form("fAODEventIn if %p",fAODEventIn));
    AliDebug(AliLog::kDebug,Form("AOD(in) tracks %d",fAODEventIn->GetNumberOfTracks()));
  }

  if (fAODEventOut) {
    AliDebug(AliLog::kDebug+1,Form("fAODEventOut if %p",fAODEventOut));
    AliDebug(AliLog::kDebug,Form("AOD(out) tracks %d",fAODEventOut->GetNumberOfTracks()));
  }

  /*
  // assign event
  if (fAODEventOut)
    fEvent.SetRef(fAODEventOut);
  else if (fESDEvent)
    fEvent.SetRef(fESDEvent, fMCEvent);
  else
    return;
  */
  if (fESDEvent)
    fEvent.SetRef(fESDEvent, fMCEvent);
  else if (fAODEventOut)
    fEvent.SetRef(fAODEventOut);
  else
    return;
  if (fEvent.GetMultiplicity()<2) return;

  // sort tracks w.r. to PID
  fPIDIndex.ResetAll(fEvent.GetMultiplicity());
  fPIDIndex.SetPriorProbability(fPrior);
  fPIDIndex.FillFromEvent(&fEvent, fESDCuts);
  fPIDIndex.SetCorrectIndexSize();

  fRsnAnalysisManager.ProcessAllPairMgrs(&fPIDIndex, &fEvent);

  AliDebug(AliLog::kDebug+2,"->");
}


//_____________________________________________________________________________
void AliRsnAnalysisSE::RsnTerminate(Option_t* )
{
  AliDebug(AliLog::kDebug+2,"<-");
  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
AliRsnAnalysisManager* AliRsnAnalysisSE::GetAnalysisManager(TString name)
{
  if (!name.IsNull()) {
    SetAnalysisManagerName(name.Data());
  }

  return &fRsnAnalysisManager;
}


//_____________________________________________________________________________
void AliRsnAnalysisSE::SetPriorProbability(AliPID::EParticleType type, Double_t p)
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
void AliRsnAnalysisSE::DumpPriors()
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
void AliRsnAnalysisSE::GetPriorProbability(Double_t *out)
{

  Int_t i;
  for (i=0;i<AliPID::kSPECIES;i++) {
    out[i] = fPrior[i];
  }
}
