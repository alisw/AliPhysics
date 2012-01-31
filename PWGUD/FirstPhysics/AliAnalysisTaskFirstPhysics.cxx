
#include "AliAnalysisTaskFirstPhysics.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TList.h"
#include "TString.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliESDHeader.h"
#include "AliMCEvent.h"
#include "AliMultiplicity.h"

#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliStack.h"

ClassImp(AliAnalysisTaskFirstPhysics)

//________________________________________________________________________
AliAnalysisTaskFirstPhysics::AliAnalysisTaskFirstPhysics(const char *name) :
  AliAnalysisTaskSE(name),
  fESD(0),
  fMCEvent(0),
  fOutput(0),
  fbReadMC(0),
  fMCProcessType(kProcND),
  fTrigger(0),
  fCutTrackPtMin(0.15),
  fCutTrackPtMax(100),
  fCutEta(0.8),
  fCutVertexZ(10)
{
  DefineOutput(1, TList::Class());
  for (Int_t i = 0; i < knTrackCuts; i ++) {
    fTrackCuts[i] = 0;
  }
}

//________________________________________________________________________
AliAnalysisTaskFirstPhysics::~AliAnalysisTaskFirstPhysics()
{
  // Destructor. Clean-up the output list, but not the histograms that are put inside
  // (the list is owner and will clean-up these histograms). Protect in PROOF case.
  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutput;
  }
  for (Int_t i = 0; i < knTrackCuts; i ++) {
    delete fTrackCuts[i];
    fTrackCuts[i] = 0;
  }
  if (fTrigger) {
    delete fTrigger;
  }
}

//________________________________________________________________________
void AliAnalysisTaskFirstPhysics::PrepareOutputList()
{
  if (fOutput) {
    AliError("fOutput already initialised.");
    return;
  }
  fOutput = new TList();
  fOutput->SetOwner();
  TH1::SetDefaultSumw2(kTRUE);
  if (fTrigger) {
    AliError("fTrigger is already initialised.");
    return;
  }
  fTrigger = new AliTriggerAnalysis;
}

//________________________________________________________________________
void AliAnalysisTaskFirstPhysics::PrepareDefaultTrackCuts()
{
  // quality cut on ITS+TPC tracks
  fTrackCuts[kTrackCutQGlo] = new AliESDtrackCuts();
  // TPC
  fTrackCuts[kTrackCutQGlo]->SetMinNClustersTPC(70);
  fTrackCuts[kTrackCutQGlo]->SetMaxChi2PerClusterTPC(4);
  fTrackCuts[kTrackCutQGlo]->SetAcceptKinkDaughters(kFALSE);
  fTrackCuts[kTrackCutQGlo]->SetRequireTPCRefit(kTRUE);
  // ITS
  fTrackCuts[kTrackCutQGlo]->SetRequireITSRefit(kTRUE);
  fTrackCuts[kTrackCutQGlo]->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
  fTrackCuts[kTrackCutQGlo]->SetEtaRange(-fCutEta, fCutEta);
  fTrackCuts[kTrackCutQGlo]->SetPtRange(fCutTrackPtMin, fCutTrackPtMax);
  // quality cut on ITS_SA tracks (complementary to ITS+TPC)
  fTrackCuts[kTrackCutQITS] = new AliESDtrackCuts();
  fTrackCuts[kTrackCutQITS]->SetRequireITSRefit(kTRUE);
  fTrackCuts[kTrackCutQITS]->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
  fTrackCuts[kTrackCutQITS]->SetEtaRange(-fCutEta, fCutEta);
  // primary selection for tracks with SPD hits
  fTrackCuts[kTrackCutDCAwSPD] = new AliESDtrackCuts();
  fTrackCuts[kTrackCutDCAwSPD]->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  fTrackCuts[kTrackCutDCAwSPD]->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
  fTrackCuts[kTrackCutDCAwSPD]->SetMaxDCAToVertexZ(0.5);
  fTrackCuts[kTrackCutDCAwSPD]->SetEtaRange(-fCutEta, fCutEta);
  // primary selection for tracks w/o SPD hits
  fTrackCuts[kTrackCutDCAwoSPD] = new AliESDtrackCuts();
  fTrackCuts[kTrackCutDCAwoSPD]->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kNone);
  fTrackCuts[kTrackCutDCAwoSPD]->SetMaxDCAToVertexXYPtDep("1.5*(0.0182+0.0350/pt^1.01)");
  fTrackCuts[kTrackCutDCAwoSPD]->SetMaxDCAToVertexZ(0.5);
  fTrackCuts[kTrackCutDCAwoSPD]->SetEtaRange(-fCutEta, fCutEta);

  // tracks without SPD hits
  fTrackCuts[kTrackCutNoSPD] = new AliESDtrackCuts();
  fTrackCuts[kTrackCutNoSPD]->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kNone);
  // tracks from TPC SA
  fTrackCuts[kTrackCutTPConly] = new AliESDtrackCuts();
  fTrackCuts[kTrackCutTPConly]->SetRequireTPCStandAlone(kTRUE);
}

//________________________________________________________________________
void AliAnalysisTaskFirstPhysics::UserCreateOutputObjects()
{
  PrepareOutputList();
  // Define cuts
  PrepareDefaultTrackCuts();

  // create more histograms here using UserHisto1d and UserHisto2d

  PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
bool AliAnalysisTaskFirstPhysics::GetESDEvent()
{
  fESD = 0;
  // Get a pointer to the reconstructed event
  AliVEvent *event = InputEvent();
  if (!event) {
    AliError("ERROR: Could not retrieve event.");
    return false;
  }
  // try to access the ESD information
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event);
  if (!esd) {
    AliError("InputEvent() is not AliESDEvent.");
    return false;
  }
  fESD = esd;
  return true;
}

//________________________________________________________________________
bool AliAnalysisTaskFirstPhysics::CheckVertex()
{
  // check for good reconstructed vertex
  if (!(fESD->GetPrimaryVertex()->GetStatus()) ||
      !(fESD->GetPrimaryVertexSPD()->GetStatus())) {
    return false;
  }
  // if vertex is from spd vertexZ, require more stringent cut
  if (fESD->GetPrimaryVertex()->IsFromVertexerZ() &&
      (fESD->GetPrimaryVertex()->GetDispersion() > 0.02 ||
       fESD->GetPrimaryVertex()->GetZRes() > 0.25 )) {
    return false;
  }
  if (TMath::Abs(fESD->GetPrimaryVertex()->GetZ()) > fCutVertexZ) {
    return false;
  }
  return true;
}

//________________________________________________________________________
bool AliAnalysisTaskFirstPhysics::CheckVertexMC() {
  if (!fbReadMC) {
    return false;
  }
  if (TMath::Abs(fMCEvent->GetPrimaryVertex()->GetZ()) > fCutVertexZ) {
    return false;
  }
  return true;
}

//________________________________________________________________________
bool AliAnalysisTaskFirstPhysics::PrepareMCInfo()
{
  fbReadMC = true;
  fMCEvent = MCEvent();
  if (!fMCEvent) {
    fbReadMC = false;
  }
  if (fbReadMC) {
    AliDebug(4, TString::Format("MC particles: %d", fMCEvent->GetNumberOfTracks()));
    AliDebug(4, TString::Format("number of particles in the stack: %d", fMCEvent->Stack()->GetNtrack()));
    AliGenPythiaEventHeader* pythiaGenHeader = 0;
    AliGenDPMjetEventHeader* dpmHeader = 0;
    pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(fMCEvent->GenEventHeader());
    dpmHeader = dynamic_cast<AliGenDPMjetEventHeader*>(fMCEvent->GenEventHeader());
    // MC event classification
    Int_t iProcessType = 0;
    Int_t iSD1flag = 0, iSD2flag = 0, iDDflag = 0,
      iCDflag = 0, iELflag = 0, iNDflag = 0;
    if (pythiaGenHeader) {
      iProcessType = pythiaGenHeader->ProcessType();
      iSD1flag = 92;
      iSD2flag = 93;
      iDDflag = 94;
      iELflag = 91;
      iCDflag = -2;
      iNDflag = -1;
      AliDebug(5, TString::Format("Pythia event, %d type", iProcessType));
    } else if (dpmHeader) {
      iProcessType = dpmHeader->ProcessType();
      iSD1flag = 5;
      iSD2flag = 6;
      iDDflag = 7;
      iELflag = 2;
      iCDflag = 4;
      iNDflag = 1;
      AliDebug(5, TString::Format("Phojet event, %d type", iProcessType));
    } else {
      AliError("neither pythia nor phojet event");
    }
    if (iProcessType == iSD1flag) {
      fMCProcessType = kProcSD1;
    } else if (iProcessType == iSD2flag) {
      fMCProcessType = kProcSD2;
    } else if (iProcessType == iDDflag) {
      fMCProcessType = kProcDD;
    } else if (iProcessType == iELflag) {
      fMCProcessType = kProcEL;
    } else if (iProcessType == iCDflag) {
      fMCProcessType = kProcCD;
    } else if (iProcessType == iNDflag) {
      fMCProcessType = kProcND;
    } else {
      fMCProcessType = kProcIndef;
    }
  }
  return fbReadMC;
}

//________________________________________________________________________
void AliAnalysisTaskFirstPhysics::UserExec(Option_t *)
{
  if (!GetESDEvent()) {
    AliError("Unable to read the ESD");
    return;
  }
  if (!CheckVertex()) {
    return;
  }

  // custom code comes here

  PostData(1, fOutput);
}

//________________________________________________________________________
TH1D* AliAnalysisTaskFirstPhysics::UserHisto1d(const char *name, const char *title, const char *xlabel, Int_t nbinsx, Double_t xlow, Double_t xup)
{
  TH1D *hist =  new TH1D(name, title, nbinsx, xlow, xup);
  hist->GetXaxis()->SetTitle(xlabel);
  fOutput->Add(hist);
  return hist;
}

//________________________________________________________________________
TH2D* AliAnalysisTaskFirstPhysics::UserHisto2d(const char *name, const char *title, const char *xlabel, Int_t nbinsx, Double_t xlow, Double_t xup, const char *ylabel, Int_t nbinsy, Double_t ylow, Double_t yup)
{
  TH2D *hist =  new TH2D(name, title, nbinsx, xlow, xup, nbinsy, ylow, yup);
  hist->GetXaxis()->SetTitle(xlabel);
  hist->GetYaxis()->SetTitle(ylabel);
  fOutput->Add(hist);
  return hist;
}

//________________________________________________________________________
bool AliAnalysisTaskFirstPhysics::GetHisto1FromOutput(const char *name, TH1D *&h) const
{
  h = dynamic_cast<TH1D*> (fOutput->FindObject(name));
  if (!h) {
    AliError(TString::Format("Unable to load histogram from output: %s", name));
    return false;
  } else {
    return true;
  }
}

//________________________________________________________________________
bool AliAnalysisTaskFirstPhysics::GetHisto2FromOutput(const char *name, TH2D *&h) const
{
  h = dynamic_cast<TH2D*> (fOutput->FindObject(name));
  if (!h) {
    AliError(TString::Format("Unable to load histogram from output: %s", name));
    return false;
  } else {
    return true;
  }
}

//________________________________________________________________________
void AliAnalysisTaskFirstPhysics::Terminate(Option_t *)
{
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if(!fOutput) {
    AliError("Could not retrieve TList fOutput.");
    return;
  }
}
