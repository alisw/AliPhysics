// $Id: AliAnalysisTaskEmcal.cxx 56756 2012-05-30 05:03:02Z loizides $
//
// Emcal base analysis task.
//
// Author: S.Aiola

#include "AliAnalysisTaskEmcal.h"

#include <TObject.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TList.h>

#include "AliESDEvent.h"
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliVCluster.h"
#include "AliVParticle.h"
#include "AliMCParticle.h"
#include "AliEmcalJet.h"
#include "AliVEventHandler.h"
#include "AliLog.h"
#include "AliEMCALGeometry.h"

ClassImp(AliAnalysisTaskEmcal)

//________________________________________________________________________
AliAnalysisTaskEmcal::AliAnalysisTaskEmcal() : 
  AliAnalysisTaskSE("AliAnalysisTaskEmcal"),
  fAnaType(kTPC),
  fInitialized(kFALSE),
  fCreateHisto(kTRUE),
  fTracksName("Tracks"),
  fCaloName("CaloClusters"),
  fNbins(500),
  fMinBinPt(0),
  fMaxBinPt(250),
  fPtCut(0.15),
  fTracks(0),
  fCaloClusters(0),
  fCent(0),
  fCentBin(-1),
  fBeamType(kNA),
  fOutput(0)
{
  // Default constructor.

  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;

  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class()); 
}

//________________________________________________________________________
AliAnalysisTaskEmcal::AliAnalysisTaskEmcal(const char *name) : 
  AliAnalysisTaskSE(name),
  fAnaType(kTPC),
  fInitialized(kFALSE),
  fCreateHisto(kTRUE),
  fTracksName("Tracks"),
  fCaloName("CaloClusters"),
  fNbins(500),
  fMinBinPt(0),
  fMaxBinPt(250),
  fPtCut(0.15),
  fTracks(0),
  fCaloClusters(0),
  fCent(0),
  fCentBin(-1),
  fBeamType(kNA),
  fOutput(0)
{
  // Standard constructor.

  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;

  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class()); 
}

//________________________________________________________________________
AliAnalysisTaskEmcal::AliAnalysisTaskEmcal(const char *name, Bool_t histo) : 
  AliAnalysisTaskSE(name),
  fAnaType(kTPC),
  fInitialized(kFALSE),
  fCreateHisto(histo),
  fTracksName("Tracks"),
  fCaloName("CaloClusters"),
  fNbins(500),
  fMinBinPt(0),
  fMaxBinPt(250),
  fPtCut(0.15),
  fTracks(0),
  fCaloClusters(0),
  fCent(0),
  fCentBin(-1),
  fBeamType(kNA),
  fOutput(0)
{
  // Standard constructor.

  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;

  if (fCreateHisto) {
    // Output slot #1 writes into a TH1 container
    DefineOutput(1, TList::Class()); 
  }
}

//________________________________________________________________________
AliAnalysisTaskEmcal::~AliAnalysisTaskEmcal()
{
  // Destructor
}

//________________________________________________________________________
void AliAnalysisTaskEmcal::UserCreateOutputObjects()
{
  // User create outputs.
}

//_____________________________________________________
AliAnalysisTaskEmcal::BeamType AliAnalysisTaskEmcal::GetBeamType()
{
  // Get beam type : pp-AA-pA
  // ESDs have it directly, AODs get it from hardcoded run number ranges

  AliESDEvent *esd = dynamic_cast<AliESDEvent*>(InputEvent());
  if (esd) {
    const AliESDRun *run = esd->GetESDRun();
    TString beamType = run->GetBeamType();
    if (beamType == "p-p")
      return kpp;
    else if (beamType == "A-A")
      return kAA;
    else if (beamType == "p-A")
      return kpA;
    else
      return kNA;
  }
  else
  {
    Int_t runNumber = InputEvent()->GetRunNumber();
    if ((runNumber >= 136851 && runNumber <= 139517) ||  // LHC10h
	(runNumber >= 166529 && runNumber <= 170593))    // LHC11h
    {
      return kAA;
    }
    else 
    {
      return kpp;
    }
  }  
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcal::RetrieveEventObjects()
{
  // Retrieve objects from event.

  if (!InputEvent()) {
    AliError("Could not retrieve event! Returning...");
    return kFALSE;
  }

  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;
  InputEvent()->GetPrimaryVertex()->GetXYZ(fVertex);

  fBeamType = GetBeamType();

  if (fBeamType == kAA) {
    AliCentrality *aliCent = InputEvent()->GetCentrality();
    if (aliCent) {
      fCent = aliCent->GetCentralityPercentile("V0M");
      if      (fCent >=  0 && fCent <   10) fCentBin = 0;
      else if (fCent >= 10 && fCent <   30) fCentBin = 1;
      else if (fCent >= 30 && fCent <   50) fCentBin = 2;
      else if (fCent >= 50 && fCent <= 100) fCentBin = 3; 
      else {
	AliWarning(Form("Negative centrality: %f. Assuming 99", fCent));
	fCentBin = 3;
      }
    }
    else {
      AliWarning(Form("Could not retrieve centrality information! Assuming 99"));
      fCentBin = 3;
    }
  }
  else {
    fCent = 99;
    fCentBin = 0;
  }

  if ((!fCaloName.IsNull()) && (fAnaType == kEMCAL)) {
    fCaloClusters =  dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fCaloName));
    if (!fCaloClusters) {
      AliWarning(Form("Could not retrieve clusters %s!", fCaloName.Data())); 
    }
  }

  if (!fTracksName.IsNull()) {
    fTracks = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTracksName));
    if (!fTracks) {
      AliWarning(Form("Could not retrieve tracks %s!", fTracksName.Data())); 
    }
  }

  return kTRUE;
}


//________________________________________________________________________
Bool_t AliAnalysisTaskEmcal::AcceptCluster(AliVCluster* clus, Bool_t acceptMC) const
{
  // Return true if cluster is accepted.
  if (!clus->IsEMCAL())
    return kFALSE;

  if (!acceptMC && clus->Chi2() == 100)
    return kFALSE;

  TLorentzVector nPart;
  clus->GetMomentum(nPart, const_cast<Double_t*>(fVertex));

  if (nPart.Et() < fPtCut)
    return kFALSE;

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcal::AcceptTrack(AliVTrack* track, Bool_t acceptMC) const
{
  // Return true if track is accepted.
  if (!acceptMC && track->GetLabel() == 100)
    return kFALSE;

  if (track->Pt() < fPtCut)
    return kFALSE;
  
  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskEmcal::UserExec(Option_t *) 
{
  // Main loop, called for each event.

  if (!RetrieveEventObjects())
    return;

  if (!Run())
    return;

  if (!FillHistograms())
    return;
    
  if (fCreateHisto) {
    // information for this iteration of the UserExec in the container
    PostData(1, fOutput);
  }
}

//________________________________________________________________________
void AliAnalysisTaskEmcal::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}
