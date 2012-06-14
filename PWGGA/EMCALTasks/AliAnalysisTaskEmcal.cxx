// $Id: AliAnalysisTaskEmcal.cxx 56756 2012-05-30 05:03:02Z loizides $
//
// Emcal base analysis task.
//
// Author: S.Aiola

#include "AliAnalysisTaskEmcal.h"

#include <TChain.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TObject.h>

#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliEMCALGeometry.h"
#include "AliESDEvent.h"
#include "AliEmcalJet.h"
#include "AliEmcalParticle.h"
#include "AliLog.h"
#include "AliMCParticle.h"
#include "AliVCluster.h"
#include "AliVEventHandler.h"
#include "AliVParticle.h"

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
  fClusPtCut(0.15),
  fTrackPtCut(0.15),
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
  fClusPtCut(0.15),
  fTrackPtCut(0.15),
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
    DefineOutput(1, TList::Class()); 
  }
}

//________________________________________________________________________
AliAnalysisTaskEmcal::~AliAnalysisTaskEmcal()
{
  // Destructor
}

//________________________________________________________________________
void AliAnalysisTaskEmcal::UserExec(Option_t *) 
{
  // Main loop, called for each event.

  if (!fInitialized)
    Init();

  if (!fInitialized)
    return;

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
Bool_t AliAnalysisTaskEmcal::AcceptCluster(AliVCluster *clus, Bool_t acceptMC) const
{
  // Return true if cluster is accepted.

  if (!clus)
    return kFALSE;

  if (!clus->IsEMCAL())
    return kFALSE;

  if (!acceptMC && clus->Chi2() == 100)
    return kFALSE;

  TLorentzVector nPart;
  clus->GetMomentum(nPart, const_cast<Double_t*>(fVertex));

  if (nPart.Et() < fClusPtCut)
    return kFALSE;

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcal::AcceptEmcalPart(AliEmcalParticle *part, Bool_t acceptMC) const
{
  // Return true if EMCal particle is accepted.

  if (!part)
    return kFALSE;

  if (fAnaType == kEMCAL && !part->IsEMCAL())
    return kFALSE;

  if ((part->IsTrack() && part->Pt() < fTrackPtCut) || (part->IsCluster() && part->Pt() < fClusPtCut))
    return kFALSE;

  if (!acceptMC && part->IsMC())
    return kFALSE;

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcal::AcceptTrack(AliVTrack *track, Bool_t acceptMC) const
{
  // Return true if track is accepted.

  if (!track)
    return kFALSE;

  if (!acceptMC && track->GetLabel() == 100)
    return kFALSE;

  if (track->Pt() < fTrackPtCut)
    return kFALSE;
  
  return kTRUE;
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
  } else {
    Int_t runNumber = InputEvent()->GetRunNumber();
    if ((runNumber >= 136851 && runNumber <= 139517) ||  // LHC10h
	(runNumber >= 166529 && runNumber <= 170593))    // LHC11h
    {
      return kAA;
    } else {
      return kpp;
    }
  }  
}

//________________________________________________________________________
void AliAnalysisTaskEmcal::Init()
{
  // Init the analysis.

  SetInitialized();
}

//________________________________________________________________________
TClonesArray *AliAnalysisTaskEmcal::GetArrayFromEvent(const char *name, const char *clname)
{
  // Get array from event.

  TClonesArray *arr = 0;
  TString sname(name);
  if (!sname.IsNull()) {
    arr = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(sname));
    if (!arr) {
      AliWarning(Form("%s: Could not retrieve array with name %s!", GetName(), name)); 
      return 0;
    }
  }
  if (!clname)
    return arr;

  TString objname(arr->GetClass()->GetName());
  TClass cls(objname);
  if (!cls.InheritsFrom(clname)) {
    AliWarning(Form("%s: Objects of type %s in %s are not inherited from %s!", 
                    GetName(), cls.GetName(), name, clname)); 
    return 0;
  }
  return arr;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcal::RetrieveEventObjects()
{
  // Retrieve objects from event.

  if (!InputEvent()) {
    AliError(Form("%s: Could not retrieve event! Returning!", GetName()));
    return kFALSE;
  }

  if (!fCaloName.IsNull() && (fAnaType == kEMCAL) && !fCaloClusters) {
    fCaloClusters =  dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fCaloName));
    if (!fCaloClusters) {
      AliError(Form("%s: Could not retrieve clusters %s!", GetName(), fCaloName.Data())); 
      return kFALSE;
    }
    else {
      TClass *cl = fCaloClusters->GetClass();
      if (!cl->GetBaseClass("AliVCluster") && !cl->GetBaseClass("AliEmcalParticle")) {
	AliError(Form("%s: Collection %s does not contain AliVCluster nor AliEmcalParticle objects!", GetName(), fCaloName.Data())); 
	fCaloClusters = 0;
	return kFALSE;
      }
    }
  }

  if (!fTracksName.IsNull() && !fTracks) {
    fTracks = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTracksName));
    if (!fTracks) {
      AliError(Form("%s: Could not retrieve tracks %s!", GetName(), fTracksName.Data())); 
      return kFALSE;
    }
    else {
      TClass *cl = fTracks->GetClass();
      if (!cl->GetBaseClass("AliVParticle") && !cl->GetBaseClass("AliEmcalParticle")) {
	AliError(Form("%s: Collection %s does not contain AliVParticle nor AliEmcalParticle objects!", GetName(), fTracksName.Data())); 
	fTracks = 0;
	return kFALSE;
      }
    }
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
	AliWarning(Form("%s: Negative centrality: %f. Assuming 99", GetName(), fCent));
	fCentBin = 3;
      }
    } else {
      AliWarning(Form("%s: Could not retrieve centrality information! Assuming 99", GetName()));
      fCentBin = 3;
    }
  } else {
    fCent = 99;
    fCentBin = 0;
  }

  return kTRUE;
}

