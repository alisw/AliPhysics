// $Id$
//
// Jet modelling task.
//
// Author: S.Aiola, C.Loizides

#include "AliJetModelBaseTask.h"

#include <TClonesArray.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TRandom3.h>

#include "AliAODCaloCluster.h"
#include "AliAnalysisManager.h"
#include "AliEMCALDigit.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecPoint.h"
#include "AliESDCaloCluster.h"
#include "AliLog.h"
#include "AliPicoTrack.h"
#include "AliVCluster.h"
#include "AliVEvent.h"

ClassImp(AliJetModelBaseTask)

//________________________________________________________________________
AliJetModelBaseTask::AliJetModelBaseTask() : 
  AliAnalysisTaskSE("AliJetModelBaseTask"),
  fGeomName(),
  fTracksName(),
  fOutTracksName(),
  fCaloName(),
  fOutCaloName(),
  fSuffix(),
  fEtaMin(-1),
  fEtaMax(1),
  fPhiMin(0),
  fPhiMax(TMath::Pi() * 2),
  fPtMin(0),
  fPtMax(0),
  fCopyArray(kTRUE),
  fNClusters(0),
  fNTracks(0),
  fMarkMC(kTRUE),
  fPtSpectrum(0),
  fIsInit(0),
  fGeom(0),
  fClusters(0),
  fOutClusters(0),
  fTracks(0),
  fOutTracks(0)
{
  // Default constructor.
}

//________________________________________________________________________
AliJetModelBaseTask::AliJetModelBaseTask(const char *name) : 
  AliAnalysisTaskSE(name),
  fGeomName(""),
  fTracksName("PicoTracks"),
  fOutTracksName("PicoTracksEmbedded"),
  fCaloName("CaloClustersCorr"),
  fOutCaloName("CaloClustersCorrEmbedded"),
  fSuffix("Processed"),
  fEtaMin(-1),
  fEtaMax(1),
  fPhiMin(0),
  fPhiMax(TMath::Pi() * 2),
  fPtMin(50),
  fPtMax(60),
  fCopyArray(kTRUE),
  fNClusters(0),
  fNTracks(1),
  fMarkMC(kTRUE),
  fPtSpectrum(0),
  fIsInit(0),
  fGeom(0),
  fClusters(0),
  fOutClusters(0),
  fTracks(0),
  fOutTracks(0)
{
  // Standard constructor.
}

//________________________________________________________________________
AliJetModelBaseTask::~AliJetModelBaseTask()
{
  // Destructor
}

//________________________________________________________________________
void AliJetModelBaseTask::UserExec(Option_t *) 
{
  // Execute per event.

  if (!fIsInit) {
    ExecOnce();
    fIsInit = 1;
  }
  Run();
}

//________________________________________________________________________
AliVCluster* AliJetModelBaseTask::AddCluster(Double_t e, Double_t eta, Double_t phi)
{
  // Add a cluster to the event.

  Int_t absId = 0;
  if (eta < -100 || phi < 0) {
    GetRandomCell(eta, phi, absId);
  }
  else {
    fGeom->EtaPhiFromIndex(absId, eta, phi);
  }

  if (absId == -1) {
    AliWarning(Form("Unable to embed cluster in eta = %f, phi = %f!"
		    " Maybe the eta-phi range is not inside the EMCal acceptance (eta = [%f, %f], phi = [%f, %f])", 
		    eta, phi, fEtaMin, fEtaMax, fPhiMin, fPhiMax));
    return 0;
  } 

  if (e < 0) {
    Double_t pt = GetRandomPt();
    TLorentzVector nPart;
    nPart.SetPtEtaPhiM(pt, eta, phi, 0);
    e = nPart.E();
  }

  return AddCluster(e, absId);
}
      
//________________________________________________________________________
AliVCluster* AliJetModelBaseTask::AddCluster(Double_t e, Int_t absId)
{
  // Add a cluster to the event.

  const Int_t nClusters = fOutClusters->GetEntriesFast();

  TClonesArray digits("AliEMCALDigit", 1);

  AliEMCALDigit *digit = static_cast<AliEMCALDigit*>(digits.New(0));
  digit->SetId(absId);
  digit->SetIndexInList(0);
  digit->SetType(AliEMCALDigit::kHG);
  digit->SetAmplitude(e);
      
  AliEMCALRecPoint recPoint("");
  recPoint.AddDigit(*digit, e, kFALSE);
  recPoint.EvalGlobalPosition(0, &digits);

  TVector3 gpos;
  recPoint.GetGlobalPosition(gpos);
  Float_t g[3];
  gpos.GetXYZ(g);
      
  AliVCluster *cluster = static_cast<AliVCluster*>(fOutClusters->New(nClusters));
  cluster->SetType(AliVCluster::kEMCALClusterv1);
  cluster->SetE(recPoint.GetEnergy());
  cluster->SetPosition(g);
  cluster->SetNCells(1);
  UShort_t shortAbsId = absId;
  cluster->SetCellsAbsId(&shortAbsId);
  Double32_t fract = 1;
  cluster->SetCellsAmplitudeFraction(&fract);
  cluster->SetID(nClusters);
  cluster->SetEmcCpvDistance(-1);
  if (fMarkMC)
    cluster->SetChi2(100); // MC flag!

  return cluster;
}

//________________________________________________________________________
AliPicoTrack* AliJetModelBaseTask::AddTrack(Double_t pt, Double_t eta, Double_t phi)
{
  // Add a track to the event.

  const Int_t nTracks = fOutTracks->GetEntriesFast();
  
  if (pt < 0) 
    pt = GetRandomPt();
  if (eta < -100) 
    eta = GetRandomEta();
  if (phi < 0) 
    phi = GetRandomPhi();

  Int_t label = fMarkMC ? 100 : 0;

  AliPicoTrack *track = new ((*fOutTracks)[nTracks]) AliPicoTrack(pt, 
								  eta, 
								  phi, 
								  1, 
								  label,    // MC flag!      
								  0, 
								  0, 
								  kFALSE);
  return track;
}

//________________________________________________________________________
void AliJetModelBaseTask::CopyClusters()
{
  // Copy all the clusters in the new collection

  Bool_t esdMode = (Bool_t)(fClusters->GetClass()->GetBaseClass("AliESDCaloCluster") != 0);
  const Int_t nClusters = fClusters->GetEntriesFast();
  
  if (esdMode) {
    for (Int_t i = 0; i < nClusters; ++i) {
      AliESDCaloCluster *esdcluster = static_cast<AliESDCaloCluster*>(fClusters->At(i));
      if (!esdcluster)
	continue;
      if (!esdcluster->IsEMCAL())
	continue;
      new ((*fOutClusters)[i]) AliESDCaloCluster(*esdcluster);
    }
  }
  else {
    for (Int_t i = 0; i < nClusters; ++i) {
      AliAODCaloCluster *aodcluster = static_cast<AliAODCaloCluster*>(fClusters->At(i));
      if (!aodcluster)
	continue;
      if (!aodcluster->IsEMCAL())
	continue;
      new ((*fOutClusters)[i]) AliAODCaloCluster(*aodcluster);
    }
  }
}

//________________________________________________________________________
void AliJetModelBaseTask::CopyTracks()
{
  const Int_t nTracks = fTracks->GetEntriesFast();
  for (Int_t i = 0; i < nTracks; ++i) {
    AliPicoTrack *track = static_cast<AliPicoTrack*>(fTracks->At(i));
    if (!track)
      continue;
    new ((*fOutTracks)[i]) AliPicoTrack(*track);
  }
}

//________________________________________________________________________
void AliJetModelBaseTask::ExecOnce()
{
  // Init task.

  if (fPtMax < fPtMin) {
    AliWarning (Form("PtMax (%f) < PtMin (%f), setting PtMax = PtMin = %f", fPtMax, fPtMin, fPtMin));
    fPtMax = fPtMin;
  }

  if (fEtaMax < fEtaMin) {
    AliWarning (Form("EtaMax (%f) < EtaMin (%f), setting EtaMax = EtaMin = %f", fEtaMax, fEtaMin, fEtaMin));
    fEtaMax = fEtaMin;
  }

  if (fPhiMax < fPhiMin) {
    AliWarning (Form("PhiMax (%f) < PhiMin (%f), setting PhiMax = PhiMin = %f", fPhiMax, fPhiMin, fPhiMin));
    fPhiMax = fPhiMin;
  }

  if (fNTracks > 0 && !fTracks && !fTracksName.IsNull()) {
    fTracks = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTracksName));
    if (!fTracks) {
      AliError(Form("%s: Couldn't retrieve tracks with name %s!", GetName(), fTracksName.Data()));
      return;
    }
    
    if (!fTracks->GetClass()->GetBaseClass("AliPicoTrack")) {
      AliError(Form("%s: Collection %s does not contain AliPicoTrack objects!", GetName(), fTracksName.Data())); 
      fTracks = 0;
      return;
    }

    if (!fOutTracks) {
      fOutTracksName = fTracksName;
      if (fCopyArray) {
	fOutTracksName += fSuffix;
	fOutTracks = new TClonesArray("AliPicoTrack", fTracks->GetSize());
	fOutTracks->SetName(fOutTracksName);
      }
      else {
	fOutTracks = fTracks;
      }
    }

    if (fCopyArray) {
      if (!(InputEvent()->FindListObject(fOutTracksName)))
	InputEvent()->AddObject(fOutTracks);
    }
  }

  if (fNClusters > 0 && !fClusters && !fCaloName.IsNull()) {
    fClusters = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fCaloName));
 
    if (!fClusters) {
      AliError(Form("%s: Couldn't retrieve clusters with name %s!", GetName(), fCaloName.Data()));
      return;
    }

    if (!fClusters->GetClass()->GetBaseClass("AliVCluster")) {
      AliError(Form("%s: Collection %s does not contain AliVCluster objects!", GetName(), fCaloName.Data())); 
      fClusters = 0;
      return;
    }

    if (!fOutClusters) {
      fOutCaloName = fCaloName;
      if (fCopyArray) {
	fOutCaloName += fSuffix;
	fOutClusters = new TClonesArray(fClusters->GetClass()->GetName(), fClusters->GetSize());
	fOutClusters->SetName(fOutCaloName);
      }
      else {
	fOutClusters = fClusters;
      }
    }

    if (fCopyArray) {
      if (!(InputEvent()->FindListObject(fOutCaloName)))
	InputEvent()->AddObject(fOutClusters);
    }

    if (!fGeom) {
      if (fGeomName.Length() > 0) {
	fGeom = AliEMCALGeometry::GetInstance(fGeomName);
	if (!fGeom)
	  AliError(Form("Could not get geometry with name %s!", fGeomName.Data()));
      } else {
	fGeom = AliEMCALGeometry::GetInstance();
	if (!fGeom) 
	  AliError("Could not get default geometry!");
      }
    }
    
    const Double_t EmcalMinEta = fGeom->GetArm1EtaMin();
    const Double_t EmcalMaxEta = fGeom->GetArm1EtaMax();
    const Double_t EmcalMinPhi = fGeom->GetArm1PhiMin() * TMath::DegToRad();
    const Double_t EmcalMaxPhi = fGeom->GetArm1PhiMax() * TMath::DegToRad();
    
    if (fEtaMax > EmcalMaxEta) fEtaMax = EmcalMaxEta;
    if (fEtaMax < EmcalMinEta) fEtaMax = EmcalMinEta;
    if (fEtaMin > EmcalMaxEta) fEtaMin = EmcalMaxEta;
    if (fEtaMin < EmcalMinEta) fEtaMin = EmcalMinEta;
  
    if (fPhiMax > EmcalMaxPhi) fPhiMax = EmcalMaxPhi;
    if (fPhiMax < EmcalMinPhi) fPhiMax = EmcalMinPhi;
    if (fPhiMin > EmcalMaxPhi) fPhiMin = EmcalMaxPhi;
    if (fPhiMin < EmcalMinPhi) fPhiMin = EmcalMinPhi;
  }
}

//________________________________________________________________________
void AliJetModelBaseTask::GetRandomCell(Double_t &eta, Double_t &phi, Int_t &absId)
{
  // Get random cell.

  Int_t repeats = 0;
  Double_t rndEta = eta;
  Double_t rndPhi = phi;
  do {
    if (eta < -100)
      rndEta = GetRandomEta();
    if (phi < 0)
      rndPhi = GetRandomPhi();
    fGeom->GetAbsCellIdFromEtaPhi(rndEta, rndPhi, absId);  
    repeats++;
  } while (absId == -1 && repeats < 100);
  
  if (!(absId > -1)) {
    AliWarning(Form("Could not extract random cluster! Random eta-phi extracted more than 100 times!\n"
		    "eta [%f, %f], phi [%f, %f]\n", fEtaMin, fEtaMax, fPhiMin, fPhiMax));
  }
  else {
    eta = rndEta;
    phi = rndPhi;
  }
}

//________________________________________________________________________
Double_t AliJetModelBaseTask::GetRandomEta()
{
  // Get random eta.

  return gRandom->Rndm() * (fEtaMax - fEtaMin) + fEtaMin;
}

//________________________________________________________________________
Double_t AliJetModelBaseTask::GetRandomPhi()
{
  // Get random phi.

  return gRandom->Rndm() * (fPhiMax - fPhiMin) + fPhiMin;
}

//________________________________________________________________________
Double_t AliJetModelBaseTask::GetRandomPt()
{
  // Get random pt.

  if (fPtSpectrum)
    return fPtSpectrum->GetRandom();
  else
    return gRandom->Rndm() * (fPtMax - fPtMin) + fPtMin;
}

//________________________________________________________________________
void AliJetModelBaseTask::Run() 
{
  // Run.
}
