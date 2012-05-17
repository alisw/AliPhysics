// $Id$
//
// Jet modelling task.
//
// Author: Salvatore Aiola, Constantin Loizides

#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TRandom3.h>

#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliVCluster.h"
#include "AliEMCALDigit.h"
#include "AliEMCALRecPoint.h"
#include "AliPicoTrack.h"
#include "AliEMCALGeometry.h"
#include "AliLog.h"

#include "AliJetModelBaseTask.h"

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
void AliJetModelBaseTask::Init()
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

  if (fNTracks > 0) {
    fTracks = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTracksName));

    if (strcmp(fTracks->GetClass()->GetName(), "AliPicoTrack")) {
      AliError("Can only embed PicoTracks!");
      return;
    }
 
    if (!fTracks) {
      AliError(Form("Couldn't retrieve tracks with name %s!", fTracksName.Data()));
      return;
    }

    if (!fOutTracks) {
      fOutTracksName = fTracksName;
      if (fCopyArray) {
	fOutTracksName += fSuffix;
	fOutTracks = new TClonesArray(*fTracks);
	fOutTracks->SetName(fOutTracksName);
      }
      else {
	fOutTracks = fTracks;
      }
    }

    if (fCopyArray) {
      if (!(InputEvent()->FindListObject(fOutTracksName)))
	InputEvent()->AddObject(fOutTracks);
      fOutTracks->Clear();
    }
  }

  if (fNClusters > 0) {
    fClusters = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fCaloName));
 
    if (!fClusters) {
      AliError(Form("Couldn't retrieve clusters with name %s!", fCaloName.Data()));
      return;
    }

    if (!fOutClusters) {
      fOutCaloName = fCaloName;
      if (fCopyArray) {
	fOutCaloName += fSuffix;
	fOutClusters = new TClonesArray(*fClusters);
	fOutClusters->SetName(fOutCaloName);
      }
      else {
	fOutClusters = fClusters;
      }
    }

    if (fCopyArray) {
      if (!(InputEvent()->FindListObject(fOutCaloName)))
	InputEvent()->AddObject(fOutClusters);
      fOutClusters->Clear();
    }
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

    // hard-coded Emcal boundaries
    const Float_t EmcalEtaMin = -0.7;
    const Float_t EmcalEtaMax = 0.7;
    const Float_t EmcalPhiMin = 80 * TMath::DegToRad();
    const Float_t EmcalPhiMax = 180 * TMath::DegToRad();

    if (fEtaMax > EmcalEtaMax) fEtaMax = EmcalEtaMax;
    if (fEtaMax < EmcalEtaMin) fEtaMax = EmcalEtaMin;
    if (fEtaMin > EmcalEtaMax) fEtaMin = EmcalEtaMax;
    if (fEtaMin < EmcalEtaMin) fEtaMin = EmcalEtaMin;
  
    if (fPhiMax > EmcalPhiMax) fPhiMax = EmcalPhiMax;
    if (fPhiMax < EmcalPhiMin) fPhiMax = EmcalPhiMin;
    if (fPhiMin > EmcalPhiMax) fPhiMin = EmcalPhiMax;
    if (fPhiMin < EmcalPhiMin) fPhiMin = EmcalPhiMin;
  }
}

//________________________________________________________________________
void AliJetModelBaseTask::GetRandomCell(Double_t &eta, Double_t &phi, Int_t &absId)
{
  Int_t repeats = 0;

  do {
    eta = gRandom->Rndm() * (fEtaMax - fEtaMin) + fEtaMin;
    phi = gRandom->Rndm() * (fPhiMax - fPhiMin) + fPhiMin;
    fGeom->GetAbsCellIdFromEtaPhi(eta, phi, absId);  
    repeats++;
  } while (absId == -1 && repeats < 100);
  
  if (!(absId > -1)) {
    AliWarning(Form("Could not extract random cluster! Random eta-phi extracted more than 100 times!\n"
		    "eta [%f, %f], phi [%f, %f]\n", fEtaMin, fEtaMax, fPhiMin, fPhiMax));
  }
}

//________________________________________________________________________
Double_t AliJetModelBaseTask::GetRandomEta()
{
  return gRandom->Rndm() * (fEtaMax - fEtaMin) + fEtaMin;
}

//________________________________________________________________________
Double_t AliJetModelBaseTask::GetRandomPhi()
{
  return gRandom->Rndm() * (fPhiMax - fPhiMin) + fPhiMin;
}

//________________________________________________________________________
Double_t AliJetModelBaseTask::GetRandomPt()
{
  return gRandom->Rndm() * (fPtMax - fPtMin) + fPtMin;
}

//________________________________________________________________________
AliVCluster* AliJetModelBaseTask::AddCluster(Double_t e, Double_t eta, Double_t phi)
{
  Int_t absId = 0;
  if (eta < 0 || phi < 0) {
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
  cluster->SetChi2(100); // MC flag!

  return cluster;
}

//________________________________________________________________________
AliPicoTrack* AliJetModelBaseTask::AddTrack(Double_t pt, Double_t eta, Double_t phi)
{
  const Int_t nTracks = fOutTracks->GetEntriesFast();
  
  if (pt < 0) 
    pt = GetRandomPt();
  if (eta < 0) 
    eta = GetRandomEta();
  if (phi < 0) 
    phi = GetRandomPhi();
  
  AliPicoTrack *track = new ((*fOutTracks)[nTracks]) AliPicoTrack(pt, 
						eta, 
						phi, 
						1, 
						100,    // MC flag!      
						0, 
						0, 
						kFALSE);
  return track;
}

//________________________________________________________________________
void AliJetModelBaseTask::Run() 
{

}

//________________________________________________________________________
void AliJetModelBaseTask::UserExec(Option_t *) 
{
  // Execute per event.

  Init();

  Run();

}
