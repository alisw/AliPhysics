// $Id$
//
// Jet embedding task.
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

#include "AliJetEmbeddingTask.h"

ClassImp(AliJetEmbeddingTask)

//________________________________________________________________________
AliJetEmbeddingTask::AliJetEmbeddingTask() : 
  AliAnalysisTaskSE("AliJetEmbeddingTask"),
  fGeomName(),
  fTracksName(),
  fOutTracksName(),
  fCaloName(),
  fOutCaloName(),
  fEtaMin(-1),
  fEtaMax(1),
  fPhiMin(0),
  fPhiMax(TMath::Pi() * 2),
  fPtMin(0),
  fPtMax(0),
  fCopyArray(kTRUE),
  fNEmbClusters(0),
  fNEmbTracks(0),
  fGeom(0),
  fClusters(0),
  fOutClusters(0),
  fTracks(0),
  fOutTracks(0)
{
  // Default constructor.
}

//________________________________________________________________________
AliJetEmbeddingTask::AliJetEmbeddingTask(const char *name) : 
  AliAnalysisTaskSE(name),
  fGeomName(""),
  fTracksName("PicoTracks"),
  fOutTracksName("PicoTracksEmbedded"),
  fCaloName("CaloClustersCorr"),
  fOutCaloName("CaloClustersCorrEmbedded"),
  fEtaMin(-1),
  fEtaMax(1),
  fPhiMin(0),
  fPhiMax(TMath::Pi() * 2),
  fPtMin(50),
  fPtMax(60),
  fCopyArray(kTRUE),
  fNEmbClusters(0),
  fNEmbTracks(1),
  fGeom(0),
  fClusters(0),
  fOutClusters(0),
  fTracks(0),
  fOutTracks(0)
{
  // Standard constructor.

}

//________________________________________________________________________
AliJetEmbeddingTask::~AliJetEmbeddingTask()
{
  // Destructor
}

//________________________________________________________________________
void AliJetEmbeddingTask::Init()
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

  if (fNEmbTracks > 0) {
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
	fOutTracksName += "Embedded";
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

  if (fNEmbClusters > 0) {
    fClusters = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fCaloName));
 
    if (!fClusters) {
      AliError(Form("Couldn't retrieve clusters with name %s!", fCaloName.Data()));
      return;
    }

    if (!fOutClusters) {
      fOutCaloName = fCaloName;
      if (fCopyArray) {
	fOutCaloName += "Embedded";
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
    if (fGeomName.Length()>0) {
      fGeom = AliEMCALGeometry::GetInstance(fGeomName);
      if (!fGeom)
        AliError(Form("Could not get geometry with name %s!", fGeomName.Data()));
    } else {
      fGeom = AliEMCALGeometry::GetInstance();
      if (!fGeom) 
        AliError("Could not get default geometry!");
    }
  }
}

//________________________________________________________________________
void AliJetEmbeddingTask::Embed() 
{
  // Embed particles.
  
  if (fNEmbClusters > 0 && fOutClusters) {
    Float_t emEtaMax = fEtaMax;
    Float_t emEtaMin = fEtaMin;
    Float_t emPhiMax = fPhiMax;
    Float_t emPhiMin = fPhiMin;

    // hard-coded Emcal boundaries
    const Float_t EmcalEtaMin = -0.7;
    const Float_t EmcalEtaMax = 0.7;
    const Float_t EmcalPhiMin = 80 * TMath::DegToRad();
    const Float_t EmcalPhiMax = 180 * TMath::DegToRad();

    if (emEtaMax > EmcalEtaMax) emEtaMax = EmcalEtaMax;
    if (emEtaMax < EmcalEtaMin) emEtaMax = EmcalEtaMin;
    if (emEtaMin > EmcalEtaMax) emEtaMin = EmcalEtaMax;
    if (emEtaMin < EmcalEtaMin) emEtaMin = EmcalEtaMin;

    if (emPhiMax > EmcalPhiMax) emPhiMax = EmcalPhiMax;
    if (emPhiMax < EmcalPhiMin) emPhiMax = EmcalPhiMin;
    if (emPhiMin > EmcalPhiMax) emPhiMin = EmcalPhiMax;
    if (emPhiMin < EmcalPhiMin) emPhiMin = EmcalPhiMin;

    const Int_t nClusters = fOutClusters->GetEntriesFast();
    TClonesArray digits("AliEMCALDigit", 1);
    for (Int_t i = 0; i < fNEmbClusters; ++i) {

      Double_t eta = 0;
      Double_t phi = 0;
      Int_t absId = 0;

      Int_t repeats = 0;
      do {
	eta = gRandom->Rndm() * (emEtaMax - emEtaMin) + emEtaMin;
	phi = gRandom->Rndm() * (emPhiMax - emPhiMin) + emPhiMin;
	fGeom->GetAbsCellIdFromEtaPhi(eta, phi, absId);  
	repeats++;
      } while (absId == -1 && repeats < 100);

      if (!(absId > -1)) {
	AliWarning(Form("Could not embed cluster! Random eta-phi extracted more than 100 times!\n"
			"eta [%f, %f], phi [%f, %f]\n", emEtaMin, emEtaMax, emPhiMin, emPhiMax));
	continue;
      }

      Double_t pt  = gRandom->Rndm() * (fPtMax - fPtMin) + fPtMin;

      TLorentzVector nPart;
      nPart.SetPtEtaPhiM(pt, eta, phi, 0);
      Double_t e = nPart.E();
      
      if (absId == -1) {
	AliWarning(Form("Unable to embed cluster in eta = %f, phi = %f!"
                        " Maybe the eta-phi range is not inside the EMCal acceptance (eta = [%f, %f], phi = [%f, %f])", 
			eta, phi, fEtaMin, fEtaMax, fPhiMin, fPhiMax));
	continue;
      }
      
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
      
      AliVCluster *cluster = static_cast<AliVCluster*>(fOutClusters->New(nClusters + i));
      cluster->SetType(AliVCluster::kEMCALClusterv1);
      cluster->SetE(recPoint.GetEnergy());
      cluster->SetPosition(g);
      cluster->SetNCells(1);
      UShort_t shortAbsId = absId;
      cluster->SetCellsAbsId(&shortAbsId);
      Double32_t fract = 1;
      cluster->SetCellsAmplitudeFraction(&fract);
      cluster->SetID(nClusters + i);
      cluster->SetEmcCpvDistance(-1);
      cluster->SetChi2(100); // MC flag!
    }
  }
 
  if (fNEmbTracks > 0 && fOutTracks) {
    Int_t nTracks = fOutTracks->GetEntriesFast();
    for (Int_t i = 0; i < fNEmbTracks; ++i) {
      Double_t pt  = gRandom->Rndm() * (fPtMax - fPtMin) + fPtMin;
      Double_t eta = gRandom->Rndm() * (fEtaMax - fEtaMin) + fEtaMin;
      Double_t phi = gRandom->Rndm() * (fPhiMax - fPhiMin) + fPhiMin;
      
      new ((*fOutTracks)[nTracks + i]) AliPicoTrack(pt, 
						    eta, 
						    phi, 
						    1, 
						    100,    // MC flag!      
						    0, 
						    0, 
						    kFALSE);
    }
  }
}

//________________________________________________________________________
void AliJetEmbeddingTask::UserExec(Option_t *) 
{
  // Execute per event.

  Init();

  Embed();
}
