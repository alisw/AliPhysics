// $Id$
//
// Jet modelling task.
//
// Author: S.Aiola, C.Loizides

#include "AliJetModelBaseTask.h"

#include <TClonesArray.h>
#include <TH1I.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TList.h>

#include "AliVEvent.h"
#include "AliAODCaloCluster.h"
#include "AliESDCaloCluster.h"
#include "AliVCluster.h"
#include "AliEMCALDigit.h"
#include "AliEMCALRecPoint.h"
#include "AliESDCaloCells.h"
#include "AliAODCaloCells.h"
#include "AliAODMCParticle.h"
#include "AliVCaloCells.h"
#include "AliPicoTrack.h"
#include "AliEMCALGeometry.h"
#include "AliLog.h"

ClassImp(AliJetModelBaseTask)

//________________________________________________________________________
AliJetModelBaseTask::AliJetModelBaseTask() : 
  AliAnalysisTaskSE("AliJetModelBaseTask"),
  fGeomName(),
  fTracksName(),
  fOutTracksName(),
  fCaloName(),
  fOutCaloName(),
  fCellsName(),
  fOutCellsName(),
  fMCParticlesName(),
  fOutMCParticlesName(),
  fSuffix(),
  fEtaMin(-1),
  fEtaMax(1),
  fPhiMin(0),
  fPhiMax(TMath::Pi() * 2),
  fPtMin(0),
  fPtMax(0),
  fCopyArray(kTRUE),
  fNClusters(0),
  fNCells(0),
  fNTracks(0),
  fMarkMC(99999),
  fPtSpectrum(0),
  fQAhistos(kFALSE),
  fIsInit(0),
  fGeom(0),
  fClusters(0),
  fOutClusters(0),
  fTracks(0),
  fOutTracks(0),
  fCaloCells(0),
  fOutCaloCells(0),
  fAddedCells(0),
  fMCParticles(0),
  fMCParticlesMap(0),
  fOutMCParticles(0),
  fOutMCParticlesMap(0),
  fMCLabelShift(0),
  fEsdMode(kFALSE),
  fOutput(0)
{
  // Default constructor.
}

//________________________________________________________________________
AliJetModelBaseTask::AliJetModelBaseTask(const char *name, Bool_t drawqa) : 
  AliAnalysisTaskSE(name),
  fGeomName(""),
  fTracksName("PicoTracks"),
  fOutTracksName("PicoTracksEmbedded"),
  fCaloName("CaloClustersCorr"),
  fOutCaloName("CaloClustersCorrEmbedded"),
  fCellsName(""),
  fOutCellsName(""),
  fMCParticlesName(""),
  fOutMCParticlesName(""),
  fSuffix("Processed"),
  fEtaMin(-1),
  fEtaMax(1),
  fPhiMin(0),
  fPhiMax(TMath::Pi() * 2),
  fPtMin(50),
  fPtMax(60),
  fCopyArray(kTRUE),
  fNClusters(0),
  fNCells(0),
  fNTracks(1),
  fMarkMC(99999),
  fPtSpectrum(0),
  fQAhistos(drawqa),
  fIsInit(0),
  fGeom(0),
  fClusters(0),
  fOutClusters(0),
  fTracks(0),
  fOutTracks(0),
  fCaloCells(0),
  fOutCaloCells(0),
  fAddedCells(0),
  fMCParticles(0),
  fMCParticlesMap(0),
  fOutMCParticles(0),
  fOutMCParticlesMap(0),
  fMCLabelShift(0),
  fEsdMode(kFALSE),
  fOutput(0)
{
  // Standard constructor.

  if (fQAhistos) {
    DefineOutput(1, TList::Class()); 
  }
}

//________________________________________________________________________
AliJetModelBaseTask::~AliJetModelBaseTask()
{
  // Destructor
}

//________________________________________________________________________
void AliJetModelBaseTask::UserCreateOutputObjects()
{
  // Create user output.
  if (!fQAhistos)
    return;

  OpenFile(1);
  fOutput = new TList();
  fOutput->SetOwner();

  PostData(1, fOutput);
}

//________________________________________________________________________
void AliJetModelBaseTask::UserExec(Option_t *) 
{
  // Execute per event.

  if (!fIsInit)
    fIsInit = ExecOnce();

  if (!fIsInit)
    return;

  if (fCopyArray) {
    if (fOutTracks)
      fOutTracks->Delete();
    if (fOutClusters)
      fOutClusters->Delete();
    if (fOutMCParticles)
      fOutMCParticles->Delete();
  }

  AliVCaloCells *tempCaloCells = 0;

  if (fCaloCells) {
    fAddedCells = 0;
    if (!fCopyArray) {
      tempCaloCells = fCaloCells;
      fCaloCells = static_cast<AliVCaloCells*>(tempCaloCells->Clone(Form("%s_old",fCaloCells->GetName())));
    }
  }

  Run();

  if (fCaloCells) {
    delete fCaloCells;
    fCaloCells = tempCaloCells;
  }
}

//________________________________________________________________________
Bool_t AliJetModelBaseTask::ExecOnce()
{
  // Init task.

  delete gRandom;
  gRandom = new TRandom3(0);

  fEsdMode = InputEvent()->InheritsFrom("AliESDEvent");

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

  if (!fCellsName.IsNull()) {
    fCaloCells = dynamic_cast<AliVCaloCells*>(InputEvent()->FindListObject(fCellsName));
    if (!fCaloCells) {
      AliWarning(Form("%s: Couldn't retrieve calo cells with name %s!", GetName(), fCellsName.Data()));
    }
    else if (!fCaloCells->InheritsFrom("AliVCaloCells")) {
      AliError(Form("%s: Collection %s does not contain a AliVCaloCells object!", GetName(), fCellsName.Data())); 
      fCaloCells = 0;
      return kFALSE;
    }

    if (!fOutCaloCells) {
      fOutCellsName = fCellsName;
      if (fCopyArray) 
	fOutCellsName += fSuffix;
      if (fCopyArray || !fCaloCells) {
	if (fEsdMode) 
	  fOutCaloCells = new AliESDCaloCells(fOutCellsName,fOutCellsName);
	else
	  fOutCaloCells = new AliAODCaloCells(fOutCellsName,fOutCellsName);

	if (InputEvent()->FindListObject(fOutCellsName)) {
	  AliFatal(Form("%s: Collection %s is already present in the event!", GetName(), fOutCellsName.Data())); 
	  return kFALSE;
	}
	else {
	  InputEvent()->AddObject(fOutCaloCells);
	}
      }
      else {
	fOutCaloCells = fCaloCells;
      }
    }
  }

  if (fNTracks > 0 && !fTracksName.IsNull()) {
    fTracks = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTracksName));
    if (!fTracks) {
      AliWarning(Form("%s: Couldn't retrieve tracks with name %s!", GetName(), fTracksName.Data()));
    }
    else if (!fTracks->GetClass()->GetBaseClass("AliPicoTrack")) {
      AliError(Form("%s: Collection %s does not contain AliPicoTrack objects!", GetName(), fTracksName.Data())); 
      fTracks = 0;
      return kFALSE;
    }

    if (!fOutTracks) {
      fOutTracksName = fTracksName;
      if (fCopyArray)
	fOutTracksName += fSuffix;
      if (fCopyArray || !fTracks) {
	fOutTracks = new TClonesArray("AliPicoTrack");
	fOutTracks->SetName(fOutTracksName);
	if (InputEvent()->FindListObject(fOutTracksName)) {
	  AliFatal(Form("%s: Collection %s is already present in the event!", GetName(), fOutTracksName.Data())); 
	  return kFALSE;
	}
	else {
	  InputEvent()->AddObject(fOutTracks);
	}
      }
      else {
	fOutTracks = fTracks;
      }
    }
  }

  if (fNClusters > 0 && !fCaloName.IsNull()) {
    fClusters = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fCaloName));
 
    if (!fClusters) {
      AliWarning(Form("%s: Couldn't retrieve clusters with name %s!", GetName(), fCaloName.Data()));
    }
    else if (!fClusters->GetClass()->GetBaseClass("AliVCluster")) {
      AliError(Form("%s: Collection %s does not contain AliVCluster objects!", GetName(), fCaloName.Data())); 
      fClusters = 0;
      return kFALSE;
    }

    if (!fOutClusters) {
      fOutCaloName = fCaloName;
      if (fCopyArray) 
	fOutCaloName += fSuffix;
      if (fCopyArray || !fClusters) {
	fOutClusters = new TClonesArray(fClusters->GetClass()->GetName(), fClusters->GetSize());
	fOutClusters->SetName(fOutCaloName);
	if (InputEvent()->FindListObject(fOutCaloName)) {
	  AliFatal(Form("%s: Collection %s is already present in the event!", GetName(), fOutCaloName.Data())); 
	  return kFALSE;
	}
	else {
	  InputEvent()->AddObject(fOutClusters);
	}
      }
      else {
	fOutClusters = fClusters;
      }
    }
  }

  if (!fMCParticlesName.IsNull()) {
    fMCParticles = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fMCParticlesName));
    if (!fMCParticles) {
      AliWarning(Form("%s: Couldn't retrieve MC particles with name %s!", GetName(), fMCParticlesName.Data()));
    }
    else {
      if (!fMCParticles->GetClass()->GetBaseClass("AliAODMCParticle")) {
	AliError(Form("%s: Collection %s does not contain AliAODMCParticle objects!", GetName(), fMCParticlesName.Data())); 
	fMCParticles = 0;
	return kFALSE;
      }
      
      fMCParticlesMap = dynamic_cast<TH1I*>(InputEvent()->FindListObject(fMCParticlesName + "_Map"));

      if (!fMCParticlesMap) {
	AliWarning(Form("%s: Could not retrieve map for MC particles %s! Will assume MC labels consistent with indexes...", GetName(), fMCParticlesName.Data())); 
	fMCParticlesMap = new TH1I(fMCParticlesName + "_Map", fMCParticlesName + "_Map",9999,0,1);
	for (Int_t i = 0; i < 9999; i++) {
	  fMCParticlesMap->SetBinContent(i,i);
	}
      }
    }

    if (!fOutMCParticles) {
      fOutMCParticlesName = fMCParticlesName;
      if (fCopyArray)
	fOutMCParticlesName += fSuffix;
      if (fCopyArray || !fMCParticles) {

	fOutMCParticles = new TClonesArray("AliAODMCParticle");
	fOutMCParticles->SetName(fOutMCParticlesName);
	if (InputEvent()->FindListObject(fOutMCParticlesName)) {
	  AliFatal(Form("%s: Collection %s is already present in the event!", GetName(), fOutMCParticlesName.Data())); 
	  return kFALSE;
	}
	else {
	  InputEvent()->AddObject(fOutMCParticles);
	}

	fOutMCParticlesMap = new TH1I(fOutMCParticlesName + "_Map", fOutMCParticlesName + "_Map",9999,0,1);
	if (InputEvent()->FindListObject(fOutMCParticlesName + "_Map")) {
	  AliFatal(Form("%s: Map %s_Map is already present in the event!", GetName(), fOutMCParticlesName.Data())); 
	  return kFALSE;
	}
	else {
	  InputEvent()->AddObject(fOutMCParticlesMap);
	}
      }
      else {
	fOutMCParticles = fMCParticles;
	fOutMCParticlesMap = fMCParticlesMap;
      }
    }
  }

  if (!fCaloName.IsNull() || !fCellsName.IsNull()) {
    if (!fGeom) {
      if (fGeomName.Length() > 0) {
	fGeom = AliEMCALGeometry::GetInstance(fGeomName);
	if (!fGeom) {
	  AliFatal(Form("Could not get geometry with name %s!", fGeomName.Data()));
	  return kFALSE;
	}
      } else {
	fGeom = AliEMCALGeometry::GetInstance();
	if (!fGeom) {
	  AliFatal("Could not get default geometry!");
	  return kFALSE;
	}
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

  return kTRUE;
}

//________________________________________________________________________
Int_t AliJetModelBaseTask::SetNumberOfOutCells(Int_t n)
{
  if (fOutCaloCells->GetNumberOfCells() < n) {
    fOutCaloCells->DeleteContainer();
    fOutCaloCells->CreateContainer(n);
  }
  else {
    fOutCaloCells->SetNumberOfCells(n);
  }

  fAddedCells = 0;

  return n;
}

//________________________________________________________________________
void AliJetModelBaseTask::CopyCells()
{
  if (!fCaloCells)
    return;

  for (Short_t i = 0; i < fCaloCells->GetNumberOfCells(); i++) {
    Int_t mclabel = 0;
    Double_t efrac = 0.;
    Double_t time = -1;
    Short_t cellNum = -1;
    Double_t amp = -1;

    fCaloCells->GetCell(i, cellNum, amp, time, mclabel, efrac);
    fOutCaloCells->SetCell(i, cellNum, amp, time, mclabel, efrac);
  }

  fAddedCells = fCaloCells->GetNumberOfCells();

  AliDebug(2, Form("%d cells from the current event", fAddedCells));
}

//________________________________________________________________________
Int_t AliJetModelBaseTask::AddCell(Double_t e, Double_t eta, Double_t phi)
{
  // Add a cell to the event.

  Int_t absId = 0;
  if (eta < -100 || phi < 0) {
    GetRandomCell(eta, phi, absId);
  }
  else {
    fGeom->EtaPhiFromIndex(absId, eta, phi);
  }

  if (absId == -1) {
    AliWarning(Form("Unable to embed cell in eta = %f, phi = %f!"
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

  return AddCell(e, absId);
}

//________________________________________________________________________
Int_t AliJetModelBaseTask::AddCell(Double_t e, Int_t absId, Double_t time, Int_t label)
{
  // Add a cell to the event.

  if (label == 0)
    label = fMarkMC;
  else
    label += fMCLabelShift;

  Bool_t r = fOutCaloCells->SetCell(fAddedCells, absId, e, time, label, 0);

  if (r) {
    fAddedCells++;
    return fAddedCells;
  }
  else {
    return 0;
  }
}


//________________________________________________________________________
AliVCluster* AliJetModelBaseTask::AddCluster(Double_t e, Double_t eta, Double_t phi, Int_t label)
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

  return AddCluster(e, absId, label);
}
      
//________________________________________________________________________
AliVCluster* AliJetModelBaseTask::AddCluster(Double_t e, Int_t absId, Int_t label)
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

  //MC label
  if (label == 0)
    label = fMarkMC;
  else
    label += fMCLabelShift;

  if (fEsdMode) {
    AliESDCaloCluster *esdClus = static_cast<AliESDCaloCluster*>(cluster);
    TArrayI parents(1, &label);
    esdClus->AddLabels(parents); 
  }
  else {
    AliAODCaloCluster *aodClus = static_cast<AliAODCaloCluster*>(cluster);
    aodClus->SetLabel(&label, 1); 
  }
  
  return cluster;
}

//________________________________________________________________________
AliPicoTrack* AliJetModelBaseTask::AddTrack(Double_t pt, Double_t eta, Double_t phi, Byte_t type, Double_t etaemc, Double_t phiemc, Bool_t ise, Int_t label)
{
  // Add a track to the event.

  const Int_t nTracks = fOutTracks->GetEntriesFast();
  
  if (pt < 0) 
    pt = GetRandomPt();
  if (eta < -100) 
    eta = GetRandomEta();
  if (phi < 0) 
    phi = GetRandomPhi();

  if (label == 0)
    label = fMarkMC;
  else
    label += fMCLabelShift;

  AliPicoTrack *track = new ((*fOutTracks)[nTracks]) AliPicoTrack(pt, 
								  eta, 
								  phi, 
								  1,
								  label,
								  type, 
								  etaemc, 
								  phiemc, 
								  ise);

  return track;
}

//________________________________________________________________________
AliAODMCParticle* AliJetModelBaseTask::AddMCParticle(AliAODMCParticle *part, Int_t origIndex)
{
  const Int_t nPart = fOutMCParticles->GetEntriesFast();

  AliAODMCParticle *aodpart = new ((*fOutMCParticles)[nPart]) AliAODMCParticle(*part);

  fOutMCParticlesMap->SetBinContent(origIndex + fMCLabelShift, nPart);

  return aodpart;
}

//________________________________________________________________________
void AliJetModelBaseTask::CopyClusters()
{
  // Copy all the clusters in the new collection
  if (!fClusters)
    return;

  const Int_t nClusters = fClusters->GetEntriesFast();
  Int_t nCopiedClusters = 0;
  
  if (fEsdMode) {
    for (Int_t i = 0; i < nClusters; ++i) {
      AliESDCaloCluster *esdcluster = static_cast<AliESDCaloCluster*>(fClusters->At(i));
      if (!esdcluster || !esdcluster->IsEMCAL())
	continue;
      new ((*fOutClusters)[nCopiedClusters]) AliESDCaloCluster(*esdcluster);
      nCopiedClusters++;
    }
  }
  else {
    for (Int_t i = 0; i < nClusters; ++i) {
      AliAODCaloCluster *aodcluster = static_cast<AliAODCaloCluster*>(fClusters->At(i));
      if (!aodcluster || !aodcluster->IsEMCAL())
	continue;
      new ((*fOutClusters)[nCopiedClusters]) AliAODCaloCluster(*aodcluster);
      nCopiedClusters++;
    }
  }
}

//________________________________________________________________________
void AliJetModelBaseTask::CopyTracks()
{
  // Copy all the tracks in the new collection

  if (!fTracks)
    return;

  const Int_t nTracks = fTracks->GetEntriesFast();
  Int_t nCopiedTracks = 0;
  for (Int_t i = 0; i < nTracks; ++i) {
    AliPicoTrack *track = static_cast<AliPicoTrack*>(fTracks->At(i));
    if (!track)
      continue;
    new ((*fOutTracks)[nCopiedTracks]) AliPicoTrack(*track);
    nCopiedTracks++;
  }
}

//________________________________________________________________________
void AliJetModelBaseTask::CopyMCParticles()
{
  // Copy all the MC particles in the new collection

  if (!fMCParticles)
    return;

  const Int_t nPart = fMCParticles->GetEntriesFast();
  Int_t nCopiedPart = 0;
  for (Int_t i = 0; i < nPart; ++i) {
    AliAODMCParticle *part = static_cast<AliAODMCParticle*>(fMCParticles->At(i));
    if (!part)
      continue;
    new ((*fOutMCParticles)[nCopiedPart]) AliAODMCParticle(*part);

    nCopiedPart++;
  }

  if (!fMCParticlesMap)
    return;

  Int_t shift = 0;

  for (Int_t i = 0; i < fMCParticlesMap->GetNbinsX()+2; i++) {
    fOutMCParticlesMap->SetBinContent(i, fMCParticlesMap->GetBinContent(i));
    if (fMCParticlesMap->GetBinContent(i) != 0)
      shift = i;
  }

  fMCLabelShift = shift;
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
