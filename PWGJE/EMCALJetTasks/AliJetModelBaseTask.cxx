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
#include "AliNamedArrayI.h"

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
  fIsMC(kFALSE),
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
  fDensitySpectrum(0),
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

  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;
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
  fIsMC(kFALSE),
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
  fDensitySpectrum(0),
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

  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;
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

  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;

  const AliVVertex *vert = InputEvent()->GetPrimaryVertex();
  if (vert)
    vert->GetXYZ(fVertex);

  if (fCopyArray) {
    if (fOutTracks)
      fOutTracks->Delete();
    if (fOutClusters)
      fOutClusters->Delete();
    if (fOutMCParticles)
      fOutMCParticles->Delete();
  }

  if (fDensitySpectrum) {
    fNTracks = fDensitySpectrum->GetRandom();
    fNCells = fDensitySpectrum->GetRandom();
    fNClusters = fDensitySpectrum->GetRandom();
  }
  
  // Clear map
  if (fOutMCParticlesMap)
    fOutMCParticlesMap->Clear();

  AliVCaloCells *tempCaloCells = 0;

  if (fCaloCells) {
    fAddedCells = 0;
    if (!fCopyArray) {
      tempCaloCells = fCaloCells;
      fCaloCells = static_cast<AliVCaloCells*>(tempCaloCells->Clone(Form("%s_old",fCaloCells->GetName())));
    }
  }

  Run();

  if (fCaloCells && !fCopyArray) {
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

  if (!fTracksName.IsNull()) {
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

  if (!fCaloName.IsNull()) {
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
      TString className;
      if (fClusters)
	className = fClusters->GetClass()->GetName();
      else if (fEsdMode)
	className = "AliESDCaloCluster";
      else
	className = "AliAODCaloCluster";
	
      if (fCopyArray || !fClusters) {
	fOutClusters = new TClonesArray(className.Data());
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
      
      fMCParticlesMap = dynamic_cast<AliNamedArrayI*>(InputEvent()->FindListObject(fMCParticlesName + "_Map"));

      if (!fMCParticlesMap) {
	AliWarning(Form("%s: Could not retrieve map for MC particles %s! Will assume MC labels consistent with indexes...", GetName(), fMCParticlesName.Data())); 
	fMCParticlesMap = new AliNamedArrayI(fMCParticlesName + "_Map", 99999);
	for (Int_t i = 0; i < 99999; i++) {
	  fMCParticlesMap->AddAt(i,i);
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

	fOutMCParticlesMap = new AliNamedArrayI(fOutMCParticlesName + "_Map",99999);
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

  if (label < 0)
    label = 0;

  label += fMarkMC + fMCLabelShift;

  Short_t pos = -1;
  if (fCaloCells)  
    pos = fCaloCells->GetCellPosition(absId);

  Double_t efrac = 1;
  Bool_t increaseOnSuccess = kFALSE;

  if (pos < 0) {
    increaseOnSuccess = kTRUE;
    pos = fAddedCells;
  }
  else {
    Short_t cellNumber = -1;
    Double_t old_e = 0;
    Double_t old_time = 0;
    Int_t old_label = 0;
    Double_t old_efrac = 0;
    fOutCaloCells->GetCell(pos, cellNumber, old_e, old_time, old_label, old_efrac);
    
    efrac = e / (old_e + e);

    if (old_label > 0 && e < old_efrac * old_e) {
      label = old_label;
      efrac = old_efrac;
      time = old_time;
    }
    
    e += old_e;
  }

  Bool_t r = fOutCaloCells->SetCell(pos, absId, e, time, label, efrac);

  if (r) {
    if (increaseOnSuccess)
      fAddedCells++;
    return fAddedCells;
  }
  else 
    return 0;
}

//________________________________________________________________________
AliVCluster* AliJetModelBaseTask::AddCluster(AliVCluster *oc)
{
  // Add a cluster to the event.

  const Int_t nClusters = fOutClusters->GetEntriesFast();
  AliVCluster *dc = static_cast<AliVCluster*>(fOutClusters->New(nClusters));
  dc->SetType(AliVCluster::kEMCALClusterv1);
  dc->SetE(oc->E());
  Float_t pos[3] = {0};
  oc->GetPosition(pos);
  dc->SetPosition(pos);
  dc->SetNCells(oc->GetNCells());
  dc->SetCellsAbsId(oc->GetCellsAbsId());
  dc->SetCellsAmplitudeFraction(oc->GetCellsAmplitudeFraction());
  dc->SetID(oc->GetID());
  dc->SetDispersion(oc->GetDispersion());
  dc->SetEmcCpvDistance(-1);
  dc->SetChi2(-1);
  dc->SetTOF(oc->GetTOF());     //time-of-flight
  dc->SetNExMax(oc->GetNExMax()); //number of local maxima
  dc->SetM02(oc->GetM02());
  dc->SetM20(oc->GetM20());
  dc->SetDistanceToBadChannel(oc->GetDistanceToBadChannel()); 

  //MC
  UInt_t nlabels = oc->GetNLabels();
  Int_t *labels = oc->GetLabels();
  TArrayI parents;
  if (nlabels != 0 && labels) 
    parents.Set(nlabels, labels);
  else {
    nlabels = 1;
    parents.Set(1);
    parents[0] = 0;
  }

  if (fMarkMC+fMCLabelShift != 0) {
    for (UInt_t i = 0; i < nlabels; i++) {
      parents[i] += fMarkMC+fMCLabelShift;
    }
  }

  AliESDCaloCluster *esdClus = dynamic_cast<AliESDCaloCluster*>(dc);
  if (esdClus) {
    esdClus->AddLabels(parents); 
  }
  else {
    AliAODCaloCluster *aodClus = dynamic_cast<AliAODCaloCluster*>(dc);
    if (aodClus) 
      aodClus->SetLabel(parents.GetArray(), nlabels); 
  }

  return dc;
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
  if (label < 0)
    label = 0;
 
  label += fMarkMC+fMCLabelShift;

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
AliPicoTrack* AliJetModelBaseTask::AddTrack(Double_t pt, Double_t eta, Double_t phi, Byte_t type, Double_t etaemc, Double_t phiemc, Double_t ptemc, Bool_t ise, Int_t label, Short_t charge)
{
  // Add a track to the event.

  const Int_t nTracks = fOutTracks->GetEntriesFast();
  
  if (pt < 0) 
    pt = GetRandomPt();
  if (eta < -100) 
    eta = GetRandomEta();
  if (phi < 0) 
    phi = GetRandomPhi();

  if (label >= 0)
    label += fMarkMC+fMCLabelShift;
  else if (label < 0)
    label -= fMarkMC+fMCLabelShift;

  AliPicoTrack *track = new ((*fOutTracks)[nTracks]) AliPicoTrack(pt, 
								  eta, 
								  phi, 
								  charge,
								  label,
								  type, 
								  etaemc, 
								  phiemc,
								  ptemc,
								  ise);

  return track;
}

//________________________________________________________________________
AliAODMCParticle* AliJetModelBaseTask::AddMCParticle(AliAODMCParticle *part, Int_t origIndex)
{
  const Int_t nPart = fOutMCParticles->GetEntriesFast();

  AliAODMCParticle *aodpart = new ((*fOutMCParticles)[nPart]) AliAODMCParticle(*part);

  if (origIndex + fMCLabelShift >= fOutMCParticlesMap->GetSize())
    fOutMCParticlesMap->Set((origIndex + fMCLabelShift)*2);

  fOutMCParticlesMap->AddAt(nPart, origIndex + fMCLabelShift);
  AliDebug(2, Form("Setting bin %d to %d (fMCLabelShift=%d, origIndex=%d)", 
		   origIndex + fMCLabelShift, fOutMCParticlesMap->At(origIndex + fMCLabelShift), fMCLabelShift, origIndex));

  return aodpart;
}

//________________________________________________________________________
void AliJetModelBaseTask::CopyCells()
{
  if (!fCaloCells)
    return;

  fAddedCells = 0;
  fCaloCells->Sort();
  for (Short_t i = 0; i < fCaloCells->GetNumberOfCells(); i++) {
    Int_t mclabel = 0;
    Double_t efrac = 0.;
    Double_t time = -1;
    Short_t cellNum = -1;
    Double_t amp = -1;

    fCaloCells->GetCell(i, cellNum, amp, time, mclabel, efrac);

    if (!fIsMC) 
      mclabel = 0;

    fOutCaloCells->SetCell(i, cellNum, amp, time, mclabel, efrac);
    fAddedCells++;
  }

  AliDebug(2, Form("%d cells from the current event", fAddedCells));
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
      AliESDCaloCluster *clus = new ((*fOutClusters)[nCopiedClusters]) AliESDCaloCluster(*esdcluster);
      if (!fIsMC) {
	TArrayI *labels = clus->GetLabelsArray();
	if (labels)
	  labels->Reset();
      }
      nCopiedClusters++;
    }
  }
  else {
    for (Int_t i = 0; i < nClusters; ++i) {
      AliAODCaloCluster *aodcluster = static_cast<AliAODCaloCluster*>(fClusters->At(i));
      if (!aodcluster || !aodcluster->IsEMCAL())
	continue;
      AliAODCaloCluster *clus = new ((*fOutClusters)[nCopiedClusters]) AliAODCaloCluster(*aodcluster);
      if (!fIsMC) 
	clus->SetLabel(0,0);
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
    AliPicoTrack *picotrack = static_cast<AliPicoTrack*>(fTracks->At(i));
    if (!picotrack)
      continue;
    AliPicoTrack *track = new ((*fOutTracks)[nCopiedTracks]) AliPicoTrack(*picotrack);
    if (!fIsMC && track->GetLabel() != 0) 
      track->SetLabel(0);
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

  if (!fMCParticlesMap || !fOutMCParticlesMap)
    return;

  if (fOutMCParticlesMap->GetSize() < fMCParticlesMap->GetSize())
    fOutMCParticlesMap->Set(fMCParticlesMap->GetSize() * 2);

  for (Int_t i = 0; i < fMCParticlesMap->GetSize(); i++) {
    fOutMCParticlesMap->AddAt(fMCParticlesMap->At(i), i);
    if (fMCParticlesMap->At(i) >= 0)
      fMCLabelShift = i;
  }

  AliDebug(2,Form("MC particles copied. fMCLabelShift=%d",fMCLabelShift));
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
      rndEta = GetRandomEta(kTRUE);
    if (phi < 0)
      rndPhi = GetRandomPhi(kTRUE);
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
Double_t AliJetModelBaseTask::GetRandomEta(Bool_t emcal)
{
  // Get random eta.

  Double_t etamax = fEtaMax;
  Double_t etamin = fEtaMin;

  if (emcal) {
    const Double_t EmcalMinEta = fGeom->GetArm1EtaMin();
    const Double_t EmcalMaxEta = fGeom->GetArm1EtaMax();
    
    if (etamax > EmcalMaxEta) etamax = EmcalMaxEta;
    if (etamax < EmcalMinEta) etamax = EmcalMinEta;
    if (etamin > EmcalMaxEta) etamin = EmcalMaxEta;
    if (etamin < EmcalMinEta) etamin = EmcalMinEta;
  }

  return gRandom->Rndm() * (etamax - etamin) + etamin;
}

//________________________________________________________________________
Double_t AliJetModelBaseTask::GetRandomPhi(Bool_t emcal)
{
  // Get random phi.
  
  Double_t phimax = fPhiMax;
  Double_t phimin = fPhiMin;

  if (emcal) {
    const Double_t EmcalMinPhi = fGeom->GetArm1PhiMin() * TMath::DegToRad();
    const Double_t EmcalMaxPhi = fGeom->GetArm1PhiMax() * TMath::DegToRad();
    
    if (phimax > EmcalMaxPhi) phimax = EmcalMaxPhi;
    if (phimax < EmcalMinPhi) phimax = EmcalMinPhi;
    if (phimin > EmcalMaxPhi) phimin = EmcalMaxPhi;
    if (phimin < EmcalMinPhi) phimin = EmcalMinPhi;
  }

  return gRandom->Rndm() * (phimax - phimin) + phimin;
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
