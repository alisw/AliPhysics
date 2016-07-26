//
// Jet modeling task.
//
// Author: S.Aiola, C.Loizides

#include "AliJetModelBaseTask.h"

#include <TClonesArray.h>
#include <TH1I.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TList.h>
#include <TF2.h>
#include <TGrid.h>
#include <TFile.h>

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
  fPythiaInfoName(""),
  fIsMC(kFALSE),
  fSuffix(),
  fEtaMin(-1),
  fEtaMax(1),
  fPhiMin(0),
  fPhiMax(TMath::Pi() * 2),
  fPtMin(0),
  fPtMax(0),
  fGenType(0),
  fCopyArray(kTRUE),
  fNClusters(0),
  fNCells(0),
  fNTracks(0),
  fMarkMC(99999),
  fPtSpectrum(0),
  fPtPhiEvPlDistribution(0),
  fDensitySpectrum(0),
  fDifferentialV2(0),
  fAddV2(kFALSE),
  fFlowFluctuations(kFALSE),
  fQAhistos(kFALSE),
  fPsi(0),
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
  fOutput(0),
  fPythiaInfo(0x0),
  fhpTEmb(0),
  fhMEmb(0),
  fhEtaEmb(0),
  fhPhiEmb(0),
  fMassFromDistr(kFALSE),
  fHMassDistrib(0),
  fHMassPtDistrib(0)
{
  // Default constructor.

  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;
}

//________________________________________________________________________
AliJetModelBaseTask::AliJetModelBaseTask(const char *name, Bool_t drawqa) : 
  AliAnalysisTaskSE(name),
  fGeomName("EMCAL_COMPLETE12SMV1"),
  fTracksName("PicoTracks"),
  fOutTracksName("PicoTracksEmbedded"),
  fCaloName("CaloClustersCorr"),
  fOutCaloName("CaloClustersCorrEmbedded"),
  fCellsName(""),
  fOutCellsName(""),
  fMCParticlesName(""),
  fOutMCParticlesName(""),
  fPythiaInfoName(""),
  fIsMC(kFALSE),
  fSuffix("Processed"),
  fEtaMin(-1),
  fEtaMax(1),
  fPhiMin(0),
  fPhiMax(TMath::Pi() * 2),
  fPtMin(50),
  fPtMax(60),
  fGenType(0),
  fCopyArray(kTRUE),
  fNClusters(0),
  fNCells(0),
  fNTracks(1),
  fMarkMC(99999),
  fPtSpectrum(0),
  fPtPhiEvPlDistribution(0),
  fDensitySpectrum(0),
  fDifferentialV2(0),
  fAddV2(kFALSE),
  fFlowFluctuations(kFALSE),
  fQAhistos(drawqa),
  fPsi(0),
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
  fOutput(0),
  fPythiaInfo(0x0),
  fhpTEmb(0),
  fhMEmb(0),
  fhEtaEmb(0),
  fhPhiEmb(0),
  fMassFromDistr(kFALSE),
  fHMassDistrib(0),
  fHMassPtDistrib(0)
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
TString AliJetModelBaseTask::GetOutTrackName() const{
   Printf("Note: If the code changes this method could give a wrong result");
   TString futurenamefOutputTrack;
   if(fCopyArray) futurenamefOutputTrack = Form("%s%s", fTracksName.Data(), fSuffix.Data());
   else futurenamefOutputTrack = fTracksName;
   return futurenamefOutputTrack;

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

  fhpTEmb  = new TH1F("fhpTEmb","#it{p}_{T} distribution; #it{p}_{T}(GeV/c)", 120, 0., 120);
  fhpTEmb->Sumw2();
  fOutput->Add(fhpTEmb);
  
  fhMEmb   = new TH1F("fhMEmb","Mass distribution; #it{M} (GeV)", 80, 0, 80.);
  fhMEmb->Sumw2();
  fOutput->Add(fhMEmb);
  
  fhEtaEmb = new TH1F("fhEtaEmb","#eta distribution; #eta", 100, -1, 1);
  fhEtaEmb->Sumw2();
  fOutput->Add(fhEtaEmb);
  
  fhPhiEmb = new TH1F("fhPhiEmb","#varphi distribution; #varphi", 100, (-1)*TMath::Pi(), TMath::Pi());
  fhPhiEmb->Sumw2();
  fOutput->Add(fhPhiEmb);
  
  fhEvents = new TH1I("fhEvents", "Number of events", 3, 0, 2);
  fOutput->Add(fhEvents);

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
    fNTracks = TMath::Nint(fDensitySpectrum->GetRandom());
    fNCells = TMath::Nint(fDensitySpectrum->GetRandom());
    fNClusters = TMath::Nint(fDensitySpectrum->GetRandom());
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

  if (fPtPhiEvPlDistribution || fAddV2)
    fPsi = gRandom->Rndm() * TMath::Pi();
  
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
	//fOutTracks->SetOwner(kTRUE);
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
	InputEvent()->AddObject(fOutTracks);
      }
    }
    InputEvent()->ls();
  }

  if(fAddV2 && (!fDifferentialV2)) {
    AliWarning(Form("%s: Cannot add v2 without diffential v2!", GetName()));
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
AliPicoTrack* AliJetModelBaseTask::AddTrack(Double_t pt, Double_t eta, Double_t phi, Byte_t type, Double_t etaemc, Double_t phiemc, Double_t ptemc, Bool_t ise, Int_t label, Short_t charge, Double_t mass)
{
  // Add a track to the event.
  if (pt < 0 && eta < -100 && phi < -100 && mass < -100) {
     GetRandomMassiveParticle(pt,eta,phi, kFALSE, mass);
  } else {
     if (pt < 0 && eta < -100 && phi < -100) {
     	GetRandomParticle(pt,eta,phi);
     }
     else {
     	if (pt < -100) 
     	   pt = GetRandomPt();
     	if (eta < -100) 
     	   eta = GetRandomEta();
     	if (phi < -100) 
     	   phi = GetRandomPhi(pt);
     }
  }
  //Printf("Adding LABEL %d", label);
  if (label >= 0)
    label += fMarkMC+fMCLabelShift;
  else if (label < 0)
    label -= fMarkMC+fMCLabelShift;
  if(fAddV2) AddV2(phi, pt);

  const Int_t nTracks = fOutTracks->GetEntriesFast();
  //Printf("+ %d = %d", fMarkMC, label);
  AliPicoTrack *track = new ((*fOutTracks)[nTracks]) AliPicoTrack(pt, 
								  eta, 
								  phi, 
								  charge,
								  label,
								  type, 
								  etaemc, 
								  phiemc,
								  ptemc,
								  ise,
								  mass);

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

//_____________________________________________________________________________
void AliJetModelBaseTask::AddV2(Double_t &phi, Double_t &pt) const
{
    // similar to AliFlowTrackSimple::AddV2, except for the flow fluctuations
    Double_t phi0(phi), v2(0.), f(0.), fp(0.), phiprev(0.);
    if(fDifferentialV2) v2 = fDifferentialV2->Eval(pt);
    if(TMath::AreEqualAbs(v2, 0, 1e-5)) return; 
    // introduce flow fluctuations (gaussian)
    if(fFlowFluctuations) v2 += TMath::Sqrt(2*(v2*.25)*(v2*.25))*TMath::ErfInverse(2*(gRandom->Uniform(0, fFlowFluctuations))-1);
    for (Int_t i(0); i < 100; i++) {
        phiprev=phi; //store last value for comparison
        f =  phi-phi0+v2*TMath::Sin(2.*(phi-fPsi));
        fp = 1.0+2.0*v2*TMath::Cos(2.*(phi-fPsi)); //first derivative
        phi -= f/fp;
        if (TMath::AreEqualAbs(phiprev, phi, 1e-10)) break; 
    }
}

//_____________________________________________________________________________
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

  Double_t result = gRandom->Rndm() * (phimax - phimin) + phimin;

  return result;
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

Double_t AliJetModelBaseTask::GetRandomM(){
   
   Double_t m = 0;
   
   if(fMassFromDistr) {
      if(fHMassDistrib) m = fHMassDistrib->GetRandom();
      else {
      	 AliError("Template distribution for mass of track embedding not found, use 0");
      	 m = 0;
      }
   }
   return m;
}

//________________________________________________________________________

void AliJetModelBaseTask::GetRandomMvsPt(Double_t &m, Double_t &pt){
   
   Double_t maxedgex = fHMassPtDistrib->GetXaxis()->GetBinLowEdge(fHMassPtDistrib->GetNbinsX());
   Double_t maxedgey = fHMassPtDistrib->GetYaxis()->GetBinLowEdge(fHMassPtDistrib->GetNbinsY());
   
   if(maxedgex < maxedgey) //condition assuming that the mass axis has smaller range than the pT axis!!
      fHMassPtDistrib->GetRandom2(m, pt);
   else fHMassPtDistrib->GetRandom2(pt, m);

}

//________________________________________________________________________
void AliJetModelBaseTask::GetRandomParticle(Double_t &pt, Double_t &eta, Double_t &phi, Bool_t emcal)
{
  // Get a random particle.
  
  eta = GetRandomEta(emcal);

  if (fPtPhiEvPlDistribution) {
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
    
    if (fPtPhiEvPlDistribution->GetXmin() > phimax || fPtPhiEvPlDistribution->GetXmax() < phimin) {
      AliWarning(Form("The hisogram %s does not overlap with the EMCal acceptance limits. It will be ignored.",fPtPhiEvPlDistribution->GetName()));
      pt = GetRandomPt();
      phi = GetRandomPhi(emcal);
    }
    else {
      do {
	fPtPhiEvPlDistribution->GetRandom2(pt,phi);
	phi += fPsi;
	if (phi > TMath::Pi() * 2) phi -= TMath::Pi() * 2;
      } while (phi > phimax || phi < phimin);
    }
  }
  else {
    pt = GetRandomPt();
    phi = GetRandomPhi(emcal);
  }
}

//________________________________________________________________________
void AliJetModelBaseTask::GetRandomMvsPtParticle(Double_t &pt, Double_t &m, Double_t &eta, Double_t &phi, Bool_t emcal){
   
   /// Get random particle from 2D mass vs pt distribution
   /// the event plane evolution in not implemented
   
   phi = GetRandomPhi(emcal);
   eta = GetRandomEta(emcal);
   GetRandomMvsPt(m, pt);

}
//________________________________________________________________________
void AliJetModelBaseTask::GetRandomMassiveParticle(Double_t &pt, Double_t &eta, Double_t &phi, Bool_t emcal, Double_t& m){
   
   if(!fHMassPtDistrib) {
   GetRandomParticle(pt,eta,phi,emcal);
   m = GetRandomM();
   
   } else GetRandomMvsPtParticle(pt,m,eta,phi,emcal);

}

//________________________________________________________________________
void AliJetModelBaseTask::Run() 
{
  // Run.
}
//________________________________________________________________________
void AliJetModelBaseTask::FillHistograms(){
   
   if(!fhpTEmb || !fhMEmb || !fhEtaEmb || !fhPhiEmb) {
      AliError("Histograms not found, are the QA histograms active?");
   }
   fhEvents->Fill(0);
   // fill histograms
   Int_t nentries = fOutTracks->GetEntries();
   for(Int_t it = 0; it<fNTracks; it++){
      AliPicoTrack *trackEmb = dynamic_cast<AliPicoTrack*>(fOutTracks->At(nentries-it-1));
      if(!trackEmb) continue;
      if(trackEmb->GetLabel() >= fMarkMC+fMCLabelShift){
      	 fhpTEmb ->Fill(trackEmb->Pt());
      	 fhMEmb  ->Fill(trackEmb->M());
      	 fhEtaEmb->Fill(trackEmb->Eta());
      	 fhPhiEmb->Fill(trackEmb->Phi());
      }
   }
   
   PostData(1, fOutput);

}

//________________________________________________________________________

void AliJetModelBaseTask::SetMassDistribution(TH1F *hM)  {
   if(!hM){
      AliError("Null histogram for mass distribution");
      return;
   }
   fMassFromDistr = kTRUE; 
   fHMassDistrib = hM;
   AliInfo("Input mass distribution set");
   
   return;
}

//________________________________________________________________________

void AliJetModelBaseTask::SetMassVsPtDistribution(TH2F *hmasspt)  {
   if(!hmasspt){
      AliError("Null histogram for mass vs pt distribution");
      return;
   }
   fMassFromDistr = kTRUE; 
   fHMassPtDistrib = hmasspt;
   AliInfo("Input mass vs pt distribution set");
   
   return;
}

//________________________________________________________________________
void AliJetModelBaseTask::SetpTDistributionFromFile(TString filename, TString histoname){
   
 SetDistributionFromFile(filename, histoname, 1);
}

//________________________________________________________________________
void AliJetModelBaseTask::SetMassDistributionFromFile(TString filename, TString histoname){
   SetDistributionFromFile(filename, histoname, 2);
}

//________________________________________________________________________
void AliJetModelBaseTask::SetMassVsPtDistributionFromFile(TString filename, TString histoname){
   SetDistributionFromFile(filename, histoname, 3);
}

//________________________________________________________________________
void AliJetModelBaseTask::SetDistributionFromFile(TString filename, TString histoname, Int_t type){
   
   if(filename.Contains("alien")) {
      TGrid::Connect("alien://");
   }
   TFile *f = TFile::Open(filename);
   if(!f){
      AliFatal(Form("File %s not found, cannot SetMassDistribution", filename.Data()));
      return;
   }
   
   TH1F* h = 0x0;
   TH2F* g = 0x0;
   if(type < 3){
      h = dynamic_cast<TH1F*> (f->Get(histoname));
      if(!h) {
      	 AliError("Input file for Mass not found");
      	 f->ls();
      }
   }
   if(type == 3){
      g = dynamic_cast<TH2F*> (f->Get(histoname));
      if(!g) {
      	 AliError("Input file for Mass not found");
      	 f->ls();
      }
   }
   
   if(type == 1) SetPtSpectrum(h);
   if(type == 2) SetMassDistribution(h);
   if(type == 3) SetMassVsPtDistribution(g);
   
   //f->Close();
   //delete f;
   
   return;

}

//________________________________________________________________________

void AliJetModelBaseTask::SetMassAndPtDistributionFromFile(TString filenameM, TString filenamepT, TString histonameM, TString histonamepT){
   SetMassDistributionFromFile(filenameM, histonameM);
   SetpTDistributionFromFile(filenamepT, histonamepT);
   return;
}


