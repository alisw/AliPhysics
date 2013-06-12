// $Id$
//
// Emcal base analysis task.
//
// Author: S.Aiola, M. Verweij

#include "AliAnalysisTaskEmcalDev.h"

#include <TClonesArray.h>
#include <TList.h>
#include <TObject.h>
#include <TH1F.h>

#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliEMCALGeometry.h"
#include "AliESDEvent.h"
#include "AliEmcalParticle.h"
#include "AliEventplane.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMCParticle.h"
#include "AliVCluster.h"
#include "AliVEventHandler.h"
#include "AliVParticle.h"

#include "AliParticleContainer.h"
#include "AliClusterContainer.h"

ClassImp(AliAnalysisTaskEmcalDev)

//________________________________________________________________________
AliAnalysisTaskEmcalDev::AliAnalysisTaskEmcalDev() : 
  AliAnalysisTaskSE("AliAnalysisTaskEmcalDev"),
  fAnaType(kTPC),
  fForceBeamType(kNA),
  fGeneralHistograms(kFALSE),
  fInitialized(kFALSE),
  fCreateHisto(kTRUE),
  fTracksName(),
  fCaloName(),
  fCaloCellsName(),
  fMinCent(-999),
  fMaxCent(-999),
  fMinVz(-999),
  fMaxVz(-999),
  fOffTrigger(AliVEvent::kAny),
  fTrigClass(),
  fNbins(500),
  fMinBinPt(0),
  fMaxBinPt(250),
  fClusPtCut(0.15),
  fTrackPtCut(0.15),
  fTrackMinEta(-0.9),
  fTrackMaxEta(0.9),
  fTrackMinPhi(-10),
  fTrackMaxPhi(10),
  fClusTimeCutLow(-10),
  fClusTimeCutUp(10),
  fMinPtTrackInEmcal(0),
  fEventPlaneVsEmcal(-1),
  fMinEventPlane(-10),
  fMaxEventPlane(10),
  fCentEst("V0M"),
  fTrackBitMap(0),
  fClusterBitMap(0),
  fMCTrackBitMap(0),
  fMCClusterBitMap(0),
  fMinMCLabel(0),
  fNcentBins(4),
  fGeom(0),
  fTracks(0),
  fCaloClusters(0),
  fCaloCells(0),
  fCent(0),
  fCentBin(-1),
  fEPV0(-1.0),
  fEPV0A(-1.0),
  fEPV0C(-1.0),
  fNVertCont(0),
  fBeamType(kNA),
  fParticleCollArray(),
  fClusterCollArray(),
  fOutput(0),
  fHistCentrality(0),
  fHistZVertex(0),
  fHistEventPlane(0)
{
  // Default constructor.

  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;

  fParticleCollArray.SetOwner(kTRUE);
  fClusterCollArray.SetOwner(kTRUE);

}

//________________________________________________________________________
AliAnalysisTaskEmcalDev::AliAnalysisTaskEmcalDev(const char *name, Bool_t histo) : 
  AliAnalysisTaskSE(name),
  fAnaType(kTPC),
  fForceBeamType(kNA),
  fGeneralHistograms(kFALSE),
  fInitialized(kFALSE),
  fCreateHisto(histo),
  fTracksName(),
  fCaloName(),
  fCaloCellsName(),
  fMinCent(-999),
  fMaxCent(-999),
  fMinVz(-999),
  fMaxVz(-999),
  fOffTrigger(AliVEvent::kAny),
  fTrigClass(),
  fNbins(500),
  fMinBinPt(0),
  fMaxBinPt(250),
  fClusPtCut(0.15),
  fTrackPtCut(0.15),
  fTrackMinEta(-0.9),
  fTrackMaxEta(0.9),
  fTrackMinPhi(-10),
  fTrackMaxPhi(10),
  fClusTimeCutLow(-10),
  fClusTimeCutUp(10),
  fMinPtTrackInEmcal(0),
  fEventPlaneVsEmcal(-1),
  fMinEventPlane(-10),
  fMaxEventPlane(10),
  fCentEst("V0M"),
  fTrackBitMap(0),
  fClusterBitMap(0),
  fMCTrackBitMap(0),
  fMCClusterBitMap(0),
  fMinMCLabel(0),
  fNcentBins(4),
  fGeom(0),
  fTracks(0),
  fCaloClusters(0),
  fCaloCells(0),
  fCent(0),
  fCentBin(-1),
  fEPV0(-1.0),
  fEPV0A(-1.0),
  fEPV0C(-1.0),
  fNVertCont(0),
  fBeamType(kNA),
  fParticleCollArray(),
  fClusterCollArray(),
  fOutput(0),
  fHistCentrality(0),
  fHistZVertex(0),
  fHistEventPlane(0)
{
  // Standard constructor.

  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;

  fParticleCollArray.SetOwner(kTRUE);
  fClusterCollArray.SetOwner(kTRUE);

  if (fCreateHisto) {
    DefineOutput(1, TList::Class()); 
  }
}

//________________________________________________________________________
AliAnalysisTaskEmcalDev::~AliAnalysisTaskEmcalDev()
{
  // Destructor
}

//________________________________________________________________________
void AliAnalysisTaskEmcalDev::UserCreateOutputObjects()
{
  // Create user output.
  if (!fCreateHisto)
    return;

  OpenFile(1);
  fOutput = new TList();
  fOutput->SetOwner();

  if (fForceBeamType == kpp)
    fNcentBins = 1;

  if (!fGeneralHistograms)
    return;

  fHistCentrality = new TH1F("fHistCentrality","Event centrality distribution", 200, 0, 100);
  fHistCentrality->GetXaxis()->SetTitle("Centrality (%)");
  fHistCentrality->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistCentrality);

  fHistZVertex = new TH1F("fHistZVertex","Z vertex position", 60, -30, 30);
  fHistZVertex->GetXaxis()->SetTitle("z");
  fHistZVertex->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistZVertex);

  fHistEventPlane = new TH1F("fHistEventPlane","Event plane", 120, -TMath::Pi(), TMath::Pi());
  fHistEventPlane->GetXaxis()->SetTitle("event plane");
  fHistEventPlane->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistEventPlane);

  PostData(1, fOutput);
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalDev::FillGeneralHistograms()
{
  fHistCentrality->Fill(fCent);
  fHistZVertex->Fill(fVertex[2]);
  fHistEventPlane->Fill(fEPV0);

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalDev::UserExec(Option_t *) 
{
  // Main loop, called for each event.

  if (!fInitialized)
    ExecOnce();

  if (!fInitialized)
    return;

  if (!RetrieveEventObjects())
    return;

  if (!IsEventSelected()) 
    return;

  if (fGeneralHistograms && fCreateHisto) {
    if (!FillGeneralHistograms())
      return;
  }

  if (!Run())
    return;

  if (fCreateHisto) {
    if (!FillHistograms())
      return;
  }
    
  if (fCreateHisto && fOutput) {
    // information for this iteration of the UserExec in the container
    PostData(1, fOutput);
  }
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalDev::AcceptCluster(AliVCluster *clus, const Int_t c) const
{
  // Return true if cluster is accepted.

  if (!clus)
    return kFALSE;

  AliClusterContainer *cont = GetClusterContainer(c);
  if(!cont) {
    AliError(Form("%s:Container %d not found",GetName(),c));
    return 0;
  }

  return cont->AcceptCluster(clus);

}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalDev::AcceptTrack(AliVParticle *track, Int_t c) const
{
  // Return true if track is accepted.

  if (!track)
    return kFALSE;

  AliParticleContainer *cont = GetParticleContainer(c);
  if(!cont) {
    AliError(Form("%s:Container %d not found",GetName(),c));
    return 0;
  }

  return cont->AcceptParticle(track);

  if (!track)
    return kFALSE;

}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalDev::AcceptEmcalPart(AliEmcalParticle *part) const
{
  // Return true if EMCal particle is accepted.

  if (!part)
    return kFALSE;

  if (part->IsTrack()) {
    if (part->IsMC()) { 
      if (part->TestBits(fMCTrackBitMap) != (Int_t)fMCTrackBitMap)
	return kFALSE;
    }
    else {
      if (part->TestBits(fTrackBitMap) != (Int_t)fTrackBitMap)
	return kFALSE;
    }

    if (part->Pt() < fTrackPtCut)
      return kFALSE;

    if (part->Eta() < fTrackMinEta || part->Eta() > fTrackMaxEta || 
	part->Phi() < fTrackMinPhi || part->Phi() > fTrackMaxPhi)
      return kFALSE;
  }

  if (part->IsCluster()) {
    if (part->IsMC(fMinMCLabel)) { 
      if (part->TestBits(fMCClusterBitMap) != (Int_t)fMCClusterBitMap)
	return kFALSE;
    }
    else {
      if (part->TestBits(fClusterBitMap) != (Int_t)fClusterBitMap)
	return kFALSE;
    }

    if (!part->IsEMCAL())
      return kFALSE;

    if (part->Pt() < fClusPtCut)
      return kFALSE;
  }

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalDev::ExecOnce()
{
  // Init the analysis.

  if (!InputEvent()) {
    AliError(Form("%s: Could not retrieve event! Returning!", GetName()));
    return;
  }

  fGeom = AliEMCALGeometry::GetInstance();
  if (!fGeom) {
    AliError(Form("%s: Can not create geometry", GetName()));
    return;
  }

  if (fEventPlaneVsEmcal >= 0) {
    Double_t ep = (fGeom->GetArm1PhiMax() + fGeom->GetArm1PhiMin()) / 2 * TMath::DegToRad() + fEventPlaneVsEmcal - TMath::Pi();
    fMinEventPlane = ep - TMath::Pi() / 4;
    fMaxEventPlane = ep + TMath::Pi() / 4;
  }

  //Load all requested track branches - each container knows name already
  for(Int_t i =0; i<fParticleCollArray.GetEntriesFast(); i++) {
    AliParticleContainer *cont = static_cast<AliParticleContainer*>(fParticleCollArray.At(i));
    cont->SetParticleArray(InputEvent());
  }

  fTracks = GetParticleArray(0);
  if(!fTracks && fParticleCollArray.GetEntriesFast()>0) {
    AliError(Form("%s: Could not retrieve first track branch!", GetName()));
    return;
  }

  //Load all requested cluster branches - each container knows name already
  for(Int_t i =0; i<fClusterCollArray.GetEntriesFast(); i++) {
    AliClusterContainer *cont = static_cast<AliClusterContainer*>(fClusterCollArray.At(i));
    cont->SetClusterArray(InputEvent());
  }
  fCaloClusters = GetClusterArray(0);
  if(!fCaloClusters && fClusterCollArray.GetEntriesFast()>0) {
    AliError(Form("%s: Could not retrieve first cluster branch!", GetName()));
    return;
  }

  if (!fCaloCellsName.IsNull() && !fCaloCells) {
    fCaloCells =  dynamic_cast<AliVCaloCells*>(InputEvent()->FindListObject(fCaloCellsName));
    if (!fCaloCells) {
      AliError(Form("%s: Could not retrieve cells %s!", GetName(), fCaloCellsName.Data())); 
      return;
    }
  }


  fInitialized = kTRUE;
}

//_____________________________________________________
AliAnalysisTaskEmcalDev::BeamType AliAnalysisTaskEmcalDev::GetBeamType()
{
  // Get beam type : pp-AA-pA
  // ESDs have it directly, AODs get it from hardcoded run number ranges

  if (fForceBeamType != kNA)
    return fForceBeamType;

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
    } 
    else if ((runNumber>=188365 && runNumber <= 188366) || // LHC12g
	     (runNumber >= 195344 && runNumber <= 196608))  // LHC13b-f
    {
      return kpA;
    } else {
      return kpp;
    }
  }  
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalDev::IsEventSelected()
{
  // Check if event is selected

  if (fOffTrigger != AliVEvent::kAny) {
    UInt_t res = 0;
    const AliESDEvent *eev = dynamic_cast<const AliESDEvent*>(InputEvent());
    if (eev) {
      res = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    } else {
      const AliAODEvent *aev = dynamic_cast<const AliAODEvent*>(InputEvent());
      if (aev) {
        res = aev->GetHeader()->GetOfflineTrigger();
      }
    }
    if ((res & fOffTrigger) == 0)
      return kFALSE;
  }

  if (!fTrigClass.IsNull()) {
    TString fired;
    const AliESDEvent *eev = dynamic_cast<const AliESDEvent*>(InputEvent());
    if (eev) {
      fired = eev->GetFiredTriggerClasses();
    } else {
      const AliAODEvent *aev = dynamic_cast<const AliAODEvent*>(InputEvent());
      if (aev) {
        fired = aev->GetFiredTriggerClasses();
      }
    }
    if (!fired.Contains("-B-"))
      return kFALSE;
    TObjArray *arr = fTrigClass.Tokenize("|");
    if (!arr)
      return kFALSE;
    Bool_t match = 0;
    for (Int_t i=0;i<arr->GetEntriesFast();++i) {
      TObject *obj = arr->At(i);
      if (!obj)
        continue;
      if (fired.Contains(obj->GetName())) {
        match = 1;
        break;
      }
    }
    delete arr;
    if (!match)
      return kFALSE;
  }

  if ((fMinCent != -999) && (fMaxCent != -999)) {
    if (fCent<fMinCent)
      return kFALSE;
    if (fCent>fMaxCent)
      return kFALSE;
  }

  if ((fMinVz != -999) && (fMaxVz != -999)) {
    if (fNVertCont == 0 )
      return kFALSE;
    Double_t vz = fVertex[2];
    if (vz<fMinVz)
      return kFALSE;
    if (vz>fMaxVz)
      return kFALSE;
  }

  if (fMinPtTrackInEmcal > 0 && fTracks && fGeom) { //what is fTracks for? remove?
    Bool_t trackInEmcalOk = kFALSE;
    Int_t ntracks = GetNParticles(0);
    for (Int_t i = 0; i < ntracks; i++) {
      AliVParticle *track = GetAcceptParticleFromArray(i,0);
      if(!track)
	continue;
      if (track->Eta() < fGeom->GetArm1EtaMin() || track->Eta() > fGeom->GetArm1EtaMax() ||
	  track->Phi() < fGeom->GetArm1PhiMin() * TMath::DegToRad() || track->Phi() > fGeom->GetArm1PhiMax() * TMath::DegToRad())
	continue;
      if (track->Pt() > fMinPtTrackInEmcal) {
	trackInEmcalOk = kTRUE;
	break;
      }
    }
    if (!trackInEmcalOk)
      return kFALSE;
  }


  if (!(fEPV0 > fMinEventPlane && fEPV0 <= fMaxEventPlane) &&
      !(fEPV0 + TMath::Pi() > fMinEventPlane && fEPV0 + TMath::Pi() <= fMaxEventPlane) &&
      !(fEPV0 - TMath::Pi() > fMinEventPlane && fEPV0 - TMath::Pi() <= fMaxEventPlane)) 
    return kFALSE;


  return kTRUE;
}

//________________________________________________________________________
TClonesArray *AliAnalysisTaskEmcalDev::GetArrayFromEvent(const char *name, const char *clname)
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
  } else {
    return 0;
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
Bool_t AliAnalysisTaskEmcalDev::RetrieveEventObjects()
{
  // Retrieve objects from event.

  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;
  fNVertCont = 0;

  const AliVVertex *vert = InputEvent()->GetPrimaryVertex();
  if (vert) {
    vert->GetXYZ(fVertex);
    fNVertCont = vert->GetNContributors();
  }

  fBeamType = GetBeamType();

  if (fBeamType == kAA || fBeamType == kpA ) {
    AliCentrality *aliCent = InputEvent()->GetCentrality();
    if (aliCent) {
      fCent = aliCent->GetCentralityPercentile(fCentEst.Data()); 
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
    AliEventplane *aliEP = InputEvent()->GetEventplane();
    if (aliEP) {
      fEPV0  = aliEP->GetEventplane("V0" ,InputEvent());
      fEPV0A = aliEP->GetEventplane("V0A",InputEvent());
      fEPV0C = aliEP->GetEventplane("V0C",InputEvent());
    } else {
      AliWarning(Form("%s: Could not retrieve event plane information!", GetName()));
    }
  } else {
    fCent = 99;
    fCentBin = 0;
  }

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalDev::AddParticleContainer(const char *n) {

  // Add particle container
  // will be called in AddTask macro

  AliParticleContainer *cont = 0x0;
  cont = new AliParticleContainer();
  cont->SetArrayName(n);
  TString contName = cont->GetArrayName();
 
  fParticleCollArray.Add(cont);

}

//________________________________________________________________________
void AliAnalysisTaskEmcalDev::AddClusterContainer(const char *n) {

  // Add cluster container
  // will be called in AddTask macro

  AliClusterContainer *cont = 0x0;
  cont = new AliClusterContainer();
  cont->SetArrayName(n);

  fClusterCollArray.Add(cont);

}

//________________________________________________________________________
AliParticleContainer* AliAnalysisTaskEmcalDev::GetParticleContainer(Int_t i) const {
  // Get i^th particle container

  if(i<0 || i>fParticleCollArray.GetEntriesFast()) return 0;
  AliParticleContainer *cont = static_cast<AliParticleContainer*>(fParticleCollArray.At(i));
  return cont;
}

//________________________________________________________________________
AliClusterContainer* AliAnalysisTaskEmcalDev::GetClusterContainer(Int_t i) const {
  // Get i^th cluster container

  if(i<0 || i>fClusterCollArray.GetEntriesFast()) return 0;
  AliClusterContainer *cont = static_cast<AliClusterContainer*>(fClusterCollArray.At(i));
  return cont;
}

//________________________________________________________________________
TClonesArray* AliAnalysisTaskEmcalDev::GetParticleArray(Int_t i) const {
  // Get i^th TClonesArray with AliVParticle

  AliParticleContainer *cont = GetParticleContainer(i);
  if(!cont) {
    AliError(Form("%s: Particle container %d not found",GetName(),i));
    return 0;
  }
  TString contName = cont->GetArrayName();
  return cont->GetArray();

}

//________________________________________________________________________
TClonesArray* AliAnalysisTaskEmcalDev::GetClusterArray(Int_t i) const {
  // Get i^th TClonesArray with AliVCluster

  AliClusterContainer *cont = GetClusterContainer(i);
  if(!cont) {
    AliError(Form("%s:Cluster container %d not found",GetName(),i));
    return 0;
  }
  return cont->GetArray();

}

//________________________________________________________________________
AliVParticle* AliAnalysisTaskEmcalDev::GetAcceptParticleFromArray(Int_t p, Int_t c) const {
  // Get particle p if accepted from  container c
  // If particle not accepted return 0

  AliParticleContainer *cont = GetParticleContainer(c);
  if(!cont) {
    AliError(Form("%s: Particle container %d not found",GetName(),c));
    return 0;
  }
  AliVParticle *vp = cont->GetAcceptParticle(p);

  return vp;
  
}

//________________________________________________________________________
AliVCluster* AliAnalysisTaskEmcalDev::GetAcceptClusterFromArray(Int_t cl, Int_t c) const {
  // Get particle p if accepted from  container c
  // If particle not accepted return 0

  AliClusterContainer *cont = GetClusterContainer(c);
  if(!cont) {
    AliError(Form("%s: Cluster container %d not found",GetName(),c));
    return 0;
  }
  AliVCluster *vc = cont->GetAcceptCluster(cl);

  return vc;
  
}

//________________________________________________________________________
Int_t AliAnalysisTaskEmcalDev::GetNParticles(Int_t i) const {
  // Get number of entries in particle array i

  AliParticleContainer *cont = GetParticleContainer(i);
  if(!cont) {
    AliError(Form("%s: Particle container %d not found",GetName(),i));
    return 0;
  }
  return cont->GetNEntries();

}

//________________________________________________________________________
Int_t AliAnalysisTaskEmcalDev::GetNClusters(Int_t i) const {
  // Get number of entries in cluster array i

  AliClusterContainer *cont = GetClusterContainer(i);
  if(!cont) {
    AliError(Form("%s: Cluster container %d not found",GetName(),i));
    return 0;
  }
  return cont->GetNEntries();

}

