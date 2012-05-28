// $Id: AliAnalysisTaskEmcal.cxx $
//
// Emcal base analysis task.
//
// Author: S.Aiola

#include "AliAnalysisTaskEmcal.h"

#include <TObject.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TList.h>

#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliVCluster.h"
#include "AliVParticle.h"
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
  fMinEta(-0.9),
  fMaxEta(0.9),
  fMinPhi(-10),
  fMaxPhi(10),
  fJetRadius(0.4),
  fTracksName("Tracks"),
  fCaloName("CaloClusters"),
  fJetsName("Jets"),
  fNbins(500),
  fMinPt(0),
  fMaxPt(250),
  fPtCut(0.15),
  fPtBiasJetTrack(10),
  fPtBiasJetClus(10),
  fJetPtCut(1),
  fJetAreaCut(0.2),
  fTracks(0),
  fCaloClusters(0),
  fJets(0),
  fCent(0),
  fCentBin(-1),
  fOutput(0)
{
  // Default constructor.

  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;
}

//________________________________________________________________________
AliAnalysisTaskEmcal::AliAnalysisTaskEmcal(const char *name) : 
  AliAnalysisTaskSE(name),
  fAnaType(kTPC),
  fInitialized(kFALSE),
  fMinEta(-0.9),
  fMaxEta(0.9),
  fMinPhi(-10),
  fMaxPhi(10),
  fJetRadius(0.4),
  fTracksName("Tracks"),
  fCaloName("CaloClusters"),
  fJetsName("Jets"),
  fNbins(500),
  fMinPt(0),
  fMaxPt(250),
  fPtCut(0.15),
  fPtBiasJetTrack(10),
  fPtBiasJetClus(10),
  fJetPtCut(1),
  fJetAreaCut(0.2),
  fTracks(0),
  fCaloClusters(0),
  fJets(0),
  fCent(0),
  fCentBin(-1),
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
AliAnalysisTaskEmcal::~AliAnalysisTaskEmcal()
{
  // Destructor
}

//________________________________________________________________________
void AliAnalysisTaskEmcal::UserCreateOutputObjects()
{
  // User create outputs.
}

//________________________________________________________________________
void AliAnalysisTaskEmcal::RetrieveEventObjects()
{
  // Retrieve objects from event.

  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;
  InputEvent()->GetPrimaryVertex()->GetXYZ(fVertex);

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

  if (!fJetsName.IsNull()) {
    fJets = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fJetsName));
    if (!fJets) {
      AliWarning(Form("Could not retrieve jets %s!", fJetsName.Data())); 
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskEmcal::Init()
{
  // Init the analysis.

  if (fAnaType == kTPC) {
    SetEtaLimits(-0.9, 0.9);
    SetPhiLimits(-10, 10);
  }
  else if (fAnaType == kEMCAL) {
    AliEMCALGeometry *geom = AliEMCALGeometry::GetInstance();
    if (!geom) {
      AliFatal("Can not create geometry");
      return;
    }
    SetEtaLimits(geom->GetArm1EtaMin(), geom->GetArm1EtaMax());
    SetPhiLimits(geom->GetArm1PhiMin() * TMath::DegToRad(), geom->GetArm1PhiMax() * TMath::DegToRad());
  }
  else {
    AliWarning("Analysis type not recognized! Assuming kTPC...");
    SetAnaType(kTPC);
    Init();
  }

  SetInitialized();
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcal::IsJetTrack(AliEmcalJet* jet, Int_t itrack, Bool_t sorted) const
{
  // Return true if track is in jet.

  for (Int_t i = 0; i < jet->GetNumberOfTracks(); i++) {
    Int_t ijettrack = jet->TrackAt(i);
    if (sorted && ijettrack > itrack)
      return kFALSE;
    if (ijettrack == itrack)
      return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcal::IsJetCluster(AliEmcalJet* jet, Int_t iclus, Bool_t sorted) const
{
  // Return true if cluster is in jet.

  for (Int_t i = 0; i < jet->GetNumberOfClusters(); i++) {
    Int_t ijetclus = jet->ClusterAt(i);
    if (sorted && ijetclus > iclus)
      return kFALSE;
    if (ijetclus == iclus)
      return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcal::AcceptJet(AliEmcalJet *jet) const
{   
  // Return true if jet is accepted.

  if (jet->Pt() <= fJetPtCut)
    return kFALSE;
  if (jet->Area() <= fJetAreaCut)
    return kFALSE;
  if (fAnaType == kEMCAL && !jet->IsInsideEmcal())
    return kFALSE;

  return (Bool_t)(jet->Eta() > fMinEta && jet->Eta() < fMaxEta && jet->Phi() > fMinPhi && jet->Phi() < fMaxPhi);
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcal::AcceptCluster(AliVCluster* clus, Bool_t acceptMC) const
{
  // Return true if cluster is accepted.

  if (!acceptMC && clus->Chi2() == 100)
    return kFALSE;

  TLorentzVector nPart;
  clus->GetMomentum(nPart, const_cast<Double_t*>(fVertex));

  if (nPart.Et() < fPtCut)
    return kFALSE;

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcal::AcceptTrack(AliVParticle* track, Bool_t acceptMC) const
{
  // Return true if track is accepted.

  if (!acceptMC && track->GetLabel() == 100)
    return kFALSE;

  if (track->Pt() < fPtCut)
    return kFALSE;
  
  return (Bool_t)(track->Eta() > fMinEta && track->Eta() < fMaxEta && track->Phi() > fMinPhi && track->Phi() < fMaxPhi);
}

//________________________________________________________________________
void AliAnalysisTaskEmcal::UserExec(Option_t *) 
{
  // Main loop, called for each event.

  if (!fInitialized) 
    Init();

  RetrieveEventObjects();

  FillHistograms();
    
  // information for this iteration of the UserExec in the container
  PostData(1, fOutput);
}

//________________________________________________________________________
void AliAnalysisTaskEmcal::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}
