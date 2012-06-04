// $Id: AliAnalysisTaskEmcalJet.cxx 56756 2012-05-30 05:03:02Z loizides $
//
// Emcal jet analysis base task.
//
// Author: S.Aiola

#include "AliAnalysisTaskEmcalJet.h"

#include <TObject.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TList.h>

#include "AliESDEvent.h"
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliVCluster.h"
#include "AliVParticle.h"
#include "AliEmcalJet.h"
#include "AliVEventHandler.h"
#include "AliLog.h"
#include "AliEMCALGeometry.h"

ClassImp(AliAnalysisTaskEmcalJet)

//________________________________________________________________________
AliAnalysisTaskEmcalJet::AliAnalysisTaskEmcalJet() : 
  AliAnalysisTaskEmcal("AliAnalysisTaskEmcalJet"),
  fJetRadius(0.4),
  fJetsName("Jets"),
  fPtBiasJetTrack(10),
  fPtBiasJetClus(10),
  fJetPtCut(1),
  fJetAreaCut(0.2),
  fMinEta(-0.9),
  fMaxEta(0.9),
  fMinPhi(-10),
  fMaxPhi(10),
  fJets(0)
{
  // Default constructor.
}

//________________________________________________________________________
AliAnalysisTaskEmcalJet::AliAnalysisTaskEmcalJet(const char *name) : 
  AliAnalysisTaskEmcal(name),
  fJetRadius(0.4),
  fJetsName("Jets"),
  fPtBiasJetTrack(10),
  fPtBiasJetClus(10),
  fJetPtCut(1),
  fJetAreaCut(0.2),
  fMinEta(-0.9),
  fMaxEta(0.9),
  fMinPhi(-10),
  fMaxPhi(10),
  fJets(0)
{
  // Standard constructor.
}

//________________________________________________________________________
AliAnalysisTaskEmcalJet::AliAnalysisTaskEmcalJet(const char *name, Bool_t histo) : 
  AliAnalysisTaskEmcal(name, histo),
  fJetRadius(0.4),
  fJetsName("Jets"),
  fPtBiasJetTrack(10),
  fPtBiasJetClus(10),
  fJetPtCut(1),
  fJetAreaCut(0.2),
  fMinEta(-0.9),
  fMaxEta(0.9),
  fMinPhi(-10),
  fMaxPhi(10),
  fJets(0)
{
  // Standard constructor.
}

//________________________________________________________________________
AliAnalysisTaskEmcalJet::~AliAnalysisTaskEmcalJet()
{
  // Destructor
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJet::RetrieveEventObjects()
{
  // Retrieve objects from event.

  if (!AliAnalysisTaskEmcal::RetrieveEventObjects())
    return kFALSE;

  if (!fJetsName.IsNull()) {
    fJets = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fJetsName));
    if (!fJets) {
      AliWarning(Form("Could not retrieve jets %s!", fJetsName.Data())); 
    }
  }

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJet::IsJetTrack(AliEmcalJet* jet, Int_t itrack, Bool_t sorted) const
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
Bool_t AliAnalysisTaskEmcalJet::IsJetCluster(AliEmcalJet* jet, Int_t iclus, Bool_t sorted) const
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
Bool_t AliAnalysisTaskEmcalJet::AcceptJet(AliEmcalJet *jet, Bool_t bias) const
{   
  // Return true if jet is accepted.

  if (jet->Pt() <= fJetPtCut)
    return kFALSE;
  if (jet->Area() <= fJetAreaCut)
    return kFALSE;
  if (bias && jet->MaxTrackPt() < fPtBiasJetTrack && (fAnaType == kTPC || jet->MaxClusterPt() < fPtBiasJetClus))
    return kFALSE;

  return (Bool_t)(jet->Eta() > fMinEta && jet->Eta() < fMaxEta && jet->Phi() > fMinPhi && jet->Phi() < fMaxPhi);
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJet::Init()
{
  // Init the analysis.

  if (fAnaType == kTPC) {
    SetEtaLimits(-0.9 + fJetRadius, 0.9 - fJetRadius);
    SetPhiLimits(-10, 10);
  }
  else if (fAnaType == kEMCAL) {
    AliEMCALGeometry *geom = AliEMCALGeometry::GetInstance();
    if (!geom) {
      AliFatal("Can not create geometry");
      return;
    }
    SetEtaLimits(geom->GetArm1EtaMin() + fJetRadius, geom->GetArm1EtaMax() - fJetRadius);
    SetPhiLimits(geom->GetArm1PhiMin() * TMath::DegToRad() + fJetRadius, geom->GetArm1PhiMax() * TMath::DegToRad() - fJetRadius);
  }
  else {
    AliWarning("Analysis type not recognized! Assuming kTPC...");
    SetAnaType(kTPC);
    Init();
  }

  SetInitialized();
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJet::UserExec(Option_t *) 
{
  if (!fInitialized) 
    Init();

  AliAnalysisTaskEmcal::UserExec("");
}
