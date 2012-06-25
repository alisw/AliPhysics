// $Id: AliAnalysisTaskEmcalJet.cxx 56756 2012-05-30 05:03:02Z loizides $
//
// Emcal jet analysis base task.
//
// Author: S.Aiola

#include "AliAnalysisTaskEmcalJet.h"

#include <TChain.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TObject.h>

#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliEMCALGeometry.h"
#include "AliESDEvent.h"
#include "AliEmcalJet.h"
#include "AliLog.h"
#include "AliRhoParameter.h"
#include "AliVCluster.h"
#include "AliVEventHandler.h"
#include "AliVParticle.h"

ClassImp(AliAnalysisTaskEmcalJet)

//________________________________________________________________________
AliAnalysisTaskEmcalJet::AliAnalysisTaskEmcalJet() : 
  AliAnalysisTaskEmcal("AliAnalysisTaskEmcalJet"),
  fJetRadius(0.4),
  fJetsName(),
  fPtBiasJetTrack(5),
  fPtBiasJetClus(5),
  fJetPtCut(1),
  fJetAreaCut(0.4),
  fPercAreaCut(0.8),
  fMinEta(-0.9),
  fMaxEta(0.9),
  fMinPhi(-10),
  fMaxPhi(10),
  fMaxClusterPt(100),
  fMaxTrackPt(100),
  fJets(0)
{
  // Default constructor.
}

//________________________________________________________________________
AliAnalysisTaskEmcalJet::AliAnalysisTaskEmcalJet(const char *name, Bool_t histo) : 
  AliAnalysisTaskEmcal(name, histo),
  fJetRadius(0.4),
  fJetsName(),
  fPtBiasJetTrack(5),
  fPtBiasJetClus(5),
  fJetPtCut(1),
  fJetAreaCut(0.4),
  fPercAreaCut(0.8),
  fMinEta(-0.9),
  fMaxEta(0.9),
  fMinPhi(-10),
  fMaxPhi(10),
  fMaxClusterPt(100),
  fMaxTrackPt(100),
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
Bool_t AliAnalysisTaskEmcalJet::AcceptBiasJet(AliEmcalJet *jet) const
{ 
  // Accept jet with a bias.

  if (jet->MaxTrackPt() < fPtBiasJetTrack && (fAnaType == kTPC || jet->MaxClusterPt() < fPtBiasJetClus))
    return kFALSE;
  else
    return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJet::AcceptJet(AliEmcalJet *jet, Bool_t bias, Bool_t upCut) const
{   
  // Return true if jet is accepted.

  if (jet->Pt() <= fJetPtCut)
    return kFALSE;
  if (jet->Area() <= fJetAreaCut)
    return kFALSE;
  if (bias && !AcceptBiasJet(jet))
    return kFALSE;
  if (upCut && (jet->MaxTrackPt() > fMaxTrackPt || jet->MaxClusterPt() > fMaxClusterPt))
    return kFALSE;

  return (Bool_t)(jet->Eta() > fMinEta && jet->Eta() < fMaxEta && jet->Phi() > fMinPhi && jet->Phi() < fMaxPhi);
}

//________________________________________________________________________
AliRhoParameter *AliAnalysisTaskEmcalJet::GetRhoFromEvent(const char *name)
{
  // Get rho from event.

  AliRhoParameter *rho = 0;
  TString sname(name);
  if (!sname.IsNull()) {
    rho = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(sname));
    if (!rho) {
      AliWarning(Form("%s: Could not retrieve rho with name %s!", GetName(), name)); 
      return 0;
    }
  }
  return rho;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJet::ExecOnce()
{
  // Init the analysis.

  if (fPercAreaCut >= 0) {
    AliInfo(Form("%s: jet area cut will be calculated as a percentage of the average area, given value will be overwritten", GetName()));
    fJetAreaCut = fPercAreaCut * fJetRadius * fJetRadius * TMath::Pi();
  }

  if (fAnaType == kTPC) {
    SetEtaLimits(-0.9 + fJetRadius, 0.9 - fJetRadius);
    SetPhiLimits(-10, 10);
  } else if (fAnaType == kEMCAL || fAnaType == kTPCSmall || fAnaType == kEMCALOnly) {
    AliEMCALGeometry *geom = AliEMCALGeometry::GetInstance();
    if (geom) {
      SetEtaLimits(geom->GetArm1EtaMin() + fJetRadius, geom->GetArm1EtaMax() - fJetRadius);
      SetPhiLimits(geom->GetArm1PhiMin() * TMath::DegToRad() + fJetRadius, geom->GetArm1PhiMax() * TMath::DegToRad() - fJetRadius);
    }
    else {
      AliWarning(Form("%s: Can not create geometry", GetName()));
    }
  } else {
    AliWarning(Form("%s: Analysis type not recognized! Assuming kTPC!", GetName()));
    SetAnaType(kTPC);
    ExecOnce();
    return;
  }

  if (fAnaType == kTPCSmall)
    fAnaType = kTPC;

  AliAnalysisTaskEmcal::ExecOnce();
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJet::IsJetCluster(AliEmcalJet* jet, Int_t iclus, Bool_t sorted) const
{
  // Return true if cluster is in jet.

  for (Int_t i = 0; i < jet->GetNumberOfClusters(); ++i) {
    Int_t ijetclus = jet->ClusterAt(i);
    if (sorted && ijetclus > iclus)
      return kFALSE;
    if (ijetclus == iclus)
      return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJet::IsJetTrack(AliEmcalJet* jet, Int_t itrack, Bool_t sorted) const
{
  // Return true if track is in jet.

  for (Int_t i = 0; i < jet->GetNumberOfTracks(); ++i) {
    Int_t ijettrack = jet->TrackAt(i);
    if (sorted && ijettrack > itrack)
      return kFALSE;
    if (ijettrack == itrack)
      return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJet::RetrieveEventObjects()
{
  // Retrieve objects from event.

  if (!AliAnalysisTaskEmcal::RetrieveEventObjects())
    return kFALSE;

  if (!fJetsName.IsNull() && !fJets) {
    fJets = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fJetsName));
    if (!fJets) {
      AliError(Form("%s: Could not retrieve jets %s!", GetName(), fJetsName.Data()));
      return kFALSE;
    }
    else if (!fJets->GetClass()->GetBaseClass("AliEmcalJet")) {
      AliError(Form("%s: Collection %s does not contain AliEmcalJet objects!", GetName(), fJetsName.Data())); 
      fJets = 0;
      return kFALSE;
    }
  }

  return kTRUE;
}
