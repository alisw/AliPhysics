// $Id$
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
  fRhoName(),
  fPtBiasJetTrack(0),
  fPtBiasJetClus(0),
  fJetPtCut(1),
  fJetAreaCut(-1),
  fPercAreaCut(-1),
  fAreaEmcCut(0),
  fJetMinEta(-0.9),
  fJetMaxEta(0.9),
  fJetMinPhi(-10),
  fJetMaxPhi(10),
  fMaxClusterPt(100),
  fMaxTrackPt(100),
  fLeadingHadronType(0),
  fNLeadingJets(2),
  fJets(0),
  fRho(0),
  fRhoVal(0)
{
  // Default constructor.
}

//________________________________________________________________________
AliAnalysisTaskEmcalJet::AliAnalysisTaskEmcalJet(const char *name, Bool_t histo) : 
  AliAnalysisTaskEmcal(name, histo),
  fJetRadius(0.4),
  fJetsName(),
  fRhoName(),
  fPtBiasJetTrack(0),
  fPtBiasJetClus(0),
  fJetPtCut(1),
  fJetAreaCut(-1),
  fPercAreaCut(-1),
  fAreaEmcCut(0),
  fJetMinEta(-0.9),
  fJetMaxEta(0.9),
  fJetMinPhi(-10),
  fJetMaxPhi(10),
  fMaxClusterPt(100),
  fMaxTrackPt(100),
  fLeadingHadronType(0),
  fNLeadingJets(2),
  fJets(0),
  fRho(0),
  fRhoVal(0)
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

  if (fLeadingHadronType == 0) {
    if (jet->MaxTrackPt() < fPtBiasJetTrack) return kFALSE;
  }
  else if (fLeadingHadronType == 1) {
    if (jet->MaxClusterPt() < fPtBiasJetClus) return kFALSE;
  }
  else {
    if (jet->MaxTrackPt() < fPtBiasJetTrack && jet->MaxClusterPt() < fPtBiasJetClus) return kFALSE;
  }

  return kTRUE;
}

//________________________________________________________________________
Float_t* AliAnalysisTaskEmcalJet::GenerateFixedBinArray(Int_t n, Float_t min, Float_t max) const
{
  Float_t *bins = new Float_t[n+1];

  Float_t binWidth = (max-min)/n;
  bins[0] = min;
  for (Int_t i = 1; i <= n; i++) {
    bins[i] = bins[i-1]+binWidth;
  }

  return bins;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalJet::GetLeadingHadronPt(AliEmcalJet *jet) const
{
  if (fLeadingHadronType == 0)       // charged leading hadron
    return jet->MaxTrackPt();
  else if (fLeadingHadronType == 1)  // neutral leading hadron
    return jet->MaxClusterPt();
  else                               // charged or neutral
    return jet->MaxPartPt();
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJet::AcceptJet(AliEmcalJet *jet) const
{   
  // Return true if jet is accepted.

  if (jet->Pt() <= fJetPtCut) 
    return kFALSE;
 
  if (jet->Area() <= fJetAreaCut) 
    return kFALSE;

  if (jet->AreaEmc()<fAreaEmcCut)
    return kFALSE;

  if (!AcceptBiasJet(jet))
    return kFALSE;
  
  if (jet->MaxTrackPt() > fMaxTrackPt || jet->MaxClusterPt() > fMaxClusterPt)
    return kFALSE;

  Double_t jetPhi = jet->Phi();
  Double_t jetEta = jet->Eta();
  
  if (fJetMinPhi < 0) // if limits are given in (-pi, pi) range
    jetPhi -= TMath::Pi() * 2;

  return (Bool_t)(jetEta > fJetMinEta && jetEta < fJetMaxEta && jetPhi > fJetMinPhi && jetPhi < fJetMaxPhi);
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

  AliAnalysisTaskEmcal::ExecOnce();

  if (fPercAreaCut >= 0) {
    if (fJetAreaCut >= 0)
      AliInfo(Form("%s: jet area cut will be calculated as a percentage of the average area, given value will be overwritten", GetName()));
    fJetAreaCut = fPercAreaCut * fJetRadius * fJetRadius * TMath::Pi();
  }
  if (fJetAreaCut < 0)
    fJetAreaCut = 0;

  if (fAnaType == kTPC) {
    SetJetEtaLimits(-0.5, 0.5);
    SetJetPhiLimits(-10, 10);
  } 
  else if (fAnaType == kEMCAL && fGeom) {
    SetJetEtaLimits(fGeom->GetArm1EtaMin() + fJetRadius, fGeom->GetArm1EtaMax() - fJetRadius);
    SetJetPhiLimits(fGeom->GetArm1PhiMin() * TMath::DegToRad() + fJetRadius, fGeom->GetArm1PhiMax() * TMath::DegToRad() - fJetRadius);
  }

  if (!fRhoName.IsNull() && !fRho) {
    fRho = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(fRhoName));
    if (!fRho) {
      AliError(Form("%s: Could not retrieve rho %s!", GetName(), fRhoName.Data()));
      fInitialized = kFALSE;
      return;
    }
  }

  if (!fJetsName.IsNull() && !fJets) {
    fJets = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fJetsName));
    if (!fJets) {
      AliError(Form("%s: Could not retrieve jets %s!", GetName(), fJetsName.Data()));
      fInitialized = kFALSE;
      return;
    }
    else if (!fJets->GetClass()->GetBaseClass("AliEmcalJet")) {
      AliError(Form("%s: Collection %s does not contain AliEmcalJet objects!", GetName(), fJetsName.Data())); 
      fJets = 0;
      fInitialized = kFALSE;
      return;
    }
  }
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJet::GetSortedArray(Int_t indexes[], TClonesArray *array, Double_t rho) const
{
  // Get the leading jets.

  static Float_t pt[9999] = {0};

  if (!array)
    return 0;

  const Int_t n = array->GetEntriesFast();

  if (n < 1)
    return kFALSE;

  if (array->GetClass()->GetBaseClass("AliEmcalJet")) {

    for (Int_t i = 0; i < n; i++) {

      pt[i] = -FLT_MAX;

      AliEmcalJet* jet = static_cast<AliEmcalJet*>(array->At(i));
      
      if (!jet) {
	AliError(Form("Could not receive jet %d", i));
	continue;
      }
      
      if (!AcceptJet(jet))
	continue;
      
      pt[i] = jet->Pt() - rho * jet->Area();
    }
  }

  else if (array->GetClass()->GetBaseClass("AliVTrack")) {

    for (Int_t i = 0; i < n; i++) {

      pt[i] = -FLT_MAX;

      AliVTrack* track = static_cast<AliVTrack*>(array->At(i));
      
      if (!track) {
	AliError(Form("Could not receive track %d", i));
	continue;
      }  
      
      if (!AcceptTrack(track))
	continue;
      
      pt[i] = track->Pt();
    }
  }

  else if (array->GetClass()->GetBaseClass("AliVCluster")) {

    for (Int_t i = 0; i < n; i++) {

      pt[i] = -FLT_MAX;

      AliVCluster* cluster = static_cast<AliVCluster*>(array->At(i));
      
      if (!cluster) {
	AliError(Form("Could not receive cluster %d", i));
	continue;
      }  
      
      if (!AcceptCluster(cluster))
	continue;

      TLorentzVector nPart;
      cluster->GetMomentum(nPart, const_cast<Double_t*>(fVertex));
      
      pt[i] = nPart.Pt();
    }
  }

  TMath::Sort(n, pt, indexes);

  if (pt[indexes[0]] == -FLT_MAX) 
    return 0;

  return kTRUE;
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

  if (fRho)
    fRhoVal = fRho->GetVal();

  return kTRUE;
}
