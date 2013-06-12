//
// Jet Container
//
// Author: M. Verweij

#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>

#include <TChain.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TObject.h>
#include "AliEmcalJet.h"
#include "AliVEvent.h"
#include "AliLog.h"
#include "AliRhoParameter.h"
#include "AliEMCALGeometry.h"

#include "AliJetContainer.h"

ClassImp(AliJetContainer)

//________________________________________________________________________
AliJetContainer::AliJetContainer():
  AliEmcalContainer("AliJetContainer"),
  fJetAcceptanceType(kUser),
  fJetRadius(0),
  fRhoName(),
  fPtBiasJetTrack(0),
  fPtBiasJetClus(0),
  fJetPtCut(1),
  fJetAreaCut(-1),
  fAreaEmcCut(0),
  fJetMinEta(-0.9),
  fJetMaxEta(0.9),
  fJetMinPhi(-10),
  fJetMaxPhi(10),
  fMaxClusterPt(1000),
  fMaxTrackPt(100),
  fLeadingHadronType(0),
  fNLeadingJets(1),
  fJetBitMap(0),
  fRho(0),
  fGeom(0)
{
  // Default constructor.

}

//________________________________________________________________________
AliJetContainer::AliJetContainer(const char *name):
  AliEmcalContainer(name),
  fJetAcceptanceType(kUser),
  fJetRadius(0),
  fRhoName(),
  fPtBiasJetTrack(0),
  fPtBiasJetClus(0),
  fJetPtCut(1),
  fJetAreaCut(-1),
  fAreaEmcCut(0),
  fJetMinEta(-0.9),
  fJetMaxEta(0.9),
  fJetMinPhi(-10),
  fJetMaxPhi(10),
  fMaxClusterPt(1000),
  fMaxTrackPt(100),
  fLeadingHadronType(0),
  fNLeadingJets(1),
  fJetBitMap(0),
  fRho(0),
  fGeom(0)
{
  // Standard constructor.

}

//________________________________________________________________________
AliJetContainer::~AliJetContainer()
{
  // Destructor.
}

//________________________________________________________________________
void AliJetContainer::SetJetArray(AliVEvent *event) 
{
  // Set jet array

  //  SetArray(event, fClArrayName.GetName(), "AliEmcalJet");
  SetArray(event, "AliEmcalJet");

  if(fJetAcceptanceType==kTPC)
    SetJetEtaPhiTPC();
  else if(fJetAcceptanceType==kEMCAL)
    SetJetEtaPhiEMCAL();
}

//________________________________________________________________________
void AliJetContainer::SetEMCALGeometry() {
  fGeom = AliEMCALGeometry::GetInstance();
  if (!fGeom) {
    AliError(Form("%s: Can not create geometry", GetName()));
    return;
  }
}

//________________________________________________________________________
void AliJetContainer::LoadRho(AliVEvent *event)
{
  //Load rho

  if (!fRhoName.IsNull() && !fRho) {
    fRho = dynamic_cast<AliRhoParameter*>(event->FindListObject(fRhoName));
    if (!fRho) {
      AliError(Form("%s: Could not retrieve rho %s!", GetName(), fRhoName.Data()));
      return;
    }
  }
}

//________________________________________________________________________
AliEmcalJet* AliJetContainer::GetJet(Int_t i) const {

  //Get i^th jet in array

  if(i<0 || i>fClArray->GetEntriesFast()) return 0;
  AliEmcalJet *jet = static_cast<AliEmcalJet*>(fClArray->At(i));
  return jet;

}


//________________________________________________________________________
Double_t AliJetContainer::GetJetPtCorr(Int_t i) const {
  AliEmcalJet *jet = GetJet(i);

  return jet->Pt() - fRho->GetVal()*jet->Area();
}

//________________________________________________________________________
AliEmcalJet* AliJetContainer::GetAcceptJet(Int_t i) const {

  //Only return jet if is accepted

  AliEmcalJet *jet = GetJet(i);
  if(!AcceptJet(jet)) return 0;

  return jet;
}

//________________________________________________________________________
Bool_t AliJetContainer::AcceptBiasJet(AliEmcalJet *jet) const
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
Bool_t AliJetContainer::AcceptJet(AliEmcalJet *jet) const
{   
  // Return true if jet is accepted.
  if (!jet)
    return kFALSE;

  if (jet->TestBits(fJetBitMap) != (Int_t)fJetBitMap)
    return kFALSE;

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
Double_t AliJetContainer::GetLeadingHadronPt(AliEmcalJet *jet) const
{
  if (fLeadingHadronType == 0)       // charged leading hadron
    return jet->MaxTrackPt();
  else if (fLeadingHadronType == 1)  // neutral leading hadron
    return jet->MaxClusterPt();
  else                               // charged or neutral
    return jet->MaxPartPt();
}


//________________________________________________________________________
void AliJetContainer::SetJetEtaPhiEMCAL()
{
  //Set default cuts for full jets

  if(!fGeom) SetEMCALGeometry();
  if(fGeom) {
    SetJetEtaLimits(fGeom->GetArm1EtaMin() + fJetRadius, fGeom->GetArm1EtaMax() - fJetRadius);
    SetJetPhiLimits(fGeom->GetArm1PhiMin() * TMath::DegToRad() + fJetRadius, fGeom->GetArm1PhiMax() * TMath::DegToRad() - fJetRadius);
  }
  else {
    AliWarning("Could not get instance of AliEMCALGeometry. Using manual settings");
    SetJetEtaLimits(-0.7+fJetRadius,0.7-fJetRadius);
    SetJetPhiLimits(1.4+fJetRadius,TMath::Pi()-fJetRadius);
  }

}

//________________________________________________________________________
void AliJetContainer::SetJetEtaPhiTPC()
{
  //Set default cuts for full jets

  SetJetEtaLimits(-0.5, 0.5);
  SetJetPhiLimits(-10, 10);

}

//________________________________________________________________________
void AliJetContainer::ResetCuts() 
{
  // Reset cuts to default values

  fPtBiasJetTrack = 0;
  fPtBiasJetClus = 0;
  fJetPtCut = 1;
  fJetAreaCut = -1;
  fAreaEmcCut = 0;
  fJetMinEta = -0.9;
  fJetMaxEta = 0.9;
  fJetMinPhi = -10;
  fJetMaxPhi = 10;
  fMaxClusterPt = 1000;
  fMaxTrackPt = 100;
  fLeadingHadronType = 0;

}
