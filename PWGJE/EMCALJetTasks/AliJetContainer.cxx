// $Id$
//
// Container with name, TClonesArray and cuts for jets
//
// Author: M. Verweij

#include <TClonesArray.h>

#include "AliEmcalJet.h"
#include "AliVEvent.h"
#include "AliLog.h"
#include "AliEMCALGeometry.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliLocalRhoParameter.h"

#include "AliJetContainer.h"

ClassImp(AliJetContainer)

//________________________________________________________________________
AliJetContainer::AliJetContainer():
  AliEmcalContainer("AliJetContainer"),
  fJetAcceptanceType(kUser),
  fJetRadius(0),
  fRhoName(),
  fLocalRhoName(),
  fFlavourSelection(0),
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
  fZLeadingEmcCut(10.),
  fZLeadingChCut(10.),
  fNEFMinCut(-10.),
  fNEFMaxCut(10.),
  fLeadingHadronType(0),
  fNLeadingJets(1),
  fJetBitMap(0),
  fJetTrigger(0),
  fParticleContainer(0),
  fClusterContainer(0),
  fRho(0),
  fLocalRho(0),
  fGeom(0),
  fRunNumber(0)
{
  // Default constructor.

  fClassName = "AliEmcalJet";
}

//________________________________________________________________________
AliJetContainer::AliJetContainer(const char *name):
  AliEmcalContainer(name),
  fJetAcceptanceType(kUser),
  fJetRadius(0),
  fRhoName(),
  fLocalRhoName(),
  fFlavourSelection(0),
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
  fZLeadingEmcCut(10.),
  fZLeadingChCut(10.),
  fNEFMinCut(-10.),
  fNEFMaxCut(10.),
  fLeadingHadronType(0),
  fNLeadingJets(1),
  fJetBitMap(0),
  fJetTrigger(0),
  fParticleContainer(0),
  fClusterContainer(0),
  fRho(0),
  fLocalRho(0),
  fGeom(0),
  fRunNumber(0)
{
  // Standard constructor.

  fClassName = "AliEmcalJet";
}

//________________________________________________________________________
void AliJetContainer::SetArray(AliVEvent *event) 
{
  // Set jet array

  AliEmcalContainer::SetArray(event);

  if(fJetAcceptanceType==kTPC) {
    AliDebug(2,Form("%s: set TPC acceptance cuts",GetName()));
    SetJetEtaPhiTPC();
  }
  else if(fJetAcceptanceType==kEMCAL) {
    AliDebug(2,Form("%s: set EMCAL acceptance cuts",GetName()));
    SetJetEtaPhiEMCAL();
 }
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
  // Load rho

  if (!fRhoName.IsNull() && !fRho) {
    fRho = dynamic_cast<AliRhoParameter*>(event->FindListObject(fRhoName));
    if (!fRho) {
      AliError(Form("%s: Could not retrieve rho %s!", GetName(), fRhoName.Data()));
      return;
    }
  }
}

//________________________________________________________________________
void AliJetContainer::LoadLocalRho(AliVEvent *event)
{
  // Load local rho

  if (!fLocalRhoName.IsNull() && !fLocalRho) {
    fLocalRho = dynamic_cast<AliLocalRhoParameter*>(event->FindListObject(fLocalRhoName));
    if (!fLocalRho) {
      AliError(Form("%s: Could not retrieve rho %s!", GetName(), fLocalRhoName.Data()));
      return;
    }
  }
}

//________________________________________________________________________
AliEmcalJet* AliJetContainer::GetLeadingJet(const char* opt)
{
  // Get the leading jet; if opt contains "rho" the sorting is according to pt-A*rho

  TString option(opt);
  option.ToLower();

  Int_t tempID = fCurrentID;

  AliEmcalJet *jetMax = GetNextAcceptJet(0);
  AliEmcalJet *jet = 0;

  if (option.Contains("rho")) {
    while ((jet = GetNextAcceptJet())) {
      if ( (jet->Pt()-jet->Area()*GetRhoVal()) > (jetMax->Pt()-jetMax->Area()*GetRhoVal()) )
	jetMax = jet;
    }
  }
  else {
    while ((jet = GetNextAcceptJet())) {
      if (jet->Pt() > jetMax->Pt()) jetMax = jet;
    }
  }

  fCurrentID = tempID;

  return jetMax;
}

//________________________________________________________________________
AliEmcalJet* AliJetContainer::GetJet(Int_t i) const {

  //Get i^th jet in array

  if(i<0 || i>fClArray->GetEntriesFast()) return 0;
  AliEmcalJet *jet = static_cast<AliEmcalJet*>(fClArray->At(i));
  return jet;

}

//________________________________________________________________________
AliEmcalJet* AliJetContainer::GetAcceptJet(Int_t i) const {

  //Only return jet if is accepted

  AliEmcalJet *jet = GetJet(i);
  if(!AcceptJet(jet)) return 0;

  return jet;
}

//________________________________________________________________________
AliEmcalJet* AliJetContainer::GetJetWithLabel(Int_t lab) const {

  //Get particle with label lab in array
  
  Int_t i = GetIndexFromLabel(lab);
  return GetJet(i);
}

//________________________________________________________________________
AliEmcalJet* AliJetContainer::GetAcceptJetWithLabel(Int_t lab) const {

  //Get particle with label lab in array
  
  Int_t i = GetIndexFromLabel(lab);
  return GetAcceptJet(i);
}

//________________________________________________________________________
AliEmcalJet* AliJetContainer::GetNextAcceptJet(Int_t i) {

  //Get next accepted jet; if i >= 0 (re)start counter from i; return 0 if no accepted jet could be found

  if (i>=0) fCurrentID = i;

  const Int_t njets = GetNEntries();
  AliEmcalJet *jet = 0;
  while (fCurrentID < njets && !jet) { 
    jet = GetAcceptJet(fCurrentID);
    fCurrentID++;
  }

  return jet;
}

//________________________________________________________________________
AliEmcalJet* AliJetContainer::GetNextJet(Int_t i) {

  //Get next jet; if i >= 0 (re)start counter from i; return 0 if no jet could be found

  if (i>=0) fCurrentID = i;

  const Int_t njets = GetNEntries();
  AliEmcalJet *jet = 0;
  while (fCurrentID < njets && !jet) { 
    jet = GetJet(fCurrentID);
    fCurrentID++;
  }

  return jet;
}

//________________________________________________________________________
Double_t AliJetContainer::GetJetPtCorr(Int_t i) const {
  AliEmcalJet *jet = GetJet(i);

  return jet->Pt() - fRho->GetVal()*jet->Area();
}

//________________________________________________________________________
Double_t AliJetContainer::GetJetPtCorrLocal(Int_t i) const {
  AliEmcalJet *jet = GetJet(i);

  return jet->Pt() - fLocalRho->GetLocalVal(jet->Phi(), fJetRadius)*jet->Area();
}

//________________________________________________________________________
void AliJetContainer::GetMomentum(TLorentzVector &mom, Int_t i) const
{
  //Get momentum of the i^th jet in array

  AliEmcalJet *jet = GetJet(i);
  if(jet) jet->GetMom(mom);
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

   if (jet->AreaEmc() < fAreaEmcCut)
      return kFALSE;
   
   if (fZLeadingChCut < 1 && GetZLeadingCharged(jet) > fZLeadingChCut)
      return kFALSE;
   
   if (fZLeadingEmcCut < 1 && GetZLeadingEmc(jet) > fZLeadingEmcCut)
      return kFALSE;

   if (jet->NEF() < fNEFMinCut || jet->NEF() > fNEFMaxCut)
      return kFALSE;
   
   if (!AcceptBiasJet(jet))
      return kFALSE;
   
   if (jet->MaxTrackPt() > fMaxTrackPt || jet->MaxClusterPt() > fMaxClusterPt)
      return kFALSE;
   
   if (fFlavourSelection != 0 && !jet->TestFlavourTag(fFlavourSelection))
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
void AliJetContainer::GetLeadingHadronMomentum(TLorentzVector &mom, AliEmcalJet *jet) const
{
  Double_t maxClusterPt = 0;
  Double_t maxClusterEta = 0;
  Double_t maxClusterPhi = 0;
  
  Double_t maxTrackPt = 0;
  Double_t maxTrackEta = 0;
  Double_t maxTrackPhi = 0;
      
  if (fClusterContainer && fClusterContainer->GetArray() && (fLeadingHadronType == 1 || fLeadingHadronType == 2)) {
    AliVCluster *cluster = jet->GetLeadingCluster(fClusterContainer->GetArray());
    if (cluster) {
      TLorentzVector nPart;
      cluster->GetMomentum(nPart, const_cast<Double_t*>(fVertex));
      
      maxClusterEta = nPart.Eta();
      maxClusterPhi = nPart.Phi();
      maxClusterPt = nPart.Pt();
    }
  }
      
  if (fParticleContainer && fParticleContainer->GetArray() && (fLeadingHadronType == 0 || fLeadingHadronType == 2)) {
    AliVParticle *track = jet->GetLeadingTrack(fParticleContainer->GetArray());
    if (track) {
      maxTrackEta = track->Eta();
      maxTrackPhi = track->Phi();
      maxTrackPt = track->Pt();
    }
  }
      
  if (maxTrackPt > maxClusterPt) 
    mom.SetPtEtaPhiM(maxTrackPt,maxTrackEta,maxTrackPhi,0.139);
  else 
    mom.SetPtEtaPhiM(maxClusterPt,maxClusterEta,maxClusterPhi,0.139);
}

//________________________________________________________________________
Double_t AliJetContainer::GetZLeadingEmc(AliEmcalJet *jet) const
{

  if (fClusterContainer && fClusterContainer->GetArray()) {
    TLorentzVector mom;
    
    AliVCluster *cluster = jet->GetLeadingCluster(fClusterContainer->GetArray());
    if (cluster) {
      cluster->GetMomentum(mom, const_cast<Double_t*>(fVertex));
      
      return GetZ(jet,mom);
    }
    else
      return -1;
  }
  else
    return -1;
}

//________________________________________________________________________
Double_t AliJetContainer::GetZLeadingCharged(AliEmcalJet *jet) const
{

  if (fParticleContainer && fParticleContainer->GetArray() ) {
    TLorentzVector mom;
    
    AliVParticle *track = jet->GetLeadingTrack(fParticleContainer->GetArray());
    if (track) {
      mom.SetPtEtaPhiM(track->Pt(),track->Eta(),track->Phi(),0.139);
      
      return GetZ(jet,mom);
    }
    else
      return -1;
  }
  else
    return -1;
}

//________________________________________________________________________
Double_t AliJetContainer::GetZ(AliEmcalJet *jet, TLorentzVector mom) const
{

  Double_t pJetSq = jet->Px()*jet->Px() + jet->Py()*jet->Py() + jet->Pz()*jet->Pz();

  if(pJetSq>1e-6)
    return (mom.Px()*jet->Px() + mom.Py()*jet->Py() + mom.Pz()*jet->Pz())/pJetSq;
  else {
    AliWarning(Form("%s: strange, pjet*pjet seems to be zero pJetSq: %f",GetName(), pJetSq));
    return -1;
  }

}

//________________________________________________________________________
void AliJetContainer::SetJetEtaPhiEMCAL()
{
  //Set default cuts for full jets

  if(!fGeom) SetEMCALGeometry();
  if(fGeom) {
    SetJetEtaLimits(fGeom->GetArm1EtaMin() + fJetRadius, fGeom->GetArm1EtaMax() - fJetRadius);

    if(fRunNumber>=177295 && fRunNumber<=197470) //small SM masked in 2012 and 2013
      SetJetPhiLimits(1.4+fJetRadius,TMath::Pi()-fJetRadius);
    else
      SetJetPhiLimits(fGeom->GetArm1PhiMin() * TMath::DegToRad() + fJetRadius, fGeom->GetArm1PhiMax() * TMath::DegToRad() - fJetRadius);

  }
  else {
    AliWarning("Could not get instance of AliEMCALGeometry. Using manual settings for EMCAL year 2011!!");
    SetJetEtaLimits(-0.7+fJetRadius,0.7-fJetRadius);
    SetJetPhiLimits(1.405+fJetRadius,3.135-fJetRadius);
  }
}

//________________________________________________________________________
void AliJetContainer::SetJetEtaPhiTPC()
{
  //Set default cuts for charged jets

  SetJetEtaLimits(-0.9+fJetRadius, 0.9-fJetRadius);
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
  fZLeadingEmcCut = 10.;
  fZLeadingChCut = 10.;
}

//________________________________________________________________________
void AliJetContainer::SetClassName(const char *clname)
{
  // Set the class name

  TClass cls(clname);
  if (cls.InheritsFrom("AliEmcalJet")) fClassName = clname;
  else AliError(Form("Unable to set class name %s for a AliJetContainer, it must inherits from AliEmcalJet!",clname));
}

//________________________________________________________________________
Double_t AliJetContainer::GetFractionSharedPt(AliEmcalJet *jet1) const
{
  //
  // Get fraction of shared pT between matched full and charged jet
  // Uses charged jet pT as baseline: fraction = \Sum_{const,full jet} pT,const,i / pT,jet,ch
  // Only works if tracks array of both jets is the same
  //

  AliEmcalJet *jet2 = jet1->ClosestJet();
  if(!jet2) return -1;

  Double_t fraction = 0.;
  Double_t jetPt2 = jet2->Pt();
 
  if(jetPt2>0) {
    Double_t sumPt = 0.;
    AliVParticle *vpf = 0x0;
    Int_t iFound = 0;
    for(Int_t icc=0; icc<jet2->GetNumberOfTracks(); icc++) {
      Int_t idx = (Int_t)jet2->TrackAt(icc);
      iFound = 0;
      for(Int_t icf=0; icf<jet1->GetNumberOfTracks(); icf++) {
	if(idx == jet1->TrackAt(icf) && iFound==0 ) {
	  iFound=1;
	  vpf = static_cast<AliVParticle*>(jet1->TrackAt(icf, fParticleContainer->GetArray()));
	  if(vpf) sumPt += vpf->Pt();
	  continue;
	}
      }
    }
    fraction = sumPt/jetPt2;
  } else 
    fraction = -1;
  
  return fraction;
}

