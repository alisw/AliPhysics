//
// Container with name, TClonesArray and cuts for jets
//
// Author: M. Verweij, S. Aiola

#include <TClonesArray.h>

#include "AliVEvent.h"
#include "AliLog.h"
#include "AliEMCALGeometry.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliLocalRhoParameter.h"
#include "AliTLorentzVector.h"

#include "AliJetContainer.h"

ClassImp(AliJetContainer)

//________________________________________________________________________
AliJetContainer::AliJetContainer():
  AliParticleContainer(),
  fJetAcceptanceType(kUser),
  fJetRadius(0),
  fRhoName(),
  fLocalRhoName(),
  fRhoMassName(),
  fFlavourSelection(0),
  fJetAreaCut(-1),
  fAreaEmcCut(-1),
  fMinClusterPt(-1),
  fMaxClusterPt(1000),
  fMinTrackPt(-1),
  fMaxTrackPt(100),
  fZLeadingEmcCut(10.),
  fZLeadingChCut(10.),
  fNEFMinCut(-10.),
  fNEFMaxCut(10.),
  fLeadingHadronType(0),
  fNLeadingJets(1),
  fMinNConstituents(-1),
  fJetTrigger(0),
  fTagStatus(-1),
  fParticleContainer(0),
  fClusterContainer(0),
  fRho(0),
  fLocalRho(0),
  fRhoMass(0),
  fGeom(0),
  fRunNumber(0)
{
  // Default constructor.

  fBaseClassName = "AliEmcalJet";
  SetClassName("AliEmcalJet");
}

//________________________________________________________________________
AliJetContainer::AliJetContainer(const char *name):
  AliParticleContainer(name),
  fJetAcceptanceType(kUser),
  fJetRadius(0),
  fRhoName(),
  fLocalRhoName(),
  fRhoMassName(),
  fFlavourSelection(0),
  fJetAreaCut(-1),
  fAreaEmcCut(-1),
  fMinClusterPt(-1),
  fMaxClusterPt(1000),
  fMinTrackPt(-1),
  fMaxTrackPt(100),
  fZLeadingEmcCut(10.),
  fZLeadingChCut(10.),
  fNEFMinCut(-10.),
  fNEFMaxCut(10.),
  fLeadingHadronType(0),
  fNLeadingJets(1),
  fMinNConstituents(-1),
  fJetTrigger(0),
  fTagStatus(-1),
  fParticleContainer(0),
  fClusterContainer(0),
  fRho(0),
  fLocalRho(0),
  fRhoMass(0),
  fGeom(0),
  fRunNumber(0)
{
  // Standard constructor.

  fBaseClassName = "AliEmcalJet";
  SetClassName("AliEmcalJet");
  SetMinPt(1);
}

//________________________________________________________________________
AliJetContainer::AliJetContainer(EJetType_t jetType, EJetAlgo_t jetAlgo, ERecoScheme_t recoScheme, Double_t radius,
    AliParticleContainer* partCont, AliClusterContainer* clusCont, TString tag):
  AliParticleContainer(GenerateJetName(jetType, jetAlgo, recoScheme, radius, partCont, clusCont, tag)),
  fJetAcceptanceType(kUser),
  fJetRadius(radius),
  fRhoName(),
  fLocalRhoName(),
  fRhoMassName(),
  fFlavourSelection(0),
  fJetAreaCut(-1),
  fAreaEmcCut(-1),
  fMinClusterPt(-1),
  fMaxClusterPt(1000),
  fMinTrackPt(-1),
  fMaxTrackPt(100),
  fZLeadingEmcCut(10.),
  fZLeadingChCut(10.),
  fNEFMinCut(-10.),
  fNEFMaxCut(10.),
  fLeadingHadronType(0),
  fNLeadingJets(1),
  fMinNConstituents(-1),
  fJetTrigger(0),
  fTagStatus(-1),
  fParticleContainer(partCont),
  fClusterContainer(clusCont),
  fRho(0),
  fLocalRho(0),
  fRhoMass(0),
  fGeom(0),
  fRunNumber(0)
{
  // Constructor.

  fBaseClassName = "AliEmcalJet";
  SetClassName("AliEmcalJet");
  SetMinPt(1);
}

//________________________________________________________________________
void AliJetContainer::SetArray(AliVEvent *event) 
{
  // Set jet array

  AliEmcalContainer::SetArray(event);

  SetAcceptanceCuts();
}

//________________________________________________________________________
void AliJetContainer::SetAcceptanceCuts()
{
  // Set acceptance

  switch (fJetAcceptanceType) {
  case kTPC:
    AliDebug(2,Form("%s: set TPC acceptance cuts",GetName()));
    SetJetEtaPhiTPC();
    break;
  case kTPCfid:
    AliDebug(2,Form("%s: set TPC acceptance cuts",GetName()));
    SetJetEtaPhiTPC(fJetRadius);
    break;
  case kEMCAL:
    AliDebug(2,Form("%s: set EMCAL acceptance cuts",GetName()));
    SetJetEtaPhiEMCAL();
  case kEMCALfid:
    AliDebug(2,Form("%s: set EMCAL acceptance cuts",GetName()));
    SetJetEtaPhiEMCAL(fJetRadius);
    break;
  case kDCALfid:
    AliDebug(2,Form("%s: set EMCAL acceptance cuts",GetName()));
    SetJetEtaPhiDCAL();
    break;
  case kDCAL:
    AliDebug(2,Form("%s: set EMCAL acceptance cuts",GetName()));
    SetJetEtaPhiDCAL(fJetRadius);
    break;
  case kUser:
    break;
  }
}


//________________________________________________________________________
void AliJetContainer::SetEMCALGeometry()
{
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
void AliJetContainer::LoadRhoMass(AliVEvent *event)
{
  // Load rho

  if (!fRhoMassName.IsNull() && !fRhoMass) {
    fRhoMass = dynamic_cast<AliRhoParameter*>(event->FindListObject(fRhoMassName));
    if (!fRhoMass) {
      AliError(Form("%s: Could not retrieve rho_mass %s!", GetName(), fRhoMassName.Data()));
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
  ResetCurrentID();

  AliEmcalJet *jetMax = GetNextAcceptJet();
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

  UInt_t rejectionReason = 0;
  AliEmcalJet *jet = GetJet(i);
  if(!AcceptJet(jet, rejectionReason)) return 0;

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
AliEmcalJet* AliJetContainer::GetNextAcceptJet() {

  //Get next accepted jet;

  const Int_t njets = GetNEntries();
  AliEmcalJet *jet = 0;
  do {
    fCurrentID++;
    if (fCurrentID >= njets) break;
    jet = GetAcceptJet(fCurrentID);
  } while (!jet);

  return jet;
}

//________________________________________________________________________
AliEmcalJet* AliJetContainer::GetNextJet() {

  //Get next jet;

  const Int_t njets = GetNEntries();
  AliEmcalJet *jet = 0;
  do {
    fCurrentID++;
    if (fCurrentID >= njets) break;
    jet = GetJet(fCurrentID);
  } while (!jet);


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
Bool_t AliJetContainer::GetMomentumFromJet(TLorentzVector &mom, const AliEmcalJet* jet, Double_t mass) const
{
  Double_t p = jet->P();
  Double_t e = TMath::Sqrt(mass*mass + p*p);

  mom.SetPtEtaPhiE(jet->Pt(), jet->Eta(), jet->Phi(), e);

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliJetContainer::GetMomentumFromJet(TLorentzVector &mom, const AliEmcalJet* jet) const
{
  if (jet) {
    if (fMassHypothesis >= 0) {
      GetMomentumFromJet(mom, jet, fMassHypothesis);
    }
    else {
      jet->GetMomentum(mom);
    }
    return kTRUE;
  }
  else {
    mom.SetPtEtaPhiM(0, 0, 0, 0);
    return kFALSE;
  }
}

//________________________________________________________________________
Bool_t AliJetContainer::GetMomentum(TLorentzVector &mom, Int_t i) const
{
  //Get momentum of the i^th particle in array

  AliEmcalJet *jet = GetJet(i);
  return GetMomentumFromJet(mom, jet);
}

//________________________________________________________________________
Bool_t AliJetContainer::GetNextMomentum(TLorentzVector &mom)
{
  //Get momentum of the next jet in array

  AliEmcalJet *jet = GetNextJet();
  return GetMomentumFromJet(mom, jet);
}

//________________________________________________________________________
Bool_t AliJetContainer::GetAcceptMomentum(TLorentzVector &mom, Int_t i) const
{
  //Get momentum of the i^th jet in array

  AliEmcalJet *jet = GetAcceptJet(i);
  return GetMomentumFromJet(mom, jet);
}

//________________________________________________________________________
Bool_t AliJetContainer::GetNextAcceptMomentum(TLorentzVector &mom)
{
  //Get momentum of the next accepted jet in array

  AliEmcalJet *jet = GetNextAcceptJet();
  return GetMomentumFromJet(mom, jet);
}

//________________________________________________________________________
Bool_t AliJetContainer::AcceptJet(const AliEmcalJet *jet, UInt_t &rejectionReason) const
{
  // Return true if jet is accepted.

  Bool_t r = ApplyJetCuts(jet, rejectionReason);
  if (!r) return kFALSE;

  AliTLorentzVector mom;
  GetMomentumFromJet(mom, jet);

  return ApplyKinematicCuts(mom, rejectionReason);
}

//________________________________________________________________________
Bool_t AliJetContainer::AcceptJet(Int_t i, UInt_t &rejectionReason) const
{
  // Return true if jet is accepted.

  Bool_t r = ApplyJetCuts(GetJet(i), rejectionReason);
  if (!r) return kFALSE;

  AliTLorentzVector mom;
  GetMomentum(mom, i);

  return ApplyKinematicCuts(mom, rejectionReason);
}


//________________________________________________________________________
Bool_t AliJetContainer::ApplyJetCuts(const AliEmcalJet *jet, UInt_t &rejectionReason) const
{   
  // Return true if jet is accepted.

  if (!jet) {
    AliDebug(11,"No jet found");
    rejectionReason |= kNullObject;
    return kFALSE;
  }

  if (jet->TestBits(fBitMap) != (Int_t)fBitMap) {
    AliDebug(11,"Cut rejecting jet: Bit map");
    rejectionReason |= kBitMapCut;
    return kFALSE;
  }

  if (jet->Area() <= fJetAreaCut)  {
    AliDebug(11,"Cut rejecting jet: Area");
    rejectionReason |= kAreaCut;
    return kFALSE;
  }

  if (jet->AreaEmc() < fAreaEmcCut) {
    AliDebug(11,"Cut rejecting jet: AreaEmc");
    rejectionReason |= kAreaEmcCut;
    return kFALSE;
  }

  if (fZLeadingChCut < 1 && GetZLeadingCharged(jet) > fZLeadingChCut) {
    AliDebug(11,"Cut rejecting jet: ZLeading");
    rejectionReason |= kZLeadingChCut;
    return kFALSE;
  }

  if (fZLeadingEmcCut < 1 && GetZLeadingEmc(jet) > fZLeadingEmcCut) {
    AliDebug(11,"Cut rejecting jet: ZLeadEmc");
    rejectionReason |= kZLeadingEmcCut;
    return kFALSE;
  }

  if (jet->NEF() < fNEFMinCut || jet->NEF() > fNEFMaxCut) {
    AliDebug(11,"Cut rejecting jet: NEF");
    rejectionReason |= kNEFCut;
    return kFALSE;
  }

  if(fMinNConstituents>0 && jet->GetNumberOfConstituents()<fMinNConstituents) {
    AliDebug(11,"Cut rejecting jet: minimum number of constituents");
    rejectionReason |= kMinNConstituents;
    return kFALSE;
  }

  if (fLeadingHadronType == 0) {
    if (jet->MaxTrackPt() < fMinTrackPt) {
      AliDebug(11,"Cut rejecting jet: Bias");
      rejectionReason |= kMinLeadPtCut;
      return kFALSE;
    }
  }
  else if (fLeadingHadronType == 1) {
    if (jet->MaxClusterPt() < fMinClusterPt) {
      AliDebug(11,"Cut rejecting jet: Bias");
      rejectionReason |= kMinLeadPtCut;
      return kFALSE;
    }
  }
  else {
    if (jet->MaxTrackPt() < fMinTrackPt && jet->MaxClusterPt() < fMinClusterPt) {
      AliDebug(11,"Cut rejecting jet: Bias");
      rejectionReason |= kMinLeadPtCut;
      return kFALSE;
    }
  }

  if (jet->MaxTrackPt() > fMaxTrackPt) {
    AliDebug(11,"Cut rejecting jet: MaxTrackPt");
    rejectionReason |= kMaxTrackPtCut;
    return kFALSE;

  }

  if (jet->MaxClusterPt() > fMaxClusterPt) {
    AliDebug(11,"Cut rejecting jet: MaxClusPt");
    rejectionReason |= kMaxClusterPtCut;
    return kFALSE;
  }

  if (fFlavourSelection != 0 && !jet->TestFlavourTag(fFlavourSelection)) {
    AliDebug(11,"Cut rejecting jet: Flavour");
    rejectionReason |= kFlavourCut;
    return kFALSE;
  }

  if (fTagStatus>-1 && jet->GetTagStatus()!=fTagStatus) {
    AliDebug(11,"Cut rejecting jet: tag status");
    rejectionReason |= kTagStatus;
    return kFALSE;
  }

  return kTRUE;
}

//________________________________________________________________________
Double_t AliJetContainer::GetLeadingHadronPt(const AliEmcalJet *jet) const
{
  if (fLeadingHadronType == 0)       // charged leading hadron
    return jet->MaxTrackPt();
  else if (fLeadingHadronType == 1)  // neutral leading hadron
    return jet->MaxClusterPt();
  else                               // charged or neutral
    return jet->MaxPartPt();
}

//________________________________________________________________________
void AliJetContainer::GetLeadingHadronMomentum(TLorentzVector &mom, const AliEmcalJet *jet) const
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
Double_t AliJetContainer::GetZLeadingEmc(const AliEmcalJet *jet) const
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
Double_t AliJetContainer::GetZLeadingCharged(const AliEmcalJet *jet) const
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
Double_t AliJetContainer::GetZ(const AliEmcalJet *jet, TLorentzVector mom) const
{
  Double_t pJetSq = jet->Px()*jet->Px() + jet->Py()*jet->Py() + jet->Pz()*jet->Pz();

  if (pJetSq < 1e-6) {
    AliWarning(Form("%s: strange, pjet*pjet seems to be zero pJetSq: %.3f",GetName(), pJetSq));
    return 0;
  }

  Double_t z = (mom.Px()*jet->Px() + mom.Py()*jet->Py() + mom.Pz()*jet->Pz()) / pJetSq;

  if (z < 0) {
    AliWarning(Form("%s: z  = %.3ff < 0, returning 0...",GetName(), z));
    z = 0;
  }

  return z;
}

//________________________________________________________________________
void AliJetContainer::SetJetEtaPhiEMCAL(Double_t r)
{
  //Set default cuts for full jets in EMCal

  if (!fGeom) SetEMCALGeometry();
  if (fGeom) {
    SetEtaLimits(fGeom->GetArm1EtaMin() + r, fGeom->GetArm1EtaMax() - r);

    if(fRunNumber>=177295 && fRunNumber<=197470) {//small SM masked in 2012 and 2013
      SetPhiLimits(1.405 + r,3.135 - r);
    }
    else {
      SetPhiLimits(fGeom->GetArm1PhiMin() * TMath::DegToRad() + r, fGeom->GetEMCALPhiMax() * TMath::DegToRad() - r);
    }
  }
  else {
    AliWarning("Could not get instance of AliEMCALGeometry. Using manual settings for EMCAL year 2011!!");
    SetEtaLimits(-0.7 + r, 0.7 - r);
    SetPhiLimits(1.405 + r, 3.135 - r);
  }
}

//________________________________________________________________________
void AliJetContainer::SetJetEtaPhiDCAL(Double_t r)
{
  //Set default cuts for full jets in DCal

  if (!fGeom) SetEMCALGeometry();
  if (fGeom) {
    SetEtaLimits(fGeom->GetArm1EtaMin() + r, fGeom->GetArm1EtaMax() - r);
    SetPhiLimits(fGeom->GetDCALPhiMin() * TMath::DegToRad() + r, fGeom->GetDCALPhiMax() * TMath::DegToRad() - r);
  }
  else {
    AliWarning("Could not get instance of AliEMCALGeometry. Using manual settings for DCAL year 2015!!");
    SetEtaLimits(-0.7 + r, 0.7 - r);
    SetPhiLimits(4.538 + r, 5.727 - r);
  }
}

//________________________________________________________________________
void AliJetContainer::SetJetEtaPhiTPC(Double_t r)
{
  //Set default cuts for charged jets

  SetEtaLimits(-0.9 + r, 0.9 - r);
  SetPhiLimits(0, 0);  // No cut on phi
}

//________________________________________________________________________
void AliJetContainer::PrintCuts() 
{
  // Reset cuts to default values
  TString arrName = GetArrayName();
  Printf("Print jet cuts for %s",arrName.Data());
  Printf("PtBiasJetTrack: %f",fMinTrackPt);
  Printf("PtBiasJetClus: %f",fMinClusterPt);
  Printf("JetPtCut: %f", fMinPt);
  Printf("JetPtCutMax: %f", fMaxPt);
  Printf("JetAreaCut: %f",fJetAreaCut);
  Printf("AreaEmcCut: %f",fAreaEmcCut);
  Printf("JetMinEta: %f", fMinEta);
  Printf("JetMaxEta: %f", fMaxEta);
  Printf("JetMinPhi: %f", fMinPhi);
  Printf("JetMaxPhi: %f", fMaxPhi);
  Printf("MaxClusterPt: %f",fMaxClusterPt);
  Printf("MaxTrackPt: %f",fMaxTrackPt);
  Printf("LeadingHadronType: %d",fLeadingHadronType);
  Printf("ZLeadingEmcCut: %f",fZLeadingEmcCut);
  Printf("ZLeadingChCut: %f",fZLeadingChCut);

}

//________________________________________________________________________
void AliJetContainer::ResetCuts() 
{
  // Reset cuts to default values

  fMinTrackPt     = 0;
  fMinClusterPt   = 0;
  fMinPt          = 0;
  fJetAreaCut     = -1;
  fAreaEmcCut     = -1;
  fMinEta         = -0.9;
  fMaxEta         = 0.9;
  fMinPhi         = 0;
  fMaxPhi         = 10;
  fMaxClusterPt   = 1000;
  fMaxTrackPt     = 100;
  fLeadingHadronType = 0;
  fZLeadingEmcCut = 10.;
  fZLeadingChCut  = 10.;
}

//________________________________________________________________________
Int_t AliJetContainer::GetNAcceptedJets()
{
  // Get number of accepted jets

  Int_t nJet = 0;
  Int_t tempID = fCurrentID;
  ResetCurrentID();

  AliEmcalJet *jet = GetNextAcceptJet();
  if(jet) nJet = 1;
  while (GetNextAcceptJet())
    nJet++;

  fCurrentID = tempID;

  return nJet;
}

//________________________________________________________________________
Double_t AliJetContainer::GetFractionSharedPt(const AliEmcalJet *jet1, AliParticleContainer *cont2) const
{
  //
  // Get fraction of shared pT between matched jets
  // Uses ClosestJet() jet pT as baseline: fraction = \Sum_{const,jet1} pT,const,i / pT,jet,closest
  // Only works if tracks array of both jets is the same -> modifying this. if no container is given than is like this, otherwise the geomerical matching is applied
  //

  AliEmcalJet *jet2 = jet1->ClosestJet();
  if(!jet2) return -1;

  Double_t fraction = 0.;
  Double_t jetPt2 = jet2->Pt();
  Int_t bgeom = kTRUE;
  if(!cont2) bgeom = kFALSE;
  if(jetPt2>0) {
    Double_t sumPt = 0.;
    AliVParticle *vpf = 0x0;
    Int_t iFound = 0;
    for(Int_t icc=0; icc<jet2->GetNumberOfTracks(); icc++) {
      Int_t idx = (Int_t)jet2->TrackAt(icc);
      //get particle
      AliVParticle *p2 = 0x0;
      if(bgeom) p2 = static_cast<AliVParticle*>(jet2->TrackAt(icc, cont2->GetArray()));
      iFound = 0;
      for(Int_t icf=0; icf<jet1->GetNumberOfTracks(); icf++) {
        if(!bgeom && idx == jet1->TrackAt(icf) && iFound==0 ) {
          iFound=1;
          vpf = static_cast<AliVParticle*>(jet1->TrackAt(icf, fParticleContainer->GetArray()));
          if(vpf) sumPt += vpf->Pt();
          continue;
        }
        if(bgeom){
          vpf = static_cast<AliVParticle*>(jet1->TrackAt(icf, fParticleContainer->GetArray()));
          if(!vpf) continue;
          if(!SamePart(vpf, p2, 1.e-4)) continue; //not the same particle
          sumPt += vpf->Pt();
        }
      }
    }
    fraction = sumPt/jetPt2;
  } else 
    fraction = -1;
  return fraction;
}

//________________________________________________________________________
TString AliJetContainer::GenerateJetName(EJetType_t jetType, EJetAlgo_t jetAlgo, ERecoScheme_t recoScheme, Double_t radius, AliParticleContainer* partCont, AliClusterContainer* clusCont, TString tag)
{
  TString algoString;
  switch (jetAlgo)
  {
  case kt_algorithm:
    algoString = "KT";
    break;
  case antikt_algorithm:
    algoString = "AKT";
    break;
  default:
    ::Warning("AliJetContainer::GenerateJetName", "Unknown jet finding algorithm '%d'!", jetAlgo);
    algoString = "";
  }

  TString typeString;
  switch (jetType) {
  case kFullJet:
    typeString = "Full";
    break;
  case kChargedJet:
    typeString = "Charged";
    break;
  case kNeutralJet:
    typeString = "Neutral";
    break;
  }

  TString radiusString = TString::Format("R%03.0f", radius*100.0);

  TString trackString;
  if (jetType != kNeutralJet && partCont) {
    trackString = "_" + TString(partCont->GetTitle());
  }

  TString clusterString;
  if (jetType != kChargedJet && clusCont) {
    clusterString = "_" + TString(clusCont->GetTitle());
  }

  TString recombSchemeString;
  switch (recoScheme) {
  case E_scheme:
    recombSchemeString = "E_scheme";
    break;
  case pt_scheme:
    recombSchemeString = "pt_scheme";
    break;
  case pt2_scheme:
    recombSchemeString = "pt2_scheme";
    break;
  case Et_scheme:
    recombSchemeString = "Et_scheme";
    break;
  case Et2_scheme:
    recombSchemeString = "Et2_scheme";
    break;
  case BIpt_scheme:
    recombSchemeString = "BIpt_scheme";
    break;
  case BIpt2_scheme:
    recombSchemeString = "BIpt2_scheme";
    break;
  case external_scheme:
    recombSchemeString = "ext_scheme";
    break;
  default:
    ::Error("AliJetContainer::GenerateJetName", "Recombination %d scheme not recognized.", recoScheme);
  }

  TString name = TString::Format("%s_%s%s%s%s%s_%s",
      tag.Data(), algoString.Data(), typeString.Data(), radiusString.Data(), trackString.Data(), clusterString.Data(), recombSchemeString.Data());

  return name;
}

//________________________________________________________________________
const char* AliJetContainer::GetTitle() const
{
  static TString jetString;

  if (GetMinPt() == 0) {
    jetString = TString::Format("_%s_pT0000", GetArrayName().Data());
  }
  else if (GetMinPt() < 1.0) {
    jetString = TString::Format("_%s_pT0%3.0f", GetArrayName().Data(), GetMinPt()*1000.0);
  }
  else {
    jetString = TString::Format("_%s_pT%4.0f", GetArrayName().Data(), GetMinPt()*1000.0);
  }

  return jetString.Data();
}
