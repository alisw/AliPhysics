/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include <TClonesArray.h>

#include "AliVEvent.h"
#include "AliLog.h"
#include "AliEMCALGeometry.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliLocalRhoParameter.h"
#include "AliTLorentzVector.h"

#include "AliJetContainer.h"

/// \cond CLASSIMP
ClassImp(AliJetContainer);
/// \endcond

/**
 * Default constructor.
 */
AliJetContainer::AliJetContainer():
  AliParticleContainer(),
  fJetAcceptanceType(0),
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
  fRunNumber(0),
  fTpcHolePos(0),
  fTpcHoleWidth(0)
{
  fBaseClassName = "AliEmcalJet";
  SetClassName("AliEmcalJet");
}

/**
 * Standard constructor.
 * @param name Name of the jet branch (TClonesArray)
 */
AliJetContainer::AliJetContainer(const char *name):
  AliParticleContainer(name),
  fJetAcceptanceType(0),
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
  fRunNumber(0),
  fTpcHolePos(0),
  fTpcHoleWidth(0)
{
  fBaseClassName = "AliEmcalJet";
  SetClassName("AliEmcalJet");
  SetMinPt(1);
}

/**
 * Jet definition constructor.
 *
 * This constructor takes a jet definition as an input and builds automatically the
 * jet branch name for itself.
 * @param jetType Type of the jet (full, charged, neutral)
 * @param jetAlgo Jet algorithm used to reconstruct the jets (anti-kt, kt, etc.)
 * @param recoScheme Jet recombination scheme (E-scheme, pt-scheme)
 * @param radius Jet resolution parameter
 * @param partCont Particle container used to feed the jet finder
 * @param clustCont Cluster container used to feed the jet finder
 * @param tag Additional tag for the jet branch name (default is "Jet")
 */
AliJetContainer::AliJetContainer(EJetType_t jetType, EJetAlgo_t jetAlgo, ERecoScheme_t recoScheme, Double_t radius,
    AliParticleContainer* partCont, AliClusterContainer* clusCont, TString tag):
  AliParticleContainer(GenerateJetName(jetType, jetAlgo, recoScheme, radius, partCont, clusCont, tag)),
  fJetAcceptanceType(0),
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
  fBaseClassName = "AliEmcalJet";
  SetClassName("AliEmcalJet");
  SetMinPt(1);
}

/**
 * Calls the base class method, then set the acceptance cuts.
 * @param event Event pointer used to retrieve the jet branch
 */
void AliJetContainer::SetArray(const AliVEvent *event)
{
  // Set jet array

  AliEmcalContainer::SetArray(event);
}

/**
 * Loads the rho object from the provided event.
 * The rho object contains information about the event average
 * energy density, used to subtract diffuse background in jet reconstructed
 * in pA or A-A collisions.
 * @param Valid pointer to a AliVEvent object from which the object is to be retrieved
 */
void AliJetContainer::LoadRho(const AliVEvent *event)
{
  if (!fRhoName.IsNull() && !fRho) {
    fRho = dynamic_cast<AliRhoParameter*>(event->FindListObject(fRhoName));
    if (!fRho) {
      AliError(Form("%s: Could not retrieve rho %s!", GetName(), fRhoName.Data()));
      return;
    }
  }
}

/**
 * Loads the local rho object from the provided event.
 * The local rho object contains information about the event
 * energy density as a function of the event geometry.
 * It is used to subtract background in jet reconstructed
 * in pA or A-A collisions.
 * @param Valid pointer to a AliVEvent object from which the object is to be retrieved
 */
void AliJetContainer::LoadLocalRho(const AliVEvent *event)
{
  if (!fLocalRhoName.IsNull() && !fLocalRho) {
    fLocalRho = dynamic_cast<AliLocalRhoParameter*>(event->FindListObject(fLocalRhoName));
    if (!fLocalRho) {
      AliError(Form("%s: Could not retrieve rho %s!", GetName(), fLocalRhoName.Data()));
      return;
    }
  }
}

/**
 * Loads the mass rho object from the provided event.
 * The mass rho object contains information about the average event
 * mass density.
 * It is used to subtract background in jet reconstructed
 * in pA or A-A collisions.
 * @param Valid pointer to a AliVEvent object from which the object is to be retrieved
 */
void AliJetContainer::LoadRhoMass(const AliVEvent *event)
{
  if (!fRhoMassName.IsNull() && !fRhoMass) {
    fRhoMass = dynamic_cast<AliRhoParameter*>(event->FindListObject(fRhoMassName));
    if (!fRhoMass) {
      AliError(Form("%s: Could not retrieve rho_mass %s!", GetName(), fRhoMassName.Data()));
      return;
    }
  }
}

/**
 * Finds the leading jet in the container according to pT.
 * If opt contains "rho" the sorting is done according to pT - area x rho instead of pT
 * @param opt Options concerning the sorting. The only implemented option is "rho".
 * @return A pointer to the leading jet in the container
 */
AliEmcalJet* AliJetContainer::GetLeadingJet(const char* opt)
{
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

/**
 * Finds the jet at position i in the container.
 * @param i Index position of the jet
 * @return A pointer to the jet object at position i
 */
AliEmcalJet* AliJetContainer::GetJet(Int_t i) const
{
  if (i < 0 || i > fClArray->GetEntriesFast()) return 0;
  AliEmcalJet *jet = static_cast<AliEmcalJet*>(fClArray->At(i));
  return jet;

}

/**
 * Finds the jet at position i in the container and checks if the jet passes the cuts.
 * @param i Index position of the jet
 * @return A pointer to the jet object at position i if the jet passes the cuts; if the jet
 * does not pass the cuts, a NULL pointer
 */
AliEmcalJet* AliJetContainer::GetAcceptJet(Int_t i) const
{
  UInt_t rejectionReason = 0;
  AliEmcalJet *jet = GetJet(i);
  if(!AcceptJet(jet, rejectionReason)) return 0;

  return jet;
}

/**
 * Iterator over accepted jets in the container. Get the next accepted
 * jet in the array. If the end is reached, NULL is returned.
 * @deprecated Only for backward compatibility - use AliJetIterableContainer instead
 * @return Next accepted jet in the array (NULL if the end is reached)
 */
AliEmcalJet* AliJetContainer::GetNextAcceptJet()
{
  const Int_t njets = GetNEntries();
  AliEmcalJet *jet = 0;
  do {
    fCurrentID++;
    if (fCurrentID >= njets) break;
    jet = GetAcceptJet(fCurrentID);
  } while (!jet);

  return jet;
}

/**
 * Iterator over jets in the container. Get the next
 * jet in the array. If the end is reached, NULL is returned.
 * @deprecated Only for backward compatibility - use AliJetIterableContainer instead
 * @return Next jet in the array (NULL if the end is reached)
 */
AliEmcalJet* AliJetContainer::GetNextJet()
{
  const Int_t njets = GetNEntries();
  AliEmcalJet *jet = 0;
  do {
    fCurrentID++;
    if (fCurrentID >= njets) break;
    jet = GetJet(fCurrentID);
  } while (!jet);


  return jet;
}

/**
 * Finds the jet at position i in the container
 * and subtracts the average background from the jet pT.
 * @param i Index position of the jet
 * @return The subtracted jet pT
 */
Double_t AliJetContainer::GetJetPtCorr(Int_t i) const
{
  AliEmcalJet *jet = GetJet(i);
  return jet->Pt() - fRho->GetVal()*jet->Area();
}

/**
 * Finds the jet at position i in the container
 * and subtracts the event geometry modulated background from the jet pT.
 * @param i Index position of the jet
 * @return The subtracted jet pT
 */
Double_t AliJetContainer::GetJetPtCorrLocal(Int_t i) const
{
  AliEmcalJet *jet = GetJet(i);

  return jet->Pt() - fLocalRho->GetLocalVal(jet->Phi(), fJetRadius)*jet->Area();
}

/**
 * Calculates the jet 4-momentum assuming a given mass.
 * @param[out] mom Reference to a TLorentzVector object where the 4-momentum is returned
 * @param[in] jet A valid pointer to an AliEmcalJet object
 * @param[in] mass Mass assumption for the jet
 * @return Always kTRUE
 */
Bool_t AliJetContainer::GetMomentumFromJet(TLorentzVector &mom, const AliEmcalJet* jet, Double_t mass) const
{
  Double_t p = jet->P();
  Double_t e = TMath::Sqrt(mass*mass + p*p);

  mom.SetPtEtaPhiE(jet->Pt(), jet->Eta(), jet->Phi(), e);

  return kTRUE;
}

/**
 * Calculates the jet 4-momentum using the default mass hypothesis of the container.
 * If the mass hypothesis is not set, it uses the mass set inside the AliEmcalJet object.
 * @param[out] mom Reference to a TLorentzVector object where the 4-momentum is returned
 * @param[in] jet A pointer to an AliEmcalJet object
 * @return kTRUE if successful, kFALSE if jet is not a valid pointer
 */
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

/**
 * Finds the jet at position i in the container and calculates the jet 4-momentum
 * using the default mass hypothesis of the container.
 * If the mass hypothesis is not set, it uses the mass set inside the AliEmcalJet object.
 * @param[out] mom Reference to a TLorentzVector object where the 4-momentum is returned
 * @param[in] i Index position of the jet in the container
 * @return kTRUE if successful, kFALSE if no jet is found at position i
 */
Bool_t AliJetContainer::GetMomentum(TLorentzVector &mom, Int_t i) const
{
  AliEmcalJet *jet = GetJet(i);
  return GetMomentumFromJet(mom, jet);
}

/**
 * Iterator over jets in the container. Get the next
 * jet in the array. If the end is reached it will return kFALSE.
 * Calculates the 4-momentum using the default mass hypothesis.
 * @deprecated Only for backward compatibility - use AliJetIterableContainer instead
 * @param[out] mom Reference to a TLorentzVector object where the 4-momentum is returned
 * @return kFALSE if the end of the container is not yet reached
 */
Bool_t AliJetContainer::GetNextMomentum(TLorentzVector &mom)
{
  AliEmcalJet *jet = GetNextJet();
  return GetMomentumFromJet(mom, jet);
}

/**
 * Finds the jet at position i in the container, check if it passes the cuts
 * and then calculates the jet 4-momentum
 * using the default mass hypothesis of the container.
 * If the mass hypothesis is not set, it uses the mass set inside the AliEmcalJet object.
 * @param[out] mom Reference to a TLorentzVector object where the 4-momentum is returned
 * @param[in] i Index position of the jet in the container
 * @return kTRUE if successful, kFALSE if no jet is found at position i or the jet is not accepted
 */
Bool_t AliJetContainer::GetAcceptMomentum(TLorentzVector &mom, Int_t i) const
{
  AliEmcalJet *jet = GetAcceptJet(i);
  return GetMomentumFromJet(mom, jet);
}

/**
 * Iterator over jets in the container. Get the next accepted
 * jet in the array. If the end is reached it will return kFALSE.
 * Calculates the 4-momentum using the default mass hypothesis.
 * @deprecated Only for backward compatibility - use AliJetIterableContainer instead
 * @param[out] mom Reference to a TLorentzVector object where the 4-momentum is returned
 * @return kFALSE if the end of the container is not yet reached
 */
Bool_t AliJetContainer::GetNextAcceptMomentum(TLorentzVector &mom)
{
  AliEmcalJet *jet = GetNextAcceptJet();
  return GetMomentumFromJet(mom, jet);
}

/**
 * Checks if a jet passes the cuts.
 * @param[in] jet Pointer to a AliEmcalJet object
 * @param[out] Rejection reason bit in case the jet does not pass the cuts
 * @return kTRUE if jet passes the cuts, kFALSE otherwise
 */
Bool_t AliJetContainer::AcceptJet(const AliEmcalJet *jet, UInt_t &rejectionReason) const
{
  if (fTpcHolePos>0) {
    Bool_t s = CheckTpcHolesOverlap(jet,rejectionReason);
    if (!s) return kFALSE; 
  }
  
  Bool_t r = ApplyJetCuts(jet, rejectionReason);
  if (!r) return kFALSE;

  AliTLorentzVector mom;
  GetMomentumFromJet(mom, jet);

  return ApplyKinematicCuts(mom, rejectionReason);
}

/**
 * Find the jet at position i in the container and checks if it passes the cuts.
 * @param[in] i Index position in the container
 * @param[out] Rejection reason bit in case the jet does not pass the cuts
 * @return kTRUE if jet passes the cuts, kFALSE otherwise
 */
Bool_t AliJetContainer::AcceptJet(Int_t i, UInt_t &rejectionReason) const
{
  if (fTpcHolePos>0) {
    Bool_t s = CheckTpcHolesOverlap(GetJet(i),rejectionReason);
    if (!s) return kFALSE; 
   }
  
  Bool_t r = ApplyJetCuts(GetJet(i), rejectionReason);
  if (!r) return kFALSE;

  AliTLorentzVector mom;
  GetMomentum(mom, i);

  return ApplyKinematicCuts(mom, rejectionReason);
}

/**
 * Apply the jet specific cuts to a jet object
 * @param[in] jet Pointer to a AliEmcalJet object
 * @param[out] Rejection reason bit in case the jet does not pass the cuts
 * @return kTRUE if jet passes the cuts, kFALSE otherwise
 */
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

  if (fMinNConstituents > 0 && jet->GetNumberOfConstituents() < fMinNConstituents) {
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
  
  // If user sets acceptance selection type, compare it to jet's acceptance type, bitwise
  if (fJetAcceptanceType != 0) {
    UInt_t jetAccType = jet->GetJetAcceptanceType();
    UInt_t isAccepted = jetAccType & fJetAcceptanceType;
    if (!isAccepted)
      return kFALSE;
  }

  return kTRUE;
}

/**
 * Retrieve the transverse momentum of leading hadron of the jet.
 * Depending on the internal settings, the leading hadron can be
 * restricted to be charged, neutral or any.
 * @param jet Pointer to a AliEmcalJet object
 * @return The transverse momentum of the leading hadron
 */
Double_t AliJetContainer::GetLeadingHadronPt(const AliEmcalJet *jet) const
{
  if (fLeadingHadronType == 0)       // charged leading hadron
    return jet->MaxTrackPt();
  else if (fLeadingHadronType == 1)  // neutral leading hadron
    return jet->MaxClusterPt();
  else                               // charged or neutral
    return jet->MaxPartPt();
}

/**
 * Retrieve the 4-momentum of leading hadron of the jet.
 * The mass hypothesis is always set to the pion mass (0.139 GeV/c^2).
 * @param[out] mom Reference to a TLorentzVector object where the result is returned
 * @param[in] jet Pointer to a AliEmcalJet object
 */
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

/**
 * Calculates the momentum fraction of the leading calorimeter cluster
 * that belongs to the jet.
 * @param jet A valid pointer to an AliEmcalJet object
 * @return The momentum fraction of the leading calorimeter cluster
 */
Double_t AliJetContainer::GetZLeadingEmc(const AliEmcalJet *jet) const
{
  if (fClusterContainer && fClusterContainer->GetArray()) {
    TLorentzVector mom;

    AliVCluster *cluster = jet->GetLeadingCluster(fClusterContainer->GetArray());
    if (cluster) {
      cluster->GetMomentum(mom, fVertex);

      return GetZ(jet,mom);
    }
    else {
      return -1;
    }
  }
  else {
    return -1;
  }
}

/**
 * Calculates the momentum fraction of the leading track
 * that belongs to the jet.
 * @param jet A valid pointer to an AliEmcalJet object
 * @return The momentum fraction of the leading track
 */
Double_t AliJetContainer::GetZLeadingCharged(const AliEmcalJet *jet) const
{

  if (fParticleContainer && fParticleContainer->GetArray() ) {
    TLorentzVector mom;

    AliVParticle *track = jet->GetLeadingTrack(fParticleContainer->GetArray());
    if (track) {
      mom.SetPtEtaPhiM(track->Pt(),track->Eta(),track->Phi(),0.139);

      return GetZ(jet,mom);
    }
    else {
      return -1;
    }
  }
  else {
    return -1;
  }
}

/**
 * Calculates the momentum fraction carried by a 4-momentum
 * with respect to the jet.
 * @param jet A valid pointer to an AliEmcalJet object
 * @param mom Constant reference to a TLorentz objetc
 * @return The momentum fraction of 4-momentum
 */
Double_t AliJetContainer::GetZ(const AliEmcalJet *jet, const TLorentzVector& mom) const
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

/**
 * Prints the current cuts to the standard output, for debug purposes.
 */
void AliJetContainer::PrintCuts() 
{
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

/**
 * Resets the cuts to the default values.
 */
void AliJetContainer::ResetCuts() 
{
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

/**
 * Calculates the total number of accepted jets.
 * @return Number of accepted jets.
 */
Int_t AliJetContainer::GetNAcceptedJets()
{
  return accepted().GetEntries();
}

/**
 * Get fraction of shared pT between matched jets.
 * Uses ClosestJet() jet pT as baseline: fraction = \Sum_{const,jet1} pT,const,i / pT,jet,closest
 * If no particle container is given, it assumes that the matched jet's constituent
 * come from the same particle container.
 * @param jet1 Pointer to an AliEmcalJet object
 * @param cont2 Pointer to an AliParticleContainer object in which the constituents of the matched jet are to be found
 * @return The momentum fraction of jet1 that is shared by its closest matched jet
 */
Double_t AliJetContainer::GetFractionSharedPt(const AliEmcalJet *jet1, AliParticleContainer *cont2) const
{
  AliEmcalJet *jet2 = jet1->ClosestJet();
  if (!jet2) return -1;

  Double_t jetPt2 = jet2->Pt();
  if (jetPt2 <= 0) return -1;

  Int_t bgeom = kTRUE;
  if (!cont2) bgeom = kFALSE;
  Double_t sumPt = 0.;
  AliVParticle *vpf = 0x0;
  Int_t iFound = 0;
  for (Int_t icc = 0; icc < jet2->GetNumberOfTracks(); icc++) {
    Int_t idx = (Int_t)jet2->TrackAt(icc);
    //get particle
    AliVParticle *p2 = 0x0;
    if (bgeom) p2 = static_cast<AliVParticle*>(jet2->TrackAt(icc, cont2->GetArray()));
    iFound = 0;
    for (Int_t icf = 0; icf < jet1->GetNumberOfTracks(); icf++) {
      if (!bgeom && idx == jet1->TrackAt(icf) && iFound == 0 ) {
        iFound = 1;
        vpf = jet1->Track(icf);
        if (vpf) sumPt += vpf->Pt();
        continue;
      }
      if (bgeom){
        vpf = jet1->Track(icf);
        if (!vpf) continue;
        if (!SamePart(vpf, p2, 1.e-4)) continue; //not the same particle
        sumPt += vpf->Pt();
      }
    }
  }

  Double_t fraction = sumPt / jetPt2;

  return fraction;
}

/**
 * Generate the jet branch name according to a given jet definition.
 * @param jetType Type of the jet (full, charged, neutral)
 * @param jetAlgo Jet algorithm used to reconstruct the jets (anti-kt, kt, etc.)
 * @param recoScheme Jet recombination scheme (E-scheme, pt-scheme)
 * @param radius Jet resolution parameter
 * @param partCont Particle container used to feed the jet finder
 * @param clustCont Cluster container used to feed the jet finder
 * @param tag Additional tag for the jet branch name (default is "Jet")
 * @return A string containing the jet branch name
 */
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

/**
 * Create an iterable container interface over all objects in the
 * EMCAL container.
 * @return iterable container over all objects in the EMCAL container
 */
const AliJetIterableContainer AliJetContainer::all() const {
  return AliJetIterableContainer(this, false);
}

/**
 * Create an iterable container interface over accepted objects in the
 * EMCAL container.
 * @return iterable container over accepted objects in the EMCAL container
 */
const AliJetIterableContainer AliJetContainer::accepted() const {
  return AliJetIterableContainer(this, true);
}

/**
 * Create an iterable container interface over all objects in the
 * EMCAL container.
 * @return iterable container over all objects in the EMCAL container
 */
const AliJetIterableMomentumContainer AliJetContainer::all_momentum() const {
  return AliJetIterableMomentumContainer(this, false);
}

/**
 * Create an iterable container interface over accepted objects in the
 * EMCAL container.
 * @return iterable container over accepted objects in the EMCAL container
 */
const AliJetIterableMomentumContainer AliJetContainer::accepted_momentum() const {
  return AliJetIterableMomentumContainer(this, true);
}

/**
 * Generates a title for this container. The title is generated putting together the jet branch
 * name and the pT cut.
 * @return A pointer to a statically allocated string which contains the title.
 */
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

/**
 * Checks that the jet is far enough (axis + radius) from a hole in the TPC acceptance.
 * @param[in] jet A pointer to an AliEmcalJet object
 * @param[out] rejectionReason If jet does not pass the cut the corresponding rejection reason bit is set.
 */
Bool_t AliJetContainer::CheckTpcHolesOverlap(const AliEmcalJet *jet, UInt_t &rejectionReason) const
{   
  if (!jet) {
    AliDebug(11,"No jet found");
    rejectionReason |= kNullObject;
    return kFALSE;
  } 

  Double_t disthole = RelativePhi(jet->Phi(), fTpcHolePos);
  if (TMath::Abs(disthole) < (fTpcHoleWidth + fJetRadius)){
    AliDebug(11,"Jet overlaps with TPC hole");
    rejectionReason |= kOverlapTpcHole;
    return kFALSE;
  }

  return kTRUE;
}
