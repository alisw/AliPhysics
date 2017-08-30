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
#include <algorithm>
#include <iostream>
#include <vector>
#include <TClonesArray.h>

#include "AliVEvent.h"
#include "AliLog.h"

#include "AliTLorentzVector.h"
#include "AliParticleContainer.h"

/// \cond CLASSIMP
ClassImp(AliParticleContainer);
/// \endcond

// Properly instantiate the object
AliEmcalContainerIndexMap <TClonesArray, AliVParticle> AliParticleContainer::fgEmcalContainerIndexMap;

/**
 * Default constructor.
 */
AliParticleContainer::AliParticleContainer():
  AliEmcalContainer(),
  fMinDistanceTPCSectorEdge(-1),
  fChargeCut(kNoChargeCut),
  fGeneratorIndex(-1)
{
  fBaseClassName = "AliVParticle";
  SetClassName("AliVParticle");
}

/**
 * Standard constructor.
 * @param name Name of the particle branch (TClonesArray)
 */
AliParticleContainer::AliParticleContainer(const char *name) :
  AliEmcalContainer(name),
  fMinDistanceTPCSectorEdge(-1),
  fChargeCut(kNoChargeCut),
  fGeneratorIndex(-1)
{
  fBaseClassName = "AliVParticle";
  SetClassName("AliVParticle");
}

/**
 * Get the leading particle in the container. If "p" is contained in the parameter opt,
 * then the absolute momentum is use instead of the transverse momentum.
 * @param[in] opt Options for the selection of the leading particle
 * @return Leading particle in the container
 */
AliVParticle* AliParticleContainer::GetLeadingParticle(const char* opt) 
{
  TString option(opt);
  option.ToLower();

  Int_t tempID = fCurrentID;
  ResetCurrentID();

  AliVParticle *partMax = GetNextAcceptParticle();
  AliVParticle *part = 0;

  if (option.Contains("p")) {
    while ((part = GetNextAcceptParticle())) {
      if (part->P() > partMax->P()) partMax = part;
    }
  }
  else {
    while ((part = GetNextAcceptParticle())) {
      if (part->Pt() > partMax->Pt()) partMax = part;
    }
  }

  fCurrentID = tempID;

  return partMax;
}

/**
 * Get \f$ i^{th} \f$ particle in the container.
 * @param[in] i Index of the particle to access
 * @return Parrticle at the given index (NULL if the index is out of range)
 */
AliVParticle* AliParticleContainer::GetParticle(Int_t i) const 
{
  //
  if (i == -1) i = fCurrentID;
  if (i < 0 || i >= this->fClArray->GetEntriesFast()) return 0;
  AliVParticle *vp = static_cast<AliVParticle*>(fClArray->At(i));
  return vp;
}

/**
 * Get \f$ i^{th} \f$ particle in the container if it is accepted.
 * In case it is not accepted a nullpointer is returned.
 * @param[in] i Index of the particle
 * @return Particle at the index if it is accepted, NULL otherwise
 */
AliVParticle* AliParticleContainer::GetAcceptParticle(Int_t i) const
{
  UInt_t rejectionReason = 0;
  if (i == -1) i = fCurrentID;
  if (AcceptParticle(i, rejectionReason)) {
      return GetParticle(i);
  }
  else {
    AliDebug(2,"Particle not accepted.");
    return 0;
  }
}

/**
 * Iterator over accepted particles in the container. Get the next accepted
 * particle in the array. If the end is reached, NULL is returned.
 * @deprecated Only for backward compatibility - use accept_iterator instead
 * @return Next accepted particle in the array (NULL if the end is reached)
 */
AliVParticle* AliParticleContainer::GetNextAcceptParticle()
{
  const Int_t n = GetNEntries();
  AliVParticle *p = 0;
  do {
    fCurrentID++;
    if (fCurrentID >= n) break;
    p = GetAcceptParticle(fCurrentID);
  } while (!p);

  return p;
}

/**
 * Iterator over all particles in the container. Get the next particle in
 * the array. If the end is reached, NULL is returned.
 * @deprecated Only for backward compatibility - use all_iterator instead.
 * @return Next particle in the array (NULL if the end is reached)
 */
AliVParticle* AliParticleContainer::GetNextParticle()
{
  const Int_t n = GetNEntries();
  AliVParticle *p = 0;
  do {
    fCurrentID++;
    if (fCurrentID >= n) break;
    p = GetParticle(fCurrentID);
  } while (!p);

  return p;
}

/**
 * Retrieve momentum information of a particle (part) and fill a TLorentzVector
 * with it. In case the optional parameter mass is provided, it is used as mass
 * hypothesis, otherwise the mass hypothesis from the particle itself is used.
 * @param[out] mom Momentum vector to be filled
 * @param[in] part Particle from which the momentum information is obtained.
 * @param[in] mass (Optional) Mass hypothesis
 * @return
 */
Bool_t AliParticleContainer::GetMomentumFromParticle(TLorentzVector &mom, const AliVParticle* part, Double_t mass) const
{
  if (part && part->Eta() < 1e6 && part->Eta() > -1e6) { // protection against FPE in sinh(eta)
    if (mass < 0) mass = part->M();
    mom.SetPtEtaPhiM(part->Pt(), part->Eta(), part->Phi(), mass);
    return kTRUE;
  }
  else {
    mom.SetPtEtaPhiM(0, 0, 0, 0);
    return kFALSE;
  }
}

/**
 * Fills a TLorentzVector with the momentum information of the particle provided
 * under a global mass hypothesis.
 * @param[out] mom Momentum vector of the particle provided
 * @param[in] part Particle from which to obtain the momentum information
 * @return Always true
 */
Bool_t AliParticleContainer::GetMomentumFromParticle(TLorentzVector &mom, const AliVParticle* part) const
{
  return GetMomentumFromParticle(mom,part,fMassHypothesis);
}

/**
 * Fills a TLorentzVector with the monentum infomation of the
 * \f$ i^{th} \f$ particle in the container, using a global
 * mass hypothesis. In case the provided index is out of
 * range, false is returned as return value.
 * @param[out] mom Momentum vector of the \f$ i^{th} \f$ particle in the array
 * @param[in] i Index of th particle to check
 * @return True if the request was successfull, false otherwise
 */
Bool_t AliParticleContainer::GetMomentum(TLorentzVector &mom, Int_t i) const
{
  if (i == -1) i = fCurrentID;
  AliVParticle *vp = GetParticle(i);
  return GetMomentumFromParticle(mom, vp);
}

/**
 * Fills a TLorentzVector with the monentum infomation of the
 * next particle in the container, using a global mass hypothesis.
 * In case the iterator reached the end of the array, false
 * is returned as return value.
 * @deprecated Old style iterator - use all_iterator instead
 * @param[out] mom Momentum vector of the next particle
 * @return True if the request was successfull, false otherwise
 */
Bool_t AliParticleContainer::GetNextMomentum(TLorentzVector &mom)
{
  AliVParticle *vp = GetNextParticle();
  return GetMomentumFromParticle(mom, vp);
}

/**
 * Fills a TLorentzVector with the monentum infomation of the
 * \f$ i^{th} \f$ accepted particle in the container, using a
 * global mass hypothesis. In case the provided index is out of
 * range, or the particle under the index is not accepted, false
 * is returned as return value.
 * @param[out] mom Momentum vector of the accepted particle
 * @param[in] i Index to check
 * @return True if the request was successfull, false otherwise
 */
Bool_t AliParticleContainer::GetAcceptMomentum(TLorentzVector &mom, Int_t i) const
{
  if (i == -1) i = fCurrentID;
  AliVParticle *vp = GetAcceptParticle(i);
  return GetMomentumFromParticle(mom, vp);
}

/**
 * Fills a TLorentzVector with the monentum infomation of the
 * next accepted particle in the container, using a global
 * mass hypothesis. In case the iteration reached the end of
 * the array, false is returned as return value.
 * @deprecated Old style iterator - use accept_iterator instead
 * @param[out] mom Momentum vector of the next particle in the array
 * @return True if the request was successfull, false (no more entries) otherwise
 */
Bool_t AliParticleContainer::GetNextAcceptMomentum(TLorentzVector &mom)
{
  AliVParticle *vp = GetNextAcceptParticle();
  return GetMomentumFromParticle(mom, vp);
}

/**
 * Perform full particle selection consisting of kinematical and particle property
 * selection on the particle provided. In case the particle is rejected, the
 * reason for the rejection is encoded in the bitmap rejectionReason.
 * @param vp Particle for which to perform the selection
 * @param rejectionReason Bitmap with the reason why the particle was accepted.
 * Note: The value is not set to 0 in the function in order to combine the information
 * with other selection steps.
 * @return True if the particle is accepted, false otherwise.
 */
Bool_t AliParticleContainer::AcceptParticle(const AliVParticle *vp, UInt_t &rejectionReason) const
{
  Bool_t r = ApplyParticleCuts(vp, rejectionReason);
  if (!r) return kFALSE;

  AliTLorentzVector mom;

  Int_t id = fClArray->IndexOf(vp);
  if (id >= 0) {
    GetMomentum(mom, id);
  }
  else {
    GetMomentumFromParticle(mom, vp);
  }

  return ApplyKinematicCuts(mom, rejectionReason);
}

/**
 * Perform full particle selection consisting of kinematical and particle property
 * selection on the \f$ i^{th} \f$ particle in the array. In case the particle is
 * rejected, the reason for the rejection is encoded in the bitmap rejectionReason.
 * @param[in] i Index of the particle to select.
 * @param[out] rejectionReason Bitmap with the reason why the particle was accepted.
 * Note: The value is not set to 0 in the function in order to combine the information
 * with other selection steps.
 * @return True if the particle was accepted, false otherwise
 */
Bool_t AliParticleContainer::AcceptParticle(Int_t i, UInt_t &rejectionReason) const
{
  Bool_t r = ApplyParticleCuts(GetParticle(i), rejectionReason);
  if (!r) return kFALSE;

  AliTLorentzVector mom;
  GetMomentum(mom, i);

  return ApplyKinematicCuts(mom, rejectionReason);
}

/**
 * Apply cuts on the particle properties. Implemented are
 * - Monte-Carlo label
 * - Charge
 * - Generator index
 *
 * @param[in] vp Particle for which the selection is performed
 * @param[out] rejectionReason Bitmap with the reason why the particle was accepted.
 * Note: The value is not set to 0 in the function in order to combine the information
 * with other selection steps.
 * @return True if the particle is accepted, False otherwise
 */
Bool_t AliParticleContainer::ApplyParticleCuts(const AliVParticle* vp, UInt_t &rejectionReason) const
{
  if (!vp) {
    rejectionReason |= kNullObject;
    return kFALSE;
  }

  if (vp->TestBits(fBitMap) != (Int_t)fBitMap) {
    rejectionReason |= kBitMapCut;
    return kFALSE;
  }

  if (fMinMCLabel >= 0 && TMath::Abs(vp->GetLabel()) < fMinMCLabel) {
    rejectionReason |= kMCLabelCut;
    return kFALSE;
  }

  if (fMaxMCLabel >= 0 && TMath::Abs(vp->GetLabel()) > fMaxMCLabel) {
    rejectionReason |= kMCLabelCut;
    return kFALSE;
  }

  switch (fChargeCut) {
  case kCharged:
    if (vp->Charge() == 0) {
      rejectionReason |= kChargeCut;
      return kFALSE;
    }
    break;

  case kNeutral:
    if (vp->Charge() != 0) {
      rejectionReason |= kChargeCut;
      return kFALSE;
    }
    break;

  case kPositiveCharge:
    if (vp->Charge() <= 0) {
      rejectionReason |= kChargeCut;
      return kFALSE;
    }
    break;

  case kNegativeCharge:
    if (vp->Charge() >= 0) {
      rejectionReason |= kChargeCut;
      return kFALSE;
    }
    break;

  default:
    break;
  }

  if (fGeneratorIndex >= 0 && fGeneratorIndex != vp->GetGeneratorIndex()) {
    rejectionReason |= kMCGeneratorCut;
    return kFALSE;
  }

  return kTRUE;
}

/**
 * Apply kinematical cuts to the momentum vector provided. In addition
 * to the standard kinematical cuts in \f$ p_{t} \f$, \f$ \eta \f$ and
 * \f$ \phi \f$, implemented in the AliEmcalContainer::ApplyKinematicCuts,
 * also the distance to the TPC sector boundary is checked.
 * @param[in] mom Momentum vector for which
 * @param[out] rejectionReason  Bitmap with the reason why the particle was accepted.
 * Note: The value is not set to 0 in the function in order to combine the information
 * with other selection steps.
 * @return True if the momentum vector was selected under the given cuts, false otherwise
 */
Bool_t AliParticleContainer::ApplyKinematicCuts(const AliTLorentzVector& mom, UInt_t &rejectionReason) const
{
  if(fMinDistanceTPCSectorEdge>0.) {
    const Double_t pi = TMath::Pi();
    const Double_t kSector = pi/9;
    Double_t phiDist = TMath::Abs(mom.Phi() - TMath::FloorNint(mom.Phi()/kSector)*kSector);
    if(phiDist<fMinDistanceTPCSectorEdge) {
      rejectionReason |= kMinDistanceTPCSectorEdgeCut;
      return kFALSE;
    }
  }

  return AliEmcalContainer::ApplyKinematicCuts(mom, rejectionReason);
}

/**
 * Get number of accepted particles. In order to get this number,
 * the selection has to be applied to each particle within this
 * function.
 * @return Number of selected particles under the given particle selection
 */
Int_t AliParticleContainer::GetNAcceptedParticles() const
{
  Int_t nPart = 0;
  for(int ipart = 0; ipart < this->GetNParticles(); ipart++){
    UInt_t rejectionReason = 0;
    if(this->AcceptParticle(ipart, rejectionReason)) nPart++;
  }
  return nPart;
}

/**
 * Make a title of the container name based on the min \f$ p_{t} \f$ used
 * in the particle selection process.
 * @return Title of the container
 */
const char* AliParticleContainer::GetTitle() const
{
  static TString trackString;

  if (GetMinPt() == 0) {
    trackString = TString::Format("%s_pT0000", GetArrayName().Data());
  }
  else if (GetMinPt() < 1.0) {
    trackString = TString::Format("%s_pT0%3.0f", GetArrayName().Data(), GetMinPt()*1000.0);
  }
  else {
    trackString = TString::Format("%s_pT%4.0f", GetArrayName().Data(), GetMinPt()*1000.0);
  }

  return trackString.Data();
}

/**
 * Connect the container to the array with content stored inside the virtual event.
 * The object name in the event must match the name given in the constructor.
 *
 * Additionally register the array into the index map.
 *
 * @param event Input event containing the array with content.
 */
void AliParticleContainer::SetArray(const AliVEvent * event)
{
  AliEmcalContainer::SetArray(event);

  // Register TClonesArray in index map
  fgEmcalContainerIndexMap.RegisterArray(GetArray());
}

/**
 * Create an iterable container interface over all objects in the
 * EMCAL container.
 * @return iterable container over all objects in the EMCAL container
 */
const AliParticleIterableContainer AliParticleContainer::all() const {
  return AliParticleIterableContainer(this, false);
}

/**
 * Create an iterable container interface over accepted objects in the
 * EMCAL container.
 * @return iterable container over accepted objects in the EMCAL container
 */
const AliParticleIterableContainer AliParticleContainer::accepted() const {
  return AliParticleIterableContainer(this, true);
}

/**
 * Create an iterable container interface over all objects in the
 * EMCAL container.
 * @return iterable container over all objects in the EMCAL container
 */
const AliParticleIterableMomentumContainer AliParticleContainer::all_momentum() const {
  return AliParticleIterableMomentumContainer(this, false);
}

/**
 * Create an iterable container interface over accepted objects in the
 * EMCAL container.
 * @return iterable container over accepted objects in the EMCAL container
 */
const AliParticleIterableMomentumContainer AliParticleContainer::accepted_momentum() const {
  return AliParticleIterableMomentumContainer(this, true);
}

/******************************************
 * Unit tests                             *
 ******************************************/

int TestParticleContainerIterator(const AliParticleContainer *const cont, int iteratorType, bool verbose){
  std::vector<AliVParticle *> reference, variation;
  AliVParticle *test = NULL;
  for(int iclust = 0; iclust < cont->GetNParticles(); iclust++){
    test = cont->GetParticle(iclust);
    if(!iteratorType){
      UInt_t rejectionReason = 0;
      if(!cont->AcceptParticle(test, rejectionReason)) continue;
    }
    reference.push_back(test);
  }

  if(!iteratorType){
    // test accept iterator
    for(auto part : cont->accepted()){
      variation.push_back(part);
    }
  } else {
    // test all iterator
    for(auto part : cont->all()){
      variation.push_back(part);
    }
  }

  if(verbose) {
    // Printing cluster addresses:
    if(reference.size() < 30){
      std::cout << "Paritcles in reference container: " << std::endl;
      std::cout << "===========================================" << std::endl;
      for(std::vector<AliVParticle *>::iterator refit = reference.begin(); refit != reference.end(); ++refit){
        std::cout << "Address: " << *refit << std::endl;
      }
    }
    if(variation.size() < 30){
      std::cout << "Paritcles in test container: " << std::endl;
      std::cout << "===========================================" << std::endl;
      for(std::vector<AliVParticle *>::iterator varit = variation.begin(); varit != variation.end(); ++varit){
        std::cout << "Address: " << *varit << std::endl;
      }
    }
  }

  int testresult = 0;
  // compare distributions: all particles in one vector needs to be found in the other vector and vice versa
  bool failure = false;
  // first test: cleck if all particles are found by the iterator
  for(std::vector<AliVParticle *>::iterator clit = reference.begin(); clit != reference.end(); ++clit){
    if(std::find(variation.begin(), variation.end(), *clit) == variation.end()) {
      if(verbose)
        std::cout << "Could not find particle with address " << *clit << " in test container" << std::endl;
      failure = true;
      break;
    }
  }
  if(!failure){
    // second test: check if there are no extra particles found
    for(std::vector<AliVParticle *>::iterator clit = variation.begin(); clit != variation.end(); ++clit){
      if(std::find(reference.begin(), reference.end(), *clit) == reference.end()) {
       if(verbose)
          std::cout << "Could not find particle with address " << *clit << " in reference container" << std::endl;
        failure = true;
        break;
      }
    }
    if(failure) testresult = 2;
  } else testresult = 1;
  if(verbose) {
    std::cout << "Unit test particle container, iterator type " << iteratorType << std::endl;
    std::cout << "Number of expected particles: " << reference.size() << ", number of found particles: " << variation.size() << std::endl;
    std::cout << "Test result: " << testresult << std::endl;
    std::cout << "-----------------------------------------------------------------------" << std::endl;
  }
  return testresult;
}
