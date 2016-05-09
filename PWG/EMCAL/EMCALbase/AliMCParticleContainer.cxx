//
// Container with name, TClonesArray and cuts for particles
//
// Author: M. Verweij, S. Aiola

#include <TClonesArray.h>

#include "AliVEvent.h"
#include "AliLog.h"

#include "AliTLorentzVector.h"
#include "AliMCParticleContainer.h"

/// \cond CLASSIMP
ClassImp(AliMCParticleContainer);
/// \endcond

/**
 * Default constructor.
 */
AliMCParticleContainer::AliMCParticleContainer():
  AliParticleContainer(),
  fMCFlag(AliAODMCParticle::kPhysicalPrim)
{
  fBaseClassName = "AliAODMCParticle";
  SetClassName("AliAODMCParticle");
}

/**
 * Standard constructor.
 * @param[in] name Name of the container (= name of the array operated on)
 */
AliMCParticleContainer::AliMCParticleContainer(const char *name):
  AliParticleContainer(name),
  fMCFlag(AliAODMCParticle::kPhysicalPrim)
{
  fBaseClassName = "AliAODMCParticle";
  SetClassName("AliAODMCParticle");
}

/**
 * Get MC particle using the MC label
 * @param[in] lab Label of the particle
 * @return pointer to particle if particle is found, NULL otherwise
 */
AliAODMCParticle* AliMCParticleContainer::GetMCParticleWithLabel(Int_t lab) const
{
  Int_t i = GetIndexFromLabel(lab);
  return GetMCParticle(i);
}

/**
 * Get MC particle using the MC label
 * @param[in] lab Label of the particle
 * @return pointer to particle if particle is found and accepted, NULL otherwise
 */
AliAODMCParticle* AliMCParticleContainer::GetAcceptMCParticleWithLabel(Int_t lab)
{
  Int_t i = GetIndexFromLabel(lab);
  return GetAcceptMCParticle(i);
}

/**
 * Get track at index in the container
 * @param[in] i Index of the particle in the container
 * @return pointer to particle if particle is found, NULL otherwise
 */
AliAODMCParticle* AliMCParticleContainer::GetMCParticle(Int_t i) const
{
  if (i == -1) i = fCurrentID;
  if (i < 0 || i >= fClArray->GetEntriesFast()) return 0;
  AliAODMCParticle *vp = static_cast<AliAODMCParticle*>(fClArray->At(i));
  return vp;
}

/**
 * Get track at index in the container
 * @param[in] i Index of the particle in the container
 * @return pointer to particle if particle is accepted, NULL otherwise
 */
AliAODMCParticle* AliMCParticleContainer::GetAcceptMCParticle(Int_t i) const
{
  //return pointer to particle if particle is accepted

  UInt_t rejectionReason = 0;
  if (i == -1) i = fCurrentID;
  if (AcceptMCParticle(i, rejectionReason)) {
      return GetMCParticle(i);
  }
  else {
    AliDebug(2,"Particle not accepted.");
    return 0;
  }
}

/**
 * Get next accepted particle in the container selected using the track cuts provided.
 * @deprecated Old style iterator - for compatibility reasons, use AliParticleContainer::accept_iterator instead
 * @return Next accepted particle (NULL if the end of the array is reached)
 */
AliAODMCParticle* AliMCParticleContainer::GetNextAcceptMCParticle()
{
  const Int_t n = GetNEntries();
  AliAODMCParticle *p = 0;
  do {
    fCurrentID++;
    if (fCurrentID >= n) break;
    p = GetAcceptMCParticle(fCurrentID);
  } while (!p);

  return p;
}

/**
 * Get next particle in the container
 * @deprecated Old style iterator - for compatibility reasons, use AliParticleContainer::all_iterator instead
 * @return Next track in the container (NULL if end of the container is reached)
 */
AliAODMCParticle* AliMCParticleContainer::GetNextMCParticle()
{
  //Get next particle

  const Int_t n = GetNEntries();
  AliAODMCParticle *p = 0;
  do {
    fCurrentID++;
    if (fCurrentID >= n) break;
    p = GetMCParticle(fCurrentID);
  } while (!p);

  return p;
}

/**
 * Perform full MC particle selection for the particle vp, consisting
 * of kinematical particle selection and MC-specific cuts
 * @param[in] vp Particle to be checked
 * @param[in] rejectionReason Bitmap encoding the reason why the
 * particle was rejected. Note: The variable is not set to NULL
 * inside this function before changing its value.
 * @return True if the particle is accepted, false otherwise
 */
Bool_t AliMCParticleContainer::AcceptMCParticle(const AliAODMCParticle *vp, UInt_t &rejectionReason) const
{
  // Return true if vp is accepted.
  Bool_t r = ApplyMCParticleCuts(vp, rejectionReason);
  if (!r) return kFALSE;

  AliTLorentzVector mom;
  GetMomentumFromParticle(mom, vp);

  return ApplyKinematicCuts(mom, rejectionReason);
}

/**
 * Perform full MC particle selection for the particle vp, consisting
 * of kinematical particle selection and MC-specific cuts
 * @param[in] i Index of the particle to check
 * @param[in] rejectionReason Bitmap encoding the reason why the
 * particle was rejected. Note: The variable is not set to NULL
 * inside this function before changing its value.
 * @return True if the particle is accepted, false otherwise
 */
Bool_t AliMCParticleContainer::AcceptMCParticle(Int_t i, UInt_t &rejectionReason) const
{
  // Return true if vp is accepted.
  Bool_t r = ApplyMCParticleCuts(GetMCParticle(i), rejectionReason);
  if (!r) return kFALSE;

  AliTLorentzVector mom;
  GetMomentum(mom, i);

  return ApplyKinematicCuts(mom, rejectionReason);
}

/**
 * Apply MC particle cuts, e.g. primary particle selection.
 * @param[in] vp Particle to be checked
 * @param[in] rejectionReason Bitmap encoding the reason why the
 * particle was rejected. Note: The variable is not set to NULL
 * inside this function before changing its value.
 * @return True if the particle is accepted, false otherwise
 */
Bool_t AliMCParticleContainer::ApplyMCParticleCuts(const AliAODMCParticle* vp, UInt_t &rejectionReason) const
{
  // Return true if i^th particle is accepted.
  // Cuts on the particle properties

  if ((vp->GetFlag() & fMCFlag) != fMCFlag) {
    rejectionReason |= kMCFlag;
    return kFALSE;
  }

  return ApplyParticleCuts(vp, rejectionReason);
}

/**
 * Create an iterable container interface over all objects in the
 * EMCAL container.
 * @return iterable container over all objects in the EMCAL container
 */
const AliMCParticleIterableContainer AliMCParticleContainer::all() const {
  return AliMCParticleIterableContainer(this, false);
}

/**
 * Create an iterable container interface over accepted objects in the
 * EMCAL container.
 * @return iterable container over accepted objects in the EMCAL container
 */
const AliMCParticleIterableContainer AliMCParticleContainer::accepted() const {
  return AliMCParticleIterableContainer(this, true);
}

/**
 * Build title of the container consisting of the container name
 * and a string encoding the minimum \f$ p_{t} \f$ cut applied
 * in the kinematic particle selection.
 * @return Title of the container
 */
const char* AliMCParticleContainer::GetTitle() const
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
