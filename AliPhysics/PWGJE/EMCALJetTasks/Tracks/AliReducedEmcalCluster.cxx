/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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
#include <TArrayD.h>
#include <TLorentzVector.h>
#include <TObjArray.h>

#include "AliReducedEmcalCluster.h"

/// \cond CLASSIMP
ClassImp(HighPtTracks::AliReducedClusterParticle)
ClassImp(HighPtTracks::AliReducedEmcalCluster)
/// \endcond

namespace HighPtTracks {

/**
 * Constructor, initialising with default values
 */
AliReducedEmcalCluster::AliReducedEmcalCluster() :
  fClusterID(-1),
  fEnergy(-1),
  fEta(-100),
  fPhi(-100),
  fM02(-1),
  fM20(-1),
  fContributors(NULL)
{
  memset(fCellEnergies, 0, sizeof(Float_t) * 3);
}

/**
 * Constructor, initialising cluster with all parameters
 * \param id ID of the cluster
 * \param energy Energy of the cluster
 * \param eta \f$ \eta \f$ position of the cluster
 * \param phi \f$ \phi \f$ position of the cluster
 * \param m02 The m02 shower shape parameter
 * \param m20 The m20 shower shape parameter
 */
AliReducedEmcalCluster::AliReducedEmcalCluster(Int_t id, Float_t energy,Float_t eta, Float_t phi, Float_t m02, Float_t m20):
  TObject(),
  fClusterID(id),
  fEnergy(energy),
  fEta(eta),
  fPhi(phi),
  fM02(m02),
  fM20(m20),
  fContributors(NULL)
{
  memset(fCellEnergies, 0, sizeof(Float_t)*3);
}

/**
 * Copy constructor
 * \param ref Reference for the copy
 */
AliReducedEmcalCluster::AliReducedEmcalCluster(const AliReducedEmcalCluster& ref) :
  TObject(ref),
  fClusterID(ref.fClusterID),
  fEnergy(ref.fEnergy),
  fEta(ref.fEta),
  fPhi(ref.fPhi),
  fM02(ref.fM02),
  fM20(ref.fM20),
  fContributors(NULL)
{
  ref.Copy(*this);
}

/**
 * Assignment operator
 * \param ref Reference for the copy
 * \return This object
 */
AliReducedEmcalCluster& AliReducedEmcalCluster::operator=(const AliReducedEmcalCluster& ref) {
  if(&ref != this){
    this->~AliReducedEmcalCluster();
    TObject::operator=(ref);
    ref.Copy(*this);
  }
  return *this;
}

/**
 * Destructor
 */
AliReducedEmcalCluster::~AliReducedEmcalCluster() {
  if(fContributors) delete fContributors;
}

/**
 * Copy content of this cluster into target
 * \param target Target for the copy
 */
void AliReducedEmcalCluster::Copy(TObject &target) const {
  AliReducedEmcalCluster *targetcluster = dynamic_cast<AliReducedEmcalCluster *>(&target);
  if(!targetcluster) return;
  targetcluster->fClusterID = fClusterID;
  targetcluster->fEnergy = fEnergy;
  targetcluster->fEta = fEta;
  targetcluster->fPhi = fPhi;
  targetcluster->fM02 = fM02;
  targetcluster->fM20 = fM20;
  targetcluster->fContributors = NULL;
  memcpy(targetcluster->fCellEnergies, fCellEnergies, sizeof(Float_t)*3);
  if(fContributors){
    targetcluster->fContributors = new TObjArray;
    targetcluster->fContributors->SetOwner(kTRUE);
    for(TIter contiter = TIter(fContributors).Begin(); contiter != TIter::End(); ++contiter){
      targetcluster->fContributors->Add(new AliReducedClusterParticle(*(static_cast<AliReducedClusterParticle *>(*contiter))));
    }
  }
}

/**
 * Fill cell energies into target array
 * \param target
 */
void AliReducedEmcalCluster::FillCellEnergies(TArrayD& target) {
  target.Set(3);
  for(Int_t icell = 0; icell < 3; icell++) target[icell] = fCellEnergies[icell];
}

/**
 * Add new true cluster contributor to the cluster
 * \param pdg PDG code of the true particle
 * \param px x-component of the momentum vector of the true particle
 * \param py y-component of the momentum vector of the true particle
 * \param pz z-component of the momentum vector of the true particle
 * \param energy Energy of the true particle
 */
void AliReducedEmcalCluster::AddTrueContributor(Int_t pdg, Double_t px, Double_t py, Double_t pz, Double_t energy) {
  if(!fContributors){
    fContributors = new TObjArray;
    fContributors->SetOwner(kTRUE);
  }
  fContributors->Add(new AliReducedClusterParticle(pdg, px, py, pz, energy));
}

/**
 * Dummy constructor
 */
HighPtTracks::AliReducedClusterParticle::AliReducedClusterParticle():
  fPdg(0),
  fEnergy(0)
{
  memset(fPvec, 0, sizeof(Float_t) * 3);
}

/**
 * Default constructor, initializing params
 * \param pdg Particle PDG code
 * \param px x-component of the momentum vector
 * \param py y-component of the momentum vector
 * \param pz z-component of the momentum vector
 * \param energy Particle energy
 */
HighPtTracks::AliReducedClusterParticle::AliReducedClusterParticle(Int_t pdg,
    Double_t px, Double_t py, Double_t pz, Double_t energy) :
  fPdg(pdg),
  fEnergy(energy)
{
  fPvec[0] = px;
  fPvec[1] = py;
  fPvec[2] = pz;
}

/**
 * Destructor, does nothing
 */
HighPtTracks::AliReducedClusterParticle::~AliReducedClusterParticle() { }

/**
 * Accessor to particle kinematical information via a TLorentzVector
 * \param target Lorentzvector to fill with particle kinematical information
 */
void HighPtTracks::AliReducedClusterParticle::FillLorentzVector(TLorentzVector& target) const {
  target.SetPxPyPzE(fPvec[0],fPvec[1],fPvec[2],fEnergy);
}

} /* namespace HighPtTracks */

