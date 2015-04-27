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
#include "AliReducedGeneratedParticle.h"

/// \cond CLASSIMP
ClassImp(HighPtTracks::AliReducedGeneratedParticle)
/// \endcond

namespace HighPtTracks {

/**
 * Default (dummy) constructor
 */
AliReducedGeneratedParticle::AliReducedGeneratedParticle():
  TObject(),
  fParticleID(-1),
  fPDGCode(-10000),
  fEnergy(-1.)
{
  memset(fPVec, 0, sizeof(Double_t) * 3);
}

/**
 * Constructor, setting parameters of the particle
 * @param id Unique ID of the particle
 * @param pdg PDG code
 * \param px x-component of the momentum vector
 * \param py y-component of the momentum vector
 * \param pz z-component of the momentum vector
 * @param energy Particle energy
 */
AliReducedGeneratedParticle::AliReducedGeneratedParticle(Int_t id, Int_t pdg,
    Double_t px, Double_t py, Double_t pz, Double_t energy):
  TObject(),
  fParticleID(id),
  fPDGCode(pdg),
  fEnergy(energy)
{
  fPVec[0] = px;
  fPVec[1] = py;
  fPVec[2] = px;
}

/**
 * Destructor, nothing to do
 */
AliReducedGeneratedParticle::~AliReducedGeneratedParticle() {}

} /* namespace HighPtTracks */
