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
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TParticlePDG.h>

#include "AliReducedJetConstituent.h"

/// \cond CLASSIMP
ClassImp(HighPtTracks::AliReducedJetConstituent)
/// \endcond

namespace HighPtTracks {

/**
 * Dummy (I/O) constructor, not to be used.
 */
AliReducedJetConstituent::AliReducedJetConstituent() :
	TObject(),
	fPx(0),
	fPy(0),
	fPz(0),
	fE(0),
	fPdgCode(0)
{
}

/**
 * Main constructor, initializing constituent with the 4-momentum of the particle and the PDG code.
 *
 * \param px x-component of the particle 3-momentum vector
 * \param py y-component of the particle 3-momentum vector
 * \param pz z-component of the particle 3-momentum vector
 * \param e Particle energy
 * \param pdg PDG code
 */
AliReducedJetConstituent::AliReducedJetConstituent(double px, double py, double pz, double e, int pdg) :
	TObject(),
	fPx(px),
	fPy(py),
	fPz(pz),
	fE(e),
	fPdgCode(pdg)
{
}

/**
 * Access to the particle 4-momentum vector via filling a TLorentzVector
 *
 * \param target Lorentz-vector filled by this function
 */
void AliReducedJetConstituent::FillLorentzVector(TLorentzVector& target) const {
	target.SetPxPyPzE(fPx, fPy, fPz, fE);
}

/**
 * Obtain more information about the particle (type, charge, quantum numbers) via the information
 * stored in the ROOT particle database, provided via a TParticlePDG. The information is linked via
 * the PDG code of the particle.
 *
 * \return Extended particle information via TParticlePDG.
 */
TParticlePDG* AliReducedJetConstituent::GetPDGParticle() const {
	return TDatabasePDG::Instance()->GetParticle(fPdgCode);
}

} /* namespace HighPtTracks */
