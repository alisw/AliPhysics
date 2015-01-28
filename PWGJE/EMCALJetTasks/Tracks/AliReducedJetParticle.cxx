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
/*
 * Object storing particle based information for tracks associated with jets
 *
 * Author:
 *   Markus Fasel <markus.fasel@cern.ch>
 */
#include <TLorentzVector.h>

#include "AliReducedJetParticle.h"

ClassImp(HighPtTracks::AliReducedJetParticle)

namespace HighPtTracks {

AliReducedJetParticle::AliReducedJetParticle() :
	TObject(),
	fPx(0),
	fPy(0),
	fPz(0),
	fE(0),
	fDr(0),
	fPdgCode(0),
	fIsReconstructed(false)
{
}

AliReducedJetParticle::AliReducedJetParticle(double px, double py, double pz, double e, int pdgcode, bool rec) :
	TObject(),
	fPx(px),
	fPy(py),
	fPz(pz),
	fE(e),
	fDr(0),
	fPdgCode(pdgcode),
	fIsReconstructed(rec)
{
}

void AliReducedJetParticle::FillLorentzVector(TLorentzVector& ref) const {
	/*
	 * Access particle kine information
	 */
	ref.SetPxPyPzE(fPx, fPy, fPz, fE);
}

} /* namespace HighPtTracks */
