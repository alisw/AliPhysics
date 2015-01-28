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
 * Object storing reduced jet constituent information
 *
 * Author:
 *   Markus Fasel <markus.fasel@cern.ch>
 */
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TParticlePDG.h>

#include "AliReducedJetConstituent.h"

ClassImp(HighPtTracks::AliReducedJetConstituent)

namespace HighPtTracks {

AliReducedJetConstituent::AliReducedJetConstituent() :
	TObject(),
	fPx(0),
	fPy(0),
	fPz(0),
	fE(0),
	fPdgCode(0)
{
}

AliReducedJetConstituent::AliReducedJetConstituent(double px, double py, double pz, double e, int pdg) :
	TObject(),
	fPx(px),
	fPy(py),
	fPz(pz),
	fE(e),
	fPdgCode(pdg)
{
}

void AliReducedJetConstituent::FillLorentzVector(TLorentzVector& target) const {
	target.SetPxPyPzE(fPx, fPy, fPz, fE);
}

TParticlePDG* AliReducedJetConstituent::GetPDGParticle() const {
	return TDatabasePDG::Instance()->GetParticle(fPdgCode);
}

} /* namespace HighPtTracks */
