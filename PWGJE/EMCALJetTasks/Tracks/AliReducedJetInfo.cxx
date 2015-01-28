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
 * Class with reduced jet information for the reduced jet tree
 *
 * Author:
 *   Markus Fasel <markus.fasel@cern.ch>
 */

#include <TLorentzVector.h>
#include <TObjArray.h>

#include "AliReducedJetConstituent.h"
#include "AliReducedJetInfo.h"
#include "AliReducedJetParticle.h"

ClassImp(HighPtTracks::AliReducedJetInfo)

namespace HighPtTracks {

AliReducedJetInfo::AliReducedJetInfo() :
	TObject(),
	fPx(0.),
	fPy(0.),
	fPz(0.),
	fE(0),
	fConstituents(NULL),
	fParticlesInCone(NULL)
{
}

AliReducedJetInfo::AliReducedJetInfo(double px, double py, double pz, double e) :
	TObject(),
	fPx(px),
	fPy(py),
	fPz(pz),
	fE(e),
	fConstituents(NULL),
	fParticlesInCone(NULL)
{
	fParticlesInCone = new TObjArray;
	fParticlesInCone->SetOwner(true);
	fConstituents = new TObjArray;
	fConstituents->SetOwner(true);
}

AliReducedJetInfo::AliReducedJetInfo(const AliReducedJetInfo& ref) :
	TObject(ref),
	fPx(ref.fPx),
	fPy(ref.fPy),
	fPz(ref.fPz),
	fE(ref.fE),
	fConstituents(NULL),
	fParticlesInCone(NULL)
{
	fParticlesInCone = new TObjArray();
	fParticlesInCone->SetOwner(true);
	fConstituents = new TObjArray;
	fConstituents->SetOwner(true);

	AliReducedJetParticle *mypart(NULL);
	TIter particleIter(ref.fParticlesInCone);
	while((mypart = dynamic_cast<AliReducedJetParticle *>(particleIter()))){
		fParticlesInCone->Add(new AliReducedJetParticle(*mypart));
	}
	AliReducedJetConstituent *myconst(NULL);
	TIter constIter(ref.fConstituents);
	while((myconst = dynamic_cast<AliReducedJetConstituent *>(constIter()))){
		fConstituents->Add(new AliReducedJetConstituent(*myconst));
	}
}

AliReducedJetInfo& AliReducedJetInfo::operator=(const AliReducedJetInfo& ref) {
	TObject::operator=(ref);
	if(&ref != this){
		fPx = ref.fPx;
		fPy = ref.fPy;
		fPz = ref.fPz;
		fE = ref.fE;

		fParticlesInCone->Clear();
		AliReducedJetParticle *mypart(NULL);
		TIter particleIter(ref.fParticlesInCone);
		while((mypart = dynamic_cast<AliReducedJetParticle *>(particleIter()))){
			fParticlesInCone->Add(new AliReducedJetParticle(*mypart));
		}

		fConstituents->Clear();
		AliReducedJetConstituent *myconst(NULL);
		TIter constIter(ref.fConstituents);
		while((myconst = dynamic_cast<AliReducedJetConstituent *>(constIter()))){
			fConstituents->Add(new AliReducedJetConstituent(*myconst));
		}
	}
	return *this;
}

AliReducedJetInfo::~AliReducedJetInfo() {
	delete fParticlesInCone;
}

void AliReducedJetInfo::AddParticleInCone(AliReducedJetParticle* part) {
	if(fParticlesInCone) fParticlesInCone->Add(part);
}

void AliReducedJetInfo::FillLorentzVector(TLorentzVector& vec) const {
	vec.SetPxPyPzE(fPx, fPy, fPx, fE);
}

int AliReducedJetInfo::GetNumberOfMatchedParticles() const {
	if(!fParticlesInCone) return 0;
	return fParticlesInCone->GetEntries();
}

void AliReducedJetInfo::AddConstituent(AliReducedJetConstituent* con) {
	if(fConstituents) fConstituents->Add(con);
}

AliReducedJetParticle* AliReducedJetInfo::GetMatchedParticle(int ipart) const {
	if(!fParticlesInCone || ipart >= fParticlesInCone->GetEntries()) return NULL;
	return dynamic_cast<AliReducedJetParticle *>(fParticlesInCone->At(ipart));
}

} /* namespace HighPtTracks */
