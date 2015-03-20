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
#include <TLorentzVector.h>
#include <TObjArray.h>

#include "AliReducedJetConstituent.h"
#include "AliReducedJetInfo.h"
#include "AliReducedJetParticle.h"

/// \cond CLASSIMP
ClassImp(HighPtTracks::AliReducedJetInfo)
/// \endcond

namespace HighPtTracks {

/***
 * Dummy (I/O) constructor, not to be used by the user. Containers are not allocated.
 */
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

/**
 * Main constructor of the reduced jet information of a reconstuced jet. Initializes the
 * 4-momentum vector. Also containers for associated particles and jet constituents are
 * allocated.
 *
 * \param px x-component of the 4-momentum vector
 * \param py y-component of the 4-momentum vector
 * \param pz z-component of the 4-momentum vector
 * \param e Reconstructed jet energy
 */
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

/**
 * Copy constructor, initializes a new jet info from a reference object. For constituents and
 * associated particles a deep copy is performed, meaning the new object keeps ownership over
 * its pointers.
 *
 * \param ref Reference object for the copy
 */
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

/**
 * Assignment operator, copying status of the reference reduced jet into this jet. For constituents and
 * associated particles a deep copy is performed, meaning the new object keeps ownership over
 * its pointers.
 *
 * \param ref Reference object assinged to this reduced jet
 * \return This object
 */
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

/**
 * Destructor, cleaning object. All containers are deleted.
 */
AliReducedJetInfo::~AliReducedJetInfo() {
	delete fParticlesInCone;
	delete fConstituents;
}

/**
 * Add new associated particle (expanded information) to the jet info object. Jet info takes
 * ownership over the object.
 *
 * \param part Particle to be added to the jet container.
 */
void AliReducedJetInfo::AddParticleInCone(AliReducedJetParticle* part) {
	if(fParticlesInCone) fParticlesInCone->Add(part);
}

/**
 * Fills a TLorentzVector with jet 4-momentum.
 *
 * \param vec The vector to be filled.
 */
void AliReducedJetInfo::FillLorentzVector(TLorentzVector& vec) const {
	vec.SetPxPyPzE(fPx, fPy, fPx, fE);
}

/**
 * Get the number of particles matched to this jet.
 *
 * \return Number of matched particles
 */
int AliReducedJetInfo::GetNumberOfMatchedParticles() const {
	if(!fParticlesInCone) return 0;
	return fParticlesInCone->GetEntries();
}

/**
 * Adds a constituent (charged or neutral) as associated by the clustering algorithm to the
 * reduced jet.
 *
 * \param con Jet constituent to be added
 */
void AliReducedJetInfo::AddConstituent(AliReducedJetConstituent* con) {
	if(fConstituents) fConstituents->Add(con);
}

/**
 * Provides access to a single particle associated to the jet at a given index in then container.
 *
 * \param ipart Index of the particle
 * \return Particle at a given index (NULL if index is out of range)
 */
AliReducedJetParticle* AliReducedJetInfo::GetMatchedParticle(int ipart) const {
	if(!fParticlesInCone || ipart >= fParticlesInCone->GetEntries()) return NULL;
	return dynamic_cast<AliReducedJetParticle *>(fParticlesInCone->At(ipart));
}

} /* namespace HighPtTracks */
