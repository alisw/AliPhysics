/**
 * \file AliReducedJetConstituent.h
 * \brief Definition of class AliReducedJetConstituent, a minimal stucture for jet
 * constituents associated to a jet by the jet clustering algorithm.
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Jan 28, 2015
 */
#ifndef ALIREDUCEDJETCONSTITUENT_H
#define ALIREDUCEDJETCONSTITUENT_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>

class TLorentzVector;
class TParticlePDG;

/**
 * \namespace HighPtTracks
 * \brief Namespace for classes creating trees of events with jets
 *
 * This namespace contains classes descibing reduced events with high
 * jets. A jet event consists of the following classes.
 *    - AliReducedJetEvent
 *    - AliReducedJetParticle
 *    - AliReducedJetConstituent
 *    - AliReducedMatchedTrack
 *  The task AliHighPtReconstructionEfficiency produces the reduced jet events
 *  and is included in the namespace as well. Also the helper class AliParticleMap
 *  is part of this namespace.
 */
namespace HighPtTracks {

/**
 * \class AliReducedJetConstituent
 * \brief Minimal stucture for jet constituents associated to a jet by the jet clustering algorithm.
 *
 * This class represents the minimal structure for constituent found by the clusterizing
 * algorithm. It extends the information of a PseudoJet in fastjet with particle type
 * information coming from the ROOT particle database, obtained by the PDG code. Jet constituents
 * are allowed to be charged or neutral.
 */
class AliReducedJetConstituent : public TObject {
public:
	AliReducedJetConstituent();
	AliReducedJetConstituent(double px, double py, double pz, double e, int pdg);
	virtual ~AliReducedJetConstituent() {}

	/**
	 * Set the particle 4-momentum vector
	 *
	 * \param px x-component of the momentum vector
	 * \param py y-component of the momentum vector
	 * \param pz z-component of the momentum vector
	 * \param e Particle energy
	 */
	void Set(double px, double py, double pz, double e){
		fPx = px;
		fPy = py;
		fPz = pz;
		fE = e;
	}

	/**
	 * Set the particle PDG code
	 *
	 * \param pdg PDG code of the particle
	 */
	void SetPdgCode(int pdg){ fPdgCode = pdg; }

	void FillLorentzVector(TLorentzVector &target) const;

	/**
	 * Get the PDG code of the particle
	 *
	 * \return PDG code of the particle
	 */
	int GetPdgCode() const { return fPdgCode; }

	TParticlePDG *GetPDGParticle() const;

private:
	double 			fPx;          ///< x-compoent of the constituent 3-momentum vector
	double 			fPy;          ///< y-compoent of the constituent 3-momentum vector
	double 			fPz;          ///< z-compoent of the constituent 3-momentum vector
	double 			fE;           ///< Energy of the constituent
	int 			  fPdgCode;     ///< PDG code of the constituent

	/// \cond CLASSIMP
	ClassDef(AliReducedJetConstituent, 1);
	/// \endcond
};

} /* namespace HighPtTracks */

#endif /* PWGJE_EMCALJETTASKS_TRACKS_ALIREDUCEDJETCONSTITUENT_H_ */
