/**
 * \file AliReducedJetInfo.h
 * \brief Definition of class AliReducedJetInfo, a structure for reduced information
 * about a reconstructed jet
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Jan 28, 2015
 */
#ifndef ALIREDUCEDJETINFO_H
#define ALIREDUCEDJETINFO_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>

class TLorentzVector;
class TObjArray;

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

class AliReducedJetParticle;
class AliReducedJetConstituent;

/**
 * \class AliReducedJetInfo
 * \brief Reduced information about a reconstructed jet
 *
 * Class with reduced jet information for the reduced jet tree. For a jet, the following information is
 * stored:
 *  -# 4-momentum vector of the jet
 *  -# list of reduced constituents found by the jet clustering algorithm
 *  -# list of all charged particles with reduced information within a certain distance
 *     to the main jet axis.
 */
class AliReducedJetInfo  : public TObject {
public:
	AliReducedJetInfo();
	AliReducedJetInfo(double px, double py, double pz, double e);
	AliReducedJetInfo(const AliReducedJetInfo &ref);
	AliReducedJetInfo &operator=(const AliReducedJetInfo &ref);
	virtual ~AliReducedJetInfo();

	/**
	 * Set the jet 4-momentum vector
	 *
	 * \param px x-component of the 4-momentum vector
	 * \param py y-component of the 4-momentum vector
	 * \param pz z-component of the 4-momentum vector
	 * \param e jet energy
	 */
	void Set(double px, double py, double pz, double e){
		fPx = px;
		fPy = py;
		fPz = pz;
		fE = e;
	}
	void AddParticleInCone(AliReducedJetParticle *part);
	void AddConstituent(AliReducedJetConstituent *con);

	/**
	 * Get the jet 4-vector and fill it into the function parameter
	 *
	 * \param px x-component of the jet 4-momentum vector, filled by the function
	 * \param py y-component of the jet 4-momentum vector, filled by the function
	 * \param pz z-component of the jet 4-momentum vector, filled by the function
	 * \param e jet energy, filled by the function
	 */
	void GetPxPyPxE(double &px, double &py, double &pz, double &e){
		px = fPx;
		py = fPy;
		pz = fPz;
		e = fE;
	}
	void FillLorentzVector(TLorentzVector &vec ) const;
	int GetNumberOfMatchedParticles() const;
	TObjArray *GetListOfMatchedParticles() const { return fParticlesInCone; }
	TObjArray *GetListOfConstituents() const { return fConstituents; }
	AliReducedJetParticle *GetMatchedParticle(int ipart) const;

private:
	double					fPx;                      ///< x-component of the 4-momentum vector
	double 					fPy;                      ///< y-component of the 4-momentum vector
	double					fPz;                      ///< z-component of the 4-momentum vector
	double					fE;                       ///< reconstructed jet energy

	TObjArray				*fConstituents;           ///< Constituents found by the jet clustering algorithm
	TObjArray				*fParticlesInCone;        ///< Particles associated to this jet via distance to the main jet axis

	/// \cond CLASSIMP
	ClassDef(AliReducedJetInfo, 1);
	/// \endcond
};

} /* namespace HighPtTracks */

#endif /* ALIREDUCEDJETINFO_H */
