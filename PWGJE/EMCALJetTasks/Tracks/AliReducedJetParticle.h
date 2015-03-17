/**
 * \class AliReducedJetParticle
 * \brief Reduced information set of particles associated with a jet
 *
 * Object storing particle based information for tracks associated with jets. The object
 * contains the 4-momentum vector of the particle, the PDG code, the distance to the main jet
 * axis of its associated jet, the number of TPC track references of the Monte-Carlo true
 * particle, and a list of reduced reconstructed tracks associated to the true particle.
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Jan 28, 2015
 */
#ifndef ALIREDUCEDJETPARTICLE_H
#define ALIREDUCEDJETPARTICLE_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObjArray.h>
#include "AliReducedMatchedTrack.h"

class TLorentzVector;

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

class AliReducedJetParticle : public TObject {
public:
	AliReducedJetParticle();
	AliReducedJetParticle(double px, double py, double pz, double e, int pdgcode);
	AliReducedJetParticle(const AliReducedJetParticle &ref);
	AliReducedJetParticle &operator=(const AliReducedJetParticle &ref);
	virtual ~AliReducedJetParticle();

	void FillLorentzVector(TLorentzVector &ref) const;

	/**
	 * \brief Get distance to the main jet axis.
	 *
	 * Get distance to the main jet axis.
	 *
	 * \return The reconstructed distance
	 */
	double GetDistanceToJetMainAxis() const { return fDr; }

	/**
	 * \brief Get the PDG code of the particle.
	 *
	 * Get the PDG code of the particle.
	 *
	 * \return The PDF code
	 */
	int GetPdgCode() const { return fPdgCode; }

	/**
	 * \brief Checks if the track is reconstructed.
	 *
	 * Checks if the track is reconstructed. Reconstructed tracks are defined as tracks
	 * with at least one associated track.
	 *
	 * \return True if at least one reconstructed track is found, false otherwise
	 */
	bool IsReconstructed() const { return fMatchedTrackContainer->GetEntries() > 0; }

	/**
	 * \brief Get the number of matched reconstructed tracks
	 *
	 * Get the number of matched reconstruced tracks associated to the particle.
	 *
	 * \return Number of matched reconstructed tracks
	 */
	int GetNumberOfMatchedTracks() const { return fMatchedTrackContainer->GetEntries(); }

	double GetDeltaPt(int itrk = 0) const;

	/**
	 * \brief Get the number of track references in the TPC for the particle
	 *
	 * Get the number of track references in the TPC for the particle
	 *
	 * \return The number of TPC track references
	 */
	unsigned short GetNumberOfTPCTrackReferences() const { return fNTPCTrackReferences; }

	TObjArray *GetMatchedTracks() const { return fMatchedTrackContainer; }

	/**
	 * \brief Get matched track with given index
	 *
	 * Get matched track with given index. In case the index is out of bounds,
	 * return a nullpointer.
	 *
	 * \param itrk Index to check
	 * \return The track with this index (NULL if out of bounds)
	 */
	AliReducedMatchedTrack *GetMatchedTrack(int itrk){
	  if(itrk >= fMatchedTrackContainer->GetEntries()) return NULL;
	  return static_cast<AliReducedMatchedTrack *>(fMatchedTrackContainer->At(itrk));
	}

	/**
	 * \brief Set the particle 4-momentum vector
	 *
	 * Set the particle 4-momentum vector
	 *
	 * \param px x-component of the momentum vector
	 * \param py y-component of the momentum vector
	 * \param pz z-component of the momentum vector
	 * \param e Particle energy
	 */
	void SetKine(double px, double py, double pz, double e){
		fPx = px;
		fPy = py;
		fPz = pz;
		fE = e;
	}

	/**
	 * \brief Set the distance to the main jet axis
	 *
	 * Set the distance to the main jet axis
	 *
	 * \param dr Distance to the main jet axis
	 */
	void SetDistanceToMainJetAxis(double dr) { fDr = dr; }

	/**
	 * \brief Set the particle PDG code
	 *
	 * Set the particle PDG code
	 *
	 * \param pdg The particle PDG code
	 */
	void SetPdgCode(int pdg) { fPdgCode = pdg; }

	void AddMatchedTrack(AliReducedMatchedTrack *trk);

	/**
	 * \brief Set the number of track references associated to the particle
	 *
	 * Set the number of track references associated to the particle
	 *
	 * \param nref Number of track references associated to the particle
	 */
	void SetNumberOfTPCtrackReferences(unsigned short nref) { fNTPCTrackReferences = nref; }

private:
	double 			fPx;              /// x-component of the momentum vector
	double 			fPy;              /// y-component of the momentum vector
	double			fPz;              /// z-component of the momentum vector
	double			fE;               /// Particle energy
	double			fDr;              /// Distance to the main jet axis
	int				fPdgCode;           /// PDG code of the particle
	unsigned short fNTPCTrackReferences;  /// Number of TPC track references associated to the particle
	TObjArray *fMatchedTrackContainer;    /// Container for matched tracks at reconstruction level

	/// \cond CLASSIMP
	ClassDef(AliReducedJetParticle, 3);
	/// \endcond
};

} /* namespace HighPtTracks */

#endif /* ALIREDUCEDJETPARTICLE_H */
