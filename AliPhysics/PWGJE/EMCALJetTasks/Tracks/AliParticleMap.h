/**
 * \file AliParticleMap.h
 * \brief Helper objects for mapping MC-particles to reconstructed particles.
 *
 * This file contains helper classes mapping MC-particles to reconstructed particels.
 *  - AliParticleList contains a list of AliVParticles with the same features
 *  - AliParticlePair contains a pair of Monte-Carlo particle and a list of all reconstructed
 *    particles with the same label
 *  - AliParticleMap maps all reconstructed particles to one label
 *
 *  \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 *  \date Jan 28, 2015
 */
#ifndef ALIPARTICLEMAP_H
#define ALIPARTICLEMAP_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <cstdlib>
#include <vector>
#include <map>

class AliVParticle;
class AliVTrack;

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
 * \class AliParticleList
 * \brief Container of reconstructed particles
 *
 * This class is a container of pointers of reconstructed particles. It can be used in the paritcle map associating
 * reconstructed particles to Monte-Carlo particles. The class does not take ownership over pointers.
 */
class AliParticleList{
public:
  /**
   * Constructor
   */
	AliParticleList():
		fParticles()
	{}

	/**
	 * Destructor
	 */
	~AliParticleList(){}

	/**
	 * Add new particle to the list
	 * \param track Track to be added
	 */
	void AddParticle(AliVTrack *track) { fParticles.push_back(track); }

	/**
	 * Access to particle in the list with a given index
	 * \param itrack Index of the particle
	 * \return Track at the given index
	 */
	AliVTrack *GetParticle(int itrack) const { return fParticles[itrack]; }

	/**
	 * Get the number of particles stored in the list
	 * \return Number of particles
	 */
	int GetNumberOfParticles() const { return fParticles.size(); }

private:

	std::vector<AliVTrack *> 			fParticles;       ///< Vector of reconstructed particles
};

/**
 * \class AliReconstructedParticlePair
 * \brief Pair of a Monte-Carlo true particle and the associated reconstructed information
 *
 * Helper class mapping combining all reconstructed particles which can be associated to one generated
 * particle.
 */
class AliReconstructedParticlePair{
public:
  /**
   *  Constructor
   */
	AliReconstructedParticlePair():
		fTrueParticle(NULL),
		fRecParticles()
	{}

	/**
	 * Copy constructor. Creates new pair from a reference object. Does not take ownership of pointers.
	 *
	 * \param ref Reference for the copy
	 */
	AliReconstructedParticlePair(const AliReconstructedParticlePair &ref):
		fTrueParticle(ref.fTrueParticle),
		fRecParticles(ref.fRecParticles)
	{}

	/**
	 * Assignment operator. Assigns values stored in reference object to this object.
	 *
	 * \param ref Reference for the copy
	 * \return This object
	 */
	AliReconstructedParticlePair &operator=(const AliReconstructedParticlePair &ref){
		if(this != &ref){
			fTrueParticle = ref.fTrueParticle;
			fRecParticles = ref.fRecParticles;
		}
		return *this;
	}

	/**
	 * Destructor
	 */
	~AliReconstructedParticlePair(){}

	/**
	 * Access to true particle
	 * \return The true MC particle
	 */
	AliVParticle *GetMCTrueParticle() const { return fTrueParticle; }

	/**
	 * Access to reconstructed particles
	 * \return List of reconstructed particles matched to this particle
	 */
	const AliParticleList &GetRecTracks() const { return fRecParticles; }

	/**
	 * Set the particle at generator level
	 * \param part MC-true particle
	 */
	void SetMCTrueParticle(AliVParticle *const part) { fTrueParticle = part; }

	/**
	 * Set the list of reconstructed particles associated to this particle
	 * \param tracks List of reconstructed tracks
	 */
	void SetRecParticles(const AliParticleList &tracks) { fRecParticles = tracks; }

private:

	AliVParticle 			       *fTrueParticle;            ///< True selected particle
	AliParticleList 	        fRecParticles;            ///< List of all matched particles according to the Monte-Carlo label
};

/**
 * \class AliParticleMap
 * \brief Map of reconstructed particles which share the same Monte-Carlo label
 *
 * Class connecting all reconstructing particles sharing the same Monte-Carlo label. Used for
 * a fast search of reconstructed tracks for given generated tracks.
 */
class AliParticleMap{
public:
  /**
   * Constructor
   */
	AliParticleMap():
		fParticles()
	{}
	~AliParticleMap();

	void AddParticle(AliVTrack *track);
	AliParticleList *GetParticles(int label) const;

	/**
	 * Get the number of true particles (according to stored labels)
	 *
	 * \return Number of entries
	 */
	int GetNumberOfParticles() const { return fParticles.size(); }

	void Print() const;

private:
	std::map<int, AliParticleList *> 				fParticles;         ///< connection of particles to labels
};

} /* namespace HighPtTracks */

#endif /* ALIPARTICLEMAP_H */
