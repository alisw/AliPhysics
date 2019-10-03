/**
 * \file AliHighPtReconstructionEfficiency.h
 * \brief Definition of the analysis task producing filtered trees with reconstructed
 * jets at generator level
 *
 * \author: Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Jan 28, 2015
 */
#ifndef ALIHIGHPTRECONSTRUCTIONEFFICIENCY_H
#define ALIHIGHPTRECONSTRUCTIONEFFICIENCY_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h"

#include <map>
#include <vector>
#include <fastjet/PseudoJet.hh>
#include "AliParticleMap.h"

class AliESDtrackCuts;
class AliVParticle;
class AliVTrack;
class TList;
class TNtuple;

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

class AliReducedJetEvent;
class AliReducedJetInfo;

/**
 * \class AliHighPtReconstructionEfficiency
 * \brief Analysis task producing filtered trees with reconstructed jets at generator level
 *
 * Analysis task creating reduced jet events. Each event contains reduced jets, which associate
 * constituents and reconstructed track. Constituents are particles associated by the jet finder
 * while tracks are all charged particles in a cone within a maximum radius (can be larger than R set
 * in the jet finder).
 */
class AliHighPtReconstructionEfficiency : public AliAnalysisTaskSE {
public:
  /**
   * \enum CutType_t
   * \brief Declaration of cut types
   */
  enum CutType_t{
    kRJStandardCuts = 0,      ///< Standard track cuts
    kRJHybridCuts = 1         ///< Hybrid track cuts
  };
	AliHighPtReconstructionEfficiency();
	AliHighPtReconstructionEfficiency(const char *name);
	virtual ~AliHighPtReconstructionEfficiency();

	virtual void UserCreateOutputObjects();
	virtual bool UserNotify();
	virtual void UserExec(Option_t* /*option*/);
	void Terminate(Option_t* /*option*/) {}

	/**
	 * Set maximum \f$ \eta \f$ for particles and reconstructed tracks
	 *
	 * \param maxeta The maximum allowed \f$ \eta \f$
	 */
	void SetMaxEtaParticles(double maxeta) { fMaxEtaParticles = maxeta; }

	/**
	 * Set maximum \f$ \eta \f$ for jets
	 *
	 * \param maxeta The maximum allowed \f$ \eta \f$
	 */
	void SetMaxEtaJets(double maxeta) { fMaxEtaJets = maxeta; }

	/**
	 * Set minimum allowed \f$ p_{t} \f$ for particles and tracks
	 *
	 * \param minpt
	 */
	void SetMinPtTracks(double minpt) { fMinPtParticles = minpt; }

	/**
	 * \brief Set Maximum allowed distance for particles associated with a jet to the main jet axis
	 *
	 * \param maxdr Maximum allowed distance
	 */
	void SetMaxDR(double maxdr) { fMaxDR = maxdr; }

	/**
	 * Set standard track cuts
	 *
	 * \param cuts Standard track cuts
	 */
	void SetStandardTrackCuts(AliESDtrackCuts *const cuts) { SetTrackCuts(cuts, kRJStandardCuts); }

	/**
	 * Set hybrid track cuts
	 *
	 * \param cuts Hybrid track cuts
	 */
	void SetHybridTrackCuts(AliESDtrackCuts *const cuts) { SetTrackCuts(cuts, kRJHybridCuts); }

	/**
	 * Set task into debug mode
	 */
	void SetTaskDebugMode() { fTaskDebugMode = true; }

protected:
	void SetTrackCuts(AliESDtrackCuts *const cuts, CutType_t cuttype) { fTrackCuts[cuttype] = cuts; }
	bool IsSelected(const AliVTrack * const track, CutType_t type) const;
	bool IsTrueSelected(const AliVParticle *const track) const;
	void SelectParticlesForJetfinding(TList &particles) const;
	void CreateRectrackLookup();
	std::vector<AliReconstructedParticlePair> SelectParticles() const;
	double GetDR(const fastjet::PseudoJet &recjet, const AliVParticle *const inputtrack) const;
	const AliParticleList * FindReconstructedParticleFast(int label) const;
	void ProcessJet(AliReducedJetInfo *const jet, const std::vector<AliReconstructedParticlePair> &particles) const;
	void ConvertConstituents(AliReducedJetInfo * const recjet, const fastjet::PseudoJet &inputjet);
	bool IsPhysicalPrimary(const AliVParticle *const part) const;
	bool PythiaInfoFromFile(const char* currFile, double &fXsec, double &fTrials, int &pthard) const ;
	unsigned short GetNumberOfTPCTrackReferences(AliVParticle *const trk) const;

private:
	AliHighPtReconstructionEfficiency(const AliHighPtReconstructionEfficiency &);
	AliHighPtReconstructionEfficiency &operator=(const AliHighPtReconstructionEfficiency &);

	AliESDtrackCuts					*fTrackCuts[2];				///< List of track cuts
	/// Map of reconstructed particles associate to a Monte-Carlo label
	AliParticleMap					*fParticleMap;				//!

	// Output objects
	TTree							        *fJetTree;          ///< Output tree
	AliReducedJetEvent				*fJetEvent;         ///< Output jet event

	double 							fMaxEtaJets;              ///< \f$ \eta \f$ cut for jets
	double							fMaxEtaParticles;         ///< \f$ \eta \f$ cut for particles
	double 							fMinPtParticles;          ///< minimium \f$ p_{t} \f$ cut for particles
	double							fMaxDR;                   ///< maximum distance of a particle to the main jet axis

	double							fCrossSection;            ///< Cross section from PYTHIA
	double							fNtrials;                 ///< Number of trials from PYTHIA
	int								fPtHardBin;                 ///< \f$ p_{t} \f$-hard bin of the event

	bool 							fTaskDebugMode;             ///< Switch for debug mode

	/// \cond CLASSIMP
	ClassDef(AliHighPtReconstructionEfficiency, 1);
	/// \endcond
};

}
#endif /* ALIHIGHPTRECONSTRUCTIONEFFICIENCY_H */
