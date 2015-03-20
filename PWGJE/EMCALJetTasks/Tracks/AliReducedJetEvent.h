/**
 * \file AliReducedJetEvent.h
 * \brief Definition of class AliReducedJetEvent, an event structure containing reduced jet information.
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Jan 28, 2015
 */
#ifndef ALIREDUCEDJETEVENT_H
#define ALIREDUCEDJETEVENT_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>

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

class AliReducedJetInfo;

/**
 * \class AliReducedJetEvent
 * \brief Event structure containing reduced jet information.
 *
 * This class represents an event structure for jets at generator level, which links true jets
 * to particles within the jet cone in order to study their properties. The event contains
 *  -# A list of reconstructed reduced jets at generator level
 *  -# Calculated cross section of the event
 *  -# Number of trials
 *  -# Generated \f$ p_{t} \f$ of the hard interaction
 */
class AliReducedJetEvent : public TObject {
public:
	AliReducedJetEvent();
	AliReducedJetEvent(double crosssection, int ntrials, double pthard);
	AliReducedJetEvent(const AliReducedJetEvent &ref);
	AliReducedJetEvent &operator=(const AliReducedJetEvent &);
	virtual ~AliReducedJetEvent();

	/**
	 * Set values calculated by pythia (Cross section, number of trials, generated \f$ p_{t} \f$
	 * of the hard interaction)
	 *
	 * \param crosssection Event cross section
	 * \param ntrials Number of trials
	 * \param pthard Generated \f$ p_{t} \f$ of the hard interaction
	 */
	void SetPythiaHardInfo(double crosssection, int ntrials, double pthard){
		fCrossSection = crosssection;
		fTrials = ntrials;
		fPtHard = pthard;
	}
	void AddReconstructedJet(AliReducedJetInfo *jet);

	/**
	 * Get the event cross section.
	 *
	 * \return The event cross section
	 */
	double GetCrossSection() const { return fCrossSection; }

	/**
	 * Get the number of trials.
	 *
	 * \return The number of trials
	 */
	int GetNumberOfTrials() const { return fTrials; }

	/**
	 * Get the \f$ p_{t} \f$ of the hard interaction.
	 * \return
	 */
	double GetPtHard() const { return fPtHard; }

	int GetNumberOfJets() const;

	AliReducedJetInfo *GetReconstructedJet(int ijet) const;

	/**
	 * Get a list of reconstructed jets found in this event.
	 *
	 * \return List of jets found in this event.
	 */
	TObjArray *GetListOfJets() const { return fReconstructedJets; }

private:
	double					fCrossSection;            ///< Event cross section
	int 					  fTrials;                  ///< Number of trials
	double					fPtHard;                  ///< Generated \f$ p_{t} \f$ of the hard interaction
	TObjArray				*fReconstructedJets;      ///< List of reconstructed jets at generator level

	/// \cond CLASSIMP
	ClassDef(AliReducedJetEvent, 1);
	/// \endcond
};

} /* namespace HighPtTracks */

#endif /* ALIREDUCEDJETEVENT_H */
