/**
 * \file AliEMCalPtTaskTrackSelectionAOD.h
 * \brief Implement virtual track selection for AOD analysis
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Dec 12, 2014
 */
#ifndef ALIEMCALPTTRACKSELECTIONAOD_H_
#define ALIEMCALPTTRACKSELECTIONAOD_H_
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliESDtrackCuts.h"
#include "AliEMCalPtTaskVTrackSelection.h"

class AliVTrack;

/**
 * \namespace EMCalTriggerPtAnalysis
 * \brief Analysis of high-\f$ p_{t} \f$ tracks in triggered events
 *
 * This namespace contains classes for the analysis of high-\f$ p_{t} \f$ tracks in
 * triggered events.
 */
namespace EMCalTriggerPtAnalysis {

/**
 * \class AliEMCalPtTaskTrackSelectionAOD
 * \brief Implement virtual track selection for AOD analysis
 *
 * Implementation of track selection in case the analysis runs on AODs
 * For the moment it uses the AliESDtrackCuts and converts AOD tracks to
 * ESD tracks, which might change in the future when an AOD track selection
 * framework becomes available.
 */
class AliEMCalPtTaskTrackSelectionAOD: public AliEMCalPtTaskVTrackSelection {
public:
	AliEMCalPtTaskTrackSelectionAOD();
	AliEMCalPtTaskTrackSelectionAOD(AliESDtrackCuts *cuts, UInt_t filterbits);
	virtual ~AliEMCalPtTaskTrackSelectionAOD() {}

	virtual TObjArray *GetAcceptedTracks(const TClonesArray * const tracks);
	virtual TObjArray *GetAcceptedTracks(const AliVEvent *const event);
	virtual bool IsTrackAccepted(AliVTrack * const trk);

	/**
	 * Add a new filter bit to the track selection. Multiple filter bits can be set
	 * at the same time (via the bitwise or operator |).
	 *
	 * \param filterbits
	 */
	void AddFilterBit(UInt_t filterbits) { fFilterBits |= filterbits; }


private:
	UInt_t			fFilterBits;				    ///< Track filter bits

	/// \cond CLASSIMP
	ClassDef(AliEMCalPtTaskTrackSelectionAOD, 1);		// Track selection class for AOD analysis
	/// \endcond
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALPTTRACKSELECTIONAOD_H_ */
