/**
 * \file AliAODTrackSelection.h
 * \brief Implement virtual track selection for AOD analysis
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Jul 24, 2015
 */
#ifndef ALIAODTRACKSELECTIONAOD_H_
#define ALIAODTRACKSELECTIONAOD_H_
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <AliVTrackSelection.h>

class AliVCuts;
class AliVTrack;

/**
 * \class AliAODTrackSelection
 * \brief Implement virtual track selection for AOD analysis
 *
 * Implementation of track selection in case the analysis runs on AODs
 * For the moment it uses the AliESDtrackCuts and converts AOD tracks to
 * ESD tracks, which might change in the future when an AOD track selection
 * framework becomes available.
 */
class AliAODTrackSelection: public AliVTrackSelection {
public:
	AliAODTrackSelection();
	AliAODTrackSelection(AliVCuts *cuts, UInt_t filterbits);
	virtual ~AliAODTrackSelection() {}

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
	ClassDef(AliAODTrackSelection, 2);
	/// \endcond
};

#endif /* ALIAODTRACKSELECTIONAOD_H_ */
