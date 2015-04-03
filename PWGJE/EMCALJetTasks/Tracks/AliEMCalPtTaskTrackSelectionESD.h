/**
 * \file AliEMCalPtTaskTrackSelectionESD.h
 * \brief Declaration of class AliEMCalPtTaskTrackSelectionESD
 *
 * In this header file the class AliEMCalPtTaskTrackSelectionESD, which implements
 * the virtual track selection for ESD tracks, is declared
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Dec 12, 2014
 */
#ifndef ALIEMCALPTTASKTRACKSELECTIONESD_H_
#define ALIEMCALPTTASKTRACKSELECTIONESD_H_
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <AliESDtrackCuts.h>
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
 * \class AliEMCalPtTaskTrackSelectionESD
 * \brief Implementation of virtual track selection for ESDs
 *
 * Implementation of the track selection for the analysis on ESDs using
 * AliESDtrackCuts as underlying structure
 */
class AliEMCalPtTaskTrackSelectionESD: public AliEMCalPtTaskVTrackSelection {
public:
	AliEMCalPtTaskTrackSelectionESD();
	AliEMCalPtTaskTrackSelectionESD(AliESDtrackCuts *cuts);
	virtual ~AliEMCalPtTaskTrackSelectionESD() {}

	virtual TObjArray *GetAcceptedTracks(const TClonesArray * const tracks);
	virtual TObjArray *GetAcceptedTracks(const AliVEvent *const event);
	virtual bool IsTrackAccepted(AliVTrack * const trk);

	/// \cond CLASSIMP
	ClassDef(AliEMCalPtTaskTrackSelectionESD,1);	// Selection of ESD tracks for analysis
	/// \endcond
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALPTTASKTRACKSELECTIONESD_H_ */
