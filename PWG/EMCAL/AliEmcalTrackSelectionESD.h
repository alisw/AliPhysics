/**
 * \file AliEmcalTrackSelectionESD.h
 * \brief Declaration of class AliEmcalTrackSelectionESD
 *
 * In this header file the class AliEmcalTrackSelectionESD, which implements
 * the virtual track selection for ESD tracks, is declared
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Jul 24, 2015
 */
#ifndef ALIEMCALTASKTRACKSELECTIONESD_H_
#define ALIEMCALTASKTRACKSELECTIONESD_H_
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <AliEmcalTrackSelection.h>

class AliVCuts;
class AliVTrack;

/**
 * \class AliEmcalTrackSelectionESD
 * \brief Implementation of virtual track selection for ESDs
 *
 * Implementation of the track selection for the analysis on ESDs using
 * AliESDtrackCuts as underlying structure
 */
class AliEmcalTrackSelectionESD: public AliEmcalTrackSelection {
public:
	AliEmcalTrackSelectionESD();
	AliEmcalTrackSelectionESD(AliVCuts *cuts);
	virtual ~AliEmcalTrackSelectionESD() {}

	virtual TObjArray *GetAcceptedTracks(const TClonesArray * const tracks);
	virtual TObjArray *GetAcceptedTracks(const AliVEvent *const event);
	virtual bool IsTrackAccepted(AliVTrack * const trk);

	/// \cond CLASSIMP
	ClassDef(AliEmcalTrackSelectionESD,1);
	/// \endcond
};

#endif /* ALIEMCALPTTASKTRACKSELECTIONESD_H_ */
