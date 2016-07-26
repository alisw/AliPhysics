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
 * \ingroup EMCALCOREFW
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Jul 24, 2015
 *
 * Implementation of the track selection for the analysis on ESDs using
 * AliESDtrackCuts as underlying structure
 */
class AliEmcalTrackSelectionESD: public AliEmcalTrackSelection {
public:
	AliEmcalTrackSelectionESD();
	AliEmcalTrackSelectionESD(AliVCuts *cuts);
	AliEmcalTrackSelectionESD(ETrackFilterType_t type, const char* period = "");
	virtual ~AliEmcalTrackSelectionESD() {}

	virtual void GenerateTrackCuts(ETrackFilterType_t type, const char* period = "");

	virtual bool IsTrackAccepted(AliVTrack * const trk);

	/// \cond CLASSIMP
	ClassDef(AliEmcalTrackSelectionESD,1);
	/// \endcond
};

#endif /* ALIEMCALPTTASKTRACKSELECTIONESD_H_ */
