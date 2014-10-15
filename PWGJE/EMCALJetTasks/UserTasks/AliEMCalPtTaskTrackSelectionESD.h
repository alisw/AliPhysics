#ifndef ALIEMCALPTTASKTRACKSELECTIONESD_H_
#define ALIEMCALPTTASKTRACKSELECTIONESD_H_
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel

#include <AliEMCalPtTaskVTrackSelection.h>
#include <AliESDtrackCuts.h>

namespace EMCalTriggerPtAnalysis {

class AliEMCalPtTaskTrackSelectionESD: public AliEMCalPtTaskVTrackSelection {
public:
	AliEMCalPtTaskTrackSelectionESD();
	AliEMCalPtTaskTrackSelectionESD(AliESDtrackCuts *cuts);
	AliEMCalPtTaskTrackSelectionESD(const AliEMCalPtTaskTrackSelectionESD &ref);
	AliEMCalPtTaskTrackSelectionESD &operator=(const AliEMCalPtTaskTrackSelectionESD &ref);
	virtual ~AliEMCalPtTaskTrackSelectionESD();

	virtual TObjArray *GetAcceptedTracks(const TClonesArray * const tracks);
	virtual TObjArray *GetAcceptedTracks(const AliVEvent *const event);

	void SetTrackCuts(AliESDtrackCuts * cuts) { fTrackCuts = cuts; }
	virtual TObject *GetTrackCuts() { return fTrackCuts; }

private:
	AliESDtrackCuts *fTrackCuts;				// Track cuts

	ClassDef(AliEMCalPtTaskTrackSelectionESD,1);	// Selection of ESD tracks for analysis
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALPTTASKTRACKSELECTIONESD_H_ */
