#ifndef ALIEMCALPTTRACKSELECTIONAOD_H_
#define ALIEMCALPTTRACKSELECTIONAOD_H_
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel

#include "AliESDtrackCuts.h"
#include "AliEMCalPtTaskVTrackSelection.h"

class AliVTrack;

namespace EMCalTriggerPtAnalysis {

class AliEMCalPtTaskTrackSelectionAOD: public AliEMCalPtTaskVTrackSelection {
public:
	AliEMCalPtTaskTrackSelectionAOD();
	AliEMCalPtTaskTrackSelectionAOD(AliESDtrackCuts *cuts, UInt_t filterbits);
	AliEMCalPtTaskTrackSelectionAOD(const AliEMCalPtTaskTrackSelectionAOD &ref);
	AliEMCalPtTaskTrackSelectionAOD &operator=(const AliEMCalPtTaskTrackSelectionAOD &ref);
	virtual ~AliEMCalPtTaskTrackSelectionAOD();

	virtual TObjArray *GetAcceptedTracks(const TClonesArray * const tracks);
	virtual TObjArray *GetAcceptedTracks(const AliVEvent *const event);
	virtual bool IsTrackAccepted(AliVTrack * const trk);

	void AddFilterBit(UInt_t filterbits) { fFilterBits |= filterbits; }
	void SetTrackCuts(AliESDtrackCuts *trackCuts) { fTrackCuts = trackCuts; }
	virtual TObject *GetTrackCuts() { return fTrackCuts; }

private:
	AliESDtrackCuts *fTrackCuts;				// Track cuts
	UInt_t			fFilterBits;				// Track filter bits

	ClassDef(AliEMCalPtTaskTrackSelectionAOD, 1);		// Track selection class for AOD analysis

};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALPTTRACKSELECTIONAOD_H_ */
