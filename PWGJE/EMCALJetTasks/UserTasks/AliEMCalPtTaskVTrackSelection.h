#ifndef ALIEMCALPTTASKVTRACKSELECTION_H_
#define ALIEMCALPTTASKVTRACKSELECTION_H_
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel

#include <TObject.h>

class TClonesArray;
class TObjArray;
class AliVEvent;

namespace EMCalTriggerPtAnalysis {

class AliEMCalPtTaskVTrackSelection : public TObject {
public:
	AliEMCalPtTaskVTrackSelection();
	AliEMCalPtTaskVTrackSelection(const AliEMCalPtTaskVTrackSelection &ref);
	AliEMCalPtTaskVTrackSelection &operator=(const AliEMCalPtTaskVTrackSelection &ref);
	virtual ~AliEMCalPtTaskVTrackSelection();

	virtual TObjArray *GetAcceptedTracks(const TClonesArray * const tracks) = 0;
	virtual TObjArray *GetAcceptedTracks(const AliVEvent *const event) = 0;

	virtual TObject *GetTrackCuts() = 0;

protected:
	TObjArray *fListOfTracks;		// TObjArray with accepted tracks

	ClassDef(AliEMCalPtTaskVTrackSelection, 1); // Track selection for the EMCal pt analysis
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALPTTASKVTRACKSELECTION_H_ */
