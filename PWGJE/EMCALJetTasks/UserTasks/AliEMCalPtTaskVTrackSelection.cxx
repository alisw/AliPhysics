/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/*
 * Interface for track selection for the analysis of charged hadrons in
 * EMCal-triggered events
 *
 * Author:
 * 		Markus Fasel
 */
#include <TObjArray.h>
#include <AliEMCalPtTaskVTrackSelection.h>

ClassImp(EMCalTriggerPtAnalysis::AliEMCalPtTaskVTrackSelection)

namespace EMCalTriggerPtAnalysis {

AliEMCalPtTaskVTrackSelection::AliEMCalPtTaskVTrackSelection() :
	TObject(),
	fListOfTracks(NULL)
{
	/*
	 * Default constructor
	 */
}

AliEMCalPtTaskVTrackSelection::AliEMCalPtTaskVTrackSelection(const AliEMCalPtTaskVTrackSelection& ref):
	TObject(ref),
	fListOfTracks(NULL)
{
	if(ref.fListOfTracks) fListOfTracks = new TObjArray(*(ref.fListOfTracks));
}

AliEMCalPtTaskVTrackSelection& AliEMCalPtTaskVTrackSelection::operator=(const AliEMCalPtTaskVTrackSelection& ref) {
	TObject::operator=(ref);
	if(this != &ref){
		this->~AliEMCalPtTaskVTrackSelection();
		if(ref.fListOfTracks) fListOfTracks = new TObjArray(*(ref.fListOfTracks));
	}
	return *this;
}

AliEMCalPtTaskVTrackSelection::~AliEMCalPtTaskVTrackSelection() {
	if(fListOfTracks) delete fListOfTracks;
}

} /* namespace EMCalTriggerPtAnalysis */
