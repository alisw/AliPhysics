/// \class AliESDTrackSelection
/// \brief Declaration of class AliESDTrackSelection
///
/// In this header file the class AliESDTrackSelection, which implements
/// the virtual track selection for ESD tracks, is declared
///
/// \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
/// \date Jul 24, 2015

#ifndef ALIESDTASKTRACKSELECTION_H_
#define ALIESDTASKTRACKSELECTION_H_
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <AliVTrackSelection.h>

class AliVCuts;
class AliVTrack;

class AliESDTrackSelection: public AliVTrackSelection {
public:
	AliESDTrackSelection();
	AliESDTrackSelection(AliVCuts *cuts);
	virtual ~AliESDTrackSelection() {}

	virtual bool IsTrackAccepted(AliVTrack * const trk);

	ClassDef(AliESDTrackSelection,1);
};

#endif /* ALIESDTRACKSELECTION_H_ */
