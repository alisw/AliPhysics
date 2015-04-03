/**
 * \file AliEMCalPtTaskVTrackSelection.h
 * \brief Declartion of class AliEMCalPtTaskVTrackSelection
 *
 * In this header file the class AliEMCalPtTaskVTrackSelection, which handles the
 * track selection in a transparent way for ESD, AOD and Pico Tracks, is declared.
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Dec 12, 2014
 */
#ifndef ALIEMCALPTTASKVTRACKSELECTION_H_
#define ALIEMCALPTTASKVTRACKSELECTION_H_
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include <TObject.h>

class TClonesArray;
class TObjArray;
class AliVCuts;
class AliVEvent;
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
 * \class AliEMCalPtTaskVTrackSelection
 * \brief Interface for virtual track selection
 *
 * Interface for track selection for the analysis of charged hadrons in
 * EMCal-triggered events. The following abstract functions need to be implemented
 * by inheriting classes:
 * - GetAcceptedTracks (with TClonesArray and AliVEvent as parameters)
 * - IsTrackAccepted (with AliVTrackCuts)
 */
class AliEMCalPtTaskVTrackSelection : public TObject {
public:
	AliEMCalPtTaskVTrackSelection();
	AliEMCalPtTaskVTrackSelection(const AliEMCalPtTaskVTrackSelection &ref);
	AliEMCalPtTaskVTrackSelection &operator=(const AliEMCalPtTaskVTrackSelection &ref);
	virtual ~AliEMCalPtTaskVTrackSelection();

	virtual TObjArray *GetAcceptedTracks(const TClonesArray * const tracks) = 0;
	virtual TObjArray *GetAcceptedTracks(const AliVEvent *const event) = 0;
	virtual bool IsTrackAccepted(AliVTrack * const trk) = 0;

	void AddTrackCuts(AliVCuts *cuts);
	Int_t GetNumberOfCutObjects() const;
	AliVCuts *GetTrackCuts(Int_t icut);

protected:
	TObjArray *fListOfTracks;		///< TObjArray with accepted tracks
	TObjArray *fListOfCuts;     ///< List of track cut objects

	/// \cond CLASSIMP
	ClassDef(AliEMCalPtTaskVTrackSelection, 1); // Track selection for the EMCal pt analysis
	/// \endcond
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALPTTASKVTRACKSELECTION_H_ */
