/**
 * \file AliEmcalTrackSelection.h
 * \brief Declartion of class AliEmcalTrackSelection
 *
 * In this header file the class AliEmcalTrackSelection, which handles the
 * track selection in a transparent way for ESD, AOD and Pico Tracks, is declared.
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Jul 24, 2015
 */
#ifndef ALIEMCALTRACKSELECTION_H_
#define ALIEMCALTRACKSELECTION_H_
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include <TObject.h>

class TClonesArray;
class TObjArray;
class AliVCuts;
class AliVEvent;
class AliVTrack;

/**
 * \class AliEmcalTrackSelection
 * \brief Interface for virtual track selection
 *
 * Interface for track selection within the EMCAL framework. Enables transparent track selection
 * for ESDs and AODs by implementing a wrapper derived from this class. The following abstract
 * functions need to be implemented by inheriting classes:
 * - GetAcceptedTracks (with TClonesArray and AliVEvent as parameters)
 * - IsTrackAccepted (with AliVTrackCuts)
 */
class AliEmcalTrackSelection : public TObject {
public:
	AliEmcalTrackSelection();
	AliEmcalTrackSelection(const AliEmcalTrackSelection &ref);
	AliEmcalTrackSelection &operator=(const AliEmcalTrackSelection &ref);
	virtual ~AliEmcalTrackSelection();

	virtual TObjArray *GetAcceptedTracks(const TClonesArray * const tracks) = 0;
	virtual TObjArray *GetAcceptedTracks(const AliVEvent *const event) = 0;
	virtual bool IsTrackAccepted(AliVTrack * const trk) = 0;

	void AddTrackCuts(AliVCuts *cuts);
	Int_t GetNumberOfCutObjects() const;
	AliVCuts *GetTrackCuts(Int_t icut);

	void SetSelectionModeAny() { fSelectionModeAny = kTRUE; }
	void SetSelectionModeAll() { fSelectionModeAny = kFALSE; }

protected:
	TObjArray *fListOfTracks;		      ///< TObjArray with accepted tracks
	TObjArray *fListOfCuts;           ///< List of track cut objects
	Bool_t     fSelectionModeAny;     ///< Accept track if any of the cuts is fulfilled

	/// \cond CLASSIMP

	/// \endcond
};

#endif /* ALIEMCALTRACKSELECTION_H_ */
