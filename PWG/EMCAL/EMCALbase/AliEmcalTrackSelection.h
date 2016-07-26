/**
 * \file AliEmcalTrackSelection.h
 * \brief Declartion of class AliEmcalTrackSelection
 *
 * In this header file the class AliEmcalTrackSelection, which handles the
 * track selection in a transparent way for ESD, AOD and Pico Tracks, is declared.
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 * \date Jan 30, 2016
 */
#ifndef ALIEMCALTRACKSELECTION_H_
#define ALIEMCALTRACKSELECTION_H_
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>
#include <TBits.h>

class TClonesArray;
class TObjArray;
class AliVCuts;
class AliVEvent;
class AliVTrack;

/**
 * \class AliEmcalTrackSelection
 * \brief Interface for virtual track selection
 * \ingroup EMCALCOREFW
 *
 * Interface for track selection within the EMCAL framework. Enables transparent track selection
 * for ESDs and AODs by implementing a wrapper derived from this class. The following abstract
 * functions need to be implemented by inheriting classes:
 * - GetAcceptedTracks (with TClonesArray and AliVEvent as parameters)
 * - IsTrackAccepted (with AliVTrackCuts)
 * - GenerateTrackCuts
 */
class AliEmcalTrackSelection : public TObject {
public:
  enum ETrackFilterType_t {
    kNoTrackFilter = 0,
    kCustomTrackFilter,
    kHybridTracks,
    kTPCOnlyTracks
  };

	AliEmcalTrackSelection();
	AliEmcalTrackSelection(const AliEmcalTrackSelection &ref);
	AliEmcalTrackSelection &operator=(const AliEmcalTrackSelection &ref);
	virtual ~AliEmcalTrackSelection();

	TObjArray *GetAcceptedTracks(const TClonesArray * const tracks);
	TObjArray *GetAcceptedTracks(const AliVEvent *const event);
	virtual bool IsTrackAccepted(AliVTrack * const trk) = 0;

	virtual void GenerateTrackCuts(ETrackFilterType_t type, const char* period = "") = 0;
	void AddTrackCuts(AliVCuts *cuts);
	void AddTrackCuts(TObjArray *cuts);
	Int_t GetNumberOfCutObjects() const;
	AliVCuts *GetTrackCuts(Int_t icut);

	const TBits& GetTrackBitmap() const { return fTrackBitmap; }
	const TClonesArray* GetAcceptedTrackBitmaps() const { return fListOfTrackBitmaps; }

	void SetSelectionModeAny() { fSelectionModeAny = kTRUE ; }
	void SetSelectionModeAll() { fSelectionModeAny = kFALSE; }

protected:
	TObjArray    *fListOfTracks;         ///< TObjArray with accepted tracks
	TClonesArray *fListOfTrackBitmaps;   ///< TClonesArray with accepted tracks' bit maps
	TBits         fTrackBitmap;          ///< Bitmap of last accepted/rejected track
	TObjArray    *fListOfCuts;           ///< List of track cut objects
	Bool_t        fSelectionModeAny;     ///< Accept track if any of the cuts is fulfilled

	/// \cond CLASSIMP

	/// \endcond
};

#endif /* ALIEMCALTRACKSELECTION_H_ */
