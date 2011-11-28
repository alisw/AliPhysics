// -*- Mode: C++ -*-
#ifndef ALIHLTMUONDECISION_H
#define ALIHLTMUONDECISION_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

// $Id$

///
/// @file   AliHLTMUONDecision.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   12 May 2008
/// @brief  Declaration of a dHLT decision object in ROOT object format.
///

#include "TObject.h"
#include "TClonesArray.h"

class AliHLTMUONTrack;
class AliHLTMUONMansoTrack;

/**
 * AliHLTMUONDecision stores converted dHLT raw trigger decision data as a ROOT object.
 * Both AliHLTMUONSinglesDecisionBlockStruct and AliHLTMUONPairsDecisionBlockStruct
 * data blocks are converted into this class by the AliHLTMUONRootifierComponent.
 * This class is mainly for testing or as a helper object for dHLT specific analysis,
 * since it is sometimes easier to store and handle ROOT objects.
 */
class AliHLTMUONDecision : public TObject
{
	/**
	 * Stream operator for usage with std::ostream classes.
	 * Allows usage such as:
	 *   AliHLTMUONDecision t; std::cout << t;
	 */
	friend std::ostream& operator << (
			std::ostream& stream,
			const AliHLTMUONDecision& decision
		);

public:

	/**
	 * The AliTrackDecision class stores per track trigger information.
	 */
	class AliTrackDecision : public TObject
	{
		/**
		 * Stream operator for usage with std::ostream classes.
		 */
		friend std::ostream& operator << (std::ostream& stream, const AliTrackDecision& decision);
	
	public:
		
		/**
		 * Constructor for new single track trigger decision object.
		 * \param pt  The calculated pT value used for the trigger decision.
		 * \param passedLowCut  Flag indicating if the track passed the low pT cut.
		 * \param passedHighCut  Flag indicating if the track passed the high pT cut.
		 * \param track  Pointer to the associated track object.
		 */
		AliTrackDecision(
				Float_t pt = -1,
				Bool_t passedLowCut = kFALSE,
				Bool_t passedHighCut = kFALSE,
				const TObject* track = NULL
			) :
			TObject(), fTrack(track), fPt(pt),
			fPassedLowCut(passedLowCut), fPassedHighCut(passedHighCut)
		{}
		
		/**
		 * Copy constructor performs shallow copy of object since we
		 * do not take ownership of the track object.
		 */
		AliTrackDecision(const AliTrackDecision& obj) :
			TObject(obj), fTrack(obj.fTrack), fPt(obj.fPt),
			fPassedLowCut(obj.fPassedLowCut), fPassedHighCut(obj.fPassedHighCut)
		{}
		
		/**
		 * Asignment operators performs shallow copy of object since we
		 * do not take ownership of the track object.
		 */
		AliTrackDecision& operator = (const AliTrackDecision& obj)
		{
		        if (this==&obj) return *this;
			TObject::operator = (obj);
			fTrack = obj.fTrack; fPt = obj.fPt;
			fPassedLowCut = obj.fPassedLowCut; fPassedHighCut = obj.fPassedHighCut;
			return *this;
		}
		
		/**
		 * Default destructor.
		 */
		virtual ~AliTrackDecision() {}
		
		/**
		 * Returns the track associated with the trigger decision or NULL if none found.
		 */
		const TObject* Track() const { return fTrack; }
		
		/**
		 * Returns the track associated with the trigger decision as a Manso track object.
		 * NULL is returned if no track is found or the track object is not a Manso track.
		 */
		const AliHLTMUONMansoTrack* MansoTrack() const;
		
		/**
		 * Returns the track associated with the trigger decision as a full track object.
		 * NULL is returned if no track is found or the track object is not a full track.
		 */
		const AliHLTMUONTrack* FullTrack() const;
		
		/**
		 * Returns the calculated pT value used for the trigger decision.
		 */
		Float_t Pt() const { return fPt; }
		
		/**
		 * Returns kTRUE if the track passed the low pT cut, else kFALSE.
		 */
		Bool_t PassedLowPtCut() const { return fPassedLowCut; }
		
		/**
		 * Returns kTRUE if the track passed the high pT cut, else kFALSE.
		 */
		Bool_t PassedHighPtCut() const { return fPassedHighCut; }
		
		/// Print method inherited from TObject.
		virtual void Print(Option_t* option = NULL) const;
	
		// Methods inherited from TObject
		virtual Bool_t IsSortable() const { return kTRUE; }
		Int_t Compare(const TObject* obj) const;

		// Implement comparison operators.
		bool operator == (const AliTrackDecision& d) const
		{
			return fTrack == d.fTrack and fPt == d.fPt
				and fPassedLowCut == d.fPassedLowCut
				and fPassedHighCut == d.fPassedHighCut;
		}

		bool operator != (const AliTrackDecision& d) const
		{
			return not this->operator == (d);
		}
	
	private:
		
		const TObject* fTrack; ///< Track associated with this decision.
		Float_t fPt; ///< Calculated pT value used for decision (GeV/c).
		Bool_t fPassedLowCut; ///< Indicates if the track passed the low pT cut.
		Bool_t fPassedHighCut; ///< Indicates if the track passed the high pT cut.
		
		ClassDef(AliHLTMUONDecision::AliTrackDecision, 4); // A single track dHLT trigger decision object.
	};
	
	/**
	 * The AliPairDecision class stores trigger information about a track pair.
	 */
	class AliPairDecision : public TObject
	{
		/**
		 * Stream operator for usage with std::ostream classes.
		 */
		friend std::ostream& operator << (std::ostream& stream, const AliPairDecision& decision);
	
	public:
		
		/**
		 * Constructor for new track pair trigger decision object.
		 * \param mass   The invariant mass of the track pair.
		 * \param passedLowCut  Indicates if the pair passed the low mass cut.
		 * \param passedHighCut  Indicates if the pair passed the high mass cut.
		 * \param unlike   Indicates if the tracks have opposite sign.
		 * \param lowPtCount  The number of tracks in the pair that passed the low pT cut.
		 *            Should be in the range [0..2].
		 * \param highPtCount  The number of tracks in the pair that passed the high pT cut.
		 *            Should be in the range [0..2].
		 * \param trackA  Pointer to the first associated track object.
		 * \param trackB  Pointer to the second associated track object.
		 */
		AliPairDecision(
				Float_t mass = -1, Bool_t passedLowCut = kFALSE,
				Bool_t passedHighCut = kFALSE, Bool_t unlike = kFALSE,
				UChar_t lowPtCount = 0, UChar_t highPtCount = 0,
				const TObject* trackA = NULL,
				const TObject* trackB = NULL
			) :
			TObject(), fTrackA(trackA), fTrackB(trackB), fMass(mass),
			fPassedLowCut(passedLowCut), fPassedHighCut(passedHighCut),
			fUnlike(unlike), fLowPtCount(lowPtCount), fHighPtCount(highPtCount)
		{}
		
		/**
		 * Copy constructor performs shallow copy of object since we
		 * do not take ownership of the track objects.
		 */
		AliPairDecision(const AliPairDecision& obj) :
			TObject(obj), fTrackA(obj.fTrackA), fTrackB(obj.fTrackB), fMass(obj.fMass),
			fPassedLowCut(obj.fPassedLowCut), fPassedHighCut(obj.fPassedHighCut),
			fUnlike(obj.fUnlike), fLowPtCount(obj.fLowPtCount), fHighPtCount(obj.fHighPtCount)
		{}
		
		/**
		 * Asignment operators performs shallow copy of object since we
		 * do not take ownership of the track objects.
		 */
		AliPairDecision& operator = (const AliPairDecision& obj)
		{
		        if (this==&obj) return *this;
			TObject::operator = (obj);
			fTrackA = obj.fTrackA; fTrackB = obj.fTrackB; fMass = obj.fMass; fPassedLowCut = obj.fPassedLowCut;
			fPassedHighCut = obj.fPassedHighCut; fUnlike = obj.fUnlike; fLowPtCount = obj.fLowPtCount; fHighPtCount = obj.fHighPtCount;
			return *this;
		}
		
		/**
		 * Default destructor.
		 */
		virtual ~AliPairDecision() {}
		
		/**
		 * Returns the first track associated with the track pair trigger decision
		 * or NULL if none found.
		 */
		const TObject* TrackA() const { return fTrackA; }
		
		/**
		 * Returns the first track associated with the pair decision as a Manso track object.
		 * NULL is returned if no track is found or the track object is not a Manso track.
		 */
		const AliHLTMUONMansoTrack* MansoTrackA() const;
		
		/**
		 * Returns the first track associated with the pair decision as a full track object.
		 * NULL is returned if no track is found or the track object is not a full track.
		 */
		const AliHLTMUONTrack* FullTrackA() const;
		
		/**
		 * Returns the second track associated with the track pair trigger decision
		 * or NULL if none found.
		 */
		const TObject* TrackB() const { return fTrackB; }
		
		/**
		 * Returns the second track associated with the pair decision as a Manso track object.
		 * NULL is returned if no track is found or the track object is not a Manso track.
		 */
		const AliHLTMUONMansoTrack* MansoTrackB() const;
		
		/**
		 * Returns the second track associated with the pair decision as a full track object.
		 * NULL is returned if no track is found or the track object is not a full track.
		 */
		const AliHLTMUONTrack* FullTrackB() const;
		
		/**
		 * Returns the calculated invariant mass value used for the trigger decision.
		 */
		Float_t Mass() const { return fMass; }
		
		/**
		 * Returns kTRUE if the track pair passed the low invariant mass cut, else kFALSE.
		 */
		Bool_t PassedLowMassCut() const { return fPassedLowCut; }
		
		/**
		 * Returns kTRUE if the track pair passed the high invariant mass cut, else kFALSE.
		 */
		Bool_t PassedHighMassCut() const { return fPassedHighCut; }
		
		/**
		 * Returns kTRUE if the track pair has unlike sign, else kFALSE.
		 */
		Bool_t UnlikeSign() const { return fUnlike; }
		
		/**
		 * Returns kTRUE if the track pair has like sign, else kFALSE.
		 */
		Bool_t LikeSign() const { return not fUnlike; }
		
		/**
		 * Returns the number of tracks in the pair that passed the low pT cut.
		 * Can be one of 0, 1 or 2.
		 */
		UChar_t NumberPassedLowPtCut() const { return fLowPtCount; }
		
		/**
		 * Returns the number of tracks in the pair that passed the high pT cut.
		 * Can be one of 0, 1 or 2.
		 */
		UChar_t NumberPassedHighPtCut() const { return fHighPtCount; }
		
		/**
		 * Returns kTRUE if both tracks passed the low pT cut, else kFALSE.
		 */
		Bool_t BothPassedLowPtCut() const { return NumberPassedLowPtCut() == 2; }
		
		/**
		 * Returns kTRUE if both tracks passed the high pT cut, else kFALSE.
		 */
		Bool_t BothPassedHighPtCut() const { return NumberPassedHighPtCut() == 2; }
		
		/// Print method inherited from TObject.
		virtual void Print(Option_t* option = NULL) const;
	
		// Methods inherited from TObject
		virtual Bool_t IsSortable() const { return kTRUE; }
		Int_t Compare(const TObject* obj) const;

		// Implement comparison operators.
		bool operator == (const AliPairDecision& d) const
		{
			return fTrackA == d.fTrackA and fTrackB == d.fTrackB and fMass == d.fMass
				and fPassedLowCut == d.fPassedLowCut and fPassedHighCut == d.fPassedHighCut
				and fUnlike == d.fUnlike and fLowPtCount == d.fLowPtCount and fHighPtCount == d.fHighPtCount;
		}

		bool operator != (const AliPairDecision& d) const
		{
			return not this->operator == (d);
		}
	
	private:
		
		const TObject* fTrackA; ///< The first track associated with this pair decision.
		const TObject* fTrackB; ///< The second track associated with this pair decision.
		Float_t fMass;  ///< The invariant mass used for the trigger decision. (GeV/c^2)
		Bool_t fPassedLowCut; ///< Indicates if the track passed the low mass cut.
		Bool_t fPassedHighCut; ///< Indicates if the track passed the high mass cut.
		Bool_t fUnlike; ///< Indicates if the track pair has unlike sign.
		UChar_t fLowPtCount; ///< The number of tracks in the pair that passed the low pT cut.
		UChar_t fHighPtCount; ///< The number of tracks in the pair that passed the high pT cut.
		
		ClassDef(AliHLTMUONDecision::AliPairDecision, 4); // A track pair dHLT trigger decision object.
	};

	/**
	 * Constructor for creating a dHLT decision object.
	 * \param nLowPt  Number of tracks above low pT cut.
	 * \param nHiPt   Number of tracks above high pT cut.
	 * \param nUnlikeAnyPt  Number of track pairs with unlike sign.
	 * \param nUnlikeLowPt  Number of unlike sign track pairs with pT > low cut.
	 * \param nUnlikeHighPt Number of unlike sign track pairs with pT > high cut.
	 * \param nLikeAnyPt   Number of track pairs with like sign.
	 * \param nLikeLowPt   Number of like sign track pairs with pT > low cut.
	 * \param nLikeHighPt  Number of like sign track pairs with pT > high cut.
	 * \param nMassAny   Number of unlike sign track pairs with invariant mass > low cut.
	 * \param nMassLow   Number of unlike sign track pairs with invariant mass > low mass cut and pT > low pT cut.
	 * \param nMassHigh  Number of unlike sign track pairs with invariant mass > high mass cut and pT > high pT cut.
	 */
	AliHLTMUONDecision(
			UInt_t nLowPt = 0, UInt_t nHiPt = 0,
			UInt_t nUnlikeAnyPt = 0, UInt_t nUnlikeLowPt = 0, UInt_t nUnlikeHighPt = 0,
			UInt_t nLikeAnyPt = 0, UInt_t nLikeLowPt = 0, UInt_t nLikeHighPt = 0,
			UInt_t nMassAny = 0, UInt_t nMassLow = 0, UInt_t nMassHigh = 0
		);
	
	/**
	 * Default destructor.
	 */
	virtual ~AliHLTMUONDecision() {}

	/**
	 * Returns the number of single low pT triggers.
	 */
	UInt_t NumberOfLowPtTriggers() const { return fNlowPt; }

	/**
	 * Returns the number of single high pT triggers.
	 */
	UInt_t NumberOfHighPtTriggers() const { return fNhighPt; }
	
	/**
	 * Returns the number of unlike sign pt triggers.
	 */
	UInt_t NumberOfUnlikePairs() const { return fNunlikeAnyPt; }
	
	/**
	 * Returns the number of unlike sign pT triggers, where both tracks have pT > low cut.
	 */
	UInt_t NumberOfUnlikeLowPtPairs() const { return fNunlikeLowPt; }
	
	/**
	 * Returns the number of unlike sign pT triggers, where both tracks have pT > high cut.
	 */
	UInt_t NumberOfUnlikeHighPtPairs() const { return fNunlikeHighPt; }
	
	/**
	 * Returns the number of like sign pt triggers.
	 */
	UInt_t NumberOfLikePairs() const { return fNlikeAnyPt; }
	
	/**
	 * Returns the number of like sign pT triggers, where both tracks have pT > low cut.
	 */
	UInt_t NumberOfLikeLowPtPairs() const { return fNlikeLowPt; }
	
	/**
	 * Returns the number of like sign pT triggers, where both tracks have pT > high cut.
	 */
	UInt_t NumberOfLikeHighPtPairs() const { return fNlikeHighPt; }
	
	/**
	 * Returns the number of invariant mass triggers.
	 */
	UInt_t NumberOfMassTriggers() const { return fNmassAny; }
	
	/**
	 * Returns the number of invariant mass triggers,
	 * where invariant mass > low cut and both tracks have pT > low cut.
	 */
	UInt_t NumberOfLowMassTriggers() const { return fNmassLow; }
	
	/**
	 * Returns the number of invariant mass triggers,
	 * where invariant mass > high cut and both tracks have pT > high cut.
	 */
	UInt_t NumberOfHighMassTriggers() const { return fNmassHigh; }
	
	/**
	 * Returns the total number of single tracks that formed part of this trigger decision.
	 */
	Int_t NumberOfTracks() const { return fTrackDecisions.GetEntriesFast(); }
	
	/**
	 * Returns the i'th trigger decision for a single track.
	 * \param i  Should be in the range [0 .. NumberOfTracks()-1]
	 */
	const AliTrackDecision* SingleTrackDecision(Int_t i) const
	{
		return static_cast<const AliTrackDecision*>(fTrackDecisions.At(i));
	}
	
	/**
	 * Returns the total number of track pairs that formed part of this trigger decision.
	 */
	Int_t NumberOfPairs() const { return fPairDecisions.GetEntriesFast(); }
	
	/**
	 * Returns the i'th trigger decision for a track pair.
	 * \param i  Should be in the range [0 .. NumberOfPairs()-1]
	 */
	const AliPairDecision* TrackPairDecision(Int_t i) const
	{
		return static_cast<const AliPairDecision*>(fPairDecisions.At(i));
	}
	
	/// Add a single track decision to the dHLT trigger.
	void AddDecision(const AliTrackDecision* decision);
	
	/// Add a single track decision to the dHLT trigger.
	void AddDecision(
			Float_t pt, Bool_t passedLowCut, Bool_t passedHighCut,
			const TObject* track
		);
	
	/// Add a track pair decision to the dHLT trigger.
	void AddDecision(const AliPairDecision* decision);
	
	/// Add a track pair decision to the dHLT trigger.
	void AddDecision(
			Float_t mass, Bool_t passedLowCut,
			Bool_t passedHighCut, Bool_t unlike,
			UChar_t lowPtCount, UChar_t highPtCount,
			const TObject* trackA, const TObject* trackB
		);
	
	/**
	 * Prints the details of the dHLT decision to screen.
	 * @param option  A case sensitive string that can contain one of the
	 *     following strings:
	 *       "compact" - Prints just the trigger scalars.
	 *       "detail" - Prints also the track decision list in a compact format.
	 *       "all" - Prints all decision information and also the track decision details.
	 *     If the string contains an empty option or NULL then the default is
	 *     to print compactly.
	 */
	virtual void Print(Option_t* option = NULL) const;
	
	// Methods inherited from TObject
	virtual Bool_t IsSortable() const { return kTRUE; }
	Int_t Compare(const TObject* obj) const;

	// Implement comparison operators.
	bool operator == (const AliHLTMUONDecision& track) const;

	bool operator != (const AliHLTMUONDecision& track) const
	{
		return not this->operator == (track);
	}

private:

	// Do not allow copying of this class.
	AliHLTMUONDecision(const AliHLTMUONDecision& track);
	AliHLTMUONDecision& operator = (const AliHLTMUONDecision& track);
	
	UInt_t fNlowPt;  ///< Number of low pt triggers.
	UInt_t fNhighPt; ///< Number of high pt triggers.
	
	/// Number of unlike sign pt triggers for both tracks having any pt.
	UInt_t fNunlikeAnyPt;
	
	/// Number of unlike sign low pt triggers where both tracks have pt > low cut.
	UInt_t fNunlikeLowPt;
	
	/// Number of unlike sign high pt triggers where both tracks have pt > high cut.
	UInt_t fNunlikeHighPt;
	
	/// Number of like sign pt triggers where both tracks have any pt.
	UInt_t fNlikeAnyPt;
	
	/// Number of like sign low pt triggers where both tracks have pt > low cut.
	UInt_t fNlikeLowPt;
	
	/// Number of like sign high pt triggers where both tracks have pt > high cut.
	UInt_t fNlikeHighPt;
	
	/// Number of pairs that have invariant mass > low mass cut, any pt and unlike sign.
	UInt_t fNmassAny;
	
	/// Number of pairs that have invariant mass > low mass cut, pt > low pt cut and unlike sign.
	UInt_t fNmassLow;
	
	/// Number of pairs that have invariant mass > high mass cut, pt > high pt cut and unlike sign.
	UInt_t fNmassHigh;
	
	TClonesArray fTrackDecisions;  ///< Array of single track decision objects.
	TClonesArray fPairDecisions;  ///< Array of track pair decision objects.

	ClassDef(AliHLTMUONDecision, 4); // Decision object containing data converted from raw internal dHLT data structures.
};

#endif // ALIHLTMUONDECISION_H
