/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors:                                                       *
 *   Artur Szostak <artursz@iafrica.com>                                  *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// $Id: $

///
/// @file   AliHLTMUONDecision.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   12 May 2008
/// @brief  Implementation of the AliHLTMUONDecision class.
///
/// The class is used to store the dHLT decision in a more convenient ROOT
/// object format for testing and analysis.
///

#include "AliHLTMUONDecision.h"
#include "AliHLTMUONTrack.h"
#include "AliHLTMUONMansoTrack.h"
#include "AliLog.h"
#include <cstring>
#include <iostream>
#include <iomanip>

ClassImp(AliHLTMUONDecision);
ClassImp(AliHLTMUONDecision::AliTrackDecision);
ClassImp(AliHLTMUONDecision::AliPairDecision);


std::ostream& operator << (
		std::ostream& stream,
		const AliHLTMUONDecision& decision
	)
{
/// Stream operator for std::ostream classes.
/// \param stream  The output stream object being written to.
/// \param track  The dHLT decision object to print to the stream.
/// \returns  Returns 'stream'.

	stream	<< "No. low/high pT: [" << decision.fNlowPt
		<< ", " << decision.fNhighPt
		<< "]; No. any/low/high unlike: [" << decision.fNunlikeAnyPt
		<< ", " << decision.fNunlikeLowPt
		<< ", " << decision.fNunlikeHighPt
		<< "]; No. any/low/high like = " << decision.fNlikeAnyPt
		<< ", " << decision.fNlikeLowPt
		<< ", " << decision.fNlikeHighPt
		<< "]; No. any/low/high mass = " << decision.fNmassAny
		<< ", " << decision.fNmassLow
		<< ", " << decision.fNmassHigh
		<< "]";
	return stream;
}


AliHLTMUONDecision::AliHLTMUONDecision(
		UInt_t nLowPt, UInt_t nHiPt,
		UInt_t nUnlikeAnyPt, UInt_t nUnlikeLowPt, UInt_t nUnlikeHighPt,
		UInt_t nLikeAnyPt, UInt_t nLikeLowPt, UInt_t nLikeHighPt,
		UInt_t nMassAny, UInt_t nMassLow, UInt_t nMassHigh
	) :
	TObject(),
	fNlowPt(nLowPt),
	fNhighPt(nHiPt),
	fNunlikeAnyPt(nUnlikeAnyPt),
	fNunlikeLowPt(nUnlikeLowPt),
	fNunlikeHighPt(nUnlikeHighPt),
	fNlikeAnyPt(nLikeAnyPt),
	fNlikeLowPt(nLikeLowPt),
	fNlikeHighPt(nLikeHighPt),
	fNmassAny(nMassAny),
	fNmassLow(nMassLow),
	fNmassHigh(nMassHigh),
	fTrackDecisions("AliHLTMUONDecision::AliTrackDecision"),
	fPairDecisions("AliHLTMUONDecision::AliPairDecision")
{
/// Constructor for creating a dHLT decision object.
/// \param nLowPt  Number of tracks above low pT cut.
/// \param nHiPt   Number of tracks above high pT cut.
/// \param nUnlikeAnyPt  Number of track pairs with unlike sign.
/// \param nUnlikeLowPt  Number of unlike sign track pairs with pT > low cut.
/// \param nUnlikeHighPt Number of unlike sign track pairs with pT > high cut.
/// \param nLikeAnyPt   Number of track pairs with like sign.
/// \param nLikeLowPt   Number of like sign track pairs with pT > low cut.
/// \param nLikeHighPt  Number of like sign track pairs with pT > high cut.
/// \param nMassAny   Number of unlike sign track pairs with invariant mass > low cut.
/// \param nMassLow   Number of unlike sign track pairs with invariant mass > low mass cut and pT > low pT cut.
/// \param nMassHigh  Number of unlike sign track pairs with invariant mass > high mass cut and pT > high pT cut.
}


void AliHLTMUONDecision::AddDecision(
		Float_t pt, Bool_t passedLowCut, Bool_t passedHighCut,
		const TObject* track
	)
{
/// Add a single track decision to the dHLT trigger.
/// \param pt  The calculated pT value used for the trigger decision.
/// \param passedLowCut  Flag indicating if the track passed the low pT cut.
/// \param passedHighCut  Flag indicating if the track passed the high pT cut.
/// \param track  Pointer to the associated track object.

	new (fTrackDecisions[fTrackDecisions.GetEntriesFast()])
		AliTrackDecision(
			pt, passedLowCut, passedHighCut, track
		);
}


void AliHLTMUONDecision::AddDecision(const AliTrackDecision* decision)
{
/// Add a single track decision to the dHLT trigger.

	new (fTrackDecisions[fTrackDecisions.GetEntriesFast()])
		AliTrackDecision(*decision);
}


void AliHLTMUONDecision::AddDecision(
		Float_t mass, Bool_t passedLowCut,
		Bool_t passedHighCut, Bool_t unlike,
		UChar_t lowPtCount, UChar_t highPtCount,
		const TObject* trackA, const TObject* trackB
	)
{
/// Add a track pair decision to the dHLT trigger.
/// \param mass   The invariant mass of the track pair.
/// \param passedLowCut  Indicates if the pair passed the low mass cut.
/// \param passedHighCut  Indicates if the pair passed the high mass cut.
/// \param unlike   Indicates if the tracks have opposite sign.
/// \param lowPtCount  The number of tracks in the pair that passed the low pT cut.
///            Should be in the range [0..2].
/// \param highPtCount  The number of tracks in the pair that passed the high pT cut.
///            Should be in the range [0..2].
/// \param trackA  Pointer to the first associated track object.
/// \param trackB  Pointer to the second associated track object.

	new (fPairDecisions[fPairDecisions.GetEntriesFast()])
		AliPairDecision(
			mass, passedLowCut, passedHighCut, unlike,
			lowPtCount, highPtCount, trackA, trackB
		);
}


void AliHLTMUONDecision::AddDecision(const AliPairDecision* decision)
{
/// Add a track pair decision to the dHLT trigger.

	new (fPairDecisions[fPairDecisions.GetEntriesFast()])
		AliPairDecision(*decision);
}


void AliHLTMUONDecision::Print(Option_t* option) const
{
/// Prints the trigger decision information to standard output (screen).
/// \param option  Can be one of the following:
///      - "compact" - prints in a compact format.
///      - "detail" - prints trigger information in a more detailed format.
///      - "all" - prints a full dump of the trigger object.

	using namespace std;
	
	if (option == NULL or strcmp(option, "") == 0 or
	    strcmp(option, "compact") == 0
	   )
	{
		cout << *this << endl;
	}
	else if (strcmp(option, "detail") == 0)
	{
		cout << "dHLT trigger decision scalars:" << endl;
		cout << "No. tracks passed pT cut," << endl;
		cout << "  low = " << fNlowPt << endl;
		cout << " high = " << fNhighPt << endl;
		cout << "No. unlike sign pairs," << endl;
		cout << "          total = " << fNunlikeAnyPt << endl;
		cout << "   pT > low cut = " << fNunlikeLowPt << endl;
		cout << "  pT > high cut = " << fNunlikeHighPt << endl;
		cout << "No. like sign pairs," << endl;
		cout << "          total = " << fNlikeAnyPt << endl;
		cout << "   pT > low cut = " << fNlikeLowPt << endl;
		cout << "  pT > high cut = " << fNlikeHighPt << endl;
		cout << "No. pairs with," << endl;
		cout << "          invariant mass > low cut = " << fNmassAny << endl;
		cout << "   invariant mass and pT > low cut = " << fNmassLow << endl;
		cout << "  invariant mass and pT > high cut = " << fNmassHigh << endl;
		
		streamsize w = cout.width();
		ios::fmtflags f = cout.flags();
		cout << "Triggers for single tracks:" << endl;
		cout << setw(10) << "" << setw(12) << "pT  " << setw(15) << "Passed pT cut" << endl;
		cout	<< setw(10) << "Track" << setw(12) << "(GeV/c)" << setw(6) << "low "
			<< setw(3) << " | " << setw(6) << "high" << endl;
		for (Int_t i = 0; i < NumberOfTracks(); i++)
		{
			const AliTrackDecision* decision = SingleTrackDecision(i);
			if (decision == NULL) continue;
			
			if (decision->MansoTrack() != NULL)
			{
				cout << setw(10) << decision->MansoTrack()->Id();
			}
			else if (decision->FullTrack() != NULL)
			{
				cout << setw(10) << decision->FullTrack()->Id();
			}
			else
			{
				cout << setw(10) << "-";
			}
			
			cout	<< setw(12) << decision->Pt()
				<< setw(6) << (decision->PassedLowPtCut() ? "yes" : "no")
				<< setw(3) << "   "
				<< setw(6) << (decision->PassedHighPtCut() ? "yes" : "no")
				<< endl;
		}
		cout << "Triggers for track pairs:" << endl;
		cout	<< setw(20) << "Track pair" << setw(6) << "Like" << setw(12) << "mass  "
			<< setw(17) << "Passed mass cut"
			<< setw(20) << "No. with pT > than" << endl;
		cout	<< setw(10) << "track A" << setw(10) << "track B"
			<< setw(6) << "sign" << setw(12) << "(GeV/c^2)"
			<< setw(8) << "low " << setw(3) << " | " << setw(6) << "high"
			<< setw(11) << "low " << setw(3) << " | " << setw(6) << "high" << endl;
		for (Int_t j = 0; j < NumberOfPairs(); j++)
		{
			const AliPairDecision* decision = TrackPairDecision(j);
			if (decision == NULL) continue;
			
			if (decision->MansoTrackA() != NULL)
			{
				cout << setw(10) << decision->MansoTrackA()->Id();
			}
			else if (decision->FullTrackA() != NULL)
			{
				cout << setw(10) << decision->FullTrackA()->Id();
			}
			else
			{
				cout << setw(10) << "-";
			}
			
			if (decision->MansoTrackB() != NULL)
			{
				cout << setw(10) << decision->MansoTrackB()->Id();
			}
			else if (decision->FullTrackB() != NULL)
			{
				cout << setw(10) << decision->FullTrackB()->Id();
			}
			else
			{
				cout << setw(10) << "-";
			}
			
			cout	<< setw(6) << (decision->LikeSign() ? "yes" : "no")
				<< setw(12) << decision->Mass()
				<< setw(8) << (decision->PassedLowMassCut() ? "yes" : "no")
				<< setw(3) << "   "
				<< setw(6) << (decision->PassedHighMassCut() ? "yes" : "no")
				<< setw(11) << Int_t(decision->NumberPassedLowPtCut())
				<< setw(3) << "   "
				<< setw(6) << Int_t(decision->NumberPassedHighPtCut())
				<< endl;
		}
		cout.width(w); // reset the field width to previous value.
		cout.flags(f); // reset the flags to previous values.
	}
	else if (strcmp(option, "all") == 0)
	{
		cout << "dHLT trigger decision scalars:" << endl;
		cout << "No. tracks passed pT cut," << endl;
		cout << "  low = " << fNlowPt << endl;
		cout << " high = " << fNhighPt << endl;
		cout << "No. unlike sign pairs," << endl;
		cout << "          total = " << fNunlikeAnyPt << endl;
		cout << "   pT > low cut = " << fNunlikeLowPt << endl;
		cout << "  pT > high cut = " << fNunlikeHighPt << endl;
		cout << "No. like sign pairs," << endl;
		cout << "          total = " << fNlikeAnyPt << endl;
		cout << "   pT > low cut = " << fNlikeLowPt << endl;
		cout << "  pT > high cut = " << fNlikeHighPt << endl;
		cout << "No. pairs with," << endl;
		cout << "          invariant mass > low cut = " << fNmassAny << endl;
		cout << "   invariant mass and pT > low cut = " << fNmassLow << endl;
		cout << "  invariant mass and pT > high cut = " << fNmassHigh << endl;
		
		cout << "===============================================" << endl;
		cout << "========== Triggers for single tracks: ========" << endl;
		for (Int_t i = 0; i < NumberOfTracks(); i++)
		{
			const AliTrackDecision* decision = SingleTrackDecision(i);
			if (decision == NULL) continue;
			decision->Print("all");
		}
		if (NumberOfTracks() == 0) cout << "(None)" << endl;
		cout << "===============================================" << endl;
		cout << "========== Triggers for track pairs: ==========" << endl;
		for (Int_t j = 0; j < NumberOfPairs(); j++)
		{
			const AliPairDecision* decision = TrackPairDecision(j);
			if (decision == NULL) continue;
			decision->Print("all");
		}
		if (NumberOfPairs() == 0) cout << "(None)" << endl;
	}
	else
	{
		AliError("Unknown option specified. Can only be one of 'compact',"
			" 'detail' or 'all'."
		);
	}
}


Int_t AliHLTMUONDecision::Compare(const TObject* obj) const
{
/// We compare this object with 'obj' first by the trigger scalars and then
/// by the signle track and track pair decision lists.
/// \param obj  This is the object to compare to. It must be of type AliHLTMUONDecision.
/// \returns  -1 if 'this' is smaller than 'obj', 1 if greater and zero if both
///      objects are the same.

	if (obj->IsA() == AliHLTMUONDecision::Class())
	{
		const AliHLTMUONDecision* d =
			static_cast<const AliHLTMUONDecision*>(obj);
		if (fNlowPt < d->fNlowPt) return -1;
		if (fNlowPt > d->fNlowPt) return 1;
		if (fNhighPt < d->fNhighPt) return -1;
		if (fNhighPt > d->fNhighPt) return 1;
		if (fNunlikeAnyPt < d->fNunlikeAnyPt) return -1;
		if (fNunlikeAnyPt > d->fNunlikeAnyPt) return 1;
		if (fNunlikeLowPt < d->fNunlikeLowPt) return -1;
		if (fNunlikeLowPt > d->fNunlikeLowPt) return 1;
		if (fNunlikeHighPt < d->fNunlikeHighPt) return -1;
		if (fNunlikeHighPt > d->fNunlikeHighPt) return 1;
		if (fNlikeAnyPt < d->fNlikeAnyPt) return -1;
		if (fNlikeAnyPt > d->fNlikeAnyPt) return 1;
		if (fNlikeLowPt < d->fNlikeLowPt) return -1;
		if (fNlikeLowPt > d->fNlikeLowPt) return 1;
		if (fNlikeHighPt < d->fNlikeHighPt) return -1;
		if (fNlikeHighPt > d->fNlikeHighPt) return 1;
		if (fNmassAny < d->fNmassAny) return -1;
		if (fNmassAny > d->fNmassAny) return 1;
		if (fNmassLow < d->fNmassLow) return -1;
		if (fNmassLow > d->fNmassLow) return 1;
		if (fNmassHigh < d->fNmassHigh) return -1;
		if (fNmassHigh > d->fNmassHigh) return 1;
		
		// Now check the track decision arrays.
		if (NumberOfTracks() < d->NumberOfTracks()) return -1;
		if (NumberOfTracks() > d->NumberOfTracks()) return 1;
		for (Int_t i = 0; i < NumberOfTracks(); i++)
		{
			Int_t result = SingleTrackDecision(i)->Compare( d->SingleTrackDecision(i) );
			if (result != 0) return result;
		}
		if (NumberOfPairs() < d->NumberOfPairs()) return -1;
		if (NumberOfPairs() > d->NumberOfPairs()) return 1;
		for (Int_t j = 0; j < NumberOfPairs(); j++)
		{
			Int_t result = TrackPairDecision(j)->Compare( d->TrackPairDecision(j) );
			if (result != 0) return result;
		}
		
		// At this point everything was equal so return 0 to indicate this fact.
		return 0;
	}
	else
	{
		AliError(Form("Do not know how to compare %s to %s.",
			this->ClassName(),
			obj->ClassName()
		));
		return -999;
	}
}


bool AliHLTMUONDecision::operator == (const AliHLTMUONDecision& decision) const
{
/// Comparison operator just compares if the scalars are the same.
/// \param decision  The trigger decision object to compare to.
/// \returns  true if 'this' object has the same scalars as 'decision', else false.

	return	fNlowPt == decision.fNlowPt
		and fNhighPt == decision.fNhighPt
		and fNunlikeAnyPt == decision.fNunlikeAnyPt
		and fNunlikeLowPt == decision.fNunlikeLowPt
		and fNunlikeHighPt == decision.fNunlikeHighPt
		and fNlikeAnyPt == decision.fNlikeAnyPt
		and fNlikeLowPt == decision.fNlikeLowPt
		and fNlikeHighPt == decision.fNlikeHighPt
		and fNmassAny == decision.fNmassAny
		and fNmassLow == decision.fNmassLow
		and fNmassHigh == decision.fNmassHigh;
}


std::ostream& operator << (
		std::ostream& stream,
		const AliHLTMUONDecision::AliTrackDecision& decision
	)
{
/// Stream operator for std::ostream classes.
/// \param stream  The output stream object being written to.
/// \param track  The dHLT decision object to print to the stream.
/// \returns  Returns 'stream'.

	stream	<< "Passed low/high pT cut: [" << (decision.fPassedLowCut ? "yes" : "no")
		<< ", " << (decision.fPassedHighCut ? "yes" : "no")
		<< "]; with pT = " << decision.fPt;
	return stream;
}


const AliHLTMUONMansoTrack* AliHLTMUONDecision::AliTrackDecision::MansoTrack() const
{
/// Returns the associated track as a Manso track object and NULL if the track
/// object is missing or not a Manso track object.

	if (fTrack == NULL) return NULL;
	return dynamic_cast<const AliHLTMUONMansoTrack*>(fTrack);
}


const AliHLTMUONTrack* AliHLTMUONDecision::AliTrackDecision::FullTrack() const
{
/// Returns the associated track as a full track object and NULL if the track
/// object is missing or not a full track object.

	if (fTrack == NULL) return NULL;
	return dynamic_cast<const AliHLTMUONTrack*>(fTrack);
}


void AliHLTMUONDecision::AliTrackDecision::Print(Option_t* option) const
{
/// Prints the trigger decision to standard output (screen).
/// \param option  Can be one of the following:
///      - "compact" - prints in a compact format.
///      - "detail" - prints trigger information in a more detailed format.
///      - "all" - prints a full dump of the trigger object.

	using namespace std;
	
	if (option == NULL or strcmp(option, "") == 0 or
	    strcmp(option, "compact") == 0
	   )
	{
		cout << *this << endl;
	}
	else if (strcmp(option, "detail") == 0)
	{
		Int_t id = -1;
		if (MansoTrack() != NULL) id = MansoTrack()->Id();
		else if (FullTrack() != NULL) id = FullTrack()->Id();
		cout << "Trigger decision for track: " << id << endl;
		cout << "pT = " << fPt << " GeV/c" << endl;
		cout << "pT cut | passed" << endl;
		cout << "-------+--------" << endl;
		cout << "  low  |  " << (fPassedLowCut ? "yes" : "no") << endl;
		cout << " high  |  " << (fPassedHighCut ? "yes" : "no") << endl;
	}
	else if (strcmp(option, "all") == 0)
	{
		Int_t id = -1;
		if (MansoTrack() != NULL) id = MansoTrack()->Id();
		else if (FullTrack() != NULL) id = FullTrack()->Id();
		cout << "Trigger decision for track: " << id << endl;
		cout << "pT = " << fPt << " GeV/c" << endl;
		cout << "pT cut | passed" << endl;
		cout << "-------+--------" << endl;
		cout << "  low  |  " << (fPassedLowCut ? "yes" : "no") << endl;
		cout << " high  |  " << (fPassedHighCut ? "yes" : "no") << endl;
		cout << "===== Track details =====" << endl;
		fTrack->Print("all");
	}
	else
	{
		AliError("Unknown option specified. Can only be one of 'compact',"
			" 'detail' or 'all'."
		);
	}
}


Int_t AliHLTMUONDecision::AliTrackDecision::Compare(const TObject* obj) const
{
/// We compare this object with 'obj' field by field.
/// \param obj  This is the object to compare to. It must be of type
///      AliHLTMUONDecision::AliTrackDecision.
/// \returns  -1 if 'this' is smaller than 'obj', 1 if greater and zero if both
///      objects are the same.

	if (obj->IsA() == AliHLTMUONDecision::Class())
	{
		const AliHLTMUONDecision::AliTrackDecision* d =
			static_cast<const AliHLTMUONDecision::AliTrackDecision*>(obj);
		if (fPt < d->fPt) return -1;
		if (fPt > d->fPt) return 1;
		if (fPassedLowCut < d->fPassedLowCut) return -1;
		if (fPassedLowCut > d->fPassedLowCut) return 1;
		if (fPassedHighCut < d->fPassedHighCut) return -1;
		if (fPassedHighCut > d->fPassedHighCut) return 1;
		return fTrack->Compare(d->fTrack);
	}
	else
	{
		AliError(Form("Do not know how to compare %s to %s.",
			this->ClassName(),
			obj->ClassName()
		));
		return -999;
	}
}


std::ostream& operator << (
		std::ostream& stream,
		const AliHLTMUONDecision::AliPairDecision& decision
	)
{
/// Stream operator for std::ostream classes.
/// \param stream  The output stream object being written to.
/// \param track  The dHLT decision object to print to the stream.
/// \returns  Returns 'stream'.

	stream	<< (decision.fUnlike ? "Unlike" : "Like")
		<< " sign pair passed low/high mass cut: [" << (decision.fPassedLowCut ? "yes" : "no")
		<< ", " << (decision.fPassedHighCut ? "yes" : "no")
		<< "]; with mass = " << decision.fMass;
	return stream;
}


const AliHLTMUONMansoTrack* AliHLTMUONDecision::AliPairDecision::MansoTrackA() const
{
/// Returns the first associated track as a Manso track object and NULL if the
/// track object is missing or not a Manso track object.

	if (fTrackA == NULL) return NULL;
	return dynamic_cast<const AliHLTMUONMansoTrack*>(fTrackA);
}


const AliHLTMUONTrack* AliHLTMUONDecision::AliPairDecision::FullTrackA() const
{
/// Returns the first associated track as a full track object and NULL if the
/// track object is missing or not a full track object.

	if (fTrackA == NULL) return NULL;
	return dynamic_cast<const AliHLTMUONTrack*>(fTrackA);
}


const AliHLTMUONMansoTrack* AliHLTMUONDecision::AliPairDecision::MansoTrackB() const
{
/// Returns the second associated track as a Manso track object and NULL if the
/// track object is missing or not a Manso track object.

	if (fTrackB == NULL) return NULL;
	return dynamic_cast<const AliHLTMUONMansoTrack*>(fTrackB);
}


const AliHLTMUONTrack* AliHLTMUONDecision::AliPairDecision::FullTrackB() const
{
/// Returns the second associated track as a full track object and NULL if the
/// track object is missing or not a full track object.

	if (fTrackB == NULL) return NULL;
	return dynamic_cast<const AliHLTMUONTrack*>(fTrackB);
}


void AliHLTMUONDecision::AliPairDecision::Print(Option_t* option) const
{
/// Prints the trigger decision to standard output (screen).
/// \param option  Can be one of the following:
///      - "compact" - prints in a compact format.
///      - "detail" - prints trigger information in a more detailed format.
///      - "all" - prints a full dump of the trigger object.

	using namespace std;
	
	if (option == NULL or strcmp(option, "") == 0 or
	    strcmp(option, "compact") == 0
	   )
	{
		cout << *this << endl;
	}
	else if (strcmp(option, "detail") == 0)
	{
		Int_t id1 = -1;
		if (MansoTrackA() != NULL) id1 = MansoTrackA()->Id();
		else if (FullTrackA() != NULL) id1 = FullTrackA()->Id();
		Int_t id2 = -1;
		if (MansoTrackB() != NULL) id2 = MansoTrackB()->Id();
		else if (FullTrackB() != NULL) id2 = FullTrackB()->Id();
		cout	<< "Trigger decision for track pair: {" << id1
			<< ", " << id2 << "}" << endl;
		cout << "Invariant mass = " << fMass << " GeV/c^2" << endl;
		cout << "mass cut | passed" << endl;
		cout << "---------+--------" << endl;
		cout << "   low   |  " << (fPassedLowCut ? "yes" : "no") << endl;
		cout << "  high   |  " << (fPassedHighCut ? "yes" : "no") << endl;
		cout << "Number of tracks in pair that passed," << endl;
		cout << "   low pT cut = " << Int_t(fLowPtCount) << endl;
		cout << "  high pT cut = " << Int_t(fHighPtCount) << endl;
	}
	else if (strcmp(option, "all") == 0)
	{
		Int_t id1 = -1;
		if (MansoTrackA() != NULL) id1 = MansoTrackA()->Id();
		else if (FullTrackA() != NULL) id1 = FullTrackA()->Id();
		Int_t id2 = -1;
		if (MansoTrackB() != NULL) id2 = MansoTrackB()->Id();
		else if (FullTrackB() != NULL) id2 = FullTrackB()->Id();
		cout	<< "Trigger decision for track pair: {" << id1
			<< ", " << id2 << "}" << endl;
		cout << "Invariant mass = " << fMass << " GeV/c^2" << endl;
		cout << "mass cut | passed" << endl;
		cout << "---------+--------" << endl;
		cout << "   low   |  " << (fPassedLowCut ? "yes" : "no") << endl;
		cout << "  high   |  " << (fPassedHighCut ? "yes" : "no") << endl;
		cout << "Number of tracks in pair that passed," << endl;
		cout << "   low pT cut = " << Int_t(fLowPtCount) << endl;
		cout << "  high pT cut = " << Int_t(fHighPtCount) << endl;
		cout << "===== First track details =====" << endl;
		fTrackA->Print("all");
		cout << "===== Second track details =====" << endl;
		fTrackB->Print("all");
	}
	else
	{
		AliError("Unknown option specified. Can only be one of 'compact',"
			" 'detail' or 'all'."
		);
	}
}


Int_t AliHLTMUONDecision::AliPairDecision::Compare(const TObject* obj) const
{
/// We compare this object with 'obj' field by field.
/// \param obj  This is the object to compare to. It must be of type
///      AliHLTMUONDecision::AliPairDecision.
/// \returns  -1 if 'this' is smaller than 'obj', 1 if greater and zero if both
///      objects are the same.

	if (obj->IsA() == AliHLTMUONDecision::Class())
	{
		const AliHLTMUONDecision::AliPairDecision* d =
			static_cast<const AliHLTMUONDecision::AliPairDecision*>(obj);
		if (fMass < d->fMass) return -1;
		if (fMass > d->fMass) return 1;
		if (fPassedLowCut < d->fPassedLowCut) return -1;
		if (fPassedLowCut > d->fPassedLowCut) return 1;
		if (fPassedHighCut < d->fPassedHighCut) return -1;
		if (fPassedHighCut > d->fPassedHighCut) return 1;
		if (fUnlike < d->fUnlike) return -1;
		if (fUnlike > d->fUnlike) return 1;
		if (fLowPtCount < d->fLowPtCount) return -1;
		if (fLowPtCount > d->fLowPtCount) return 1;
		if (fHighPtCount < d->fHighPtCount) return -1;
		if (fHighPtCount > d->fHighPtCount) return 1;
		Int_t result = fTrackA->Compare(d->fTrackA);
		if (result != 0) return result;
		return fTrackB->Compare(d->fTrackB);
	}
	else
	{
		AliError(Form("Do not know how to compare %s to %s.",
			this->ClassName(),
			obj->ClassName()
		));
		return -999;
	}
}

