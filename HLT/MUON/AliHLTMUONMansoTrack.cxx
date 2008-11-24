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

/* $Id$ */

///
/// @file   AliHLTMUONMansoTrack.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   29 Sep 2007
/// @brief  Implementation of the AliHLTMUONMansoTrack class.
///
/// The Manso track class is used to store converted track data from dHLT raw
/// internal data blocks.
///

#include "AliHLTMUONMansoTrack.h"
#include "AliHLTMUONRecHit.h"
#include "AliHLTMUONTriggerRecord.h"
#include "AliLog.h"
#include "mapping/AliMpDEManager.h"
#include <cstring>
#include <iostream>
#include <iomanip>

ClassImp(AliHLTMUONMansoTrack);


std::ostream& operator << (
		std::ostream& stream,
		const AliHLTMUONMansoTrack& track
	)
{
/// Stream operator for std::ostream classes.
/// \param stream  The output stream object being written to.
/// \param track  The track object to print to the stream.
/// \returns  Returns 'stream'.

	stream	<< "ID: " << track.fId
		<< "; sign: " << track.fSign
		<< "; p = (" << track.Px()
		<< ", " << track.Py()
		<< ", " << track.Pz()
		<< "); chi^2: " << track.fChi2;
	return stream;
}


AliHLTMUONMansoTrack::AliHLTMUONMansoTrack(
		Int_t id, Int_t sign,
		Float_t px, Float_t py, Float_t pz,
		Float_t chi2,
		const AliHLTMUONTriggerRecord* trigrec,
		const AliHLTMUONRecHit* hit7,
		const AliHLTMUONRecHit* hit8,
		const AliHLTMUONRecHit* hit9,
		const AliHLTMUONRecHit* hit10,
		Float_t zf, Float_t qbl
	) :
	TObject(),
	fId(id),
	fSign(sign),
	fMomentum(px, py, pz),
	fChi2(chi2),
	fTrigRec(trigrec),
	fZmiddle(zf),
	fQBL(qbl)
{
/// Default constructor.
/// @param id       The track ID number which must be unique for any event.
/// @param sign     The particle's sign: -1, 1 or 0 if unknown.
/// @param px       X component of the particle's momentum (GeV/c).
/// @param py       Y component of the particle's momentum (GeV/c).
/// @param pz       Z component of the particle's momentum (GeV/c).
/// @param chi2     The chi squared of the track fit.
/// @param trigrec  Corresponding trigger record used as a seed to find
///                 this track.
/// @param hit7     Hit on chamber 7, tracking station 4.
/// @param hit8     Hit on chamber 8, tracking station 4.
/// @param hit9     Hit on chamber 9, tracking station 5.
/// @param hit10    Hit on chamber 10, tracking station 5.
/// @param zf    The Z coordinate of the middle of the magnetic field assumed
///              during momentum calculation.
/// @param qbl   The integrated magnetic field strength assumed during momentum
///              calculation.

	if (sign < -1 or 1 < sign)
	{
		AliError(Form("Trying to set the sign to %d. This is outside the"
			" valid range of [-1..1]", sign
		));
		fSign = 0;
	}
	
	fHit[0] = hit7;
	fHit[1] = hit8;
	fHit[2] = hit9;
	fHit[3] = hit10;
}


const AliHLTMUONRecHit* AliHLTMUONMansoTrack::Hit(Int_t chamber) const
{
/// Returns the hit on the specified chamber.
/// \param chamber  The chamber to return the hit for. Must be a value in the range [7..10].
/// \returns  A pointer to the hit object else NULL if it does not exist.

	if (7 <= chamber and chamber <= 10) return fHit[chamber - 7];
	
	AliError(Form(
		"Chamber number %d is not in the valid range [7..10].",
		int(chamber)
	));
	return NULL;
}


void AliHLTMUONMansoTrack::Print(Option_t* option) const
{
/// Prints the track information to standard output (screen).
/// \param option  Can be one of the following:
///           - "compact" - prints in a compact format.
///           - "detail" - prints track information in a more detailed format.
///           - "all" - prints a full dump of the track object.

	using namespace std;
	
	if (	option == NULL or strcmp(option, "") == 0 or
		strcmp(option, "compact") == 0
	   )
	{
		cout << *this << endl;
	}
	else if (strcmp(option, "detail") == 0)
	{
		cout << "Track ID = " << fId << "; sign = ";
		if (fSign != 0)
			cout << fSign;
		else
			cout << "unknown";
		cout << "; momentum: (px = " << Px()
			<< " GeV/c, py = " << Py()
			<< " GeV/c, pz = " << Pz()
			<< " GeV/c); chi^2 = " << fChi2 << endl;
		cout << "Used Zmiddle = " << fZmiddle << " cm and QBL = "
			<< fQBL << " T.m for the momentum calculation." << endl;
		
		streamsize w = cout.width();
		ios::fmtflags f = cout.flags();
		cout << setw(9) << "Chamber" << setw(0) << "  Hit:"
			<< setw(8) << "X (cm)"
			<< setw(12) << "Y (cm)"
			<< setw(12) << "Z (cm)" << endl;
		for (int i = 0; i < 4; i++)
		{
			cout << setw(9) << i+11;
			if (fHit[i] != NULL)
			{
				cout << setw(14) << fHit[i]->X()
					<< setw(12) << fHit[i]->Y()
					<< setw(12) << fHit[i]->Z();
			}
			else
			{
				cout << setw(14) << "-"
					<< setw(12) << "-"
					<< setw(12) << "-";
			}
			cout << endl;
		}
		cout.width(w); // reset the field width to previous value.
		cout.flags(f); // reset the flags to previous values.
	}
	else if (strcmp(option, "all") == 0)
	{
		cout << "Track ID = " << fId << "; sign = ";
		if (fSign != 0)
			cout << fSign;
		else
			cout << "unknown";
		cout << "; momentum: (px = " << Px()
			<< " GeV/c, py = " << Py()
			<< " GeV/c, pz = " << Pz()
			<< " GeV/c); chi^2 = " << fChi2 << endl;
		cout << "Used Zmiddle = " << fZmiddle << " cm and QBL = "
			<< fQBL << " T.m for the momentum calculation." << endl;
		
		for (int i = 0; i < 4; i++)
		{
			cout << "===== Hit on chamber: " << i+7 << " =====" << endl;
			if (fHit[i] != NULL)
				fHit[i]->Print("all");
			else
				cout << "No hit found." << endl;
		}
		
		cout << "===== Trigger Record =====" << endl;
		if (fTrigRec != NULL)
			fTrigRec->Print("all");
		else
			cout << "No trigger record associated with track." << endl;
	}
	else
	{
		AliError("Unknown option specified. Can only be one of 'compact',"
			" 'detail' or 'all'."
		);
	}
}


Int_t AliHLTMUONMansoTrack::Compare(const TObject* obj) const
{
/// We compare this object with 'obj' first by track ID, then by sign, then
/// by momentum and finally by chi2.
/// \param obj  This is the object to compare to. It must be of type AliHLTMUONMansoTrack.
/// \returns  -1 if 'this' is smaller than 'obj', 1 if greater and zero if both
///      objects are the same.

	if (obj->IsA() == AliHLTMUONMansoTrack::Class())
	{
		const AliHLTMUONMansoTrack* t =
			static_cast<const AliHLTMUONMansoTrack*>(obj);
		if (fId < t->fId) return -1;
		if (fId > t->fId) return 1;
		if (fSign < t->fSign) return -1;
		if (fSign > t->fSign) return 1;
		if (Px() < t->Px()) return -1;
		if (Px() > t->Px()) return 1;
		if (Py() < t->Py()) return -1;
		if (Py() > t->Py()) return 1;
		if (Pz() < t->Pz()) return -1;
		if (Pz() > t->Pz()) return 1;
		if (fChi2 < t->fChi2) return -1;
		if (fChi2 > t->fChi2) return 1;
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


bool AliHLTMUONMansoTrack::operator == (const AliHLTMUONMansoTrack& track) const
{
/// Comparison operator.
/// \param track  The track object to compare to.
/// \returns  true if 'this' object is identical to 'track' else false.

	return	fId == track.fId and fSign == track.fSign
		and fMomentum == track.fMomentum
		and fChi2 == track.fChi2 and fTrigRec == track.fTrigRec
		and fHit[0] == track.fHit[0] and fHit[1] == track.fHit[1]
		and fHit[2] == track.fHit[2] and fHit[3] == track.fHit[3];
}
