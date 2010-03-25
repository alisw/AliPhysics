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

// $Id: AliHLTMUONTrack.cxx 37131 2009-11-23 11:08:09Z aszostak $

///
/// @file   AliHLTMUONTrack.cxx
/// @author Indranil Das <indra.ehep@gmail.com> and Artur Szostak <artursz@iafrica.com>
/// @date   10 March 2010
/// @brief  Implementation of the AliHLTMUONTrack class.
///
/// The track class is used to store converted track data from dHLT raw internal
/// data blocks as a ROOT object. This allows storing these in ROOT files.
///

#include "AliHLTMUONTrack.h"
#include "AliHLTMUONRecHit.h"
#include "AliHLTMUONTriggerRecord.h"
#include "AliLog.h"
#include "mapping/AliMpDEManager.h"
#include <cstring>
#include <iostream>
#include <iomanip>

ClassImp(AliHLTMUONTrack);


std::ostream& operator << (
		std::ostream& stream,
		const AliHLTMUONTrack& track
	)
{
/// Stream operator for std::ostream classes.
/// \param stream  The output stream object being written to.
/// \param track  The track object to print to the stream.
/// \returns  the stream object that was written to, <i>stream</i>.

	stream	<< "ID: " << track.fId
		<< "; sign: " << track.fSign
		<< "; p = (" << track.Px()
		<< ", " << track.Py()
		<< ", " << track.Pz()
		<< "); chi^2: " << track.fChi2;
	return stream;
}


AliHLTMUONTrack::AliHLTMUONTrack(
		Int_t id, Int_t sign,
		Float_t px, Float_t py, Float_t pz,
		Float_t invmom, Float_t thetaX, Float_t thetaY,
		Float_t x, Float_t y, Float_t z,
		Float_t chi2,
		const AliHLTMUONTriggerRecord* trigrec,
		const AliHLTMUONRecHit* hits[16]
	) :
	TObject(),
	fId(id),
	fSign(sign),
	fInverseBendingMomentum(invmom),
	fThetaX(thetaX),
	fThetaY(thetaY),
	fMomentum(px, py, pz),
	fVertexDCA(x, y, z),
	fChi2(chi2),
	fTrigRec(trigrec)
{
/// Default constructor.
/// @param id       The track ID number which must be unique for any event.
/// @param sign     The particle's sign: -1, 1 or 0 if unknown.
/// @param px       X component of the particle's momentum (GeV/c).
/// @param py       Y component of the particle's momentum (GeV/c).
/// @param pz       Z component of the particle's momentum (GeV/c).
/// @param invmom   Inverse bending momentum (GeV/c).
/// @param thetaX   The non-bending plane slope of the fitted track.
/// @param thetaY   The bending plane slope of the fitted track.
/// @param x        X coordinate of the particle's distance of closest
///                 approach (DCA) position (cm).
/// @param y        Y coordinate of the particle's DCA position (cm).
/// @param z        Z coordinate of the particle's DCA position (cm).
/// @param chi2     The chi squared of the track fit.
/// @param trigrec  Corresponding trigger record used as a seed to find
///                 this track.
/// @param hits     The array of 16 hit coordinates found for the track.
///                 If NULL then then all hit coordinates are set to empty.

	if (sign < -1 or 1 < sign)
	{
		AliError(Form("Trying to set the sign to %d. This is outside the"
			" valid range of [-1..1]", sign
		));
		fSign = 0;
	}
	
	// Copy individual hit pointers if hits were set.
	if (hits != NULL)
	{
		for (int i = 0; i < 16; ++i)
		{
			fHit[i] = hits[i];
		}
	}
	else
	{
		for (int i = 0; i < 16; ++i)
		{
			fHit[i] = NULL;
		}
	}
}


const AliHLTMUONRecHit* AliHLTMUONTrack::Hit(Int_t i) const
{
/// Returns the i'th hit.
/// \param i  The number of the hit to return. Must be a value in the range [0..15].
/// \returns  A pointer to the hit object else NULL if it does not exist.

	if (0 <= i and i <= 15) return fHit[i];
	AliError(Form("Hit index number %d is not in the valid range [0..15].", int(i)));
	return NULL;
}


const AliHLTMUONRecHit* AliHLTMUONTrack::HitByChamber(Int_t chamber) const
{
/// Returns the first hit found for the specified chamber.
/// \param chamber  The chamber to return the hit for. Must be a value in the range [1..14].
/// \returns  A pointer to the hit object else NULL if it does not exist.

	if (not (1 <= chamber and chamber <= 14))
	{
		AliError(Form(
			"Chamber number %d is not in the valid range [1..14].",
			int(chamber)
		));
		return NULL;
	}
	
	const AliHLTMUONRecHit* hit = NULL;
	for (int i = 0; i < 16; ++i)
	{
		if (fHit[i] == NULL) continue;
		if (fHit[i]->Chamber(false) == chamber)
		{
			hit = fHit[i];
			break;
		}
	}
	return hit;
}


void AliHLTMUONTrack::Print(Option_t* option) const
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
			<< " GeV/c); vertex DCA: (x = " << X()
			<< " cm, y = " << Y()
			<< " cm, z = " << Z()
			<< " cm); chi^2 = " << fChi2 << endl;
		
		streamsize w = cout.width();
		ios::fmtflags f = cout.flags();
		cout << setw(9) << "Hit no."
			<< setw(9) << "Chamber"
			<< setw(0) << "  Pos:"
			<< setw(8) << "X (cm)"
			<< setw(12) << "Y (cm)"
			<< setw(12) << "Z (cm)" << endl;
		for (int i = 0; i < 16; i++)
		{
			cout << setw(9) << i;
			if (fHit[i] != NULL)
			{
				cout << setw(9) << fHit[i]->Chamber(false)
					<< setw(14) << fHit[i]->X()
					<< setw(12) << fHit[i]->Y()
					<< setw(12) << fHit[i]->Z();
			}
			else
			{
				cout << setw(9) << "-"
					<< setw(14) << "-"
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
			<< " GeV/c); vertex DCA: (x = " << X()
			<< " cm, y = " << Y()
			<< " cm, z = " << Z()
			<< " cm); chi^2 = " << fChi2 << endl;
		cout << "Inverse bending momentum = " << fInverseBendingMomentum
			<< "; non-bending plane slope = " << fThetaX
			<< "; bending plane slope = " << fThetaY << endl;
		
		streamsize w = cout.width();
		ios::fmtflags f = cout.flags();
		cout.width(w); // reset the field width to previous value.
		cout.flags(f); // reset the flags to previous values.
		
		for (int i = 0; i < 16; i++)
		{
			cout << "======== Hit " << i << " ========" << endl;
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


Int_t AliHLTMUONTrack::Compare(const TObject* obj) const
{
/// We compare this object with 'obj' first by track ID, then by sign, then
/// by momentum, then by DCA vertex and finally by chi2.
/// \param obj  This is the object to compare to. It must be of type AliHLTMUONTrack.
/// \returns  -1 if 'this' is smaller than 'obj', 1 if greater and zero if both
///      objects are the same.

	if (obj->IsA() == AliHLTMUONTrack::Class())
	{
		const AliHLTMUONTrack* t =
			static_cast<const AliHLTMUONTrack*>(obj);
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
		if (X() < t->X()) return -1;
		if (X() > t->X()) return 1;
		if (Y() < t->Y()) return -1;
		if (Y() > t->Y()) return 1;
		if (Z() < t->Z()) return -1;
		if (Z() > t->Z()) return 1;
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


bool AliHLTMUONTrack::operator == (const AliHLTMUONTrack& track) const
{
/// Comparison operator.
/// \param track  The track object to compare to.
/// \returns  true if 'this' object is identical to 'track' else false.

	for (int i = 0; i < 16; i++)
	{
		if (fHit[i] != track.fHit[i]) return false;
	}
	
	return	fId == track.fId and fSign == track.fSign
		and fInverseBendingMomentum == track.fInverseBendingMomentum
		and fThetaX == track.fThetaX and fThetaY == track.fThetaY
		and fMomentum == track.fMomentum and fVertexDCA == track.fVertexDCA
		and fChi2 == track.fChi2 and fTrigRec == track.fTrigRec;
}
