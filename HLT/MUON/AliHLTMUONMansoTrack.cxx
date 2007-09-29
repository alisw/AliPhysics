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

/**
 * @file   AliHLTMUONMansoTrack.cxx
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   
 * @brief  Implementation of the AliHLTMUONMansoTrack class.
 */

#include "AliHLTMUONMansoTrack.h"
#include "AliLog.h"
#include "mapping/AliMpDEManager.h"
#include <cstring>
#include <iostream>
#include <iomanip>
using namespace std;

ClassImp(AliHLTMUONMansoTrack);


std::ostream& operator << (
		std::ostream& stream,
		const AliHLTMUONMansoTrack& track
	)
{
// Stream operator.

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
	fId(id), fSign(sign), fMomentum(px, py, pz),
	fChi2(chi2), fTrigRec(trigrec), fZmiddle(zf), fQBL(qbl)
{
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
// Returns the hit on the specified chamber.

	if (7 <= chamber and chamber <= 10) return fHit[chamber - 7];
	
	AliError(Form(
		"Chamber number %d is not in the valid range [7..10].",
		int(chamber)
	));
	return NULL;
}


void AliHLTMUONMansoTrack::Print(Option_t* option) const
{
// Prints the track information to standard output (screen).

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
// We compare this object with 'obj' first by track ID, then by sign, then
// by momentum and finally by chi2.

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
	return	fId == track.fId and fSign == track.fSign
		and fMomentum == track.fMomentum
		and fChi2 == track.fChi2 and fTrigRec == track.fTrigRec
		and fHit[0] == track.fHit[0] and fHit[1] == track.fHit[1]
		and fHit[2] == track.fHit[2] and fHit[3] == track.fHit[3];
}
