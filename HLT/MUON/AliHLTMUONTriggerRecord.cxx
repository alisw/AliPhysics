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
/// @file   AliHLTMUONTriggerRecord.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   29 Sep 2007
/// @brief  Implementation of the AliHLTMUONTriggerRecord class.
///
/// Code for the trigger record structure containing data corresponding to the
/// L0 trigger local board output.
///

#include "AliHLTMUONTriggerRecord.h"
#include "AliLog.h"
#include "mapping/AliMpDEManager.h"
#include <cstring>
#include <iostream>
#include <iomanip>
using namespace std;

ClassImp(AliHLTMUONTriggerRecord);


std::ostream& operator << (
		std::ostream& stream,
		const AliHLTMUONTriggerRecord& trigrec
	)
{
/// Stream operator for std::ostream classes.
/// \param stream  The output stream object being written to.
/// \param trigrec  The trigger record object to print to the stream.
/// \returns  Returns 'stream'.

	stream	<< "ID: " << trigrec.fId
		<< "; sign: " << trigrec.fSign
		<< "; p = (" << trigrec.Px()
		<< ", " << trigrec.Py()
		<< ", " << trigrec.Pz()
		<< ")";
	return stream;
}


AliHLTMUONTriggerRecord::AliHLTMUONTriggerRecord(
		Int_t id, Int_t sign,
		Float_t px, Float_t py, Float_t pz,
		Int_t sourceDDL, Float_t zf, Float_t qbl
	) :
	TObject(),
	fId(id),
	fSign(sign),
	fMomentum(px, py, pz),
	fSourceDDL(sourceDDL),
	fZmiddle(zf), fQBL(qbl)
{
/// Constructor for creating a new trigger record.
/// @param id    The trigger record ID number unique for an event.
/// @param sign  The particle's sign. Must be -1, 1 or 0 if the sign is unknown.
/// @param px    X component of the particle's momentum.
/// @param py    Y component of the particle's momentum.
/// @param pz    Z component of the particle's momentum.
/// @param sourceDDL  The DDL from which this trigger record originates.
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

	// Fill the debugging information to invalid values by default.
	for (int i = 0; i < 4; i++)
	{
		fDetElemId[i] = -1;
		fPatternX[i] = -1;
		fPatternY[i] = -1;
	}
}


const TVector3& AliHLTMUONTriggerRecord::Hit(Int_t chamber) const
{
/// Returns the hit on the specified chamber.
/// \param chamber  The chamber for which to fetch the hit. Valid values
///                 are in the range [11..14].
/// \returns  The 3D corrdinate of the hit.

	if (11 <= chamber and chamber <= 14) return fHit[chamber - 11];
	
	AliError(Form(
		"Chamber number %d is not in the valid range [11..14].",
		int(chamber)
	));
	return fHit[0];
}


Int_t AliHLTMUONTriggerRecord::DetElemId(Int_t chamber) const
{
/// Returns the detector element ID for the specified chamber associated
/// to the hit on that chamber.
/// @param chamber  The chamber for which to fetch the detector element ID.
///                Valid values are in the range [11..14].
/// \returns  The detector element ID or -1 if not known.

	if (11 <= chamber and chamber <= 14) return fDetElemId[chamber - 11];
	
	AliError(Form(
		"Chamber number %d is not in the valid range [11..14].",
		int(chamber)
	));
	return fDetElemId[0];
}


Int_t AliHLTMUONTriggerRecord::PatternX(Int_t chamber) const
{
/// Returns the raw data X pattern of the hit on the specified chamber.
/// \param chamber  The chamber for which to fetch the bit pattern.
///                 Valid values are in the range [11..14].
/// \returns X bit pattern of the hit.

	if (11 <= chamber and chamber <= 14) return fPatternX[chamber - 11];
	
	AliError(Form(
		"Chamber number %d is not in the valid range [11..14].",
		int(chamber)
	));
	return fPatternX[0];
}


Int_t AliHLTMUONTriggerRecord::PatternY(Int_t chamber) const
{
/// Returns the raw data Y pattern of the hit on the specified chamber.
/// \param chamber  The chamber for which to fetch the bit pattern.
///                 Valid values are in the range [11..14].
/// \returns Y bit pattern of the hit.

	if (11 <= chamber and chamber <= 14) return fPatternY[chamber - 11];
	
	AliError(Form(
		"Chamber number %d is not in the valid range [11..14].",
		int(chamber)
	));
	return fPatternY[0];
}


void AliHLTMUONTriggerRecord::SetHit(Int_t chamber, Float_t x, Float_t y, Float_t z)
{
/// Fills the hit coordinate on the specified chamber.
/// @param chamber  The chamber for which to set the hit. Valid values
///                 are in the range [11..14].
/// @param x  The X coordinate of the hit in centimetres.
/// @param y  The Y coordinate of the hit in centimetres.
/// @param z  The Z coordinate of the hit in centimetres.

	if (11 <= chamber and chamber <= 14)
	{
		fHit[chamber - 11].SetXYZ(x, y, z);
	}
	else
	{
		AliError(Form(
			"Chamber number %d is not in the valid range [11..14].",
			int(chamber)
		));
	}
}


void AliHLTMUONTriggerRecord::SetHitDebugInfo(
		Int_t chamber,
		Int_t detElemId, UShort_t patternX, UShort_t patternY
	)
{
/// Fills the debugging information corresponding to the hit on the specified chamber.
/// Sets the debugging information for the hit on the specified chamber.
/// @param chamber  The chamber for which to set the debugging information.
///                Valid values are in the range [11..14].
/// @param detElemId  The detector element ID.
/// @param patterX    The X bit pattern from the local board.
/// @param patterY    The Y bit pattern from the local board.

	if (11 <= chamber and chamber <= 14)
	{
		fDetElemId[chamber - 11] = detElemId;
		fPatternX[chamber - 11] = patternX;
		fPatternY[chamber - 11] = patternY;
	}
	else
	{
		AliError(Form(
			"Chamber number %d is not in the valid range [11..14].",
			int(chamber)
		));
	}
}


void AliHLTMUONTriggerRecord::Print(Option_t* option) const
{
/// Prints the trigger record to standard output (screen).
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
		cout << "Trigger record ID = " << fId << "; sign = ";
		if (fSign != 0)
			cout << fSign;
		else
			cout << "unknown";
		cout << "; momentum: (px = " << Px()
			<< " GeV/c, py = " << Py()
			<< " GeV/c, pz = " << Pz()
			<< " GeV/c)" << endl;
		cout << "Used Zmiddle = " << fZmiddle << " cm and QBL = "
			<< fQBL << " T.m for the momentum calculation." << endl;
	}
	else if (strcmp(option, "all") == 0)
	{
		cout << "Trigger record ID = " << fId << "; sign = ";
		if (fSign != 0)
			cout << fSign;
		else
			cout << "unknown";
		cout << "; momentum: (px = " << Px()
			<< " GeV/c, py = " << Py()
			<< " GeV/c, pz = " << Pz()
			<< " GeV/c)" << endl;
		cout << "Used Zmiddle = " << fZmiddle << " cm and QBL = "
			<< fQBL << " T.m for the momentum calculation." << endl;
		
		streamsize w = cout.width();
		ios::fmtflags f = cout.flags();
		cout << setw(9) << "Chamber" << setw(0) << "  Hit:"
			<< setw(8) << "X (cm)"
			<< setw(12) << "Y (cm)"
			<< setw(12) << "Z (cm)"
			<< setw(12) << "DetElemID"
			<< setw(18) << "X bit pattern"
			<< setw(18) << "Y bit pattern" << endl;
		for (int i = 0; i < 4; i++)
		{
			cout << setw(9) << i+11
				<< setw(14) << fHit[i].X()
				<< setw(12) << fHit[i].Y()
				<< setw(12) << fHit[i].Z()
				<< setw(12) << fDetElemId[i];
			if (fPatternX[i] != -1)
			{
				// Print the X pattern as a bit pattern.
				cout << "  ";
				for (Int_t j = 15; j >= 0; j--)
					cout << (((fPatternX[i] & (1 << j)) > 0) ? "1" : "0");
			}
			else
			{
				cout << "  ----------------";
			}
			if (fPatternY[i] != -1)
			{
				// Print the Y pattern as a bit pattern.
				cout << "  ";
				for (Int_t j = 15; j >= 0; j--)
					cout << (((fPatternY[i] & (1 << j)) > 0) ? "1" : "0");
			}
			else
			{
				cout << "  ----------------";
			}
			cout << endl;
		}
		cout.width(w); // reset the field width to previous value.
		cout.flags(f); // reset the flags to previous values.
	}
	else
	{
		AliError("Unknown option specified. Can only be one of 'compact',"
			" 'detail' or 'all'."
		);
	}
}


Int_t AliHLTMUONTriggerRecord::Compare(const TObject* obj) const
{
/// We compare this object with 'obj' first by trigger record ID, then
/// by sign and finally by momentum.
/// \param obj  This is the object to compare to. It must be of type AliHLTMUONTriggerRecord.
/// \returns  -1 if 'this' is smaller than 'obj', 1 if greater and zero if both
///      objects are the same.

	if (obj->IsA() == AliHLTMUONTriggerRecord::Class())
	{
		const AliHLTMUONTriggerRecord* tr =
			static_cast<const AliHLTMUONTriggerRecord*>(obj);
		if (fId < tr->fId) return -1;
		if (fId > tr->fId) return 1;
		if (fSign < tr->fSign) return -1;
		if (fSign > tr->fSign) return 1;
		if (Px() < tr->Px()) return -1;
		if (Px() > tr->Px()) return 1;
		if (Py() < tr->Py()) return -1;
		if (Py() > tr->Py()) return 1;
		if (Pz() < tr->Pz()) return -1;
		if (Pz() > tr->Pz()) return 1;
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


bool AliHLTMUONTriggerRecord::operator == (const AliHLTMUONTriggerRecord& trigrec) const
{
/// Compares the trigger record 'trigrec' to this one.
/// \param trigrec  The trigger record object to compare to.
/// \returns  true if 'this' object is identical to 'trigrec' else false.

	return	fId == trigrec.fId and fSign == trigrec.fSign
		and fMomentum == trigrec.fMomentum
		and fHit[0] == trigrec.fHit[0] and fHit[1] == trigrec.fHit[1]
		and fHit[2] == trigrec.fHit[2] and fHit[3] == trigrec.fHit[3];
}
