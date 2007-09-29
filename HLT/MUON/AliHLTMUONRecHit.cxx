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
 * @file   AliHLTMUONRecHit.cxx
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   
 * @brief  Implementation of the AliHLTMUONRecHit class.
 */

#include "AliHLTMUONRecHit.h"
#include "AliLog.h"
#include "mapping/AliMpDEManager.h"
#include <cstring>
#include <iostream>
#include <iomanip>
using namespace std;

ClassImp(AliHLTMUONRecHit);
ClassImp(AliHLTMUONRecHit::Channel);


std::ostream& operator << (std::ostream& stream, const AliHLTMUONRecHit& hit)
{
// Stream operator for std::ostream classes.

	stream << "(" << hit.X() << ", " << hit.Y() << ", " << hit.Z() << ")";
	return stream;
}


void AliHLTMUONRecHit::SetDebugInfo(
		Int_t detElemId, Int_t clusterId, UInt_t nChExp, Int_t sourceDDL
	)
{
// Sets the debugging information.

	fSourceDDL = sourceDDL;
	fDetElemId = detElemId;
	fClusterId = clusterId;
	fNchExp = nChExp;
}
	

Int_t AliHLTMUONRecHit::Chamber(bool warn) const
{
// Returns the chamber ID for this hit.

	if (fSourceDDL != -1)
	{
		if (fSourceDDL < 1 or fSourceDDL > 20)
		{
			return ((fSourceDDL-1) / 2) + 1;
		}
		else if (warn)
		{
			AliError(Form("The DDL source number: %d is out of range."
				" Valid values are [1..20]", fSourceDDL
			));
		}
	}
	if (fDetElemId != -1) return AliMpDEManager::GetChamberId(fDetElemId, warn);
	
	if (warn)
	{
		AliWarning("Neither the DDL source nor the detector element ID was not set,"
			" so we do not know on which chamber this hit was reconstructed."
		);
	}
	return -1;
}


void AliHLTMUONRecHit::AddChannel(
		Short_t manu, Short_t channel, Short_t signal,
		UInt_t rawDataWord
	)
{
// Adds a new channel to the channels list forming this hit's cluster.

	Int_t index = fChannels.GetEntriesFast();
	new (fChannels[index]) Channel(manu, channel, signal, rawDataWord);
}


void AliHLTMUONRecHit::Print(Option_t* option) const
{
// Prints the coordinates of this hit to standard output (screen).

	if (	option == NULL or strcmp(option, "") == 0 or
		strcmp(option, "compact") == 0
	   )
	{
		cout << *this << endl;
	}
	else if (strcmp(option, "detail") == 0)
	{
		cout << "(x = " << X() << " cm, y = " << Y()
			<< " cm, z = " << Z()
			<< " cm); source DDL = " << fSourceDDL
			<< "; DetElemID = " << fDetElemId
			<< "; cluster ID = " << fClusterId
			<< "; expected #ch = " << fNchExp << endl;
	}
	else if (strcmp(option, "all") == 0)
	{
		cout << "(x = " << X() << " cm, y = " << Y()
			<< " cm, z = " << Z()
			<< " cm); source DDL = " << fSourceDDL
			<< "; DetElemID = " << fDetElemId
			<< "; cluster ID = " << fClusterId
			<< "; expected #ch = " << fNchExp << endl;
		if (fChannels.GetEntriesFast() == 0)
		{
			cout << "No channels found for this hit." << endl;
		}
		else
		{
			streamsize w = cout.width();
			ios::fmtflags f = cout.flags();
			cout << setw(12) << "MANU"
				<< setw(12) << "Channel"
				<< setw(12) << "Signal"
				<< setw(15) << "Raw data word" << endl;
			cout << showbase;
			for (Int_t i = 0; i < fChannels.GetEntriesFast(); i++)
			{
				const AliHLTMUONRecHit::Channel* c =
					static_cast<const AliHLTMUONRecHit::Channel*>(fChannels[i]);
				cout << dec << setw(12) << c->Manu()
					<< setw(12) << c->Address()
					<< setw(12) << c->Signal()
					<< "     " << hex << setw(10) << internal;
				ios::char_type fc = cout.fill('0');
				cout << c->RawDataWord() << right << endl;
				cout.fill(fc);
			}
			cout.width(w); // reset the field width to previous value.
			cout.flags(f); // reset the flags to previous values.
		}
	}
	else
	{
		AliError("Unknown option specified. Can only be one of 'compact',"
			" 'detail' or 'all'."
		);
	}
}


Int_t AliHLTMUONRecHit::Compare(const TObject* obj) const
{
// We compare this object with 'obj' first by X, then Y, then Z.

	if (obj->IsA() == AliHLTMUONRecHit::Class())
	{
		const AliHLTMUONRecHit* h = static_cast<const AliHLTMUONRecHit*>(obj);
		if (X() < h->X()) return -1;
		if (X() > h->X()) return 1;
		if (Y() < h->Y()) return -1;
		if (Y() > h->Y()) return 1;
		if (Z() < h->Z()) return -1;
		if (Z() > h->Z()) return 1;
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


std::ostream& operator << (std::ostream& stream, const AliHLTMUONRecHit::Channel& c)
{
// Stream operator for std::ostream classes.

	stream << "Channel: " << c.fManu << " , " << c.fAddress
		<< "; ADC: " << c.fSignal;
	return stream;
}


void AliHLTMUONRecHit::Channel::Print(Option_t* option) const
{
// Prints the details of this channel to standard output (screen).

	if (	option == NULL or strcmp(option, "") == 0 or
		strcmp(option, "compact") == 0
	   )
	{
		cout << *this << endl;
	}
	else if (strcmp(option, "detail") == 0)
	{
		streamsize w = cout.width();
		ios::fmtflags f = cout.flags();
		cout << "MANU = " << fManu << ", Channel address = " << fAddress
			<< ", Signal = " << fSignal
			<< "; Raw data word = " << hex << setw(10) << internal;
		ios::char_type fc = cout.fill('0');
		cout << fRawDataWord << endl;
		cout.fill(fc); // reset fill character
		cout.width(w); // reset the field width to previous value.
		cout.flags(f); // reset the flags to previous values.
	}
	else
	{
		AliError("Unknown option specified. Can only be one of"
			" 'compact' or 'detail'."
		);
	}
}


Int_t AliHLTMUONRecHit::Channel::Compare(const TObject* obj) const
{
// We compare this object with 'obj' first by MANU number, then by MANU channel
// address, then ADC signal.

	if (obj->IsA() == Channel::Class())
	{
		const Channel* c = static_cast<const Channel*>(obj);
		if (fManu < c->Manu()) return -1;
		if (fManu > c->Manu()) return 1;
		if (fAddress < c->Address()) return -1;
		if (fAddress > c->Address()) return 1;
		if (fSignal < c->Signal()) return -1;
		if (fSignal > c->Signal()) return 1;
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
