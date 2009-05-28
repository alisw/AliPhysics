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
/// @file   AliHLTMUONRecHit.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   29 Sep 2007
/// @brief  Implementation of the AliHLTMUONRecHit class.
///
/// The AliHLTMUONRecHit object is used to store 3D hit coordinates translated
/// from dHLT raw data.
///

#include "AliHLTMUONRecHit.h"
#include "AliLog.h"
#include "mapping/AliMpDEManager.h"
#include <cstring>
#include <iostream>
#include <iomanip>

ClassImp(AliHLTMUONRecHit);
ClassImp(AliHLTMUONRecHit::AliChannel);


std::ostream& operator << (std::ostream& stream, const AliHLTMUONRecHit& hit)
{
/// Stream operator for std::ostream classes.
/// \param stream  The output stream object being written to.
/// \param track  The hit object to print to the stream.
/// \returns  Returns 'stream'.

	stream << "(" << hit.X() << ", " << hit.Y() << ", " << hit.Z() << ")";
	return stream;
}


void AliHLTMUONRecHit::SetDebugInfo(
		Int_t detElemId, Int_t clusterId, UInt_t nChExp,
		Float_t charge, Int_t sourceDDL
	)
{
/// Sets the debugging information.
/// @param detElemId  The detector element ID.
/// @param clusterId  Cluster ID of the hit's cluster.
/// @param nChExp     Number of expected channels forming the cluster.
/// @param charge  The total charge of the cluster.
/// @param sourceDDL  The source DDL of this hit.

	fSourceDDL = sourceDDL;
	fDetElemId = detElemId;
	fClusterId = clusterId;
	fNchExp = nChExp;
	fCharge = charge;
}
	

Int_t AliHLTMUONRecHit::Chamber(bool warn) const
{
/// Returns the chamber ID for this hit.
/// \param warn  Indicates if any warning should be printed in case of problems.
/// \returns The chamber number of this hit in the range [1..14] or -1 if not known.

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
/// Adds a new channel to the channels list forming this hit's cluster.
/// @param manu    The MANU number
/// @param channel The MANU channel address.
/// @param signal  The ADC signal value measured on the channel.
/// @param rawDataWord This is the raw data word as read from the DDL.

	Int_t index = fChannels.GetEntriesFast();
	new (fChannels[index]) AliChannel(manu, channel, signal, rawDataWord);
}


void AliHLTMUONRecHit::Print(Option_t* option) const
{
/// Prints the coordinates of this hit to standard output (screen).
/// \param option  Can be one of the following:
///           - "compact" - prints in a compact format.
///           - "detail" - prints hit information in a more detailed format.
///           - "all" - prints a full dump of the hit object.

	using namespace std;

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
			<< "; total charge = " << fCharge
			<< "; expected #ch = " << fNchExp << endl;
	}
	else if (strcmp(option, "all") == 0)
	{
		cout << "(x = " << X() << " cm, y = " << Y()
			<< " cm, z = " << Z()
			<< " cm); source DDL = " << fSourceDDL
			<< "; DetElemID = " << fDetElemId
			<< "; cluster ID = " << fClusterId
			<< "; total charge = " << fCharge
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
				const AliHLTMUONRecHit::AliChannel* c =
					static_cast<const AliHLTMUONRecHit::AliChannel*>(fChannels[i]);
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
/// We compare this object with 'obj' first by X, then Y, then Z.
/// \param obj  This is the object to compare to. It must be of type AliHLTMUONRecHit.
/// \returns  -1 if 'this' is smaller than 'obj', 1 if greater and zero if both
///      objects are the same.

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


std::ostream& operator << (std::ostream& stream, const AliHLTMUONRecHit::AliChannel& c)
{
/// Stream operator for std::ostream classes.
/// \param stream  The output stream object being written to.
/// \param c  The channel object to print to the stream.
/// \returns  Returns 'stream'.

	stream << "Channel: " << c.fManu << " , " << c.fAddress
		<< "; ADC: " << c.fSignal;
	return stream;
}


void AliHLTMUONRecHit::AliChannel::Print(Option_t* option) const
{
/// Prints the details of this channel to standard output (screen).
/// \param option  Can be one of the following:
///           - "compact" - prints in a compact format.
///           - "detail" - prints channel information in a more detailed format.

	using namespace std;
	
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


Int_t AliHLTMUONRecHit::AliChannel::Compare(const TObject* obj) const
{
/// We compare this object with 'obj' first by MANU number, then by MANU channel
/// address, then ADC signal.
/// \param obj  This is the object to compare to. It must be of type AliHLTMUONRecHit::Channel.
/// \returns  -1 if 'this' is smaller than 'obj', 1 if greater and zero if both
///      objects are the same.

	if (obj->IsA() == AliChannel::Class())
	{
		const AliChannel* c = static_cast<const AliChannel*>(obj);
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
