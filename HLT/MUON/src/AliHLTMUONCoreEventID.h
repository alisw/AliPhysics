#ifndef ALIHLTMUONCOREEVENTID_H
#define ALIHLTMUONCOREEVENTID_H
/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
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

////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "AliHLTMUONBasicTypes.h"
#include "AliHLTMUONUtils.h"


// Must correct this definition!!!!
struct AliHLTMUONCoreEventID
{
	UInt fBunch;
	UInt fTimeStamp;


	AliHLTMUONCoreEventID(UInt bunch = 0, UInt timestamp = 0)
		: fBunch(bunch), fTimeStamp(timestamp)
	{};


	bool operator == (const AliHLTMUONCoreEventID& rhs)
	{
		return fBunch == rhs.fBunch and fTimeStamp == rhs.fTimeStamp;
	};

	bool operator != (const AliHLTMUONCoreEventID& rhs)
	{
		return not (*this == rhs);
	};
};


#endif // ALIHLTMUONCOREEVENTID_H
