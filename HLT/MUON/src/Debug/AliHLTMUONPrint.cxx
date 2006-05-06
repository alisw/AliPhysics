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
// Print routines to display internal dHLT data on the console.
// The << operators are overloaded to write to ostreams.
//
////////////////////////////////////////////////////////////////////////////////

#include "Debug/AliHLTMUONPrint.h"
#include "AliHLTMUONCoreEventID.h"
#include "AliHLTMUONCorePoint.h"
#include "AliHLTMUONCoreTriggerRecord.h"
#include "AliHLTMUONCoreRegionOfInterest.h"

#include <iostream>
using std::endl;
using std::cout;
using std::ostream;


ostream& operator << (ostream& os, const AliHLTMUONCoreEventID& id)
{
	os << "<" << id.fBunch << ":" << id.fTimeStamp << ">";
	return os;
}


ostream& operator << (ostream& os, const AliHLTMUONCorePoint& p)
{
	os << "[" << p.X() << ", " << p.Y() << "]";
	return os;
}


ostream& operator << (ostream& os, const AliHLTMUONCoreParticleSign s)
{
	switch (s)
	{
	case kSignMinus:   os << "Minus";   break;
	case kSignPlus:    os << "Plus";    break;
	case kUnknownSign: os << "Unknown"; break;
	default:           os << "FAULT!";
	}
	return os;
}


ostream& operator << (ostream& os, const AliHLTMUONCoreTriggerRecord& rec)
{
	char* signstr;
	switch (rec.fSign)
	{
	case kSignMinus:   signstr = "Minus  "; break;
	case kSignPlus:    signstr = "Plus   "; break;
	case kUnknownSign: signstr = "Unknown"; break;
	default:           signstr = "FAULT!!";
	}
	
	os << "{ sign = " << signstr << ", pt = " 
	   << rec.fPt << ", impact on station 1: "
	   << rec.fStation1impact << " , station 2: "
	   << rec.fStation2impact << " }";
	return os;
}


ostream& operator << (ostream& os, const AliHLTMUONCoreChamberID chamber)
{
	Int ch = (Int)chamber;
	if (0 <= ch && ch <= 9)
		os << ch + 1;
	else
		os << "FAULT!!";
	return os;
}

