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

/*
$Log$
Revision 1.1.1.1  2005/09/12 22:11:40  byordano
SHUTTLE package

Revision 1.2  2005/08/30 10:53:23  byordano
some more descriptions added

*/

//
// This class represents the main value structure
// which forms so called 'historical data' in any SCADA system.
// When a value (which represent a parameter of some real world object)
// is measured in the time, a value serie (called value set) is formed.
// Each element of this value series has two fields:
// fValue - primitive value which represents the real measured value
// fTimestamp - timestamp when the measurement was made
//

#include "AliDCSValue.h"

#include "TTimeStamp.h"

ClassImp(AliDCSValue)

AliDCSValue::AliDCSValue() {

}

AliDCSValue::AliDCSValue(const AliSimpleValue& value, UInt_t timeStamp):
	fValue(value), fTimeStamp(timeStamp)
{

}


TString AliDCSValue::ToString() const {

	return fValue.ToString() + ", Timestmap: " +
		 TTimeStamp(fTimeStamp).AsString();
}


