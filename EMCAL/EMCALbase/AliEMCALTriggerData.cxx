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
 * about the suitability of this software for any purpeateose. It is      *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*
 
 
EMCal trigger data container
for data (both raw & rec) persistency
Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/

#include "AliEMCALTriggerData.h"
//#include "AliEMCALTriggerPatch.h"
#include "AliLog.h"
#include "TIterator.h"
#include "Riostream.h"

ClassImp(AliEMCALTriggerData)

//_____________
AliEMCALTriggerData::AliEMCALTriggerData() : TObject(),
fMode(0),
fL1GammaThreshold(),
fL1JetThreshold(),
fL1V0(),
fL1FrameMask(0),
fL1TriggerType(),
fL1DataDecoded(0),
fL1RawData(0),
fMedian(0)
{  
	// Ctor
		
	fL1GammaThreshold[0] = fL1GammaThreshold[1] = 0;
	fL1JetThreshold[0] = fL1JetThreshold[1] = 0;
	
	fL1V0[0] = fL1V0[1] = 0;
	for (Int_t i = 0; i < 19; i++) fL1TriggerType[i] = 0;	
}

//_____________
AliEMCALTriggerData::~AliEMCALTriggerData()
{
	// Dtor	
}

//_____________
void AliEMCALTriggerData::Scan() const
{
	// Dump

	for (int i = 0; i < 2; i++){
		printf("\tL1 thresholds[%d]: gamma %d\tjet %d\n", i, fL1GammaThreshold[i], fL1JetThreshold[i]);
	}
}

//_____________
void AliEMCALTriggerData::Reset()
{
	// Reset
	
	fL1DataDecoded = 0;
}



