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

 


 
Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/

#include "AliEMCALTriggerDCSConfig.h"
#include "AliEMCALTriggerSTUDCSConfig.h"
#include "AliEMCALTriggerTRUDCSConfig.h"

ClassImp(AliEMCALTriggerDCSConfig)

//_____________________________________________________________________________
AliEMCALTriggerDCSConfig::AliEMCALTriggerDCSConfig() : TObject()
,fTRUArr(new TClonesArray("AliEMCALTriggerTRUDCSConfig",32))
,fSTUObj(new AliEMCALTriggerSTUDCSConfig())
{
	//
	// AliEMCALTriggerDCSConfig default constructor
	//
}

//_____________________________________________________________________________
AliEMCALTriggerDCSConfig::~AliEMCALTriggerDCSConfig()
{
	//
	//
	//
	delete fTRUArr; fTRUArr = 0x0;
	delete fSTUObj; fSTUObj = 0x0;
}
