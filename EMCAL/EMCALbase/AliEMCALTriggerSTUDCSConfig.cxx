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

#include "AliEMCALTriggerSTUDCSConfig.h"
#include "TVector2.h"

ClassImp(AliEMCALTriggerSTUDCSConfig)
  
//_____________________________________________________________________________
AliEMCALTriggerSTUDCSConfig::AliEMCALTriggerSTUDCSConfig() : TObject(),
fGetRawData(1),
fRegion(0xFFFFFFFF),
fFw(0x2A012)
{
	//
	// AliEMCALTriggerSTUDCSConfig default constructor
	//
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 2; j++) {
			fG[i][j] = 0;
			fJ[i][j] = 0;
		}	
	}
	memset(fPHOSScale, 0, sizeof(Int_t) * 32);
}

//_____________________________________________________________________________
void AliEMCALTriggerSTUDCSConfig::GetSegmentation(TVector2& v1, TVector2& v2, TVector2& v3, TVector2& v4) const
{
	// Get Segmentation
	
	v1.Set(1., 1.);
	v2.Set(2., 2.);
	v3.Set(4., 4.);
	
	Double_t js = 2 + (fFw >> 16);
	v4.Set(js, js);
}

