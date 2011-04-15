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
 Author: Jiri Kral, JYU
*/

#include "AliEMCALTriggerTRUDCSConfig.h"

ClassImp(AliEMCALTriggerTRUDCSConfig)
  
//_____________________________________________________________________________
AliEMCALTriggerTRUDCSConfig::AliEMCALTriggerTRUDCSConfig() : TObject()
,fSELPF(0)
,fL0SEL(0)
,fL0COSM(0)
,fGTHRL0(0)
{
	//
	// AliEMCALTriggerTRUDCSConfig default constructor
	//
	for (Int_t i=0;i<6;i++) fMaskReg[i] = 0;
}


