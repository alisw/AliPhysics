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

ClassImp(AliEMCALTriggerSTUDCSConfig)
  
//_____________________________________________________________________________
AliEMCALTriggerSTUDCSConfig::AliEMCALTriggerSTUDCSConfig() : TObject()
  ,fGA(0)
  ,fGB(1)
  ,fGC(0)
  ,fJA(0)
  ,fJB(1)
  ,fJC(0)
  ,fGetRawData(0)
  ,fRegion(0xFFFFFFFF)
  ,fFw(2223)
{
  //
  // AliEMCALTriggerSTUDCSConfig default constructor
  //
}


