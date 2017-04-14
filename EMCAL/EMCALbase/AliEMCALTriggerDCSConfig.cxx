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

#include "AliEMCALTriggerDCSConfig.h"
#include "AliEMCALTriggerSTUDCSConfig.h"
#include "AliEMCALTriggerTRUDCSConfig.h"

/// \cond CLASSIMP
ClassImp(AliEMCALTriggerDCSConfig) ;
/// \endcond

///
/// Default constructor
//_____________________________________________________________________________
AliEMCALTriggerDCSConfig::AliEMCALTriggerDCSConfig() : TObject()
,fTRUArr(0x0)
,fSTUObj(0x0)
,fSTUDCAL(0x0)
{
  /*
	fTRUArr = new TClonesArray("AliEMCALTriggerTRUDCSConfig",62);
	fSTUObj = new AliEMCALTriggerSTUDCSConfig();
	*/
}

///
/// Destructor
//_____________________________________________________________________________
AliEMCALTriggerDCSConfig::~AliEMCALTriggerDCSConfig()
{	
  delete fTRUArr; fTRUArr = 0x0;
  delete fSTUObj; fSTUObj = 0x0;
}
