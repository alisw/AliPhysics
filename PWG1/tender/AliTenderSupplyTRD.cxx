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

////////////////////////////////////////////////////////////////////
// Correct ESD-TRD info based on better calibration
// - PID recalculation for new gain calibrations
//
// Author : Alexandru Bercuci <A.Bercuci@gsi.de>
//
////////////////////////////////////////////////////////////////////

#include "AliTenderSupplyTRD.h"


ClassImp(AliTenderSupplyTRD)

//_________________________________________________________________
AliTenderSupplyTRD::AliTenderSupplyTRD()
  : AliTenderSupply("TRD", NULL)
{
// Default constructor
}

//_________________________________________________________________
void AliTenderSupplyTRD::Init()
{
// Link TRD specific info

}

//_________________________________________________________________
void AliTenderSupplyTRD::ProcessEvent()
{
// Apply the following corrections
//  - PID
//  - ...


}

