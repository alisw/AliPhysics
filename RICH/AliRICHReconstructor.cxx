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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for RICH reconstruction                                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliRICHReconstructor.h"
#include "AliRICHClusterFinder.h"
#include "AliRICHHelix.h"
#include <AliRunLoader.h>
#include <AliRun.h>
#include <AliESD.h>

ClassImp(AliRICHReconstructor)

//__________________________________________________________________________________________________
void AliRICHReconstructor::Reconstruct(AliRunLoader* pAL) const
{
//Finds clusters out of digits
  AliDebug(1,"Start cluster finder.");AliRICHClusterFinder clus(GetRICH(pAL));  clus.Exec();
}
//__________________________________________________________________________________________________
void AliRICHReconstructor::FillESD(AliRunLoader* /*pAL*/, AliESD* /*pESD*/) const
{
//All is done in AliRICHTracker. Nothing to do here.
}//FillESD
//__________________________________________________________________________________________________
AliRICH* AliRICHReconstructor::GetRICH(AliRunLoader* pAL) const
{
// get the RICH detector

  if (!pAL->GetAliRun()) pAL->LoadgAlice();
  if (!pAL->GetAliRun()) {AliError("couldn't get AliRun object"); return NULL;  }
  AliRICH* pRich = (AliRICH*) pAL->GetAliRun()->GetDetector("RICH");
  if (!pRich) {AliError("couldn't get RICH detector");    return NULL;  }
  return pRich;
}
//__________________________________________________________________________________________________
