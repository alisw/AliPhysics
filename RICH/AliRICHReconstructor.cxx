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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for RICH reconstruction                                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliRICHReconstructor.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliRICHClusterFinder.h"


ClassImp(AliRICHReconstructor)


//_____________________________________________________________________________
void AliRICHReconstructor::Reconstruct(AliRunLoader* runLoader) const
{
// reconstruct clusters

  AliRICH* rich = GetRICH(runLoader);
  if (!rich) return;
  AliRICHClusterFinder clusterer(rich);
  clusterer.Exec();
}

//_____________________________________________________________________________
void AliRICHReconstructor::FillESD(AliRunLoader* /*runLoader*/, 
				   AliESD* /*esd*/) const
{
// nothing to be done

}


//_____________________________________________________________________________
AliRICH* AliRICHReconstructor::GetRICH(AliRunLoader* runLoader) const
{
// get the RICH detector

  if (!runLoader->GetAliRun()) runLoader->LoadgAlice();
  if (!runLoader->GetAliRun()) {
    Error("GetRICH", "couldn't get AliRun object");
    return NULL;
  }
  AliRICH* rich = (AliRICH*) runLoader->GetAliRun()->GetDetector("RICH");
  if (!rich) {
    Error("GetRICH", "couldn't get RICH detector");
    return NULL;
  }
  return rich;
}
