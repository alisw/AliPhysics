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
// class for TOF reconstruction                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliTOFReconstructor.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliTOF.h"
#include "AliTOFtracker.h"


ClassImp(AliTOFReconstructor)


//_____________________________________________________________________________
void AliTOFReconstructor::Reconstruct(AliRunLoader* /*runLoader*/) const
{
// nothing to be done

}

//_____________________________________________________________________________
AliTracker* AliTOFReconstructor::CreateTracker(AliRunLoader* runLoader) const
{
// create a TOF tracker

  AliTOFGeometry* geom = GetTOFGeometry(runLoader);
  if (!geom) return NULL;
  Double_t parPID[] = {130., 5.};
  return new AliTOFtracker(geom, parPID);
}

//_____________________________________________________________________________
void AliTOFReconstructor::FillESD(AliRunLoader* /*runLoader*/, 
				  AliESD* /*esd*/) const
{
// nothing to be done

}


//_____________________________________________________________________________
AliTOFGeometry* AliTOFReconstructor::GetTOFGeometry(AliRunLoader* runLoader) const
{
// get the TOF parameters

  if (!runLoader->GetAliRun()) runLoader->LoadgAlice();
  if (!runLoader->GetAliRun()) {
    Error("GetTOFGeometry", "couldn't get AliRun object");
    return NULL;
  }
  AliTOF* tof = (AliTOF*) runLoader->GetAliRun()->GetDetector("TOF");
  if (!tof) {
    Error("GetTOFGeometry", "couldn't get TOF detector");
    return NULL;
  }
  if (!tof->GetGeometry()) {
    Error("GetTOFGeometry", "no TOF geometry available");
    return NULL;
  }
  return tof->GetGeometry();
}
