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
//_________________________________________________________________________
// Algorithm Base class to construct EMCAL track segments
// Associates EMC and PPSD clusters
// Unfolds the EMC cluster   
//*-- 
//*-- Author: Dmitri Peressounko (RRC Ki & SUBATECH)
//             Adapted from PHOS by Y. Schutz (SUBATECH)


// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---
#include "AliEMCALTrackSegmentMaker.h"
#include "AliEMCALGetter.h"

ClassImp( AliEMCALTrackSegmentMaker) 


//____________________________________________________________________________
  AliEMCALTrackSegmentMaker:: AliEMCALTrackSegmentMaker() : TTask("","")
{
  // ctor
  fEventFolderName = "" ; 

}

//____________________________________________________________________________
AliEMCALTrackSegmentMaker::AliEMCALTrackSegmentMaker(const TString alirunFileName,
						     const TString eventFolderName):
  TTask("EMCAL"+AliConfig::Instance()->GetTrackerTaskName(), alirunFileName),
  fEventFolderName(eventFolderName)
{
  // ctor

}

//____________________________________________________________________________
AliEMCALTrackSegmentMaker::~AliEMCALTrackSegmentMaker()
{
 //Remove this from the parental task before destroying
  AliEMCALGetter::Instance()->EmcalLoader()->CleanTracker();
}

