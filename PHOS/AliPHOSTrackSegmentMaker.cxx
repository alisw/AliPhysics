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

//_________________________________________________________________________
// A brief description of the class
//*-- Author : Yves Schutz  SUBATECH 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

#include "TObjArray.h"
#include "TClonesArray.h"

// --- Standard library ---

#include <iostream>

// --- AliRoot header files ---

#include "AliPHOSTrackSegmentMaker.h"
#include "AliPHOSTrackSegment.h"
#include "AliPHOSLink.h"
#include "AliPHOSv0.h"
#include "AliRun.h"

ClassImp( AliPHOSTrackSegmentMaker) 


//____________________________________________________________________________
 AliPHOSTrackSegmentMaker:: AliPHOSTrackSegmentMaker() // ctor
{
  fR0 = 4. ;   
}

//____________________________________________________________________________
void  AliPHOSTrackSegmentMaker::MakeTrackSegments(DigitsList * DL, RecPointsList * emcl, RecPointsList * ppsdl,
                                        TrackSegmentsList * trsl)
{
 
}

