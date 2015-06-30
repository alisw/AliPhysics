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

//-----------------------------------------------------------------------------
// Class AliMFTTrack
//-------------------
// Description of an ALICE Standalone MFT track
//-------------------
// Contact author: raphael.tieulent@cern.ch
//-----------------------------------------------------------------------------


#include "AliMFTTrack.h"

#include "AliMFTCATrack.h"
#include "AliMFTTrackParam.h"

/// \cond CLASSIMP
ClassImp(AliMFTTrack); // Class implementation in ROOT context
/// \endcond


//=============================================================================================

AliMFTTrack::AliMFTTrack():TObject(),
fChi2(0.),
fTrackParamAtVertex(NULL),
fCATrack(NULL)
{
  /// Default constructor
  
}

//=============================================================================================

AliMFTTrack::AliMFTTrack(AliMFTCATrack *catrack):TObject(),
fChi2(0.),
fTrackParamAtVertex(NULL),
fCATrack(catrack)
{
  /// Constructor from a AliMFTCATrack

  
}

//=============================================================================================


AliMFTTrack::~AliMFTTrack() {
  delete fTrackParamAtVertex;
  delete fCATrack;
  
}

