/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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
#include <TArrayI.h>
#include <TList.h>
#include <TVector3.h>

#include "AliReducedReconstructedTrack.h"

/// \cond CLASSIMP
ClassImp(HighPtTracks::AliReducedReconstructedTrack)
/// \endcond

namespace HighPtTracks {

/**
 * Dummy constructor
 */
AliReducedReconstructedTrack::AliReducedReconstructedTrack() :
   fTrackCutsMap(0),
   fClusterIndex(-1),
   fParticleIndex(-1),
   fGoodMCTrack(kFALSE),
   fClustersTPC(0),
   fCrossedRowsTPC(0),
   fSharedClustersTPC(0),
   fFindableClustersTPC(0)
{
  memset(fPVec, 0, sizeof(Double_t)*3);
}

/**
 * Destructor
 */
AliReducedReconstructedTrack::~AliReducedReconstructedTrack() {
}

/**
 * Fill a TVector3 with the track momentum information
 * \param pvec The vector to be filled
 */
void AliReducedReconstructedTrack::FillMomentumVector(TVector3& pvec) {
  pvec.SetXYZ(pvec[0], pvec[1], pvec[2]);
}

} /* namespace HighPtTracks */


