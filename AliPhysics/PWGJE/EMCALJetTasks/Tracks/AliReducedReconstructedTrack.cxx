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
   TObject(),
   fCharge(0),
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
void AliReducedReconstructedTrack::FillMomentumVector(TVector3& pvec) const {
  pvec.SetXYZ(fPVec[0], fPVec[1], fPVec[2]);
}
/**
 * Get reconstructed track \f$ p_{t} \f$
 * \return The reconstructed \f$ p_{t} \f$
 */
Double_t AliReducedReconstructedTrack::Pt() const {
  TVector3 pvec;
  FillMomentumVector(pvec);
  return pvec.Pt();
}

/**
 * Get reconstructed track \f$ \eta \f$
 * \return The reconstructed \f$ \eta \f$
 */
Double_t AliReducedReconstructedTrack::Eta() const {
  TVector3 pvec;
  FillMomentumVector(pvec);
  return pvec.Eta();
}

/**
 * Get reconstructed track \f$ \phi \f$
 * \return The reconstructed \f$ \phi \f$
 */
Double_t AliReducedReconstructedTrack::Phi() const {
  TVector3 pvec;
  FillMomentumVector(pvec);
  return pvec.Phi();
}

} /* namespace HighPtTracks */


