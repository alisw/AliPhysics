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
#include <TVector3.h>
#include "AliReducedMatchedTrack.h"

/// \cond CLASSIMP
ClassImp(HighPtTracks::AliReducedMatchedTrack);
/// \endcond

namespace HighPtTracks {

/**
 * Dummy (I/O) constructor, not to be used
 */
AliReducedMatchedTrack::AliReducedMatchedTrack() :
    TObject(),
    fPx(0),
    fPy(0),
    fPz(0),
    fGoodTrackLabel(false),
    fNclustersTPC(0),
    fTrackCuts(0)
{
}

/**
 * Main constructor, initializing 3-vector of the momentum
 *
 * \param px x-component of the 3-vector
 * \param py y-component of the 3-vector
 * \param pz z-component of the 3-vector
 */
AliReducedMatchedTrack::AliReducedMatchedTrack(double px, double py, double pz):
    TObject(),
    fPx(px),
    fPy(py),
    fPz(pz),
    fGoodTrackLabel(false),
    fNclustersTPC(0),
    fTrackCuts(0)
{
}

/**
 *  Get track \f$ p_{t} \f$ from the track 3-momentum
 *
 * \return track \f$ p_{t} \f$
 */
double AliReducedMatchedTrack::Pt() const {
  TVector3 pvec(fPx, fPy, fPz);
  return pvec.Pt();
}

/**
 *  Get track \f$ \eta \f$ from the track 3-momentum
 *
 * \return track \f$ \eta \f$
 */
double AliReducedMatchedTrack::Eta() const {
  TVector3 pvec(fPx, fPy, fPz);
  return pvec.Eta();
}

/**
 *  Get track \f$ \phi \f$ from the track 3-momentum
 *
 * \return track \f$ \phi \f$
 */
double AliReducedMatchedTrack::Phi() const {
  TVector3 pvec(fPx, fPy, fPz);
  return pvec.Phi();
}

/**
 * Fill 3-vector with momentum information
 *
 * \param vec The vector to be filled
 */
void AliReducedMatchedTrack::FillVector(TVector3& vec) {
  vec.SetXYZ(fPx, fPy, fPz);
}

} /* namespace HighPtTracks */
