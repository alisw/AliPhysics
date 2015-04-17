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
#include "AliReducedPatchInfo.h"

/// \cond CLASSIMP
ClassImp(HighPtTracks::AliReducedPatchInfo)
/// \endcond

namespace HighPtTracks {

  /**
   * Constructor, initialising with default values
   */
  AliReducedPatchInfo::AliReducedPatchInfo():
    TObject(),
    fEnergy(-1),
    fAmplitude(-1),
    fEta(-100),
    fPhi(-1)
  {}
  /**
   * Constructor, initialising patch energy and position
   * \param energy Patch energy
   * \param amplitude Patch amplitude
   * \param eta Patch \f$ \eta \f$ (center)
   * \param phi Patch \f$ \phi \f$ (center)
   */
  AliReducedPatchInfo::AliReducedPatchInfo(Float_t energy, Float_t amplitude, Float_t eta, Float_t phi):
    TObject(),
    fEnergy(energy),
    fAmplitude(amplitude),
    fEta(eta),
    fPhi(phi)
  {}

} /* namespace HighPtTracks */
