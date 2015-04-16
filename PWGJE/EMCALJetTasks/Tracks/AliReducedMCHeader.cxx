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
#include "AliReducedMCHeader.h"

/// \cond CLASSIMP
ClassImp(HighPtTracks::AliReducedMCHeader)
/// \endcond

namespace HighPtTracks {

/**
 * Dummy constructor
 */
AliReducedMCHeader::AliReducedMCHeader():
  TObject(),
  fCrossSection(-1),
  fNumberOfTrials(-1),
  fPtHard(-1)
{
}

/**
 * Constructor, initialising also the parameter of the event header
 * \param crosssection The event cross section
 * \param numberOfTrials The number of trials
 * \param pthard The \f$ p_{t} \f$ of the hard interaction
 */
AliReducedMCHeader::AliReducedMCHeader(Double_t crosssection,Double_t numberOfTrials, Double_t pthard):
  TObject(),
  fCrossSection(crosssection),
  fNumberOfTrials(numberOfTrials),
  fPtHard(pthard)
{
}

/**
 * Destructor, nothing to do
 */
AliReducedMCHeader::~AliReducedMCHeader() {}

} /* namespace HighPtTracks */
