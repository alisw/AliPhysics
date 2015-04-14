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
#include "AliReducedEmcalCluster.h"

/// \cond CLASSIMP
ClassImp(HighPtTracks::AliReducedEmcalCluster)
/// \endcond

namespace HighPtTracks {

/**
 * Constructor, initialising with default values
 */
AliReducedEmcalCluster::AliReducedEmcalCluster() :
  fClusterID(-1),
  fEnergy(-1),
  fEta(-100),
  fPhi(-100),
  fM02(-1),
  fM20(-1)
{
}

/**
 * Constructor, initialising cluster with all parameters
 * \param id ID of the cluster
 * \param energy Energy of the cluster
 * \param eta \f$ \eta \f$ position of the cluster
 * \param phi \f$ \phi \f$ position of the cluster
 * \param m02 The m02 shower shape parameter
 * \param m20 The m20 shower shape parameter
 */
AliReducedEmcalCluster::AliReducedEmcalCluster(Int_t id, Float_t energy,Float_t eta, Float_t phi, Float_t m02, Float_t m20):
  fClusterID(id),
  fEnergy(energy),
  fEta(eta),
  fPhi(phi),
  fM02(m02),
  fM20(m20)
{
}

/**
 * Destructor, nothing to do
 */
AliReducedEmcalCluster::~AliReducedEmcalCluster() { }

} /* namespace HighPtTracks */
