

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Oystein Djuvsland <oysteind@ift.uib.no>                       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTCaloClusterizerNbyN.h"
#include "AliHLTLogging.h"

/**
 * @file   AliHLTCaloClusterizerNbyN.cxx
 * @author Oystein Djuvsland
 * @date
 * @brief  Clusterizer for CALO HLT
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
 

AliHLTCaloClusterizerNbyN::AliHLTCaloClusterizerNbyN(TString det) : AliHLTCaloClusterizer(det)
,fN(3)
{
// Constructor
}
AliHLTCaloClusterizerNbyN::~AliHLTCaloClusterizerNbyN()
{
// Destructor
}

Int_t AliHLTCaloClusterizerNbyN::ClusterizeEvent(Int_t nDigits)
{
  HLTFatal("Clusterizer deprecated");
  return 0;
}

