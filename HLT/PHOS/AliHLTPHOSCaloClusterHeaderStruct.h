
#ifndef ALIHLTCALOCLUSTERHEADERSTRUCT_H
#define ALIHLTCALOCLUSTERHEADERSTRUCT_H

/**************************************************************************
 * Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved.      *
 *                                                                        *
 * Author: Per Thomas Hille <perthi@fys.uio.no> for the ALICE HLT Project.*
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

/**
 * Calo cluster struct for  HLT
 *
 * @file   AliHLTCaloClusterHeaderStruct
 * @author Oystein Djuvsland
 * @date
 * @brief  Calorimeter cluster header for HLT
 */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

/**
 * @struct AliHLTCaloClusterHeaderStruct
 * Calorimeter cluster header describing the number of 
 * clusters in the following block
 *
 * @ingroup alihlt_phos
 */

#include "Rtypes.h"

struct AliHLTCaloClusterHeaderStruct
{
  Short_t fNClusters;
};

#endif
