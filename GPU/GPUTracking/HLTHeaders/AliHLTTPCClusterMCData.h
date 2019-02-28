//**************************************************************************\
//* This file is property of and copyright by the ALICE Project            *\
//* ALICE Experiment at CERN, All rights reserved.                         *\
//*                                                                        *\
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *\
//*                  for The ALICE HLT Project.                            *\
//*                                                                        *\
//* Permission to use, copy, modify and distribute this software and its   *\
//* documentation strictly for non-commercial purposes is hereby granted   *\
//* without fee, provided that the above copyright notice appears in all   *\
//* copies and that both the copyright notice and this permission notice   *\
//* appear in the supporting documentation. The authors make no claims     *\
//* about the suitability of this software for any purpose. It is          *\
//* provided "as is" without express or implied warranty.                  *\
//**************************************************************************

/// \file AliHLTTPCClusterMCData.h
/// \author ALICE HLT Project

#ifndef _ALIHLTTPCCLUSTERMCDATA_H_
#define _ALIHLTTPCCLUSTERMCDATA_H_

/**
 * @struct AliHLTTPCClusterMCWeight
 * This in a struct for MC weights
 * @ingroup alihlt_tpc
 */
struct AliHLTTPCClusterMCWeight {
  //* constructor **/
  AliHLTTPCClusterMCWeight() : fMCID(-1), fWeight(0) {}

  int fMCID;     // MC track ID
  float fWeight; // weight of the track ID
};

typedef struct AliHLTTPCClusterMCWeight AliHLTTPCClusterMCWeight;

/**
 * @struct AliHLTTPCClusterMCLabel
 * This in a struct for MC labels
 * @ingroup alihlt_tpc
 */
struct AliHLTTPCClusterMCLabel {
  AliHLTTPCClusterMCWeight fClusterID[3]; // three most relevant MC labels
};

typedef struct AliHLTTPCClusterMCLabel AliHLTTPCClusterMCLabel;

/**
 * @struct AliHLTTPCClusterMCData
 * This in a container for MC labels
 * @ingroup alihlt_tpc
 */
struct AliHLTTPCClusterMCData {
  int fCount;
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
  AliHLTTPCClusterMCLabel fLabels[1];
#else
  AliHLTTPCClusterMCLabel fLabels[0];
#endif
};

typedef struct AliHLTTPCClusterMCData AliHLTTPCClusterMCData;

#endif
