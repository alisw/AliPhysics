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

/// \file ClusterNativeAccessExt.h
/// \author David Rohr

#ifndef CLUSTERNATIVEACCESSEXT_H
#define CLUSTERNATIVEACCESSEXT_H

#include "GPUTPCDef.h"

#ifdef HAVE_O2HEADERS
#include "DataFormatsTPC/ClusterNative.h"
#else
namespace o2
{
namespace tpc
{
struct ClusterNative {
};
struct ClusterNativeAccessFullTPC {
  const ClusterNative* clusters[GPUCA_NSLICES][GPUCA_ROW_COUNT];
  unsigned int nClusters[GPUCA_NSLICES][GPUCA_ROW_COUNT];
};
struct Constants {
  static constexpr int MAXSECTOR = GPUCA_NSLICES;
  static constexpr int MAXGLOBALPADROW = GPUCA_ROW_COUNT;
};
} // namespace tpc
} // namespace o2
#endif

namespace GPUCA_NAMESPACE
{
namespace gpu
{
struct ClusterNativeAccessExt : public o2::tpc::ClusterNativeAccessFullTPC {
  unsigned int clusterOffset[GPUCA_NSLICES][GPUCA_ROW_COUNT];
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
