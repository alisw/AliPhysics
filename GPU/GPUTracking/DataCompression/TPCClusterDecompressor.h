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

/// \file TPCClusterDecompressor.h
/// \author David Rohr

#ifndef TPCCLUSTERDECOMPRESSOR_H
#define TPCCLUSTERDECOMPRESSOR_H

#include "GPUTPCCompression.h"
#include <vector>

namespace o2
{
namespace tpc
{
struct ClusterNativeAccess;
struct ClusterNative;
} // namespace tpc
} // namespace o2

namespace GPUCA_NAMESPACE
{
namespace gpu
{
struct GPUParam;

class TPCClusterDecompressor
{
 public:
  static constexpr unsigned int NSLICES = GPUCA_NSLICES;
  int decompress(const o2::tpc::CompressedClustersFlat* clustersCompressed, o2::tpc::ClusterNativeAccess& clustersNative, std::vector<o2::tpc::ClusterNative>& clusterBuffer, const GPUParam& param);
  int decompress(const o2::tpc::CompressedClusters* clustersCompressed, o2::tpc::ClusterNativeAccess& clustersNative, std::vector<o2::tpc::ClusterNative>& clusterBuffer, const GPUParam& param);

 protected:
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
