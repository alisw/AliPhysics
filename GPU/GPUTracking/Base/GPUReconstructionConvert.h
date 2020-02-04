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

/// \file GPUReconstructionConvert.h
/// \author David Rohr

#ifndef GPURECONSTRUCTIONCONVERT_H
#define GPURECONSTRUCTIONCONVERT_H

#include <memory>
#include "GPUDef.h"

namespace o2
{
namespace tpc
{
struct ClusterNative;
struct ClusterNativeAccess;
} // namespace tpc
} // namespace o2

class AliHLTTPCRawCluster;

namespace GPUCA_NAMESPACE
{
namespace gpu
{
class GPUParam;
class GPUTPCClusterData;
class TPCFastTransform;
struct GPUTrackingInOutDigits;
struct GPUTrackingInOutZS;
namespace deprecated
{
struct PackedDigit;
}

class GPUReconstructionConvert
{
 public:
  constexpr static unsigned int NSLICES = GPUCA_NSLICES;
  static void ConvertNativeToClusterData(o2::tpc::ClusterNativeAccess* native, std::unique_ptr<GPUTPCClusterData[]>* clusters, unsigned int* nClusters, const TPCFastTransform* transform, int continuousMaxTimeBin = 0);
  static void ConvertRun2RawToNative(o2::tpc::ClusterNativeAccess& native, std::unique_ptr<o2::tpc::ClusterNative[]>& nativeBuffer, const AliHLTTPCRawCluster** rawClusters, unsigned int* nRawClusters);
  static void RunZSEncoder(const GPUTrackingInOutDigits* in, GPUTrackingInOutZS*& out, const GPUParam& param, bool zs12bit);
  static void RunZSFilter(std::unique_ptr<deprecated::PackedDigit[]>* buffers, const deprecated::PackedDigit** ptrs, size_t* ns, const GPUParam& param, bool zs12bit);
  static int GetMaxTimeBin(const o2::tpc::ClusterNativeAccess& native);
  static int GetMaxTimeBin(const GPUTrackingInOutDigits& digits);
  static int GetMaxTimeBin(const GPUTrackingInOutZS& zspages);

 private:
  static void ZSstreamOut(unsigned short* bufIn, unsigned int& lenIn, unsigned char* bufOut, unsigned int& lenOut, unsigned int nBits);
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
