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

/// \file GPUMemorySizeScalers.h
/// \author David Rohr

#ifndef O2_GPU_GPUMEMORYSIZESCALERS_H
#define O2_GPU_GPUMEMORYSIZESCALERS_H

#include "GPUDef.h"

namespace GPUCA_NAMESPACE::gpu
{

struct GPUMemorySizeScalers {
  // Input sizes
  double nTPCdigits = 0;
  double nTPCHits = 0;
  double nTRDTracklets = 0;
  double nITSTracks = 0;

  // Offset
  static constexpr double offset = 1000.;

  // Scaling Factors
  static constexpr double tpcClustersPerDigit = 0.2;
  static constexpr double tpcClustersNoiseOffset = 20000;
  static constexpr double tpcTrackletsPerHit = 0.08;
  static constexpr double tpcSectorTracksPerHit = 0.02;
  static constexpr double tpcSectorTrackHitsPerHit = 0.8f;
  static constexpr double tpcTracksPerHit = 0.012;
  static constexpr double tpcTrackHitsPerHit = 0.7;

  double NTPCClusters(double tpcDigits) { return offset + tpcDigits * tpcClustersPerDigit + tpcClustersNoiseOffset; }
  double NTPCTracklets(double tpcHits) { return offset + tpcHits * tpcTrackletsPerHit; }
  double NTPCSectorTracks(double tpcHits) { return offset + tpcHits * tpcSectorTracksPerHit; }
  double NTPCSectorTrackHits(double tpcHits) { return offset + tpcHits * tpcSectorTrackHitsPerHit; }
  double NTPCTracks(double tpcHits) { return offset + tpcHits * tpcTracksPerHit; }
  double NTPCTrackHits(double tpcHits) { return offset + tpcHits * tpcTrackHitsPerHit; }
  double NTPCMaxRowStartHits(double tpcHits) { return offset + NTPCTracklets(tpcHits) / GPUCA_ROW_COUNT * 3.; }
};

} // namespace GPUCA_NAMESPACE::gpu

#endif
