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
  double factor = 1;

  // Offset
  static constexpr double offset = 1000.;
  static constexpr double hitOffset = 20000;

  // Scaling Factors
  static constexpr double tpcPeaksPerDigit = 0.2;
  static constexpr double tpcClustersPerPeak = 0.9;
  static constexpr double tpcStartHitsPerHit = 0.08;
  static constexpr double tpcTrackletsPerStartHit = 0.8;
  static constexpr double tpcTrackletHitsPerHit = 5;
  static constexpr double tpcSectorTracksPerHit = 0.02;
  static constexpr double tpcSectorTrackHitsPerHit = 0.8f;
  static constexpr double tpcTracksPerHit = 0.012;
  static constexpr double tpcTrackHitsPerHit = 0.7;

  double NTPCPeaks(double tpcDigits) { return hitOffset + tpcDigits * tpcPeaksPerDigit * factor; }
  double NTPCClusters(double tpcDigits) { return tpcClustersPerPeak * NTPCPeaks(tpcDigits) * factor; }
  double NTPCStartHits(double tpcHits) { return offset + tpcHits * tpcStartHitsPerHit * factor; }
  double NTPCTracklets(double tpcHits) { return NTPCStartHits(tpcHits) * tpcTrackletsPerStartHit * factor; }
  double NTPCTrackletHits(double tpcHits) { return hitOffset + tpcHits * tpcTrackletHitsPerHit * factor; }
  double NTPCSectorTracks(double tpcHits) { return offset + tpcHits * tpcSectorTracksPerHit * factor; }
  double NTPCSectorTrackHits(double tpcHits) { return offset + tpcHits * tpcSectorTrackHitsPerHit * factor; }
  double NTPCTracks(double tpcHits) { return offset + tpcHits * tpcTracksPerHit * factor; }
  double NTPCTrackHits(double tpcHits) { return offset + tpcHits * tpcTrackHitsPerHit * factor; }
  double NTPCMaxRowStartHits(double tpcHits) { return offset + NTPCStartHits(tpcHits) / GPUCA_ROW_COUNT * 3. * factor; }
};

} // namespace GPUCA_NAMESPACE::gpu

#endif
