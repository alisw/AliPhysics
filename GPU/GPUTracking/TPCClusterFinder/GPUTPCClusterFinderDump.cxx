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

/// \file GPUTPCClusterFinderDump.cxx
/// \author David Rohr

#include "GPUTPCClusterFinder.h"
#include "GPUReconstruction.h"
#include "ClusterNative.h"
#include "Digit.h"

using namespace GPUCA_NAMESPACE::gpu;

void GPUTPCClusterFinder::DumpDigits(std::ostream& out)
{
  out << "Clusterer - Digits - Slice " << mISlice << ": " << mPmemory->counters.nDigits << "\n";
  for (size_t i = 0; i < mPmemory->counters.nDigits; i++) {
    out << i << ": " << mPdigits[i].charge << ", " << mPdigits[i].time << ", " << (int)mPdigits[i].pad << ", " << (int)mPdigits[i].row << "\n";
  }
}

void GPUTPCClusterFinder::DumpChargeMap(std::ostream& out, std::string_view title)
{
  out << "Clusterer - " << title << " - Slice " << mISlice << "\n";
  for (unsigned int i = 0; i < TPC_MAX_TIME_PADDED; i++) {
    out << "Line " << i;
    for (unsigned int j = 0; j < TPC_NUM_OF_PADS; j++) {
      if (mPchargeMap[i * TPC_NUM_OF_PADS + j]) {
        out << " " << std::hex << mPchargeMap[i * TPC_NUM_OF_PADS + j] << std::dec;
      }
    }
    out << "\n";
  }
}

void GPUTPCClusterFinder::DumpPeaks(std::ostream& out)
{
  out << "Clusterer - Peaks - Slice " << mISlice << "\n";
  for (unsigned int i = 0; i < mPmemory->counters.nDigits; i++) {
    out << (int)mPisPeak[i] << " ";
    if ((i + 1) % 100 == 0) {
      out << "\n";
    }
  }
  out << "\n";
}

void GPUTPCClusterFinder::DumpPeaksCompacted(std::ostream& out)
{
  out << "Clusterer - Compacted Peaks - Slice " << mISlice << ": " << mPmemory->counters.nPeaks << "\n";
  for (size_t i = 0; i < mPmemory->counters.nPeaks; i++) {
    out << i << ": " << mPpeaks[i].charge << ", " << mPpeaks[i].time << ", " << (int)mPpeaks[i].pad << ", " << (int)mPpeaks[i].row << "\n";
  }
}

void GPUTPCClusterFinder::DumpSuppressedPeaks(std::ostream& out)
{
  out << "Clusterer - NoiseSuppression - Slice " << mISlice << "\n";
  for (unsigned int i = 0; i < mPmemory->counters.nPeaks; i++) {
    out << (int)mPisPeak[i] << " ";
    if ((i + 1) % 100 == 0) {
      out << "\n";
    }
  }
  out << "\n";
}

void GPUTPCClusterFinder::DumpSuppressedPeaksCompacted(std::ostream& out)
{
  out << "Clusterer - Noise Suppression Peaks Compacted - Slice " << mISlice << ": " << mPmemory->counters.nClusters << "\n";
  for (size_t i = 0; i < mPmemory->counters.nClusters; i++) {
    out << i << ": " << mPfilteredPeaks[i].charge << ", " << mPfilteredPeaks[i].time << ", " << (int)mPfilteredPeaks[i].pad << ", " << (int)mPfilteredPeaks[i].row << "\n";
  }
}

void GPUTPCClusterFinder::DumpCountedPeaks(std::ostream& out)
{
  out << "Clusterer - Peak Counts - Slice " << mISlice << "\n";
  for (int i = 0; i < GPUCA_ROW_COUNT; i++) {
    out << i << ": " << mPclusterInRow[i] << "\n";
  }
}

void GPUTPCClusterFinder::DumpClusters(std::ostream& out)
{
  out << "Clusterer - Clusters - Slice " << mISlice << "\n";

  for (int i = 0; i < GPUCA_ROW_COUNT; i++) {
    size_t N = mPclusterInRow[i];
    out << "Row: " << i << ": " << N << "\n";
    std::vector<deprecated::ClusterNative> sortedCluster(N);
    deprecated::ClusterNative* row = &mPclusterByRow[i * mNMaxClusterPerRow];
    std::copy(row, &row[N], sortedCluster.begin());

    std::sort(sortedCluster.begin(), sortedCluster.end(), [](const auto& c1, const auto& c2) {
      float t1 = deprecated::cnGetTime(&c1);
      float t2 = deprecated::cnGetTime(&c2);
      float p1 = deprecated::cnGetPad(&c1);
      float p2 = deprecated::cnGetPad(&c2);
      uchar f1 = deprecated::cnGetFlags(&c1);
      uchar f2 = deprecated::cnGetFlags(&c2);

      if (t1 < t2) {
        return true;
      } else if (t1 == t2) {
        if (p1 < p2) {
          return true;
        } else if (p1 == p2) {
          return f1 < f2;
        } else {
          return false;
        }
      } else {
        return false;
      }
    });

    for (const auto& cl : sortedCluster) {
      out << std::hex << cl.timeFlagsPacked << std::dec << ", " << cl.padPacked << ", " << (int)cl.sigmaTimePacked << ", " << (int)cl.sigmaPadPacked << ", " << cl.qmax << ", " << cl.qtot << "\n";
    }
  }
}
