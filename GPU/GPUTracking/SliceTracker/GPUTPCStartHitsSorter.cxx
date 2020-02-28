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

/// \file GPUTPCStartHitsSorter.cxx
/// \author David Rohr

#include "GPUTPCStartHitsSorter.h"
#include "GPUTPCTracker.h"

using namespace GPUCA_NAMESPACE::gpu;

template <>
GPUdii() void GPUTPCStartHitsSorter::Thread<0>(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() MEM_LOCAL(GPUSharedMemory) & GPUrestrict() s, processorType& GPUrestrict() tracker)
{
  // Sorts the Start Hits by Row Index
  if (iThread == 0) {
    const int tmpNRows = GPUCA_ROW_COUNT - 6;
    const int nRows = iBlock == (nBlocks - 1) ? (tmpNRows - (tmpNRows / nBlocks) * (nBlocks - 1)) : (tmpNRows / nBlocks);
    const int nStartRow = (tmpNRows / nBlocks) * iBlock + 1;
    int startOffset2 = 0;

    for (int ir = 1; ir < GPUCA_ROW_COUNT - 5; ir++) {
      if (ir < nStartRow) {
        startOffset2 += tracker.RowStartHitCountOffset()[ir];
      }
    }
    s.mStartOffset = startOffset2;
    s.mNRows = nRows;
    s.mStartRow = nStartRow;
  }
  GPUbarrier();

  int startOffset = s.mStartOffset;
  for (int ir = 0; ir < s.mNRows; ir++) {
    GPUglobalref() GPUTPCHitId* const GPUrestrict() startHits = tracker.TrackletStartHits();
    GPUglobalref() GPUTPCHitId* const GPUrestrict() tmpStartHits = tracker.TrackletTmpStartHits() + (s.mStartRow + ir) * tracker.NMaxRowStartHits();
    const int tmpLen = tracker.RowStartHitCountOffset()[ir + s.mStartRow]; // Length of hits in row stored by StartHitsFinder

    for (int j = iThread; j < tmpLen; j += nThreads) {
      startHits[startOffset + j] = tmpStartHits[j];
    }
    startOffset += tmpLen;
  }
}
