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

/// \file GPUTPCStartHitsFinder.cxx
/// \author Sergey Gorbunov, Ivan Kisel, David Rohr

#include "GPUTPCStartHitsFinder.h"
#include "GPUTPCTracker.h"
#include "GPUCommonMath.h"

using namespace GPUCA_NAMESPACE::gpu;

template <>
GPUdii() void GPUTPCStartHitsFinder::Thread<0>(int /*nBlocks*/, int nThreads, int iBlock, int iThread, GPUsharedref() MEM_LOCAL(GPUSharedMemory) & GPUrestrict() s, processorType& GPUrestrict() tracker)
{
  // find start hits for tracklets

  if (iThread == 0) {
    s.mIRow = iBlock + 1;
    s.mNRowStartHits = 0;
    if (s.mIRow <= GPUCA_ROW_COUNT - 4) {
      s.mNHits = tracker.Row(s.mIRow).NHits();
    } else {
      s.mNHits = -1;
    }
  }
  GPUbarrier();
  GPUglobalref() const MEM_GLOBAL(GPUTPCRow) & GPUrestrict() row = tracker.Row(s.mIRow);
  GPUglobalref() const MEM_GLOBAL(GPUTPCRow) & GPUrestrict() rowUp = tracker.Row(s.mIRow + 2);
  for (int ih = iThread; ih < s.mNHits; ih += nThreads) {
    if (tracker.HitLinkDownData(row, ih) == CALINK_INVAL && tracker.HitLinkUpData(row, ih) != CALINK_INVAL && tracker.HitLinkUpData(rowUp, tracker.HitLinkUpData(row, ih)) != CALINK_INVAL) {
#ifdef GPUCA_SORT_STARTHITS
      GPUglobalref() GPUTPCHitId* const GPUrestrict() startHits = tracker.TrackletTmpStartHits() + s.mIRow * tracker.NMaxRowStartHits();
      unsigned int nextRowStartHits = CAMath::AtomicAddShared(&s.mNRowStartHits, 1);
      CONSTEXPR int errCode = GPUCA_ERROR_ROWSTARTHIT_OVERFLOW;
      if (nextRowStartHits + 1 >= tracker.NMaxRowStartHits())
#else
      GPUglobalref() GPUTPCHitId* const GPUrestrict() startHits = tracker.TrackletStartHits();
      unsigned int nextRowStartHits = CAMath::AtomicAdd(tracker.NTracklets(), 1);
      CONSTEXPR int errCode = GPUCA_ERROR_TRACKLET_OVERFLOW;
      if (nextRowStartHits + 1 >= tracker.NMaxStartHits())
#endif
      {
        tracker.CommonMemory()->kernelError = errCode;
        CAMath::AtomicExch(tracker.NTracklets(), 0);
        break;
      }
      startHits[nextRowStartHits].Set(s.mIRow, ih);
    }
  }
  GPUbarrier();

#ifdef GPUCA_SORT_STARTHITS
  if (iThread == 0) {
    unsigned int nOffset = CAMath::AtomicAdd(tracker.NTracklets(), s.mNRowStartHits);
    tracker.RowStartHitCountOffset()[s.mIRow] = s.mNRowStartHits;
    if (nOffset + s.mNRowStartHits >= tracker.NMaxStartHits()) {
      tracker.CommonMemory()->kernelError = GPUCA_ERROR_TRACKLET_OVERFLOW;
      CAMath::AtomicExch(tracker.NTracklets(), 0);
    }
  }
#endif
}
