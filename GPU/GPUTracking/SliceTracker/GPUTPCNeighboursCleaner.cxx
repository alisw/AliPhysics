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

/// \file GPUTPCNeighboursCleaner.cxx
/// \author Sergey Gorbunov, Ivan Kisel, David Rohr

#include "GPUTPCNeighboursCleaner.h"
#include "GPUTPCTracker.h"
#include "GPUCommonMath.h"
using namespace GPUCA_NAMESPACE::gpu;

template <>
GPUdii() void GPUTPCNeighboursCleaner::Thread<0>(int /*nBlocks*/, int nThreads, int iBlock, int iThread, GPUsharedref() MEM_LOCAL(GPUSharedMemory) & s, processorType& tracker)
{
  // *
  // * kill link to the neighbour if the neighbour is not pointed to the cluster
  // *

  if (iThread == 0) {
    s.mIRow = iBlock + 2;
    if (s.mIRow <= GPUCA_ROW_COUNT - 3) {
      s.mIRowUp = s.mIRow + 2;
      s.mIRowDn = s.mIRow - 2;
      s.mNHits = tracker.Row(s.mIRow).NHits();
    }
  }
  GPUbarrier();

  if (s.mIRow <= GPUCA_ROW_COUNT - 3) {
#ifdef GPUCA_GPUCODE
    int Up = s.mIRowUp;
    int Dn = s.mIRowDn;
    GPUglobalref() const MEM_GLOBAL(GPUTPCRow)& row = tracker.Row(s.mIRow);
    GPUglobalref() const MEM_GLOBAL(GPUTPCRow)& rowUp = tracker.Row(Up);
    GPUglobalref() const MEM_GLOBAL(GPUTPCRow)& rowDn = tracker.Row(Dn);
#else
    const GPUTPCRow& row = tracker.Row(s.mIRow);
    const GPUTPCRow& rowUp = tracker.Row(s.mIRowUp);
    const GPUTPCRow& rowDn = tracker.Row(s.mIRowDn);
#endif

    // - look at up link, if it's valid but the down link in the row above doesn't link to us remove
    //   the link
    // - look at down link, if it's valid but the up link in the row below doesn't link to us remove
    //   the link
    for (int ih = iThread; ih < s.mNHits; ih += nThreads) {
      calink up = tracker.HitLinkUpData(row, ih);
      if (up != CALINK_INVAL) {
        calink upDn = tracker.HitLinkDownData(rowUp, up);
        if ((upDn != (calink)ih)) {
          tracker.SetHitLinkUpData(row, ih, CALINK_INVAL);
        }
      }
      calink dn = tracker.HitLinkDownData(row, ih);
      if (dn != CALINK_INVAL) {
        calink dnUp = tracker.HitLinkUpData(rowDn, dn);
        if (dnUp != (calink)ih) {
          tracker.SetHitLinkDownData(row, ih, CALINK_INVAL);
        }
      }
    }
  }
}
