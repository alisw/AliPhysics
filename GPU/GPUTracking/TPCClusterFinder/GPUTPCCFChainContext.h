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

/// \file GPUTPCCFChainContext.h
/// \author David Rohr

#ifndef O2_GPU_TPCCFCHAINCONTEXT_H
#define O2_GPU_TPCCFCHAINCONTEXT_H

#include "clusterFinderDefs.h"
#include "GPUDataTypes.h"
#include "GPUTPCClusterFinder.h"
#include "CfFragment.h"
#include <vector>
#include <utility>

namespace GPUCA_NAMESPACE
{
namespace gpu
{

struct GPUTPCCFChainContext {
  struct FragmentData {
    unsigned int nDigits[GPUCA_NSLICES][GPUTrackingInOutZS::NENDPOINTS];
    unsigned int nPages[GPUCA_NSLICES][GPUTrackingInOutZS::NENDPOINTS];
    std::vector<unsigned short> pageDigits[GPUCA_NSLICES][GPUTrackingInOutZS::NENDPOINTS];
    GPUTPCClusterFinder::MinMaxCN minMaxCN[GPUCA_NSLICES][GPUTrackingInOutZS::NENDPOINTS];
  };

  struct PtrSave {
    GPUTPCClusterFinder::ZSOffset* zsOffsetHost;
    GPUTPCClusterFinder::ZSOffset* zsOffsetDevice;
    unsigned char* zsDevice;
  };

  std::vector<FragmentData> fragmentData;
  unsigned int nPagesTotal;
  unsigned int nPagesSectorMax;
  unsigned int nPagesFragmentMax;
  unsigned int nPagesSector[GPUCA_NSLICES];
  size_t nMaxDigitsFragment[GPUCA_NSLICES];
  unsigned int tpcMaxTimeBin;
  unsigned int nFragments;
  CfFragment fragmentFirst;
  std::pair<unsigned int, unsigned int> nextPos[GPUCA_NSLICES];
  PtrSave ptrSave[GPUCA_NSLICES];
  const o2::tpc::ClusterNativeAccess* ptrClusterNativeSave;

  void prepare(bool tpcZS, const CfFragment& fragmentMax)
  {
    nPagesTotal = nPagesSectorMax = nPagesFragmentMax = 0;
    for (unsigned int i = 0; i < GPUCA_NSLICES; i++) {
      nPagesSector[i] = 0;
      nMaxDigitsFragment[i] = 0;
    }

    if (tpcZS) {
      tpcMaxTimeBin = 0;
      CfFragment f = fragmentMax;
      while (!f.isEnd()) {
        f = f.next();
      }
      nFragments = f.index;
      if (fragmentData.size() < nFragments) {
        fragmentData.resize(nFragments);
      }

      for (unsigned int i = 0; i < nFragments; i++) {
        for (unsigned int j = 0; j < GPUCA_NSLICES; j++) {
          for (unsigned int k = 0; k < GPUTrackingInOutZS::NENDPOINTS; k++) {
            fragmentData[i].nDigits[j][k] = fragmentData[i].nPages[j][k] = 0;
            fragmentData[i].pageDigits[j][k].clear();
          }
        }
      }
    }
  }
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
