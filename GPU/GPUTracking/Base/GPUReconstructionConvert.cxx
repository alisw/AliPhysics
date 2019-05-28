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

/// \file GPUReconstructionConvert.cxx
/// \author David Rohr

#include "GPUReconstructionConvert.h"
#include "TPCFastTransform.h"
#include "GPUTPCClusterData.h"
#include "ClusterNativeAccessExt.h"
#include "AliHLTTPCRawCluster.h"

using namespace GPUCA_NAMESPACE::gpu;
using namespace o2::tpc;

void GPUReconstructionConvert::ConvertNativeToClusterData(ClusterNativeAccessExt* native, std::unique_ptr<GPUTPCClusterData[]>* clusters, unsigned int* nClusters, const TPCFastTransform* transform, int continuousMaxTimeBin)
{
#ifdef HAVE_O2HEADERS
  memset(nClusters, 0, NSLICES * sizeof(nClusters[0]));
  unsigned int offset = 0;
  for (unsigned int i = 0; i < NSLICES; i++) {
    unsigned int nClSlice = 0;
    for (int j = 0; j < Constants::MAXGLOBALPADROW; j++) {
      nClSlice += native->nClusters[i][j];
    }
    nClusters[i] = nClSlice;
    clusters[i].reset(new GPUTPCClusterData[nClSlice]);
    nClSlice = 0;
    for (int j = 0; j < Constants::MAXGLOBALPADROW; j++) {
      for (unsigned int k = 0; k < native->nClusters[i][j]; k++) {
        const auto& cin = native->clusters[i][j][k];
        float x = 0, y = 0, z = 0;
        if (continuousMaxTimeBin == 0) {
          transform->Transform(i, j, cin.getPad(), cin.getTime(), x, y, z);
        } else {
          transform->TransformInTimeFrame(i, j, cin.getPad(), cin.getTime(), x, y, z, continuousMaxTimeBin);
        }
        auto& cout = clusters[i].get()[nClSlice];
        cout.x = x;
        cout.y = y;
        cout.z = z;
        cout.row = j;
        cout.amp = cin.qMax;
        cout.flags = cin.getFlags();
        cout.id = offset + k;
        nClSlice++;
      }
      native->clusterOffset[i][j] = offset;
      offset += native->nClusters[i][j];
    }
  }
#endif
}

void GPUReconstructionConvert::ConvertRun2RawToNative(ClusterNativeAccessExt* native, std::unique_ptr<ClusterNative[]>* nativeBuffers, const AliHLTTPCRawCluster** rawClusters, unsigned int* nRawClusters)
{
#ifdef HAVE_O2HEADERS
  memset((void*)native, 0, sizeof(*native));
  for (int i = 0; i < Constants::MAXSECTOR; i++) {
    for (unsigned int j = 0; j < nRawClusters[i]; j++) {
      native->nClusters[i][rawClusters[i][j].GetPadRow()]++;
    }
    for (int j = 0; j < Constants::MAXGLOBALPADROW; j++) {
      nativeBuffers[i * Constants::MAXGLOBALPADROW + j].reset(new ClusterNative[native->nClusters[i][j]]);
      native->clusters[i][j] = nativeBuffers[i * Constants::MAXGLOBALPADROW + j].get();
      native->nClusters[i][j] = 0;
    }
    for (unsigned int j = 0; j < nRawClusters[i]; j++) {
      const AliHLTTPCRawCluster& org = rawClusters[i][j];
      int row = org.GetPadRow();
      ClusterNative& c = nativeBuffers[i * Constants::MAXGLOBALPADROW + row][native->nClusters[i][row]++];
      c.setTimeFlags(org.GetTime(), org.GetFlags());
      c.setPad(org.GetPad());
      c.setSigmaTime(std::sqrt(org.GetSigmaTime2()));
      c.setSigmaPad(std::sqrt(org.GetSigmaPad2()));
      c.qMax = org.GetQMax();
      c.qTot = org.GetCharge();
    }
  }
#endif
}
