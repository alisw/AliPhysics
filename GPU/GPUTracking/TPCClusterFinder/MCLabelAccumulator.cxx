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

/// \file MCLabelAccumulator.cxx
/// \author Felix Weiglhofer

#include "MCLabelAccumulator.h"

#include "GPUHostDataTypes.h"
#include "GPUTPCClusterFinder.h"
#include "SimulationDataFormat/MCTruthContainer.h"

using namespace GPUCA_NAMESPACE::gpu;

MCLabelAccumulator::MCLabelAccumulator(GPUTPCClusterFinder& clusterer)
  : mIndexMap(clusterer.mPindexMap), mLabels(clusterer.mPinputLabels), mOutput(clusterer.mPlabelsByRow)
{
  mClusterLabels.reserve(32);
}

void MCLabelAccumulator::collect(const ChargePos& pos, Charge q)
{
  if (q == 0 || !engaged()) {
    return;
  }

  uint index = mIndexMap[pos];

  auto labels = mLabels->getLabels(index);

  for (const auto& label : labels) {
    int h = label.getRawValue() % mMaybeHasLabel.size();

    if (mMaybeHasLabel[h]) {
      auto lookup = std::find(mClusterLabels.begin(), mClusterLabels.end(), label);
      if (lookup != mClusterLabels.end()) {
        continue;
      }
    }

    mMaybeHasLabel[h] = true;
    mClusterLabels.emplace_back(label);
  }
}

void MCLabelAccumulator::commit(Row row, uint indexInRow, uint maxElemsPerBucket)
{
  if (indexInRow > maxElemsPerBucket || !engaged()) {
    return;
  }

  auto& out = mOutput[row * maxElemsPerBucket + indexInRow];
  out.labels = std::move(mClusterLabels);
}
