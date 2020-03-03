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

/// \file MCLabelAccumulator.h
/// \author Felix Weiglhofer

#ifndef O2_GPU_MC_LABEL_ACCUMULATOR_H
#define O2_GPU_MC_LABEL_ACCUMULATOR_H

#include "clusterFinderDefs.h"
#include "Array2D.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include <bitset>
#include <vector>

namespace GPUCA_NAMESPACE
{
namespace dataformats
{
template <typename T>
class MCTruthContainer;
}

namespace gpu
{

class GPUTPCClusterFinder;
struct GPUTPCClusterMCInterim;

class MCLabelAccumulator
{

 public:
  MCLabelAccumulator(GPUTPCClusterFinder&);

  void collect(const ChargePos&, Charge);

  bool engaged() const { return mLabels != nullptr && mOutput != nullptr; }

  void commit(Row, uint, uint);

 private:
  using MCLabelContainer = o2::dataformats::MCTruthContainer<o2::MCCompLabel>;

  Array2D<const uint> mIndexMap;
  const MCLabelContainer* mLabels = nullptr;
  GPUTPCClusterMCInterim* mOutput = nullptr;

  std::bitset<64> mMaybeHasLabel;
  std::vector<uint64_t> mClusterLabels;
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
