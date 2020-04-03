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

/// \file GPUHostDataTypes.h
/// \author David Rohr

#ifndef GPUHOSTDATATYPES_H
#define GPUHOSTDATATYPES_H

#include "GPUCommonDef.h"

// These are complex data types wrapped in simple structs, which can be forward declared.
// Structures used on the GPU can have pointers to these wrappers, when the wrappers are forward declared.
// These wrapped complex types are not meant for usage on GPU

#if defined(GPUCA_GPUCODE)
#error "GPUHostDataTypes.h should never be included on GPU."
#endif

#include <vector>
#include <array>
#include <memory>
#include <mutex>
#include "DataFormatsTPC/Constants.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{

struct GPUTPCDigitsMCInput {
  std::array<const o2::dataformats::MCTruthContainer<o2::MCCompLabel>*, o2::tpc::Constants::MAXSECTOR> v;
};

struct GPUTPCClusterMCInterim {
  std::vector<o2::MCCompLabel> labels;
  uint offset;
};

struct GPUTPCLinearLabels {
  std::vector<o2::dataformats::MCTruthHeaderElement> header;
  std::vector<o2::MCCompLabel> data;
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
