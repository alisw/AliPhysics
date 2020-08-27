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

/// \file GPUSettings.h
/// \author David Rohr

#ifndef GPUSETTINGS_H
#define GPUSETTINGS_H

#include "GPUCommonDef.h"
#ifndef __OPENCL__
#include <vector>
#include <string>
#endif

namespace GPUCA_NAMESPACE
{
namespace gpu
{
class GPUDisplayBackend;
class GPUReconstruction;

class GPUSettings
{
 public:
  enum CompressionModes { CompressionTruncate = 1,
                          CompressionDifferences = 2,
                          CompressionTrackModel = 4,
                          CompressionFull = 7 };
  enum CompressionSort { SortTime = 0,
                         SortPad = 1,
                         SortZTimePad = 2,
                         SortZPadTime = 3,
                         SortNoSort = 4 };
  enum CompressionRejection { RejectionNone = 0,
                              RejectionStrategyA = 1,
                              RejectionStrategyB = 2 };

#if !defined(__OPENCL__) || defined(__OPENCLCPP__)
  static CONSTEXPR unsigned int TPC_MAX_TF_TIME_BIN = ((256 * 3564 + 2 * 8 - 2) / 8);
#endif
};

// Settings describing the events / time frames
struct GPUSettingsEvent {
  // All new members must be sizeof(int) resp. sizeof(float) for alignment reasons!
  GPUSettingsEvent();
  float solenoidBz;         // solenoid field strength
  int constBz;              // for test-MC events with constant Bz
  int homemadeEvents;       // Toy-MC events
  int continuousMaxTimeBin; // 0 for triggered events, -1 for default of 23ms
  int needsClusterer;       // Set to true if the data requires the clusterizer
};

// Settings defining the setup of the GPUReconstruction processing (basically selecting the device / class instance)
struct GPUSettingsDeviceBackend {
  GPUSettingsDeviceBackend();
  unsigned int deviceType;   // Device type, shall use GPUDataTypes::DEVICE_TYPE constants, e.g. CPU / CUDA
  char forceDeviceType;      // Fail if device initialization fails, otherwise falls back to CPU
  GPUReconstruction* master; // GPUReconstruction master object
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#ifdef GPUCA_GPUCODE_DEVICE
#define QCONFIG_GPU
#endif
// See GPUSettingsList.h for a list of all available settings of GPU Reconstruction
#include "utils/qconfig.h"

#endif
