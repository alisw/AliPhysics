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

/// \file GPUSettings.cxx
/// \author David Rohr

#include "GPUSettings.h"
#include "GPUDef.h"
#include "GPUDataTypes.h"
#include <cstring>

using namespace GPUCA_NAMESPACE::gpu;

GPUSettingsDeviceBackend::GPUSettingsDeviceBackend()
{
  deviceType = GPUDataTypes::DeviceType::CPU;
  forceDeviceType = true;
  master = nullptr;
}

GPUSettingsEvent::GPUSettingsEvent()
{
  solenoidBz = -5.00668;
  constBz = 0;
  homemadeEvents = 0;
  continuousMaxTimeBin = 0;
  needsClusterer = 0;
}
