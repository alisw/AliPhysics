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

/// \file GPUO2InterfaceConfiguration.h
/// \author David Rohr

#ifndef GPUO2INTERFACECONFIGURATION_H
#define GPUO2INTERFACECONFIGURATION_H

#ifndef GPUCA_O2_LIB
#define GPUCA_O2_LIB
#endif
#ifndef HAVE_O2HEADERS
#define HAVE_O2HEADERS
#endif
#ifndef GPUCA_TPC_GEOMETRY_O2
#define GPUCA_TPC_GEOMETRY_O2
#endif

#include <memory>
#include "GPUSettings.h"
#include "GPUDisplayConfig.h"
#include "GPUQAConfig.h"
#include "TPCFastTransform.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{
// Full configuration structure with all available settings of GPU...
struct GPUO2InterfaceConfiguration {
  GPUO2InterfaceConfiguration() = default;
  ~GPUO2InterfaceConfiguration() = default;
  GPUO2InterfaceConfiguration(const GPUO2InterfaceConfiguration&) = default;

  // Settings for the Interface class
  struct GPUInterfaceSettings {
    bool dumpEvents = false;
  };

  GPUSettingsProcessing configProcessing;
  GPUSettingsDeviceProcessing configDeviceProcessing;
  GPUSettingsEvent configEvent;
  GPUSettingsRec configReconstruction;
  GPUDisplayConfig configDisplay;
  GPUQAConfig configQA;
  GPUInterfaceSettings configInterface;
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
