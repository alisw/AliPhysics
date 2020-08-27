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

/// \file GPUO2InterfaceConfigurableParam.cxx
/// \author David Rohr

#include "GPUO2InterfaceConfigurableParam.h"
#include "GPUO2InterfaceConfiguration.h"
#include "GPUDataTypes.h"

using namespace o2::gpu;
#define BeginNamespace(name)
#define EndNamespace()
#define AddOption(name, type, default, optname, optnameshort, help, ...)
#define AddVariable(name, type, default)
#define AddOptionSet(name, type, value, optname, optnameshort, help, ...)
#define AddOptionVec(name, type, optname, optnameshort, help, ...)
#define AddOptionArray(name, type, count, default, optname, optnameshort, help, ...)
#define AddSubConfig(name, instance)
#define BeginSubConfig(name, instance, parent, preoptname, preoptnameshort, descr) O2ParamImpl(GPUCA_M_CAT(GPUConfigurableParam, name))
#define EndConfig()
#define AddCustomCPP(...)
#define AddHelp(...)
#define AddShortcut(...)
#include "GPUSettingsList.h"
#undef BeginNamespace
#undef EndNamespace
#undef AddOption
#undef AddVariable
#undef AddOptionSet
#undef AddOptionVec
#undef AddOptionArray
#undef AddSubConfig
#undef BeginSubConfig
#undef EndConfig
#undef AddCustomCPP
#undef AddHelp
#undef AddShortcut

GPUSettingsO2 GPUO2InterfaceConfiguration::ReadConfigurableParam()
{
#define BeginNamespace(name)
#define EndNamespace()
#define AddOption(name, type, default, optname, optnameshort, help, ...) dst.name = src.name;
#define AddVariable(name, type, default)
#define AddOptionSet(name, type, value, optname, optnameshort, help, ...)
#define AddOptionVec(name, type, optname, optnameshort, help, ...)
#define AddOptionArray(name, type, count, default, optname, optnameshort, help, ...) \
  for (int i = 0; i < count; i++) {                                                  \
    dst.name[i] = src.name[i];                                                       \
  }
#define AddSubConfig(name, instance)
#define BeginSubConfig(name, instance, parent, preoptname, preoptnameshort, descr) \
  name instance;                                                                   \
  {                                                                                \
    auto& src = GPUCA_M_CAT(GPUConfigurableParam, name)::Instance();               \
    name& dst = instance;
#define EndConfig() }
#define AddCustomCPP(...)
#define AddHelp(...)
#define AddShortcut(...)
#include "GPUSettingsList.h"
#undef BeginNamespace
#undef EndNamespace
#undef AddOption
#undef AddVariable
#undef AddOptionSet
#undef AddOptionVec
#undef AddOptionArray
#undef AddSubConfig
#undef BeginSubConfig
#undef EndConfig
#undef AddCustomCPP
#undef AddHelp
#undef AddShortcut

  configProcessing = proc;
  configReconstruction = rec;
  configDisplay = GL;
  configQA = QA;
  if (global.continuousMaxTimeBin) {
    configEvent.continuousMaxTimeBin = global.continuousMaxTimeBin;
  }
  if (global.solenoidBz > -1000.f) {
    configEvent.solenoidBz = global.solenoidBz;
  }
  if (global.constBz) {
    configEvent.constBz = global.constBz;
  }
  if (configReconstruction.TrackReferenceX == 1000.f) {
    configReconstruction.TrackReferenceX = 83.f;
  }
  configDeviceBackend.deviceType = GPUDataTypes::GetDeviceType(global.deviceType.c_str());
  configDeviceBackend.forceDeviceType = global.forceDeviceType;
  return global;
}
