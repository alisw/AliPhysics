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

/// \file GPUO2InterfaceConfigurableParam.h
/// \author David Rohr

// This file auto-generates a ConfigurableParam object from the GPU parameter macros.
// Set via:
// --configKeyValues "GPU_global.[x]=[y]" : for global GPU run configurations, like solenoidBz, gpuType, configuration object files.
// --configKeyValues "GPU_rec.[x]=[y]" : for GPU reconstruction related settings used on the GPU, like pt threshold for track rejection.
// --configKeyValues "GPU_proc.[x]=[y]" : for processing options steering GPU reconstruction like GPU device ID, debug output level, number of CPU threads.
// Check GPUSettingsList.h for all options

#ifndef GPUO2INTERFACECONFIGURABLEPARAM_H
#define GPUO2INTERFACECONFIGURABLEPARAM_H

// Some defines denoting that we are compiling for O2
#ifndef HAVE_O2HEADERS
#define HAVE_O2HEADERS
#endif
#ifndef GPUCA_TPC_GEOMETRY_O2
#define GPUCA_TPC_GEOMETRY_O2
#endif
#ifndef GPUCA_O2_INTERFACE
#define GPUCA_O2_INTERFACE
#endif

#include "CommonUtils/ConfigurableParam.h"
#include "CommonUtils/ConfigurableParamHelper.h"
#include "GPUSettings.h"
#include "GPUDefMacros.h"
#include <vector>

#define BeginNamespace(name) \
  namespace name             \
  {
#define EndNamespace() }
#define AddOption(name, type, default, optname, optnameshort, help, ...) type name = default;
#define AddVariable(name, type, default)
#define AddOptionSet(name, type, value, optname, optnameshort, help, ...)
#define AddOptionVec(name, type, optname, optnameshort, help, ...)
#define AddOptionArray(name, type, count, default, optname, optnameshort, help, ...) type name[count] = {default};
#define AddSubConfig(name, instance)
#define BeginSubConfig(name, instance, parent, preoptname, preoptnameshort, descr)                                                     \
  struct GPUCA_M_CAT(GPUConfigurableParam, name) : public o2::conf::ConfigurableParamHelper<GPUCA_M_CAT(GPUConfigurableParam, name)> { \
    O2ParamDef(GPUCA_M_CAT(GPUConfigurableParam, name), GPUCA_M_STR(GPUCA_M_CAT(GPU_, instance))) public:
#define EndConfig() \
  }                 \
  ;
#define AddCustomCPP(...) __VA_ARGS__
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

#endif
