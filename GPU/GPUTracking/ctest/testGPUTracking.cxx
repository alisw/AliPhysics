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

/// \file testGPUTracking.cxx
/// \author David Rohr

#define BOOST_TEST_MODULE Test TPC CA GPU Tracking
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include "GPUO2Interface.h"
#include "TPCFastTransform.h"
#include "GPUO2InterfaceConfiguration.h"
#include "GPUDataTypes.h"

using namespace o2::gpu;

/// @brief Basic test if we can create the interface
BOOST_AUTO_TEST_CASE(CATracking_test1)
{
  auto interface = new GPUTPCO2Interface;

  float solenoidBz = -5.00668; // B-field
  float refX = 1000.;          // transport tracks to this x after tracking, >500 for disabling
  bool continuous = true;      // time frame data v.s. triggered events

  GPUO2InterfaceConfiguration config;
  config.configProcessing.deviceType = GPUDataTypes::DeviceType::CPU;
  config.configProcessing.forceDeviceType = true;

  config.configDeviceProcessing.nThreads = 4;           // 4 threads if we run on the CPU, 1 = default, 0 = auto-detect
  config.configDeviceProcessing.runQA = true;           // Run QA after tracking
  config.configDeviceProcessing.eventDisplay = nullptr; // Ptr to event display backend, for running standalone OpenGL event display
                                                        // config.configDeviceProcessing.eventDisplay = new GPUDisplayBackendX11;

  config.configEvent.solenoidBz = solenoidBz;
  config.configEvent.continuousMaxTimeBin = continuous ? 0.023 * 5e6 : 0; // Number of timebins in timeframe if continuous, 0 otherwise

  config.configReconstruction.NWays = 3;               // Should always be 3!
  config.configReconstruction.NWaysOuter = true;       // Will create outer param for TRD
  config.configReconstruction.SearchWindowDZDR = 2.5f; // Should always be 2.5 for looper-finding and/or continuous tracking
  config.configReconstruction.TrackReferenceX = refX;

  interface->Initialize(config);
  delete interface;
}
