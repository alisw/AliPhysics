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

/// \file GPUO2DataTypes.h
/// \author David Rohr

#ifndef O2_GPU_GPUO2DATATYPES_H
#define O2_GPU_GPUO2DATATYPES_H

// Pull in several O2 headers with basic data types, or load a header with empty fake classes if O2 headers not available

#if defined(HAVE_O2HEADERS) && (!defined(__OPENCL__) || defined(__OPENCLCPP__))
#include "DataFormatsTPC/ClusterNative.h"
#include "DetectorsBase/MatLayerCylSet.h"
#include "TRDBase/TRDGeometryFlat.h"
#else
#include "GPUO2FakeClasses.h"
#endif

#if !defined(__OPENCL__) || defined(__OPENCLCPP__)
#include "GPUdEdxInfo.h"
#endif

#endif
