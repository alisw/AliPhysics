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

/// \file GPUDisplayExt.h
/// \author David Rohr

#ifndef GPUDISPLAYEXT_H
#define GPUDISPLAYEXT_H
#ifdef GPUCA_BUILD_EVENT_DISPLAY

#include "GPUCommonDef.h"

#if defined(GPUCA_DISPLAY_GL3W) && !defined(GPUCA_DISPLAY_OPENGL_CORE)
#define GPUCA_DISPLAY_OPENGL_CORE
#endif

#ifdef GPUCA_DISPLAY_GL3W
#include "GL/gl3w.h"
#else
#include <GL/glew.h>
#endif

namespace GPUCA_NAMESPACE
{
namespace gpu
{
#ifdef GPUCA_DISPLAY_GL3W
static int GPUDisplayExtInit()
{
  return gl3wInit();
}
#else
static int GPUDisplayExtInit()
{
  return glewInit();
}
#endif
#ifdef GPUCA_DISPLAY_OPENGL_CORE
static constexpr bool GPUCA_DISPLAY_OPENGL_CORE_FLAGS = true;
#else
static constexpr bool GPUCA_DISPLAY_OPENGL_CORE_FLAGS = false;
#endif
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
#endif
