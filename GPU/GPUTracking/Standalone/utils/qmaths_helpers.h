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

/// \file qmaths_helpers.h
/// \author David Rohr

#ifndef QMATH_HELPERS_H
#define QMATH_HELPERS_H

#if defined __has_include
#if __has_include(<xmmintrin.h>) && __has_include(<pmmintrin.h>)
#include <xmmintrin.h>
#include <pmmintrin.h>
#if defined(_MM_FLUSH_ZERO_OFF) && defined(_MM_DENORMALS_ZERO_ON)
static void disable_denormals()
{
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
}
#define XMM_HAS_DENORMAL_DEACTIVATE
#endif
#endif
#endif
#ifdef XMM_HAS_DENORMAL_DEACTIVATE
#undef XMM_HAS_DENORMAL_DEACTIVATE
#else
static void disable_denormals() {}
#endif

#endif
