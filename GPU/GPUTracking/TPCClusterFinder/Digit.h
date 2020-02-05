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

/// \file digit.h
/// \author Felix Weiglhofer
//
#if !defined(SHARED_DIGIT_H)
#define SHARED_DIGIT_H

#include "clusterFinderDefs.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{
namespace deprecated
{

struct PackedDigit {
  float charge;
  Timestamp time;
  Pad pad;
  Row row;
};

using Digit = PackedDigit;

} // namespace deprecated
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif //!defined(SHARED_DIGIT_H)

// vim: set ts=4 sw=4 sts=4 expandtab:
