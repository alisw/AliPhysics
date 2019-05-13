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

/// \file GPUTPCMCPoint.cxx
/// \author Sergey Gorbunov, Ivan Kisel, David Rohr

#include "GPUTPCMCPoint.h"

GPUTPCMCPoint::GPUTPCMCPoint() : fX(0), fY(0), fZ(0), fSx(0), fSy(0), fSz(0), fTime(0), mISlice(0), fTrackID(0)
{
  //* Default constructor
}
