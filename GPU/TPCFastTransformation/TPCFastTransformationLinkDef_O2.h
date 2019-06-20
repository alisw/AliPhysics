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

/// \file TPCFastTransformationLinkDef_O2.h
/// \author Sergey Gorbunov

#ifdef __CLING__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class o2::gpu::RegularSpline1D+;
#pragma link C++ class o2::gpu::IrregularSpline1D+;
#pragma link C++ class o2::gpu::IrregularSpline2D3D+;
#pragma link C++ class o2::gpu::SemiregularSpline2D3D+;
#pragma link C++ class o2::gpu::IrregularSpline2D3DCalibrator+;
#pragma link C++ class o2::gpu::TPCFastTransformGeo+;
#pragma link C++ class o2::gpu::TPCFastTransform+;
#pragma link C++ class o2::gpu::TPCDistortionIRS+;

#endif
