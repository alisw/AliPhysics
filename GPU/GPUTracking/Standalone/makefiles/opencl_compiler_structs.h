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

/// \file opencl_compiler_structs.h
/// \author David Rohr

struct _makefiles_opencl_platform_info {
  char platform_profile[64];
  char platform_version[64];
  char platform_name[64];
  char platform_vendor[64];
  cl_uint count;
};

struct _makefiles_opencl_device_info {
  char device_name[64];
  char device_vendor[64];
  cl_uint nbits;
  size_t binary_size;
};
