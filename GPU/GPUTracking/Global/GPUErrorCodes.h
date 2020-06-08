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

/// \file GPUErrorCodes.h
/// \author David Rohr

// Error Codes for GPU Tracker
GPUCA_ERROR_CODE(0, ERROR_NONE)
GPUCA_ERROR_CODE(1, ERROR_ROWSTARTHIT_OVERFLOW)
GPUCA_ERROR_CODE(2, ERROR_STARTHIT_OVERFLOW)
GPUCA_ERROR_CODE(3, ERROR_TRACKLET_OVERFLOW)
GPUCA_ERROR_CODE(4, ERROR_TRACKLET_HIT_OVERFLOW)
GPUCA_ERROR_CODE(5, ERROR_TRACK_OVERFLOW)
GPUCA_ERROR_CODE(6, ERROR_TRACK_HIT_OVERFLOW)
GPUCA_ERROR_CODE(7, ERROR_GLOBAL_TRACKING_TRACK_OVERFLOW)
GPUCA_ERROR_CODE(8, ERROR_GLOBAL_TRACKING_TRACK_HIT_OVERFLOW)
GPUCA_ERROR_CODE(9, ERROR_LOOPER_OVERFLOW)
GPUCA_ERROR_CODE(10, ERROR_MERGER_CE_HIT_OVERFLOW)
GPUCA_ERROR_CODE(11, ERROR_MERGER_LOOPER_OVERFLOW)
GPUCA_ERROR_CODE(12, ERROR_SLICEDATA_FIRSTHITINBIN_OVERFLOW)
GPUCA_ERROR_CODE(13, ERROR_SLICEDATA_HITINROW_OVERFLOW)
GPUCA_ERROR_CODE(14, ERROR_SLICEDATA_BIN_OVERFLOW)
GPUCA_ERROR_CODE(15, ERROR_SLICEDATA_Z_OVERFLOW)
