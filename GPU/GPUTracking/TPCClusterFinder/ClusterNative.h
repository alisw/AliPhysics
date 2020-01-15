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
#if !defined(CLUSTER_NATIVE_H)
#define CLUSTER_NATIVE_H

#define CN_SCALE_TIME_PACKED 64
#define CN_SCALE_PAD_PACKED 64
#define CN_SCALE_SIGMA_TIME_PACKED 32
#define CN_SCALE_SIGMA_PAD_PACKED 32

#define CN_TIME_MASK 0xFFFFFF
#define CN_FLAG_MASK 0xFF000000

#if defined(__OPENCL__) && !defined(__OPENCL_CPP__)
#define CAST(type, val) ((type)(val))
#else
#define CAST(type, val) static_cast<type>(val)
#endif

namespace GPUCA_NAMESPACE
{
namespace gpu
{
namespace deprecated
{

struct ClusterNative {
  uint timeFlagsPacked;
  ushort padPacked;
  unsigned char sigmaTimePacked;
  unsigned char sigmaPadPacked;
  ushort qmax;
  ushort qtot;
};

enum CnFlagPos {
  CN_FLAG_POS_IS_EDGE_CLUSTER = 0,
  CN_FLAG_POS_SPLIT_IN_TIME = 1,
  CN_FLAG_POS_SPLIT_IN_PAD = 2,
};

enum CnFlag {
  CN_FLAG_IS_EDGE_CLUSTER = (1 << CN_FLAG_POS_IS_EDGE_CLUSTER),
  CN_FLAG_SPLIT_IN_TIME = (1 << CN_FLAG_POS_SPLIT_IN_TIME),
  CN_FLAG_SPLIT_IN_PAD = (1 << CN_FLAG_POS_SPLIT_IN_PAD),
};

GPUdi() ushort cnPackPad(float pad)
{
  return CAST(ushort, pad * CN_SCALE_PAD_PACKED + 0.5f);
}

GPUdi() uint cnPackTime(float time)
{
  return CAST(uint, time * CN_SCALE_TIME_PACKED + 0.5f);
}

GPUdi() float cnUnpackPad(ushort pad)
{
  return CAST(float, pad) * (1.f / CN_SCALE_PAD_PACKED);
}

GPUdi() float cnUnpackTime(uint time)
{
  return CAST(float, time) * (1.f / CN_SCALE_TIME_PACKED);
}

GPUdi() unsigned char cnPackSigma(float sigma, float scale)
{
  uint tmp = sigma * scale + 0.5f;
  return (tmp > 0xFF) ? 0xFF : tmp;
}

GPUdi() float cnUnpackSigma(
  unsigned char sigmaPacked,
  float scale)
{
  return CAST(float, sigmaPacked) * (1.f / scale);
}

GPUdi() unsigned char cnGetFlags(const ClusterNative* c)
{
  return c->timeFlagsPacked >> 24;
}

GPUdi() uint cnGetTimePacked(const ClusterNative* c)
{
  return c->timeFlagsPacked & CN_TIME_MASK;
}

GPUdi() void cnSetTimePackedFlags(
  ClusterNative* c,
  uint timePacked,
  unsigned char flags)
{
  c->timeFlagsPacked = (timePacked & CN_TIME_MASK) | CAST(uint, flags) << 24;
}

GPUdi() void cnSetTimePacked(ClusterNative* c, uint timePacked)
{
  c->timeFlagsPacked = (timePacked & CN_TIME_MASK) | (c->timeFlagsPacked & CN_FLAG_MASK);
}

GPUdi() void cnSetFlags(ClusterNative* c, unsigned char flags)
{
  c->timeFlagsPacked = (c->timeFlagsPacked & CN_TIME_MASK) | (CAST(uint, flags) << 24);
}

GPUdi() float cnGetTime(const ClusterNative* c)
{
  return cnUnpackTime(c->timeFlagsPacked & CN_TIME_MASK);
}

GPUdi() void cnSetTime(ClusterNative* c, float time)
{
  c->timeFlagsPacked = (cnPackTime(time) & CN_TIME_MASK) | (c->timeFlagsPacked & CN_FLAG_MASK);
}

GPUdi() void cnSetTimeFlags(
  ClusterNative* c,
  float time,
  unsigned char flags)
{
  c->timeFlagsPacked = (cnPackTime(time) & CN_TIME_MASK) | (CAST(uint, flags) << 24);
}

GPUdi() float cnGetPad(const ClusterNative* c)
{
  return cnUnpackPad(c->padPacked);
}

GPUdi() void cnSetPad(ClusterNative* c, float pad)
{
  c->padPacked = cnPackPad(pad);
}

GPUdi() float cnGetSigmaTime(const ClusterNative* c)
{
  return cnUnpackSigma(c->sigmaTimePacked, CN_SCALE_SIGMA_TIME_PACKED);
}

GPUdi() void cnSetSigmaTime(ClusterNative* c, float sigmaTime)
{
  c->sigmaTimePacked = cnPackSigma(sigmaTime, CN_SCALE_SIGMA_TIME_PACKED);
}

GPUdi() float cnGetSigmaPad(const ClusterNative* c)
{
  return cnUnpackSigma(c->sigmaPadPacked, CN_SCALE_SIGMA_PAD_PACKED);
}

GPUdi() void cnSetSigmaPad(ClusterNative* c, float sigmaPad)
{
  c->sigmaPadPacked = cnPackSigma(sigmaPad, CN_SCALE_SIGMA_PAD_PACKED);
}

} // namespace deprecated
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif // !defined(CLUSTER_NATIVE_H)
// vim: set ts=4 sw=4 sts=4 expandtab:
