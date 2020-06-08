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

/// \file PackedCharge.h
/// \author Felix Weiglhofer

#ifndef O2_GPU_ARRAY2D_H
#define O2_GPU_ARRAY2D_H

#include "clusterFinderDefs.h"
#include "ChargePos.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{

template <typename T, typename Layout>
class AbstractArray2D
{

 public:
  GPUdi() explicit AbstractArray2D(T* d) : data(d) {}

  GPUdi() T& operator[](const ChargePos& p) { return data[Layout::idx(p)]; }
  GPUdi() const T& operator[](const ChargePos& p) const { return data[Layout::idx(p)]; }

  GPUdi() void safeWrite(const ChargePos& p, const T& v)
  {
    if (data != nullptr) {
      (*this)[p] = v;
    }
  }

 private:
  T* data;
};

template <typename Grid>
class TilingLayout
{
 public:
  enum {
    Height = Grid::Height,
    Width = Grid::Width,
    WidthInTiles = (TPC_NUM_OF_PADS + Width - 1) / Width,
  };

  GPUdi() static tpccf::SizeT idx(const ChargePos& p)
  {
    const tpccf::SizeT tilePad = p.gpad / Width;
    const tpccf::SizeT tileTime = p.timePadded / Height;

    const tpccf::SizeT inTilePad = p.gpad % Width;
    const tpccf::SizeT inTileTime = p.timePadded % Height;

    return (tileTime * WidthInTiles + tilePad) * (Width * Height) + inTileTime * Width + inTilePad;
  }

  GPUd() static size_t items()
  {
    return (TPC_NUM_OF_PADS + Width - 1) / Width * Width * (TPC_MAX_FRAGMENT_LEN_PADDED + Height - 1) / Height * Height;
  }
};

class LinearLayout
{
 public:
  GPUdi() static tpccf::SizeT idx(const ChargePos& p)
  {
    return TPC_NUM_OF_PADS * p.timePadded + p.gpad;
  }

  GPUd() static size_t items()
  {
    return TPC_NUM_OF_PADS * TPC_MAX_FRAGMENT_LEN_PADDED;
  }
};

template <tpccf::SizeT S>
struct GridSize;

template <>
struct GridSize<1> {
  enum {
    Width = 8,
    Height = 8,
  };
};

template <>
struct GridSize<2> {
  enum {
    Width = 4,
    Height = 4,
  };
};

template <>
struct GridSize<4> {
  enum {
    Width = 4,
    Height = 4,
  };
};

#if defined(CHARGEMAP_TILING_LAYOUT)
template <typename T>
using TPCMapMemoryLayout = TilingLayout<GridSize<sizeof(T)>>;
#else
template <typename T>
using TPCMapMemoryLayout = LinearLayout;
#endif

template <typename T>
using Array2D = AbstractArray2D<T, TPCMapMemoryLayout<T>>;

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
