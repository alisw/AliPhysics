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

namespace GPUCA_NAMESPACE
{
namespace gpu
{

class Array2D
{

 public:
  // Maps the position of a pad given as row and index in that row to a unique
  // index between 0 and TPC_NUM_OF_PADS.
  static GPUdi() GlobalPad tpcGlobalPadIdx(Row row, Pad pad)
  {
    return TPC_PADS_PER_ROW_PADDED * row + pad + PADDING_PAD;
  }

  static GPUdi() size_t idxSquareTiling(GlobalPad gpad, Timestamp time, size_t N)
  {
    /* time += PADDING; */

    /* const size_t C = TPC_NUM_OF_PADS + N - 1; */

    /* const size_t inTilePad = gpad % N; */
    /* const size_t inTileTime = time % N; */

    /* return N * (time * C + gpad + inTileTime) + inTilePad; */

    time += PADDING_TIME;

    const size_t tileW = N;
    const size_t tileH = N;
    const size_t widthInTiles = (TPC_NUM_OF_PADS + tileW - 1) / tileW;

    const size_t tilePad = gpad / tileW;
    const size_t tileTime = time / tileH;

    const size_t inTilePad = gpad % tileW;
    const size_t inTileTime = time % tileH;

    return (tileTime * widthInTiles + tilePad) * (tileW * tileH) + inTileTime * tileW + inTilePad;
  }

  static GPUdi() size_t idxTiling4x4(GlobalPad gpad, Timestamp time)
  {
    return idxSquareTiling(gpad, time, 4);
  }

  static GPUdi() size_t idxTiling8x8(GlobalPad gpad, Timestamp time)
  {
    return idxSquareTiling(gpad, time, 8);
  }

  static GPUdi() size_t chargemapIdx(GlobalPad gpad, Timestamp time)
  {

#if defined(CHARGEMAP_4x4_TILING_LAYOUT) || defined(CHARGEMAP_4x8_TILING_LAYOUT) || defined(CHARGEMAP_8x4_TILING_LAYOUT)
#define CHARGEMAP_TILING_LAYOUT
#endif

#if defined(CHARGEMAP_4x4_TILING_LAYOUT)

    return idxTiling4x4(gpad, time);

#elif defined(CHARGEMAP_4x8_TILING_LAYOUT)
#define TILE_WIDTH 4
#define TILE_HEIGHT 8

    time += PADDING_TIME;

    const size_t tileW = TILE_WIDTH;
    const size_t tileH = TILE_HEIGHT;
    const size_t widthInTiles = (TPC_NUM_OF_PADS + tileW - 1) / tileW;

    const size_t tilePad = gpad / tileW;
    const size_t tileTime = time / tileH;

    const size_t inTilePad = gpad % tileW;
    const size_t inTileTime = time % tileH;

    return (tileTime * widthInTiles + tilePad) * (tileW * tileH) + inTileTime * tileW + inTilePad;

#undef TILE_WIDTH
#undef TILE_HEIGHT

#elif defined(CHARGEMAP_8x4_TILING_LAYOUT)
#define TILE_WIDTH 8
#define TILE_HEIGHT 4

    time += PADDING_TIME;

    const size_t tileW = TILE_WIDTH;
    const size_t tileH = TILE_HEIGHT;
    const size_t widthInTiles = (TPC_NUM_OF_PADS + tileW - 1) / tileW;

    const size_t tilePad = gpad / tileW;
    const size_t tileTime = time / tileH;

    const size_t inTilePad = gpad % tileW;
    const size_t inTileTime = time % tileH;

    return (tileTime * widthInTiles + tilePad) * (tileW * tileH) + inTileTime * tileW + inTilePad;

#undef TILE_WIDTH
#undef TILE_HEIGHT

#elif defined(CHARGEMAP_PAD_MAJOR_LAYOUT)

    time += PADDING_TIME;
    return TPC_MAX_TIME_PADDED * gpad + time;

#else // Use row time-major layout

    time += PADDING_TIME;
    return TPC_NUM_OF_PADS * time + gpad;

#endif
  }
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#define ACCESS_2D(map, idxFunc, gpad, time) map[idxFunc(gpad, time)]

#define CHARGE(map, gpad, time) ACCESS_2D(map, Array2D::chargemapIdx, gpad, time)

#if defined(CHARGEMAP_TILING_LAYOUT)
/* #if 0 */
#define IS_PEAK(map, gpad, time) ACCESS_2D(map, Array2D::idxTiling8x8, gpad, time)
#else
#define IS_PEAK(map, gpad, time) ACCESS_2D(map, Array2D::chargemapIdx, gpad, time)
#endif

#endif
