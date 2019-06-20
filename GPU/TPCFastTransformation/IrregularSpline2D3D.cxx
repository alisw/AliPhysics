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

/// \file  IrregularSpline2D3D.cxx
/// \brief Implementation of IrregularSpline2D3D class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#include "IrregularSpline2D3D.h"

#if !defined(GPUCA_GPUCODE)
#include <iostream>
#endif

using namespace GPUCA_NAMESPACE::gpu;

IrregularSpline2D3D::IrregularSpline2D3D() : FlatObject(), mGridU(), mGridV()
{
  /// Default constructor. Creates an empty uninitialised object
}

void IrregularSpline2D3D::destroy()
{
  /// See FlatObject for description
  mGridU.destroy();
  mGridV.destroy();
  FlatObject::destroy();
}

void IrregularSpline2D3D::relocateBufferPointers(const char* oldBuffer, char* actualBuffer)
{
  /// relocate pointers from old to new buffer location

  char* bufferU = FlatObject::relocatePointer(oldBuffer, actualBuffer, mGridU.getFlatBufferPtr());
  mGridU.setActualBufferAddress(bufferU);

  char* bufferV = FlatObject::relocatePointer(oldBuffer, actualBuffer, mGridV.getFlatBufferPtr());
  mGridV.setActualBufferAddress(bufferV);
}

void IrregularSpline2D3D::cloneFromObject(const IrregularSpline2D3D& obj, char* newFlatBufferPtr)
{
  /// See FlatObject for description

  const char* oldFlatBufferPtr = obj.mFlatBufferPtr;

  FlatObject::cloneFromObject(obj, newFlatBufferPtr);

  char* bufferU = FlatObject::relocatePointer(oldFlatBufferPtr, mFlatBufferPtr, obj.mGridU.getFlatBufferPtr());
  mGridU.cloneFromObject(obj.mGridU, bufferU);

  char* bufferV = FlatObject::relocatePointer(oldFlatBufferPtr, mFlatBufferPtr, obj.mGridV.getFlatBufferPtr());
  mGridV.cloneFromObject(obj.mGridV, bufferV);
}

void IrregularSpline2D3D::moveBufferTo(char* newFlatBufferPtr)
{
  /// See FlatObject for description
  const char* oldFlatBufferPtr = mFlatBufferPtr;
  FlatObject::moveBufferTo(newFlatBufferPtr);
  relocateBufferPointers(oldFlatBufferPtr, mFlatBufferPtr);
}

void IrregularSpline2D3D::setActualBufferAddress(char* actualFlatBufferPtr)
{
  /// See FlatObject for description
  const char* oldFlatBufferPtr = mFlatBufferPtr;
  FlatObject::setActualBufferAddress(actualFlatBufferPtr);
  relocateBufferPointers(oldFlatBufferPtr, mFlatBufferPtr);
}

void IrregularSpline2D3D::setFutureBufferAddress(char* futureFlatBufferPtr)
{
  /// See FlatObject for description
  const char* oldFlatBufferPtr = mFlatBufferPtr;

  char* bufferU = relocatePointer(oldFlatBufferPtr, futureFlatBufferPtr, mGridU.getFlatBufferPtr());
  mGridU.setFutureBufferAddress(bufferU);

  char* bufferV = relocatePointer(oldFlatBufferPtr, futureFlatBufferPtr, mGridV.getFlatBufferPtr());
  mGridV.setFutureBufferAddress(bufferV);

  FlatObject::setFutureBufferAddress(futureFlatBufferPtr);
}

void IrregularSpline2D3D::construct(int numberOfKnotsU, const float knotsU[], int numberOfAxisBinsU, int numberOfKnotsV, const float knotsV[], int numberOfAxisBinsV)
{
  /// Constructor
  ///
  /// Number of knots created and their values may differ from the input values:
  /// - Edge knots 0.f and 1.f will be added if they are not present.
  /// - Knot values are rounded to closest axis bins: k*1./numberOfAxisBins.
  /// - Knots which are too close to each other will be merged
  /// - At least 5 knots and at least 4 axis bins will be created for consistency reason
  ///
  /// \param numberOfKnotsU     U axis: Number of knots in knots[] array
  /// \param knotsU             U axis: Array of knots.
  /// \param numberOfAxisBinsU  U axis: Number of axis bins to map U coordinate to
  ///                           an appropriate [knot(i),knot(i+1)] interval.
  ///                           The knot positions have a "granularity" of 1./numberOfAxisBins
  ///
  /// \param numberOfKnotsV     V axis: Number of knots in knots[] array
  /// \param knotsV             V axis: Array of knots.
  /// \param numberOfAxisBinsV  V axis: Number of axis bins to map U coordinate to
  ///                           an appropriate [knot(i),knot(i+1)] interval.
  ///                           The knot positions have a "granularity" of 1./numberOfAxisBins
  ///

  FlatObject::startConstruction();

  mGridU.construct(numberOfKnotsU, knotsU, numberOfAxisBinsU);
  mGridV.construct(numberOfKnotsV, knotsV, numberOfAxisBinsV);

  size_t vOffset = alignSize(mGridU.getFlatBufferSize(), mGridV.getBufferAlignmentBytes());

  FlatObject::finishConstruction(vOffset + mGridV.getFlatBufferSize());

  mGridU.moveBufferTo(mFlatBufferPtr);
  mGridV.moveBufferTo(mFlatBufferPtr + vOffset);
}

void IrregularSpline2D3D::constructRegular(int numberOfKnotsU, int numberOfKnotsV)
{
  /// Constructor for a regular spline
  /// \param numberOfKnotsU     U axis: Number of knots in knots[] array
  /// \param numberOfKnotsV     V axis: Number of knots in knots[] array
  ///

  FlatObject::startConstruction();

  mGridU.constructRegular(numberOfKnotsU);
  mGridV.constructRegular(numberOfKnotsV);

  size_t vOffset = alignSize(mGridU.getFlatBufferSize(), mGridV.getBufferAlignmentBytes());

  FlatObject::finishConstruction(vOffset + mGridV.getFlatBufferSize());

  mGridU.moveBufferTo(mFlatBufferPtr);
  mGridV.moveBufferTo(mFlatBufferPtr + vOffset);
}

void IrregularSpline2D3D::print() const
{
#if !defined(GPUCA_GPUCODE)
  std::cout << " Irregular Spline 2D3D: " << std::endl;
  std::cout << " grid U: " << std::endl;
  mGridU.print();
  std::cout << " grid V: " << std::endl;
  mGridV.print();
#endif
}
