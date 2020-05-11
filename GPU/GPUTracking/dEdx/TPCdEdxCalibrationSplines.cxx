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

/// \file  dEdxCalibrationSplines.cxx
/// \brief Definition of dEdxCalibrationSplines class
///
/// \author  Matthias Kleiner <matthias.kleiner@cern.ch>

#include "TPCdEdxCalibrationSplines.h"

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE) // code invisible on GPU and in the standalone compilation
#include "TFile.h"
#endif

using namespace GPUCA_NAMESPACE::gpu;

#if !defined(GPUCA_GPUCODE) && defined(GPUCA_STANDALONE)
TPCdEdxCalibrationSplines::TPCdEdxCalibrationSplines()
  : FlatObject()
{
  /// Empty constructor
}
#elif !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE)
TPCdEdxCalibrationSplines::TPCdEdxCalibrationSplines()
  : FlatObject()
{
  /// Default constructor
  setDefaultSplines();
}

TPCdEdxCalibrationSplines::TPCdEdxCalibrationSplines(const TPCdEdxCalibrationSplines& obj)
  : FlatObject()
{
  /// Copy constructor
  this->cloneFromObject(obj, nullptr);
}

TPCdEdxCalibrationSplines& TPCdEdxCalibrationSplines::operator=(const TPCdEdxCalibrationSplines& obj)
{
  /// Assignment operator
  this->cloneFromObject(obj, nullptr);
  return *this;
}

void TPCdEdxCalibrationSplines::recreate(int nKnotsU1[], int nKnotsU2[])
{
  /// Default constructor

  FlatObject::startConstruction();

  int buffSize = 0;
  int offsets1[FSplines];
  int offsets2[FSplines];
  for (unsigned int i = 0; i < FSplines; i++) {
    mCalibSplinesqMax[i].recreate(nKnotsU1[i], nKnotsU2[i]);
    buffSize = alignSize(buffSize, mCalibSplinesqMax[i].getBufferAlignmentBytes());
    offsets1[i] = buffSize;
    buffSize += mCalibSplinesqMax[i].getFlatBufferSize();
  }
  for (unsigned int i = 0; i < FSplines; i++) {
    mCalibSplinesqTot[i].recreate(nKnotsU1[i], nKnotsU2[i]);
    buffSize = alignSize(buffSize, mCalibSplinesqTot[i].getBufferAlignmentBytes());
    offsets2[i] = buffSize;
    buffSize += mCalibSplinesqTot[i].getFlatBufferSize();
  }

  FlatObject::finishConstruction(buffSize);

  for (unsigned int i = 0; i < FSplines; i++) {
    mCalibSplinesqMax[i].moveBufferTo(mFlatBufferPtr + offsets1[i]);
  }
  for (unsigned int i = 0; i < FSplines; i++) {
    mCalibSplinesqTot[i].moveBufferTo(mFlatBufferPtr + offsets2[i]);
  }
}

void TPCdEdxCalibrationSplines::cloneFromObject(const TPCdEdxCalibrationSplines& obj, char* newFlatBufferPtr)
{
  /// See FlatObject for description

  const char* oldFlatBufferPtr = obj.mFlatBufferPtr;
  FlatObject::cloneFromObject(obj, newFlatBufferPtr);

  for (unsigned int i = 0; i < FSplines; i++) {
    char* buffer = FlatObject::relocatePointer(oldFlatBufferPtr, mFlatBufferPtr, obj.mCalibSplinesqMax[i].getFlatBufferPtr());
    mCalibSplinesqMax[i].cloneFromObject(obj.mCalibSplinesqMax[i], buffer);
  }

  for (unsigned int i = 0; i < FSplines; i++) {
    char* buffer = FlatObject::relocatePointer(oldFlatBufferPtr, mFlatBufferPtr, obj.mCalibSplinesqTot[i].getFlatBufferPtr());
    mCalibSplinesqTot[i].cloneFromObject(obj.mCalibSplinesqTot[i], buffer);
  }
}

void TPCdEdxCalibrationSplines::moveBufferTo(char* newFlatBufferPtr)
{
  /// See FlatObject for description
  char* oldFlatBufferPtr = mFlatBufferPtr;
  FlatObject::moveBufferTo(newFlatBufferPtr);
  char* currFlatBufferPtr = mFlatBufferPtr;
  mFlatBufferPtr = oldFlatBufferPtr;
  setActualBufferAddress(currFlatBufferPtr);
}
#endif

void TPCdEdxCalibrationSplines::destroy()
{
  /// See FlatObject for description
  for (unsigned int i = 0; i < FSplines; i++) {
    mCalibSplinesqMax[i].destroy();
    mCalibSplinesqTot[i].destroy();
  }
  FlatObject::destroy();
}

void TPCdEdxCalibrationSplines::setActualBufferAddress(char* actualFlatBufferPtr)
{
  /// See FlatObject for description

  FlatObject::setActualBufferAddress(actualFlatBufferPtr);
  int offset = 0;
  for (unsigned int i = 0; i < FSplines; i++) {
    offset = alignSize(offset, mCalibSplinesqMax[i].getBufferAlignmentBytes());
    mCalibSplinesqMax[i].setActualBufferAddress(mFlatBufferPtr + offset);
    offset += mCalibSplinesqMax[i].getFlatBufferSize();
  }
  for (unsigned int i = 0; i < FSplines; i++) {
    offset = alignSize(offset, mCalibSplinesqTot[i].getBufferAlignmentBytes());
    mCalibSplinesqTot[i].setActualBufferAddress(mFlatBufferPtr + offset);
    offset += mCalibSplinesqTot[i].getFlatBufferSize();
  }
}

void TPCdEdxCalibrationSplines::setFutureBufferAddress(char* futureFlatBufferPtr)
{
  /// See FlatObject for description

  for (unsigned int i = 0; i < FSplines; i++) {
    char* buffer = relocatePointer(mFlatBufferPtr, futureFlatBufferPtr, mCalibSplinesqMax[i].getFlatBufferPtr());
    mCalibSplinesqMax[i].setFutureBufferAddress(buffer);
  }
  for (unsigned int i = 0; i < FSplines; i++) {
    char* buffer = relocatePointer(mFlatBufferPtr, futureFlatBufferPtr, mCalibSplinesqTot[i].getFlatBufferPtr());
    mCalibSplinesqTot[i].setFutureBufferAddress(buffer);
  }
  FlatObject::setFutureBufferAddress(futureFlatBufferPtr);
}

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE)

TPCdEdxCalibrationSplines* TPCdEdxCalibrationSplines::readFromFile(
  TFile& inpf, const char* name)
{
  /// read a class object from the file
  return FlatObject::readFromFile<TPCdEdxCalibrationSplines>(inpf, name);
}

int TPCdEdxCalibrationSplines::writeToFile(TFile& outf, const char* name)
{
  /// write a class object to the file
  return FlatObject::writeToFile(*this, outf, name);
}

void TPCdEdxCalibrationSplines::setDefaultSplines()
{
  FlatObject::startConstruction();

  int buffSize = 0;
  int offsets1[FSplines];
  int offsets2[FSplines];

  auto defaultFnd2D = [&](float x1, float x2, float f[]) {
    f[0] = 1.f;
  };

  for (unsigned int ireg = 0; ireg < FSplines; ++ireg) {
    o2::gpu::Spline2D<float, 1> splineTmpqMax;
    splineTmpqMax.approximateFunction(0., 1., 0., 1., defaultFnd2D);
    mCalibSplinesqMax[ireg] = splineTmpqMax;
    buffSize = alignSize(buffSize, mCalibSplinesqMax[ireg].getBufferAlignmentBytes());
    offsets1[ireg] = buffSize;
    buffSize += mCalibSplinesqMax[ireg].getFlatBufferSize();
  }

  for (unsigned int ireg = 0; ireg < FSplines; ++ireg) {
    o2::gpu::Spline2D<float, 1> splineTmpqTot;
    splineTmpqTot.approximateFunction(0., 1., 0., 1., defaultFnd2D);
    mCalibSplinesqTot[ireg] = splineTmpqTot;
    buffSize = alignSize(buffSize, mCalibSplinesqTot[ireg].getBufferAlignmentBytes());
    offsets2[ireg] = buffSize;
    buffSize += mCalibSplinesqTot[ireg].getFlatBufferSize();
  }

  FlatObject::finishConstruction(buffSize);

  for (unsigned int i = 0; i < FSplines; i++) {
    mCalibSplinesqMax[i].moveBufferTo(mFlatBufferPtr + offsets1[i]);
  }
  for (unsigned int i = 0; i < FSplines; i++) {
    mCalibSplinesqTot[i].moveBufferTo(mFlatBufferPtr + offsets2[i]);
  }
}
#endif
