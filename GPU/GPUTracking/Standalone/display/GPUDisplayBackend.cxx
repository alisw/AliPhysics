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

/// \file GPUDisplayBackend.cxx
/// \author David Rohr

#include "GPUDisplayBackend.h"
#include "GPUDisplay.h"

using namespace GPUCA_NAMESPACE::gpu;

void* GPUDisplayBackend::OpenGLWrapper(void* ptr)
{
  GPUDisplayBackend* me = reinterpret_cast<GPUDisplayBackend*>(ptr);
  int retVal = me->OpenGLMain();
  if (retVal == -1) {
    me->InitGL(true);
  }
  return ((void*)(size_t)retVal);
}

void GPUDisplayBackend::HandleSendKey()
{
  if (mSendKey) {
    mDisplay->HandleSendKey(mSendKey);
    mSendKey = 0;
  }
}

void GPUDisplayBackend::HandleKeyRelease(unsigned char key) { mDisplay->HandleKeyRelease(key); }
int GPUDisplayBackend::DrawGLScene(bool mixAnimation, float animateTime) { return mDisplay->DrawGLScene(mixAnimation, animateTime); }
void GPUDisplayBackend::ReSizeGLScene(int width, int height)
{
  mDisplayHeight = height;
  mDisplayWidth = width;
  mDisplay->ReSizeGLScene(width, height);
}
int GPUDisplayBackend::InitGL(bool initFailure) { return mDisplay->InitGL(initFailure); }
void GPUDisplayBackend::ExitGL() { return mDisplay->ExitGL(); }
bool GPUDisplayBackend::EnableSendKey() { return true; }
