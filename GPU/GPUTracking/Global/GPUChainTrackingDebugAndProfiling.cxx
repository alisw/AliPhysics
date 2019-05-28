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

/// \file GPUChainTrackingDebugAndProfiling.cxx
/// \author David Rohr

#include "GPUChainTracking.h"
#ifdef GPUCA_TRACKLET_CONSTRUCTOR_DO_PROFILE
#include "bitmapfile.h"
#include <memory>
#endif
#define PROFILE_MAX_SIZE (100 * 1024 * 1024)

using namespace GPUCA_NAMESPACE::gpu;

static inline unsigned int RGB(unsigned char r, unsigned char g, unsigned char b) { return (unsigned int)r | ((unsigned int)g << 8) | ((unsigned int)b << 16); }

int GPUChainTracking::PrepareProfile()
{
#ifdef GPUCA_TRACKLET_CONSTRUCTOR_DO_PROFILE
  char* tmpMem = (char*)mRec->AllocateUnmanagedMemory(PROFILE_MAX_SIZE, GPUMemoryResource::MEMORY_GPU);
  processorsShadow()->tpcTrackers[0].mStageAtSync = tmpMem;
  runKernel<GPUMemClean16>({ BlockCount(), ThreadCount(), -1 }, nullptr, krnlRunRangeNone, krnlEventNone, tmpMem, PROFILE_MAX_SIZE);
#endif
  return 0;
}

int GPUChainTracking::DoProfile()
{
#ifdef GPUCA_TRACKLET_CONSTRUCTOR_DO_PROFILE
  std::unique_ptr<char[]> stageAtSync{ new char[PROFILE_MAX_SIZE] };
  mRec->GPUMemCpy(stageAtSync.get(), processorsShadow()->tpcTrackers[0].mStageAtSync, PROFILE_MAX_SIZE, -1, false);

  FILE* fp = fopen("profile.txt", "w+");
  FILE* fp2 = fopen("profile.bmp", "w+b");

  const int bmpheight = 8192;
  BITMAPFILEHEADER bmpFH;
  BITMAPINFOHEADER bmpIH;
  memset(&bmpFH, 0, sizeof(bmpFH));
  memset(&bmpIH, 0, sizeof(bmpIH));

  bmpFH.bfType = 19778; //"BM"
  bmpFH.bfSize = sizeof(bmpFH) + sizeof(bmpIH) + (ConstructorBlockCount() * ConstructorThreadCount() / 32 * 33 - 1) * bmpheight;
  bmpFH.bfOffBits = sizeof(bmpFH) + sizeof(bmpIH);

  bmpIH.biSize = sizeof(bmpIH);
  bmpIH.biWidth = ConstructorBlockCount() * ConstructorThreadCount() / 32 * 33 - 1;
  bmpIH.biHeight = bmpheight;
  bmpIH.biPlanes = 1;
  bmpIH.biBitCount = 32;

  fwrite(&bmpFH, 1, sizeof(bmpFH), fp2);
  fwrite(&bmpIH, 1, sizeof(bmpIH), fp2);

  int nEmptySync = 0;
  for (unsigned int i = 0; i < bmpheight * ConstructorBlockCount() * ConstructorThreadCount(); i += ConstructorBlockCount() * ConstructorThreadCount()) {
    int fEmpty = 1;
    for (unsigned int j = 0; j < ConstructorBlockCount() * ConstructorThreadCount(); j++) {
      fprintf(fp, "%d\t", stageAtSync[i + j]);
      int color = 0;
      if (stageAtSync[i + j] == 1) {
        color = RGB(255, 0, 0);
      }
      if (stageAtSync[i + j] == 2) {
        color = RGB(0, 255, 0);
      }
      if (stageAtSync[i + j] == 3) {
        color = RGB(0, 0, 255);
      }
      if (stageAtSync[i + j] == 4) {
        color = RGB(255, 255, 0);
      }
      fwrite(&color, 1, sizeof(int), fp2);
      if (j > 0 && j % 32 == 0) {
        color = RGB(255, 255, 255);
        fwrite(&color, 1, 4, fp2);
      }
      if (stageAtSync[i + j]) {
        fEmpty = 0;
      }
    }
    fprintf(fp, "\n");
    if (fEmpty) {
      nEmptySync++;
    } else {
      nEmptySync = 0;
    }
    (void)nEmptySync;
    // if (nEmptySync == GPUCA_SCHED_ROW_STEP + 2) break;
  }

  fclose(fp);
  fclose(fp2);
#endif
  return 0;
}
