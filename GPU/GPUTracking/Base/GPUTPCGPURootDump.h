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

/// \file GPUTPCGPURootDump.h
/// \author David Rohr

#ifndef GPUTPCGPUROOTDUMP_H
#define GPUTPCGPUROOTDUMP_H

#if (!defined(GPUCA_STANDALONE) || defined(BUILD_QA)) && !defined(GPUCA_GPUCODE)
#include <TTree.h>
#include <TFile.h>
#include <TNtuple.h>

namespace GPUCA_NAMESPACE
{
namespace gpu
{
namespace
{
template <class S>
struct internal_Branch {
  template <typename... Args>
  static void Branch(S* p, Args... args)
  {
  }
};
template <>
struct internal_Branch<TTree> {
  template <typename... Args>
  static void Branch(TTree* p, Args... args)
  {
    p->Branch(args...);
  }
};
} // namespace

template <class T>
class GPUTPCGPURootDump
{
 public:
  GPUTPCGPURootDump() = delete;
  GPUTPCGPURootDump(const GPUTPCGPURootDump<T>&) = delete;
  GPUTPCGPURootDump<T> operator=(const GPUTPCGPURootDump<T>&) = delete;
  template <typename... Args>
  GPUTPCGPURootDump(const char* filename, Args... args)
  {
    fFile = new TFile(filename, "recreate");
    fTree = new T(args...);
  }

  ~GPUTPCGPURootDump()
  {
    fTree->Write();
    fFile->Write();
    fFile->Close();
    delete fFile;
  }

  template <typename... Args>
  void Fill(Args... args)
  {
    fTree->Fill(args...);
  }
  template <typename... Args>
  void Branch(Args... args)
  {
    internal_Branch<T>::Branch(fTree, args...);
  }

 private:
  TFile* fFile = nullptr;
  T* fTree = nullptr;
};
#else
template <class T>
class GPUTPCGPURootDump
{
 public:
  GPUTPCGPURootDump() = delete;
  template <typename... Args>
  GPUTPCGPURootDump(const char* filename, Args... args)
  {
  }
  template <typename... Args>
  void Fill(Args... args)
  {
  }
  template <typename... Args>
  void Branch(Args... args)
  {
  }

 private:
  void *a, *b;
};
}
}
#endif

#endif
