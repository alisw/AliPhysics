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

/// \file vecpod.h
/// \author David Rohr

#include <vector>

template <class T>
struct vecpod_allocator {
  typedef T value_type;
  vecpod_allocator() noexcept : stdalloc() {}
  T* allocate(std::size_t n) { return stdalloc.allocate(n); }
  void deallocate(T* p, std::size_t n) { stdalloc.deallocate(p, n); }
  static void construct(T*) {}
  std::allocator<T> stdalloc;
};

template <class T>
using vecpod = typename std::vector<T, vecpod_allocator<T>>;
// template <class T> using vecpod = typename std::vector<T>;
