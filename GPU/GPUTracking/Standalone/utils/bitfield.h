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

/// \file bitfield.h
/// \author David Rohr

#ifndef Q_BITFIELD_H
#define Q_BITFIELD_H

#ifdef GPUCA_NOCOMPAT_ALLOPENCL
#include <type_traits>
#endif

template <class T, class S>
class bitfield
{
 public:
  bitfield(T v) : bits((S)v) {}
  bitfield(S v = 0) : bits(v) {}
  bitfield(const bitfield&) = default;
  bitfield& operator=(const bitfield&) = default;
  bitfield operator|(const bitfield v) const { return bits | v.bits; }
  bitfield operator|(const T v) const { return bits | static_cast<S>(v); }
  bitfield& operator|=(const bitfield v)
  {
    bits |= v.bits;
    return *this;
  }
  bitfield operator&(const bitfield v) const { return bits & v.bits; }
  bitfield operator&(const T v) const { return bits & static_cast<S>(v); }
  bitfield& operator&=(const bitfield v)
  {
    bits &= v.bits;
    return *this;
  }
  bitfield operator~() const { return ~bits; }
  bool operator==(const bitfield v) { return bits == v.bits; }
  bool operator==(const T v) { return bits == static_cast<S>(v); }
  bool operator!=(const bitfield v) { return bits != v.bits; }
  bool operator!=(const T v) { return bits != static_cast<S>(v); }
  bitfield& setBits(const bitfield v, bool w)
  {
    if (w) {
      bits |= v.bits;
    } else {
      bits &= ~v.bits;
    }
    return *this;
  }
  void clear() { bits = 0; }
  void set(S v) { bits = v; }
  void set(T v) { bits = static_cast<S>(v); }
  template <typename... Args>
  void set(T v, Args... args)
  {
    this->set(args...);
    bits |= static_cast<S>(v);
  }
  S get() const { return bits; }
  operator bool() const { return bits; }
  operator S() const { return bits; }
  bool isSet(const bitfield& v) const { return *this & v; }
  bool isSet(const S v) const { return bits & v; }

#ifdef GPUCA_NOCOMPAT_ALLOPENCL
  static_assert(std::is_integral<S>::value, "Storage type non integral");
  static_assert(sizeof(S) >= sizeof(T), "Storage type has insufficient capacity");
#endif

 private:
  S bits;
};

#endif
