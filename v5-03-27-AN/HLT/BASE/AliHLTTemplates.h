//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTTEMPLATES_H
#define ALIHLTTTEMPLATES_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTemplates.h
/// @author Matthias Richter
/// @date   2011-04-29
/// @brief  A collection of HLT template definitions
///

namespace HLT
{
  // unary predicate
  // checks if return value (type T) of a specified member function of a class C
  // instance has a certain value
  // used together with stl algorithms to extract specific predicate from the class
  template <class C, typename T>
  class AliHLTUnaryPredicate {
  public:
    /// constructor
    AliHLTUnaryPredicate(T (C::*pFct)() const, T value)
      : fpFct(pFct), fValue(value) { };

    /// copy contructor
    AliHLTUnaryPredicate(const AliHLTUnaryPredicate& p)
      : fpFct(p.fpFct), fValue(p.fValue) { };
    /// assignment operator
    AliHLTUnaryPredicate& operator=(const AliHLTUnaryPredicate& p) {
      fpFct=p.fpFct; fValue=p.fValue; return *this;
    }

    /// override operator== execute member function and compare result to value
    bool operator==(const C& object) const {
      if (!fpFct) return false;
      return (object.*fpFct)()==fValue;
    }

  private:
    /// standard contructor prohibited
    AliHLTUnaryPredicate();

    T (C::*fpFct)() const;   //! pointer to member function
    T  fValue;               //! value to match
  };

  template <class C, typename T>
  class AliHLTGetValue {
  public:
    /// constructor
    AliHLTGetValue(const C& object, T (C::*pFct)() const)
      : fObject(object), fpFct(pFct) { };

    /// override operator== execute member function and compare result to value
    operator T() const {
      return (fObject.*fpFct)();
    }

  private:
    /// standard contructor prohibited
    AliHLTGetValue();
    /// copy contructor prohibited
    AliHLTGetValue(const AliHLTGetValue&);
    /// assignment operator prohibited
    AliHLTGetValue& operator=(const AliHLTGetValue&);

    const T& fObject;       //! object
    T (C::*fpFct)() const;   //! pointer to member function
  };

  // operator== to be used as predicates for condition classes in stl algorithms
  // need to change the order of the to parameters in order to get
  // into the operator== of the condition class 
  template <class T, class C>
  bool operator==(const T& p, const C& c)
  {
    return c==p;
  }

  // operator== to be used as predicates for condition classes in stl algorithms
  // template for maps, use value of the pair in the condition
  // need to change the order of the to parameters in order to get
  // into the operator== of the condition class 
  template <class K, class T, class C>
  bool operator==(const std::pair<K, T>& p, const C& c)
  {
    return c==p.second;
  }

  // copy function for maps
  // stl algorithms can not be used here because pair.first can not be assigned
  // in maps but has to be used as index
  template <class InputIterator, class K, class V, class C>
  int copy_map_if ( InputIterator first, InputIterator last,
		    std::map<K, V>& result, const C& value )
  {
    for ( ; first != last; ++first)
      if ((*first == value)) result[first->first] = first->second;
    return 0;
  }

  // get the key of a pair
  class AliGetKey {
  public:
    template <typename T>
    typename T::first_type operator()(T pair) const {
      return pair.first;
    }
  };

} // end of namespace
#endif
