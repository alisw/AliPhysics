///
/// \file AliFemto/AliFemtoConfigObject.h
///

#pragma once

#ifndef ALIFEMTOCONFIGOBJECT_H
#define ALIFEMTOCONFIGOBJECT_H

#include <map>
#include <list>
#include <string>
#include <vector>
#include <ostream>
#include <utility>
#include <iostream>
#include <iterator>

#include <TBuffer.h>
#include <TString.h>
#include <TObjString.h>
#include <TText.h>
#include <TPaveText.h>


#if __cplusplus >= 201103L
#define ENABLE_MOVE_SEMANTICS 1
#endif


#define FORWARD_STANDARD_TYPES(__FUNC)          \
  __FUNC(BoolValue_t, kBOOL, fValueBool);       \
  __FUNC(IntValue_t, kINT, fValueInt);          \
  __FUNC(FloatValue_t, kFLOAT, fValueFloat);    \
  __FUNC(StringValue_t, kSTRING, fValueString); \
  __FUNC(ArrayValue_t, kARRAY, fValueArray);    \
  __FUNC(MapValue_t, kMAP, fValueMap);          \
  __FUNC(RangeValue_t, kRANGE, fValueRange);    \
  __FUNC(RangeListValue_t, kRANGELIST, fValueRangeList);


/// \class AliFemtoConfigObject
/// \brief Polymorphic type for easily specifying and storing parameters
///
/// ConfigObject is a JSON-like data format used
///
/// #### Numbers
///
/// Basic ints and floats are supported.
/// Only base-10 integers are supported (no hex-numbers, yet)
///
///
/// #### String
///
/// Strings start and stop with a single quote, for easily embedding with
/// within normal C++ strings.
/// To prevent the hassle of dealing with escaping quotes, these
/// config-objects have no concept of escape symbols. This means there
/// is no way to have a single or double quote.
///
///
/// #### Range & RangeList
///
/// It is often useful to specify two numerical limits on a value; to
/// simplify handling this special-case, ConfigObject has a 'Range'
/// type (``std::pair<double,double>``) that can be specified by simply
/// separating two numbres with a colon (ex, -1:1 or 2.71:3.14).
///
/// Groups of ranges can be set using multiple colon-separated numbers
/// surrounded by parentheses "()". These groups define non-overlapping
/// continuous ranges, so `(0:5:10:15:20)` would be a list of 4 ranges:
/// `0-5, 5-10, 10-15, 15-20`.
/// Non-continuous ranges are suppored using comma syntax:
///   (0:5:10,15:20) would be the same as above without 10-15
///
///
/// #### Containers
///
/// The strength of this class lies in the ability to store other ConfigObjects
/// within.
/// There are two container types: lists "[1,'a',3:4]" and maps ("{ a: 1, b: ['x', 2.3] }")
/// Map keys are NOT surrounded by quotes - use only their names.
/// Map keys provide a shorthand for making maps within maps using '.' or '/' as
/// key separators:
///   `{a.b: "c"} == {a/b: "c"} == {a: {b: "c"}}; // true`
///
/// Slash notation might be removed...
///
///
/// Usage
/// -----
///
///
///
/// Examples
/// --------
///
/// ```cpp
/// AliFemtoConfigObject hist_cfg = AliFemtoConfigObject::Parse(
///     "{ name: 'h', title: 'Example histogram', bincount: 120, maxmin: -3:5}");
///
/// hist_cfg.find<std::string>("name", "DEFAULT_NAME"); // "h"
/// hist_cfg.find<std::string>("hist_title"); // "" (not found)
/// hist_cfg.find_int("bincount");  // 120
///
/// double start = 0, stop = 1; // default values
/// // only sets start & stop if specified and the correct type!
/// if (!hist_cfg.load_range("binlimits", start, stop)) {
///
//
///    std::cerr << "Warning - Histogram limits not set"
///                 " - using defaults (" << start << ", " << stop << ")\n";
/// }

///
/// ```
///
///
/// \author Andrew Kubera, The Ohio State University, andrew.kubera@cern.ch
///
class AliFemtoConfigObject : public TObject {
public:

  /// \defgroup Associated Types
  /// @{
#if __cplusplus < 201103L
  typedef Bool_t BoolValue_t;
  typedef double FloatValue_t;
  typedef Long64_t IntValue_t;
  typedef std::string StringValue_t;
  typedef std::vector<AliFemtoConfigObject> ArrayValue_t;
  typedef std::string Key_t;
  typedef std::map<Key_t, AliFemtoConfigObject> MapValue_t;
  typedef std::pair<FloatValue_t, FloatValue_t> RangeValue_t;
  typedef std::vector<RangeValue_t> RangeListValue_t;
#else
  using BoolValue_t = Bool_t;
  using FloatValue_t = double;
  using IntValue_t = Long64_t;
  using StringValue_t = std::string;
  using ArrayValue_t = std::vector<AliFemtoConfigObject>;
  using Key_t = std::string;
  using MapValue_t = std::map<Key_t, AliFemtoConfigObject>;
  using RangeValue_t = std::pair<FloatValue_t, FloatValue_t>;
  using RangeListValue_t = std::vector<RangeValue_t>;
#endif
    /// @}

  /// class for visualizing on a canvas (TPad)
  class Painter;

  static const ULong_t EMPTY_ARRAY_HASH = 0x742e2b13ab659L,
                       EMPTY_MAP_HASH = 0x66fa474be7ae5L,
                       EMPTY_RANGELIST_HASH = 0x905193fc5bfc52L;

protected:

  /// finite number of posslbe types
  enum TypeTagEnum_t {
    kEMPTY,
    kBOOL,
    kINT,
    kFLOAT,
    kSTRING,
    kARRAY,
    kMAP,
    kRANGE,
    kRANGELIST
  };

  /// static mapping from enum to type
  // template <TypeTagEnum_t> struct type_from_enum { };

public:
  /// Create object by parsing cstring
  static AliFemtoConfigObject Parse(const char*);

  /// Create object by parsing string
  static AliFemtoConfigObject Parse(const std::string &);

  /// Create object by parsing TString
  static AliFemtoConfigObject Parse(const TString &);

  /// Create object by parsing TString
  static AliFemtoConfigObject Parse(const TObjString &);

  /// Create object from string; if a map, also parse defaults and
  /// insert any values missing in object with defaults.
  static AliFemtoConfigObject ParseWithDefaults(const std::string &, const std::string &defaults);

  /// Create object from string; if a map, also parse defaults and
  /// insert any values missing in object with defaults.
  static AliFemtoConfigObject ParseWithDefaults(const std::string &, const AliFemtoConfigObject &defaults);

  /// Create object in empty state
  AliFemtoConfigObject():
    TObject()
    , fTypeTag(kEMPTY)
    , fPainter(NULL)
    {
    }

  /// Construct from TObject
  AliFemtoConfigObject(const TObject *ptr);
  AliFemtoConfigObject(const TObject &orig);

  /// Calls destructor on internal data
  virtual ~AliFemtoConfigObject();

  /// Copy value
  AliFemtoConfigObject(AliFemtoConfigObject const &orig):
    TObject(orig)
    , fTypeTag(kEMPTY)
    , fPainter(NULL)
    {
      _CopyConstructValue(orig);
    }

  /// Assign value
  AliFemtoConfigObject& operator=(AliFemtoConfigObject const &);

  // I don't know why this needs to be specified explicitly
  // instead of the macro-generated ones below
  AliFemtoConfigObject(const char * v):
    fTypeTag(kSTRING), fValueString(v), fPainter(nullptr) { }

  // define constructors
  #define ConstructorDef(__type, __tag, __dest) \
    AliFemtoConfigObject(const __type &v):      \
      fTypeTag(__tag), __dest(v), fPainter(nullptr) { }

    // ConstructorDef(char*, kSTRING, fValueString);

    // constructors with all associated value-types (IntValue_t, MapValue_t, etc..)
    FORWARD_STANDARD_TYPES(ConstructorDef)

    //ConstructorDef(TString, kSTRING, fValueString);   // implicit cast - not allowed by clang 9.0, must be declared explicitly
    ConstructorDef(int, kINT, fValueInt);
    ConstructorDef(long, kINT, fValueInt);
    ConstructorDef(unsigned, kINT, fValueInt);
    ConstructorDef(unsigned long, kINT, fValueInt);
    ConstructorDef(float, kFLOAT, fValueFloat);
  #undef ConstructorDef
  // Explicit constructor for TString - needed for clang 9.0
  AliFemtoConfigObject(const TString &v) : fTypeTag(kSTRING), fValueString(v.Data()), fPainter(nullptr) {}

  // iterator-based constructors
  #define ConstructorDefItor(__type, __tag, __dest) \
    AliFemtoConfigObject(const __type::iterator a, const __type::iterator b): fTypeTag(__tag), __dest(a,b), fPainter(nullptr) { }

    ConstructorDefItor(std::string, kSTRING, fValueString);
    ConstructorDefItor(std::list<AliFemtoConfigObject>, kARRAY, fValueArray);
  #undef ConstructorDefItor

  /// Construct Range from two floating points
  AliFemtoConfigObject(const Float_t a, const Float_t b): fTypeTag(kRANGE), fValueRange(a,b), fPainter(nullptr) { }


// moveable heap-allocated data
#ifdef ENABLE_MOVE_SEMANTICS

  /// Move constructor
  AliFemtoConfigObject(AliFemtoConfigObject &&orig):
    TObject(orig)
    , fTypeTag(kEMPTY)
    , fPainter(nullptr)
    {
      _MoveConstructValue(std::move(orig));
      /*
      std::swap(fPainter, orig.fPainter);
      if (fPainter) {
        fPainter->ResetData(this);
      }
      */
    }

  /// Move assignment
  AliFemtoConfigObject& operator=(AliFemtoConfigObject &&rhs)
    {
      // not same type - destroy data and move
      if (rhs.fTypeTag != fTypeTag) {
        _DeleteValue();
        _MoveConstructValue(std::move(rhs));
      }
      // same type - can move-assign data
      else {
        switch (rhs.fTypeTag) {
          case kEMPTY: break;
          case kBOOL: fValueBool = std::move(rhs.fValueBool); break;
          case kFLOAT: fValueFloat = std::move(rhs.fValueFloat); break;
          case kINT: fValueInt = std::move(rhs.fValueInt); break;
          case kSTRING: fValueString = std::move(rhs.fValueString); break;
          case kARRAY: fValueArray = std::move(rhs.fValueArray); break;
          case kMAP: fValueMap = std::move(rhs.fValueMap); break;
          case kRANGE: fValueRange = std::move(rhs.fValueRange); break;
          case kRANGELIST: fValueRangeList = std::move(rhs.fValueRangeList); break;
        }

        rhs._DeleteValue();
      }

      // unsure if this is proper behavior
      /*
      delete fPainter;
      fPainter = nullptr;
      std::swap(fPainter, rhs.fPainter);
      if (fPainter) {
        fPainter->ResetData(this);
      }
      */

      return *this;
    }

  #define ConstructorDef(__type, __tag, __dest) AliFemtoConfigObject(__type &&v): fTypeTag(__tag), __dest(std::move(v)), fPainter(nullptr) { }
    FORWARD_STANDARD_TYPES(ConstructorDef)
  #undef ConstructorDef

#endif // move-semantics

  typedef std::pair<float,float> pair_of_floats;
  typedef std::pair<int,int> pair_of_ints;

  /// \class BuildStruct
  /// \brief helper class for building a mapping object
  class BuildMap {
  public:
    BuildMap(): fMap() {}

#ifdef ENABLE_MOVE_SEMANTICS

    #define IMPL_BUILDITEM(__type, __a, __b) \
      BuildMap& operator()(const Key_t &key, const __type& val) { fMap.emplace(key, val); return *this; }
    #define IMPL_BUILDITEM_CASTED(__type, __savedtype) \
      BuildMap& operator()(const Key_t &key, const __type& val) { fMap.emplace(key, static_cast<__savedtype>(val)); return *this; } \
      BuildMap& operator()(const Key_t &key, __type&& val) { fMap.emplace(key, std::move(static_cast<__savedtype>(val))); return *this; } \

#else

    #define IMPL_BUILDITEM(__type, __a, __b) \
      BuildMap& operator()(const Key_t &key, const __type& val) { fMap[key] = val; return *this; }
    #define IMPL_BUILDITEM_CASTED(__type, __savedtype) \
      BuildMap& operator()(const Key_t &key, const __type& val) { fMap[key] = static_cast<__savedtype>(val); return *this; }

#endif

    FORWARD_STANDARD_TYPES(IMPL_BUILDITEM)

    IMPL_BUILDITEM_CASTED(Float_t, FloatValue_t);
    IMPL_BUILDITEM_CASTED(Int_t, IntValue_t);
    IMPL_BUILDITEM_CASTED(long, IntValue_t);
    IMPL_BUILDITEM_CASTED(UInt_t, IntValue_t);
    IMPL_BUILDITEM_CASTED(ULong64_t, IntValue_t);
    IMPL_BUILDITEM_CASTED(pair_of_ints, RangeValue_t);
    IMPL_BUILDITEM_CASTED(pair_of_floats, RangeValue_t);

    IMPL_BUILDITEM(AliFemtoConfigObject, 0, 0);

    #undef IMPL_BUILDITEM
    #undef IMPL_BUILDITEM_CASTED


    // -- custom operator() methods --

#ifdef ENABLE_MOVE_SEMANTICS
    BuildMap& operator()(const Key_t &key, const TString &val) { fMap.emplace(key, val.Data()); return *this; };
    BuildMap& operator()(const Key_t &key, const char* val) { fMap.emplace(key, StringValue_t(val)); return *this; }
#else
    BuildMap& operator()(const Key_t &key, const TString &val) { fMap[key] = val.Data(); return *this; };
    BuildMap& operator()(const Key_t &key, const char* val) { fMap[key] = StringValue_t(val); return *this; }
#endif

    operator AliFemtoConfigObject()
      { return into(); }

    AliFemtoConfigObject into()
      { return AliFemtoConfigObject(std::move(fMap)); }

    AliFemtoConfigObject into() const
      { return AliFemtoConfigObject(fMap); }


  protected:
    MapValue_t fMap;
  };


  // ========================
  //      C++ Operators
  // ========================

  /// Comparison based on the hash of the object.
  ///
  /// Note this guarantees transativity of config objects, but holds
  /// NO gurantees on the ordering of the values within the config
  /// object; if two objects containing integers, cfg1 & cfg2, are
  /// compared and cfg1 < cfg2 is true, the integer value of cfg1 is
  /// not necessarily less than the one in cfg2 (but is guaranteed to
  /// be not equal)
  ///
  /// This operator should mostly be used for ordering purposes (such
  /// as storage classes), and not for
  /// values extracted for comparison.
  ///
  /// Note that this may be an expensive operation for complex
  /// structures like Map and List, as a hash is computed for every
  /// contained value.
  ///
  bool operator<(const AliFemtoConfigObject &rhs) const { return Hash() < rhs.Hash(); };

  bool operator==(const AliFemtoConfigObject &rhs) const;
  bool operator!=(const AliFemtoConfigObject &rhs) const { return !(*this == rhs); }
  #define IMPL_EQUALITY(__type, __tag, __target) \
    bool operator==(const __type &rhs) const { return fTypeTag == __tag && __target == rhs; } \
    bool operator!=(const __type &rhs) const { return !operator==(rhs); }

    FORWARD_STANDARD_TYPES(IMPL_EQUALITY)
    IMPL_EQUALITY(TString, kSTRING, fValueInt); // fValueString);

  #undef IMPL_EQUALITY


  // ========================
  //   ConfigObject Methods
  // ========================

  bool is_empty() const { return fTypeTag == kEMPTY; }
  bool is_bool() const { return fTypeTag == kBOOL; }
  bool is_float() const { return fTypeTag == kFLOAT; }
  bool is_int() const { return fTypeTag == kINT; }
  bool is_num() const { return fTypeTag == kINT || fTypeTag == kFLOAT; }
  bool is_str() const { return fTypeTag == kSTRING; }
  bool is_array() const { return fTypeTag == kARRAY; }
  bool is_map() const { return fTypeTag == kMAP; }
  bool is_range() const { return fTypeTag == kRANGE; }
  bool is_rangelist() const { return fTypeTag == kRANGELIST; }

  /// Return integer value - no typecheck
  BoolValue_t as_bool() const { return fValueBool; }
  IntValue_t as_int() const { return fValueInt; }
  FloatValue_t as_float() const { return fValueFloat; }
  RangeValue_t as_range() const { return fValueRange; }
  RangeListValue_t as_rangelist() const
    {
      return is_rangelist() ? fValueRangeList
           : is_range() ? RangeListValue_t(1, fValueRange)
           : RangeListValue_t();
    }

  // template <typename BoolType> bool load_bool(BoolType &v) const { return is_bool() ? v = fValueBool, true : false; }
  // template <typename StringType> bool load_str(StringType &v) const { return is_str() ? v = fValueString, true : false; }
  // template <typename FloatType> bool load_float(FloatType &v) const { return is_float() ? v = fValueFloat, true : false; }
  // template <typename IntType> bool load_int(IntType &v) const { return is_int() ? v = fValueInt, true : false; }
  // template <typename RangeType> bool load_range(RangeType &v) const { return is_range() ? v = fValueRange, true : false; }
  // template <typename ArrayType> bool load_array(ArrayType &v) const { return is_array() ? v = fValueArray, true : false; }
  // template <typename ArrayType> bool load_into_array(ArrayType &v) const { return is_array() ? std::copy(fValueArray.begin(), fValueArray.end(), v), true : false; }

  bool load_bool(BoolValue_t &b) const { return is_bool() ? b = fValueBool, true : false; }
  bool load_float(FloatValue_t &f) const { return is_float() ? f = fValueFloat, true : false; }
  bool load_float(float &f) const { return is_float() ? f = fValueFloat, true : false; }
  bool load_int(IntValue_t &i) const { return is_int() ? i = fValueInt, true : false; }
  bool load_int(int &i) const { return is_int() ? i = static_cast<int>(fValueInt), true : false; }
  bool load_int(unsigned int &i) const { return is_int() ? i = static_cast<unsigned int>(fValueInt), true : false; }
  bool load_uint(ULong_t &i) const { return is_int() ? i = static_cast<ULong_t>(fValueInt), true : false; }
  bool load_uint(ULong64_t &i) const { return is_int() ? i = static_cast<ULong64_t>(fValueInt), true : false; }
  bool load_num(FloatValue_t &f) const { return is_int() ? f = fValueInt, true : load_float(f); }
  bool load_str(std::string &s) const { return is_str() ? s = fValueString, true : false; }
  bool load_str(TString &s) const { return is_str() ? s = fValueString, true : false; }
  bool load_array(ArrayValue_t &a) const { return is_array() ? a = fValueArray, true : false; }
  bool load_map(MapValue_t &m) const { return is_map() ? m = fValueMap, true : false; }

  bool load_range(RangeValue_t &r) const { return is_range() ? r = fValueRange, true : false; }
  bool load_range(pair_of_floats &r) const { return is_range() ? r = fValueRange, true : false; }
  bool load_range(pair_of_ints &r) const { return is_range() ? r = fValueRange, true : false; }
  bool load_range(typename RangeValue_t::first_type &a, typename RangeValue_t::second_type &b) const
    { return is_range() ? a = fValueRange.first, b = fValueRange.second, true : false; }
  bool load_range(Float_t &a, Float_t &b) const
    { return is_range() ? a = fValueRange.first, b = fValueRange.second, true : false; }
  bool load_range(int &a, int &b) const
    { return is_range() ? a = fValueRange.first, b = fValueRange.second, true : false; }

  bool load_rangelist(RangeListValue_t &r) const
    { return is_rangelist() ? r = fValueRangeList, true : false; }

  /// same as load_rangelist but interprets single range as rangelist
  /// of length 1
  bool load_ranges(RangeListValue_t &r) const
    { return is_range() ? r.clear(), r.push_back(fValueRange), true : load_rangelist(r); }


  /// Return the commonname of this object's contained type; usefull for debugging
  TString name_of_type() const
    { return NameFromtype(fTypeTag); }

  /// \defgroup Find&Load Methods
  /// @{
  /// Copies item identified by *key*, returns true if found
  #define IMPL_FINDANDLOAD(__dest_type, __load, __tag, __source)    \
    bool find_and_load(const Key_t &key, __dest_type &dest) const { \
      if (!is_map()) { return false; }                              \
      const AliFemtoConfigObject* found = find(key);                \
      if (!found) { return false; }                                 \
      return found-> __load (dest); }

    // FORWARD_STANDARD_TYPES(IMPL_FINDANDLOAD)

    IMPL_FINDANDLOAD(BoolValue_t, load_bool, kBOOL, fValueBool);
    IMPL_FINDANDLOAD(IntValue_t, load_int, kINT, fValueInt);
    IMPL_FINDANDLOAD(FloatValue_t, load_float, kFLOAT, fValueFloat);
    IMPL_FINDANDLOAD(StringValue_t, load_str, kSTRING, fValueString);
    IMPL_FINDANDLOAD(ArrayValue_t, load_array, kARRAY, fValueArray);
    IMPL_FINDANDLOAD(MapValue_t, load_map, kMAP, fValueMap);
    IMPL_FINDANDLOAD(RangeValue_t, load_range, kRANGE, fValueRange);
    IMPL_FINDANDLOAD(RangeListValue_t, load_rangelist, kRANGELIST, fValueRangeList);

    IMPL_FINDANDLOAD(pair_of_floats, load_range, kRANGE, fValueRange);
    IMPL_FINDANDLOAD(pair_of_ints, load_range, kRANGE, fValueRange);
    // IMPL_FINDANDLOAD(unsigned int, load_int, kINT, fValueInt);
    IMPL_FINDANDLOAD(ULong64_t, load_uint, kINT, fValueInt);
    // IMPL_FINDANDLOAD(int, load_int, kINT, fValueInt);
    IMPL_FINDANDLOAD(Float_t, load_float, kFLOAT, fValueFloat);
    IMPL_FINDANDLOAD(TString, load_str, kSTRING, fValueString);

  #undef IMPL_FINDANDLOAD

  bool find_and_load(const Key_t &key, AliFemtoConfigObject &dest) const
    {
      if (is_map()) {
        if (const AliFemtoConfigObject *found = find(key)) {
          dest = *found;
          return true;
        }
      }
      return false;
    }

  #define IMPL_GET(__type, _ignored, __ignored) \
    __type get(const Key_t &key, const __type &defval) const \
      { __type result(defval); find_and_load(key, result);  \
        return result; }

  FORWARD_STANDARD_TYPES(IMPL_GET)
  IMPL_GET(pair_of_floats, kRANGE, fValueRange);
  IMPL_GET(pair_of_ints, kRANGE, fValueRange);
  // IMPL_GET(int, kINT, fValueInt);
  IMPL_GET(ULong64_t, kINT, fValueInt);
  IMPL_GET(Float_t, kFLOAT, fValueFloat);
  IMPL_GET(TString, kSTRING, fValueString);

  #undef IMPL_GET

  AliFemtoConfigObject get_object(const Key_t &key) const
    {
      AliFemtoConfigObject result;
      find_and_load(key, result);
      return result;
    }

  double get_num(const Key_t &key, const double defval=NAN) const
    {
      double result = defval;

      if (is_map()) {
        const AliFemtoConfigObject *found = find(key);
        if (found) {
          if (found->is_float()) {
            result = found->fValueFloat;
          }
          else if (found->is_int()) {
            result = found->fValueInt;
          }
        }
      }

      return result;
    }

  TString get_str(const Key_t &key, const TString &defval="") const
    {
      TString result = defval;

      if (is_map()) {
        const AliFemtoConfigObject *found = find(key);
        if (found && found->is_str()) {
          result = found->fValueString;
        }
      }

      return result;
    }

  #if __cplusplus < 201103L

  #define IMPL_INSERT(__value_type, _ignored, __ignored)       \
    void insert(const Key_t &key, const __value_type &value) { \
      if (!is_map()) { return; }                               \
      fValueMap[key] = AliFemtoConfigObject(value); }

  #define IMPL_CASTED_INSERT(__value_type, __stored_type, _)     \
    void insert(const Key_t &key, const __value_type &value) {   \
      if (!is_map()) { return; }                                 \
      fValueMap[key] = static_cast<__stored_type>(value); }

  #else

  #define IMPL_INSERT(__value_type, _ignored, __ignored)       \
    void insert(const Key_t &key, const __value_type &value) { \
      if (!is_map()) { return; }                               \
      fValueMap.emplace(key, value); }

  #define IMPL_CASTED_INSERT(__value_type, __stored_type, _)     \
    void insert(const Key_t &key, const __value_type &value) {   \
      if (!is_map()) { return; }                                 \
      fValueMap.emplace(key, static_cast<__stored_type>(value)); }

  #endif

    FORWARD_STANDARD_TYPES(IMPL_INSERT)

    IMPL_INSERT(pair_of_floats, kRANGE, fValueRange);
    IMPL_INSERT(pair_of_ints, kRANGE, fValueRange);
    IMPL_INSERT(int, kINT, fValueInt);
    IMPL_INSERT(Float_t, kFLOAT, fValueFloat);
    IMPL_INSERT(TString, kSTRING, fValueString);
    IMPL_INSERT(AliFemtoConfigObject, void, void);

    IMPL_CASTED_INSERT(ULong_t, IntValue_t, 0);
    IMPL_CASTED_INSERT(ULong64_t, IntValue_t, 0);

  #undef IMPL_INSERT
  #undef IMPL_CASTED_INSERT

  #ifdef ENABLE_MOVE_SEMANTICS
  #define IMPL_INSERT(__value_type, _i, __i)              \
    void insert(const Key_t &key, __value_type &&value) { \
      if (!is_map()) { return; }                          \
      fValueMap.emplace(key, std::move(value)); }

  FORWARD_STANDARD_TYPES(IMPL_INSERT)

  #undef IMPL_INSERT
  #endif

  /// Checks for existence of key in this object.
  /// Returns false if this is not a map or key is not found.
  /// Does not check sub-objects (but this SHOULD be done at
  /// some point)
  bool has_key(const Key_t key) const
  {
    return fTypeTag == kMAP
           && fValueMap.find(key) != fValueMap.end();
  }

  /// Remove and returns pointer to object at *index* if an array
  ///
  /// Delete the pointer after use!
  ///
  /// Index starts at 0. Negative numbers wrap around to the back, so '-1'
  /// removes the last element (default
  ///
  /// If not array, or index out of bounds returns nullptr
  AliFemtoConfigObject* pop(Int_t index=-1);

  /// Removes and returns pointer to object associated with *key*
  ///
  /// If not a struct nor has *key*, return nullptr.
  AliFemtoConfigObject* pop(const Key_t &key);

  /// Removes and returns pointer to object associated with *key*,
  /// or pointer to ConfigObject constructed from default if not
  /// found or not a struct.
  template <typename ReturnType>
  AliFemtoConfigObject* pop(const Key_t &key, const ReturnType &default_)
    {
      if (auto *obj = pop(key)) {
        return obj;
      }
      return new AliFemtoConfigObject(default_);
    }


  /// If object is map and key is present, returns pointer to subobject
  /// otherwise return NULL
  ///
  /// Returned pointer is still owned by this object; do NOT delete.
  ///
  /// This should return std::weak_ptr<AliFemtoConfigObject> whenever
  /// AliPhysics fully supports C++11
  ///
  AliFemtoConfigObject* find(const Key_t &key);

  /// Apply find, but return a pointer to a constant ConfigObject
  const AliFemtoConfigObject* find(const Key_t &key) const;


  /// Return copy of object with a key removed from map.
  /// If not a map, or key not present, object is simply copied
  AliFemtoConfigObject WithoutKey(const Key_t &key) const;

  /// Return copy of object with given keys removed
  AliFemtoConfigObject WithoutKeys(const std::vector<Key_t>&) const;

  #define IMPL_POP_ITEM(__dest_type, __name)           \
    __dest_type __name(const Key_t &key, const __dest_type &_def) \
      { __dest_type res(_def); pop_and_load(key, res); return res; }

    IMPL_POP_ITEM(StringValue_t, pop_str);
    IMPL_POP_ITEM(TString, pop_str);
    IMPL_POP_ITEM(IntValue_t, pop_int);
    IMPL_POP_ITEM(unsigned int, pop_uint);
    IMPL_POP_ITEM(FloatValue_t, pop_float);
    IMPL_POP_ITEM(MapValue_t, pop_map);
    IMPL_POP_ITEM(ArrayValue_t, pop_array);
    IMPL_POP_ITEM(BoolValue_t, pop_bool);
    IMPL_POP_ITEM(RangeValue_t, pop_range);
    IMPL_POP_ITEM(pair_of_floats, pop_range);
    IMPL_POP_ITEM(pair_of_ints, pop_range);

    StringValue_t pop_str(const Key_t &key, const char* _default)
      { StringValue_t res(_default); pop_and_load(key, res); return res; }

    /// Pop generic number (int or float) out of object
    FloatValue_t pop_num(const Key_t &key, FloatValue_t _default=NAN)
      {
        double res(_default);
        if (const AliFemtoConfigObject *found = find(key)) {
          if (found->is_float()) {
            res = found->as_float();
          }
          else if (found->is_int()) {
            res = found->as_int();
          }
          else {
            return _default;
          }
          pop(key);
        }
        return res;
      }

  #undef IMPL_POP_ITEM


  /// \defgroup Pop&Load Methods
  /// @{
  /// Removes item identified by *key*

  #define IMPL_POPANDLOAD(__dest_type, __tag, __source)      \
    bool pop_and_load(const Key_t &key, __dest_type &dest) { \
      if (!is_map()) { return false; }                       \
      MapValue_t::const_iterator found = fValueMap.find(key);\
      if (found == fValueMap.cend()) { return false; }       \
      if (found->second.fTypeTag != __tag) { return false; } \
      dest = found->second. __source;                        \
      fValueMap.erase(found); return true; }

    FORWARD_STANDARD_TYPES(IMPL_POPANDLOAD)

    IMPL_POPANDLOAD(pair_of_floats, kRANGE, fValueRange);
    IMPL_POPANDLOAD(pair_of_ints, kRANGE, fValueRange);
    IMPL_POPANDLOAD(int, kINT, fValueInt);
    IMPL_POPANDLOAD(unsigned int, kINT, fValueInt);
    IMPL_POPANDLOAD(Float_t, kFLOAT, fValueFloat);
    IMPL_POPANDLOAD(TString, kSTRING, fValueString);

  #undef IMPL_POPANDLOAD

  /// Find object by key and move value into another object
  ///
  /// Returns true if the object was found; false if not found, or
  /// this object is not a map.
  ///
  bool pop_and_load(const Key_t &key, AliFemtoConfigObject &dest)
  {
    if (is_map()) {
      MapValue_t::const_iterator found = fValueMap.find(key);
      if (found != fValueMap.cend()) {
        dest = found->second;
        fValueMap.erase(found);
        return true;
      }
    }
    return false;
  }


  /// @}

  /// \class list_iterator
  /// \brief Iterates over object if it is a list
  ///
  /// Use this class by calling `list_begin()` and `list_end()` on a
  /// config object.
  /// If the object is not a list, then list_begin() will return the
  /// same value as list_end(), and there is no iteration (using
  /// standard iteration rules)
  ///
  class list_iterator : public std::iterator<std::bidirectional_iterator_tag, ArrayValue_t::value_type> {
    AliFemtoConfigObject *fParent;
    bool fIsArray;
    ArrayValue_t::iterator fInternal;
    friend class AliFemtoConfigObject;

    list_iterator(AliFemtoConfigObject *obj, ArrayValue_t::iterator it):
      fParent(obj), fIsArray(true), fInternal(it) {}

  public:
    /// construct from config object
    list_iterator(AliFemtoConfigObject &obj)
      : fParent(&obj), fIsArray(obj.is_array()), fInternal()
    {
      if (fIsArray) {
        fInternal = obj.fValueArray.begin();
      }
    }

    list_iterator(const list_iterator &orig)
      : fParent(orig.fParent), fIsArray(orig.fParent), fInternal()
    {
      if (fIsArray) {
        fInternal = orig.fInternal;
      }
    }

    value_type& operator*()
      { return *fInternal; }

    list_iterator& operator++()
      { ++fInternal; return *this; }

    list_iterator operator++(int)
    {
      list_iterator temp(*this);
      ++fInternal;
      return temp;
    }

    bool operator!=(const list_iterator &rhs) const {
      return rhs.fParent != fParent                  // if we don't have same parent - different
          || fIsArray ? (fInternal != rhs.fInternal) // if array, compare internal iterator
                      : false;                       // if not array, we are equal
    }
  };

  /// Iterator to first element of list
  ///
  /// If this value is not a list the value returned is equal to
  /// `list_end`, making it safe to use in a for look without
  /// doing a type check.
  ///
  list_iterator list_begin()
    { return list_iterator(*this); }

  /// Iterator to 'end' of list
  ///
  /// If this config object is a list, the iterator returned is
  /// equivalent to the value returned by a standard `std::end()`
  /// function call; otherwise it returns a unique value which should
  /// be compared to the value returned by `list_begin`.
  ///
  list_iterator list_end()
  {
    if (is_array()) {
      return list_iterator(this, fValueArray.end());
    }
    return list_iterator(*this);
  }

  /// \class iterator_over_list
  /// \brief Returned by method `items_in_list` for c++11 foreach
  ///        loop syntax
  ///
  /// ```cpp
  /// for (const auto &obj : container.items_in_list()) {
  ///   ...
  /// }
  /// ```
  ///
  class iterator_over_list {
    AliFemtoConfigObject *fParent;
  public:

    iterator_over_list(AliFemtoConfigObject &obj): fParent(&obj) {}
    list_iterator begin() { return fParent->list_begin(); }
    list_iterator end() { return fParent->list_end(); }
  };

  iterator_over_list items_in_list()
    { return iterator_over_list(*this); }

  /// \class map_iterator
  /// \brief Iterates over object if it is a map
  ///
  /// Use this class by calling `map_begin()` and `map_end()` on a
  /// config object.
  /// If the object is not a list, then map_begin() will return the
  /// same value as map_end(), and there is no iteration (using
  /// standard iteration rules)
  ///
  class map_iterator : public std::iterator<std::bidirectional_iterator_tag, MapValue_t::value_type> {
    AliFemtoConfigObject *fParent;
    bool fIsMap;
    MapValue_t::iterator fInternal;
    friend class AliFemtoConfigObject;

    map_iterator(AliFemtoConfigObject *obj, MapValue_t::iterator it)
      : fParent(obj), fIsMap(true), fInternal(it)
    { }

  public:
    /// construct from config object
    map_iterator(AliFemtoConfigObject &obj)
      : fParent(&obj), fIsMap(obj.is_map()), fInternal()
    {
      if (fIsMap) {
        fInternal = obj.fValueMap.begin();
      }
    }

    value_type& operator*()
      { return *fInternal; }

    // prefix-increment (++it)
    map_iterator& operator++()
      { ++fInternal; return *this; }

    // postfix-increment (it++)
    map_iterator operator++(int)
      {
        map_iterator tmp(*this);
        ++fInternal;
        return tmp;
      }

    bool operator!=(map_iterator const &rhs) const
      {
        return rhs.fParent != fParent                // if we don't have same parent - different
            || fIsMap ? (fInternal != rhs.fInternal) // if array, compare internal iterator
                      : false;                       // if not array, we are equal
      }
  };

  /// Begin looping over map pairs
  map_iterator map_begin()
    { return map_iterator(*this); }

  /// Ending of map iteration
  map_iterator map_end()
    {
      if (is_map()) {
        return map_iterator(this, fValueMap.end());
      }
      return map_iterator(*this);
    }

  /// \class iterator_over_map
  /// \brief Returned by method `items_in_map` for c++11 foreach
  ///        loop syntax
  ///
  /// ```cpp
  /// for (const auto &pair : container.items_in_map()) {
  ///   ...
  /// }
  /// ```
  ///
  class iterator_over_map {
    AliFemtoConfigObject *fParent;
  public:
    iterator_over_map(AliFemtoConfigObject &obj): fParent(&obj) {}
    map_iterator begin() { return fParent->map_begin(); }
    map_iterator end() { return fParent->map_end(); }

    /// Python-interop method
    ///
    /// Use as iterator
    ///
    /// Does not check if type is map! MAY CRASH!
    ///
    MapValue_t __iter__() const
      { return fParent->is_map() ? fParent->fValueMap : MapValue_t(); }
  };

  /// Simplified loop-over-map syntax function
  iterator_over_map items_in_map()
    { return iterator_over_map(*this); }

  /// \class Popper
  /// \brief Struct used for 'poping' many values from object
  struct Popper {
    AliFemtoConfigObject* src;
    Popper(const AliFemtoConfigObject &obj);
    virtual ~Popper();

    template <typename T>
    Popper& operator()(const Key_t &key, T &dest) {
      src->pop_and_load(key, dest); return *this; }

    Popper& WarnOfRemainingItems() {
      src->WarnOfRemainingItems();
      return *this;
    }

    Popper(const Popper &orig): src(orig.src) {}
    Popper& operator=(const Popper &rhs) { src = rhs.src; return *this; }
  };

  Popper pop_all() const
    { return Popper(*this); }

  /// Return a pointer to new config object with values copied
  virtual AliFemtoConfigObject* Clone(const char *newname="") const
    { return new AliFemtoConfigObject(*this); }

  /// Return a new config object with values copied
  AliFemtoConfigObject DeepCopy() const
    { return AliFemtoConfigObject(*this); }

  /// return a string of valid JSON that may be parsed by other
  /// languages into their equivalent data types.
  ///
  /// Similar to Stringify method, this output is valid JSON.
  ///
  std::string as_JSON_string() const;

  /// Pretty-print the value
  TString Stringify(const bool pretty=false, int deep = 0) const;

  /// Consume this config object while constructing new class
  ///
  /// \param warn Will print warning if there are remaining elements
  ///   in this config object. This is the only way to be notified of
  ///   key "typos" in the config definition
  ///
  /// \return New object of type whatever
  template <typename T>
  T* Into(bool warn=true);

  /// Attempt to construct a configuration object from arbitrary source
  template <typename T>
  static AliFemtoConfigObject From(const T&);

  /// Perform a "deep" update of object into this map object.
  /// If either is not a map, do nothing.
  /// Loop through keys in "default" object, if any key is missing
  /// in `this`, copy value over. If both keys point to maps, recursively
  /// call SetDefault on the sub-objects.
  ///
  void SetDefault(const AliFemtoConfigObject &default_mapobj);
#ifdef ENABLE_MOVE_SEMANTICS
  void SetDefault(AliFemtoConfigObject &&default_mapobj);
#endif

  /// Update map with key/value pairs from other dict.
  /// If either is not a map, do nothing.
  ///
  /// \param source     Copy values from this object
  /// \param all_keys   If true, all keys are read in from source,
  ///                   otherwise only copy in the keys already present
  ///
  void Update(const AliFemtoConfigObject &source, bool all_keys=true);

  /// Print warning that multiple objects remain in this map or
  /// array. This is useful for ensuring all fields are used after
  /// "popping" all expected key-value pairs.
  void WarnOfRemainingItems(std::ostream& out=std::cerr) const;

  // ========================
  //     TObject Methods
  // ========================

  /// Maps values of ConfigObjects to an unsigned integer.
  /// If obj1 == obj2, then obj1.Hash() == obj2.Hash().
  /// No particular algorithm is used std::hash, and there is no
  /// guarantee that the same object will yield the same hash between
  /// different versions of AliPhysics.
  virtual ULong_t Hash() const;

  /// (static) title of "Configuration Object"
  virtual const char* GetTitle() const
    { return "Configuration Object"; }

  /// Called by merging actions like `hadd`
  virtual Long64_t Merge(TCollection *);

  // virtual int DistanceToPrimitive(int x, int y) { return TObject::DistanceToPrimitive(x, y); }
  virtual void Draw(Option_t* opt="");
  // virtual void ExecuteEvent(int event, int px, int py);

  /// Called when pad requests update.
  ///
  /// Calls to paint are forwarded to the painter object
  ///
  virtual void Paint(Option_t *opt=nullptr);

  /// Double clicking on object in TBrowser simply draws the object
  /// to the current pad
  virtual void Browse(TBrowser *);


protected:

  /// Current type
  TypeTagEnum_t fTypeTag;

  /// Union of all possible types
  union {
    BoolValue_t fValueBool;
    FloatValue_t fValueFloat;
    IntValue_t fValueInt;
    StringValue_t fValueString;
    ArrayValue_t fValueArray;
    MapValue_t fValueMap;
    RangeValue_t fValueRange;
    RangeListValue_t fValueRangeList;
  };

  /// Object used to do stuff
  Painter* fPainter; //!

  // template <typename OutStreamValue_t>
  // friend OutStreamValue_t& operator<<(OutStreamValue_t &, const AliFemtoConfigObject &);

  // template <typename InStreamValue_t>
  // friend InStreamValue_t& operator>>(InStreamValue_t &, AliFemtoConfigObject &);

  friend std::ostream& operator<<(std::ostream&, const AliFemtoConfigObject&);
  friend std::istream& operator>>(std::istream&, AliFemtoConfigObject&);

  friend TBuffer& operator<<(TBuffer&, const AliFemtoConfigObject&);
  friend TBuffer& operator>>(TBuffer&, AliFemtoConfigObject&);


  typedef std::string::const_iterator StringIter_t;

  // internal parsing methods - throw unhelpful errors to public one
  static AliFemtoConfigObject Parse(StringIter_t &, const StringIter_t);
  static AliFemtoConfigObject ParseMap(StringIter_t &, const StringIter_t);
  static AliFemtoConfigObject ParseArray(StringIter_t &, const StringIter_t);

  /// Parse a string 'into' destination
  // static void Parse(const std::string &, AliFemtoConfigObject &);

  static TString NameFromtype(TypeTagEnum_t tag)
  {
    switch (tag) {
      case kEMPTY: return "empty";
      case kBOOL: return "bool";
      case kINT: return "int";
      case kFLOAT: return "float";
      case kSTRING: return "str";
      case kRANGE: return "range";
      case kRANGELIST: return "rangelist";
      case kARRAY: return "array";
      case kMAP: return "map";
    }
    return "--";
  }

  static std::vector<std::string> split_key(Key_t key);
  MapValue_t::iterator find_item(Key_t key);

private:

  void _DeleteValue();
  void _CopyConstructValue(const AliFemtoConfigObject &);
#ifdef ENABLE_MOVE_SEMANTICS
  void _MoveConstructValue(AliFemtoConfigObject &&);
#endif

  /// \cond CLASSIMP
  ClassDef(AliFemtoConfigObject, 2);
  /// \endcond
};


/// \class AliFemtoConfigObject::Painter
/// \brief Responsible for painting routines
///
class AliFemtoConfigObject::Painter {

  friend class AliFemtoConfigObject;

public:
  /// Construct with parent configuration object
  Painter(AliFemtoConfigObject &data);

  /// Copy constructor
  Painter(const Painter &o)
    : fData(o.fData)
    , fTitle(o.fTitle)
    , fBody(o.fBody)
  {}

  Painter& operator=(const Painter &r)
  {
    if (this != &r) {
      fData = r.fData;
      fTitle = r.fTitle;
      fBody = r.fBody;
    }
    return *this;
  }

  void Paint();
  void ExecuteEvent(int, int, int);

protected:

  void ResetData(AliFemtoConfigObject *data)
    { fData = data; }

  AliFemtoConfigObject *fData; //!

  TText fTitle;
  TPaveText fBody;
};

inline
AliFemtoConfigObject::Popper::Popper(const AliFemtoConfigObject &obj)
: src(new AliFemtoConfigObject(obj))
{}

inline
AliFemtoConfigObject::Popper::~Popper()
{
  delete src;
}


inline
void
AliFemtoConfigObject::Paint(Option_t*)
{
  if (!fPainter) {
    fPainter = new Painter(*this);
  }

  return fPainter->Paint();
}



// /// map enums to types - useful for templates?
// #define MAP_ENUM_TO_TYPE(__enum, __type)
//   template <> struct AliFemtoConfigObject::type_from_enum<AliFemtoConfigObject::__enum > { typedef typename AliFemtoConfigObject::__type type; };
//   // template <> struct AliFemtoConfigObject::enum_from_type<AliFemtoConfigObject::__type > { static constexpr TypeTagEnum_t value { __enum }; };
//
//   MAP_ENUM_TO_TYPE(kBOOL, BoolValue_t);
//   MAP_ENUM_TO_TYPE(kINT, IntValue_t);
//   MAP_ENUM_TO_TYPE(kFLOAT, FloatValue_t);
//   MAP_ENUM_TO_TYPE(kSTRING, StringValue_t);
//   MAP_ENUM_TO_TYPE(kARRAY, ArrayValue_t);
//   MAP_ENUM_TO_TYPE(kMAP, MapValue_t);
//   MAP_ENUM_TO_TYPE(kRANGE, RangeValue_t);
//   MAP_ENUM_TO_TYPE(kRANGELIST, RangeListValue_t);
//
//  // MAP_ENUM_TO_TYPE(kSTRING, const char *);
//
// #undef MAP_ENUM_TO_TYPE

// #define MAP_ENUM_TO_TYPE(__type, __enum, _)
//   template <> struct AliFemtoConfigObject::type_from_enum<AliFemtoConfigObject::__enum > { typedef typename AliFemtoConfigObject::__type type; };
//
//   FORWARD_STANDARD_TYPES(MAP_ENUM_TO_TYPE);
//
// #undef MAP_ENUM_TO_TYPE

inline
AliFemtoConfigObject::~AliFemtoConfigObject()
{
  _DeleteValue();
  delete fPainter;
}

inline
void AliFemtoConfigObject::_DeleteValue()
{
  TypeTagEnum_t typetag = fTypeTag;
  fTypeTag = kEMPTY;
  switch (typetag) {
    case kEMPTY: break;
    case kBOOL: return fValueBool.~BoolValue_t();
    case kINT: return fValueInt.~IntValue_t();
    case kFLOAT: return fValueFloat.~FloatValue_t();
    case kSTRING: return fValueString.~StringValue_t();
    case kRANGE: return fValueRange.~RangeValue_t();
    case kRANGELIST: return fValueRangeList.~RangeListValue_t();
    case kARRAY: return fValueArray.~ArrayValue_t();
    case kMAP: return fValueMap.~MapValue_t();
  }
}

inline
void AliFemtoConfigObject::_CopyConstructValue(const AliFemtoConfigObject &src)
{
  #define COPY(__type, __dest) new (& __dest) __type (src. __dest)

  switch (src.fTypeTag) {
    case kEMPTY: break;
    case kBOOL: COPY(BoolValue_t, fValueBool); break;
    case kINT: COPY(IntValue_t, fValueInt); break;
    case kFLOAT: COPY(FloatValue_t, fValueFloat); break;
    case kSTRING: COPY(StringValue_t, fValueString); break;
    case kARRAY: COPY(ArrayValue_t, fValueArray); break;
    case kMAP: COPY(MapValue_t, fValueMap); break;
    case kRANGE: COPY(RangeValue_t, fValueRange); break;
    case kRANGELIST: COPY(RangeListValue_t, fValueRangeList); break;
  }

  #undef COPY

  fTypeTag = src.fTypeTag;
}

#ifdef ENABLE_MOVE_SEMANTICS

inline
void AliFemtoConfigObject::_MoveConstructValue(AliFemtoConfigObject &&src)
{
  #define MOVE(__type, __dest) new (& __dest) __type (std::move(src. __dest))

  switch (src.fTypeTag) {
    case kEMPTY: break;
    case kBOOL: MOVE(BoolValue_t, fValueBool); break;
    case kINT: MOVE(IntValue_t, fValueInt); break;
    case kFLOAT: MOVE(FloatValue_t, fValueFloat); break;
    case kSTRING: MOVE(StringValue_t, fValueString); break;
    case kARRAY: MOVE(ArrayValue_t, fValueArray); break;
    case kMAP: MOVE(MapValue_t, fValueMap); break;
    case kRANGE: MOVE(RangeValue_t, fValueRange); break;
    case kRANGELIST: MOVE(RangeListValue_t, fValueRangeList); break;
  }

  fTypeTag = src.fTypeTag;
  src._DeleteValue();

  #undef MOVE
}

#endif

inline
AliFemtoConfigObject&
AliFemtoConfigObject::operator=(AliFemtoConfigObject const &rhs)
{
  if (this == &rhs) {
    return *this;
  }

  TObject::operator=(rhs);

  if (rhs.fTypeTag != fTypeTag) {
    _DeleteValue();
    _CopyConstructValue(rhs);
  } else {
    #define COPY(__dest) __dest = rhs.__dest;

    switch (rhs.fTypeTag) {
      case kEMPTY: break;
      case kBOOL: COPY(fValueBool); break;
      case kINT: COPY(fValueInt); break;
      case kFLOAT: COPY(fValueFloat); break;
      case kSTRING: COPY(fValueString); break;
      case kARRAY: COPY(fValueArray); break;
      case kMAP: COPY(fValueMap); break;
      case kRANGE: COPY(fValueRange); break;
      case kRANGELIST: COPY(fValueRangeList); break;
    }

    #undef COPY

    fTypeTag = rhs.fTypeTag;
  }
  return *this;
}

/*
template <>
AliFemtoConfigObject::Popper&
AliFemtoConfigObject::Popper::operator()<int>(const Key_t &key, int &dest) {
  src->pop_and_load(key, (Int_t&)dest); return *this; }

template <>
AliFemtoConfigObject::Popper&
AliFemtoConfigObject::Popper::operator()<UInt_t>(const Key_t &key, UInt_t &dest) {
  src->pop_and_load(key, (Int_t&)dest); return *this; }

*/

#undef FORWARD_STANDARD_TYPES


#endif // ALIFEMTOCONFIGOBJECT_H
