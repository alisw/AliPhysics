///
/// \file AliFemtoUser/AliFemtoConfigObject.h
///

#pragma once

#ifndef ALIFEMTOCONFIGOBJECT_H
#define ALIFEMTOCONFIGOBJECT_H

#include <Rtypes.h>
#include <map>
#include <list>
#include <string>
#include <vector>
#include <ostream>

#include <TString.h>
#include <TText.h>

#if false && __cplusplus >= 201103L
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
  typedef Bool_t BoolValue_t;
  typedef double FloatValue_t;
  typedef Long64_t IntValue_t;
  typedef std::string StringValue_t;
  typedef std::vector<AliFemtoConfigObject> ArrayValue_t;
  typedef std::string Key_t;
  typedef std::map<Key_t, AliFemtoConfigObject> MapValue_t;
  typedef std::pair<FloatValue_t, FloatValue_t> RangeValue_t;
  typedef std::vector<RangeValue_t> RangeListValue_t;
  /// @}

  /// class for visualizing on a canvas (TPad)
  class Painter;

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

  /// Create object by parsing string
  static AliFemtoConfigObject Parse(const std::string &);

  /// Create object in empty state
  AliFemtoConfigObject();

  /// Construct from TObject
  AliFemtoConfigObject(const TObject *);
  AliFemtoConfigObject(const TObject &);


  /// Calls destructor on internal data
  virtual ~AliFemtoConfigObject();

  /// Copy value
  AliFemtoConfigObject(AliFemtoConfigObject const &);

  /// Assign value
  AliFemtoConfigObject& operator=(AliFemtoConfigObject const &);

  // define constructors
  #define ConstructorDef(__type, __tag, __dest) \
    AliFemtoConfigObject(const __type &v): fTypeTag(__tag), __dest(v), fPainter(nullptr) { }


    // constructors with all associated value-types (IntValue_t, MapValue_t, etc..)
    FORWARD_STANDARD_TYPES(ConstructorDef);

    ConstructorDef(TString, kSTRING, fValueString);
    ConstructorDef(char *, kSTRING, fValueString);
    ConstructorDef(int, kINT, fValueInt);
    ConstructorDef(float, kFLOAT, fValueFloat);
  #undef ConstructorDef

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
  AliFemtoConfigObject(AliFemtoConfigObject &&);

  /// Move assignment
  AliFemtoConfigObject& operator=(AliFemtoConfigObject &&rhs);

  #define ConstructorDef(__type, __tag, __dest) AliFemtoConfigObject(__type &&v): fTypeTag(__tag), __dest(std::move(v)), fPainter(nullptr) { }
    FORWARD_STANDARD_TYPES(ConstructorDef);
  #undef ConstructorDef

#endif // move-semantics


  // ========================
  //      C++ Operators
  // ========================
  bool operator==(const AliFemtoConfigObject &rhs) const;
  bool operator!=(const AliFemtoConfigObject &rhs) const { return !(*this == rhs); }
  #define IMPL_EQUALITY(__type, __tag, __target) \
    bool operator==(const __type &rhs) const { return fTypeTag == __tag && __target == rhs; } \
    bool operator!=(const __type &rhs) const { return !operator==(rhs); }

    FORWARD_STANDARD_TYPES(IMPL_EQUALITY);
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

  bool load_bool(BoolValue_t &b) const { return is_bool() ? b = fValueBool, true : false; }
  bool load_float(FloatValue_t &f) const { return is_float() ? f = fValueFloat, true : false; }
  bool load_int(IntValue_t &i) const { return is_int() ? i = fValueInt, true : false; }
  bool load_num(FloatValue_t &f) const { return is_int() ? f = fValueInt, true : load_float(f); }
  bool load_str(StringValue_t &s) const { return is_str() ? s = fValueString, true : false; }
  bool load_array(ArrayValue_t &a) const { return is_array() ? a = fValueArray, true : false; }
  bool load_map(MapValue_t &m) const { return is_map() ? m = fValueMap, true : false; }
  bool load_range(RangeValue_t &r) const { return is_range() ? r = fValueRange, true : false; }
  bool load_range(Float_t &a, Float_t &b) const {
    return is_range() ? a = fValueRange.first, b = fValueRange.second, true : false; }
  bool load_rangelist(RangeListValue_t &r) const {
    return is_rangelist() ? r = fValueRangeList, true : load_rangelist(r); }
  bool load_ranges(RangeListValue_t &r) const {
    return is_range() ? r = RangeListValue_t({fValueRange}), true : load_rangelist(r); }


  /// Pretty-print the value
  TString Stringify(const bool pretty=true) const;


  // ========================
  //     TObject Methods
  // ========================

  virtual ULong_t Hash() const;
  virtual const char* GetTitle() const { return "Configuration Object"; }

  /// Called by merging actions like `hadd`
  virtual Long64_t Merge(TCollection *);

  // virtual int DistanceToPrimitive(int x, int y) { return TObject::DistanceToPrimitive(x, y); }
  // virtual void Draw(Option_t* ="");
  // virtual void ExecuteEvent(int event, int px, int py);

  /// Called when pad requests update.
  ///
  /// Calls to paint are forwarded to the painter object
  ///
  virtual void Paint(Option_t *opt=nullptr);

  /// Double clicking on object in TBrowser simply draws the object
  /// to the current pad
  virtual void Browse(TBrowser *) { Draw(); }


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

  template <typename OutStreamValue_t>
  friend OutStreamValue_t& operator<<(OutStreamValue_t &, const AliFemtoConfigObject &);

  template <typename InStreamValue_t>
  friend InStreamValue_t& operator>>(InStreamValue_t &, AliFemtoConfigObject &);

  typedef std::string::const_iterator StringIter_t;

  // internal parsing methods - throw unhelpful errors to public one
  static AliFemtoConfigObject Parse(StringIter_t &, const StringIter_t);
  static AliFemtoConfigObject ParseMap(StringIter_t &, const StringIter_t);
  static AliFemtoConfigObject ParseArray(StringIter_t &, const StringIter_t);

  /// Parse a string 'into' destination
  // static void Parse(const std::string &, AliFemtoConfigObject &);


private:

  void _DeleteValue();
  void _CopyConstructValue(const AliFemtoConfigObject &);
#ifdef ENABLE_MOVE_SEMANTICS
  void _MoveConstructValue(const AliFemtoConfigObject &&);
#endif

  ClassDef(AliFemtoConfigObject, 1);
};


/// \class AliFemtoConfigObject::Painter
/// \brief Responsible for painting routines
///
class AliFemtoConfigObject::Painter {
protected:
  AliFemtoConfigObject *fData; //!

public:
  /// Construct with a
  Painter(AliFemtoConfigObject &data);

  void Paint();
  void ExecuteEvent(int, int, int);

protected:
  TText fTitle;
};

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
  switch (fTypeTag) {
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
}

inline
AliFemtoConfigObject::AliFemtoConfigObject():
  TObject()
  , fTypeTag(kEMPTY)
  , fPainter(NULL)
{
}

inline
AliFemtoConfigObject::AliFemtoConfigObject(AliFemtoConfigObject const &orig):
  TObject(orig)
  , fTypeTag(orig.fTypeTag)
  , fPainter(NULL)
{
  _CopyConstructValue(orig);
}

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
  }
  fTypeTag = rhs.fTypeTag;
  return *this;
}


#ifdef ENABLE_MOVE_SEMANTICS

inline
AliFemtoConfigObject::AliFemtoConfigObject(AliFemtoConfigObject &&orig):
  TObject(orig)
  , fTypeTag(orig.fTypeTag)
  , fPainter(orig.fPainter)
{
  _MoveConstructValue(std::move(orig));
  orig._DeleteValue();
  orig.fTypeTag = kEMPTY;
  orig.fPainter = nullptr;
}

/// Move assignment-operator
inline
AliFemtoConfigObject& AliFemtoConfigObject::operator=(AliFemtoConfigObject &&rhs)
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
  }

  fTypeTag = rhs.fTypeTag;

  rhs._DeleteValue();
  rhs.fTypeTag = kEMPTY;

  delete fPainter;
  fPainter = rhs.fPainter;
  rhs.fPainter = nullptr;

  return *this;
}

#endif // move-semantics


#undef FORWARD_STANDARD_TYPES





#endif // ALIFEMTOCONFIGOBJECT_H
