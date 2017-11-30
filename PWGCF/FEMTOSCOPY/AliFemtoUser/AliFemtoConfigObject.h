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
#include <iostream>

#include <TBuffer.h>
#include <TString.h>
#include <TText.h>
#include <TPaveText.h>


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

  /// Create object from string; if a map, also parse defaults and
  /// insert any values missing in object with defaults.
  static AliFemtoConfigObject ParseWithDefaults(const std::string &, const std::string &defaults);

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

    //ConstructorDef(TString, kSTRING, fValueString);   // implicit cast - not allowed by clang 9.0, must be declared explicitly
    ConstructorDef(char *, kSTRING, fValueString);
    ConstructorDef(int, kINT, fValueInt);
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
  AliFemtoConfigObject(AliFemtoConfigObject &&);

  /// Move assignment
  AliFemtoConfigObject& operator=(AliFemtoConfigObject &&rhs);

  #define ConstructorDef(__type, __tag, __dest) AliFemtoConfigObject(__type &&v): fTypeTag(__tag), __dest(std::move(v)), fPainter(nullptr) { }
    FORWARD_STANDARD_TYPES(ConstructorDef);
  #undef ConstructorDef

#endif // move-semantics

  typedef std::pair<float,float> pair_of_floats;

  /// \class BuildStruct
  /// \brief helper class for building a mapping object
  class BuildMap {
  public:
    BuildMap() {}

#ifdef ENABLE_MOVE_SEMANTICS
    #define IMPL_BUILDITEM(__type, __a, __b) \
      BuildMap& operator()(const Key_t &key, const __type& val) { fMap[key] = val; return *this; }  \
      BuildMap& operator()(const Key_t &key, __type && val) { fMap.insert(std::make_pair(key, std::move(val))); return *this; }
    #define IMPL_CASTED_BUILDITEM(__type, __savedtype) \
      BuildMap& operator()(const Key_t &key, const __type& val) { fMap[key] = static_cast<__savedtype>(val); return *this; }  \
      BuildMap& operator()(const Key_t &key, __type && val) { fMap.insert(std::make_pair(key, std::move(static_cast<__savedtype>(val)))); return *this; }
#else
    #define IMPL_BUILDITEM(__type, __a, __b) \
      BuildMap& operator()(const Key_t &key, const __type& val) { fMap.insert(std::make_pair(key, val)); return *this; }
    #define IMPL_CASTED_BUILDITEM(__type, __savedtype) \
      BuildMap& operator()(const Key_t &key, const __type& val) { fMap.insert(std::make_pair(key, static_cast<__savedtype>(val))); return *this; }
#endif

    FORWARD_STANDARD_TYPES(IMPL_BUILDITEM);

    IMPL_CASTED_BUILDITEM(Int_t, IntValue_t);
    IMPL_CASTED_BUILDITEM(UInt_t, IntValue_t);
    IMPL_CASTED_BUILDITEM(pair_of_floats, RangeValue_t);

    IMPL_BUILDITEM(AliFemtoConfigObject, 0, 0);
    #undef IMPL_BUILDITEM

    operator AliFemtoConfigObject() {
      return AliFemtoConfigObject(std::move(fMap));
    }

  protected:
    MapValue_t fMap;
  };


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
    if (is_range()) {
      r.clear();
      r.push_back(fValueRange);
      return true;
    } else {
      return load_rangelist(r);
    }
    // return is_range() ? r.clear(), r.push_back(fValueRange), true : load_rangelist(r);
  }

  bool load_range(std::pair<float, float> &r) const { return is_range() ? r = fValueRange, true : false; }



  /// \defgroup Find&Load Methods
  /// @{
  /// Copies item identified by *key*, returns true if found
  #define IMPL_FINDANDLOAD(__dest_type, __tag, __source)            \
    bool find_and_load(const Key_t &key, __dest_type &dest) const { \
      if (!is_map()) { return false; }                        \
      auto found = fValueMap.find(key);                       \
      if (found == fValueMap.end()) { return false; }         \
      if (found->second.fTypeTag != __tag) { return false; }  \
      dest = found->second. __source;                         \
      return true; }

    FORWARD_STANDARD_TYPES(IMPL_FINDANDLOAD)

    IMPL_FINDANDLOAD(pair_of_floats, kRANGE, fValueRange);
    IMPL_FINDANDLOAD(int, kINT, fValueInt);
    IMPL_FINDANDLOAD(unsigned int, kINT, fValueInt);

  #undef IMPL_FINDANDLOAD
  /// @}


  /// Remove and returns pointer to object at *index* if an array
  ///
  /// Delete the pointer after use!
  ///
  /// Index starts at 0. Negative numbers wrap around to the back, so '-1'
  /// removes the last element (default behavior).
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
  AliFemtoConfigObject* pop(const Key_t &key, const ReturnType &default_);

  /// Return copy of object with a key removed from map.
  /// If not a map, or key not present, object is simply copied
  AliFemtoConfigObject WithoutKey(const Key_t &key) const;

  /// Return copy of object with given keys removed
  AliFemtoConfigObject WithoutKeys(const std::vector<Key_t>&) const;

  /// \defgroup Pop&Load Methods
  /// @{
  /// Removes item identified by *key*

  #define IMPL_POPANDLOAD(__dest_type, __tag, __source)      \
    bool pop_and_load(const Key_t &key, __dest_type &dest) { \
      if (!is_map()) { return false; }                       \
      auto found = fValueMap.find(key);                      \
      if (found == fValueMap.end()) { return false; }        \
      if (found->second.fTypeTag != __tag) { return false; } \
      dest = found->second. __source;                        \
      fValueMap.erase(found); return true; }

    FORWARD_STANDARD_TYPES(IMPL_POPANDLOAD)

    IMPL_POPANDLOAD(pair_of_floats, kRANGE, fValueRange);
    IMPL_POPANDLOAD(int, kINT, fValueInt);
    IMPL_POPANDLOAD(unsigned int, kINT, fValueInt);

  #undef IMPL_POPANDLOAD

  /// @}

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
  };

  Popper pop_all() const { return Popper(*this); }

  /// return a string of valid JSON that may be parsed by other
  /// languages into their equivalent data types.
  ///
  /// Similar to Stringify method, this output is valid JSON.
  ///
  std::string as_JSON_string() const;

  /// A general template method for building objects with config object
  ///
  /// Default implementation simply calls the type's constructor with
  /// a const reference to this object.
  template <typename T>
  T* Construct() const;

  /// Consume this config object while constructing new class
  ///
  /// \param warn Will print warning if there are remaining elements
  ///   in this config object. This is the only way to be notified of
  ///   key "typos" in the config definition
  ///
  /// \return New object of type whatever
  template <typename T>
  T* Into(bool warn=true);

  /// Pretty-print the value
  TString Stringify(const bool pretty=false, int deep = 0) const;

  /// Perform a "deep" update of object into this map object.
  /// If either is not a map, do nothing.
  /// Loop through keys in "default" object, if any key is missing
  /// in `this`, copy value over. If both keys point to maps, recursively
  /// call SetDefault on thses sub-objects.
  ///
  void SetDefault(const AliFemtoConfigObject &default_mapobj);
#ifdef ENABLE_MOVE_SEMANTICS
  void SetDefault(AliFemtoConfigObject &&default_mapobj);
#endif

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


template <typename StreamType>
StreamType& operator<<(StreamType &stream, const AliFemtoConfigObject *ptr)
{
  return std::operator<<(stream, *ptr);
}

template <typename StreamType>
StreamType& operator>>(StreamType &stream, AliFemtoConfigObject *&ptr)
{
  if (!ptr) {
    ptr = new AliFemtoConfigObject();
  }
  return std::operator>>(stream, *ptr);
}

template <typename ReturnType>
AliFemtoConfigObject* AliFemtoConfigObject::pop(const Key_t &key, const ReturnType &default_)
{
  if (auto *obj = pop(key)) {
    return obj;
  }
  return new AliFemtoConfigObject(default_);
}

template <typename T>
T* AliFemtoConfigObject::Construct() const
{
  return new T(*this);
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
