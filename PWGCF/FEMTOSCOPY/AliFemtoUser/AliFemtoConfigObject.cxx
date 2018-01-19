///
/// \file AliFemtoUser/AliFemtoConfigObject.cxx
///

#include "AliFemtoConfigObject.h"

#include <TObjString.h>
#include <TObjArray.h>
#include <TCollection.h>

#include <regex>
#include <cctype>
#include <sstream>
#include <iostream>
#include <exception>
#include <functional>
#include <numeric>

/// \cond CLASSIMP
ClassImp(AliFemtoConfigObject);
/// \endcond


AliFemtoConfigObject::AliFemtoConfigObject(const TObject *obj)
: AliFemtoConfigObject(*obj)
{
  // std::cout << "[AliFemtoConfigObject::AliFemtoConfigObject(const TObject*)] " << obj << "\n";
}

AliFemtoConfigObject::AliFemtoConfigObject(const TObject &obj)
: TObject(obj)
, fPainter(nullptr)
{

  if (const auto *config = dynamic_cast<const AliFemtoConfigObject*>(&obj)) {
    *this = *config;
  }
  else if (const auto *str_ptr = dynamic_cast<const TObjString*>(&obj)) {
    fTypeTag = kSTRING;
    new (&fValueString) StringValue_t (str_ptr->GetString().Data());
  }
  else {
    throw std::runtime_error(std::string("Unexpected TObject initalization with point to ") + obj.ClassName());
  }
}

bool
AliFemtoConfigObject::operator==(const AliFemtoConfigObject &rhs) const
{
  if (this == &rhs) {
    return true;
  }

  if (fTypeTag != rhs.fTypeTag) {
    return false;
  }

  switch (fTypeTag) {
    case kEMPTY: return true;
    case kBOOL: return fValueBool == rhs.fValueBool;
    case kINT: return fValueInt == rhs.fValueInt;
    case kFLOAT: return fValueFloat == rhs.fValueFloat;
    case kSTRING: return fValueString == rhs.fValueString;
    case kARRAY: return fValueArray == rhs.fValueArray;
    case kMAP: return fValueMap == rhs.fValueMap;
    case kRANGE: return fValueRange == rhs.fValueRange;
    case kRANGELIST: return fValueRangeList == rhs.fValueRangeList;
  }

  return false;
}

// template<>
// struct AliFemtoConfigObject::enum_from_type< bool >
// {
//    static const TypeTagEnum_t value;
//   };

//  const AliFemtoConfigObject::TypeTagEnum_t AliFemtoConfigObject::enum_from_type< AliFemtoConfigObject::BoolValue_t >::value = static_cast<AliFemtoConfigObject::TypeTagEnum_t>(1);

TString
AliFemtoConfigObject::Stringify(bool pretty, int deep) const
{
#define INDENTSTEP 4
  switch (fTypeTag) {
    case kEMPTY: return "";
    case kBOOL: return fValueBool ? "true" : "false";
    case kINT: return TString::Format("%lld", fValueInt);
    case kFLOAT: return TString::Format("%f", fValueFloat);
    case kSTRING: return TString::Format("'%s'", fValueString.c_str());
    case kRANGE: return TString::Format("%f:%f", fValueRange.first, fValueRange.second);
    case kARRAY: {
      TString result = "[";
      auto it = fValueArray.cbegin(),
          end = fValueArray.cend();
      if (it != end) {
        result += it->Stringify(pretty);
      }
      for (++it; it != end; ++it) {
        result += ", " + it->Stringify(pretty);
      }
      result += ']';
      return result;
    }
    case kRANGELIST: {
      TString result = '(';
      auto it = fValueRangeList.cbegin(),
          end = fValueRangeList.cend();
      if (it != end) {
        result += TString::Format("%f:%f", it->first, it->second);
      }
      for (++it; it != end; ++it) {
        if (it->first == std::prev(it)->second) {
          result += TString::Format(":%f", it->second);
        } else {
          result += TString::Format(", %f:%f", it->first, it->second);
        }
      }
      result += ')';
      return result;
    }
    case kMAP: {
     TString result = '{';
     if (pretty) {
       result += TString('\n');
     }
     auto it = fValueMap.cbegin(),
         end = fValueMap.cend();
      if (it != end) {
        if (pretty) {
          result += TString(' ',(deep+1)*INDENTSTEP);
          result += it->first + ": " + it->second.Stringify(pretty,deep+1);
          result += TString('\n');
        }
        else {
          result += it->first + ": " + it->second.Stringify(pretty);
        }
      }
      for (++it; it != end; ++it) {
        if (pretty) {
          result += TString(' ',(deep+1)*INDENTSTEP);
          result += TString::Format("%s: %s", it->first.c_str(), it->second.Stringify(pretty, deep+1).Data());
          result += TString('\n');
        }
        else {
          result += TString::Format(", %s: %s", it->first.c_str(), it->second.Stringify(pretty).Data());
        }
      }
      if (pretty) {
        result += TString(' ',deep*INDENTSTEP);
      }
      result += '}';
      return result;
    }
  }
  return "";
}


AliFemtoConfigObject*
AliFemtoConfigObject::pop(Int_t idx)
{
  AliFemtoConfigObject *result = nullptr;
  if (is_array()) {
    auto N = static_cast<Int_t>(fValueArray.size());

    // out of bounds
    if (idx <= -N || N <= idx) {
      return nullptr;
    }

    auto it = fValueArray.begin();

    // wrap negative numbers around
    if (idx < 0) {
      idx += N;
    }

    // move iterator to index, move value and erase from array
    std::advance(it, idx);
    result = new AliFemtoConfigObject(std::move(*it));
    fValueArray.erase(it);
  }
  return result;
}


AliFemtoConfigObject*
AliFemtoConfigObject::pop(const Key_t &key)
{
  AliFemtoConfigObject *result = nullptr;

  if (is_map()) {
    auto found = fValueMap.find(key);
    if (found != fValueMap.end()) {
      result = new AliFemtoConfigObject(std::move(found->second));
      fValueMap.erase(found);
    }
  }
  return result;
}


void
AliFemtoConfigObject::WarnOfRemainingItems(std::ostream& out) const
{
  if (is_map() && fValueMap.size() > 0) {
    out << "Warning - AliFemtoConfigObject map has unexpected leftover objects: ";
    for (auto &it : fValueMap) {
      out << it.first << "[" << NameFromtype(it.second.fTypeTag) << "] ";
    }
    out << "\n";
  }
  else if(is_array() && fValueArray.size() > 0) {
    out << "Warning - AliFemtoConfigObject array has " << fValueArray.size() << " unexpected leftover objects: ";
    for (auto &val : fValueArray) {
      out << "[" << NameFromtype(val.fTypeTag) << "] ";
    }
    out << "\n";
  }
}


//=================================
//
//   TObject Methods
//

Long64_t
AliFemtoConfigObject::Merge(TCollection *collection)
{
  if (collection == nullptr) {
    Warning("AliFemtoConfigObjet::Merge", "Asked to merge with null collection");
    return 0;
  }

  TIter next_obj(collection);

  while (TObject *obj = next_obj()) {
    auto cfg = dynamic_cast<AliFemtoConfigObject*>(obj);
    // already in collection
    if (*cfg == *this) {
      return 0;
    }
  }

  collection->Add(this);
  return 1;
}

// template<>
TBuffer& operator<<(TBuffer &stream, const AliFemtoConfigObject &cfg)
{
  stream << cfg.fTypeTag;

  switch (cfg.fTypeTag) {
    case AliFemtoConfigObject::kEMPTY:
      break;
    case AliFemtoConfigObject::kBOOL:
      stream << cfg.fValueBool;
      break;
    case AliFemtoConfigObject::kFLOAT:
      // std::cout << "<< kFLOAT " << cfg.fValueFloat << "\n";
      stream << cfg.fValueFloat;
      break;
    case AliFemtoConfigObject::kINT:
      stream << cfg.fValueInt;
      break;
    case AliFemtoConfigObject::kRANGE:
      stream << std::get<0>(cfg.fValueRange);
      stream << std::get<1>(cfg.fValueRange);
      break;
    case AliFemtoConfigObject::kRANGELIST:
      stream << static_cast<std::size_t>(cfg.fValueRangeList.size());
      for (const auto &range : cfg.fValueRangeList) {
        stream << std::get<0>(range);
        stream << std::get<1>(range);
      }
      break;

    case AliFemtoConfigObject::kARRAY:
      stream << static_cast<std::size_t>(cfg.fValueArray.size());
      for (const auto &val : cfg.fValueArray) {
        stream << val;
      }
      break;

    case AliFemtoConfigObject::kMAP:
      stream << static_cast<std::size_t>(cfg.fValueMap.size());
      for (const auto &val : cfg.fValueMap) {
        stream << TString(std::get<0>(val));
        stream << std::get<1>(val);
      }
      break;

    case AliFemtoConfigObject::kSTRING:
      stream << TString(cfg.fValueString);
      break;
  }
  return stream;
}

// template<>
TBuffer& operator>>(TBuffer &stream, AliFemtoConfigObject &cfg)
{
  using ENUM_TYPE = decltype(AliFemtoConfigObject::kEMPTY);
  UInt_t type = 0;

  cfg._DeleteValue();

  //stream.ReadUInt(type);
  stream >> type;
  cfg.fTypeTag = static_cast<ENUM_TYPE>(type);

  std::size_t count = 1;

  switch (cfg.fTypeTag) {
    case AliFemtoConfigObject::kEMPTY:
      break;

    case AliFemtoConfigObject::kBOOL:
      new (&cfg.fValueBool) AliFemtoConfigObject::BoolValue_t ();
      stream >> cfg.fValueBool;
      break;

    case AliFemtoConfigObject::kINT:
      new (&cfg.fValueInt) AliFemtoConfigObject::IntValue_t ();
      stream >> cfg.fValueInt;
      break;

    case AliFemtoConfigObject::kFLOAT:
      new (&cfg.fValueFloat) AliFemtoConfigObject::FloatValue_t ();
      stream >> cfg.fValueFloat;
      break;

    case AliFemtoConfigObject::kRANGE:
      new (&cfg.fValueRange) AliFemtoConfigObject::RangeValue_t ();
      stream >> std::get<0>(cfg.fValueRange);
      stream >> std::get<1>(cfg.fValueRange);
      break;

    case AliFemtoConfigObject::kRANGELIST:
      stream >> count;
      new (&cfg.fValueRangeList) AliFemtoConfigObject::RangeListValue_t (count);
      for (auto &range : cfg.fValueRangeList) {
        stream >> std::get<0>(range);
        stream >> std::get<1>(range);
      }
      break;

    case AliFemtoConfigObject::kSTRING:
      {
        TString stringbuff;
        stream >> stringbuff;
        new (&cfg.fValueString) AliFemtoConfigObject::StringValue_t (stringbuff.Data());
      }
      break;

    case AliFemtoConfigObject::kARRAY:
      // new (&cfg.fValueArray) AliFemtoConfigObject::Array_t (count);
      // for (auto &val : cfg.fValueArray) {
      //   stream >> val;
      // }
      new (&cfg.fValueArray) AliFemtoConfigObject::ArrayValue_t ();
      cfg.fValueArray.reserve(count);
      for (std::size_t i = 0; i < count; ++i) {
        AliFemtoConfigObject tmp;
        stream >> tmp;
        cfg.fValueArray.emplace_back(std::move(tmp));
      }
      break;

      case AliFemtoConfigObject::kMAP:
      stream >> count;
      new (&cfg.fValueMap) AliFemtoConfigObject::MapValue_t ();
      for (std::size_t i = 0; i < count; ++i) {
        TString keybuff;
        AliFemtoConfigObject val;
        stream >> keybuff;
        stream >> val;
        cfg.fValueMap.emplace(keybuff.Data(), std::move(val));
      }
      break;
  }

  return stream;
}

void
AliFemtoConfigObject::Streamer(TBuffer &buff)
{
  AliFemtoConfigObject &value = *this;

  if (buff.IsReading()) {
    Version_t v = buff.ReadVersion();
    if (v != 1) {
      std::cerr << "W-AliFemtoConfigObject: Unknown AliFemtoConfigObject version " << v << "\n";
    }
    TObject::Streamer(buff);
    // std::cout << "[AliFemtoConfigObject::Streamer] Reading into AliFemtoConfigObject value at " << this << "\n";
    buff >> value;
    // std::cout << "    is-empty? " << value.is_empty() << " " << " (this->is_empty(): " << this->is_empty() << ")\n";
  } else {
    buff.WriteVersion(AliFemtoConfigObject::IsA());
    TObject::Streamer(buff);
    buff << value;
  }

}


//===============================================
//
//    AliFemtoConfigObject::Painter Methods
//


AliFemtoConfigObject::Painter::Painter(AliFemtoConfigObject &data):
fData(&data)
, fTitle(0.05, 0.9, "[-] AliFemtoConfigObject")
, fBody(0.05, 0.85, 0.9, 0.1, "NB")
{
  // std::cout << "[AliFemtoConfigObjectPainter::AliFemtoConfigObjectPainter]\n";
  fBody.SetTextSize(18);
}

void
AliFemtoConfigObject::Painter::ExecuteEvent(int e, int x, int y)
{
  printf("[AliFemtoConfigObject::Painter::ExecuteEvent] %d (%d %d)\n", e, x, y);
}

void
AliFemtoConfigObject::Painter::Paint()
{
  auto t = &fTitle;
  t->SetTextAlign(13);
  t->SetTextColor(kRed+2);
  t->SetTextFont(43);
  t->SetTextSize(25);
  // t->SetTextAngle(45);
  t->Draw();


  TString result(fData->Stringify(true));
/*
  if (fData->is_map()) {
    for (auto &pair : fData->fValueMap) {
      result += pair.first + ": ";
      if (pair.second.is_map()) {
        aresult += "{ [+] }\n";
      } if (pair.second.is_float()) {
        pair += TString::Format("%f", )
      }

    }

  }
*/

  fBody.Clear();
  fBody.SetTextAlign(13);
  fBody.SetTextColor(kBlue+2);
  fBody.SetTextFont(43);
  fBody.SetTextSize(18);
  fBody.SetFillColor(0);
  TObjArray *lines = result.Tokenize('\n');
  for (Int_t iline = 0; iline < lines->GetEntriesFast(); iline++)
    fBody.AddText(((TObjString*) lines->At(iline))->String().Data());
  fBody.Draw();
  delete lines;
}




//========================
//
//    Parsing Methods!
//

using StringIter_t = std::string::const_iterator;


struct ParseError : public std::runtime_error {

  StringIter_t it;
  std::string type_being_parsed;

  ParseError(StringIter_t i, const std::string &type, const std::string &msg)
  : std::runtime_error(msg)
  , it(i)
  , type_being_parsed(type)
  { }

};

#define DEFINE_PARSING_ERROR(__name) \
  struct __name : ParseError {       \
    __name(StringIter_t it): ParseError(it, "", "") { } \
    __name(StringIter_t it, const std::string &type, const std::string &str): ParseError(it, type, str) { } };

  DEFINE_PARSING_ERROR(EmptySourceString);
  DEFINE_PARSING_ERROR(UnexpectedCharacter);
  DEFINE_PARSING_ERROR(UnexpectedLeadingChar);
  DEFINE_PARSING_ERROR(UnexpectedEndOfInput);

#undef DEFINE_PARSING_ERROR

AliFemtoConfigObject
AliFemtoConfigObject::Parse(const char *src)
{
  std::string s(src);
  return AliFemtoConfigObject::Parse(s);
}

AliFemtoConfigObject
AliFemtoConfigObject::Parse(const TString &src)
{
  std::string s(src.Data());
  return AliFemtoConfigObject::Parse(s);
}

/// Create object by parsing TString
AliFemtoConfigObject
AliFemtoConfigObject::Parse(const TObjString &src)
{
  std::string s(src.GetString().Data());
  return AliFemtoConfigObject::Parse(s);
}



// public parsing method - This should handle all user-facing parsing errors
AliFemtoConfigObject
AliFemtoConfigObject::Parse(const std::string &src)
{
  StringIter_t begin = src.begin(),
                       it = begin;

  try {
    return Parse(it, src.cend());
  }

  catch (EmptySourceString &) {
    std::cerr << "Warning - Parsed empty source string\n";
    return AliFemtoConfigObject();
  }

  catch (UnexpectedLeadingChar &err) {
    std::stringstream ss;
    ss << "Unexpected leading character '" << *it << "' when parsing type:" << err.type_being_parsed << ".";
    throw std::runtime_error(ss.str());
  }

  catch (std::runtime_error &err) {
    std::cerr << "Error!" << err.what() << "\n";
    throw;
  }

}

AliFemtoConfigObject
AliFemtoConfigObject::ParseWithDefaults(const std::string &src, const std::string &defaults)
{
  AliFemtoConfigObject obj = Parse(src);
  if (obj.is_map()) {
    auto d = Parse(defaults);
    obj.SetDefault(d);
  }
  return obj;
}

void
AliFemtoConfigObject::SetDefault(const AliFemtoConfigObject &d)
{
  // only proceed if both are maps
  if (!is_map() || !d.is_map()) {
    return;
  }

  for (auto &kv_pair : d.fValueMap) {
    auto found = fValueMap.find(kv_pair.first);
    if (found == fValueMap.end()) {
      fValueMap.emplace(kv_pair.first, kv_pair.second);
    }
    else if (found->second.is_map() && kv_pair.second.is_map()) {
      found->second.SetDefault(kv_pair.second);
    }
  }
}

#define INT_PATTERN "\\-?\\d+"
#define FLT_PATTERN "\\-?(?:inf|nan|(?:\\d+\\.\\d*|\\.\\d+)(?:e[+\\-]?\\d+)?|\\d+e[+\\-]?\\d+)"
#define NUM_PATTERN "(?:" FLT_PATTERN "|" INT_PATTERN ")"
#define STR_PATTERN "'[^']*'"
#define RANGE_PATTERN "(" NUM_PATTERN "):(" NUM_PATTERN ")"
#define BOOL_PATTERN "(?:true|false)"
// numbers separated by colons
#define MULTIRANGE_PATTERN NUM_PATTERN "(?:\\s*\\:\\s*" NUM_PATTERN  ")+"
// group [1] is the contents between parens
#define RANGELIST_PATTERN   \
  "\\(\\s*"                 \
    "("                     \
    MULTIRANGE_PATTERN      \
      "(?:\\s*,\\s*"        \
        MULTIRANGE_PATTERN  \
      ")*"                  \
    ")"                     \
  ",?\\s*\\)\\s*"

    // NAME_PATTERN = "[A-Za-z_]+[A-Za-z_0-9]*",

#define IDENT_PATTERN "[a-zA-Z_]+[a-zA-Z0-9_]*"
#define KEY_PATTERN IDENT_PATTERN "(?:\\." IDENT_PATTERN "|/" IDENT_PATTERN ")*"

    // IDENT_PATTERN = "[a-zA-Z_]+[a-zA-Z0-9_]*(?:\\.[a-zA-Z_]+[a-zA-Z0-9_]*)*",

#define LBRK_PAT "^\\s*" "\\{" "\\s*"
#define RBRK_PAT  "\\s*" "\\}" "\\s*"
#define LPRN_PAT        "\\(" "\\s*"
#define RPRN_PAT         "\\)" "\\s*"

#define LSQRBRK_PAT      "\\[" "\\s*"
#define RSQRBRK_PAT      "\\]" "\\s*"

#define COMMA_PAT "\\s*" "," "\\s*"
#define COLON_PAT "\\s*" ":" "\\s*"

#if !defined(__clang__) && __GNUC__ == 4 && __GNUC_MINOR__ < 9
static const std::regex
    NUM_RX(""),
    INT_RX(""),
    FLT_RX(""),
    BOOL_RX(""),
    STR_RX(""),
    ID_RX(""),
    KEY_RX(""),
    RANGE_RX(""),
    RANGELIST_RX(""),
    LBRK_RX(""),
    RBRK_RX(""),
    LSQRBRK_RX(""),
    RSQRBRK_RX(""),
    LPAREN_RX(""),
    RPAREN_RX(""),
    COLON_RX(""),
    COMMA_RX("");

#else
static const std::regex
    NUM_RX("^" NUM_PATTERN, std::regex_constants::icase),
    INT_RX("^" INT_PATTERN),
    FLT_RX("^(" FLT_PATTERN ")", std::regex_constants::icase),
    BOOL_RX("^" BOOL_PATTERN),
    STR_RX("^" STR_PATTERN),
    ID_RX("^" IDENT_PATTERN),
    KEY_RX("^" KEY_PATTERN),
    RANGE_RX("^(?:" RANGE_PATTERN ")"),
    RANGELIST_RX("^" RANGELIST_PATTERN),
    LBRK_RX("^" LBRK_PAT),
    RBRK_RX("^" RBRK_PAT),
    LSQRBRK_RX("^" LSQRBRK_PAT),
    RSQRBRK_RX("^" RSQRBRK_PAT),
    LPAREN_RX("^" LPRN_PAT),
    RPAREN_RX("^" RPRN_PAT),
    COLON_RX("^" COLON_PAT),
    COMMA_RX("^" COMMA_PAT);
#endif

std::vector<std::string>
AliFemtoConfigObject::split_key(Key_t key)
{
  std::smatch match;
  std::vector<std::string> result;
  if (std::regex_match(key.cbegin(), key.cend(), match, KEY_RX)) {
    for (auto ptr = match[0].first, stop = match[0].second;
         std::regex_search(ptr, stop, match, ID_RX);
         ptr = match[0].second + 1) {
      result.emplace_back(match[0].first, match[0].second);
      if (match[0].second == key.cend()) {
        break;
      }
    }
  }
  return result;
}


std::vector<AliFemtoConfigObject::Key_t>
parse_map_key_unchecked(const StringIter_t& it, const StringIter_t stop)
{
  std::smatch match;
  std::vector<std::string> result;
  StringIter_t ptr = it;
  while (std::regex_search(ptr, stop, match, ID_RX)) {
    result.emplace_back(match[0].first, match[0].second);
    ptr = match[0].second + 1;

    if (match[0].second == stop) {
      break;
    }
  }

  return result;
}


AliFemtoConfigObject
AliFemtoConfigObject::ParseArray(StringIter_t& it, const StringIter_t stop)
{
  StringIter_t START = it;

  ArrayValue_t result;

  std::smatch match;

  // ensure we are starting with bracket
  if (std::regex_search(it, stop, match, LSQRBRK_RX)) {
    it = match[0].second;
  } else {
    // UnexpectedEndOfInput
    throw UnexpectedLeadingChar(it, "Array", "square-bracket '['");
    // throw std::runtime_error(std::string("Expected to match square-bracket, found '") + *it + "'");
  }

  // loop until ']'
  while (!std::regex_search(it, stop, match, RSQRBRK_RX)) {
    auto val = Parse(it, stop);
    result.emplace_back(val);

    // there must be a comma or ']'
    if (std::regex_search(it, stop, match, COMMA_RX)) {
      it = match[0].second;
    } else if (!std::regex_search(it, stop, match, RSQRBRK_RX)) {
      // throw std::runtime_error(std::string("Expected ',' or ']', found '") + *it + "'");
      throw UnexpectedCharacter(START, "Array", "',' or ']'");
    }
  }

  it = match[0].second;
  return AliFemtoConfigObject(std::move(result));
}


// template <typename Iter_t>
AliFemtoConfigObject
AliFemtoConfigObject::ParseMap(StringIter_t& it, const StringIter_t stop)
// parse_struct(Iter_t& it, const Iter_t stop)
{
  StringIter_t START = it;
  // std::cout << "Parsing struct " << *it << "\n";
  AliFemtoConfigObject::MapValue_t result_map;

  std::smatch match;

  if (!std::regex_search(it, stop, match, LBRK_RX)) {
      throw std::runtime_error(std::string("Expected to match left-bracket character '{', found '") + *it + "'");
  }

  it = match[0].second;

  // loop until we match right bracket }
  while (!std::regex_search(it, stop, match, RBRK_RX)) {
    // load key into match[0]
    if (!std::regex_search(it, stop, match, KEY_RX)) {
      throw UnexpectedCharacter(START, "Map", "Expected beginning of key or '}' char");
    }

    auto keys = parse_map_key_unchecked(match[0].first, match[0].second);

    if (!std::regex_search(match[0].second, stop, match, COLON_RX)) {
      throw UnexpectedCharacter(START, "Map", "Expected colon ':' after key " + match[0].str());
    }
    // 'it' now should point to beginning of value
    it = match[0].second;
    AliFemtoConfigObject value = Parse(it, stop);

    if (std::regex_search(it, stop, match, COMMA_RX)) {
      it = match[0].second;
    } else if (!std::regex_search(it, stop, match, RBRK_RX)) {
      throw UnexpectedCharacter(START, "Map", "Expected closing brace '}'");
    }

    // push key-value pairs into back
    if (keys.size() > 1) {
      auto key_it = keys.rbegin(),
           key_it_stop = keys.rend();

      AliFemtoConfigObject::MapValue_t tmp_map;

      // loop until last one
      for (; std::next(key_it) != key_it_stop; ++key_it) {
        tmp_map.emplace(*key_it, std::move(value)); // value should now be empty
        value = std::move(tmp_map); // tmp_map should now be empty - value is new map object
      }
    }

    result_map.emplace(keys.front(), std::move(value));
  }

  it = match[0].second;
  return AliFemtoConfigObject(std::move(result_map));
}

AliFemtoConfigObject
match_to_rangelist(const std::smatch& m)
{
  AliFemtoConfigObject::RangeListValue_t result;
  using Float_t = AliFemtoConfigObject::RangeListValue_t::value_type::first_type;

  const std::regex comma_re(COMMA_PAT),
                   colon_re(COLON_PAT);

  std::sregex_token_iterator group(m[1].first, m[1].second, comma_re, -1),
                             stop;

  for (; group != stop; group++) {

    std::sregex_token_iterator num(group->first, group->second, colon_re, -1);

    Float_t first = std::stod(*num++);

    for (num++; num != stop; ++num) {
      Float_t second = std::stod(*num);
      result.emplace_back(first, second);
      first = second;
    }
  }

  return AliFemtoConfigObject(result);
}

AliFemtoConfigObject
match_to_range(const std::smatch &match)
{
  AliFemtoConfigObject::FloatValue_t range_start = std::stod(match[1].str()),
                                      range_stop = std::stod(match[2].str());
  return AliFemtoConfigObject(range_start, range_stop);
}

AliFemtoConfigObject
AliFemtoConfigObject::Parse(StringIter_t& it, const StringIter_t stop)
{
  std::smatch match;
  // std::cout << "Parsing value:\n";

  // skip leading whitespace
  while (it != stop && std::isspace(*it)) {
    ++it;
  }

  if (it == stop) {
    // throw std::runtime_error("Unexpected end of string.");
    throw EmptySourceString(it);
  }

  else if (std::regex_search(it, stop, match, RANGE_RX)) {
    it = match[0].second;
    return match_to_range(match);
  }

  else if (std::regex_search(it, stop, match, BOOL_RX)) {
    it = match[0].second;
    return AliFemtoConfigObject(match[0].str() == "true");
  }

  else if (std::regex_search(it, stop, match, FLT_RX)) {
    it = match[0].second;
    return AliFemtoConfigObject(std::stod(match.str()));
  }

  else if (std::regex_search(it, stop, match, INT_RX)) {
    it = match[0].second;
    return AliFemtoConfigObject(std::stoi(match.str()));
  }

  else if (std::regex_search(it, stop, match, STR_RX)) {
    it = match[0].second;
    return AliFemtoConfigObject(std::string(match[0].first + 1, it - 1));
  }

  else if (std::regex_search(it, stop, match, RANGELIST_RX)) {
    it = match[0].second;
    return match_to_rangelist(match);
  }

  else if (*it == '{') {
    return ParseMap(it, stop);
  }

  else if (*it == '[') {
    return ParseArray(it, stop);
  }

  throw UnexpectedLeadingChar(it, "undetermined-type", "");
}

namespace std {
  template<>
  struct hash<AliFemtoConfigObject> {
    using Range_t = AliFemtoConfigObject::RangeValue_t;

    ULong_t operator()(Range_t const& pair) const {
      return std::hash<Range_t::first_type>{}(pair.first)
             ^ (std::hash<Range_t::second_type>{}(pair.second) << 1);
    }

    ULong_t operator()(AliFemtoConfigObject::ArrayValue_t const& array) const {
      using Item_t = AliFemtoConfigObject::ArrayValue_t::value_type;

      return std::accumulate(
        array.cbegin(), array.cend(), 0,
        [] (ULong_t hash, const Item_t &obj) { return obj.Hash() ^ (hash << 1); });
    }

    ULong_t operator()(AliFemtoConfigObject::RangeListValue_t const& list) const {
      using Item_t = AliFemtoConfigObject::RangeListValue_t::value_type;

      return std::accumulate(
        list.cbegin(), list.cend(), 0,
        [&] (ULong_t a, const Item_t &pair) {
          return operator()(pair) ^ (a << 1); });
    }
  };
}


ULong_t
AliFemtoConfigObject::Hash() const
{
  static const ULong_t MAGIC_HASH_NUMBER = 3527539;

  ULong_t result = MAGIC_HASH_NUMBER ^ std::hash<int>{}(fTypeTag);

  switch (fTypeTag) {
    case kEMPTY: return result;
    case kBOOL: return result ^ std::hash<BoolValue_t>{}(fValueBool);
    case kINT: return result ^ std::hash<IntValue_t>{}(fValueInt);
    case kFLOAT: return result ^ std::hash<FloatValue_t>{}(fValueFloat);
    case kSTRING: return result ^ std::hash<StringValue_t>{}(fValueString);
    case kRANGE: return result ^ std::hash<AliFemtoConfigObject>{}(fValueRange);
    case kRANGELIST: return result ^ std::hash<AliFemtoConfigObject>{}(fValueRangeList);

    case kARRAY:

    case kMAP:
      return std::accumulate(
        fValueMap.cbegin(), fValueMap.cend(), result,
        [] (ULong_t hash, const MapValue_t::value_type &pair) {
          return hash
                 ^ (std::hash<Key_t>{}(pair.first) << 1)
                 ^ (pair.second.Hash() << 2); });
  }

  return 0;
}


std::string AliFemtoConfigObject::as_JSON_string() const
{
  if (is_empty()) {
    return "null";
  }

  std::stringstream ss;
  if (is_map()) {
    ss << "{";
    std::string seperator = "";
    for (const auto &pair : fValueMap) {
      ss << seperator << "\"" << pair.first << "\": "
         << pair.second.as_JSON_string();
      seperator = ", ";
    }
    ss << "}";
  }
  else if (is_array()) {
    ss << "[";
    std::string seperator = "";
    for (const auto &item : fValueArray) {
      ss << seperator << item.as_JSON_string();
      seperator = ", ";
    }
    ss << "]";
  }
  else if (is_range()) {
    ss << "{\"MIN\": " << std::get<0>(fValueRange)
       << ", \"MAX\": " << std::get<1>(fValueRange)
       << "}";
  }
  else if (is_rangelist()) {
    ss << "[";
    std::string seperator = "";
    for (const auto &item : fValueRangeList) {
      ss << "{\"MIN\": " << std::get<0>(item)
         << ", \"MAX\": " << std::get<1>(item)
         << "}";
      seperator = ", ";
    }
    ss << "]";
  }
  else if (is_int()) {
    ss << fValueInt;
  }
  else if (is_float()) {
    ss << fValueFloat;
  }
  else if (is_bool()) {
    ss << (fValueBool ? "true" : "false");
  }
  else if (is_str()) {\
    ss << '"' << fValueString << '"';
  }

  return ss.str();
}

AliFemtoConfigObject::MapValue_t::iterator
AliFemtoConfigObject::find_item(const AliFemtoConfigObject::Key_t key)
{
  auto it = key.begin(),
       stop = key.end();

  std::smatch match;
  if (!is_map() || !std::regex_search(it, stop, match, KEY_RX)) {
    return AliFemtoConfigObject::MapValue_t::iterator();
  }

  auto keys = parse_map_key_unchecked(match[0].first, match[0].second);
  auto next_key = keys.begin();

  AliFemtoConfigObject::MapValue_t::iterator found = fValueMap.find(*next_key);

  while (++next_key != keys.end()) {
    if (!found->second.is_map()) {
      return AliFemtoConfigObject::MapValue_t::iterator();
    }
    found = found->second.fValueMap.find(*next_key);
  }

  return found;
}


AliFemtoConfigObject
AliFemtoConfigObject::WithoutKey(const Key_t &key) const
{
  AliFemtoConfigObject result(*this);

  AliFemtoConfigObject *ptr = &result,
                       *container = nullptr;

  auto keys = split_key(key);

  for (auto &subkey : keys) {
    if (ptr->is_map()) {
      auto found = ptr->fValueMap.find(subkey);
      if (found == ptr->fValueMap.end()) {
        container = nullptr;
        break;
      }
      container = ptr;
      ptr = &found->second;
    } else {
      container = nullptr;
      break;
    }
  }

  if (container) {
    auto it = container->fValueMap.find(keys.back());
    container->fValueMap.erase(it);
  }

  return result;
}

AliFemtoConfigObject
AliFemtoConfigObject::WithoutKeys(const std::vector<Key_t> &keys) const
{
  AliFemtoConfigObject result(*this);
  for (auto &key : keys) {
    delete result.pop(key);
  }
  return result;
}
