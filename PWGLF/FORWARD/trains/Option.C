/**
 * @file   Option.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Tue Oct 16 18:59:04 2012
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_trains_util
 */
#ifndef OPTION_C
#define OPTION_C
#include <TNamed.h>
#include <iomanip>
#ifndef __CINT__
# include <TString.h>
# include <TList.h>
# include <TObjArray.h>
# include <TObjString.h>
# include <TMath.h>
# include <iostream>
# include <iomanip>
# include <cstdarg>
# include <cstdio>
#else
class TString;
class TList;
class TObjArray;
#endif

/** 
 * An option.  The value is stored as a string 
 *
 * @ingroup pwglf_forward_trains_util
 */
struct Option /* : public TNamed */
{
  /** 
   * Constructor 
   * 
   * @param name           Name of option
   * @param arg            Dummy argument (possibly null)
   * @param description    Description
   * @param value          Default value 
   */
  Option(const TString& name, 
	 const TString& arg, 
	 const TString& description,
	 const TString& value)
    : fName(name), fDescription(description), 
      fValue(value), fArg(arg), fIsSet(false)
  {}
  /** 
   * Copy constructor 
   * 
   * @param other Object to copy from 
   */
  Option(const Option& other)
    : fName(other.fName), 
      fDescription(other.fDescription),
      fValue(other.fValue), 
      fArg(other.fArg), 
      fIsSet(other.fIsSet)
  {}
  /** 
   * Assignment operator 
   * 
   * @param other Object to assign from 
   *
   * @return Reference to this object
   */
  Option& operator=(const Option& other)
  {
    if (&other == this) return *this;
    
    fName        = other.fName;
    fDescription = other.fDescription;
    fValue       = other.fValue;
    fArg         = other.fArg;
    fIsSet       = other.fIsSet;
    
    return *this;
  }
  /** 
   * Set the value 
   * 
   * @param val New value 
   */
  void Set(const TString& val) 
  { 
    if (HasArg()) {
      fIsSet = true;
      fValue = val;
      // Info("Set", "Setting option %s with arg %s to %s", 
      // fName.Data(), fArg.Data(), fValue.Data());
      return;
    }

    // Info("Set", "Setting option %s with no arg", fName.Data());
    // Allow flags to get =true, =1, =false, =0 values
    if (!val.IsNull() && 
	(val.EqualTo("false", TString::kIgnoreCase) || 
	 (val.IsDigit() && !val.Atoi()))) {
      fIsSet = false;
      fValue = "false";
    }
    else {
      fIsSet = true;
      fValue = "true";
    }
  }
  /** 
   * Set the value
   */
  void Set() 
  { 
    if (HasArg()) { 
      Error("Option::Set", "Option %s needs an argument", fName.Data());
      return;
    }
    Set("true");
  }
  /** 
   * Reset the set flag
   * 
   */
  void Reset() 
  { 
    fIsSet = false; 
    if (!HasArg()) fValue = "false";
  }
  /**
   * @return constant reference to value 
   */
  const TString& Get() const { return fValue; }
  /**
   * @return true if this option was set externally
   */
  Bool_t IsSet() const { return fIsSet; }
  /**
   * @return true if this option needs an argument
   */
  Bool_t HasArg() const { return !fArg.IsNull(); }
  /**
   * @return value as a boolean value
   */
  Bool_t AsBool() const { return fIsSet; }
  /**
   * @return value as an integer value
   */
  Int_t    AsInt() const { return fValue.Atoi(); }
  /**
   * @return value as a long integer value
   */
  Long64_t AsLong() const { return fValue.Atoll(); } 
  /**
   * @return value as a double precision value
   */
  Double_t AsDouble() const { return fValue.Atof(); } 
  /**
   * @return value as a C string
   */
  const char* AsString() const { return fValue.Data(); } 
  /** 
   * @return Width of the name 
   */
  Int_t NameWidth() const { return fName.IsNull() ? 0 : fName.Length(); }
  /** 
   * @return the width of the dummy argument
   */
  Int_t ArgWidth() const { return fArg.IsNull() ? 0 : fArg.Length(); }
  /** 
   * Show help 
   * 
   * @param o Output stream 
   * @param w With of option plus argument 
   */
  void Help(std::ostream& o, Int_t w=-1) const
  {
    TString tmp(fName);
    if (HasArg()) { 
      tmp.Append("=");
      tmp.Append(fArg);
    }
    if (w <= 0) w = NameWidth() + ArgWidth() + 1;
    std::ios::fmtflags oldf = o.setf(std::ios::left);
    o << std::setw(w) << tmp << "  " << fDescription << " [";
    if (!HasArg()) o << (IsSet() ? "true" : "false");
    else           o << fValue;
    o << "]" << std::endl;
    o.setf(oldf);
  }
  /** 
   * Show the option
   * 
   * @param o Output stream
   * @param w With of name 
   */
  void Show(std::ostream& o, Int_t w=-1) const 
  {
    if (w <= 0) w = NameWidth();
    std::ios::fmtflags oldf = o.setf(std::ios::left);
    o << std::setw(w) << fName << ": ";
    if (!HasArg()) o << (IsSet() ? "true" : "false");
    else           o << fValue;
    o << std::endl;
    o.setf(oldf);
  }
  /** 
   * Store option and possible value 
   * 
   * @param o Output stream
   * @param quote If true, quote output 
   */
  void Store(std::ostream& o, bool quote=true) const 
  {
    o << fName;
    if (!HasArg()) return;
    o << "=" << (quote ? "\"" : "") << fValue << (quote ? "\"" : "");
  }
  TString fName;        // Name 
  TString fDescription; // Description
  TString fValue;       // Value 
  TString fArg;         // Argument dummy 
  Bool_t  fIsSet;       // True if this option was set 
  
  // ClassDef(Option,1) // Option 
};

/** 
 * A List of options 
 */
struct OptionList 
{
  // Linked list element
  struct Link 
  {
    Link*   fPrev;
    Link*   fNext;
    Option* fThis;
    Link() : fPrev(0), fNext(0), fThis(0) {}
    Link(Link* next, Option* opt)
      : fPrev(next ? next->fPrev : 0), // Set previous 
	fNext(next), // Set next link
	fThis(opt) // Set data
    {
      if (fPrev) fPrev->fNext = this; // Set forward link
      if (fNext) fNext->fPrev = this; // Set previous link
    }
    ~Link() { 
      if (fPrev) fPrev->fNext = fNext;
      if (fNext) fNext->fPrev = fPrev;
      delete fThis;
    }
    Link(const Link& o)
      : fPrev(o.fPrev), 
	fNext(o.fNext), 
	fThis(o.fThis)
    {
    }
    Link& operator=(const Link& o)
    {
      if (&o == this) return *this;
      fPrev = o.fPrev;
      fNext = o.fNext;
      fThis = o.fThis;
      return *this;
    }	
  };
  /** 
   * Constructor 
   */
  OptionList() : fList(0) { }
  /** 
   * Copy constructor 
   *
   * @param other Object to copy from 
   */
  OptionList(const OptionList& other) 
    : fList(0)
  { 
    // fList.SetOwner(); 
    // TIter next(&other.fList);
    // Option* o = 0;
    // while ((o = static_cast<Option*>(next()))) 
    //   fList.Add(new Option(*o));
    Copy(other);
  }
  /** 
   * Destructor 
   */
  ~OptionList() { Delete(); }

  /** 
   * Remove all options 
   */
  void Delete()
  { 
    Link* cur = fList;
    while (cur) {
      Link* tmp = cur->fNext;
      delete cur;
      cur = tmp;
      // Remove(cur->fThis->fName);
      // cur = tmp;
    }
    fList = 0;
  }
  /** 
   * Assignment operator 
   *
   * @param other  Object to assign from 
   *
   * @return reference to this 
   */
  OptionList& operator=(const OptionList& other) 
  { 
    if (&other == this) return *this;
    Delete();
    Copy(other);
    
    return *this; 
  }
  /** 
   * Copy list from other object
   * 
   * @param other Object to copy from 
   */
  void Copy(const OptionList& other)
  {
    Delete();
    const Link* ocur = other.fList;
    Link*       cur  = fList;
    Link*       prev = fList;
    while (ocur) { 
      cur        = new Link;
      cur->fThis = new Option(*(ocur->fThis));
      cur->fNext = 0;
      cur->fPrev = prev;
      if (fList == 0) fList = cur;
      if (prev)       prev->fNext = cur;
      prev = cur;
      ocur = ocur->fNext;
    }
  }
  void DebugLink(const Link* link) const
  {
    std::cout << "Link=" << link;
    if (link) {
      std::cout << " prev=" << link->fPrev 
		<< " next=" << link->fNext
		<< " obj=" << link->fThis;
      if (link->fThis) 
        std::cout << " name=" << link->fThis->fName;
    }
    std::cout <<std::endl;
  }
  /** 
   * Find an optio by name 
   * 
   * @param name Name of option to find
   * 
   * @return Pointer to option or null
   */
  Option* Find(const TString& name) const
  {
    const Link* cur = fList;
    // DebugLink(cur);
    while (cur && cur->fThis) { 
      if (name.EqualTo(cur->fThis->fName)) return cur->fThis;
      cur = cur->fNext;
    }
    return 0;
  }
  /** 
   * Add an option with argument 
   * 
   * @param name Name of option
   * @param arg  Dummy argument
   * @param desc Description
   * @param val  Default value 
   * 
   * @return Newly added option 
   */
  Option* Add(const TString& name, 
	      const TString& arg, 
	      const TString& desc, 
	      const TString& val="")
  {
    Option* o = Find(name);
    if (o) { 
      Warning("OptionList::Add", "Option %s already registered", name.Data());
      return o;
    }
    // Info("Add", "New option %s with arg %s (%s) and value %s",    
    //	 name.Data(), arg.Data(), desc.Data(), val.Data());
    o = new Option(name, arg, desc, val);
    Link*  cur  = fList;
    if (!cur) { 
      cur = new Link;
      cur->fThis = o;
      cur->fNext = 0;
      cur->fPrev = 0;
      fList      = cur;
    }
    else {
      Link* n = 0;
      Link* l = 0;
      while (cur) { 
	if (cur->fThis->fName.CompareTo(name) < 0) { 
	  l   = cur;
	  cur = cur->fNext;
	  continue;
	}
	n = new Link(cur, o);
	if (cur == fList) fList = n;
	break;
      }
      if (!n) { 
	n = new Link;
	l->fNext = n;
	n->fPrev = l;
	n->fNext = 0;
	n->fThis = o;
      }
    }
    return o;
  }
  /** 
   * Add an option with no argument
   * 
   * @param name Name of option
   * @param desc Description 
   * 
   * @return Newly created option 
   */
  Option* Add(const TString& name, 
	      const TString& desc)
  {
    return Add(name, "", desc, "");
  }
  /** 
   * Add an option with no argument
   * 
   * @param name Name of option
   * @param desc Description 
   * @param def  Default value (true, or false)
   * 
   * @return Newly created option 
   */
  Option* Add(const TString& name, 
	      const TString& desc,
	      Bool_t         def)
  {
    Option* o = Add(name, "", desc, "");
    if (o) o->Set(def ? "true" : "false");
    return o;
  }
  /** 
   * Add an option with argument 
   * 
   * @param name Name of option
   * @param arg  Dummy argument
   * @param desc Description
   * @param val  Default value 
   * @param asHex If true, interpret as hex number 
   * 
   * @return Newly added option 
   */
  Option* Add(const TString& name, 
	      const TString& arg, 
	      const TString& desc, 
	      Int_t          val, 
	      Bool_t         asHex=false)
  {
    if (asHex) {
      UInt_t uval = val;
      return Add(name, arg, desc, Form("0x%x", uval));
    }
    return Add(name, arg, desc, Form("%d", val));
  }
  /** 
   * Add an option with argument 
   * 
   * @param name Name of option
   * @param arg  Dummy argument
   * @param desc Description
   * @param val  Default value 
   * @param asHex If true, interpret as hex 
   * 
   * @return Newly added option 
   */
  Option* Add(const TString& name, 
	      const TString& arg, 
	      const TString& desc, 
	      Long64_t       val, 
	      Bool_t         asHex=false)
  {
    if (asHex) {
      ULong64_t uval = val;
      return Add(name, arg, desc, Form("0x%llx", uval));
    }
    return Add(name, arg, desc, Form("%lld", val));
  }
  /** 
   * Add an option with argument 
   * 
   * @param name Name of option
   * @param arg  Dummy argument
   * @param desc Description
   * @param val  Default value 
   * 
   * @return Newly added option 
   */
  Option* Add(const TString& name, 
	      const TString& arg, 
	      const TString& desc, 
	      Double_t val)
  {
    return Add(name, arg, desc, Form("%lg", val));
  }

  /** 
   * Remove an option 
   * 
   * @param name Name of option to remove 
   */
  void Remove(const TString& name)
  {
    Link* cur = fList;
    while (cur) { 
      if (!cur->fThis->fName.EqualTo(name)) {
	cur = cur->fNext;
	continue;
      }
      if (fList == cur) fList = cur->fNext;
      delete cur;
      break;
    }
  }
  /** 
   * Check if a given option was set externally
   * 
   * @param name Name of option
   * 
   * @return true if option exists and was set externally 
   */
  Bool_t Has(const TString& name) const
  {
    Option* o = Find(name);
    return (o && o->IsSet());
  }
  /** 
   * Get the value of an option
   * 
   * @param name Name of option
   * 
   * @return Value of option, or empty string 
   */
  const TString& Get(const TString& name) const
  {
    static TString null("");
    Option* o = Find(name);
    if (!o) return null;
    return o->Get();
  }
  /** 
   * Get a value using a format statement. Remember argument(s) must
   * be passed by address (as pointers)
   * 
   * @param name Name of option @param format Format statement.
   * Remeber, double and long needs the "l" modifier
   * 
   * @return true on success
   */
  Bool_t GetF(const TString& name, const Char_t* format, ...) const
  {
    Option* o = Find(name);
    if (!o) return false;
    
    va_list ap;
    va_start(ap, format);
    int ret = vsscanf(o->fValue.Data(), format, ap);
    va_end(ap);
    
    return ret > 0;
  }
  /** 
   * Get value of an option as a boolean 
   * 
   * @param name Name of 
   * 
   * @return Value or false if not found
   */
  Bool_t AsBool(const TString& name) const 
  {
    Option* o = Find(name);
    if (!o) return false;
    return o->AsBool();
  }
  /** 
   * Return value of an option as an integer
   * 
   * @param name Name of option 
   * @param def  Default value if options isn't found 
   * 
   * @return Value or default value 
   */
  Int_t AsInt(const TString& name, Int_t def=0) const
  {
    Option* o = Find(name);
    if (!o) return def;
    return o->AsInt();
  }
  /** 
   * Return value of an option as an integer
   * 
   * @param name Name of option 
   * @param def  Default value if options isn't found 
   * 
   * @return Value or default value 
   */
  Long64_t AsLong(const TString& name, Long64_t def=0) const
  {
    Option* o = Find(name);
    if (!o) return def;
    return o->AsLong();
  }
  /** 
   * Return value of an option as a double precision real number
   * 
   * @param name Name of option 
   * @param def  Default value if options isn't found 
   * 
   * @return Value or default value 
   */
  Double_t AsDouble(const TString& name, Double_t def=0) const
  {
    Option* o = Find(name);
    if (!o) return def;
    return o->AsDouble();
  }
  /** 
   * Set value using a format statement 
   * 
   * @param name   Name of option.
   * @param format Format statement
   */
  void SetF(const TString& name, const Char_t* format, ...)
  {
    Option* o = Find(name);
    if (!o) return;
    
    static char buf[1024];
    va_list ap;
    
    va_start(ap, format);
    vsnprintf(buf, 1023, format, ap);
    buf[1023] = '\0';
    va_end(ap);
    
    o->Set(buf);
  }
  /** 
   * Set an option
   * 
   * @param name   Name of option
   * @param value  Value of option
   */
  void Set(const TString& name, const TString& value)
  {
    Option* o = Find(name);
    if (!o) return;
    o->Set(value);
  }
  /** 
   * Set a flag 
   * 
   * @param name Name of flag 
   */
  void Set(const TString& name)
  {
    Option* o = Find(name);
    if (!o) return;
    o->Set();
  }   
  /** 
   * Set long integer value 
   * 
   * @param name Name of option
   * @param val  Value 
   * @param asHex If true, interpret as hex 
   */
  void Set(const TString& name, Int_t val, Bool_t asHex=false)
  {
    if (asHex) Set(name, Form("0x%x", val));
    else       Set(name, Form("%d", val));
  }
  /** 
   * Set long integer value 
   * 
   * @param name Name of option
   * @param val  Value 
   * @param asHex If true, interpret as hex value 
   */
  void Set(const TString& name, Long64_t val, Bool_t asHex=false)
  {
    ULong64_t uval = val;
    if (asHex) Set(name, Form("0x%llx", uval));
    else       Set(name, Form("%lld", val));
  }
  /** 
   * Set double precision floating point value 
   * 
   * @param name Name of option
   * @param val  Value 
   */
  void Set(const TString& name, Double_t val) 
  {
    Set(name, Form("%lg", val));
  }
  /** 
   * Parse the options given in tmp 
   * 
   * @param tmp     String to pass
   * @param delims  Delimiters 
   * 
   * @return true on success 
   */
  Bool_t Parse(const TString& tmp, const TString& delims)
  {
    TObjArray* opts = tmp.Tokenize(delims);
    // Info("OptionList::Parse", "Parsing options %s", tmp.Data());
    Bool_t ret = Parse(opts);
    opts->Delete();
    return ret;
  }
  /** 
   * Parse options given in a collection 
   * 
   * @param opts List of arguments 
   * @param ignoreUnknown If true, ignore unknown options 
   * 
   * @return true on success
   */
  Bool_t Parse(const TCollection* opts, Bool_t ignoreUnknown=false)
  {
    // Info("OptionList::Parse", "List of options");
    // opts->ls();
    TIter       next(opts);
    TObjString* o = 0;
    while ((o = static_cast<TObjString*>(next()))) { 
      TString& s   = o->String();
      TString  key = s;
      TString  val = "";
      Int_t    eq  = s.Index("=");
      if (eq != kNPOS) {
	key = s(0, eq);
	val = s(eq+1, s.Length()-eq-1);
      }

      // Info("OptionList::Parse", "Looking for key=%s", key.Data());
      Option*  opt = Find(key);
      if (!opt) { 
	if (!ignoreUnknown)
	  Warning("OptionList::Parse", "Unknown option: \"%s\"", s.Data());
	continue;
      }
      if (opt->HasArg() && val.IsNull()) { 
	Warning("OptionList::Parse", 
		"Option %s needs an argument, using default %s", 
		key.Data(), opt->fValue.Data());
	val = opt->fValue;
	// return false;
      }
      opt->Set(val);
    }
    return true;
  }
  /** 
   * Find the widest name and dummy argument
   * 
   * @param nWidth On return, the largest width of option names 
   * @param aWidth On return, the largest width of option dummy args
   */
  void Widest(Int_t& nWidth, Int_t& aWidth) const 
  {
    nWidth = 0;
    aWidth = 0;
    const Link* cur = fList;
    while (cur) { 
      Option* opt = cur->fThis;
      nWidth = TMath::Max(nWidth, opt->NameWidth());
      aWidth = TMath::Max(aWidth, opt->ArgWidth());
      cur = cur->fNext;
    }

    // while ((opt = static_cast<Option*>(next()))) {
    // }
  }
  /** 
   * Display option help
   * 
   * @param o Output stream 
   * @param prefix Prefix for each option.
   */
  void Help(std::ostream& o, const char* prefix="  ") const 
  {
    Int_t nWidth, aWidth;
    Widest(nWidth, aWidth);
    if (aWidth > 0) nWidth += aWidth+1;

    const Link* cur = fList;
    while (cur) { 
      Option* opt = cur->fThis;
      o << prefix;
      opt->Help(o, nWidth);
      cur = cur->fNext;
    }
  }
  /** 
   * Show the values of options 
   * 
   * @param o Output stream
   * @param prefix Prefix for each option
   */
  void Show(std::ostream& o, const char* prefix="  ") const 
  {
    Int_t nWidth, aWidth;
    Widest(nWidth, aWidth);


    const Link* cur = fList;
    while (cur) { 
      Option* opt = cur->fThis;
      o << prefix;
      opt->Show(o, nWidth);
      cur = cur->fNext;
    }
  }
  /** 
   * Show the values of options 
   * 
   * @param o Output stream
   * @param prefix Prefix for each option
   * @param delim Delimters 
   * @param quote Quote output 
   * @param onlySet if true, only output set options 
   */
  void Store(std::ostream& o, const char* prefix="", 
	     const char* delim=",", bool quote=true, 
	     bool onlySet=false) const 
  {
    Int_t nWidth, aWidth;
    Widest(nWidth, aWidth);

    const Link* cur = fList;
    while (cur) { 
      Option* opt = cur->fThis;
      if ((!opt->HasArg() || onlySet) && !opt->IsSet()) {
	cur = cur->fNext;
	continue;
      }
      o << prefix;
      opt->Store(o, quote);
      o << delim;
      cur = cur->fNext;
    }
  }
  // Our linked list
  Link* fList;

  static void Test(const char* opts="")
  {
    OptionList l;
    l.Add("int", "NUMBER", "Integer", "42");
    l.Add("float", "NUMBER", "Floating point", "3.14");
    l.Add("bool", "Flag");
    l.Add("hex", "NUMBER", "Hexadecimal", "0xdead");
    l.Add("string", "STRING", "A string", "Hello, world");
    l.Show(std::cout, "\t");
    
    std::cout << "Find" << std::endl;
    Option* b = l.Find("bool");
    b->Set("true");
    b->Show(std::cout);

    std::cout << "SetF" << std::endl;
    l.SetF("float", "%f", 2.17);
    l.Show(std::cout, "\t");
    
    std::cout << "GetF" << std::endl;
    Float_t f;
    l.GetF("float", "%f", &f);
    std::cout << "\tf=" << f << std::endl;
    std::cout << "Remove" << std::endl;
    l.Remove("float");
    l.Show(std::cout, "\t");

    std::cout << "Set" << std::endl;
    l.Set("int", "10");
    l.Set("hex", 0xbeef, true);
    l.Set("bool", "false");
    l.Show(std::cout, "\t");

    std::cout << "Copy" << std::endl;
    OptionList c(l);
    c.Show(std::cout, "\t");

    std::cout << "Parse" << std::endl;
    c.Parse(opts,",");
    c.Show(std::cout, "\t");
    std::cout << "End of test" << std::endl;
  }
  // TList fList;
};

#endif

