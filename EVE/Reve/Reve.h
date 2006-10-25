// $Header$

#ifndef REVE_Reve_H
#define REVE_Reve_H

// #include <string>
#include <exception>
#include <TObject.h>
#include <TString.h>
#include <TError.h>
#include <Gtypes.h>

class TVirtualPad;
class TGeoManager;

namespace Reve {


/**************************************************************************/
// Exceptions, string functions
/**************************************************************************/

bool operator==(const TString& t, const std::string& s);
bool operator==(const std::string& s, const TString& t);

class Exc_t : public std::exception, public TString
{
 public:
  Exc_t() {}
  Exc_t(const TString& s) : TString(s) {}
  Exc_t(const char* s)    : TString(s) {}
  Exc_t(const std::string& s);

  virtual ~Exc_t() throw () {}

  virtual const char* what() const throw () { return Data(); }

  ClassDef(Exc_t, 1);
};

Exc_t operator+(const Exc_t &s1, const std::string  &s2);
Exc_t operator+(const Exc_t &s1, const TString &s2);
Exc_t operator+(const Exc_t &s1, const char    *s2);

void  WarnCaller(const TString& warning);


/**************************************************************************/
// Environment, Macro functions
/**************************************************************************/

void   SetupEnvironment();

Bool_t CheckMacro(const Text_t* mac);
void   AssertMacro(const Text_t* mac);
void   Macro(const Text_t* mac);
void   LoadMacro(const Text_t* mac);


/**************************************************************************/
// Local cache for global Pad, GeoManager 
/**************************************************************************/

TVirtualPad* PushPad(TVirtualPad* new_gpad=0, Int_t subpad=0);
TVirtualPad* PopPad(Bool_t modify_update_p=false);

class PadHolder
{
private:
  Bool_t fModifyUpdateP;
public:
  PadHolder(Bool_t modify_update_p, TVirtualPad* new_gpad=0, Int_t subpad=0) :
    fModifyUpdateP(modify_update_p)
  { PushPad(new_gpad, subpad); }

  virtual ~PadHolder() { PopPad(fModifyUpdateP); }

  ClassDef(PadHolder, 0);
};

class GeoManagerHolder
{
private:
  TGeoManager* fManager;

  GeoManagerHolder(const GeoManagerHolder&);            // Not implemented
  GeoManagerHolder& operator=(const GeoManagerHolder&); // Not implemented

public:
  GeoManagerHolder(TGeoManager* new_gmgr=0);
  virtual ~GeoManagerHolder();

  ClassDef(GeoManagerHolder, 0);
};


/**************************************************************************/
// ReferenceCount base-class (interface)
/**************************************************************************/

class ReferenceCount
{
protected:
  Int_t fRefCount;

public:
  ReferenceCount() : fRefCount(0) {}
  virtual ~ReferenceCount() {}

  void IncRefCount() { ++fRefCount; }
  void DecRefCount() { if(--fRefCount <= 0) OnZeroRefCount(); }

  virtual void OnZeroRefCount() { delete this; }

  ClassDef(ReferenceCount, 0);
};


/**************************************************************************/
// Color, palette management
/**************************************************************************/

void     ColorFromIdx(Color_t ci, UChar_t* col, Bool_t alpha=kTRUE);
void     ColorFromIdx(Float_t f1, Color_t c1, Float_t f2, Color_t c2,
		      UChar_t* col, Bool_t alpha=kTRUE);
Color_t* FindColorVar(TObject* obj, const Text_t* varname);

class RGBAPalette : public TObject, public ReferenceCount
{
protected:
  Int_t     fMinVal;
  Int_t     fMaxVal;
  Int_t     fNBins;
  Bool_t    fInterpolate;
  Bool_t    fWrap;

  mutable UChar_t* fColorArray;

  void SetupColor(Int_t val, UChar_t* pix) const;
  
public:
  RGBAPalette();
  RGBAPalette(Int_t min, Int_t max);
  RGBAPalette(Int_t min, Int_t max, Bool_t interp, Bool_t wrap);
  virtual ~RGBAPalette();

  void SetupColorArray() const;
  void ClearColorArray();

  UChar_t* ColorFromArray(Int_t val) const;
  void     ColorFromArray(Int_t val, UChar_t* pix, Bool_t alpha=kTRUE) const;

  Int_t  GetMinVal() const        { return fMinVal; }
  Int_t  GetMaxVal() const        { return fMaxVal; }
  Bool_t GetInterpolate() const   { return fInterpolate; }
  Bool_t GetWrap() const          { return fWrap; }

  void   SetMinMax(Int_t min, Int_t max);
  void   SetInterpolate(Bool_t b);
  void   SetWrap(Bool_t b);

  ClassDef(RGBAPalette, 1)
};


inline UChar_t* RGBAPalette::ColorFromArray(Int_t val) const
{
  if(!fColorArray)  SetupColorArray();
  if(val < fMinVal) val = fWrap ? ((val+1-fMinVal)%fNBins + fMaxVal) : fMinVal;
  if(val > fMaxVal) val = fWrap ? ((val-1-fMaxVal)%fNBins + fMinVal) : fMaxVal;
  return fColorArray + 4 * (val - fMinVal);
}

inline void RGBAPalette::ColorFromArray(Int_t val, UChar_t* pix, Bool_t alpha) const
{
  UChar_t* c = ColorFromArray(val);
  pix[0] = c[0]; pix[1] = c[1]; pix[2] = c[2];
  if (alpha) pix[3] = c[3];
}

/**************************************************************************/

}

#endif
