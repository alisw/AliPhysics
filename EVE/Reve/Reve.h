// $Header$

#ifndef REVE_Reve_H
#define REVE_Reve_H

// #include <string>
#include <exception>
#include <TString.h>
#include <TError.h>
#include <Gtypes.h>

class TVirtualPad;
class TGeoManager;

namespace Reve {

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

void WarnCaller(const TString& warning);

void     ColorFromIdx(Color_t ci, UChar_t* col);
Color_t* FindColorVar(TObject* obj, const Text_t* varname);

/**************************************************************************/
/**************************************************************************/

void SetupEnvironment();

Bool_t CheckMacro(const Text_t* mac);
void   AssertMacro(const Text_t* mac);
void   Macro(const Text_t* mac);
void   LoadMacro(const Text_t* mac);

/**************************************************************************/
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

}

#endif
