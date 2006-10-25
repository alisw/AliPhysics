// $Header$

#include "Reve.h"

#include <TError.h>
#include <TPad.h>
#include <TGeoManager.h>

#include <TStyle.h>
#include <TColor.h>

#include <TROOT.h>
#include <TInterpreter.h>

#include <list>
#include <string>
#include <iostream>

//______________________________________________________________________
// Reve
//

/**************************************************************************/
/**************************************************************************/
namespace Reve {
/**************************************************************************/


// TString .vs. string

bool operator==(const TString& t, const std::string& s)
{ return (s == t.Data()); }

bool operator==(const std::string&  s, const TString& t)
{ return (s == t.Data()); }

// Exc

Exc_t::Exc_t(const std::string& s) : TString(s.c_str()) {}

// Exc + ops

Exc_t operator+(const Exc_t &s1, const std::string &s2)
{ return Exc_t((std::string&)s1 + s2); }

Exc_t operator+(const Exc_t &s1, const TString &s2)
{ return Exc_t((std::string&)s1 + s2.Data()); }

Exc_t operator+(const Exc_t &s1,  const char *s2)
{ return Exc_t((std::string&)s1 + s2); }

// ----------------------------------------------------------------

void WarnCaller(const TString& warning)
{
  std::cout << "WRN: " << warning << std::endl;
}

/**************************************************************************/
/**************************************************************************/

void SetupEnvironment()
{
  // Check if REVESYS exists, try fallback to $ALICE_ROOT/EVE.
  // Setup Include and Macro paths.

  static const Exc_t eH("Reve::SetupEnvironment");
  static Bool_t setupDone = kFALSE;

  if (setupDone) {
    Info(eH.Data(), "has already been run.");
    return;
  }

  if(gSystem->Getenv("REVESYS") == 0) {
    if(gSystem->Getenv("ALICE_ROOT") != 0) {
      Info(eH.Data(), "setting REVESYS from ALICE_ROOT.");
      gSystem->Setenv("REVESYS", Form("%s/EVE", gSystem->Getenv("ALICE_ROOT")));
    } else {
      Error(eH.Data(), "REVESYS not defined, neither is ALICE_ROOT.");
      gSystem->Exit(1);
    }
  }
  if(gSystem->AccessPathName(gSystem->Getenv("REVESYS")) == kTRUE) {
    Error(eH.Data(), "REVESYS '%s' does not exist.", gSystem->Getenv("REVESYS"));
    gSystem->Exit(1);
  }

  TString macPath(gROOT->GetMacroPath());
  macPath += Form(":%s/macros", gSystem->Getenv("REVESYS"));
  gInterpreter->AddIncludePath(gSystem->Getenv("REVESYS"));
  if(gSystem->Getenv("ALICE_ROOT") != 0) {
    macPath += Form(":%s/alice-macros", gSystem->Getenv("REVESYS"));
    gInterpreter->AddIncludePath(Form("%s/include", gSystem->Getenv("ALICE_ROOT")));
    gInterpreter->AddIncludePath(gSystem->Getenv("ALICE_ROOT"));
  }
  gROOT->SetMacroPath(macPath);

  setupDone = kTRUE;
}

/**************************************************************************/

namespace {
  void ChompTail(TString& s, char c='.') {
    Ssiz_t p = s.Last(c);
    if(p != kNPOS)
      s.Remove(p);
  }
}

Bool_t CheckMacro(const Text_t* mac)
{
  // Checks if macro 'mac' is loaded.

  return gROOT->GetInterpreter()->IsLoaded(mac);

  // Previous version expected function with same name and used ROOT's
  // list of global functions.
  /*
  TString foo(mac); ChompTail(foo);
  if(recreate) {
    TCollection* logf = gROOT->GetListOfGlobalFunctions(kFALSE);
    logf->SetOwner();
    logf->Clear();
  }
  return (gROOT->GetGlobalFunction(foo.Data(), 0, kTRUE) != 0);
  */
}

void AssertMacro(const Text_t* mac)
{
  // Load and execute macro 'mac' if it has not been loaded yet.

  if(CheckMacro(mac) == kFALSE) {
    gROOT->Macro(mac);
  }
}

void Macro(const Text_t* mac)
{
  // Execute macro 'mac'. Do not reload the macro.

  if(CheckMacro(mac) == kFALSE) {
    gROOT->LoadMacro(mac);
  }
  TString foo(mac); ChompTail(foo); foo += "()";
  gROOT->ProcessLine(foo.Data());
}

void LoadMacro(const Text_t* mac)
{
  // Makes sure that macro 'mac' is loaded, but do not reload it.

  if(CheckMacro(mac) == kFALSE) {
    gROOT->LoadMacro(mac);
  }
}

/**************************************************************************/
/**************************************************************************/

// Pad stack for RINT/GUI thread.
std::list<TVirtualPad*> s_Pad_Stack;

TVirtualPad* PushPad(TVirtualPad* new_gpad, Int_t subpad)
{
  // printf("Reve::PushPad old=%p, new=%p\n", gPad, new_gpad);
  s_Pad_Stack.push_back(gPad);
  if(new_gpad != 0)
    new_gpad->cd(subpad);
  else
    gPad = 0;
  return gPad;
}

TVirtualPad* PopPad(Bool_t modify_update_p)
{
  // printf("Reve::PopPad old=%p, new=%p\n", gPad, s_Pad_Stack.empty() ? 0 : s_Pad_Stack.back());
  if(s_Pad_Stack.empty()) {
    Warning("Reve::PopTPad", "stack empty.");
  } else {
    if(modify_update_p && gPad != 0) {
      gPad->Modified();
      gPad->Update();
    }
    gPad = s_Pad_Stack.back();
    s_Pad_Stack.pop_back();
  }
  return gPad;
}

/**************************************************************************/
// 
/**************************************************************************/

GeoManagerHolder::GeoManagerHolder(TGeoManager* new_gmgr) :
  fManager(gGeoManager)
{
  gGeoManager = new_gmgr;
}

GeoManagerHolder::~GeoManagerHolder()
{
  gGeoManager = fManager;
}

/**************************************************************************/
// Color, palette management
/**************************************************************************/

void ColorFromIdx(Color_t ci, UChar_t* col, Bool_t alpha)
{
  if (ci < 0) {
    col[0] = col[1] = col[2] = col[3] = 0;
    return;
  }
  TColor* c = gROOT->GetColor(ci);
  if(c) { 
    col[0] = (UChar_t)(255*c->GetRed());  
    col[1] = (UChar_t)(255*c->GetGreen());
    col[2] = (UChar_t)(255*c->GetBlue()); 
    if (alpha) col[3] = 255;
  }
}

void ColorFromIdx(Float_t f1, Color_t c1, Float_t f2, Color_t c2,
		  UChar_t* col, Bool_t alpha)
{
  TColor* t1 = gROOT->GetColor(c1);
  TColor* t2 = gROOT->GetColor(c2);
  if(t1 && t2) { 
    col[0] = (UChar_t)(255*(f1*t1->GetRed()   + f2*t2->GetRed()));
    col[1] = (UChar_t)(255*(f1*t1->GetGreen() + f2*t2->GetGreen()));
    col[2] = (UChar_t)(255*(f1*t1->GetBlue()  + f2*t2->GetBlue())); 
    if (alpha) col[3] = 255;
  }
}

Color_t* FindColorVar(TObject* obj, const Text_t* varname)
{
  static const Exc_t eH("Reve::FindColorVar");

  Int_t off = obj->IsA()->GetDataMemberOffset(varname);
  if(off == 0)
    throw(eH + "could not find member '" + varname + "' in class " + obj->IsA()->GetName() + ".");
  return (Color_t*) (((char*)obj) + off);
}

/**************************************************************************/
// RGBAPalette
/**************************************************************************/

RGBAPalette::RGBAPalette() :
  TObject(),
  Reve::ReferenceCount(),
  fInterpolate (kFALSE),
  fWrap        (kFALSE),
  fColorArray  (0)
{
  SetMinMax(0, 100);
}

RGBAPalette::RGBAPalette(Int_t min, Int_t max) :
  TObject(),
  Reve::ReferenceCount(),
  fInterpolate (kFALSE),
  fWrap        (kFALSE),
  fColorArray  (0)
{
  SetMinMax(min, max);
}

RGBAPalette::RGBAPalette(Int_t min, Int_t max, Bool_t interp, Bool_t wrap) :
  TObject(),
  Reve::ReferenceCount(),
  fInterpolate (interp),
  fWrap        (wrap),
  fColorArray  (0)
{
  SetMinMax(min, max);
}

RGBAPalette::~RGBAPalette()
{
  delete [] fColorArray;
}

/**************************************************************************/

void RGBAPalette::SetupColor(Int_t val, UChar_t* pixel) const
{
  using namespace TMath;
  Float_t div  = Max(1, fMaxVal - fMinVal);
  Int_t   nCol = gStyle->GetNumberOfColors();

  Float_t f;
  if      (val <= fMinVal) f = 0;
  else if (val >= fMaxVal) f = nCol - 1;
  else                     f = (val - fMinVal)/div*(nCol - 1);

  if (fInterpolate) {
    Int_t  bin = (Int_t) f;
    Float_t f1 = f - bin, f2 = 1.0f - f1;
    ColorFromIdx(f1, gStyle->GetColorPalette(bin),
		 f2, gStyle->GetColorPalette(Min(bin + 1, nCol - 1)),
		 pixel);
  } else {
    ColorFromIdx(gStyle->GetColorPalette((Int_t) Nint(f)), pixel);    
  }
}

void RGBAPalette::SetupColorArray() const
{
  if(fColorArray)
    return;

  fColorArray = new UChar_t [4 * fNBins];
  UChar_t* p = fColorArray;
  for(Int_t v=fMinVal; v<=fMaxVal; ++v, p+=4)
    SetupColor(v, p);
}

void RGBAPalette::ClearColorArray()
{
  if(fColorArray) {
    delete [] fColorArray;
    fColorArray = 0;
  }
}

/**************************************************************************/

void RGBAPalette::SetMinMax(Int_t min, Int_t max)
{
  fMinVal = min;
  fMaxVal = max;
  fNBins  = max - min + 1;
  ClearColorArray();
}

void RGBAPalette::SetInterpolate(Bool_t b)
{
  fInterpolate = b;
  ClearColorArray();
}

void RGBAPalette::SetWrap(Bool_t b)
{
  fWrap = b;
}

/**************************************************************************/
} // end namespace Reve
/**************************************************************************/
/**************************************************************************/
