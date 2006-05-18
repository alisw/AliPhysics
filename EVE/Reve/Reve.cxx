// $Header$

#include "Reve.h"

#include <TError.h>
#include <TPad.h>
#include <TGeoManager.h>
#include <TColor.h>

#include <TROOT.h>

#include <list>

#include <iostream>

//______________________________________________________________________
// Reve
//

namespace Reve {

Exc_t operator+(const Exc_t &s1, const std::string &s2)
{ return Exc_t((std::string&)s1 + s2); }

Exc_t operator+(const Exc_t &s1, const TString &s2)
{ return Exc_t((std::string&)s1 + s2.Data()); }

Exc_t operator+(const Exc_t &s1,  const char *s2)
{ return Exc_t((std::string&)s1 + s2); }


void WarnCaller(const TString& warning)
{
  std::cout << "WRN: " << warning << std::endl;
}

void ColorFromIdx(Color_t ci, UChar_t* col)
{
  TColor* c = gROOT->GetColor(ci);
  if(c) { 
    col[0] = (UChar_t)(255*c->GetRed());  
    col[1] = (UChar_t)(255*c->GetGreen());
    col[2] = (UChar_t)(255*c->GetBlue()); 
    col[3] = 255;
  }
}

/**************************************************************************/
/**************************************************************************/

void SetupEnvironment()
{
  // Check REVESYS exists, try fallback to $ALICE_ROOT/EVE.

  static const Exc_t eH("Reve::SetupEnvironment");

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
}

/**************************************************************************/

namespace {
  void ChompTail(TString& s, char c='.') {
    Ssiz_t p = s.Last(c);
    if(p != kNPOS)
      s.Remove(p);
  }
}

void AssertMacro(const Text_t* mac)
{
  // Load and execute macro 'mac' if it has not been loaded yet.

  TString foo(mac); ChompTail(foo);
  if(gROOT->GetGlobalFunction(foo.Data(), 0, true) == 0) {
    gROOT->Macro(mac);
  }
}

void Macro(const Text_t* mac)
{
  // Execute macro 'mac'. Do not reload the macro.

  TString foo(mac); ChompTail(foo);
  if(gROOT->GetGlobalFunction(foo.Data(), 0, true) == 0)
    gROOT->LoadMacro(mac);

  foo += "()";
  gROOT->ProcessLine(foo.Data());
}

void LoadMacro(const Text_t* mac)
{
  // Makes sure that macro 'mac' is loaded, but do not reload it.

  TString foo(mac); ChompTail(foo);
  if(gROOT->GetGlobalFunction(foo.Data(), 0, true) == 0)
    gROOT->LoadMacro(mac);
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


}
