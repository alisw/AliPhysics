// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveMacro_H
#define AliEveMacro_H

#include <TObject.h>
#include <TString.h>

//______________________________________________________________________________
// Short description of AliEveMacro
//

class AliEveMacro : public TObject
{
  friend class AliEveMacroEditor;

public:
  enum DataSource_e { kNone = 0, kRunLoader, kESD, kESDfriend, kRawReader };

  AliEveMacro(Int_t src, const TString& tags, const TString& mac, const TString& foo,
	      const TString& args="", Bool_t act=kTRUE);
  virtual ~AliEveMacro() {}

  Int_t          GetSources() const         { return fSources; }
  void           SetSources(Int_t x)        { fSources = x; }
  const TString& GetTags() const            { return fTags; }
  void           SetTags(const TString& x)  { fTags = x; }
  const TString& GetMacro() const           { return fMacro; }
  void           SetMacro(const TString& x) { fMacro = x; }
  const TString& GetFunc() const            { return fFunc; }
  void           SetFunc(const TString& x)  { fFunc = x; }
  const TString& GetArgs() const            { return fArgs; }
  void           SetArgs(const TString& x)  { fArgs = x; }
  Bool_t         GetActive() const          { return fActive; }
  void           SetActive(Bool_t x)        { fActive = x; }

  TString        FormForExec() const;
  TString        FormForDisplay() const;

protected:
  Int_t   fSources; // Source of data.
  TString fTags;    // Tags describing the macro (for selection).
  TString fMacro;   // Macro where func is defined; if null, assume it is there.
  TString fFunc;    // Function to call.
  TString fArgs;    // Arguments for the function.
  Bool_t  fActive;  // Flag if macro is active.

private:

  ClassDef(AliEveMacro, 0); // Short description.
};

#endif
