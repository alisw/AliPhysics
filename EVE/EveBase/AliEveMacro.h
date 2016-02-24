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

class TEveElement;

//______________________________________________________________________________
// Short description of AliEveMacro
//

class AliEveMacro : public TObject
{
  friend class AliEveMacroEditor;

public:
  enum ExecStatus_e { kNotRun    = -2,
                      kNoData    = -1,
                      kOK        =  0,
                      kException =  1,
                      kError     =  2 };

  AliEveMacro(const TString& tags, const TString& mac, const TString& foo,
	      const TString& args="");
  virtual ~AliEveMacro() {}

  const TString& GetTags() const            { return fTags; }
  void           SetTags(const TString& x)  { fTags = x; }
  const TString& GetMacro() const           { return fMacro; }
  void           SetMacro(const TString& x) { fMacro = x; }
  const TString& GetFunc() const            { return fFunc; }
  void           SetFunc(const TString& x)  { fFunc = x; }
  const TString& GetArgs() const            { return fArgs; }
  void           SetArgs(const TString& x)  { fArgs = x; }

  TString        FormForExec() const;

protected:
  TString fTags;    // Tags describing the macro (for selection).
  TString fMacro;   // Macro where func is defined; if null, assume it is there.
  TString fFunc;    // Function to call.
  TString fArgs;    // Arguments for the function.

private:
  AliEveMacro(const AliEveMacro&);            // Not implemented
  AliEveMacro& operator=(const AliEveMacro&); // Not implemented

  ClassDef(AliEveMacro, 0); // Encapsulation of data reqired for execution of a CINT macro and the result of its last execution.
};

#endif
