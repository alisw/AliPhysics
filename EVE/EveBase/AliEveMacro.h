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
  enum DataSource_e { kNone      = 0,
                      kRunLoader = 1,
                      kESD       = 2,
                      kESDfriend = 4,
                      kRawReader = 8,
                      kAOD       = 16 };

  enum ExecStatus_e { kNotRun    = -2,
                      kNoData    = -1,
                      kOK        =  0,
                      kException =  1,
                      kError     =  2 };

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

  Bool_t         RequiresRunLoader() const { return fSources & kRunLoader; }
  Bool_t         RequiresESD()       const { return fSources & kESD;       }
  Bool_t         RequiresESDfriend() const { return fSources & kESDfriend; }
  Bool_t         RequiresRawReader() const { return fSources & kRawReader; }
  Bool_t         RequiresAOD()       const { return fSources & kAOD;       }

  void           ResetExecState();

  void           SetExecNoData();
  void           SetExecOK(TEveElement* result);
  void           SetExecException(const TString& exception);
  void           SetExecError();

  ExecStatus_e   GetExecStatus()    const { return fExecStatus;    }
  const TString& GetExecException() const { return fExecExcString; }
  TEveElement*   GetExecResult()    const { return fExecResult;    }

  Bool_t         WasExecTried()     const { return fExecStatus >= kOK; }
  Bool_t         WasExecOK()        const { return fExecStatus == kOK; }

  TString        FormForExec() const;
  TString        FormForDisplay() const;

protected:
  Int_t   fSources; // Source of data, bitwise or of DataSource_e entries.
  TString fTags;    // Tags describing the macro (for selection).
  TString fMacro;   // Macro where func is defined; if null, assume it is there.
  TString fFunc;    // Function to call.
  TString fArgs;    // Arguments for the function.
  Bool_t  fActive;  // Flag if macro is active.

  ExecStatus_e  fExecStatus;
  TString       fExecExcString;
  TEveElement  *fExecResult;

private:
  AliEveMacro(const AliEveMacro&);            // Not implemented
  AliEveMacro& operator=(const AliEveMacro&); // Not implemented

  ClassDef(AliEveMacro, 0); // Encapsulation of data reqired for execution of a CINT macro and the result of its last execution.
};

#endif
