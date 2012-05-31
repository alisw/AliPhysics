// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveMacro.h"

//______________________________________________________________________________
//
// Member fSources is a bitfield, but we do not have a widget
// that can show/edit this (a combo-box with a check-box for each
// entry). So ... use a single value for now,

ClassImp(AliEveMacro)

//______________________________________________________________________________
AliEveMacro::AliEveMacro(Int_t src, const TString& tags, const TString& mac,
			 const TString& foo, const TString& args, Bool_t act) :
  TObject(),
  fSources(src), fTags(tags), fMacro (mac),
  fFunc   (foo), fArgs(args), fActive(act),
  fExecStatus(kNotRun), fExecExcString(), fExecResult(0)
{
  // Constructor.
}

/******************************************************************************/

void AliEveMacro::ResetExecState()
{
  // Reset exec variables into state as if the macro has not been run.

  fExecStatus      = kNotRun;
  fExecExcString = "";
  fExecResult      = 0;
}

void AliEveMacro::SetExecNoData()
{
  // Set last execution state to 'NoData'.

  fExecStatus = kNoData;
}


void AliEveMacro::SetExecOK(TEveElement* result)
{
  // Set last execution state to 'OK' and register result.

  fExecStatus = kOK;
  fExecResult = result;
}


void AliEveMacro::SetExecException(const TString& exception)
{
  // Set last execution state to 'Exception' and store the exception string.

  fExecStatus    = kException;
  fExecExcString = exception;
}

void AliEveMacro::SetExecError()
{
  // Set last execution state to 'Error'.

  fExecStatus = kError;
}

/******************************************************************************/

TString AliEveMacro::FormForExec() const
{
  // Return string suitable for execution.

  return fFunc + "(" + fArgs + ");";
}

TString AliEveMacro::FormForDisplay() const
{
  // Return string suitable for display.

  return TString::Format
    (" %c %-22s  %-30s  %-30s  %-s", fActive ? 'x' : ' ',
     fMacro.Data(), fFunc.Data(), fArgs.Data(), fTags.Data());
}
