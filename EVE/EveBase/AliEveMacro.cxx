// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveMacro.h"

//______________________________________________________________________________
// Full description of AliEveMacro
//
// !!! Eventually, fSources should be a bitfield, but we need a widget
// that can show/edit this. Like combo-box with a check-box for each
// entry.

ClassImp(AliEveMacro)

//______________________________________________________________________________
AliEveMacro::AliEveMacro(Int_t src, const TString& mac, const TString& foo,
			 const TString& args, Bool_t act) :
  TObject(),
  fSources(src), fMacro(mac), fFunc(foo), fArgs(args), fActive(act)
{
  // Constructor.
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

  TString act(fActive ? " x " : "    ");
  return act + fMacro + " :: " + fFunc + " (" + fArgs + ")";
}
