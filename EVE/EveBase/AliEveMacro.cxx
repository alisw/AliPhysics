// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveMacro.h"

ClassImp(AliEveMacro)

AliEveMacro::AliEveMacro(const TString& tags, const TString& mac,
			 const TString& foo, const TString& args) :
  TObject(),fTags(tags), fMacro (mac),
  fFunc   (foo), fArgs(args)
{
  // Constructor.
}

TString AliEveMacro::FormForExec() const
{
  // Return string suitable for execution.
  return fFunc + "(" + fArgs + ");";
}

