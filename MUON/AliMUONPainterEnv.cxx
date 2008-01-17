/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// $Id$

#include "AliMUONPainterEnv.h"

#include <TEnv.h>
#include <TSystem.h>

///\class AliMUONPainterEnv
///
/// A resource file handling class.
///
/// Used to get some things persistent between two sessions of the mchview
/// program.
///
///\author Laurent Aphecetche, Subatech

///\cond CLASSIMP
ClassImp(AliMUONPainterEnv)
///\endcond

//_____________________________________________________________________________
AliMUONPainterEnv::AliMUONPainterEnv(const char* resourceFile)
: fEnv(new TEnv(resourceFile))
{
  /// Ctor
}

//_____________________________________________________________________________
AliMUONPainterEnv::~AliMUONPainterEnv()
{
  /// dtor
}

//_____________________________________________________________________________
const char* 
AliMUONPainterEnv::String(const char* resourceName, const char* defaultValue)
{
  /// Retrieve the value associated with a given source, as a string
  
  return fEnv->GetValue(resourceName,defaultValue);
}

//_____________________________________________________________________________
Int_t 
AliMUONPainterEnv::Integer(const char* resourceName, Int_t defaultValue)
{
  /// Retrieve the value associated with a given source, as an integer

  return fEnv->GetValue(resourceName,defaultValue);
}

//_____________________________________________________________________________
Double_t 
AliMUONPainterEnv::Double(const char* resourceName, Double_t defaultValue)
{
  /// Retrieve the value associated with a given source, as a double

  return fEnv->GetValue(resourceName,defaultValue);
}

//_____________________________________________________________________________
void
AliMUONPainterEnv::Save()
{
  /// Save the resource file
  fEnv->WriteFile(gSystem->ExpandPathName(Form("$HOME/%s",fEnv->GetRcName())));
}

//_____________________________________________________________________________
void 
AliMUONPainterEnv::Set(const char* resourceName, Int_t value)
{
  /// Set an integer resource

  fEnv->SetValue(resourceName,Form("%d",value));
}

//_____________________________________________________________________________
void 
AliMUONPainterEnv::Set(const char* resourceName, const char* value)
{
  /// Set a string resource

  fEnv->SetValue(resourceName,value);
}

//_____________________________________________________________________________
void 
AliMUONPainterEnv::Set(const char* resourceName, Double_t value)
{
  /// Set a double resource

  fEnv->SetValue(resourceName,Form("%g",value));
}
