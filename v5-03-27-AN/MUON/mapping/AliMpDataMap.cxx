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
// Category: basic

//-----------------------------------------------------------------------------
// Class AliMpDataMap
// ------------------------
// TObject class containing a map of strings to strings
// Author:Ivana Hrivnacova; IPN Orsay
//-----------------------------------------------------------------------------

#include "AliMpDataMap.h"

#include "AliLog.h"

#include <TObjString.h>

/// \cond CLASSIMP
ClassImp(AliMpDataMap)
/// \endcond

//_____________________________________________________________________________
AliMpDataMap::AliMpDataMap() 
  : TObject(),
    fMap()
{
/// Standard & default constructor

}

//_____________________________________________________________________________
AliMpDataMap::~AliMpDataMap() 
{
/// Destructor 
}

//
// private methods
//

//_____________________________________________________________________________
void AliMpDataMap::Add(const TString& path, const TString& data)
{
/// Add map element

  fMap.Add(new TObjString(path), new TObjString(data));
}  

//_____________________________________________________________________________
TString  AliMpDataMap::Get(const TString& path, Bool_t warn) const
{
/// Find the data string for given path;
/// give error and return empty string if not found 

  TObject* object = fMap.GetValue(path.Data());

  if ( ! object )  {
    if ( warn ) {
      AliWarningStream()
        << path << " not found in the map." << endl;
    }    
    return "";
  }    
  
  return ((TObjString*)object)->String();
}  
  
