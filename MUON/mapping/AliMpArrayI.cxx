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
// $MpId: AliMpArrayI.cxx,v 1.5 2006/05/24 13:58:29 ivana Exp $
// Category: basic
// ------------------------
// Class AliMpArrayI
// ------------------------
// Helper class for sorted integer array
// Author:Ivana Hrivnacova; IPN Orsay

#include "AliMpArrayI.h"

#include "AliLog.h"

#include <TClass.h>
#include <TString.h>
#include <Riostream.h>

#include <stdlib.h>
#include <limits.h>

/// \cond CLASSIMP
ClassImp(AliMpArrayI)
/// \endcond

const Int_t AliMpArrayI::fgkDefaultSize = 100;

//_____________________________________________________________________________
AliMpArrayI::AliMpArrayI(Bool_t sort) 
  : TObject(),
    fSort(sort),
    fNofValues(0),
    fValues(fgkDefaultSize),
    fLimits(INT_MIN,INT_MAX)
{
/// Standard & default constructor

}

//_____________________________________________________________________________
AliMpArrayI::AliMpArrayI(TRootIOCtor* /*ioCtor*/) 
  : TObject(),
    fSort(),
    fNofValues(),
    fValues(),
    fLimits()
{
/// IO constructor
}

//_____________________________________________________________________________
AliMpArrayI::~AliMpArrayI() 
{
/// Destructor 
}

//
// private methods
//

//_____________________________________________________________________________
Int_t  AliMpArrayI::GetPosition(Int_t value) const
{
/// Return the new positon where the value should be put

  for ( Int_t i=0; i<fNofValues; i++ ) {
    if ( fValues.At(i) > value ) return i;
  }  

  return fNofValues;
}  

//
// public methods
//

//_____________________________________________________________________________
Bool_t AliMpArrayI::Add(Int_t value)
{
/// Add object with its key to the map and arrays
  
  // Resize array if needed
  if ( fValues.GetSize() == fNofValues ) {
   fValues.Set(2*fValues.GetSize());
   AliWarningStream() << "Resized array." << endl;
  }
  
  
  // The position for the new value  
  Int_t pos;
  if ( fSort ) {
    pos = GetPosition(value);

    // Move elements 
    for ( Int_t i=fNofValues; i>=pos; i-- )
      fValues.AddAt(fValues.At(i), i+1);
  }  
  else
    pos = fNofValues;     
     
  // Add the new value in freed space
  fValues.AddAt(value, pos);  
  ++fNofValues;
  
  // Update linits
  if ( value < fLimits.GetFirst() )  fLimits.SetFirst(value);
  if ( value > fLimits.GetSecond() ) fLimits.SetSecond(value);
  
  return true;
}

//_____________________________________________________________________________
Bool_t AliMpArrayI::Remove(Int_t value)
{
/// Remove value from the array
  
  // Find the position for the new value  
  Int_t pos = GetPosition(value); 
   
  // Return if value is not present
  if ( pos == fNofValues ) return false;

  // Move elements 
  for ( Int_t i=pos; i<fNofValues-1; i++ )
    fValues.AddAt(fValues.At(i+1), i);
    
  // Decrement number of values
  --fNofValues;
  
  return true;
}

//_____________________________________________________________________________
void AliMpArrayI::SetSize(Int_t size)
{
/// Set given size to the array

  fValues.Set(size);
} 

//_____________________________________________________________________________
Int_t AliMpArrayI::GetSize() const
{
/// Return the map size

  return fNofValues;
}

//_____________________________________________________________________________
Int_t AliMpArrayI::GetValue(Int_t index) const
{
/// Return the index-th value 

  if ( index < 0 || index >= fNofValues ) {
    AliErrorStream() << "Index outside limits" << endl;
    return 0;
  }
  
  return fValues.At(index);
}

//_____________________________________________________________________________
Bool_t AliMpArrayI::HasValue(Int_t value) const
{
/// Return true if contains the given value

  if ( ! fNofValues ) return false;

  if ( value < fLimits.GetFirst() || value > fLimits.GetSecond() ) 
    return false;

  for ( Int_t i=0; i<fNofValues; i++ )
    if ( fValues.At(i) == value ) return true;
    
  return false;  
}

