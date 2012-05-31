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
// $MpId: AliMpStringObjMap.cxx,v 1.4 2006/05/24 13:58:29 ivana Exp $

//-----------------------------------------------------------------------------
// Class AliMpStringObjMap
// ------------------------------------ 
// Helper class that substitutes map <string, int> 
// which ALICE does not allow to use 
// Author: Ivana Hrivnacova, IPN Orsay
//-----------------------------------------------------------------------------

#include "AliMpStringObjMap.h"

#include "AliLog.h"

#include <TObjString.h>
#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMpStringObjMap)
/// \endcond

//
// static private methods
//

//______________________________________________________________________________
const TString& AliMpStringObjMap::GetUndefinedKey()
{
  /// Undefined key string
  static const TString kUndefinedKey = "Undefined";
  return kUndefinedKey;
}  

//
// ctors, dtor
//

//______________________________________________________________________________
AliMpStringObjMap::AliMpStringObjMap(Bool_t isOwner)
 : TObject(),
   fNofItems(0),
   fFirstArray(),
   fSecondArray(),
   fCurrentIndex(0)
{
/// Standard constructor

  fFirstArray.SetOwner(true);
  fSecondArray.SetOwner(isOwner);
}

//______________________________________________________________________________
AliMpStringObjMap::~AliMpStringObjMap()
{
/// Destructor

  fFirstArray.Delete();
}

//
// public methods
//

//______________________________________________________________________________
Bool_t  AliMpStringObjMap::Add(const TString& first, TObject* second)
{
/// Add map element if first not yet present

  TObject* second2 = Get(first);
  if ( second2 ) {
    AliError(Form("%s is already present in the map", first.Data()));
    return false;
  }
  
  fFirstArray.Add(new TObjString(first)); 
  fSecondArray.Add(second);
  fNofItems++;
   
  return true;
}  

//______________________________________________________________________________
TObject*  AliMpStringObjMap::Get(const TString& first) const
{
/// Find the element with specified key (first)
  
  for (Int_t i=0; i<fNofItems; i++) {
    if ( ((TObjString*)fFirstArray.At(i))->GetString() == first )
      return fSecondArray.At(i);
  }
  
  return 0;
}      

//______________________________________________________________________________
Int_t  AliMpStringObjMap::GetNofItems() const
{
/// Return the number of elements

  return fNofItems;
}  

//______________________________________________________________________________
void  AliMpStringObjMap::Clear(Option_t* /*option*/)
{
/// Delete the elements

  fNofItems = 0;
  fFirstArray.Delete();
  fSecondArray.Delete();
}  
    
//______________________________________________________________________________
void AliMpStringObjMap::Print(const char* /*option*/) const
{
/// Print the map elements

  for (Int_t i=0; i<fNofItems; i++) {
    cout << setw(4)
         << i << "  "
         << ((TObjString*)fFirstArray.At(i))->GetString()
	 << "  "
	 << setw(5)
	 << fSecondArray.At(i)
	 << endl;
  }
}  	 

//______________________________________________________________________________
void AliMpStringObjMap::Print(const TString& key, ofstream& out) const
{
/// Special printing 

  for (Int_t i=0; i<fNofItems; i++) {
    out  << key << "  "
         << ((TObjString*)fFirstArray.At(i))->GetString()
	 << "  "
	 << setw(5)
	 << fSecondArray.At(i)
	 << endl;
  }
}  	 

//______________________________________________________________________________
void  AliMpStringObjMap::First()
{
/// Set iterator to the first item and return its object

  fCurrentIndex = 0; 
}  
  

//______________________________________________________________________________
void  AliMpStringObjMap::Next()
{
/// Set iterator to the next item and return its object
/// Return 0 if there are no more items

  ++fCurrentIndex;
}  
    

//______________________________________________________________________________
TObject* AliMpStringObjMap::CurrentItem()
{
/// Set iterator to the first item and return its object

  if ( fCurrentIndex >= fNofItems ) return 0;

  return fSecondArray.At(fCurrentIndex);
}  
  

//______________________________________________________________________________
TString AliMpStringObjMap::CurrentKey()
{
/// Set iterator to the first item and return its object

  if ( fCurrentIndex >= fNofItems ) return GetUndefinedKey();
  
  return ((TObjString*)fFirstArray.At(fCurrentIndex))->GetString();
}  
  

//______________________________________________________________________________
Bool_t  AliMpStringObjMap::IsDone() const
{
/// Return true if the iterator reached the end of map

  return fCurrentIndex >= fNofItems;
}   
    
