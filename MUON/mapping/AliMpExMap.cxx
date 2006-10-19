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
// $MpId: AliMpExMap.cxx,v 1.5 2006/05/24 13:58:29 ivana Exp $
// Category: basic
// ------------------------
// Class AliMpExMap
// ------------------------
// Helper class making Root persistent TExMap
// Author:Ivana Hrivnacova; IPN Orsay

#include "AliMpExMap.h"
#include "AliMpIntPair.h"

#include "AliLog.h"

#include <TClass.h>
#include <TString.h>
#include <Riostream.h>

#include <stdlib.h>

/// \cond CLASSIMP
ClassImp(AliMpExMap)
/// \endcond

//
// static members
//

const Int_t   AliMpExMap::fgkDefaultSize = 300;
const Bool_t  AliMpExMap::fgkDefaultOwnership = true;

const Int_t AliMpExMap::fgkSeparator1 = 10000;
const Int_t AliMpExMap::fgkSeparator2 = 100;

const TString  AliMpExMap::fgkCharacterMap 
  = " 1234567890abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ.";

//
// static methods
//

//______________________________________________________________________________
Long_t  AliMpExMap::GetIndex(const AliMpIntPair& pair)
{
/// Convert the pair of integers to integer.

  if (pair.GetFirst() >= fgkSeparator1 || pair.GetSecond() >= fgkSeparator1) {
    AliFatalClass("Index out of limit.");
    exit(1); 
  }  
      
  return pair.GetFirst()*fgkSeparator1 + pair.GetSecond() + 1;
}  

//_____________________________________________________________________________
Long_t  AliMpExMap::GetIndex(const TString& s)
{
/// Convert the TString to integer.

  if (s.Length() > 5) {
    AliFatalClass("String too long.");
    return 0;
  }  

  Long_t index = 0;
  for (Int_t i=s.Length()-1; i>=0; --i)  
    index = index*fgkSeparator2 + fgkCharacterMap.First(s(i));
  
  return index;
}

//______________________________________________________________________________
AliMpIntPair  AliMpExMap::GetPair(Long_t index)
{
/// Convert the integer index to the pair of integers.

  return AliMpIntPair((index-1)/fgkSeparator1,(index-1)%fgkSeparator1);
}  

//_____________________________________________________________________________
TString  AliMpExMap::GetString(Long_t index)
{
/// Convert the integer index to the string.

  TString s;
  while (index >0) {
    Char_t c = fgkCharacterMap(index%fgkSeparator2);
    s += c;
    index = index/fgkSeparator2;
  }
  return s;
}

//
// constructors/destructor
//

//_____________________________________________________________________________
AliMpExMap::AliMpExMap(Bool_t /*standardConstructor*/) 
  : TObject(),
    fMap(fgkDefaultSize),
    fObjects(fgkDefaultSize),
    fKeys(fgkDefaultSize)
{
/// Standard constructor

  fObjects.SetOwner(fgkDefaultOwnership);
}

//_____________________________________________________________________________
AliMpExMap::AliMpExMap() 
  : TObject(),
    fMap(),
    fObjects(),
    fKeys()
{
/// Default constructor
}

//_____________________________________________________________________________
AliMpExMap::~AliMpExMap() 
{
/// Destructor 
}

//
// private methods
//

//_____________________________________________________________________________
void AliMpExMap::FillMap()
{
/// Fill transient map from the arrays of objects and keys

  for (Int_t i=0; i<fObjects.GetEntriesFast(); i++) 
    fMap.Add(fKeys.At(i), (Long_t)fObjects.At(i)); 
}

//_____________________________________________________________________________
void AliMpExMap::AddKey(Long_t key)
{
/// Add key in array with checking size

  // Resize array if needed
  if (fObjects.GetEntriesFast() == fKeys.GetSize()) {
   fKeys.Set(2*fKeys.GetSize());
   AliWarningStream() << "AliMpExMap::AddKey: resized Key array " << endl;
  } 
   
  fKeys.AddAt(key, fObjects.GetEntriesFast());      
}

//
// public methods
//

//_____________________________________________________________________________
void AliMpExMap::Add(const AliMpIntPair& key, TObject* object)
{
/// Add object with its key to the map and arrays
  
  fMap.Add(GetIndex(key), (Long_t)object);
  AddKey(GetIndex(key));
  fObjects.Add(object);
}

//_____________________________________________________________________________
void AliMpExMap::Add(const TString& key, TObject* object)
{
/// Add object with its key to the map and arrays
  
  fMap.Add(GetIndex(key), (Long_t)object);
  AddKey(GetIndex(key));
  fObjects.Add(object);
}

//_____________________________________________________________________________
void AliMpExMap::Add(Int_t key, TObject* object)
{
/// Add object with its key to the map and arrays
  
  fMap.Add(key, (Long_t)object);
  AddKey(key);
  fObjects.Add(object);
}

//_____________________________________________________________________________
void AliMpExMap::SetSize(Int_t size)
{
/// Set given size to the key array

  // fMap.Set(size);
  // fObjects.Set(size);
  fKeys.Set(size);
} 

//_____________________________________________________________________________
void AliMpExMap::SetOwner(Bool_t owner)
{
/// Set given ownership to object array

  fObjects.SetOwner(owner);
}  

//_____________________________________________________________________________
Int_t AliMpExMap::GetSize() const
{
/// Return the map size

  return fObjects.GetEntriesFast();
}

//_____________________________________________________________________________
TExMapIter AliMpExMap::GetIterator() const
{
/// Return TExMap iterator set to the beginning of the map

  return TExMapIter(&fMap);
}

//_____________________________________________________________________________
TObject* AliMpExMap::GetObject(Int_t index) const
{
/// Return the object via its index in the array
/// (This function makes possible looping over map as over an array)

  if ( index < 0 || index >= fObjects.GetEntriesFast() ) {
    AliErrorStream() << "Index outside limits" << endl;
    return 0;
  }
  
  return fObjects.At(index);
}      

//_____________________________________________________________________________
TObject* AliMpExMap::GetValue(const AliMpIntPair& key) const
{
/// Return the object associated with the given key if found,
/// otherwise return 0

  return reinterpret_cast<TObject*>(fMap.GetValue(GetIndex(key)));
}

//_____________________________________________________________________________
TObject*  AliMpExMap::GetValue(const TString& key) const
{
/// Return the object associated with the given key if found,
/// otherwise return 0

  return reinterpret_cast<TObject*>(fMap.GetValue(GetIndex(key)));
}

//_____________________________________________________________________________
TObject*  AliMpExMap::GetValue(Int_t key) const
{
/// Return the object associated with the given key if found,
/// otherwise return 0

  return reinterpret_cast<TObject*>(fMap.GetValue(key));
}

//_____________________________________________________________________________
void AliMpExMap::Streamer(TBuffer &R__b)
{
/// Customized streamer                                                     \n
/// After the arrays are read, fill the transient map

  if (R__b.IsReading()) {
    AliMpExMap::Class()->ReadBuffer(R__b, this);
    FillMap();
  } 
  else {
    AliMpExMap::Class()->WriteBuffer(R__b, this);
  }
}
