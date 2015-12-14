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

//-----------------------------------------------------------------------------
// Class AliMpExMap
// ------------------------
// Helper class making Root persistent TExMap
// Author:Ivana Hrivnacova; IPN Orsay
//-----------------------------------------------------------------------------

#include "AliMpExMap.h"
#include "AliMpExMapIterator.h"

#include "AliLog.h"

#include "TBuffer.h"
#include <TClass.h>
#include <TString.h>
#include <Riostream.h>

#include <stdlib.h>

using std::cout;
using std::endl;
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

//
// static methods
//

//_____________________________________________________________________________
const TString&  AliMpExMap::GetCharacterMap()
{
  /// Return the string mapping characters to integers
  static const TString kCharacterMap 
    = " 1234567890abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ.-";
  return kCharacterMap;  
}

//_____________________________________________________________________________
Long_t  AliMpExMap::GetIndex(const TString& s)
{
/// Convert the TString to integer.

  if (s.Length() > 5) {
    AliErrorClass("String too long.");
    return -1;
  }  

  Long_t index = 0;
  for (Int_t i=s.Length()-1; i>=0; --i)  
    index = index*fgkSeparator2 + GetCharacterMap().First(s(i));
  
  return index;
}

//_____________________________________________________________________________
TString  AliMpExMap::GetString(Long_t index)
{
/// Convert the integer index to the string.

  TString s;
  while (index >0) {
    Char_t c = GetCharacterMap()(index%fgkSeparator2);
    s += c;
    index = index/fgkSeparator2;
  }
  return s;
}

//
// constructors/destructor
//

//_____________________________________________________________________________
AliMpExMap::AliMpExMap() 
  : TObject(),
    fMap(fgkDefaultSize),
    fObjects(fgkDefaultSize),
    fKeys(fgkDefaultSize)
{
      /// Default constructor

  fObjects.SetOwner(fgkDefaultOwnership);
}

//_____________________________________________________________________________
AliMpExMap::AliMpExMap(TRootIOCtor*) 
  : TObject(),
    fMap(),
    fObjects(),
    fKeys()
{
      /// "Root - I/O" constructor
}


//_____________________________________________________________________________
AliMpExMap::AliMpExMap(const AliMpExMap& rhs)
  : TObject(),
    fMap(),
    fObjects(),
    fKeys()

{
  /// Copy ctor
  rhs.Copy(*this);
}

//_____________________________________________________________________________
AliMpExMap&
AliMpExMap::operator=(const AliMpExMap& rhs)
{
  /// Assignment operator

  // check assignment to self
  if (this == &rhs) return *this;

  rhs.Copy(*this);
  return *this;
}

//_____________________________________________________________________________
AliMpExMap::~AliMpExMap() 
{
/// Destructor 
}

//
// private static methods
//

//______________________________________________________________________________
Long_t  AliMpExMap::GetIndex(Int_t first, Int_t second)
{
/// Convert the pair of integers to integer.

  if ( first >= 0xFFFF || second >= 0xFFFF ) 
  {
    AliFatalClass("Index out of limit");
    return 0;
  }
  
  return 1 + ( first | ( second << 16 ) );
           
//  if (pair.GetFirst() >= fgkSeparator1 || pair.GetSecond() >= fgkSeparator1) {
//    AliFatalClass("Index out of limit.");
//    exit(1); 
//  }  
//      
//  return pair.GetFirst()*fgkSeparator1 + pair.GetSecond() + 1;
}  

//______________________________________________________________________________
Int_t  AliMpExMap::GetPairFirst(Long_t index) 
{
/// Return first integer from index (encoded pair)

  return (index-1) & 0xFFFF ;
}  

//______________________________________________________________________________
Int_t  AliMpExMap::GetPairSecond(Long_t index)
{
/// Return second integer from index (encoded pair)

  return ( (index-1) & 0xFFFF0000 ) >> 16 ;
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
   AliDebugStream(1) << "AliMpExMap::AddKey: resized Key array " << endl;
  } 
   
  fKeys.AddAt(key, fObjects.GetEntriesFast());      
}

//_____________________________________________________________________________
void
AliMpExMap::Copy(TObject& dest) const
{
  /// Copy this to dest
  /// Copy implies that dest will become owner of its objects, whatever
  /// the ownership of (*this) is.
  
  AliDebug(1,"");
  
  TObject::Copy(dest);
  AliMpExMap& m = static_cast<AliMpExMap&>(dest);
  m.fKeys = fKeys;
  m.fMap.Delete();
  m.fObjects.Clear();
  
  for ( Int_t i = 0; i <= fObjects.GetLast(); ++i ) 
  {
    TObject* o = fObjects.At(i)->Clone();
    if (!o)
    {
      AliError("Object was not cloned properly ! Please investigate...");
    }
    m.fObjects.AddLast(o);
  }
  m.FillMap();
  m.fObjects.SetOwner(kTRUE);
}

//
// public methods
//

//_____________________________________________________________________________
void AliMpExMap::Clear(Option_t* option)
{
/// Clear memory

  fMap.Delete();
  fObjects.Clear(option);
  fKeys.Reset();
}

//_____________________________________________________________________________
void AliMpExMap::Print(Option_t* opt) const
{
/// Print out

  cout << Form("fMap size/capacity %d/%d",fMap.GetSize(),fMap.Capacity()) 
       << Form(" fObjects.GetSize/Entries %d/%d",fObjects.GetSize(),fObjects.GetEntries()) 
       << Form(" fKeys.GetSize %d",fKeys.GetSize()) << endl;
  
  TString sopt(opt);
  sopt.ToUpper();
  
  if ( sopt.Contains("FULL") ) 
  {
    TIter next(CreateIterator());
    TObject* o;
    while ( ( o = next() ) )
    {
      o->Print();
    }
  }
}

//_____________________________________________________________________________
void AliMpExMap::Add(Int_t keyFirst, Int_t keySecond, TObject* object)
{
/// Add object with its key to the map and arrays
  
  fMap.Add(GetIndex(keyFirst, keySecond), (Long_t)object);
  AddKey(GetIndex(keyFirst, keySecond));
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
Int_t AliMpExMap::GetCapacity() const
{
  /// Return the map capacity
  
  return fObjects.GetSize();
}

//_____________________________________________________________________________
AliMpExMapIterator*
AliMpExMap::CreateIterator() const
{
/// Return iterator set to the beginning of the map

  return new AliMpExMapIterator(*this);
}

//_____________________________________________________________________________
TObject* AliMpExMap::GetValue(Int_t keyFirst, Int_t keySecond) const
{
/// Return the object associated with the given key if found,
/// otherwise return 0

  return reinterpret_cast<TObject*>(fMap.GetValue(GetIndex(keyFirst, keySecond)));
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
// Customized streamer                                                     \n
// After the arrays are read, fill the transient map

  if (R__b.IsReading()) {
    AliMpExMap::Class()->ReadBuffer(R__b, this);
    FillMap();
  } 
  else {
    AliMpExMap::Class()->WriteBuffer(R__b, this);
  }
}
