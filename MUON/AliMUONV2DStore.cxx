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

#include "AliMUONV2DStore.h"

/// \class AliMUONV2DStore
/// Defines an interface equivalent to a 2D array of TObject, indexed
/// by a pair of integers (somehow a matrix, 
/// except that indices are not necessarily sequential).
/// 
/// It's extremely simple and hopefully allow many implementations.
/// It also makes the object ownership self-evident.
///
/// \author Laurent Aphecetche

/// \cond CLASSIMP
ClassImp(AliMUONV2DStore)
/// \endcond

#include "AliLog.h"
#include "AliMUONObjectPair.h"
#include "AliMUONVDataIterator.h"
#include "AliMpIntPair.h"
#include <TMap.h>
#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <Riostream.h>

namespace
{
  //_____________________________________________________________________________
  void Decode(TMap& m, const TString& s, const char* sep)
{
    /// Fills the map m with (key,value) extracted from s
    /// where s is of the form :
    /// key1=value1;key2=value2;key3=value3
    
    TString ss(s);
    ss.ToUpper();
    
    m.SetOwner(true);
    
    TObjArray* a = ss.Tokenize(sep);
    TIter next(a);
    TObjString* o;
    
    while ( ( o = static_cast<TObjString*>(next()) ) )
    {
      TString& os(o->String());
      TObjArray* b = os.Tokenize("=");
      if (b->GetEntries()==2)
      {
        m.Add(b->At(0),b->At(1));
      }
    }
}

//_____________________________________________________________________________
Bool_t FindValue(const TMap& m, const TString& key, TString& value)
{
  /// Find value corresponding to key in map m.
  /// Return false if key not found.
  
  TString skey(key);
  skey.ToUpper();
  value = "";
  TPair* p = static_cast<TPair*>(m.FindObject(skey));
  if (p) 
  {
    value = (static_cast<TObjString*>(p->Value()))->String();
    return kTRUE;
  }
  return kFALSE;
}

}

//_____________________________________________________________________________
AliMUONV2DStore::AliMUONV2DStore()
{
/// Default constructor
}

//_____________________________________________________________________________
AliMUONV2DStore::~AliMUONV2DStore()
{
/// Destructor
}


//_____________________________________________________________________________
void
AliMUONV2DStore::Print(Option_t* opt) const
{
  /// Printout
  /// opt is used to filter which (i,j) couple you want to see
  /// e.g opt="I=12;J=1;opt=Full" to see complete values, but only for the 
  /// (12,1) pair.
  /// Warning : decoding of opt format is not really bullet-proof (yet?)

  AliMUONVDataIterator* it = this->Iterator();
  
  AliMUONObjectPair* pair;
  
  TMap m;
  Decode(m,opt,";");
  
  TString si;  
  Bool_t selectI = FindValue(m,"i",si);
  TString sj;
  Bool_t selectJ = FindValue(m,"j",sj);
  TString sopt;
  FindValue(m,"opt",sopt);
  
  while ( ( pair = static_cast<AliMUONObjectPair*>(it->Next() ) ) )
  {
    AliMpIntPair* ip = static_cast<AliMpIntPair*>(pair->First());
    Int_t i = ip->GetFirst();
    Int_t j = ip->GetSecond();
    if ( selectI && i != si.Atoi() ) continue;
    if ( selectJ && j != sj.Atoi() ) continue;
    cout << Form("[%d,%d]",i,j) << endl;
    TObject* o = pair->Second();
    if (o) 
    {
      o->Print(sopt.Data());
    }
    if ( it->IsOwner() ) delete pair;
  }
  delete it;
}  


