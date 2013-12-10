#ifndef ALIMPEXMAPITERATOR_H
#define ALIMPEXMAPITERATOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup core
/// \class AliMpExMapIterator
/// \brief Implementation of TIterator for AliMpExMap
/// 
//  Author Laurent Aphecetche

#ifndef ROOT_TIterator
#  include "TIterator.h"
#endif

class AliMpExMap;
class TString;
class TExMapIter;

//_____________________________________________________________________________
class AliMpExMapIterator : public TIterator
{
  friend class AliMpExMap;

public:
  AliMpExMapIterator(const AliMpExMap& theMap);
  AliMpExMapIterator(const AliMpExMapIterator& rhs);
  AliMpExMapIterator& operator=(const AliMpExMapIterator& rhs);
  AliMpExMapIterator& operator=(const TIterator& rhs);
  
  virtual ~AliMpExMapIterator();
  
  /// The returned object must not be deleted by the user.  

  // Iterating without retrieving a key
  virtual TObject* Next();

  // Iterating with retrieving a key
  TObject*  Next(Int_t& key);
  TObject*  Next(Int_t& keyFirst, Int_t& keySecond);
  TObject*  Next(TString& key);
  
  virtual void Reset(); 
  
  virtual const TCollection* GetCollection() const;

private:
#if (defined(ROOT_VERSION_CODE) && ROOT_VERSION_CODE >= 334082) || \
  (defined(ROOT_SVN_REVISION) && ROOT_SVN_REVISION >= 29598) 
  // ROOT_VERSION_CODE = 334082 corresponds to Root version 5.25/02
  // ROOT_SVN_REVISION is not defined in Root versions after moving in git
  typedef Long64_t Index_t;
#else    
  typedef Long_t Index_t;
#endif    
  Bool_t Next(Index_t& index, TObject*& object);
  
  TExMapIter* fIterator; ///< iterator we are wrapping

  ClassDef(AliMpExMapIterator,0) // TIterator for AliMpExMap
};


#endif
