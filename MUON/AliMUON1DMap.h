/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup calib
/// \class AliMUON1DMap
/// \brief Implementation of AliMUONVStore
/// 
//  Author Laurent Aphecetche

#ifndef ALIMUON1DMAP_H
#define ALIMUON1DMAP_H

#ifndef ALIMUONV1DSTORE_H
#  include "AliMUONVStore.h"
#endif

class AliMpExMap;

class AliMUON1DMap : public AliMUONVStore
{
public:
  AliMUON1DMap(Int_t theSize=0);
  AliMUON1DMap(const AliMUON1DMap& other);
  AliMUON1DMap& operator=(const AliMUON1DMap& other);
  virtual ~AliMUON1DMap();

  virtual Bool_t Add(TObject* object);

  virtual Bool_t CanConnect() const { return kFALSE; }
  
  virtual void Clear(Option_t* opt="");

  virtual AliMUON1DMap* Create() const;
  
  using AliMUONVStore::FindObject;
  
  virtual TObject* FindObject(UInt_t i) const;
  
  virtual TIterator* CreateIterator() const;
  
  using AliMUONVStore::GetSize;
  
  virtual Int_t GetSize() const;
  
private:
   void CopyTo(AliMUON1DMap& to) const;
  /** Set the object stored at i.
    if replace=false and there's already an object there, returns kFALSE
    */
  virtual Bool_t Set(Int_t i, TObject* object, Bool_t replace);
  
private:  
    
    AliMpExMap* fMap; ///< Internal array (map)
  
    ClassDef(AliMUON1DMap,1) // Implementation of AliMUONVStore
};

#endif
