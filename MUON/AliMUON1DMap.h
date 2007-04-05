/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup calib
/// \class AliMUON1DMap
/// \brief Implementation of AliMUONV1DStore
/// 
//  Author Laurent Aphecetche

#ifndef ALIMUON1DMAP_H
#define ALIMUON1DMAP_H

#ifndef ALIMUONV1DSTORE_H
#  include "AliMUONV1DStore.h"
#endif

class AliMpExMap;

class AliMUON1DMap : public AliMUONV1DStore
{
public:
  AliMUON1DMap(Int_t theSize=0);
  AliMUON1DMap(const AliMUON1DMap& other);
  AliMUON1DMap& operator=(const AliMUON1DMap& other);
  
  virtual ~AliMUON1DMap();
  
  /// Return the object stored at i.
  virtual TObject* Get(Int_t i) const;
  
  virtual AliMUONVDataIterator* Iterator() const;
  
  /** Set the object stored at i.
    if replace=false and there's already an object there, returns kFALSE
    */
  virtual Bool_t Set(Int_t i, TObject* object, Bool_t replace);
  
  /// Whether or not this container is the owner of its contents.
  virtual Bool_t IsOwner() const { return kTRUE; }
  
private:
   void CopyTo(AliMUON1DMap& to) const;
  
private:  
    
    AliMpExMap* fMap; ///< Internal array (map)
  
    ClassDef(AliMUON1DMap,1) // Implementation of AliMUONV1DStore
};

#endif
