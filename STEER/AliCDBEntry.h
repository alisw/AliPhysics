#ifndef ALICDBENTRY_H
#define ALICDBENTRY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///
/// class that contains an object from the data base and knows about its
/// validity range (meta data)
///

#include <TObject.h>
#include "AliCDBMetaData.h"


class AliCDBEntry: public TObject {
public:
  AliCDBEntry();
  AliCDBEntry(const TObject* object, const AliCDBMetaData& metaData);
  virtual ~AliCDBEntry();

  AliCDBEntry(const AliCDBEntry& entry);
  AliCDBEntry& operator = (const AliCDBEntry& entry);

  void                 SetVersion(Int_t version = -1)
    {fMetaData.SetVersion(version);}

  void                 SetRunRange(Int_t firstRun = -1, Int_t lastRun=-1)
    {fMetaData.SetRunRange(firstRun, lastRun);}

  virtual const char*  GetName() const;
  const TObject*       GetObject() const {return fObject;}
  const AliCDBMetaData&   GetCDBMetaData() const {return fMetaData;}

  virtual Int_t        Compare(const TObject* object) const;

private:
  TObject*             fObject;         // pointer to the data base entry obj.
  AliCDBMetaData    fMetaData;    // object's meta data

  ClassDef(AliCDBEntry, 2)   // container for a data base entry object
};

#endif
