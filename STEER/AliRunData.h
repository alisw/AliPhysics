#ifndef ALIRUNDATA_H
#define ALIRUNDATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///
/// class that contains an object from the data base and knows about its
/// validity range (meta data)
///

#include <TObject.h>
#include "AliMetaData.h"


class AliRunData: public TObject {
public:
  AliRunData();
  AliRunData(TObject* object, const AliMetaData& metaData);
  virtual ~AliRunData();

  AliRunData(const AliRunData& entry);
  AliRunData& operator = (const AliRunData& entry);

  void                 SetVersion(Int_t version = -1)
    {fMetaData.SetVersion(version);}

  virtual const char*  GetName() const;
  const TObject*       GetObject() const {return fObject;}
  const AliMetaData&   GetMetaData() const {return fMetaData;}

  virtual Int_t        Compare(const TObject* object) const;

private:
  TObject*             fObject;         // pointer to the data base entry obj.
  AliMetaData          fMetaData;       // meta data

  ClassDef(AliRunData, 1)   // container for a data base entry object
};

#endif
