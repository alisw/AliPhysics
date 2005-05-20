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
#include "AliObjectMetaData.h"


class AliRunData: public TObject {
public:
  AliRunData();
  AliRunData(TObject* object, const AliObjectMetaData& objMetaData);
  virtual ~AliRunData();

  AliRunData(const AliRunData& entry);
  AliRunData& operator = (const AliRunData& entry);

  void                 SetVersion(Int_t version = -1)
    {fObjMetaData.SetVersion(version);}

  void                 SetRunRange(Int_t firstRun = -1, Int_t lastRun=-1)
    {fObjMetaData.SetRunRange(firstRun, lastRun);}

  virtual const char*  GetName() const;
  const TObject*       GetObject() const {return fObject;}
  const AliObjectMetaData&   GetObjectMetaData() const {return fObjMetaData;}

  virtual Int_t        Compare(const TObject* object) const;

private:
  TObject*             fObject;         // pointer to the data base entry obj.
  AliObjectMetaData    fObjMetaData;    // object's meta data

  ClassDef(AliRunData, 2)   // container for a data base entry object
};

#endif
