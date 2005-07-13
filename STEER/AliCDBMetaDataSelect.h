#ifndef ALICDBMETADATASELECT_H
#define ALICDBMETADATASELECT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///
/// subsample of the CDB object metadata, used to retrieve a stored database object: name="Detector/DBType/DetSpecType", run validity, version
///

#include <TString.h>

#include "AliCDBMetaData.h"


class AliCDBMetaDataSelect: public AliCDBMetaData {
public:
  AliCDBMetaDataSelect();	// default constructor
  AliCDBMetaDataSelect(const char* name, Int_t firstRun = -1, Int_t lastRun = -1, Int_t version = -1);	// constructor
  virtual ~AliCDBMetaDataSelect() {};	// destructor

  AliCDBMetaDataSelect(const AliCDBMetaDataSelect& entry);	// copy constructor
  AliCDBMetaDataSelect(const AliCDBMetaData& entry);	// selection metadata from base class
  AliCDBMetaDataSelect& operator = (const AliCDBMetaDataSelect& entry);	// assignment operator  
  
 ClassDef(AliCDBMetaDataSelect, 1)   // metadata used to retrieve a DB object
};

#endif
