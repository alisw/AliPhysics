#ifndef ALISELECTIONMETADATA_H
#define ALISELECTIONMETADATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///
/// subsample of the object metadata, used to retrieve a stored database object: name="Detector/DBType/DetSpecType", run validity, version
///

#include <TString.h>

#include "AliMetaData.h"


class AliSelectionMetaData: public AliMetaData {
public:
  AliSelectionMetaData();	// default constructor
  AliSelectionMetaData(const char* name, Int_t firstRun = -1, Int_t lastRun = -1, Int_t version = -1);	// constructor
  virtual ~AliSelectionMetaData() {};	// destructor

  AliSelectionMetaData(const AliSelectionMetaData& entry);	// copy constructor
  AliSelectionMetaData(const AliMetaData& entry);	// selection metadata from base class
  AliSelectionMetaData& operator = (const AliSelectionMetaData& entry);	// assignment operator
  
ClassDef(AliSelectionMetaData, 1)   // metadata used to retrieve a DB object
};

#endif
