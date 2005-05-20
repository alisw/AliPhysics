#ifndef ALIOBJECTMETADATA_H
#define ALIOBJECTMETADATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///
///  Object meta data: it fully describes a run dependent database object. 
///  It is attached to the object to build an AliRunData object
/// 
///

#include <TObject.h>
#include <TString.h>

#include "AliMetaData.h"


class AliObjectMetaData: public AliMetaData {
public:
  AliObjectMetaData();	// default constructor
  AliObjectMetaData
         (const char* name, Int_t firstRun = -1, Int_t lastRun = -1, Int_t period=-1, 
	  const char* objFormat="", const char* responsible="Duck, Donald", 
	  const char* extraInfo="");	// constructor
  virtual ~AliObjectMetaData() {};	// destructor

  AliObjectMetaData(const AliObjectMetaData& entry);	// copy contructor	
  AliObjectMetaData& operator = (const AliObjectMetaData& entry);	// assignment  operator

  void         SetFormat(const char* objFormat) {fFormat = objFormat;} // infos about object's format (array o floats, histos...)
  void         SetResponsible(const char* responsible) {fResponsible = responsible;}	// who made the object?
  void         SetExtraInfo(const char* extraInfo) {fExtraInfo = extraInfo;}	// something else you would like to know
  void         SetPeriod(Int_t period) {fPeriod = period;}	// number of beam period

  const char*  GetFormat() const;
  const char*  GetResponsible() const;
  const char*  GetExtraInfo() const;
  const Int_t  GetPeriod() const;
  
private:
  Int_t        fPeriod;          // beam period 
  TString      fFormat;         // object's format
  TString      fResponsible;    // name of the person responsible for the object
  TString      fExtraInfo;      // extra info about the object
  

ClassDef(AliObjectMetaData, 1)   // Object meta data: full description of a run dependent database object
};

#endif
