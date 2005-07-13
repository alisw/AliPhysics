#ifndef ALICDBMETADATA_H
#define ALICDBMETADATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///
///  CDB meta data: it fully describes a run dependent database object. 
///  It is attached to the object to build an AliCDBEntry object
/// 
///

#include <TObject.h>
#include <TString.h>

class AliCDBMetaData: public TObject {
public:
  AliCDBMetaData();	// default constructor
  AliCDBMetaData
         (const char* name, Int_t firstRun = -1, Int_t lastRun = -1, Int_t period=-1, 
	  const char* objFormat="", const char* responsible="Duck, Donald", 
	  const char* extraInfo="");	// constructor

  virtual ~AliCDBMetaData() {};	// destructor

  AliCDBMetaData(const AliCDBMetaData& entry);	// copy contructor	
  AliCDBMetaData& operator = (const AliCDBMetaData& entry);	// assignment  operator


  void  	SetName(const char* name) {fName = name; DecodeName();} // object's name ("Detector/DBType/DetSpecType")
  void 		SetDetector(const char* Detector) {fDetector = Detector; EncodeName();} // set the detector's name (ZDC,ITS ...)
  void 		SetDBType(const char* DBType) {fDBType = DBType; EncodeName();} // set the database type (Calib, Align ...)
  void 		SetDetSpecType(const char* DetSpecType) {fDetSpecType = DetSpecType; EncodeName();} // set the detector's specific type name (Pedestals, GainConst, DeadChannelMaps...)
  void 		SetRunRange(Int_t firstRun = -1, Int_t lastRun = 1) 
    			{fFirstRun = firstRun; fLastRun = lastRun;} // set the run validity range
  void 		SetVersion(Int_t version = -1) {fVersion = version;} // object's version (automatically set during storage!)
  void 		SetPeriod(Int_t period) {fPeriod = period;}	// set number of beam period
  void 		SetFormat(const char* objFormat) {fFormat = objFormat;} // set infos about object's format (array o floats, histos...)
  void 		SetResponsible(const char* responsible) {fResponsible = responsible;}	// who made the object?
  void 		SetExtraInfo(const char* extraInfo) {fExtraInfo = extraInfo;}	// something else you would like to know
  

  const char* 	GetName() const {return fName.Data();} // get the name 
  const char* 	GetDetector() const {return fDetector.Data();} // get the detector's name 
  const char* 	GetDBType() const {return fDBType.Data();} // get the database type 
  const char* 	GetDetSpecType() const {return fDetSpecType.Data();} // get the detector's specific type name 
  
  const Int_t 	GetFirstRun() const {return fFirstRun;} // get the first valid run
  const Int_t 	GetLastRun() const {return fLastRun;} // get the last valid run
  const Int_t 	GetVersion() const {return fVersion;} // get the version
  const Int_t 	GetPeriod() const {return fPeriod;} // get the beam period
  const char* 	GetFormat() const {return fFormat.Data();} // get the object's format
  const char* 	GetResponsible() const  {return fResponsible.Data();} // get the responsible's name
  const char* 	GetExtraInfo() const {return fExtraInfo.Data();} // get the object's extra info
  
  Bool_t 	IsValid(Int_t runNumber,     
     	           	AliCDBMetaData* metaData = NULL) const; 
  Bool_t 	IsStrictlyValid(Int_t runNumber, 
     	           	AliCDBMetaData* metaData = NULL) const;
  Int_t		Compare(const TObject* object) const;
  Bool_t	Matches(const char* name, Int_t runNumber) const;

protected:
  TString 	fName;           // name of the entry
  TString 	fDetector;       // name of the detector (ZDC, TPC, etc...)
  TString 	fDBType;         // name of the database type (Calib, Align)
  TString 	fDetSpecType;    // name of the detector's specific data type (pedestals, gain coeffs...)
  Int_t  	fFirstRun;       // index of first valid run
  Int_t 	fLastRun;        // index of last valid run
  Int_t 	fVersion;        // version of the entry
  Int_t 	fPeriod;         // beam period 
  TString 	fFormat;         // object's format
  TString 	fResponsible;    // name of the person responsible for the object
  TString 	fExtraInfo;      // extra info about the object
  
  void 		EncodeName();
  void 		DecodeName();

ClassDef(AliCDBMetaData, 1)   // CDB meta data: full description of a run dependent database object
};

extern Bool_t operator == (const AliCDBMetaData& entry1, 
			   const AliCDBMetaData& entry2);

#endif
