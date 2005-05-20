#ifndef ALIMETADATA_H
#define ALIMETADATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///
/// base class of the meta data of run dependent objects
/// Derived classes: AliObjectMetaData, AliSelectionMetaData
///

#include <TObject.h>
#include <TString.h>


class AliMetaData: public TObject {
public:
  AliMetaData();
  AliMetaData(const char* name, 
	      Int_t firstRun = -1, Int_t lastRun = -1, Int_t version = -1);
  virtual ~AliMetaData() {};

  AliMetaData(const AliMetaData& entry);
  AliMetaData& operator = (const AliMetaData& entry);

  virtual void         SetName(const char* name) {fName = name; DecodeName();}
  /*virtual*/ void         SetRunRange(Int_t firstRun = -1, Int_t lastRun = 1) 
    {fFirstRun = firstRun; fLastRun = lastRun;}
  /*virtual*/ void         SetVersion(Int_t version = -1) {fVersion = version;}
  
  void         SetDetector(const char* Detector) {fDetector = Detector; EncodeName();}
  void         SetDBType(const char* DBType) {fDBType = DBType; EncodeName();}
  void         SetDetSpecType(const char* DetSpecType) {fDetSpecType = DetSpecType; EncodeName();}

  const char*  GetDetector() const;
  const char*  GetDBType() const;
  const char*  GetDetSpecType() const;
  

  /*virtual*/ const char*  GetName() const;
  /*virtual*/ Int_t        GetFirstRun() const {return fFirstRun;}
  /*virtual*/ Int_t        GetLastRun() const {return fLastRun;}
  /*virtual*/ Int_t        GetVersion() const {return fVersion;}

  /*virtual*/ Bool_t       IsValid(Int_t runNumber, 	
			       AliMetaData* metaData = NULL) const; 
  /*virtual*/ Bool_t       IsStrictlyValid(Int_t runNumber, 
			       AliMetaData* metaData = NULL) const;
  /*virtual*/ Int_t        Compare(const TObject* object) const;
  /*virtual*/ Bool_t       Matches(const char* name, Int_t runNumber) const;

protected:
  TString              fName;           // name of the entry
  TString              fDetector;       // name of the detector (ZDC, TPC, etc...)
  TString              fDBType;         // name of the database type (Calib, Align)
  TString              fDetSpecType;    // name of the detector's specific data type (pedestals, gain coeffs...)
  Int_t                fFirstRun;       // index of first valid run
  Int_t                fLastRun;        // index of last valid run
  Int_t                fVersion;        // version of the entry

  void                 EncodeName();
  void                 DecodeName();


  ClassDef(AliMetaData, 2)   // base class of the meta data of run dependent objects
};

extern Bool_t operator == (const AliMetaData& entry1, 
			   const AliMetaData& entry2);

#endif
