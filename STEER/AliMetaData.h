#ifndef ALIMETADATA_H
#define ALIMETADATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///
/// meta data of run dependent objects
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

  void                 SetName(const char* name) {fName = name;}
  void                 SetRunRange(Int_t firstRun = -1, Int_t lastRun = 1) 
    {fFirstRun = firstRun; fLastRun = lastRun;}
  void                 SetVersion(Int_t version = -1) {fVersion = version;}

  virtual const char*  GetName() const;
  Int_t                GetFirstRun() const {return fFirstRun;}
  Int_t                GetLastRun() const {return fLastRun;}
  Int_t                GetVersion() const {return fVersion;}

  Bool_t               IsValid(Int_t runNumber, 
			       AliMetaData* metaData = NULL) const;
  virtual Int_t        Compare(const TObject* object) const;
  Bool_t               Matches(const char* name, Int_t runNumber) const;

private:
  TString              fName;           // name of the entry
  Int_t                fFirstRun;       // index of first valid run
  Int_t                fLastRun;        // index of last valid run
  Int_t                fVersion;        // version of the entry

  ClassDef(AliMetaData, 1)   // container for a data base entry object
};

extern Bool_t operator == (const AliMetaData& entry1, 
			   const AliMetaData& entry2);

#endif
