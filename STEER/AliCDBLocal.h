#ifndef ALICDBLOCAL_H
#define ALICDBLOCAL_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///
///  access class to a DB file inside an organized directory structure 
///  file name = "DBFolder/detector/dbType/detSpecType/Run#firstRun-#lastRun _v#version.root"
///

#include "AliCDBStorage.h"
#include "AliCDBMetaDataSelect.h"
#include "AliCDBMetaData.h"
#include "AliCDBEntry.h"

class AliCDBLocal: public AliCDBStorage {

public:
//  AliCDBLocal();
  AliCDBLocal(const char* DBFolder = "$(ALICE_ROOT)/DB");
  virtual ~AliCDBLocal();

  TObjArray*   FindDBFiles(const char *name, Int_t runNumber); // return list of files valid for run number
  void	       TagForProduction(const AliCDBMetaDataSelect& selMetaData, UInt_t prodVers); // tag a DB file for production mode
  
protected:
  virtual AliCDBEntry*	GetEntry(AliCDBMetaDataSelect& selMetaData, Int_t runNumber);
  virtual Bool_t        PutEntry(AliCDBEntry* entry);

private:
  AliCDBLocal(const AliCDBLocal& db);
  AliCDBLocal& operator = (const AliCDBLocal& db);
  
  Bool_t     DecodeFileName(const TString strName, int *numArray, TString prefix="_v");	// Gets firstRun, lastRun, version from file name "strName"
  TString    EncodeFileName(int firstRun, int lastRun, int version, TString prefix="_v"); // returns file name from firstRun, lastRun, version	

  TString    fDBFolder;   // the DB folder

ClassDef(AliCDBLocal, 0)      // access class to a DB file in an organized directory structure (DBFolder/detector/dbType/detSpecType)
};

#endif
