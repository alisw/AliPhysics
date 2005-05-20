#ifndef ALIRUNDATAORGANIZEDFILE_H
#define ALIRUNDATAORGANIZEDFILE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///
///  access class to a DB file inside an organized directory structure 
///  (DBFolder/Detector/DBType/DetSpecType)
///

#include "AliRunDataStorage.h"
#include "AliSelectionMetaData.h"
#include "AliObjectMetaData.h"

class AliRunDataOrganizedFile: public AliRunDataStorage {

public:
//  AliRunDataOrganizedFile();
  AliRunDataOrganizedFile(const char* DBFolder = "$(ALICE_ROOT)/DB");
  virtual ~AliRunDataOrganizedFile();

  TObjArray*   FindDataBaseFile(AliSelectionMetaData& selMetaData, Int_t runNumber);
  
protected:
  virtual AliRunData*	GetEntry(AliSelectionMetaData& selMetaData, Int_t runNumber);
  virtual Bool_t        PutEntry(AliRunData* entry);

private:
  AliRunDataOrganizedFile(const AliRunDataOrganizedFile& db);
  AliRunDataOrganizedFile& operator = (const AliRunDataOrganizedFile& db);
  
  void       GetNumbers(const TString strName, int *numArray);	

  TString    fDBFolder;   // the DB folder

ClassDef(AliRunDataOrganizedFile, 0)      // access class to a DB file in an organized directory structure (DBFolder/Detector/DBType/DetSpecType)
};

#endif
