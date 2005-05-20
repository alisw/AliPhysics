#ifndef ALIRUNDATAFILE_H
#define ALIRUNDATAFILE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///
/// access classes for a data base in a LOCAL file
///

#include "AliRunDataStorage.h"
#include "AliSelectionMetaData.h"

class TFile;


class AliRunDataFile: public AliRunDataStorage {
public:
  AliRunDataFile(const char* fileName = "DB.root", Bool_t readOnly = kTRUE);
  virtual ~AliRunDataFile();

protected:
  virtual AliRunData*    GetEntry(AliSelectionMetaData& selMetaData, Int_t runNumber);

  virtual Bool_t         PutEntry(AliRunData* entry);

private:
  AliRunDataFile(const AliRunDataFile& db);
  AliRunDataFile& operator = (const AliRunDataFile& db);

  TFile*                 fFile;    //! the DB local file

  ClassDef(AliRunDataFile, 0)   // access classes for a data base in a LOCAL file
};

#endif
