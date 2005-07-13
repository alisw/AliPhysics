#ifndef ALICDBDUMP_H
#define ALICDBDUMP_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///
/// access classes for a data base in a LOCAL file
///

#include "AliCDBStorage.h"
#include "AliCDBMetaDataSelect.h"

class TFile;


class AliCDBDump: public AliCDBStorage {
public:
  AliCDBDump(const char* fileName = "DB.root", Bool_t readOnly = kTRUE);
  virtual ~AliCDBDump();
  void	       TagForProduction(const AliCDBMetaDataSelect& /* selMetaData */, UInt_t /* prodVers */);

protected:
  virtual AliCDBEntry*    GetEntry(AliCDBMetaDataSelect& selMetaData, Int_t runNumber);

  virtual Bool_t         PutEntry(AliCDBEntry* entry);

private:
  AliCDBDump(const AliCDBDump& db);
  AliCDBDump& operator = (const AliCDBDump& db);

  TFile*                 fFile;    //! the DB local file

  ClassDef(AliCDBDump, 0)   // access classes for a data base in a LOCAL file
};

#endif
