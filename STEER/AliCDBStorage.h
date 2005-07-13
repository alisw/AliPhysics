#ifndef ALICDBSTORAGE_H
#define ALICDBSTORAGE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///
/// base class for data base access classes
///

#include <TObject.h>
#include <TObjArray.h>
#include "AliCDBMetaDataSelect.h"
#include "AliCDBMetaData.h"


class TFile;
class AliCDBMetaDataSelect;
class AliCDBMetaData;
class AliCDBEntry;

class AliCDBStorage: public TObject {
public:

typedef enum {kDevelopment, kProduction} StorageMode_t;

  virtual ~AliCDBStorage();

  const TObject*         Get(const char* name, Int_t runNumber);	// Gets an object from the database

  Bool_t                 Put(const TObject* object, 
			     const AliCDBMetaData& metaData);	// Put an object into the database

  void                   Select(const AliCDBMetaDataSelect& selMetaData);	// Add a selection criterion 

  Bool_t                 DumpToFile(const char* fileName = "DB.root");	// prepares to dump the retrieved entries to a local file
  
  const AliCDBMetaData& GetCDBMetaData(const char* name);		// Gets the CDBMetaData of the retrieved entry (name=entry's name)
  
  virtual void		 TagForProduction(const AliCDBMetaDataSelect& selMetaData, UInt_t prodVers) = 0;

  static AliCDBStorage* Instance();		// Instance of the current AliCDBStorage object (AliCDBDump, AliCDBLocalFile etc...)

  void 		SetStorageMode(StorageMode_t mode)         {fStorageMode = mode;}
  StorageMode_t GetStorageMode()         {return fStorageMode;}
  

protected:
  AliCDBStorage();

  virtual AliCDBEntry*    GetEntry(AliCDBMetaDataSelect& selMetaData, Int_t runNumber) = 0;	// virtual, see the correspondent method of the derived classes

  virtual Bool_t          PutEntry(AliCDBEntry* entry);			// virtual, see the correspondent method of the derived classes

  AliCDBEntry*            GetCurrentEntry(const char* name) const
    {return (AliCDBEntry*) fEntries.FindObject(name);}
  
  StorageMode_t fStorageMode;

private:
  AliCDBStorage(const AliCDBStorage& db);
  AliCDBStorage& operator = (const AliCDBStorage& db);

  TObjArray              fSelection;   //! meta data selection

  TObjArray              fEntries;     //! array of current AliCDBEntry objects
  TFile*                 fDumpFile;  //! file for dumped entries

  static AliCDBStorage* fgInstance;   //! pointer to the DB instance

  ClassDef(AliCDBStorage, 0)     // base class for data base access classes
};

#endif
