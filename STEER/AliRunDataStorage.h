#ifndef ALIRUNDATASTORAGE_H
#define ALIRUNDATASTORAGE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///
/// base class for data base access classes
///

#include <TObject.h>
#include <TObjArray.h>

class TFile;
class AliMetaData;
class AliRunData;


class AliRunDataStorage: public TObject {
public:
  virtual ~AliRunDataStorage();

  const TObject*         Get(const char* name, Int_t runNumber);

  Bool_t                 Put(const TObject* object, 
			     const AliMetaData& metaData);

  void                   Select(const AliMetaData& metaData);

  Bool_t                 RecordToFile(const char* fileName = "DB.root");

  static AliRunDataStorage* Instance();

protected:
  AliRunDataStorage();

  virtual AliRunData*    GetEntry(AliMetaData& metaData, Int_t runNumber) = 0;

  virtual Bool_t         PutEntry(AliRunData* entry);

  AliRunData*            GetCurrentEntry(const char* name) const
    {return (AliRunData*) fEntries.FindObject(name);}

private:
  AliRunDataStorage(const AliRunDataStorage& db);
  AliRunDataStorage& operator = (const AliRunDataStorage& db);

  TObjArray              fSelection;   //! meta data selection

  TObjArray              fEntries;     //! array of current AliRunData objects
  TFile*                 fRecordFile;  //! file for recorded entries

  static AliRunDataStorage* fgInstance;   //! pointer to the DB instance

  ClassDef(AliRunDataStorage, 0)     // base class for data base access classes
};

#endif
