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
class AliSelectionMetaData;
class AliObjectMetaData;
class AliRunData;


class AliRunDataStorage: public TObject {
public:
  virtual ~AliRunDataStorage();

  const TObject*         Get(const char* name, Int_t runNumber);	// Gets an object from the database

  Bool_t                 Put(const TObject* object, 
			     const AliObjectMetaData& objMetaData);	// Put an object into the database

  void                   Select(const AliSelectionMetaData& selMetaData);	// Add a selection criterion 

  Bool_t                 RecordToFile(const char* fileName = "DB.root");	// prepares to record the retrieved entries to a local file
  
  const AliObjectMetaData& GetObjectMetaData(const char* name);		// Gets the ObjectMetaData of the retrieved entry (name=entry's name)
  
//  virtual void		 TagForProduction(const AliSelectionMetaData& selMetaData, Uint_t prodVers);

  static AliRunDataStorage* Instance();		// Instance of the current AliRunDataStorage object (AliRunDataFile, AliRunDataOrganizedFile etc...)

protected:
  AliRunDataStorage();

  virtual AliRunData*    GetEntry(AliSelectionMetaData& selMetaData, Int_t runNumber) = 0;	// virtual, see the correspondent method of the derived classes

  virtual Bool_t         PutEntry(AliRunData* entry);			// virtual, see the correspondent method of the derived classes

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
