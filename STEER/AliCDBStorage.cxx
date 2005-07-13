/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// base class for data base access classes                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <TFile.h>
#include <TKey.h>
#include <TROOT.h>

#include "AliLog.h"
#include "AliCDBEntry.h"
#include "AliCDBStorage.h"


ClassImp(AliCDBStorage)


AliCDBStorage* AliCDBStorage::fgInstance = NULL;


//_____________________________________________________________________________
AliCDBStorage::AliCDBStorage() :
  TObject(),
  fSelection(),
  fEntries(),
  fDumpFile(NULL)
{
// default constructor
  if (fgInstance) delete fgInstance;
  fgInstance = this;
  fStorageMode=kDevelopment;
}

//_____________________________________________________________________________
AliCDBStorage::~AliCDBStorage()
{
// destructor
  fSelection.Delete();
  fEntries.Delete();
  if (fDumpFile) {
    fDumpFile->Close();
    delete fDumpFile;
  }
  fgInstance = NULL;
}

//_____________________________________________________________________________
AliCDBStorage::AliCDBStorage(const AliCDBStorage& db) :
  TObject(db),
  fSelection(),
  fEntries(),
  fDumpFile(NULL)
{
// copy constructor

  AliFatal("not implemented");
}

//_____________________________________________________________________________
AliCDBStorage& AliCDBStorage::operator = (const AliCDBStorage& /*db*/)
{
// assignment operator

  AliFatal("not implemented");
  return *this;
}


//_____________________________________________________________________________
const TObject* AliCDBStorage::Get(const char* name, Int_t runNumber)
{
// get an object from the data base
// (AliCDBStorage is NOT the owner of the returned object)
// name must be in the form "Detector/DBType/DetSpecType"
// es: "ZDC/Calib/Pedestals"

  AliCDBMetaDataSelect defaultMetaData;
  AliCDBMetaDataSelect* selectedMetaData = &defaultMetaData;

  // look for a meta data selection
  for (Int_t i = 0; i < fSelection.GetEntriesFast(); i++) {
    AliCDBMetaDataSelect* selection = (AliCDBMetaDataSelect*) fSelection[i];
    if (!selection) continue;
    if (selection->Matches(name, runNumber)) {
      selectedMetaData = selection;
    }
  }

  // get the entry
  AliCDBMetaDataSelect selMetaData(*selectedMetaData);
  selMetaData.SetName(name);
  AliCDBEntry* entry = GetEntry(selMetaData, runNumber);
  if (entry) {
    AliDebug(2, "got the entry:");
    ToAliDebug(2, entry->Dump());
  } else {
    AliDebug(2, Form("got no entry for %s", name));
  }

  // update array of current entries
  if (!entry) return NULL;  
  TObject* oldEntry = fEntries.FindObject(entry->GetName());
  if (oldEntry) {
    delete fEntries.Remove(oldEntry);
  }
  fEntries.Add(entry);

  // Dump entry to a file (in the same way as AliCDBDump::PutEntry, 
  // so that the file can be opened as a AliCDBDump!)

  if (fDumpFile) {
    fDumpFile->cd();
    TDirectory* saveDir = gDirectory;

    // go to or create the directory
    TString strname(name);
    while (strname.BeginsWith("/")) strname.Remove(0);
    TDirectory* dir = fDumpFile;
    Int_t index = -1;
    while ((index = strname.Index("/")) >= 0) {
      TString dirName(strname(0, index));
      if ((index > 0) && !dir->Get(dirName)) dir->mkdir(dirName);
      dir->cd(dirName);
      dir = gDirectory;
      strname.Remove(0, index+1);
    } 

    entry->Write(strname);
    if (saveDir) saveDir->cd(); else gROOT->cd();
  
  }
  
  return (entry->GetObject())->Clone();

}

//_____________________________________________________________________________
Bool_t AliCDBStorage::Put(const TObject* object, 
			      const AliCDBMetaData& metaData)
{
// put an object into the data base
// (AliCDBStorage does not adopt the object)
// location of where the object is stored is defined by 
// the AliCDBMetaData's name ("Detector/DBType/DetSpecType")
// and run Range. Storage is handled by the PutEntry method
// of the current AliCDBStorage instance. 

  if (!object) return kFALSE;

  AliCDBEntry *entry= new AliCDBEntry(object, metaData);

  Bool_t result = PutEntry(entry);
    
  delete entry;

  return result;
}

//_____________________________________________________________________________
Bool_t AliCDBStorage::PutEntry(AliCDBEntry* entry)
{
// put an object into the data base
// Refer to the specific method of the current AliCDBStorage instance
// (AliCDBDump, AliCDBLocalFile, AliCDBGrid)

  if (!entry) return kFALSE;
  AliError(Form("This is a read only data base. "
		"The object %s was not inserted", entry->GetName()));
  return kFALSE;
}


//_____________________________________________________________________________
void AliCDBStorage::Select(const AliCDBMetaDataSelect& selMetaData)
{
// add some meta data selection criteria

  fSelection.Add(new AliCDBMetaDataSelect(selMetaData));
}


//_____________________________________________________________________________
Bool_t AliCDBStorage::DumpToFile(const char* fileName)
{
// Dump entries retrieved from the data base to a file with the given name

  if (fDumpFile) {
    fDumpFile->Close();
    delete fDumpFile;
  }

  TDirectory* dir = gDirectory;
  fDumpFile = TFile::Open(fileName, "UPDATE");
  if (dir) dir->cd(); else gROOT->cd();
  if (!fDumpFile || !fDumpFile->IsOpen()) {
    AliError(Form("could not open file %s", fileName));
    delete fDumpFile;
    fDumpFile = NULL;
    return kFALSE;
  }
  return kTRUE;
}


//_____________________________________________________________________________
AliCDBStorage* AliCDBStorage::Instance()
{
// return the current instance of the DB (AliCDBDump, AliCDBLocalFile...)
// Example of usage: after creating an istance of AliCDBStorage:
// AliCDBStorage::Instance()->Get(...)

  return fgInstance;
}


//_____________________________________________________________________________
const AliCDBMetaData& AliCDBStorage::GetCDBMetaData(const char* name)
{
// Returns the object's metadata of the already retrieved object
// (useful, for example, if you want to know the format of the object you have 
// retrieved)

AliCDBEntry *entry = (AliCDBEntry*) fEntries.FindObject(name);
 if(!entry){
    AliError(Form("Entry %s not found! You make me crash!",name));
 }
 return entry->GetCDBMetaData(); 

}

