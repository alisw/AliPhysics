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
// access classes for a data base in a LOCAL file                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <TFile.h>
#include <TKey.h>
#include <TROOT.h>
#include "AliLog.h"
#include "AliCDBEntry.h"
#include "AliCDBMetaData.h"
#include "AliCDBDump.h"


ClassImp(AliCDBDump)


//_____________________________________________________________________________
AliCDBDump::AliCDBDump(const char* fileName, Bool_t readOnly) :
  AliCDBStorage(),
  fFile(NULL)
{
// constructor

  if (!fileName) {
    AliError("no file name given");
    return;
  }
  TDirectory* saveDir = gDirectory;
  fFile = TFile::Open(fileName, ((readOnly) ? "READ" : "UPDATE"));
  if (saveDir) saveDir->cd(); else gROOT->cd();
  if (!fFile || !fFile->IsOpen()) {
    AliError(Form("could not open file %s", fileName));
    fFile = NULL;
  }
}

//_____________________________________________________________________________
AliCDBDump::~AliCDBDump()
{
// destructor

  if (fFile) {
    fFile->Close();
    delete fFile;
  }
}

//_____________________________________________________________________________
AliCDBDump::AliCDBDump(const AliCDBDump& /*db*/) :
  AliCDBStorage(),
  fFile(NULL)
{
// copy constructor

  AliFatal("not implemented");
}

//_____________________________________________________________________________
AliCDBDump& AliCDBDump::operator = (const AliCDBDump& /*db*/)
{
// assignment operator

  AliFatal("not implemented");
  return *this;
}


//_____________________________________________________________________________
AliCDBEntry* AliCDBDump::GetEntry(AliCDBMetaDataSelect& selMetaData, Int_t runNumber)
{
// get an object from the data base

  // go to the directory
  TDirectory* saveDir = gDirectory;
  TDirectory *dir = fFile;
  TString name(selMetaData.GetName());
  Int_t last = name.Last('/');
  if (last < 0) {
    fFile->cd();
  } else {
    TString dirName(name(0, last));
      if (!dir->cd(dirName)) {
      AliError(Form("no directory %s found", dirName.Data()));
      if (saveDir) saveDir->cd(); else gROOT->cd();
      return NULL;
    }
    name.Remove(0, last+1);
  }

  dir = gDirectory;
  TKey* key = dir->GetKey(name); 
  if (!key) {
    AliError(Form("no object with name %s found", selMetaData.GetName()));
    if (saveDir) saveDir->cd(); else gROOT->cd();
    return NULL;
  }
  Int_t nCycles = key->GetCycle();

  // find the closest entry
  AliCDBEntry* closestEntry = NULL;
  for (Int_t iCycle = nCycles; iCycle > 0; iCycle--) {
    key = dir->GetKey(name, iCycle);
    
    if (!key) continue;
    AliCDBEntry* entry = (AliCDBEntry*) key->ReadObj();
    if (!entry) continue;
    if (!entry->InheritsFrom(AliCDBEntry::Class())) {
      AliCDBMetaData metaData;
      entry = new AliCDBEntry(entry, metaData);
    }
    if (!entry->GetCDBMetaData().IsStrictlyValid(runNumber, &selMetaData) ||
	(entry->Compare(closestEntry) <= 0)) {
      delete entry;
      continue;
    }
    delete closestEntry;
    closestEntry = entry;
  }
  if (saveDir) saveDir->cd(); else gROOT->cd();
  if(!closestEntry) AliError(Form("No valid entry found for: name %s, version %d, run %d!!!",
            selMetaData.GetName(),selMetaData.GetVersion(),runNumber));
  if (!closestEntry) return NULL;

//   if(selMetaData.GetVersion() > -1 && (closestEntry->GetCDBMetaData()).GetVersion() != selMetaData.GetVersion()) 
//      AliWarning(Form("Warning: selected version (%d) not found, got version %d instead",
//             selMetaData.GetVersion(),(closestEntry->GetCDBMetaData()).GetVersion()));

  return closestEntry;
}

//_____________________________________________________________________________
Bool_t AliCDBDump::PutEntry(AliCDBEntry* entry)
{
// puts an object into the database

// AliCDBEntry entry is composed by the object and its MetaData
// this method takes the metaData, reads the name, runRange and Version
// creates the TDirectory structure into the file
// looks for runs with same name, if exist increment version
// (therefore version should not be put in the metadata)
// Note: the key name of the entry is "DetSpecType"
// return result 

  if (!entry || !fFile) return kFALSE;
  if (!fFile->IsWritable()) {
    AliError(Form("The data base file was opened in read only mode. "
		  "The object %s was not inserted", entry->GetName()));
    return kFALSE;
  }
  
  fFile->cd();
  TDirectory* saveDir = gDirectory;

  // go to or create the directory
  TString name(entry->GetName());
  while (name.BeginsWith("/")) name.Remove(0);
  TDirectory* dir = fFile;
  Int_t index = -1;
  while ((index = name.Index("/")) >= 0) {
    TString dirName(name(0, index));
    if ((index > 0) && !dir->Get(dirName)) dir->mkdir(dirName);
    dir->cd(dirName);
    dir = gDirectory;
    name.Remove(0, index+1);
  } 

  // determine the version number
  Int_t version = 0;
  TKey* key = dir->GetKey(name); 
  if (key) {
    Int_t nCycles = key->GetCycle();
    for (Int_t iCycle = nCycles; iCycle > 0; iCycle--) {
      key = dir->GetKey(name, iCycle); 
      if (!key) continue;
      AliCDBEntry* oldEntry = (AliCDBEntry*) key->ReadObj();
      if (!oldEntry) continue;
      if (oldEntry->InheritsFrom(AliCDBEntry::Class())) {
	if (version <= oldEntry->GetCDBMetaData().GetVersion()) {
	  version = oldEntry->GetCDBMetaData().GetVersion()+1;
	}
      }
      delete oldEntry;
    }
  }
  entry->SetVersion(version);

  Bool_t result = (entry->Write(name) != 0);
  if (saveDir) saveDir->cd(); else gROOT->cd();
  
  if(result) {
    AliInfo(Form("Run object %s",entry->GetName()));
    AliInfo(Form("was successfully written into file %s",fFile->GetName()));
  }  
  
  return result;
}

//_____________________________________________________________________________
void AliCDBDump::TagForProduction
		(const AliCDBMetaDataSelect& /* selMetaData */, UInt_t /* prodVers */){

  AliError(Form("Not implemented in this case"));
}
