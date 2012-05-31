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

#include <TROOT.h>
#include <TFile.h>
#include <TObjString.h>
#include <TTree.h>

#include "AliLog.h"
#include "AliLoader.h"
#include "AliRun.h"
#include "AliStream.h"

////////////////////////////////////////////////////////////////////////
//
// AliStream.cxx
//
// - store file names associated with a given stream
// - open and close files
// - return serial event number of the next event in the stream
// and the TFile pointer for a proper file
//  Author: Jiri Chudoba (CERN), 2001
//
////////////////////////////////////////////////////////////////////////

ClassImp(AliStream)

AliStream::AliStream():
  fLastEventSerialNr(-1),
  fLastEventNr(0),
  fCurrentFileIndex(-1),
  fEvents(0),
  fMode(0),
  fFileNames(0x0),
  fEventFolderName(0)
{
  // root requires default ctor, where no new objects can be created
  // do not use this ctor, it is supplied only for root needs
}
//_______________________________________________________________________

AliStream::AliStream(const char* foldername,Option_t *option):
  fLastEventSerialNr(-1),
  fLastEventNr(0),
  fCurrentFileIndex(-1),
  fEvents(0),
  fMode(option),
  fFileNames(new TObjArray(1)),
  fEventFolderName(foldername)
{
// ctor
}
//_______________________________________________________________________

AliStream::AliStream(const AliStream &as):
  TNamed(as),
  fLastEventSerialNr(-1),
  fLastEventNr(0),
  fCurrentFileIndex(-1),
  fEvents(0),
  fMode(0),
  fFileNames(0x0),
  fEventFolderName(" ")
{
  //
  // Copy ctor
  //
  as.Copy(*this);
}
//_______________________________________________________________________

AliStream::~AliStream()
{
// default dtor
  delete AliRunLoader::GetRunLoader(fEventFolderName); //clear the eventuall session
  if (fFileNames) delete fFileNames;
}
//_______________________________________________________________________

void AliStream::Copy(TObject &) const
{
  //
  // Copy function
  //
  AliFatal("Not implemented!");
}
//_______________________________________________________________________

void AliStream::AddFile(const char *fileName)
{
// stores the name of the file
  TObjString *name = new TObjString(fileName);
  fFileNames->Add(name);
}
//_______________________________________________________________________

Bool_t AliStream::NextEventInStream()
{
// returns kFALSE if no more events
// returns kTRUE and the serial nr of the next event
// fCurrentFile points to the file containing offered event

// no files given:
  if (fFileNames->GetLast() < 0) return kFALSE;
  
  AliRunLoader* currentloader = AliRunLoader::GetRunLoader(fEventFolderName);
  if (currentloader == 0x0) 
   {
    AliDebug(1, Form(
         "Can not get RL from folder named %s. Attempting to open next file",
         fEventFolderName.Data()));
    Int_t res = OpenNextFile();
    if ( res == 0) return kFALSE;
    currentloader = AliRunLoader::GetRunLoader(fEventFolderName);
   }
  
  if (fLastEventSerialNr+1 >= fEvents) 
   {
    if (!OpenNextFile()) return kFALSE;
   }
  AliDebug(1, Form("Trying to get event %d",fLastEventSerialNr+1));
  currentloader->GetEvent(++fLastEventSerialNr);
  return kTRUE;
}
//_______________________________________________________________________

void AliStream::ChangeMode(Option_t* option)
{
  // set the mode to READ or UPDATE, reopen file with the new mode
  // only change from UPDATE to READ have sense in the current scheme,
  // other changes are possible but not usefull

  fMode = option;
  AliRunLoader* currentloader = AliRunLoader::GetRunLoader(fEventFolderName);
  if (currentloader) {
    delete currentloader;
    fCurrentFileIndex--;
    OpenNextFile();
  }
}
//_______________________________________________________________________

Bool_t AliStream::OpenNextFile()
{
  //
  // Opens next file in the list
  //
  if (++fCurrentFileIndex > fFileNames->GetLast()) {
    AliInfo("No more files in the stream") ;
    return kFALSE;
  }

  const char* filename =   static_cast<TObjString*>(fFileNames->At(fCurrentFileIndex))->GetName();

// check if the file was already opened by some other code
  TFile *f = (TFile *)(gROOT->GetListOfFiles()->FindObject(filename));
  if (f) f->Close();

  AliRunLoader* currentloader = AliRunLoader::GetRunLoader(fEventFolderName);
  
  if (currentloader) 
   {
     delete currentloader;
   }
  
  currentloader = AliRunLoader::Open(filename,fEventFolderName,fMode);
  

  if (currentloader == 0x0) 
   {
// cannot open file specified on input. Do not skip it silently.
    AliError("Cannot open session ");
    return kFALSE;
   }
   
// find nr of events in the given file  
  
  if ( AliLoader::TestFileOption(fMode) )//tests if file is opened in read or update mode
   {
    Int_t res = currentloader->LoadHeader();
    if (res)
     {
       AliError("Problems with loading header");
       return kFALSE;
     }
    fEvents = static_cast<Int_t>(currentloader->TreeE()->GetEntries());
   }
  else
    {
     //if it is new, create or recreate there is no chance to find header in file
      fEvents = 0;
    }
   
  fLastEventSerialNr = -1;
  return kTRUE;
}
//_______________________________________________________________________

Bool_t AliStream::ImportgAlice()
{
  //
  // Imports gAlice object from file
  //
  if (fFileNames->GetLast() < 0) return kFALSE;
  
  AliRunLoader* currentloader = AliRunLoader::GetRunLoader(fEventFolderName);
  if (!currentloader) 
   {
    if (!OpenNextFile()) return kFALSE;
    currentloader = AliRunLoader::GetRunLoader(fEventFolderName);
   }
  currentloader->LoadgAlice();
  gAlice = currentloader->GetAliRun();
  if (!gAlice)  return kFALSE;
  return kTRUE;
}

//_______________________________________________________________________
TString AliStream::GetFileName(Int_t order) const
{
  // returns name of the order-th file
  // returns empty string if such file does not exist
  // first file in the input stream is 0
  TString fileName("");
  if (order > fFileNames->GetLast()) return fileName;
  TObjString *fileNameStored = dynamic_cast<TObjString*>(fFileNames->At(order));
  if (fileNameStored) fileName = fileNameStored->GetString();
  return fileName;
}

