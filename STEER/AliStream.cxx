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

/*
$Log$
Revision 1.8  2002/10/22 15:02:15  alibrary
Introducing Riostream.h

Revision 1.7  2002/10/14 14:57:33  hristov
Merging the VirtualMC branch to the main development branch (HEAD)

Revision 1.4.6.2  2002/07/24 10:08:13  alibrary
Updating VirtualMC

Revision 1.6  2002/07/16 13:48:39  jchudoba
Add methods to get access to names of files used in merging. Correct memory leak in dtor (thanks to Yves Schutz.)

Revision 1.5  2002/04/09 13:38:47  jchudoba
Add const to the filename argument

Revision 1.4  2001/12/03 07:10:13  jchudoba
Default ctor cannot create new objects, create dummy default ctor which leaves object in not well defined state - to be used only by root for I/O

Revision 1.3  2001/10/15 17:31:56  jchudoba
Bug correction

Revision 1.2  2001/10/04 15:58:52  jchudoba
Option to open the stream in READ or UPDATE mode

Revision 1.1  2001/09/19 06:20:50  jchudoba
Class to manage input filenames, used by AliRunDigitizer

*/

////////////////////////////////////////////////////////////////////////
//
// AliStream.cxx
//
// - store file names associated with a given stream
// - open and close files
// - return serial event number of the next event in the stream
// and the TFile pointer for a proper file
//
////////////////////////////////////////////////////////////////////////

#include <Riostream.h>

#include "TFile.h"
#include "TObjString.h"
#include "TROOT.h"
#include "TTree.h"

#include "AliRun.h"
#include "AliStream.h"

ClassImp(AliStream)

//_______________________________________________________________________
AliStream::AliStream():
  fLastEventSerialNr(0),
  fLastEventNr(0),
  fCurrentFileIndex(0),
  fEvents(0),
  fMode(0),
  fCurrentFile(0),
  fFileNames(0)
{
  //
  // root requires default ctor, where no new objects can be created
  // do not use this ctor, it is supplied only for root needs
  //
}

//_______________________________________________________________________
AliStream::AliStream(Option_t *option):
  fLastEventSerialNr(-1),
  fLastEventNr(0),
  fCurrentFileIndex(-1),
  fEvents(0),
  fMode(option),
  fCurrentFile(0),
  fFileNames(new TObjArray(1))
{
  //
  // Default ctor
  //
}

//_______________________________________________________________________
AliStream::AliStream(const AliStream& str):
  TNamed(str),
  fLastEventSerialNr(0),
  fLastEventNr(0),
  fCurrentFileIndex(0),
  fEvents(0),
  fMode(0),
  fCurrentFile(0),
  fFileNames(0)
{
  //
  // Copy ctor
  //
  str.Copy(*this);
}

//_______________________________________________________________________
void AliStream::Copy(AliStream& ) const
{
  Fatal("Copy","Not implemented\n");
}

//_______________________________________________________________________
AliStream::~AliStream()
{
// default dtor
  if (fFileNames) {
    fFileNames->Delete();
    delete fFileNames;
  }
}

//_______________________________________________________________________
void AliStream::AddFile(const char *fileName)
{
// stores the name of the file
  TObjString *name = new TObjString(fileName);
  fFileNames->Add(name);
}

//_______________________________________________________________________
Bool_t AliStream::NextEventInStream(Int_t &serialNr)
{
// returns kFALSE if no more events
// returns kTRUE and the serial nr of the next event
// fCurrentFile points to the file containing offered event

// no files given:
  if (fFileNames->GetLast() < 0) return kFALSE;

  if (!fCurrentFile) {
    if (!OpenNextFile()) return kFALSE;
  }
  
  if (fLastEventSerialNr+1 >= fEvents) {
    if (!OpenNextFile()) return kFALSE;
  }
  serialNr = ++fLastEventSerialNr;
  return kTRUE;
}

//_______________________________________________________________________
void AliStream::ChangeMode(Option_t* option)
// set the mode to READ or UPDATE, reopen file with the new mode
// only change from UPDATE to READ have sense in the current scheme,
// other changes are possible but not usefull
{
  fMode = option;
  if (fCurrentFile) {
    fCurrentFile->Close();
    fCurrentFileIndex--;
    OpenNextFile();
  }
}

//_______________________________________________________________________
Bool_t AliStream::OpenNextFile()
{
  if (++fCurrentFileIndex > fFileNames->GetLast()) {
    cerr<<"No more files in the stream"<<endl;
    return kFALSE;
  }

  const char * filename = 
    static_cast<TObjString*>(fFileNames->At(fCurrentFileIndex))->GetName();

// check if the file was already opened by some other code
  TFile *f = dynamic_cast<TFile*>(gROOT->GetListOfFiles()->FindObject(filename));
  if (f) f->Close();

  if (fCurrentFile) {
    if (fCurrentFile->IsOpen()) {
      fCurrentFile->Close();
    }
  }

  fCurrentFile = TFile::Open(filename,fMode.Data());
  if (!fCurrentFile) {
// cannot open file specified on input. Do not skip it silently.
    cerr<<"Cannot open file "<<filename<<endl;
    return kFALSE;
  }
// find nr of events in the given file  
  TTree * te = dynamic_cast<TTree *>(fCurrentFile->Get("TE"));
  if (!te) {
    Error("OpenNextFile", "input file does not contain TE");
    return kFALSE;
  }
  fEvents = static_cast<Int_t>(te->GetEntries());
  fLastEventSerialNr = -1;
  return kTRUE;
}

//_______________________________________________________________________
Bool_t AliStream::ImportgAlice()
{
  if (fFileNames->GetLast() < 0) return kFALSE;
  if (!fCurrentFile) {
    if (!OpenNextFile()) return kFALSE;
  }
  gAlice = dynamic_cast<AliRun*>(fCurrentFile->Get("gAlice"));
  if (!gAlice)  return kFALSE;
  return kTRUE;
}

//_______________________________________________________________________
TString AliStream::GetFileName(const Int_t order) const
// returns name of the order-th file
// returns empty string if such file does not exist
// first file in the input stream is 0
{
  TString fileName("");
  if (order > fFileNames->GetLast()) return fileName;
  TObjString *fileNameStored = dynamic_cast<TObjString*>(fFileNames->At(order));
  if (fileNameStored) fileName = fileNameStored->GetString();
  return fileName;
}

