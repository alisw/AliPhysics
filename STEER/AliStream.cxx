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

#include <iostream.h>

#include "TTree.h"

#include "AliStream.h"

#include "AliRun.h"

ClassImp(AliStream)

AliStream::AliStream()
{
// default ctor
  fLastEventSerialNr = -1;
  fLastEventNr = 0;
  fCurrentFileIndex = -1;
  fCurrentFile = 0;
  fEvents = 0;
  fFileNames = new TObjArray(1);
}

////////////////////////////////////////////////////////////////////////
AliStream::~AliStream()
{
// default dtor
  if (fFileNames) delete fFileNames;
}

////////////////////////////////////////////////////////////////////////
void AliStream::AddFile(char *fileName)
{
// stores the name of the file
  TObjString *name = new TObjString(fileName);
  fFileNames->Add(name);
}

////////////////////////////////////////////////////////////////////////
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

////////////////////////////////////////////////////////////////////////
Bool_t AliStream::OpenNextFile()
{
  if (++fCurrentFileIndex > fFileNames->GetLast()) {
    cerr<<"No more files in the stream"<<endl;
    return kFALSE;
  }

  const char * filename = 
    static_cast<TObjString*>(fFileNames->At(fCurrentFileIndex))->GetName();
  fCurrentFile = TFile::Open(filename,"READ");
  if (!fCurrentFile) {
// cannot open file specified on input. Do not skip it silently.
    cerr<<"Cannot open file "<<filename<<endl;
    return kFALSE;
  }
// find nr of events in the given file  
  TTree * te = (TTree *) fCurrentFile->Get("TE") ;
  if (!te) {
    Error("OpenNextFile", "input file does not contain TE");
    return kFALSE;
  }
  fEvents = static_cast<Int_t>(te->GetEntries());
  fLastEventSerialNr = -1;
  return kTRUE;
}

////////////////////////////////////////////////////////////////////////
Bool_t AliStream::ImportgAlice()
{
  if (fFileNames->GetLast() < 0) return kFALSE;
  if (!fCurrentFile) {
    if (!OpenNextFile()) return kFALSE;
  }
  gAlice = (AliRun*)fCurrentFile->Get("gAlice");
  if (!gAlice)  return kFALSE;
  return kTRUE;
}
