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
// Generator for vertices taken from a file                                  //
//                                                                           //
// The file name of the galice file is passed as argument to the             //
// constructor. If a second argument is given, this determines the number    //
// of events for which the same vertex is used.                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <TArrayF.h>
#include <TFile.h>
#include <TTree.h>

#include "AliVertexGenFile.h"
#include "AliLog.h"
#include "AliGenEventHeader.h"
#include "AliHeader.h"


ClassImp(AliVertexGenFile)


//_____________________________________________________________________________
AliVertexGenFile::AliVertexGenFile() :
  fFile(NULL),
  fTree(NULL),
  fHeader(NULL),
  fEventsPerEntry(0),
  fEvent(0)
{
// default constructor: initialize data members

}

//_____________________________________________________________________________
AliVertexGenFile::AliVertexGenFile(const char* fileName, 
				   Int_t eventsPerEntry) :
  fFile(NULL),
  fTree(NULL),
  fHeader(NULL),
  fEventsPerEntry(eventsPerEntry),
  fEvent(0)
{
// main constructor:
// fileName is the name of the galice file containing the vertices
// eventsPerEntry is the number of events for which the same vertex is used

  TDirectory* dir = gDirectory;

  fFile = TFile::Open(fileName);
  if (!fFile || !fFile->IsOpen()) {
    AliError(Form("could not open file %s", fileName));
    delete fFile;
    fFile = NULL;
    return;
  }
  fTree = (TTree*) fFile->Get("TE");
  if (!fTree) {
    AliError(Form("no header tree found in file %s", fileName));
    dir->cd();
    return;
  }
  fHeader = new AliHeader;
  fTree->SetBranchAddress("Header", &fHeader);

  dir->cd();
}

//_____________________________________________________________________________
AliVertexGenFile::~AliVertexGenFile()
{
// clean up

  if (fFile) fFile->Close();
  delete fFile;
  delete fHeader;
}

//_____________________________________________________________________________
time_t AliVertexGenFile::GetHeaderTimeStamp() const
{
  // get the timestamp of the last header used for the vertex
  if (fHeader) return fHeader->GetTimeStamp();
  else AliFatal("No header was loaded yet");
  return 0;
}

//_____________________________________________________________________________
TVector3 AliVertexGenFile::GetVertex()
{
// get the vertex from the event header tree

  Int_t entry = fEvent++ / fEventsPerEntry;
  if (!fTree) {
    AliError("no header tree");
    return TVector3(0,0,0);
  }

  if (fTree->GetEntry(entry) <= 0) {
    AliError(Form("error loading entry %d", entry));
    return TVector3(0,0,0);
  }

  if (!fHeader->GenEventHeader()) {
    AliError("no generator event header");
    return TVector3(0,0,0);
  }

  TArrayF vertex(3);
  fHeader->GenEventHeader()->PrimaryVertex(vertex);
  return TVector3(vertex[0], vertex[1], vertex[2]);
}


