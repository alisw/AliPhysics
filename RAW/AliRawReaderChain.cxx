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

///////////////////////////////////////////////////////////////////////////////
///
/// This is a class for reading raw data from a root chain.
/// There are two constructors available - one from a text file containing the
/// list of root raw-data files to be processed and one directly from
/// TFileCollection.
///
/// cvetan.cheshkov@cern.ch 29/07/2008
///
///////////////////////////////////////////////////////////////////////////////

#include <TChain.h>
#include <TFileCollection.h>

#include "AliRawReaderChain.h"
#include "AliRawEvent.h"

ClassImp(AliRawReaderChain)

AliRawReaderChain::AliRawReaderChain() :
  AliRawReaderRoot(),
  fChain(NULL)
{
  // default constructor
}

AliRawReaderChain::AliRawReaderChain(TString listFileName) :
  AliRawReaderRoot(),
  fChain(NULL)
{
// create raw-reader objects which takes as an input a root chain
// from the file list found in 'listFileName'

  TFileCollection collection("RAW",
			     "Collection with raw-data files",
			     listFileName.Data());

  TChain* fChain = new TChain("RAW");
  if (!fChain->AddFileInfoList((TCollection*)(collection.GetList()))) {
    Error("AliRawReaderChain","Bad file list in collection, the chain is empty");
    return;
  }

  fChain->SetBranchStatus("*",0);
  fChain->SetBranchStatus("rawevent",1);

  fEvent = new AliRawEvent;
  fChain->SetBranchAddress("rawevent", &fEvent);
}

AliRawReaderChain::AliRawReaderChain(TFileCollection *collection) :
  AliRawReaderRoot(),
  fChain(NULL)
{
// create raw-reader objects which takes as an input a root chain
// from a root file collection

  TChain* fChain = new TChain("RAW");
  if (!fChain->AddFileInfoList((TCollection*)(collection->GetList()))) {
    Error("AliRawReaderChain","Bad file list in collection, the chain is empty");
    return;
  }

  fChain->SetBranchStatus("*",0);
  fChain->SetBranchStatus("rawevent",1);

  fEvent = new AliRawEvent;
  fChain->SetBranchAddress("rawevent", &fEvent);
}

AliRawReaderChain::AliRawReaderChain(const AliRawReaderChain& rawReader) :
  AliRawReaderRoot(rawReader),
  fChain(rawReader.fChain)
{
// copy constructor
}

AliRawReaderChain& AliRawReaderChain::operator = (const AliRawReaderChain& 
						  rawReader)
{
// assignment operator

  this->~AliRawReaderChain();
  new(this) AliRawReaderChain(rawReader);
  return *this;
}

AliRawReaderChain::~AliRawReaderChain()
{
// delete objects and close root file

  if (fChain) {
    delete fChain;
    fChain = NULL;
  }
}

Bool_t AliRawReaderChain::NextEvent()
{
// go to the next event in the root file

  if (!fChain || !fChain->GetListOfFiles()->GetEntriesFast()) return kFALSE;

  do {
    delete fEvent;
    fEvent = new AliRawEvent;
    if (fChain->GetEntry(fEventIndex+1) <= 0)
      return kFALSE;
    fEventIndex++;
  } while (!IsEventSelected());
  fEventNumber++;
  return Reset();
}

Bool_t AliRawReaderChain::RewindEvents()
{
// go back to the beginning of the root file

  fEventIndex = -1;
  delete fEvent;
  fEvent = new AliRawEvent;
  fEventNumber = -1;
  return Reset();
}
