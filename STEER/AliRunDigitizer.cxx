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
Revision 1.4  2001/09/19 06:23:50  jchudoba
Move some tasks to AliStream and AliMergeCombi classes

Revision 1.3  2001/07/30 14:04:18  jchudoba
correct bug in the initialization

Revision 1.2  2001/07/28 10:44:32  hristov
Loop variable declared once; typos corrected

Revision 1.1  2001/07/27 12:59:00  jchudoba
Manager class for merging/digitization

*/

////////////////////////////////////////////////////////////////////////
//
// AliRunDigitizer.cxx
//
// Manager object for merging/digitization
//
// Instance of this class manages the digitization and/or merging of
// Sdigits into Digits. 
//
// Only one instance of this class is created in the macro:
//   AliRunDigitizer * manager = 
//      new AliRunDigitizer(nInputStreams,SPERB);
// where nInputStreams is number of input streams and SPERB is
// signals per background variable, which determines how combinations
// of signal and background events are generated.
// Then instances of specific detector digitizers are created:
//   AliMUONDigitizer *dMUON  = new AliMUONDigitizer(manager)
// and the I/O configured (you have to specify input files 
// and an output file). The manager connects appropriate trees from 
// the input files according a combination returned by AliMergeCombi 
// classcreates. It creates TreeD in the output and runs once per 
// event Digitize method of all existing AliDetDigitizers 
// (without any option). AliDetDigitizers ask manager
// for a TTree with input (manager->GetInputTreeS(Int_t i),
// merge all inputs, digitize it, and save it in the TreeD 
// obtained by manager->GetTreeD(). Output events are stored with 
// numbers from 0, this default can be changed by 
// manager->SetFirstOutputEventNr(Int_t) method. The particle numbers
// in the output are shifted by MASK, which is taken from manager.
//
// Single input file is permitted. Maximum MAXSTREAMSTOMERGE can be merged.
// Input from the memory (on-the-fly merging) is not yet 
// supported, as well as access to the input data by invoking methods
// on the output data.
//
// Access to the some data is via gAlice for now (supposing the 
// same geometry in all input files), gAlice is taken from the first 
// input file on the first stream.
//
// Example with MUON digitizer:
//
//  AliRunDigitizer * manager = new AliRunDigitizer(2,1);
//  manager->SetInputStream(0,"1track_10events_phi45_60.root");
//  manager->SetInputStream(1,"1track_10events_phi120_135.root");
//  manager->SetOutputDir("/tmp");
//  manager->SetOutputFile("digits.root");
//  AliMUONDigitizer *dMUON  = new AliMUONDigitizer(manager);
//  manager->SetNrOfEventsToWrite(1);
//  manager->Digitize();
//
//////////////////////////////////////////////////////////////////////// 

// system includes

#include <iostream.h>

// ROOT includes

#include "TFile.h"
#include "TTree.h"

// AliROOT includes

#include "AliRunDigitizer.h"
#include "AliDigitizer.h"
#include "AliRun.h"
#include "AliHeader.h"
#include "TParticle.h"
#include "AliStream.h"
#include "AliMergeCombi.h"

ClassImp(AliRunDigitizer)

////////////////////////////////////////////////////////////////////////

AliRunDigitizer::AliRunDigitizer() : TNamed("AliRunDigitizer","")
{
// default ctor
  cerr<<"Don't use"<<endl;
  fCombi = 0;
  fInputFiles = 0;
  fNDigitizers = 0;
  fNinputs = 0;
  fInputStreams = 0;
}

////////////////////////////////////////////////////////////////////////
AliRunDigitizer::AliRunDigitizer(Int_t nInputStreams, Int_t sperb) : TNamed("AliRunDigitizer","")
{
// default ctor
  if (nInputStreams == 0) {
    Error("AliRunDigitizer","Specify nr of input streams");
    return;
  }
  Int_t i;
  for (i=0;i<MAXDETECTORS;i++) fDigitizers[i]=0;
  fNDigitizers = 0;
  fNinputs = nInputStreams;
  fOutputFileName = "digits.root";
  fOutputDirName = "/tmp/";
  fCombination.Set(MAXSTREAMSTOMERGE);
  for (i=0;i<MAXSTREAMSTOMERGE;i++) {
    fArrayTreeS[i]=fArrayTreeH[i]=fArrayTreeTPCS[i]=NULL;
    fCombination[i]=-1;
  }
  fkMASKSTEP = 10000000;
  fkMASK[0] = 0;
  for (i=1;i<MAXSTREAMSTOMERGE;i++) {
    fkMASK[i] = fkMASK[i-1] + fkMASKSTEP;
  }
  fInputStreams = new TClonesArray("AliStream",nInputStreams);
  TClonesArray &lInputStreams = *fInputStreams;
  for (i=0;i<nInputStreams;i++) {
    new(lInputStreams[i]) AliStream();
  }
  fInputFiles = new TClonesArray("TFile",1);
  fOutput = 0;
  fEvent = 0;
  fNrOfEventsToWrite = 0;
  fNrOfEventsWritten = 0;
  fCopyTreesFromInput = -1;
  fCombi = new AliMergeCombi(nInputStreams,sperb);
  fDebug = 3;
  if (GetDebug()>2) 
    cerr<<"AliRunDigitizer::AliRunDigitizer() called"<<endl;
}

////////////////////////////////////////////////////////////////////////

AliRunDigitizer::~AliRunDigitizer() {
// dtor

  if (fInputFiles) {
    delete fInputFiles;
    fInputFiles = 0;
  }
  if (fInputStreams) {
    delete fInputStreams;
    fInputStreams = 0;
  }
  if (fCombi) {
    delete fCombi;
    fCombi = 0;
  }

}

////////////////////////////////////////////////////////////////////////
void AliRunDigitizer::AddDigitizer(AliDigitizer *digitizer)
{
// add digitizer to the list of active digitizers

  if (fNDigitizers >= MAXDETECTORS) {
    cerr<<"Too many detectors to digitize. Increase value of MAXDETECTORS"
	<<" constant in AliRunDigitizer.h and recompile or decrease the"
	<<" the number of detectors"<<endl;
  } else {
    fDigitizers[fNDigitizers++] = digitizer;
  }
}

////////////////////////////////////////////////////////////////////////

void AliRunDigitizer::SetInputStream(Int_t i, char *inputFile)
{
  if (i > fInputStreams->GetLast()) {
    Error("SetInputStream","Input stream number too high");
    return;
  }
  static_cast<AliStream*>(fInputStreams->At(i))->AddFile(inputFile);
}

////////////////////////////////////////////////////////////////////////
void AliRunDigitizer::Digitize()
{
// get a new combination of inputs, connect input trees and loop 
// over all digitizers

  if (!InitGlobal()) {
    cerr<<"False from InitGlobal"<<endl;
    return;
  }
// take gAlice from the first input file. It is needed to access
//  geometry data
  if (!static_cast<AliStream*>(fInputStreams->At(0))->ImportgAlice()) {
    cerr<<"gAlice object not found in the first file of "
	<<"the 1st stream"<<endl;
    return;
  }
  Int_t eventsCreated = 0;
  while (eventsCreated++ < fNrOfEventsToWrite) {
//    if (GetDebug()>2) PrintCombination();
    ConnectInputTrees();
    InitEvent();
// loop over all registered digitizers and let them do the work
    for (Int_t i=0;i<fNDigitizers; i++) {
      fDigitizers[i]->Digitize();
    }
    FinishEvent();
  }
  FinishGlobal();
}

////////////////////////////////////////////////////////////////////////
Bool_t AliRunDigitizer::ConnectInputTrees()
{
// fill arrays fArrayTreeS, fArrayTreeH and fArrayTreeTPCS with 
// pointers to the correct events according fCombination values
// null pointers can be in the output, AliDigitizer has to check it

  TTree *tree;
  char treeName[20];
  Int_t serialNr;
  Int_t eventNr[MAXSTREAMSTOMERGE], delta[MAXSTREAMSTOMERGE];
  fCombi->Combination(eventNr, delta);
  for (Int_t i=0;i<fNinputs;i++) {
    if (delta[i] == 1) {
      AliStream *iStream = static_cast<AliStream*>(fInputStreams->At(i));
      if (!iStream->NextEventInStream(serialNr)) return kFALSE;
      sprintf(treeName,"TreeS%d",serialNr);
      tree = static_cast<TTree*>(iStream->CurrentFile()->Get(treeName));
      fArrayTreeS[i] = tree;
      sprintf(treeName,"TreeH%d",serialNr);
      tree = static_cast<TTree*>(iStream->CurrentFile()->Get(treeName));
      fArrayTreeH[i] = tree;
      sprintf(treeName,"TreeS_75x40_100x60_%d",serialNr);
      tree = static_cast<TTree*>(iStream->CurrentFile()->Get(treeName));
      fArrayTreeTPCS[i] = tree;
    } else if (delta[i] != 0) {
      Error("ConnectInputTrees","Only delta 0 or 1 is implemented");
      return kFALSE;
    }
  }
  return kTRUE;
}

////////////////////////////////////////////////////////////////////////
Bool_t AliRunDigitizer::InitGlobal()
{
// called once before Digitize() is called, initialize digitizers and output

  if (!InitOutputGlobal()) return kFALSE;
  for (Int_t i=0;i<fNDigitizers; i++) {
    if (!fDigitizers[i]->Init()) return kFALSE;
  }
  return kTRUE;
}

////////////////////////////////////////////////////////////////////////
Bool_t AliRunDigitizer::InitOutputGlobal()
{
// Creates the output file, called once by InitGlobal()

  TString fn;
  fn = fOutputDirName + '/' + fOutputFileName;
  fOutput = new TFile(fn,"recreate");
  if (GetDebug()>2) {
    cerr<<"file "<<fn.Data()<<" was opened"<<endl;
    cerr<<"fOutput = "<<fOutput<<endl;
  }
  if (fOutput) return kTRUE;
  Error("InitOutputGlobal","Could not create output file.");
  return kFALSE;
}


////////////////////////////////////////////////////////////////////////
void AliRunDigitizer::InitEvent()
{
// Creates TreeDxx in the output file, called from Digitize() once for 
//  each event. xx = fEvent

  if (GetDebug()>2) 
    cerr<<"AliRunDigitizer::InitEvent: fEvent = "<<fEvent<<endl;
  fOutput->cd();
  char hname[30];
  sprintf(hname,"TreeD%d",fEvent);
  fTreeD = new TTree(hname,"Digits");
  fTreeD->Write();   // Do I have to write it here???
}

////////////////////////////////////////////////////////////////////////
void AliRunDigitizer::FinishEvent()
{
// called at the end of loop over digitizers

  fOutput->cd();
  if (fCopyTreesFromInput > -1) {
    char treeName[20];
    Int_t i = fCopyTreesFromInput; 
    sprintf(treeName,"TreeK%d",fCombination[i]);
    ((TFile*)fInputFiles->At(i))->Get(treeName)->Clone()->Write();
    sprintf(treeName,"TreeH%d",fCombination[i]);
    ((TFile*)fInputFiles->At(i))->Get(treeName)->Clone()->Write();
  }
  fEvent++;
  fNrOfEventsWritten++;
  if (fTreeD) {
    delete fTreeD;
    fTreeD = 0;
  }
}
////////////////////////////////////////////////////////////////////////
void AliRunDigitizer::FinishGlobal()
{
// called at the end of Exec
// save unique objects to the output file

  fOutput->cd();
  this->Write();
  if (fCopyTreesFromInput > -1) {
    ((TFile*)fInputFiles->At(fCopyTreesFromInput))->Get("TE")
      ->Clone()->Write();
    gAlice->Write();
  }
  fOutput->Close();
}


////////////////////////////////////////////////////////////////////////
Int_t  AliRunDigitizer::GetNParticles(Int_t event)
{
// return number of particles in all input files for a given
// event (as numbered in the output file)
// return -1 if some file cannot be accessed

  Int_t sum = 0;
  Int_t sumI;
  for (Int_t i = 0; i < fNinputs; i++) {
    sumI = GetNParticles(GetInputEventNumber(event,i), i);
    if (sumI < 0) return -1;
    sum += sumI;
  }
  return sum;
}

////////////////////////////////////////////////////////////////////////
Int_t  AliRunDigitizer::GetNParticles(Int_t event, Int_t input)
{
// return number of particles in input file input for a given
// event (as numbered in this input file)
// return -1 if some error

// Must be revised in the version with AliStream

  return -1;

/*
  TFile *file = ConnectInputFile(input);
  if (!file) {
    Error("GetNParticles","Cannot open input file");
    return -1;
  }

// find the header and get Nprimaries and Nsecondaries
  TTree* tE = (TTree *)file->Get("TE") ;
  if (!tE) {
    Error("GetNParticles","input file does not contain TE");
    return -1;
  }
  AliHeader* header;
  header = 0;
  tE->SetBranchAddress("Header", &header);
  if (!tE->GetEntry(event)) {
    Error("GetNParticles","event %d not found",event);
    return -1;
  }
  if (GetDebug()>2) {
    cerr<<"Nprimary: "<< header->GetNprimary()<<endl;
    cerr<<"Nsecondary: "<<header->GetNsecondary()<<endl;
  }
  return header->GetNprimary() + header->GetNsecondary();
*/
}

////////////////////////////////////////////////////////////////////////
Int_t* AliRunDigitizer::GetInputEventNumbers(Int_t event)
{
// return pointer to an int array with input event numbers which were
// merged in the output event event

// simplified for now, implement later
  Int_t a[MAXSTREAMSTOMERGE];
  for (Int_t i = 0; i < fNinputs; i++) {
    a[i] = event;
  }
  return a;
}
////////////////////////////////////////////////////////////////////////
Int_t AliRunDigitizer::GetInputEventNumber(Int_t event, Int_t input)
{
// return an event number of an eventInput from input file input
// which was merged to create output event event

// simplified for now, implement later
  return event;
}
////////////////////////////////////////////////////////////////////////
TParticle* AliRunDigitizer::GetParticle(Int_t i, Int_t event)
{
// return pointer to particle with index i (index with mask)

// decode the MASK
  Int_t input = i/fkMASKSTEP;
  return GetParticle(i,input,GetInputEventNumber(event,input));
}

////////////////////////////////////////////////////////////////////////
TParticle* AliRunDigitizer::GetParticle(Int_t i, Int_t input, Int_t event)
{
// return pointer to particle with index i in the input file input
// (index without mask)
// event is the event number in the file input
// return 0 if it does not exist

// Must be revised in the version with AliStream

  return 0;
/*
  TFile *file = ConnectInputFile(input);
  if (!file) {
    Error("GetParticle","Cannot open input file");
    return 0;
  }

// find the header and get Nprimaries and Nsecondaries
  TTree* tE = (TTree *)file->Get("TE") ;
  if (!tE) {
    Error("GetParticle","input file does not contain TE");
    return 0;
  }
  AliHeader* header;
  header = 0;
  tE->SetBranchAddress("Header", &header);
  if (!tE->GetEntry(event)) {
    Error("GetParticle","event %d not found",event);
    return 0;
  }
  
// connect TreeK  
  char treeName[30];
  sprintf(treeName,"TreeK%d",event);  
  TTree* tK = static_cast<TTree*>(file->Get(treeName));
  if (!tK) {
    Error("GetParticle","input file does not contain TreeK%d",event);
    return 0;
  }
  TParticle *particleBuffer;
  particleBuffer = 0;
  tK->SetBranchAddress("Particles", &particleBuffer);


// algorithmic way of getting entry index
// (primary particles are filled after secondaries)
  Int_t entry;
  if (i<header->GetNprimary())
    entry = i+header->GetNsecondary();
  else 
    entry = i-header->GetNprimary();
  Int_t bytesRead = tK->GetEntry(entry);
//  new ((*fParticles)[nentries]) TParticle(*fParticleBuffer);
  if (bytesRead)
    return particleBuffer;
  return  0;
*/
}
////////////////////////////////////////////////////////////////////////

