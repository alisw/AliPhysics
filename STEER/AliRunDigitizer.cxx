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
// AliRunDigitizer.cxx
//
// Manager object for merging/digitization
//
// Instance of this class manages the digitization and/or merging of
// Sdigits into Digits. 
//
// Only one instance of this class is created in the macro:
//   AliRunDigitizer * manager = new AliRunDigitizer();
// Then instances of specific detector digitizers are created:
//   AliMUONDigitizer *dMUON  = new AliMUONDigitizer(manager)
// and the manager is configured (you have to specify input files 
// and an output file). The manager generates a combination of input 
// events according to an option set by SetCombinationType. Then it
// connects appropriate trees from the input files, creates TreeD 
// in the output and runs once per event Digitize method of all existing 
// AliDetDigitizers (without any option). AliDetDigitizers ask manager
// for a TTree with input (manager->GetNextTreeH(TTree *) 
// or manager->GetNextTreeS(TTree *),
// merge all inputs, digitize it, and save it in the TreeD 
// obtained by manager->GetTreeD(). Output events are stored with 
// numbers from 0, this default can be changed by 
// manager->SetFirstOutputEventNr(Int_t) method. The particle numbers
// in the output are shifted by MASK, which is taken from manager.
//
// Single input file is permitted. Maximum MAXFILESTOMERGE can be merged.
// Input from the memory (on-the-fly merging) is not yet 
// supported, as well as access to the input data by invoking methods
// on the output data.
//
// Access to the geometrical data is via gAlice for now (supposing the 
// same geometry in all input files), gAlice is taken from the defined 
// input file.
//
// Example with MUON digitizer:
//
//   AliRunDigitizer * manager = new AliRunDigitizer();
//   manager->SetInput("1track_10events_phi45_60.root");
//   manager->SetInput("1track_10events_phi120_135.root");
//   manager->SetOutputDir("/home/z2/jchudoba/batch/jobtmp");
//   manager->SetOutputFile("digits.root");
//   AliMUONDigitizer *dMUON  = new AliMUONDigitizer(manager);
//   manager->SetNrOfEventsToWrite(3);
//   manager->SetCopyTreesFromInput(1);
//   manager->Digitize();
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

ClassImp(AliRunDigitizer)

////////////////////////////////////////////////////////////////////////

  AliRunDigitizer::AliRunDigitizer() : TNamed("AliRunDigitizer","")
{
// default ctor

  for (Int_t i=0;i<MAXDETECTORS;i++) fDigitizers[i]=0;
  fNDigitizers = 0;
  fNinputs = 0;
  fOutputFileName = "digits.root";
  fOutputDirName = "/tmp/";
  fCombination.Set(MAXFILESTOMERGE);
  for (Int_t i=0;i<MAXFILESTOMERGE;i++) {
    fArrayTreeS[i]=fArrayTreeH[i]=fArrayTreeTPCS[i]=NULL;
    fCombination[i]=-1;
  }
  fkMASKSTEP = 10000000;
  fkMASK[0] = 0;
  for (Int_t i=0;i<MAXFILESTOMERGE;i++) {
    fkMASK[i] = fkMASK[i-1] + fkMASKSTEP;
  }
  fInputFileNames = new TClonesArray("TObjString",1);
  fInputFiles = new TClonesArray("TFile",1);
  fMinNEvents = 99999999;
  fCombinationType=1;
  fCombinationFileName = "combination.txt";
  fOutput = 0;
  fEvent = 0;
  fNrOfEventsToWrite = 0;
  fNrOfEventsWritten = 0;
  fCopyTreesFromInput = -1;
  fDebug = 3;
  if (GetDebug()>2) 
    cerr<<"AliRunDigitizer::AliRunDigitizer() called"<<endl;
}

////////////////////////////////////////////////////////////////////////

AliRunDigitizer::~AliRunDigitizer() {
// dtor

  if (fInputFiles) delete fInputFiles;
  if (fInputFileNames) delete fInputFileNames;
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

Bool_t AliRunDigitizer::SetInput(char *inputFile)
{
// receives the name of the input file. Opens it and stores pointer 
// to it, returns kFALSE if fails
// if inputFile = MEMORY, uses pointers from gAlice - not yet implemented

// input cannot be output - open READ only
  TFile *file = new((*fInputFiles)[fNinputs]) TFile(inputFile,"READ"); 

  if (GetDebug()>2) cerr<<"AliRunDigitizer::SetInput: file = "<<file<<endl;
  if (!file->IsOpen()) {
    cerr<<"ERROR: AliRunDigitizer::SetInput: cannot open file "
	<<inputFile<<endl;
    return kFALSE;
  }

// find Header and get nr of events there
  TTree * te = (TTree *) file->Get("TE") ;
  if (!te) {
    cerr<<"ERROR: AliRunDigitizer::SetInput:  input file does "
	<<"not contain TE"<<endl;
    return kFALSE;
  }
  Int_t nEntries = (Int_t) te->GetEntries();

  if (GetDebug()>2) cerr<<"AliRunDigitizer::SetInput: nEntries = "
			<<nEntries<<endl;
  if (nEntries < 1) {
    cerr<<"ERROR: AliRunDigitizer::SetInput:  input file does "
	<<"not contain any event"<<endl;
    return kFALSE;
  }

  if (nEntries < fMinNEvents) fNrOfEventsToWrite = fMinNEvents = nEntries;

// find gAlice object if it is a first input file
//  this is only temporary solution, we need gAlice to access detector
//  geometry parameters. Unfortunately I have to include AliRun header file.
  if (fNinputs == 0) {
    gAlice = (AliRun*)file->Get("gAlice");
    if (GetDebug() > 2) cerr<<"gAlice taken from the first input: "
			    <<gAlice<<endl;
    if (!gAlice) {
      cerr<<"ERROR: AliRunDigitizer::SetInput:  first input file "
	  <<"does not contain gAlice object"<<endl;
      return kFALSE;
    }
  }
    
// store this file name if it is OK
  new((*fInputFileNames)[fNinputs])  TObjString(inputFile);
  fNinputs++;
  if (GetDebug() > 2) cerr<<"fNinputs = "<<fNinputs<<endl;
  return kTRUE;
}

////////////////////////////////////////////////////////////////////////
Bool_t AliRunDigitizer::MakeCombination()
{
// make a new combination of events from different files

  Int_t type = fCombinationType;
  if (fNrOfEventsWritten >= fNrOfEventsToWrite) return kFALSE;

  switch (type) {
	
  case 1:
// type = 1:  1-1-1 - take the same event number from each file
    if (fCombination[0]<fMinNEvents-1) {	
      for (Int_t i=0;i<fNinputs;i++) fCombination[i]++;
      return kTRUE;
    }
    if (GetDebug()>2)
      cerr<<"AliRunDigitizer::MakeCombination: Warning: "
	  <<"maximum number of Events in an input file "
	  <<"was reached."<<endl;
    break;

  case 2:
// type = 2: read from the file combinations.ascii
//  not yet implemented 100% correctly - requires 4 entries in the row
    static FILE *fp ;
    static Int_t linesRead;
    if (!fp) {
      fp = fopen(fCombinationFileName.Data(),"r");
      linesRead = 0;
    }
    if (!fp) {
      cerr<<"AliRunDigitizer::MakeCombination ERROR: "
	  <<"Cannot open input file with combinations."<<endl;
      return kFALSE;
    }
    Int_t var[4], nInputs;
    char line[80];
// since I do not close or rewind the file, the position should be correct
    if (fgets(line,80,fp)) {
      nInputs = sscanf(&line[0],"%d%d%d%d",&var[0],&var[1],&var[2],&var[3]);
      if (nInputs != fNinputs) {
	cerr<<"AliRunDigitizer::MakeCombination ERROR: "
	    <<"Nr. of input files is different from nr "
	    <<"integers defining the combination"<<endl;
	return kFALSE;
      }
      while(nInputs--) {
	fCombination[nInputs] = var[nInputs];
      }
      return kTRUE;
    } else {
      cerr<<"AliRunDigitizer::MakeCombination ERROR: "
	  <<"no more input in the file with combinations"<<endl;
      return kFALSE;
    }

  default:
    cerr<<"AliRunDigitizer::MakeCombination: ERROR: "
	<<"wrong type of required combination type: "<<type<<endl;
  }
  return kFALSE;
}

////////////////////////////////////////////////////////////////////////
void AliRunDigitizer::PrintCombination()
{
// debug method to print current combination

  cerr<<"AliRunDigitizer::PrintCombination: Current events combination:"<<endl;
  for (Int_t i=0;i<fNinputs;i++) {
    cerr<<"File: "<<((TObjString *)fInputFileNames->At(i))->GetString()<<"\tEvent: "<<fCombination[i]<<endl;
  }
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
  while (MakeCombination()) {
    if (GetDebug()>2) PrintCombination();
    ConnectInputTrees();
    InitPerEvent();
// loop over all registered digitizers and let them do the work
    for (Int_t i=0;i<fNDigitizers; i++) {
      fDigitizers[i]->Digitize();
    }
    FinishPerEvent();
  }
  FinishGlobal();
}

////////////////////////////////////////////////////////////////////////
void AliRunDigitizer::ConnectInputTrees()
{
// fill arrays fArrayTreeS, fArrayTreeH and fArrayTreeTPCS with 
// pointers to the correct events according fCombination values
// null pointers can be in the output, AliDigitizer has to check it

  TFile *file;
  TTree *tree;
  char treeName[20];
  for (Int_t i=0;i<fNinputs;i++) {
    file = (TFile*)(*fInputFiles)[i];
    sprintf(treeName,"TreeS%d",fCombination[i]);
    tree = (TTree *) file->Get(treeName);
    fArrayTreeS[i] = tree;
    sprintf(treeName,"TreeH%d",fCombination[i]);
    tree = (TTree *) file->Get(treeName);
    fArrayTreeH[i] = tree;
    sprintf(treeName,"TreeD_75x40_100x60_%d",fCombination[i]);
    tree = (TTree *) file->Get(treeName);
    fArrayTreeTPCS[i] = tree;
   }
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
void AliRunDigitizer::InitPerEvent()
{
// Creates TreeDxx in the output file, called from Digitize() once for 
//  each event. xx = fEvent

  if (GetDebug()>2) 
    cerr<<"AliRunDigitizer::InitPerEvent: fEvent = "<<fEvent<<endl;
  fOutput->cd();
  char hname[30];
  sprintf(hname,"TreeD%d",fEvent);
  fTreeD = new TTree(hname,"Digits");
//  fTreeD->Write();   // Do I have to write it here???
}


////////////////////////////////////////////////////////////////////////
TTree* AliRunDigitizer::GetNextTreeH(TTree *current) const
{
// returns next one after the current one
// returns the first if the current is NULL
// returns NULL if the current is the last one

  if (fNinputs <= 0) return 0;
  if (current == 0) return fArrayTreeH[0];
  for (Int_t i=0; i<fNinputs-1; i++) {
    if (current == fArrayTreeH[i]) return fArrayTreeH[i+1];
  }
  return 0;
}
////////////////////////////////////////////////////////////////////////
TTree* AliRunDigitizer::GetNextTreeS(TTree *current) const
{
// returns next one after the current one
// returns the first if the current is NULL
// returns NULL if the current is the last one

  if (fNinputs <= 0) return 0;
  if (current == 0) return fArrayTreeS[0];
  for (Int_t i=0; i<fNinputs-1; i++) {
    if (current == fArrayTreeS[i]) return fArrayTreeS[i+1];
  }
  return 0;
}
////////////////////////////////////////////////////////////////////////
TTree* AliRunDigitizer::GetNextTreeTPCS(TTree *current) const
{
// returns next one after the current one
// returns the first if the current is NULL
// returns NULL if the current is the last one

  if (fNinputs <= 0) return 0;
  if (current == 0) return fArrayTreeTPCS[0];
  for (Int_t i=0; i<fNinputs-1; i++) {
    if (current == fArrayTreeTPCS[i]) return fArrayTreeTPCS[i+1];
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////
Int_t AliRunDigitizer::GetNextMask(Int_t current) const
{
// returns next one after the current one
// returns the first if the current is negative
// returns negative if the current is the last one

  if (fNinputs <= 0) return -1;
  if (current < 0) return fkMASK[0]; 
  for (Int_t i=0; i<fNinputs-1; i++) {
    if (current == fkMASK[i]) return fkMASK[i+1];
  }
  return -1;
}
////////////////////////////////////////////////////////////////////////
void AliRunDigitizer::FinishPerEvent()
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
}

////////////////////////////////////////////////////////////////////////
Int_t* AliRunDigitizer::GetInputEventNumbers(Int_t event)
{
// return pointer to an int array with input event numbers which were
// merged in the output event event

// simplified for now, implement later
  Int_t a[MAXFILESTOMERGE];
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
TFile* AliRunDigitizer::ConnectInputFile(Int_t input)
{
// open input file with index input
// 1st file has index 0
// return 0x0 if fails

// check if file with index input is already open
  TFile *file = 0;
  if (fInputFiles->GetEntriesFast() > input)
    file = static_cast<TFile *>(fInputFiles->At(input));

  if (!file) {
    // find the file name and open it
    TObjString *fn = static_cast<TObjString *>(fInputFileNames->At(input));
    file = new((*fInputFiles)[input]) TFile((fn->GetString()).Data(),"READ");
    if (!file) {
      Error("ConnectInputFile","Cannot open input file");
      return 0;
    }
    if (!file->IsOpen()) {
      Error("ConnectInputFile","Cannot open input file");
      return 0;
    }
  }
  return file;
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
// return 0 i fit does not exist

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
//  tK->GetEntry(0);          // do I need this???
  Int_t bytesRead = tK->GetEntry(entry);
//  new ((*fParticles)[nentries]) TParticle(*fParticleBuffer);
  if (bytesRead)
    return particleBuffer;
  return  0;
}
////////////////////////////////////////////////////////////////////////

