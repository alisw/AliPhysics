#ifndef IceRoot_h
#define IceRoot_h

// Copyright(c) 2003, IceCube Experiment at the South Pole, All rights reserved.
// See cxx source for full Copyright notice.

// $Id$

#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"

#include "AliJob.h"
#include "AliObjMatrix.h"

#include "IceAOM.h"
#include "IceEvent.h"

class IceRoot : public AliJob
{
 public :
  IceRoot(const char* name="IceRoot",const char* title="");   // Constructor
  virtual ~IceRoot();                                         // Destructor
  void SetMaxEvents(Int_t n);       // Set maximum number of events to be processed
  void SetPrintFreq(Int_t f);       // Set printfrequency to provide info every f events
  void SetSplitLevel(Int_t split);  // Set split level for the produced ROOT data file
  void SetBufferSize(Int_t bsize);  // Set buffersize for the produced ROOT data file
  void AddInputFile(TString name);  // Add name of simple Root input file to the list
  void SetOutputFile(TFile* ofile); // Set output file for the ROOT data structures       
  void SetOutputFile(TString name); // Create output file for the ROOT data structures
  TFile* GetOutputFile();           // Provide pointer to the ROOT output file
  virtual void Exec(Option_t* opt); // Perform the format conversion

 protected :
  Int_t fSplit;        // The split level of the produced ROOT data file
  Int_t fBsize;        // The buffersize of the produced ROOT data file
  Int_t fMaxevt;       // The maximum number of events to be processed
  Int_t fPrintfreq;    // The event info printing frequency
  TObjArray* fInfiles; // Names of all the simple Root data input files
  TFile* fOutfile;     // The ROOT output file
  void PutTrigger(Int_t year);   // Put the trigger info from the simple Root data event into the IcePack structure
  void PutWaveforms(Int_t year); // Put the waveforms from the simple Root data event into the IcePack structure

  TFile* fInput;       // Simple Root data input file
  TTree* fTree;        // Tree with simple Root data

 ClassDef(IceRoot,1) // Job for conversion of simple Root data into IceEvent data structures.
};
#endif
