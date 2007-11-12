#ifndef IceRootx_h
#define IceRootx_h

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

class IceRootx : public AliJob
{
 public :
  IceRootx(const char* name="IceRootx",const char* title="");   // Constructor
  virtual ~IceRootx();                                         // Destructor
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

  TTree* fTree;        // Tree with simple Root data

 ClassDef(IceRootx,1) // Job for conversion of simple Root data into IceEvent data structures.
};
#endif
