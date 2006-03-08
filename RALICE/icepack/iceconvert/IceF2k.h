#ifndef IceF2k_h
#define IceF2k_h

// Copyright(c) 2003, IceCube Experiment at the South Pole, All rights reserved.
// See cxx source for full Copyright notice.

// $Id$

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"

#include "AliJob.h"
#include "AliObjMatrix.h"

#include "IceAOM.h"
#include "IceEvent.h"

#include "rdmc.h"

class IceF2k : public AliJob
{
 public :
  IceF2k(const char* name="IceF2k",const char* title=""); // Constructor
  virtual ~IceF2k();                                      // Destructor
  void SetMaxEvents(Int_t n);                             // Set maximum number of events to be processed
  void SetPrintFreq(Int_t f);                             // Set printfrequency to provide info every f events
  void SetSplitLevel(Int_t split);                        // Set split level for the produced ROOT data file
  void SetBufferSize(Int_t bsize);                        // Set buffersize for the produced ROO data file
  void SetInputFile(TString name);                        // Set name of F2K input file (obsolete)
  void AddInputFile(TString name);                        // Add name of F2K input file to the list
  void SetOutputFile(TFile* ofile);                       // Set output file for the ROOT data structures       
  void SetOutputFile(TString name);                       // Create output file for the ROOT data structures
  void SetMcToffset(Float_t toffset);                     // Set user defined time offset for MC data
  TFile* GetOutputFile();                                 // Provide pointer to the ROOT output file
  TDatabasePDG* GetPDG();           // Provide pointer to the PDG database
  AliObjMatrix* GetOMdbase();       // Provide pointer to the OM geometry, calib. etc... database
  AliDevice* GetFitdefs();          // Provide pointer to the Fit definition parameters
  AliDevice* GetTrigdefs();         // Provide pointer to the Fit definition parameters
  virtual void Exec(Option_t* opt); // Perform the format conversion

 protected :
  Int_t fSplit;        // The split level of the produced ROOT data file
  Int_t fBsize;        // The buffersize of the produced ROOT data file
  Int_t fMaxevt;       // The maximum number of events to be processed
  Int_t fPrintfreq;    // The event info printing frequency
  TObjArray* fInfiles; // Names of all the F2K input files
  TFile* fOutfile;     // The ROOT output file

  TDatabasePDG* fPdg;   // Database with PDG information
  AliObjMatrix* fOmdb;  // Database of all OM devices with their geometry, calib. etc... data
  AliDevice* fFitdefs;  // Fit definitions as indicated in the header of the F2000 input file
  AliDevice* fTrigdefs; // Trigger definitions as indicated in the header of the F2000 input file
  Float_t fToffset;     // Trigger time offset which might have been introduced during filtering
  Float_t fMctoffset;   // Trigger time offset which might have been introduced during MC generation

  void FillOMdbase();   // Fill geometry and calib. parameters of all devices
  void SetFitdefs();    // Set the fit definitions as used in the F2000 input file
  void SetTrigdefs();   // Set the trigger definitions as used in the F2000 input file
  void PutMcTracks();   // Put the MC tracks from the F2000 file into the IcePack structure
  void PutRecoTracks(); // Put the reconstructed tracks from the F2000 file into the IcePack structure
  void PutHits();       // Put the hits and waveforms from the F2000 file into the IcePack structure
  void PutTrigger();    // Put the trigger information from the F2000 file into the IcePack structure

  mcfile* fInput;  //! Structure holding the input file characteristics
  array   fHeader; //! Structure holding the file header info
  mevt    fEvent;  //! Structure holding the actual event data (hits, tracks, etc...)

 ClassDef(IceF2k,6) // Job for conversion of F2K data into IceEvent physics event structures.
};
#endif
