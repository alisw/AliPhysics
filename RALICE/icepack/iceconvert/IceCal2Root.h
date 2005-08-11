#ifndef IceCal2Root_h
#define IceCal2Root_h

// Copyright(c) 2003, IceCube Experiment at the South Pole, All rights reserved.
// See cxx source for full Copyright notice.

// $Id$

#include "TFile.h"
#include "TString.h"
#include "TDatabasePDG.h"

#include "AliJob.h"
#include "AliObjMatrix.h"

#include "IceAOM.h"

#include "Riostream.h"

class IceCal2Root : public AliJob
{
 public :
  IceCal2Root(const char* name="IceCal2Root",const char* title=""); // Constructor
  virtual ~IceCal2Root();                                           // Destructor
  void SetAmacalibFile(TString name); // Set name of the Amacalib input file
  void SetOutputFile(TString name);   // Set output file for the ROOT data structures           
  TDatabasePDG* GetPDG();             // Provide pointer to the PDG database
  AliObjMatrix* GetOMdbase();         // Provide pointer to the OM geometry, calib. etc... database
  virtual void Exec(Option_t* opt);   // Perform the format conversion

 protected :
  ifstream fInput;         // Input stream for generic use of reading data

  TString fAmacalFileName; // Name of the Amacalib input file
  TString fBadomFileName;  // Name of the bad OM input file
  TString fRootFileName;   // Name of the ROOT output file
  TFile* fOutfile;         // The ROOT output file

  TDatabasePDG* fPdg;      // Database with PDG information
  AliObjMatrix* fOmdb;     // Database of all OM devices with their geometry, calib. etc... data

  void GetCalibData();     // Fill geometry, calibration and Xtalk parameters of all devices

 ClassDef(IceCal2Root,1) // Job for conversion of Amacalib ascii data into an AliObjMatrix OM dbase
};
#endif
