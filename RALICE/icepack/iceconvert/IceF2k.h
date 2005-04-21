#ifndef IceF2k_h
#define IceF2k_h

// Copyright(c) 2003, IceCube Experiment at the South Pole, All rights reserved.
// See cxx source for full Copyright notice.

// $Id$

#include "TObject.h"
#include "TChain.h"
#include "TFile.h"
#include "TDatabasePDG.h"
#include "TString.h"

#include "AliObjMatrix.h"

#include "IceAOM.h"
#include "IceEvent.h"

#include "rdmc.h"

class IceF2k : public TObject
{
 public :
  IceF2k(char* fname=0,Int_t split=0,Int_t bsize=32000);         // Constructor
  virtual ~IceF2k();                                             // Destructor
  void Loop(TTree* otree=0,Int_t nentries=-1,Int_t printfreq=1); // Perform the format conversion
  TDatabasePDG* GetPDG();     // Provide pointer to the PDG database
  AliObjMatrix* GetOMdbase(); // Provide pointer to the OM geometry, calib. etc... database
  AliDevice* GetFitdefs();    // Provide pointer to the Fit definition parameters

 protected :
  Int_t fSplit; // The split level of the produced ROOT data file
  Int_t fBsize; // The buffersize of the produced ROOT data file

  TDatabasePDG* fPdg;  // Database with PDG information
  AliObjMatrix* fOmdb; // Database of all OM devices with their geometry, calib. etc... data
  AliDevice* fFitdefs; // Fit definitions as indicated in the header of the F2000 input file

  void FillOMdbase();                // Fill geometry and calib. parameters of all devices
  void SetFitdefs();                 // Set the fit definitions as used in the F2000 input file
  void PutMcTracks(IceEvent* evt);   // Put the MC tracks from the F2000 file into the IcePack structure
  void PutRecoTracks(IceEvent* evt); // Put the reconstructed tracks from the F2000 file into the IcePack structure
  void PutHits(IceEvent* evt);       // Put the hits and waveforms from the F2000 file into the IcePack structure

  mcfile* fInput;  //! Structure holding the input file characteristics
  array   fHeader; //! Structure holding the file header info
  mevt    fEvent;  //! Structure holding the actual event data (hits, tracks, etc...)

 ClassDef(IceF2k,1) // Conversion of F2K data into IceEvent physics event structures.
};
#endif
