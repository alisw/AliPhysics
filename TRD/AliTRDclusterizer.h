#ifndef TRDclusterizer_h
#define TRDclusterizer_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TNamed.h>
#include <TFile.h>

///////////////////////////////////////////////////////
//  Finds and handles cluster                        //
///////////////////////////////////////////////////////

class AliTRDclusterizer : public TNamed {

 public:

  AliTRDclusterizer();
  AliTRDclusterizer(const Text_t* name, const Text_t* title);
  ~AliTRDclusterizer();

  virtual void    Init();
  virtual Bool_t  Open(const Char_t *name, Int_t nEvent = 0);
  virtual Bool_t  MakeCluster() = 0;
  virtual Bool_t  WriteCluster();

 protected:

  TFile   *fInputFile;             //! AliROOT input file
  
  Int_t    fEvent;                 //! Event number

  ClassDef(AliTRDclusterizer,1)    // TRD-Cluster manager base class

};

//_____________________________________________________________________________
class AliTRDcluster : public TObject {

public:

  Int_t    fDetector;       // TRD detector number

  Int_t    fTimeSlice;      // Timeslice in chamber where cluster has been found
  Float_t  fEnergy;         // Charge sum of this cluster

  Float_t  fX;              // X coord in ALICE reference frame
  Float_t  fY;              // Y coord in ALICE reference frame
  Float_t  fZ;              // Z coord in ALICE reference frame

  Int_t    fTracks[3];      // Track information

public:

  AliTRDcluster() {};
  AliTRDcluster(Int_t *tracks, Int_t *cluster, Float_t energy, Float_t *pos);
  virtual ~AliTRDcluster() {};

  inline virtual Int_t *GetTracks() { return &fTracks[0]; }

  ClassDef(AliTRDcluster,1) // Cluster for Transition Radiation Detector

};

#endif
