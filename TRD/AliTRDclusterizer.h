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

#endif
