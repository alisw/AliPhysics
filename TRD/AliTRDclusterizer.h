#ifndef ALITRDCLUSTERIZER_H
#define ALITRDCLUSTERIZER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TNamed.h>

class TFile;

///////////////////////////////////////////////////////
//  Finds and handles cluster                        //
///////////////////////////////////////////////////////

class AliTRDclusterizer : public TNamed {

 public:

  AliTRDclusterizer();
  AliTRDclusterizer(const Text_t* name, const Text_t* title);
  AliTRDclusterizer(const AliTRDclusterizer &c);
  virtual ~AliTRDclusterizer();
  AliTRDclusterizer &operator=(const AliTRDclusterizer &c);

  virtual void    Copy(TObject &c);
  virtual Bool_t  Open(const Char_t *name, Int_t nEvent = 0);
  virtual Bool_t  Open(const Char_t *inname, const Char_t *outname, Int_t nEvent = 0);
  virtual Bool_t  OpenInput(const Char_t *name, Int_t nEvent = 0);
  virtual Bool_t  OpenOutput(const Char_t *name);
  virtual Bool_t  MakeClusters() = 0;
  virtual Bool_t  WriteClusters(Int_t det);

          void    SetVerbose(Int_t v = 1) { fVerbose = v; };

 protected:

  TFile   *fInputFile;             //! AliROOT input file
  TFile   *fOutputFile;            //! AliROOT output file
  TTree   *fClusterTree;           //! Tree with the cluster
  AliTRD  *fTRD;                   //! The TRD object

  Int_t    fEvent;                 //! Event number
  Int_t    fVerbose;               //  Sets the verbose level

  ClassDef(AliTRDclusterizer,2)    //  TRD-Cluster manager base class

};

#endif
