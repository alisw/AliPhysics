#ifndef ALIVERTEXGENFILE_H
#define ALIVERTEXGENFILE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliVertexGenerator.h"

class TFile;
class TTree;
class AliHeader;


class AliVertexGenFile: public AliVertexGenerator {
 public:
  AliVertexGenFile();
  AliVertexGenFile(const char* fileName, Int_t eventsPerEntry = 1);
  virtual ~AliVertexGenFile();

  virtual TVector3 GetVertex();

 private:
  TFile*           fFile;           //! galice file with vertices
  TTree*           fTree;           //! tree with headers
  AliHeader*       fHeader;         //! event header
  Int_t            fEventsPerEntry; // number of events with same vertex
  Int_t            fEvent;          //! current event number

  ClassDef(AliVertexGenFile, 1)     // generator for vertices taken from a file
};

#endif














