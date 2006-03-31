#ifndef ALITOFALIGNMENT_H
#define ALITOFALIGNMENT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////
//  class for TOF Alignment::                                   //
//////////////////////////////////////////////////////////////////

#include "TTask.h"
#include "AliAlignObj.h"

class AliTOFAlignment :public TTask{

enum {kMaxAlignObj=2000}; //maximal number of the TOF Alignable Objects

public:

 AliTOFAlignment(); 
 AliTOFAlignment(const AliTOFAlignment &t); //Copy Ctor 
  virtual ~AliTOFAlignment() {delete fTOFAlignObjArray;}
  virtual void WriteParOnCDB(Char_t *sel, Int_t minrun, Int_t maxrun);
  virtual void ReadParFromCDB(Char_t *sel, Int_t nrun);
  virtual void WriteSimParOnCDB(Char_t *sel, Int_t minrun, Int_t maxrun);
  virtual void ReadSimParFromCDB(Char_t *sel, Int_t nrun);
  virtual void Smear(Float_t *tr=0, Float_t *rot=0); // create a set of AlignObj for TOF
  virtual void Align(Float_t *tr=0, Float_t *rot=0); // create a set of AlignObj for TOF
  virtual void WriteOnCDBforDC();
  virtual void ReadFromCDBforDC();
  TObjArray * GetTOFAlignArray() const {return fTOFAlignObjArray;}

private:

  Int_t fNTOFAlignObj;      // Number of Alignable Objects
  TObjArray *fTOFAlignObjArray;
  ClassDef(AliTOFAlignment,1) // TOF Alignment 
};

#endif

