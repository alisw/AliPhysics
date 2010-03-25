#ifndef ALITOFDELTABCOFFSET_H
#define ALITOFDELTABCOFFSET_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */

// *
// *
// *
// * this class defines the DeltaBCOffset object to be stored
// * in OCDB in order to apply DeltaBC correction during 
// * reconstruction. 
// *
// *
// *

#include "TObject.h"

class AliTOFDeltaBCOffset :
public TObject
{

 public:

  AliTOFDeltaBCOffset(); // default constructor
  virtual ~AliTOFDeltaBCOffset(); // default destructor
  AliTOFDeltaBCOffset(const AliTOFDeltaBCOffset &source); // copy constructor
  AliTOFDeltaBCOffset &operator=(const AliTOFDeltaBCOffset &source); // operator=
  Int_t GetDeltaBCOffset() const {return fDeltaBCOffset;}; // getter
  void SetDeltaBCOffset(Int_t value) {fDeltaBCOffset = value;}; // setter

 private:

  Int_t fDeltaBCOffset; // deltaBC offset (BC bins)

  ClassDef(AliTOFDeltaBCOffset, 1);
};

#endif /* ALITOFDELTABCOFFSET_H */
