#ifndef ALITOFT0FILLONLINECALIB_H
#define ALITOFT0FILLONLINECALIB_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */

// *
// *
// *
// * this class defines the T0Fill online calibration object to be stored
// * in OCDB in order to obtain a corrected T0Fill from online algorithm
// * 
// *
// *
// *

#include "TObject.h"

class AliTOFT0FillOnlineCalib :
public TObject
{

 public:

  AliTOFT0FillOnlineCalib(); // default constructor
  virtual ~AliTOFT0FillOnlineCalib(); // default destructor
  AliTOFT0FillOnlineCalib(const AliTOFT0FillOnlineCalib &source); // copy constructor
  AliTOFT0FillOnlineCalib &operator=(const AliTOFT0FillOnlineCalib &source); // operator=
  Float_t GetOffset() const {return fOffset;}; // getter
  Float_t GetCoefficient() const {return fCoefficient;}; // getter
  void SetOffset(Float_t value) {fOffset = value;}; // setter
  void SetCoefficient(Float_t value) {fCoefficient = value;}; // setter

 private:

  Float_t fOffset; // offset (ps)
  Float_t fCoefficient; // coefficient

  ClassDef(AliTOFT0FillOnlineCalib, 1);
};

#endif /* ALITOFT0FILLOnlineCalib_H */
