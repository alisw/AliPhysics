#ifndef ALITOFT0FILL_H
#define ALITOFT0FILL_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */

// *
// *
// *
// * this class defines the T0Fill object to be stored
// * in OCDB in order to apply T0Fill correction during 
// * reconstruction. 
// *
// *
// *

#include "TObject.h"

class AliTOFT0Fill :
public TObject
{

 public:

  AliTOFT0Fill(); // default constructor
  virtual ~AliTOFT0Fill(); // default destructor
  AliTOFT0Fill(const AliTOFT0Fill &source); // copy constructor
  AliTOFT0Fill &operator=(const AliTOFT0Fill &source); // operator=
  Float_t GetT0Fill() const {return fT0Fill;}; // getter
  void SetT0Fill(Float_t value) {fT0Fill = value;}; // setter

 private:

  Float_t fT0Fill; // event time (ps)

  ClassDef(AliTOFT0Fill, 1);
};

#endif /* ALITOFT0FILL_H */
