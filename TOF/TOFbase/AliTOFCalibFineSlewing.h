#ifndef ALITOFCALIBFINESLEWING_H
#define ALITOFCALIBFINESLEWING_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */

// *
// *
// *
// * this class defines the TOF object to be stored
// * in OCDB for the fine time-slewing calibration
// * 
// * 
// *
// *
// *

#include "TObject.h"

class AliTOFCalibFineSlewing : public TObject
{

 public:

  AliTOFCalibFineSlewing(); // default constructor
  AliTOFCalibFineSlewing(TObjArray *oa); // standard constructor
  AliTOFCalibFineSlewing(const AliTOFCalibFineSlewing &src); // copy constructor
  AliTOFCalibFineSlewing &operator=(const AliTOFCalibFineSlewing &src); // operator=
  virtual ~AliTOFCalibFineSlewing(); // default destructor

  Bool_t IsPointValid(Float_t x, Float_t y); // is point valid
  Float_t Eval(Int_t ich, Float_t x); // eval

 protected:

  UInt_t    fSize; // total number of points
  UShort_t *fX; //[fSize] x values
  Short_t  *fY; //[fSize] y values
  UInt_t    fStart[157248]; // index of first entry 

  ClassDef(AliTOFCalibFineSlewing, 1);
};

#endif /* ALITOFCALIBFINESLEWING_H */
