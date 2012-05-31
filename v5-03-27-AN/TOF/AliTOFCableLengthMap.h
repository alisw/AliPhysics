#ifndef ALITOFCABLELENGTHMAP_H
#define ALITOFCABLELENGTHMAP_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TOF Cable Length Map class                                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "AliTOFGeometry.h"

class AliTOFCableLengthMap: public TObject{

 public:
  AliTOFCableLengthMap();
  virtual ~AliTOFCableLengthMap();
  static Float_t GetCableLength(Int_t icrate, Int_t islot, Int_t ichain, Int_t itdc);
  static Float_t GetCableTimeShift(Int_t icrate, Int_t islot, Int_t ichain, Int_t itdc);
  static Int_t GetCableTimeShiftBin(Int_t icrate, Int_t islot, Int_t ichain, Int_t itdc);
  static Float_t GetPropagationDelay() {return fgkPropagationDelay;};

 private:
  
  static const Float_t fgkCableLength[72][10][2][5];//Cable Length
  static const Float_t fgkPropagationDelay;// Propagation delay [ns/cm]

  ClassDef(AliTOFCableLengthMap,0) // TOF Cable Length Map class
    };

#endif
