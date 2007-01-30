#ifndef AliTRDCALGLOBALS_H
#define AliTRDCALGLOBALS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for global TRD parameters                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class AliTRDCalGlobals : public TNamed {

 public:

  AliTRDCalGlobals();
  AliTRDCalGlobals(const Text_t *name, const Text_t *title);
  virtual ~AliTRDCalGlobals() {};
    
  void    SetNumberOfTimeBins(Int_t value)   { fNumberOfTimeBins = value; };
  Int_t   GetNumberOfTimeBins() const        { return fNumberOfTimeBins;  };
  
 protected:

  Int_t   fNumberOfTimeBins;       // Number of timebins  
    
  ClassDef(AliTRDCalGlobals,2)     // TRD calibration class for global TRD parameters

};

#endif
