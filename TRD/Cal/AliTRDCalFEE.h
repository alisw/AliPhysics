#ifndef AliTRDCALFEE_H
#define AliTRDCALFEE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDCalFEE.h 18952 2007-06-08 11:36:12Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD FEE parameters                             //
//  Empty dummy class for backward compability                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class TString;

class AliTRDCalFEE : public TNamed {

 public:

  AliTRDCalFEE();
  AliTRDCalFEE(const Text_t *name, const Text_t *title);
  virtual ~AliTRDCalFEE() { };
    
 protected:

  ClassDef(AliTRDCalFEE,2)         //  TRD calibration class for TRD FEE parameters

};
#endif
