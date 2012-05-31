#ifndef AliTRDCALDCSPTRTlmu_H
#define AliTRDCALDCSPTRTlmu_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDCalDCSPTRTlmu.h 18952 2007-06-08 11:36:12Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD GTU configuration parameters               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class TString;

class AliTRDCalDCSPTRTlmu : public TNamed {

 public:

  AliTRDCalDCSPTRTlmu();
  AliTRDCalDCSPTRTlmu(const char *name, const char *title);
  AliTRDCalDCSPTRTlmu(const AliTRDCalDCSPTRTlmu &);
  virtual ~AliTRDCalDCSPTRTlmu() { };

 protected:


  ClassDef(AliTRDCalDCSPTRTlmu,1)      //  TRD calibration class for TRD GTU parameters

};
#endif
