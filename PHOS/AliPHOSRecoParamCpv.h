#ifndef ALIPHOSRECOPARAMCPV_H
#define ALIPHOSRECOPARAMCPV_H
/* Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

/* $Id$ */
                                              
// This class contains the PHOS CPV reconstruction parameters.
// See cxx source for use case.

#include "AliPHOSRecoParam.h"


class AliPHOSRecoParamCpv : public AliPHOSRecoParam {

public:

  AliPHOSRecoParamCpv();
  virtual ~AliPHOSRecoParamCpv() {}

  static AliPHOSRecoParam* GetCpvDefaultParameters();

  ClassDef(AliPHOSRecoParamCpv,1)
};

#endif
