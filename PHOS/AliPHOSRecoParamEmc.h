#ifndef ALIPHOSRECOPARAMEMC_H
#define ALIPHOSRECOPARAMEMC_H
/* Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

/* $Id$ */
                                              
// This class contains the PHOS EMC reconstruction parameters.
// See cxx source for use case.

#include "AliPHOSRecoParam.h"


class AliPHOSRecoParamEmc : public AliPHOSRecoParam {

public:

  AliPHOSRecoParamEmc();
  virtual ~AliPHOSRecoParamEmc() {}

  static AliPHOSRecoParam* GetEmcDefaultParameters();

  ClassDef(AliPHOSRecoParamEmc,1)
};

#endif
