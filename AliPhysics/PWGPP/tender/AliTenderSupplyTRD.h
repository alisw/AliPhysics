#ifndef ALITENDERSUPPLYTRD_H
#define ALITENDERSUPPLYTRD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTenderSupplyTRD.h 34521 2009-09-01 13:21:10Z abercuci $ */
// Author: Alexandru Bercuci, 15/11/2009

////////////////////////////////////////////////////////////
//
// Correct ESD TRD info based on up-to-date calibration
//
////////////////////////////////////////////////////////////


#ifndef ALITENDERSUPPLY_H
#include "AliTenderSupply.h"
#endif

class AliTenderSupplyTRD : public AliTenderSupply
{
public:
  AliTenderSupplyTRD();

  virtual   ~AliTenderSupplyTRD(){}
  void      Init();
  void      ProcessEvent();

private:
  AliTenderSupplyTRD(const AliTenderSupplyTRD &ref);
  AliTenderSupplyTRD& operator=(const AliTenderSupplyTRD &ref);

  ClassDef(AliTenderSupplyTRD,1)  // tender for TRD detector
};

#endif

