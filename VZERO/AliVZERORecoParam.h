#ifndef ALIVZERORECOPARAM_H
#define ALIVZERORECOPARAM_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class with VZERO reconstruction parameters                                //
// Origin: Brigitte Cheynis                                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliDetectorRecoParam.h"

class AliVZERORecoParam : public AliDetectorRecoParam
{
 public: 
  AliVZERORecoParam();
  virtual ~AliVZERORecoParam();

 private:

  AliVZERORecoParam(const AliVZERORecoParam & param);
  AliVZERORecoParam & operator=(const AliVZERORecoParam &param);

  ClassDef(AliVZERORecoParam,1) // VZERO reco parameters
};

#endif
