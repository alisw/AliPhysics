#ifndef AliMFTRecoParam_H
#define AliMFTRecoParam_H

/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//====================================================================================================================================================
//
//      Class with MFT recosntruction parameters
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "AliDetectorRecoParam.h"

//====================================================================================================================================================

class AliMFTRecoParam : public AliDetectorRecoParam {

public: 

  AliMFTRecoParam();
  virtual ~AliMFTRecoParam() {;}
  
private:

  AliMFTRecoParam(const AliMFTRecoParam &param);
  AliMFTRecoParam &operator=(const AliMFTRecoParam &param);

  ClassDef(AliMFTRecoParam,1) // MFT reco parameters

};

//====================================================================================================================================================

#endif
