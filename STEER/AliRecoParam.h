#ifndef ALIRECOPARAM_H
#define ALIRECOPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Base Class for Detector reconstruction parameters                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "TNamed.h"
class AliDetectorRecoParam;

class AliRecoParam : public TNamed
{
  enum EventType0 {kUndef=0, kPhysic=1, kCalib=2};  
 public: 
  AliRecoParam();
  virtual ~AliRecoParam();  
  static AliRecoParam * Instance();
  //
  virtual void        Print(Option_t *option="") const;
  TObjArray * GetRecoParam(const char * detType, Int_t *eventType=0);  
  void        RegisterRecoParam(AliDetectorRecoParam* param);
protected:
  TObjArray *fRecoParamArray;   //array with registerd reconstruction parameters
  static AliRecoParam* fgInstance; // Reconstruction parameters instance
  ClassDef(AliRecoParam, 1)
};


#endif
