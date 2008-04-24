#ifndef ALITPCCALIBBASE_H
#define ALITPCCALIBBASE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////
////
////

#include "TNamed.h"
class AliTPCseed;
class AliESDevent;
class TCollection;

class AliTPCcalibBase:public TNamed {
public:
  AliTPCcalibBase(); 
  virtual ~AliTPCcalibBase();
  virtual void     Process(AliESDevent */*event*/){return;}
  virtual void     Process(AliTPCseed */*track*/){return;}
  virtual Long64_t Merge(TCollection */*li*/){return 0;}
  virtual void     Analyze(){return;}
 private: 
  ClassDef(AliTPCcalibBase,1)
};

#endif
