#ifndef ALITPCCALIBLASER_H
#define ALITPCCALIBLASER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////
////
////

#include "TObject.h"
#include "TObjArray.h"
#include "TLinearFitter.h"
#include "AliTPCcalibBase.h"
#include "TH1.h"

class AliExternalTrackParam;
class AliESDtrack;
class TGraphErrors;

class AliTPCcalibLaser:public AliTPCcalibBase {
public:
  AliTPCcalibLaser();
  AliTPCcalibLaser(const Text_t *name, const Text_t *title);
  virtual ~AliTPCcalibLaser();
  virtual void Process(AliESDtrack *track);
  virtual void Analyze();
  virtual void Terminate();  
  //

private:
  ClassDef(AliTPCcalibLaser,1)
};





#endif
