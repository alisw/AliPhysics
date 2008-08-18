#ifndef ALITPCCALIBCALIB_H
#define ALITPCCALIBCALIB_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////
////
////

#include "AliTPCcalibBase.h"
class AliTPCseed;
class AliESDEvent;
class AliESDtrack;
class TCollection;
class TTreeSRedirector;
class AliExternalTrackParam;
class AliTPCclusterMI;

class AliTPCcalibCalib:public AliTPCcalibBase {
public:
  AliTPCcalibCalib(); 
  AliTPCcalibCalib(const Text_t *name, const Text_t *title);
  AliTPCcalibCalib(const AliTPCcalibCalib&calib);
  AliTPCcalibCalib &operator=(const AliTPCcalibCalib&calib);
  virtual ~AliTPCcalibCalib();
  virtual void     Process(AliESDEvent *event);
  virtual void     Analyze(){return;}
  
  Bool_t  RefitTrack(AliESDtrack * track, AliTPCseed *seed);
  Bool_t  RejectCluster(AliTPCclusterMI* cl, AliExternalTrackParam * param);
protected: 
private:
  ClassDef(AliTPCcalibCalib,1)
};

#endif
