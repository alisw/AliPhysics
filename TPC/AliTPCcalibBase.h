#ifndef ALITPCCALIBBASE_H
#define ALITPCCALIBBASE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////
////
////

#include "TNamed.h"
class AliTPCseed;
class AliESDEvent;
class TCollection;
class TTreeSRedirector;

class AliTPCcalibBase:public TNamed {
public:
  AliTPCcalibBase(); 
  AliTPCcalibBase(const AliTPCcalibBase&calib);
  AliTPCcalibBase &operator=(const AliTPCcalibBase&calib);
  virtual ~AliTPCcalibBase();
  virtual void     Process(AliESDEvent */*event*/){return;}
  virtual void     Process(AliTPCseed */*track*/){return;}
  virtual Long64_t Merge(TCollection */*li*/){return 0;}
  virtual void     Analyze(){return;}
  virtual void     Terminate();
  //
  // debug streamer support
  TTreeSRedirector *GetDebugStreamer();
  void       SetStreamLevel(Int_t streamLevel){fStreamLevel=streamLevel;}
  void       SetDebugLevel(Int_t level) {fDebugLevel = level;}
  Int_t      GetStreamLevel() const {return fStreamLevel;}
  Int_t      GetDebugLevel() const {return fDebugLevel;}
protected: 
  TTreeSRedirector *fDebugStreamer;     //! debug streamer
  Int_t  fStreamLevel;                  //  debug stream level
private:
  Int_t  fDebugLevel;                   //  debug level
  ClassDef(AliTPCcalibBase,1)
};

#endif
