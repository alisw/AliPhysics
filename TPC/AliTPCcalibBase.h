#ifndef ALITPCCALIBBASE_H
#define ALITPCCALIBBASE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////
////
////

#include "TNamed.h"
#include "TObjString.h"
class AliTPCseed;
class AliESDEvent;
class AliESDtrack;
class TCollection;
class TTreeSRedirector;

class AliTPCcalibBase:public TNamed {
public:
  AliTPCcalibBase(); 
  AliTPCcalibBase(const char * name, const char * title); 
  AliTPCcalibBase(const AliTPCcalibBase&calib);
  AliTPCcalibBase &operator=(const AliTPCcalibBase&calib);
  virtual ~AliTPCcalibBase();
  virtual void     Process(AliESDEvent */*event*/){return;}
  virtual void     Process(AliTPCseed */*track*/){return;}
  virtual void     Process(AliESDtrack */*track*/, Int_t /*runNo=-1*/){return;}
  virtual Long64_t Merge(TCollection */*li*/){return 0;}
  virtual void     Analyze(){return;}
  virtual void     Terminate();
  virtual void     UpdateEventInfo(AliESDEvent * event);
  virtual Bool_t   AcceptTrigger();
  virtual void     SetTriggerMask(Int_t accept, Int_t reject, Bool_t rejectLaser){fTriggerMaskAccept=accept;fTriggerMaskReject=reject; fRejectLaser = rejectLaser;}
 
  //
  // debug streamer support
  TTreeSRedirector *GetDebugStreamer();
  void       SetStreamLevel(Int_t streamLevel){fStreamLevel=streamLevel;}
  void       SetDebugLevel(Int_t level) {fDebugLevel = level;}
  Int_t      GetStreamLevel() const {return fStreamLevel;}
  Int_t      GetDebugLevel() const {return fDebugLevel;}
  virtual void RegisterDebugOutput(const char *path);
  static     Bool_t HasLaser(AliESDEvent *event);

protected: 
  TTreeSRedirector *fDebugStreamer;     //! debug streamer
  Int_t  fStreamLevel;                  //  debug stream level
  Int_t  fRun;                          //!  current Run number
  Int_t  fEvent;                        //! current Event number
  Int_t  fTime;                         //!  current Time
  ULong64_t  fTrigger;                  //! current trigger mask
  Float_t fMagF;                        //! current magnetic field 
  Int_t   fTriggerMaskReject;           //trigger mask - non accept trigger
  Int_t   fTriggerMaskAccept;           //trigger mask - accept
  Bool_t  fHasLaser;                    //flag the laser is overlayed with given event
  Bool_t  fRejectLaser;                 //flag- reject laser
  TObjString fTriggerClass;                // trigger class
private:
  Int_t  fDebugLevel;                   //  debug level

  ClassDef(AliTPCcalibBase,1)
};

#endif
