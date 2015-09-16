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
//class AliESDEvent;
class AliVEvent;
//class AliESDtrack;
class AliVTrack;
//class AliESDfriendTrack;
class AliVfriendTrack;
class TCollection;
class TTreeSRedirector;
class TGraph;
class TGraphErrors;
class THnSparse;
class TH1;
class TH2;
class TAxis;

class AliTPCcalibBase:public TNamed {
public:
  AliTPCcalibBase(); 
  AliTPCcalibBase(const char * name, const char * title); 
  AliTPCcalibBase(const AliTPCcalibBase&calib);
  AliTPCcalibBase &operator=(const AliTPCcalibBase&calib);
  virtual ~AliTPCcalibBase();
  virtual void     Process(AliVEvent *event){
      //Printf("AliTPCCalibBase::Process(ESD event)...");
      fCurrentEvent = event; return;}
  virtual void     Process(AliTPCseed *track){
      //Printf("AliTPCCalibBase::Process(TPC track)...");
      fCurrentSeed = track; return;}
  virtual void     Process(AliVTrack *track, Int_t /*runNo=-1*/){
      //Printf("AliTPCCalibBase::Process(ESD track)...");
      fCurrentTrack=track; return;}
  virtual Long64_t Merge(TCollection */*li*/){return 0;}
  virtual void     Analyze(){return;}
  virtual void     Terminate();
  virtual void     UpdateEventInfo(AliVEvent * event);
  virtual Bool_t   AcceptTrigger();
  virtual void     SetTriggerMask(Int_t accept, Int_t reject, Bool_t rejectLaser){fTriggerMaskAccept=accept;fTriggerMaskReject=reject; fRejectLaser = rejectLaser;}

  // full reset: discard all statistics, zero histograms, start again.
  // called in online mode (HLT) after sending output for merging.
  virtual Bool_t            ResetOutputData() {return kFALSE;}

  //
  // debug streamer support
  TTreeSRedirector *GetDebugStreamer();
  void       SetStreamLevel(Int_t streamLevel){fStreamLevel=streamLevel;}
  void       SetDebugLevel(Int_t level) {fDebugLevel = level;}
  Int_t      GetStreamLevel() const {return fStreamLevel;}
  Int_t      GetDebugLevel() const {return fDebugLevel;}
  virtual void RegisterDebugOutput(const char *path);
  static     Bool_t HasLaser(AliVEvent *event);
  static TGraphErrors *        FitSlices(THnSparse *h, Int_t axisDim1, Int_t axisDim2, Int_t minEntries, Int_t nmaxBin, Float_t fracLow=0.1, Float_t fracUp=0.9, Bool_t useMedian=kFALSE, TTreeSRedirector *cstream=0, Int_t ival=1);
  static TGraphErrors *        FitSlices(TH2* hist, Int_t minEntries, Int_t nmaxBin, Float_t fracLow=0.1, Float_t fracUp=0.9, Bool_t useMedian=kFALSE, TTreeSRedirector *cstream=0, Int_t ival=1);
  static TH2*            NormalizedProjection(THnSparse *h, Int_t axisDim1, Int_t axisDim2, Int_t normDim, Float_t minStatFrac=0.5);
  static void            BinLogX(THnSparse *h, Int_t axisDim);
  static void            BinLogX(TH1 *h);
  static void            BinLogX(TAxis * axis);
  void SetRun(Int_t run){ fRun=run;}
protected: 
  TTreeSRedirector *fDebugStreamer;     //! debug streamer
  Int_t  fStreamLevel;                  //  debug stream level
  Int_t  fRun;                          //!  current Run number
  Int_t  fEvent;                        //! current Event number
  Int_t  fTime;                         //!  current Time
  ULong64_t  fTrigger;                  //! current trigger mask
  Float_t fMagF;                        // current magnetic field 
  Int_t   fTriggerMaskReject;           //trigger mask - non accept trigger
  Int_t   fTriggerMaskAccept;           //trigger mask - accept
  Bool_t  fHasLaser;                    //flag the laser is overlayed with given event
  Bool_t  fRejectLaser;                 //flag- reject laser
  TObjString fTriggerClass;             // trigger class
  AliVEvent  *fCurrentEvent;          //! current event
  AliVTrack *fCurrentTrack;           //! current esd track
  AliVfriendTrack *fCurrentFriendTrack;     //! current friend track
  AliTPCseed   *fCurrentSeed;           //! current seed
private:
  Int_t  fDebugLevel;                   //  debug level

  ClassDef(AliTPCcalibBase,3)
};

#endif
