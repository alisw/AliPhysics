#ifndef ALITPCCALIBCALIB_H
#define ALITPCCALIBCALIB_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////
////
////

#include "AliTPCcalibBase.h"
class AliTPCseed;
//class AliESDEvent;
class AliVEvent;
//class AliESDtrack;
class AliVTrack;
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
  virtual void     Process(AliVEvent *event);
  virtual void     Analyze(){return;}
  
  Bool_t  RefitTrack(AliVTrack * track, AliTPCseed *seed, Float_t magesd);
  Bool_t  RejectCluster(AliTPCclusterMI* cl, AliExternalTrackParam * param);
  //void     Process(AliESDtrack *track, Int_t runNo=-1){AliTPCcalibBase::Process(track,runNo);}
  void     Process(AliVTrack *track, Int_t runNo=-1){AliTPCcalibBase::Process(track,runNo);}
  void     Process(AliTPCseed *track){return AliTPCcalibBase::Process(track);}

  void  SetApplyExBCorrection(Int_t flag){fApplyExBCorrection=flag;}
  void  SetApplyTOFCorrection(Int_t flag){fApplyTOFCorrection=flag;}
  void  SetApplyPositionCorrection(Int_t flag){fApplyPositionCorrection=flag;}
  void  SetApplySectorAlignment(Int_t flag){fApplySectorAlignment=flag;}
  void  SetApplyRPhiCorrection(Int_t flag){fApplyRPhiCorrection=flag;}
  void  SetApplyRCorrection(Int_t flag){fApplyRCorrection=flag;}

  //
  Int_t GetApplyExBCorrection() const {return fApplyExBCorrection;}
  Int_t GetApplyTOFCorrection() const {return fApplyTOFCorrection;}
  Int_t GetApplyPositionCorrection() const {return fApplyPositionCorrection;}
  Int_t GetApplySectorAlignment() const {return fApplySectorAlignment;}
  Int_t GetApplyRPhiCorrection() const {return fApplyRPhiCorrection;}
  Int_t GetApplyRCorrection() const {return fApplyRCorrection;}

protected: 
  Int_t fApplyExBCorrection;      // apply ExB correction (in AliTPCTransform)
  Int_t fApplyTOFCorrection;      // apply TOF correction (in AliTPCTransform)
  Int_t fApplyPositionCorrection; // apply position correction
  Int_t fApplySectorAlignment;    // apply sector alignment
  Int_t fApplyRPhiCorrection;     // apply R-Phi correction
  Int_t fApplyRCorrection;        // apply Radial correction
private:
  ClassDef(AliTPCcalibCalib,2)
};

#endif
