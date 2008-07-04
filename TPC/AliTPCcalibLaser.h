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
class AliESDEvent;
class AliESDfriend;
class TGraphErrors;


class AliTPCcalibLaser:public AliTPCcalibBase {
public:
  AliTPCcalibLaser();
  AliTPCcalibLaser(const Text_t *name, const Text_t *title);
  virtual ~AliTPCcalibLaser();
  virtual void     Process(AliESDEvent *event);
  virtual void Analyze();
  virtual void DumpMeanInfo(Float_t bfield);
  //
  //
  virtual void DumpLaser(Int_t id);
  virtual void RefitLaser(Int_t id);
  void         FitDriftV();
  void         MakeDistHisto();
  Int_t  FindMirror(AliESDtrack *track, AliTPCseed *seed);
  Bool_t AcceptLaser(Int_t id);
  
  AliESDEvent  * fESD;             //! ESD event  - not OWNER
  AliESDfriend * fESDfriend;       //! ESD event  - not OWNER
  TObjArray      fTracksMirror;    //! tracks with mirror information
  TObjArray      fTracksEsd;       //! tracks with reconstructed information - 
  //                               not owner ESD
  TObjArray      fTracksEsdParam;  //! tracks with reconstructed information - 
  //                               is owner ESD at mirror
  TObjArray      fTracksTPC;       //! tracks with reconstructed information - TPC
  //
  TObjArray      fDeltaZ;          //-> array of histograms of delta z for each track
  TObjArray      fDeltaPhi;        //-> array of histograms of delta z for each track
  TObjArray      fDeltaPhiP;       //-> array of histograms of delta z for each track
  TObjArray      fSignals;         //->Array of dedx signals
  TVectorD*      fFitAside;        //! drift fit - A side
  TVectorD*      fFitCside;        //! drift fit - C- side
  //
  Int_t          fRun;             // current run number
private:
  ClassDef(AliTPCcalibLaser,1)
};





#endif
