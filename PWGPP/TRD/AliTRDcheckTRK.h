#ifndef ALITRDCHECKTRK_H
#define ALITRDCHECKTRK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  TRD tracker systematic                                                //
//                                                                        //
//  Authors:                                                              //
//    Alexandru Bercuci <A.Bercuci@gsi.de>                                //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ALITRDRESOLUTION_H
#include "AliTRDresolution.h"
#endif

class TObjArray;
class AliTRDtrackV1;
class AliTRDcheckTRK : public AliTRDresolution
{
public:
  enum ETRDcheckTRKsteer {
     kKalmanUpdate  = 0
    ,kTrkltRefit
    ,kClRecalibrate
    ,kUseITS
  };
  enum ETRDcheckTRKconst {
     kNptBins    = 25 // no of log bins in pt spectrum
    ,kNSigmaBins = 25 // no of sigma bins
    ,kNclusters  = 19 // no of no of clusters
    ,kNdim       = 7  // no of dimensions in THnSparse
  };
  AliTRDcheckTRK();
  AliTRDcheckTRK(char* name);
  virtual ~AliTRDcheckTRK();
  static Float_t  GetKalmanStep()                      { return fgKalmanStep;}
  static Bool_t   HasClRecalibrate()                   { return TESTBIT(fgSteer, kClRecalibrate);}
  static Bool_t   HasKalmanUpdate()                    { return TESTBIT(fgSteer, kKalmanUpdate);}
  static Bool_t   HasTrkltRefit()                      { return TESTBIT(fgSteer, kTrkltRefit);}
  virtual TObjArray*  Histos();
  TH1*            PlotTrack(const AliTRDtrackV1 *t=NULL);
  TH1*            DoRoads(const AliTRDtrackV1 *t=NULL);
  static Bool_t   PropagateKalman(AliTRDtrackV1 &t, AliExternalTrackParam *ref);
  static void     SetKalmanStep(Float_t step)          { fgKalmanStep=step;}
  static void     SetClRecalibrate(Bool_t s=kTRUE)     { if(s){SETBIT(fgSteer, kClRecalibrate); SetTrkltRefit();} else CLRBIT(fgSteer, kClRecalibrate);}
  static void     SetKalmanUpdate(Bool_t s=kTRUE)      { if(s) SETBIT(fgSteer, kKalmanUpdate); else CLRBIT(fgSteer, kKalmanUpdate);}
  static void     SetTrkltRefit(Bool_t s=kTRUE)        { if(s) SETBIT(fgSteer, kTrkltRefit); else CLRBIT(fgSteer, kTrkltRefit);}
  static void     SetUseITS(Bool_t s=kTRUE)            { if(s) SETBIT(fgSteer, kUseITS); else CLRBIT(fgSteer, kUseITS);}
  static Bool_t   UseITS()                             { return TESTBIT(fgSteer, kUseITS);}
private:
  AliTRDcheckTRK(const AliTRDcheckTRK&);
  AliTRDcheckTRK& operator=(const AliTRDcheckTRK&);
  void     MakePtCalib(Float_t pt0=0.3, Float_t dpt=0.002);
  Int_t    GetPtBinCalib(Float_t pt);

  // kalman related settings
  static UChar_t fgSteer;         // steering bit map
  static Float_t fgKalmanStep;    // Kalman stepping
  Float_t        fPtBinCalib[kNptBins+1];  //! pt segmentation

  ClassDef(AliTRDcheckTRK, 1) // TRD tracker systematic
};

#endif
