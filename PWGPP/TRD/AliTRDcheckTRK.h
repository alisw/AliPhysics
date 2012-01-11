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

/*template <typename Value> class TVectorT;
typedef struct TVectorT<Double_t> TVectorD;*/
class TObjArray;
class AliTRDtrackV1;
class AliTRDcheckTRK : public AliTRDresolution
{
public:
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
  static Bool_t   HasClRecalibrate()                   { return fgClRecalibrate;}
  static Bool_t   HasKalmanUpdate()                    { return fgKalmanUpdate;}
  static Bool_t   HasTrkltRefit()                      { return fgTrkltRefit;}
  virtual TObjArray*  Histos();
  TH1*            PlotTrack(const AliTRDtrackV1 *t=NULL);
  TH1*            DoRoads(const AliTRDtrackV1 *t=NULL);
  static Bool_t   PropagateKalman(AliTRDtrackV1 &t, AliExternalTrackParam *ref);
  static void     SetKalmanStep(Float_t step)          { fgKalmanStep=step;}
  static void     SetClRecalibrate(Bool_t set=kTRUE)   { fgClRecalibrate=set;}
  static void     SetKalmanUpdate(Bool_t set=kTRUE)    { fgKalmanUpdate=set;}
  static void     SetTrkltRefit(Bool_t set=kTRUE)      { fgTrkltRefit=set;}

private:
  AliTRDcheckTRK(const AliTRDcheckTRK&);
  AliTRDcheckTRK& operator=(const AliTRDcheckTRK&);
  void     MakePtCalib(Float_t pt0=0.3, Float_t dpt=0.002);
  Int_t    GetPtBinCalib(Float_t pt);

  // kalman related settings
  static Bool_t  fgKalmanUpdate;  // update Kalman with TRD point
  static Bool_t  fgTrkltRefit;    // refit tracklet
  static Bool_t  fgClRecalibrate; // recalibrate clusters and recalculate tracklet fit
  static Float_t fgKalmanStep;    // Kalman stepping
  Float_t        fPtBinCalib[kNptBins+1];  //! pt segmentation

  ClassDef(AliTRDcheckTRK, 1) // TRD tracker systematic
};

#endif
