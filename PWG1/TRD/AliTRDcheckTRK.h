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
     kNbunchCross = 3  // no of classes for bunch crossing
    ,kNpt         = 24 // no of log bins in pt spectrum
    ,kNcharge     = 2
  };
  enum ETRDcheckTRKclasses {
     kEntry = 0
    ,kPropagation
    ,kNclasses
  };
  AliTRDcheckTRK();
  AliTRDcheckTRK(char* name);
  virtual ~AliTRDcheckTRK();
  static Float_t  GetKalmanStep()                      { return fgKalmanStep;}
  static Float_t* GetCalib(Int_t det)                  { return fgCalib[det];}
  Int_t           GetSpeciesByMass(Float_t m);
  Int_t           GetPtBin(Float_t pt);
  static Bool_t   HasClRecalibrate()                   { return fgClRecalibrate;}
  static Bool_t   HasKalmanUpdate()                    { return fgKalmanUpdate;}
  static void     LoadCalib(Float_t *calib)            { memcpy(fgCalib, calib, 540*2*sizeof(Float_t));}
  TH1*            PlotTrack(const AliTRDtrackV1 *t=NULL);
  static Bool_t   PropagateKalman(AliTRDtrackV1 &t);
  static void     SetKalmanStep(Float_t step)          { fgKalmanStep=step;}
  static void     SetClRecalibrate(Bool_t set=kTRUE)   { fgClRecalibrate=set;}
  static void     SetKalmanUpdate(Bool_t set=kTRUE)    { fgKalmanUpdate=set;}

private:
  AliTRDcheckTRK(const AliTRDcheckTRK&);
  AliTRDcheckTRK& operator=(const AliTRDcheckTRK&);

  // kalman related settings
  static Bool_t  fgKalmanUpdate;  // update Kalman with TRD point
  static Bool_t  fgClRecalibrate; // recalibrate clusters and recalculate tracklet fit
  static Float_t fgKalmanStep;    // Kalman stepping
  static Float_t fgCalib[540][2]; //! test new calibration params [method][detector][vd/exb]

  ClassDef(AliTRDcheckTRK, 1) // TRD tracker systematic
};

#endif
