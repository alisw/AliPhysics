#ifndef ALITRDTRANSFORM_H
#define ALITRDTRANSFORM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Transforms clusters into space points with calibrated positions       //
//  defined in the local tracking system                                  //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TObject.h"

class TGeoHMatrix;

class AliTRDgeometry;
class AliTRDcluster;
class AliTRDCommonParam;
class AliTRDcalibDB;
class AliTRDCalROC;
class AliTRDCalDet;
class AliTRDpadPlane;

class AliTRDtransform : public TObject {

 public:

  AliTRDtransform();
  AliTRDtransform(Int_t det);
  AliTRDtransform(const AliTRDtransform &t);
  virtual ~AliTRDtransform();
  AliTRDtransform &operator=(const AliTRDtransform &t);
  
  virtual void     Copy(TObject &t) const;
  AliTRDpadPlane*  GetPadPlane() const {return fPadPlane;}
  virtual Bool_t   Transform(AliTRDcluster *c);
  virtual void     Recalibrate(AliTRDcluster *c, Bool_t setDet = kTRUE);

          void     SetDetector(Int_t det);
  static  AliTRDgeometry& Geometry(); 

  protected:
  Int_t               fDetector;            //  Detector number

  AliTRDCommonParam  *fParam;               //  TRD common parameters

  AliTRDcalibDB      *fCalibration;         //  TRD calibration interface object
  AliTRDCalROC       *fCalVdriftROC;        //  Pad wise Vdrift calibration object
  AliTRDCalROC       *fCalT0ROC;            //  Pad wise T0 calibration object
  AliTRDCalROC       *fCalPRFROC;           //  Pad wise PRF calibration object
  const AliTRDCalDet *fkCalVdriftDet;       //  ROC wise Vdrift calibration object
  const AliTRDCalDet *fkCalT0Det;           //  ROC wise T0 calibration object
  const AliTRDCalDet *fkCalExBDet;          //  ROC wise ExB calibration object
  Double_t            fCalVdriftDetValue;   //  ROC wise Vdrift calibration value
  Double_t            fCalT0DetValue;       //  ROC wise T0 calibration value
  Double_t            fCalExBDetValue;      //  Det wise ExB calibration value

  Double_t            fSamplingFrequency;   //  ADC sampling frequency

  AliTRDpadPlane     *fPadPlane;            //  The current pad plane object
  Double_t            fZShiftIdeal;         //  Needed to define Z-position relative to middle of chamber

  TGeoHMatrix        *fMatrix;              //  Transformation matrix for a given chamber

  ClassDef(AliTRDtransform, 3)              //  Transforms clusters

};
#endif
