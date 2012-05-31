/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Transforms clusters into space points with calibrated positions       //
//  defined in the local tracking system                                  //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TGeoMatrix.h>

#include "AliLog.h"

#include "AliTRDtransform.h"
#include "AliTRDcluster.h"
#include "AliTRDgeometry.h"
#include "AliTRDpadPlane.h"
#include "AliTRDCommonParam.h"
#include "AliTRDcalibDB.h"
#include "Cal/AliTRDCalDet.h"
#include "Cal/AliTRDCalROC.h"

ClassImp(AliTRDtransform)

//_____________________________________________________________________________
AliTRDtransform::AliTRDtransform()
  :TObject()
  ,fDetector(0)
  ,fParam(0x0)
  ,fCalibration(0x0)
  ,fCalVdriftROC(0x0)
  ,fCalT0ROC(0x0)
  ,fCalPRFROC(0x0)
  ,fkCalVdriftDet(0x0)
  ,fkCalT0Det(0x0)
  ,fkCalExBDet(0x0)
  ,fCalVdriftDetValue(0)
  ,fCalT0DetValue(0)
  ,fCalExBDetValue(0)
  ,fSamplingFrequency(0)
  ,fPadPlane(0x0)
  ,fZShiftIdeal(0)
  ,fMatrix(0x0)
{
  //
  // AliTRDtransform default constructor
  //

  fParam             = AliTRDCommonParam::Instance();
  if (!fParam) {
    AliError("Could not get common parameters\n");
  }
  fSamplingFrequency = fParam->GetSamplingFrequency();

  fCalibration       = AliTRDcalibDB::Instance();
  if (!fCalibration) {
    AliError("Cannot find calibration object");
  }

  // Get the calibration objects for the global calibration
  fkCalVdriftDet     = fCalibration->GetVdriftDet();
  fkCalT0Det         = fCalibration->GetT0Det();
  fkCalExBDet        = fCalibration->GetExBDet();

}

//_____________________________________________________________________________
AliTRDtransform::AliTRDtransform(Int_t det)
  :TObject()
  ,fDetector(0)
  ,fParam(0x0)
  ,fCalibration(0x0)
  ,fCalVdriftROC(0x0)
  ,fCalT0ROC(0x0)
  ,fCalPRFROC(0x0)
  ,fkCalVdriftDet(0x0)
  ,fkCalT0Det(0x0)
  ,fkCalExBDet(0x0)
  ,fCalVdriftDetValue(0)
  ,fCalT0DetValue(0)
  ,fCalExBDetValue(0)
  ,fSamplingFrequency(0)
  ,fPadPlane(0x0)
  ,fZShiftIdeal(0)
  ,fMatrix(0x0)
{
  //
  // AliTRDtransform constructor for a given detector
  //

  fParam             = AliTRDCommonParam::Instance();
  if (!fParam) {
    AliError("Could not get common parameters\n");
  }
  fSamplingFrequency = fParam->GetSamplingFrequency();

  fCalibration       = AliTRDcalibDB::Instance();
  if (!fCalibration) {
    AliError("Cannot find calibration object");
  }

  // Get the calibration objects for the global calibration
  fkCalVdriftDet     = fCalibration->GetVdriftDet();
  fkCalT0Det         = fCalibration->GetT0Det();
  fkCalExBDet        = fCalibration->GetExBDet();

  SetDetector(det);

}

//_____________________________________________________________________________
AliTRDtransform::AliTRDtransform(const AliTRDtransform &t)
  :TObject(t)
  ,fDetector(t.fDetector)
  ,fParam(0x0)
  ,fCalibration(0x0)
  ,fCalVdriftROC(0x0)
  ,fCalT0ROC(0x0)
  ,fCalPRFROC(0x0)
  ,fkCalVdriftDet(0x0)
  ,fkCalT0Det(0x0)
  ,fkCalExBDet(0x0)
  ,fCalVdriftDetValue(0)
  ,fCalT0DetValue(0)
  ,fCalExBDetValue(0)
  ,fSamplingFrequency(0)
  ,fPadPlane(0x0)
  ,fZShiftIdeal(0)
  ,fMatrix(0x0)
{
  //
  // AliTRDtransform copy constructor
  //

  fParam             = AliTRDCommonParam::Instance();
  if (!fParam) {
    AliError("Could not get common parameters\n");
  }
  fSamplingFrequency = fParam->GetSamplingFrequency();

  fCalibration = AliTRDcalibDB::Instance();
  if (!fCalibration) {
    AliError("Cannot find calibration object");
  }
  fkCalVdriftDet     = fCalibration->GetVdriftDet();
  fkCalT0Det         = fCalibration->GetT0Det();
  fkCalExBDet        = fCalibration->GetExBDet();
}

//_____________________________________________________________________________
AliTRDtransform::~AliTRDtransform()
{
  //
  // AliTRDtransform destructor
  //

}

//_____________________________________________________________________________
AliTRDtransform &AliTRDtransform::operator=(const AliTRDtransform &t)
{
  //
  // Assignment operator
  //

  if (this != &t) {
    ((AliTRDtransform &) t).Copy(*this);
  }

  return *this;

}

//_____________________________________________________________________________
void AliTRDtransform::Copy(TObject &t) const
{
  //
  // Copy function
  //

  ((AliTRDtransform &) t).fDetector          = fDetector;
  ((AliTRDtransform &) t).fParam             = AliTRDCommonParam::Instance();
  ((AliTRDtransform &) t).fCalibration       = AliTRDcalibDB::Instance();
  if (fCalibration) {
    ((AliTRDtransform &) t).fkCalVdriftDet   = fCalibration->GetVdriftDet();
    ((AliTRDtransform &) t).fkCalT0Det       = fCalibration->GetT0Det();
    ((AliTRDtransform &) t).fkCalExBDet      = fCalibration->GetExBDet();
  }
  else {
    ((AliTRDtransform &) t).fkCalVdriftDet   = 0;
    ((AliTRDtransform &) t).fkCalT0Det       = 0;
    ((AliTRDtransform &) t).fkCalExBDet      = 0;
  }
  ((AliTRDtransform &) t).fCalVdriftROC      = 0x0;
  ((AliTRDtransform &) t).fCalT0ROC          = 0x0;
  ((AliTRDtransform &) t).fCalPRFROC         = 0x0;
  ((AliTRDtransform &) t).fCalVdriftDetValue = 0; 
  ((AliTRDtransform &) t).fCalT0DetValue     = 0;
  ((AliTRDtransform &) t).fCalExBDetValue     = 0;
  ((AliTRDtransform &) t).fSamplingFrequency = 0;
  ((AliTRDtransform &) t).fPadPlane          = 0x0;
  ((AliTRDtransform &) t).fZShiftIdeal       = 0;
  ((AliTRDtransform &) t).fMatrix            = 0x0;

}

//_____________________________________________________________________________
AliTRDgeometry& AliTRDtransform::Geometry()
{
  static AliTRDgeometry g;
  if (!g.CreateClusterMatrixArray()) {
    AliErrorGeneral("AliTRDtransform::Geometry()", "Could not get transformation matrices\n");
  }
  return g;
}

//_____________________________________________________________________________
void AliTRDtransform::SetDetector(Int_t det)
{
  //
  // Set to a new detector number and update the calibration objects
  // and values accordingly
  //

  fDetector          = det;

  // Get the calibration objects for the pad-by-pad calibration
  fCalVdriftROC      = fCalibration->GetVdriftROC(det);
  fCalT0ROC          = fCalibration->GetT0ROC(det);
  fCalPRFROC         = fCalibration->GetPRFROC(det);

  // Get the detector wise defined calibration values
  fCalVdriftDetValue = fkCalVdriftDet->GetValue(det);
  fCalT0DetValue     = fkCalT0Det->GetValue(det);
  fCalExBDetValue    = fkCalExBDet->GetValue(det);

  // Shift needed to define Z-position relative to middle of chamber
  Int_t layer        = Geometry().GetLayer(det);
  Int_t stack        = Geometry().GetStack(det);
  fPadPlane          = Geometry().GetPadPlane(layer,stack);
  fZShiftIdeal       = 0.5 * (fPadPlane->GetRow0() + fPadPlane->GetRowEnd());

  // Get the current transformation matrix
  fMatrix            = Geometry().GetClusterMatrix(det);

}

//_____________________________________________________________________________
Bool_t AliTRDtransform::Transform(AliTRDcluster *c)
{
  //
  // Transforms the local cluster coordinates into calibrated 
  // space point positions defined in the local tracking system.
  //
  // Here the calibration for T0, Vdrift and ExB is applied as well.
  //
  // Input: Cluster in the local chamber coordinates
  // Output: Tracking cluster

  if (!fMatrix) return kFALSE;


  // Parameters to adjust the X position of clusters in the alignable volume
  const Double_t kX0shift = AliTRDgeometry::AnodePos(); //[cm]

 
  // Retrieve calibration values
  Int_t col = c->GetPadCol(), row = c->GetPadRow();
  // drift velocity
  Double_t vd  = fCalVdriftDetValue * fCalVdriftROC->GetValue(col,row);
  // t0
  Double_t t0  = fCalT0DetValue     + fCalT0ROC->GetValue(col,row);
  t0 /= fSamplingFrequency;
  // ExB correction
  Double_t exb = fCalExBDetValue;//AliTRDCommonParam::Instance()->GetOmegaTau(vd);

  Float_t x = c->GetXloc(t0, vd);

  // Pad dimensions
  Double_t rs = fPadPlane->GetRowSize(row);
  Double_t cs = fPadPlane->GetColSize(col);

  // cluster error with diffusion corrections
  Double_t s2  = cs*fCalPRFROC->GetValue(col, row); 
  s2 *= s2; 
  Float_t dl, dt;
  AliTRDCommonParam::Instance()->GetDiffCoeff(dl, dt, vd);

  Double_t y0 = fPadPlane->GetColPos(col) + .5*cs;
  Double_t loc[] = {
    kX0shift-x,                    // Invert the X-position,
    c->GetYloc(y0, s2, cs) - x*exb,// apply ExB correction
    fPadPlane->GetRowPos(row) - .5*rs - fZShiftIdeal // move the Z-position relative to the middle of the chamber
  };

  // Go to tracking coordinates
  Double_t trk[3];
  fMatrix->LocalToMaster(loc, trk);

  // store tracking values
  c->SetX(trk[0]);c->SetY(trk[1]);c->SetZ(trk[2]);
  c->SetSigmaY2(s2, dt, exb, x);
  c->SetSigmaZ2(fPadPlane->GetRowSize(row)*fPadPlane->GetRowSize(row)/12.);
  
  return kTRUE;

}

//_____________________________________________________________________________
void AliTRDtransform::Recalibrate(AliTRDcluster *c, Bool_t setDet)
{
  //
  // Recalibrates the position of a given cluster
  // If <setDet> is TRUE, the detector number is set for each cluster
  // automatically. Otherwise, AliTRDtransform::SetDetector() has to
  // be used.
  //

  if (setDet) SetDetector(c->GetDetector());
  Transform(c);

}
