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
#include "AliTracker.h"
#include "AliCodeTimer.h"

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
//AliTRDtransform::AliTRDtransform()
//  :AliTransform()
AliTRDtransform::AliTRDtransform()
  :TObject()
  ,fGeo(0x0)
  ,fDetector(0)
  ,fParam(0x0)
  ,fCalibration(0x0)
  ,fCalVdriftROC(0x0)
  ,fCalT0ROC(0x0)
  ,fCalVdriftDet(0x0)
  ,fCalT0Det(0x0)
  ,fCalVdriftDetValue(0)
  ,fCalT0DetValue(0)
  ,fSamplingFrequency(0)
  ,fPadPlane(0x0)
  ,fZShiftIdeal(0)
  ,fMatrix(0x0)
{
  //
  // AliTRDtransform default constructor
  //

}

//_____________________________________________________________________________
//AliTRDtransform::AliTRDtransform(Int_t det)
//  :AliTransform()
AliTRDtransform::AliTRDtransform(Int_t det)
  :TObject()
  ,fGeo(0x0)
  ,fDetector(0)
  ,fParam(0x0)
  ,fCalibration(0x0)
  ,fCalVdriftROC(0x0)
  ,fCalT0ROC(0x0)
  ,fCalVdriftDet(0x0)
  ,fCalT0Det(0x0)
  ,fCalVdriftDetValue(0)
  ,fCalT0DetValue(0)
  ,fSamplingFrequency(0)
  ,fPadPlane(0x0)
  ,fZShiftIdeal(0)
  ,fMatrix(0x0)
{
  //
  // AliTRDtransform constructor for a given detector
  //

  fGeo               = new AliTRDgeometry();
  if (!fGeo->CreateClusterMatrixArray()) {
    AliError("Could not get transformation matrices\n");
  }

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
  fCalVdriftDet      = fCalibration->GetVdriftDet();
  fCalT0Det          = fCalibration->GetT0Det();

  SetDetector(det);

}

//_____________________________________________________________________________
//AliTRDtransform::AliTRDtransform(const AliTRDtransform &t)
//  :AliTransform(t)
AliTRDtransform::AliTRDtransform(const AliTRDtransform &t)
  :TObject(t)
  ,fGeo(0x0)
  ,fDetector(t.fDetector)
  ,fParam(0x0)
  ,fCalibration(0x0)
  ,fCalVdriftROC(0x0)
  ,fCalT0ROC(0x0)
  ,fCalVdriftDet(0x0)
  ,fCalT0Det(0x0)
  ,fCalVdriftDetValue(0)
  ,fCalT0DetValue(0)
  ,fSamplingFrequency(0)
  ,fPadPlane(0x0)
  ,fZShiftIdeal(0)
  ,fMatrix(0x0)
{
  //
  // AliTRDtransform copy constructor
  //

  if (fGeo) {
    delete fGeo;
  }
  fGeo               = new AliTRDgeometry();
  fGeo->CreateClusterMatrixArray();

  fParam             = AliTRDCommonParam::Instance();
  if (!fParam) {
    AliError("Could not get common parameters\n");
  }
  fSamplingFrequency = fParam->GetSamplingFrequency();

  fCalibration = AliTRDcalibDB::Instance();
  if (!fCalibration) {
    AliError("Cannot find calibration object");
  }
  fCalVdriftDet      = fCalibration->GetVdriftDet();
  fCalT0Det          = fCalibration->GetT0Det();

}

//_____________________________________________________________________________
AliTRDtransform::~AliTRDtransform()
{
  //
  // AliTRDtransform destructor
  //

  if (fGeo) {
    delete fGeo;
  }

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

  // Get the detector wise defined calibration values
  fCalVdriftDetValue = fCalVdriftDet->GetValue(det);
  fCalT0DetValue     = fCalT0Det->GetValue(det);

  // Shift needed to define Z-position relative to middle of chamber
  Int_t layer        = fGeo->GetLayer(det);
  Int_t stack        = fGeo->GetStack(det);
  fPadPlane          = fGeo->GetPadPlane(layer,stack);
  fZShiftIdeal       = 0.5 * (fPadPlane->GetRow0() + fPadPlane->GetRowEnd());

  // Get the current transformation matrix
  fMatrix            = fGeo->GetClusterMatrix(det);

}

//_____________________________________________________________________________
Bool_t AliTRDtransform::Transform(Double_t *x, Int_t *i, UInt_t time, Bool_t &out, Int_t  /*coordinateType*/)
{
  //
  // Transforms the local cluster coordinates into calibrated 
  // space point positions defined in the local tracking system.
  //
  // Here the calibration for T0, Vdrift and ExB is applied as well.
  //
  // Input:
  //   x[0] = COL-position relative to the center pad (pad units)
  //   x[1] = cluster signal in left pad
  //   x[2] = cluster signal in middle pad
  //   x[3] = cluster signal in right pad
  //   i[0] = ROW pad number
  //   i[1] = COL pad number
  //   time = time bin number (uncalibrated for t0)
  //
  // Output:
  //   x[0] = X-positions in tracking CS
  //   x[1] = Y-positions in tracking CS
  //   x[2] = Z-positions in tracking CS
  //   x[3] = total cluster charge
  //   x[4] = error in Y-direction
  //   x[5] = error in Z-direction
  //   i[2] = time bin number (calibrated for t0)
  //

  Double_t posLocal[3];
  Double_t posTracking[3];

  Int_t row = i[0];
  Int_t col = i[1];

  // Parameters to adjust the X position
  const Double_t kX0shift = 2.52 + 0.04273; //[cm]
  const Double_t kT0shift = 3.19e-3;        //[us]

  if (!fMatrix) {

    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
    x[4] = 0.0;
    x[5] = 0.0;
    i[2] = 0;

    return kFALSE;

  }
  else {
 
    // Calibration values
    Double_t vdrift      = fCalVdriftDetValue * fCalVdriftROC->GetValue(col,row);
    Double_t t0          = fCalT0DetValue     + fCalT0ROC->GetValue(col,row);

    // T0 correction
    Double_t timeT0Cal   = time - t0;
    // Calculate the X-position,
    Double_t xLocal      = ((timeT0Cal + 0.5) / fSamplingFrequency + kT0shift) * vdrift; 

    // Length of the amplification region
    Double_t ampLength   = (Double_t) AliTRDgeometry::CamHght();
    // The drift distance
    Double_t driftLength = TMath::Max(xLocal - 0.5*ampLength,0.0);
    // ExB correction
    Double_t exbCorr     = AliTRDCommonParam::Instance()->GetOmegaTau(vdrift);

    // Pad dimensions
    Double_t rowSize     = fPadPlane->GetRowSize(row);
    Double_t colSize     = fPadPlane->GetColSize(col);

    // Invert the X-position,
    // apply ExB correction to the Y-position
    // and move to the Z-position relative to the middle of the chamber
    posLocal[0] = -xLocal;
    posLocal[1] =  (fPadPlane->GetColPos(col) + (0.5 + x[0]) * colSize) - driftLength * exbCorr;
    posLocal[2] =  (fPadPlane->GetRowPos(row) -         0.5  * rowSize) - fZShiftIdeal;

    // Go to tracking coordinates
    fMatrix->LocalToMaster(posLocal,posTracking);

    // The total charge of the cluster
    Double_t q0             = x[1];
    Double_t q1             = x[2];
    Double_t q2             = x[3];
    Double_t clusterCharge  = q0 + q1 + q2; 
    Double_t clusterSigmaY2 = 0.0;
    if (clusterCharge > 0.0) {
      clusterSigmaY2 = (q1 * (q0 + q2) + 4.0 * q0 * q2)
    	             / (clusterCharge*clusterCharge);
    }

    // Output values
    x[0] = posTracking[0] + kX0shift;
    x[1] = posTracking[1];
    x[2] = posTracking[2];
    x[3] = clusterCharge;
    x[4] = colSize*colSize * (clusterSigmaY2 + 1.0/12.0);
    x[5] = rowSize*rowSize / 12.0;                                       
    i[2] = TMath::Nint(timeT0Cal);
		
    // For TRD tracking calibration awareness
    out  = ((i[2] < 0) || (i[2] > Int_t(3.5 * fSamplingFrequency/vdrift))) ? kTRUE : kFALSE; 

    return kTRUE;

  }

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

  if (setDet) {
    SetDetector(c->GetDetector());
  }

  // Transform the local cluster coordinates into recalibrated 
  // space point positions defined in the local tracking system.
  // Here the calibration for T0, Vdrift and ExB is applied as well.
  Double_t clusterXYZ[6];
  clusterXYZ[0] = c->GetCenter();
  clusterXYZ[1] = 0.0;
  clusterXYZ[2] = 0.0;
  clusterXYZ[3] = 0.0;
  clusterXYZ[4] = 0.0;
  clusterXYZ[5] = 0.0;
  Int_t    clusterRCT[3];
  clusterRCT[0] = c->GetPadRow();
  clusterRCT[1] = c->GetPadCol();
  clusterRCT[2] = 0;
  Int_t time    = c->GetPadTime();
  Bool_t out;
  Transform(clusterXYZ,clusterRCT,((UInt_t) time), out, 0);

  // Set the recalibrated coordinates
  c->SetX(clusterXYZ[0]);
  c->SetY(clusterXYZ[1]);
  c->SetZ(clusterXYZ[2]);
  c->SetLocalTimeBin(((Char_t) clusterRCT[2]));
  c->SetInChamber(!out);

}
