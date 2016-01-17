#ifndef ALITPCTRANSFORM_H
#define ALITPCTRANSFORM_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/// \class AliTPCTransform
/// \brief Class for tranformation of the coordinate frame
///
/// Transformation
///  local coordinate frame (sector, padrow, pad, timebine) ==>
///  rotated global (tracking) cooridnate frame (sector, lx,ly,lz)

class AliTPCRecoParam;
#include "AliTPCChebCorr.h"
#include "AliTransform.h"
#include <time.h>

class AliTPCTransform:public AliTransform {
public:
  AliTPCTransform();
  AliTPCTransform(const AliTPCTransform& transform);

  virtual ~AliTPCTransform();
  virtual void Transform(Double_t *x,Int_t *i,UInt_t time,
			 Int_t coordinateType);
  void SetPrimVertex(Double_t *vtx);
  void Local2RotatedGlobal(Int_t sec,  Double_t *x) const;
  void RotatedGlobal2Global(Int_t sector,Double_t *x) const;
  void Global2RotatedGlobal(Int_t sector,Double_t *x) const;
  void GetCosAndSin(Int_t sector,Double_t &cos,Double_t &sin) const;
  time_t GetCurrentTimeStamp() const { return fCurrentTimeStamp;}
  const AliTPCRecoParam * GetCurrentRecoParam() const {return fCurrentRecoParam;}
  AliTPCRecoParam * GetCurrentRecoParamNonConst() const {return fCurrentRecoParam;}
  UInt_t GetCurrentRunNumber() const { return fCurrentRun;}
  const AliTPCChebCorr* GetCorrMapCache0() const {return fCorrMapCache0;}
  const AliTPCChebCorr* GetCorrMapCache1() const {return fCorrMapCache1;}
  //
  TObjArray* LoadCorrectionMaps() const;
  // set current values
  //
  void SetCurrentRecoParam(AliTPCRecoParam* param){fCurrentRecoParam=param;}
  void SetCurrentRun(Int_t run){fCurrentRun=run;}
  void SetCurrentTimeStamp(time_t timeStamp);
  void ApplyTransformations(Double_t *xyz, Int_t volID);
  //
  // new correction maps
  Bool_t  UpdateTimeDependentCache();
  void    ApplyCorrectionMap(int roc, int row, double xyzSect[3]);
  void    ApplyDistortionMap(int roc, double xyzLab[3]);
  void    EvalCorrectionMap(int roc, int row, const double xyz[3], float res[3]);
  Float_t EvalCorrectionMap(int roc, int row, const double xyz[3], int dimOut);
  void    EvalDistortionMap(int roc, const double xyzSector[3], float res[3]);
  const   Float_t* GetLastMapCorrection() const {return fLastCorr;}
  //
  static void RotateToSectorUp(float *x, int& idROC);
  static void RotateToSectorDown(float *x, int& idROC);
  static void RotateToSectorUp(double *x, int& idROC);
  static void RotateToSectorDown(double *x,  int& idROC);
  static int  SectorUp(int idROC);
  static int  SectorDown(int idROC);
  static double GetMaxY2X() {return fgkMaxY2X;}
  //
private:
  AliTPCTransform& operator=(const AliTPCTransform&); // not implemented
  Float_t  fLastCorr[3]; ///!<! last correction from the map
  Double_t fCoss[18];  ///< cache the transformation
  Double_t fSins[18];  ///< cache the transformation
  Double_t fPrimVtx[3];///< position of the primary vertex - needed for TOF correction
  AliTPCRecoParam * fCurrentRecoParam; //!<! current reconstruction parameters
  const AliTPCChebCorr* fCorrMapCache0;      //!<! current correction map0 (for 1st time bin if time-dependent)
  const AliTPCChebCorr* fCorrMapCache1;      //!<! current correction map1 (for 2nd time bin if time-dependent)
  Int_t    fCurrentRun;                //!<! current run
  time_t   fCurrentTimeStamp;          //!<! current time stamp
  Bool_t   fTimeDependentUpdated;      //!<! flag successful update of time dependent stuff
  /// \cond CLASSIMP
  static const Double_t fgkSin20;       // sin(20)
  static const Double_t fgkCos20;       // sin(20)
  static const Double_t fgkMaxY2X;      // tg(10)
  //
  ClassDef(AliTPCTransform,3)
  /// \endcond
};

//______________________________________________________
inline void AliTPCTransform::ApplyCorrectionMap(int roc, int row, double xyzSect[3])
{
  // apply correction from the map to a point at given ROC and row (IROC/OROC convention)
  EvalCorrectionMap(roc, row, xyzSect, fLastCorr);
  for (int i=3;i--;) xyzSect[i] += fLastCorr[i];
  //
}

//_________________________________________________
inline void AliTPCTransform::RotateToSectorUp(float *x, int& idROC)
{
  // rotate point in sector coordinates to sector+1
  idROC += ((idROC%18)==17) ? -17 : 1;
  float tmp = x[0];
  x[0] = fgkCos20*tmp - fgkSin20*x[1];
  x[1] = fgkSin20*tmp + fgkCos20*x[1];
}

//_________________________________________________
inline void AliTPCTransform::RotateToSectorDown(float *x, int& idROC)
{
  // rotate point in sector coordinates to sector-1
  idROC += ((idROC%18)== 0) ?  17 : -1; // change to the lower sector
  float tmp = x[0];
  x[0] = fgkCos20*tmp + fgkSin20*x[1];
  x[1] =-fgkSin20*tmp + fgkCos20*x[1];
}

//_________________________________________________
inline void AliTPCTransform::RotateToSectorUp(double *x, int& idROC)
{
  // rotate point in sector coordinates to sector+1
  idROC += ((idROC%18)==17) ? -17 : 1;
  double tmp = x[0];
  x[0] = fgkCos20*tmp - fgkSin20*x[1];
  x[1] = fgkSin20*tmp + fgkCos20*x[1];
}

//_________________________________________________
inline void AliTPCTransform::RotateToSectorDown(double *x, int& idROC)
{
  // rotate point in sector coordinates to sector-1
  idROC += ((idROC%18)== 0) ?  17 : -1; // change to the lower sector
  double tmp = x[0];
  x[0] = fgkCos20*tmp + fgkSin20*x[1];
  x[1] =-fgkSin20*tmp + fgkCos20*x[1];
}

//_________________________________________________
inline int AliTPCTransform::SectorUp(int idROC)
{
  // sector+1
  return idROC + (((idROC%18)==17) ? -17 : 1);
}

//_________________________________________________
inline int AliTPCTransform::SectorDown(int idROC)
{
  // sector-1
  return idROC + (((idROC%18)== 0) ?  17 : -1); // change to the lower sector
}


#endif
