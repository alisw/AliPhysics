#ifndef TPCParamSR_H
#define TPCParamSR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for TPC parameters          //
////////////////////////////////////////////////
#include "AliTPCParam.h"

class AliTPCRF1D;
class AliTPCPRF2D;

class AliTPCParamSR : public AliTPCParam {
public:
  AliTPCParamSR();
  virtual ~AliTPCParamSR();
  Int_t  CalcResponse(Float_t* x, Int_t * index, Int_t row);
  //calculate bin response as function of the input position -x 
  //return number of valid response bin
  
  void XYZtoCRXYZ(Float_t *xyz, 
		  Int_t &sector, Int_t &padrow, Int_t option=3) const;
  //transform global position to the position relative to the sector padrow
  //if option=0  X calculate absolute            calculate sector
  //if option=1  X           absolute            use input sector
  //if option=2  X           relative to pad row calculate sector
  //if option=3  X           relative            use input sector

  void CRXYZtoXYZ(Float_t *xyz,
	    const Int_t &sector, const Int_t & padrow, Int_t option=3) const;  
  //transform relative position  to the gloabal position
  void TransformTo8(Float_t *xyz, Int_t *index) const;
  void TransformTo2(Float_t *xyz, Int_t *index) const;
  Bool_t Update();            //recalculate and check geometric parameters 
  void SetDefault();          //set default parameters
  void   SetInnerPRF(AliTPCPRF2D * prf) {fInnerPRF = prf;}
  void   SetOuterPRF(AliTPCPRF2D * prf) {fOuterPRF = prf;}
  void   SetTimeRF(AliTPCRF1D * timerf) {fTimeRF = timerf;}

  AliTPCPRF2D * GetInnerPRF() const {return fInnerPRF;}
  AliTPCPRF2D * GetOuterPRF() const {return fOuterPRF;}
  AliTPCRF1D  * GetTimeRF()   const {return fTimeRF;}
  void SetFacSigmaPadRow(Float_t fac=3.) {fFacSigmaPadRow=fac;}
  void SetFacSigmaPad(Float_t fac=3.) {fFacSigmaPad=fac;}
  void SetFacSigmaTime(Float_t fac=3.) {fFacSigmaTime=fac;}

  virtual Float_t GetPrimaryLoss(Float_t *x, Int_t *index, Float_t *angle);
  virtual Float_t GetTotalLoss(Float_t *x, Int_t *index, Float_t *angle);

  virtual void GetClusterSize(Float_t *x, Int_t *index, Float_t *angle, Int_t mode, Float_t *sigma);
  virtual void GetSpaceResolution(Float_t *x, Int_t *index, Float_t *angle, Float_t amplitude, Int_t mode,Float_t *sigma);
  virtual Float_t  GetAmp(Float_t *x, Int_t *index, Float_t *angle);
  virtual Float_t * GetAnglesAccMomentum(Float_t *x, Int_t * index, Float_t* momentum, Float_t *angle); 
 
protected:
  AliTPCPRF2D * fInnerPRF;         //pad response function object for inner sector
  AliTPCPRF2D * fOuterPRF;         //pad response function object for inner sector  
  AliTPCRF1D  * fTimeRF;           //time response function object
  Float_t      fFacSigmaPadRow;    //factor-how many sigma of response I accept
  Float_t      fFacSigmaPad;
  Float_t      fFacSigmaTime; 
  ClassDef(AliTPCParamSR,1)  //parameter  object for set:TPC
};

#endif  






