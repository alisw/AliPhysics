#ifndef TPCParamCR_H
#define TPCParamCR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for TPC parameters          //
////////////////////////////////////////////////
#include "AliTPCParam.h"

class AliTPCRF1D;
class AliTPCPRF2D;

class AliTPCParamCR : public AliTPCParam {
public:
  AliTPCParamCR();
  virtual ~AliTPCParamCR();
  Int_t  CalcResponse(Float_t* x, Int_t * index, Int_t dummy=0);
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
  Bool_t Update();            //recalculate and check geometric parameters 
  void SetDefault();          //set default parameters
  void   SetInnerPRF(AliTPCPRF2D * prf) {fInnerPRF = prf;}
  void   SetOuterPRF(AliTPCPRF2D * prf) {fOuterPRF = prf;}
  void   SetTimeRF(AliTPCRF1D * timerf) {fTimeRF = timerf;}

  AliTPCPRF2D * GetInnerPRF() const {return fInnerPRF;}
  AliTPCPRF2D * GetOuterPRF() const {return fOuterPRF;}
  AliTPCRF1D  * GetTimeRF()   const {return fTimeRF;}
protected:
  AliTPCPRF2D * fInnerPRF;         //!pad response function object for inner sector
  AliTPCPRF2D * fOuterPRF;         //!pad response function object for inner sector  
  AliTPCRF1D  * fTimeRF;           //!time response function object
  Float_t       fFacSigma;         //factor-how many sigma of response I accept
  ClassDef(AliTPCParamCR,1)  //parameter  object for set:TPC
};

#endif  
