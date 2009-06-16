#ifndef ALIZDCRECOPARAM_H
#define ALIZDCRECOPARAM_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////
//                                                      //
//   Class with ZDC reconstruction parameters           //
//   Origin: Chiara.Oppedisano@to.infn.it               //
//                                                      //
//////////////////////////////////////////////////////////

#include <TH2F.h>
#include <TH1D.h>
#include <TF1.h>
#include "AliDetectorRecoParam.h"

//class TF1;

class AliZDCRecoParam : public AliDetectorRecoParam {
 public:
  AliZDCRecoParam();
  virtual ~AliZDCRecoParam();

  virtual TH2F*   GethZDCvsZEM()  const = 0;	
  virtual TH2F*   GethZDCCvsZEM() const = 0;	
  virtual TH2F*   GethZDCAvsZEM() const = 0;	
  virtual TH1D*   GethNpartDist() const = 0;	
  virtual TH1D*   GethbDist()     const = 0;	
  virtual Float_t GetClkCenter()  const = 0;
  
  virtual void SetZDCvsZEM(TH2F* /*hCorr*/)  {printf(" AliZDCRecoParam::SetZDCvsZEM doesn't set anything!!!");}    
  virtual void SetZDCCvsZEM(TH2F* /*hCorr*/) {printf(" AliZDCRecoParam::SetZDCCvsZEM doesn't set anything!!!");}    
  virtual void SetZDCAvsZEM(TH2F* /*hCorr*/) {printf(" AliZDCRecoParam::SetZDCAvsZEM doesn't set anything!!!");}    

    
  virtual void PrintParameters() const {;} 
  
 protected:
  
  AliZDCRecoParam(const AliZDCRecoParam&);
  AliZDCRecoParam& operator =(const AliZDCRecoParam&);
   
 ClassDef(AliZDCRecoParam, 2)

};

#endif
