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

#include <TF1.h>
#include "AliDetectorRecoParam.h"

//class TF1;

class AliZDCRecoParam : public AliDetectorRecoParam {
 public:
  AliZDCRecoParam();
  virtual ~AliZDCRecoParam();

  virtual TF1* GetfZNCen()  const =0;    
  virtual TF1* GetfZNPer()  const =0;    
  virtual TF1* GetfZPCen()  const =0;    
  virtual TF1* GetfZPPer()  const =0;    
  virtual TF1* GetfZDCCen() const =0;   
  virtual TF1* GetfZDCPer() const =0;   
  virtual TF1* GetfbCen()   const =0;     
  virtual TF1* GetfbPer()   const =0;     
  virtual TF1* GetfZEMn()   const =0;     
  virtual TF1* GetfZEMp()   const =0;     
  virtual TF1* GetfZEMsp()  const =0;    
  virtual TF1* GetfZEMb()   const =0;   
  //
  virtual Float_t GetZEMEndValue()    const =0;     
  virtual Float_t GetZEMCutFraction() const =0;  
  virtual Float_t GetDZEMSup()        const =0;  	     
  virtual Float_t GetDZEMInf()        const =0;  	     
  virtual Float_t GetEZN1MaxValue()   const =0;    
  virtual Float_t GetEZP1MaxValue()   const =0;    
  virtual Float_t GetEZDC1MaxValue()  const =0;   
  virtual Float_t GetEZN2MaxValue()   const =0;    
  virtual Float_t GetEZP2MaxValue()   const =0;    
  virtual Float_t GetEZDC2MaxValue()  const =0;   
  
    
  virtual void PrintParameters() const {} 
  
 protected:
  
  AliZDCRecoParam(const AliZDCRecoParam&);
  AliZDCRecoParam& operator =(const AliZDCRecoParam&);
   
 ClassDef(AliZDCRecoParam, 2)

};

#endif
