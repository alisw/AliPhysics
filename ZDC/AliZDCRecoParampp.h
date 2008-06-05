#ifndef ALIZDCRECOPARAMPP_H
#define ALIZDCRECOPARAMPP_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////
//                                                      //
//   Class with ZDC reconstruction parameters           //
//   		proton - ptoron collisions              //
//   Origin: Chiara.Oppedisano@to.infn.it               //
//                                                      //
//////////////////////////////////////////////////////////

#include <TF1.h>
#include "AliLog.h"
#include "AliZDCRecoParam.h"

//class TF1;

class AliZDCRecoParampp : public AliZDCRecoParam {
 public:
  AliZDCRecoParampp();
  virtual ~AliZDCRecoParampp();

  // make reco parameters for p-p collisions
  static AliZDCRecoParampp *GetppRecoParam();
  
  void PrintParameters() const; 
  
  virtual TF1* GetfZNCen()  const {AliError("NO function can be defined in pp!"); return 0;}    
  virtual TF1* GetfZNPer()  const {AliError("NO function can be defined in pp!"); return 0;}    
  virtual TF1* GetfZPCen()  const {AliError("NO function can be defined in pp!"); return 0;}    
  virtual TF1* GetfZPPer()  const {AliError("NO function can be defined in pp!"); return 0;}    
  virtual TF1* GetfZDCCen() const {AliError("NO function can be defined in pp!"); return 0;}   
  virtual TF1* GetfZDCPer() const {AliError("NO function can be defined in pp!"); return 0;}   
  virtual TF1* GetfbCen()   const {AliError("NO function can be defined in pp!"); return 0;}     
  virtual TF1* GetfbPer()   const {AliError("NO function can be defined in pp!"); return 0;}     
  virtual TF1* GetfZEMn()   const {AliError("NO function can be defined in pp!"); return 0;}     
  virtual TF1* GetfZEMp()   const {AliError("NO function can be defined in pp!"); return 0;}     
  virtual TF1* GetfZEMsp()  const {AliError("NO function can be defined in pp!"); return 0;}    
  virtual TF1* GetfZEMb()   const {AliError("NO function can be defined in pp!"); return 0;}   
  //
  virtual Float_t GetZEMEndValue()    const {AliError("NO function can be defined in pp!"); return 0;}     
  virtual Float_t GetZEMCutFraction() const {AliError("NO function can be defined in pp!"); return 0;}  
  virtual Float_t GetDZEMSup()        const {AliError("NO function can be defined in pp!"); return 0;}  	     
  virtual Float_t GetDZEMInf()        const {AliError("NO function can be defined in pp!"); return 0;}  	     
  virtual Float_t GetEZN1MaxValue()   const {AliError("NO function can be defined in pp!"); return 0;}    
  virtual Float_t GetEZP1MaxValue()   const {AliError("NO function can be defined in pp!"); return 0;}    
  virtual Float_t GetEZDC1MaxValue()  const {AliError("NO function can be defined in pp!"); return 0;}   
  virtual Float_t GetEZN2MaxValue()   const {AliError("NO function can be defined in pp!"); return 0;}    
  virtual Float_t GetEZP2MaxValue()   const {AliError("NO function can be defined in pp!"); return 0;}    
  virtual Float_t GetEZDC2MaxValue()  const {AliError("NO function can be defined in pp!"); return 0;}   
  
 protected:
  
  AliZDCRecoParampp(const AliZDCRecoParampp&);
  AliZDCRecoParampp& operator =(const AliZDCRecoParampp&);
  
 ClassDef(AliZDCRecoParampp, 1)

};

#endif
