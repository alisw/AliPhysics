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

#include "AliLog.h"
#include "AliZDCRecoParam.h"

class AliZDCRecoParampp : public AliZDCRecoParam {
 public:
  AliZDCRecoParampp();
  virtual ~AliZDCRecoParampp() {;}

  // make reco parameters for p-p collisions
  static AliZDCRecoParampp *GetLowFluxParam();
  
  virtual TH1D* GethNpartDist() const {AliError("NO centrality determination in p-p!"); return 0;}
  virtual TH1D* GethbDist() 	const {AliError("NO centrality determination in p-p!"); return 0;}
  virtual Float_t GetClkCenter() const {AliError("NO centrality determination in p-p!");return 0;}
  
 protected:
  
  AliZDCRecoParampp(const AliZDCRecoParampp&);
  AliZDCRecoParampp& operator =(const AliZDCRecoParampp&);
  
 ClassDef(AliZDCRecoParampp, 1)

};

#endif
