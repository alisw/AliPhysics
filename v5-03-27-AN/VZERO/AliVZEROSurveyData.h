#ifndef ALIVZEROSurveyDATA_H
#define ALIVZEROSurveyDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//                                            // 
//  class for VZERO survey points management  //
//                                            //
////////////////////////////////////////////////

#include "TNamed.h"
#include "AliVZERO.h"

class AliVZEROSurveyData: public TNamed {

 public:
  AliVZEROSurveyData();
  AliVZEROSurveyData(const char* name);
  AliVZEROSurveyData(const AliVZEROSurveyData &calibda);
  AliVZEROSurveyData& operator= (const AliVZEROSurveyData &calibda);
  virtual ~AliVZEROSurveyData();
  void Reset();

  Float_t  GetPointA(Int_t ngA)   const {return fngA[ngA];}
  Float_t* GetPointA()   const {return (float*)fngA;}
  Float_t  GetPointB(Int_t ngB)   const {return fngB[ngB];}
  Float_t* GetPointB()   const {return (float*)fngB;}
  Float_t  GetPointC(Int_t ngC)   const {return fngC[ngC];}
  Float_t* GetPointC()   const {return (float*)fngC;}
  Float_t  GetPointD(Int_t ngD)   const {return fngD[ngD];}
  Float_t* GetPointD()   const {return (float*)fngD;}
  
  void     SetPointA(Float_t val, Int_t channel) {fngA[channel]=val;}
  void     SetPointA(Float_t* ngA);
  void     SetPointB(Float_t val, Int_t channel) {fngB[channel]=val;}
  void     SetPointB(Float_t* ngB);
  void     SetPointC(Float_t val, Int_t channel) {fngC[channel]=val;}
  void     SetPointC(Float_t* ngC);
  void     SetPointD(Float_t val, Int_t channel) {fngD[channel]=val;}
  void     SetPointD(Float_t* ngD);
  
 protected:
  Float_t  fngA[3];     //  Fiducial point A
  Float_t  fngB[3];     //  Fiducial point B
  Float_t  fngC[3];     //  Fiducial point C
  Float_t  fngD[3];     //  Fiducial point D
  
  ClassDef(AliVZEROSurveyData,1)    // VZERO Survey data
};

#endif
