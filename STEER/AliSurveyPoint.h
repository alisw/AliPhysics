#ifndef ALI_SURVEY_POINT_H
#define ALI_SURVEY_POINT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliSurveyPoint					   //
//  Retrieve and Convert survey data into ROOT Objects		   //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TString.h>

class AliSurveyPoint: public TObject {

 public:
  AliSurveyPoint();
  AliSurveyPoint(TString name, Float_t x, Float_t y, Float_t z,
		 Float_t precX, Float_t precY, Float_t precZ, 
                 Char_t type, Bool_t Target);
  virtual ~AliSurveyPoint();
  
  TString GetPointName() {return fPointName;};
  Float_t GetX() {return fX;};
  Float_t GetY() {return fY;};
  Float_t GetZ() {return fZ;};
  Float_t GetPrecisionX() {return fPrecisionX;};
  Float_t GetPrecisionY() {return fPrecisionY;};
  Float_t GetPrecisionZ() {return fPrecisionZ;};
  Char_t GetType() {return fType;};
  Bool_t GetTarget() {return fTargetUsed;};
  
  virtual const char* GetName() const {return fPointName.Data();}; 

  void PrintPoint();
  
 private:
  TString fPointName;
  Float_t fX;
  Float_t fY;
  Float_t fZ;
  Float_t fPrecisionX;
  Float_t fPrecisionY;
  Float_t fPrecisionZ;
  Char_t fType;
  Bool_t fTargetUsed;
  
  ClassDef(AliSurveyPoint, 1);
};

#endif
