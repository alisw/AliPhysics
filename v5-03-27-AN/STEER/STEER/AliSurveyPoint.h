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
  
  TString GetPointName() const {return fPointName;};
  Float_t GetX() const {return fX;};
  Float_t GetY() const {return fY;};
  Float_t GetZ() const {return fZ;};
  Float_t GetPrecisionX() const {return fPrecisionX;};
  Float_t GetPrecisionY() const {return fPrecisionY;};
  Float_t GetPrecisionZ() const {return fPrecisionZ;};
  Char_t GetType() const {return fType;};
  Bool_t GetTarget() const {return fTargetUsed;};

  // Overwritten to enable the use of TObjArray::FindObject()
  virtual const char* GetName() const {return fPointName.Data();}; 

  void PrintPoint();
  
 private:
  TString fPointName; // Point name used by the survey group (TS/SU)
  Float_t fX; // X Coordinate
  Float_t fY; // Y Coordinate
  Float_t fZ; // Z Coordinate
  Float_t fPrecisionX; // Precision in the X axis
  Float_t fPrecisionY; // Precision in the Y axis
  Float_t fPrecisionZ; // Precision in the Z axis
  Char_t fType; // Type of point: M(easured), T(ransformed)
  Bool_t fTargetUsed; // Was a survey target used to measure this point?
                      // If yes, the user must take that into account
  
  ClassDef(AliSurveyPoint, 1);
};

#endif
