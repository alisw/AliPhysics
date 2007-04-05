/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup calib
/// \class AliMUONCalibParamNI
/// \brief Implementation of AliMUONVCalibParam for tuples of ints
/// 
//  Author Laurent Aphecetche

#ifndef ALIMUONCALIBPARAMNI_H
#define ALIMUONCALIBPARAMNI_H

#ifndef ALIMUONVCALIBPARAM_H
#  include "AliMUONVCalibParam.h"
#endif

class AliMUONCalibParamNI : public AliMUONVCalibParam
{
public:
  AliMUONCalibParamNI();
  AliMUONCalibParamNI(Int_t dimension, Int_t theSize, Int_t fillWithValue=0, Int_t packingValue=0);
  AliMUONCalibParamNI(const AliMUONCalibParamNI& other);
  AliMUONCalibParamNI& operator=(const AliMUONCalibParamNI& other);
  
  virtual ~AliMUONCalibParamNI();

  /// Return dimension
  virtual Int_t Dimension() const { return fDimension; }
  
  virtual void Print(Option_t* opt="") const;
  
  virtual void SetValueAsFloat(Int_t i, Int_t j, Float_t value);
  virtual void SetValueAsInt(Int_t i, Int_t j, Int_t value);
  
  /// Return size - the number of float pair we hold
  virtual Int_t Size() const { return fSize; } 

  virtual Float_t ValueAsFloat(Int_t i, Int_t j=0) const;
  virtual Int_t ValueAsInt(Int_t i, Int_t j=0) const;
     
  virtual Bool_t UnpackValue(Int_t value, Int_t& a, Int_t& b) const;
  
  virtual Bool_t PackValues(Int_t a, Int_t b, Int_t& packedValue) const;
  
  virtual Bool_t IsPacked() const;

private:
  void CopyTo(AliMUONCalibParamNI& destination) const;
  Int_t Index(Int_t i, Int_t j) const;  
    
private:
  Int_t fDimension; ///< dimension of this object
  Int_t fSize; ///< The number of float pair we hold
  Int_t fN;    ///< The total number of floats we hold (fDimension*fSize)
  Int_t fPackingFactor; ///< packing factor, i.e. value = a*fPackingFactor + b
  
  /// The values array
  Int_t* fValues; //[fN] The values array
  
  ClassDef(AliMUONCalibParamNI,1) // Container for calibration parameters
};

#endif
