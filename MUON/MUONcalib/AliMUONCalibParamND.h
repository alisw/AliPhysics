/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup calib
/// \class AliMUONCalibParamND
/// \brief Implementation of AliMUONVCalibParam for tuples of double
/// 
//  Author Laurent Aphecetche

#ifndef ALIMUONCALIBPARAMND_H
#define ALIMUONCALIBPARAMND_H

#ifndef ALIMUONVCALIBPARAM_H
#  include "AliMUONVCalibParam.h"
#endif

class AliMUONCalibParamND : public AliMUONVCalibParam
{
public:
  AliMUONCalibParamND();
  AliMUONCalibParamND(Int_t dimension, Int_t theSize, Int_t id0, Int_t id1, Double_t fillWithValue=0);
  AliMUONCalibParamND(const AliMUONCalibParamND& other);
  AliMUONCalibParamND& operator=(const AliMUONCalibParamND& other);
  
  virtual ~AliMUONCalibParamND();

  /// Own clone methods (as the default TObject::Clone turned out to be pretty slow !)
  virtual TObject* Clone(const char* = "") const { return new AliMUONCalibParamND(*this); }
  
  /// Return dimension
  virtual Int_t Dimension() const { return fDimension; }
  
  virtual void Print(Option_t* opt="") const;
  
  virtual void SetValueAsDouble(Int_t i, Int_t j, Double_t value);
  virtual void SetValueAsDoubleFast(Int_t i, Int_t j, Double_t value);
  virtual void SetValueAsFloat(Int_t i, Int_t j, Float_t value);
  virtual void SetValueAsFloatFast(Int_t i, Int_t j, Float_t value);
  virtual void SetValueAsInt(Int_t i, Int_t j, Int_t value);
  virtual void SetValueAsIntFast(Int_t i, Int_t j, Int_t value);
  
  /// Return size - the number of double tuples we hold
  virtual Int_t Size() const { return fSize; } 

  virtual Float_t ValueAsFloat(Int_t i, Int_t j=0) const;
  virtual Float_t ValueAsFloatFast(Int_t i, Int_t j=0) const;
  virtual Double_t ValueAsDouble(Int_t i, Int_t j=0) const;
  virtual Double_t ValueAsDoubleFast(Int_t i, Int_t j=0) const;
  virtual Int_t ValueAsInt(Int_t i, Int_t j=0) const;
  virtual Int_t ValueAsIntFast(Int_t i, Int_t j=0) const;
     
  /// Advertise that we can store double precision values
  virtual Bool_t IsDoublePrecision() const { return kTRUE; }
  
private:
  void CopyTo(AliMUONCalibParamND& destination) const;
  Int_t Index(Int_t i, Int_t j) const;  
  Int_t IndexFast(Int_t i, Int_t j) const;  
    
private:
  Int_t fDimension; ///< dimension of this object
  Int_t fSize; ///< The number of double tuples we hold
  Int_t fN;    ///< The total number of floats we hold (fDimension*fSize)

  /// The values array
  Double_t* fValues; //[fN] The values array
  
  ClassDef(AliMUONCalibParamND,1) // Container for calibration parameters
};

#endif
