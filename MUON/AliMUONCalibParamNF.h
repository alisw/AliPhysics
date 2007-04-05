/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup calib
/// \class AliMUONCalibParamNF
/// \brief Implementation of AliMUONVCalibParam for tuples of floats
/// 
//  Author Laurent Aphecetche

#ifndef ALIMUONCALIBPARAMNF_H
#define ALIMUONCALIBPARAMNF_H

#ifndef ALIMUONVCALIBPARAM_H
#  include "AliMUONVCalibParam.h"
#endif

class AliMUONCalibParamNF : public AliMUONVCalibParam
{
public:
  AliMUONCalibParamNF();
  AliMUONCalibParamNF(Int_t dimension, Int_t theSize, Float_t fillWithValue=0);
  AliMUONCalibParamNF(const AliMUONCalibParamNF& other);
  AliMUONCalibParamNF& operator=(const AliMUONCalibParamNF& other);
  
  virtual ~AliMUONCalibParamNF();

  /// Return dimension
  virtual Int_t Dimension() const { return fDimension; }
  
  virtual void Print(Option_t* opt="") const;
  
  virtual void SetValueAsFloat(Int_t i, Int_t j, Float_t value);
  virtual void SetValueAsInt(Int_t i, Int_t j, Int_t value);
  
  /// Return size - the number of float pair we hold
  virtual Int_t Size() const { return fSize; } 

  virtual Float_t ValueAsFloat(Int_t i, Int_t j=0) const;
  virtual Int_t ValueAsInt(Int_t i, Int_t j=0) const;
     
private:
  void CopyTo(AliMUONCalibParamNF& destination) const;
  Int_t Index(Int_t i, Int_t j) const;  
    
private:
  Int_t fDimension; ///< dimension of this object
  Int_t fSize; ///< The number of float pair we hold
  Int_t fN;    ///< The total number of floats we hold (fDimension*fSize)

  /// The values array
  Float_t* fValues; //[fN] The values array
  
  ClassDef(AliMUONCalibParamNF,1) // Container for calibration parameters
};

#endif
