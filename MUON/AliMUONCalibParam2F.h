/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUONCalibParam2F
/// \brief Implementation of AliMUONVCalibParam for pair of floats
/// 
//  Author Laurent Aphecetche

#ifndef ALIMUONCALIBPARAM2F_H
#define ALIMUONCALIBPARAM2F_H

#ifndef ALIMUONVCALIBPARAM_H
#  include "AliMUONVCalibParam.h"
#endif

class AliMUONCalibParam2F : public AliMUONVCalibParam
{
public:
  AliMUONCalibParam2F();
  AliMUONCalibParam2F(Int_t theSize, Float_t fillWithValue=0);
  AliMUONCalibParam2F(const AliMUONCalibParam2F& other);
  AliMUONCalibParam2F& operator=(const AliMUONCalibParam2F& other);
  
  virtual ~AliMUONCalibParam2F();

  /// Return dimension (here 2)
  virtual Int_t Dimension() const { return 2; } 
  
  virtual void Print(Option_t* opt="") const;
  
  virtual void SetValueAsFloat(Int_t i, Int_t j, Float_t value);
  virtual void SetValueAsInt(Int_t i, Int_t j, Int_t value);
  
  /// Return size - the number of float pair we hold
  virtual Int_t Size() const { return fSize; } 

  virtual Float_t ValueAsFloat(Int_t i, Int_t j=0) const;
  virtual Int_t ValueAsInt(Int_t i, Int_t j=0) const;
     
private:
  void CopyTo(AliMUONCalibParam2F& destination) const;
  Int_t Index(Int_t i, Int_t j) const;  
    
private:
  Int_t fSize; ///< The number of float pair we hold
  Int_t fN;    ///< The total number of floats we hold (2*fSize)

  /// The values array
  Float_t* fValues; //[fN] The values array
  
  ClassDef(AliMUONCalibParam2F,1) // Container for calibration parameters
};

#endif
