/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUONVCalibParam
/// \brief Container of calibration values for a given number of channels.
/// 
//  Author Laurent Aphecetche

#ifndef ALIMUONVCALIBPARAM_H
#define ALIMUONVCALIBPARAM_H

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONVCalibParam : public TObject
{
public:
  AliMUONVCalibParam();
  virtual ~AliMUONVCalibParam();

  /// j indices in following methods are valid from 0 to Dimension()-1.
  virtual Int_t Dimension() const = 0;
  
  /// Set one value, for channel i, dimension j. Consider value is a float.
  virtual void SetValueAsFloat(Int_t i, Int_t j, Float_t value) = 0;
  
  /// Set one value, for channel i, dimension j. Consider value is an integer.
  virtual void SetValueAsInt(Int_t i, Int_t j, Int_t value) = 0;
  
  /// The number of channels handled by this object.
  virtual Int_t Size() const = 0;

  /// Retrieve the value for a given (channel,dim) as a float.
  virtual Float_t ValueAsFloat(Int_t i, Int_t j=0) const = 0;
  
    /// Retrieve the value for a given (channel,dim) as an integer.
  virtual Int_t ValueAsInt(Int_t i, Int_t j=0) const = 0;

  static Float_t InvalidFloatValue() { return 1E38; }
  
  ClassDef(AliMUONVCalibParam,0) // 
};

#endif
