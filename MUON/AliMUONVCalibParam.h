/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup calib
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

  /// whether or not the value we store are packed, e.g. as v = a*cste + b
  virtual Bool_t IsPacked() const { return kFALSE; }
  
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

  /// Unpack a value into a couple (a,b). Returns false if IsPacked()==kFALSE
  virtual Bool_t UnpackValue(Int_t /*value*/, Int_t& /*a*/, Int_t& /*b*/) const { return kFALSE; }
  
  /// Pack (a,b) as a single int. Returns false if IsPacked()==kFALSE
  virtual Bool_t PackValues(Int_t /*a*/, Int_t /*b*/, Int_t& /*packedValue*/) const { return kFALSE; }
  
  /// Return 1E38 as invalid float value
  static Float_t InvalidFloatValue() { return 1E38; }
  
  ClassDef(AliMUONVCalibParam,0) // 
};

#endif
