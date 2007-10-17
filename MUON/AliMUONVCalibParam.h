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
  //AliMUONVCalibParam(Int_t id0);
  AliMUONVCalibParam(Int_t id0, Int_t id1);
  virtual ~AliMUONVCalibParam();

  virtual const char* GetName() const;
  
  /// First id of this object
  virtual Int_t ID0() const;
  
  /// Second id of this object (might not be required)
  virtual Int_t ID1() const;
  
  /// whether or not the value we store are packed, e.g. as v = a*cste + b
  virtual Bool_t IsPacked() const { return kFALSE; }
  
  /// j indices in following methods are valid from 0 to Dimension()-1.
  virtual Int_t Dimension() const = 0;

  /** Set one value, for channel i, dimension j. Consider value is a double.
    Only ok to use if IsDoublePrecision() is kTRUE.
    */
  virtual void SetValueAsDouble(Int_t i, Int_t j, Double_t value);

  /// Same as above but w/o bound checking
  virtual void SetValueAsDoubleFast(Int_t i, Int_t j, Double_t value);

  /// Set one value, for channel i, dimension j. Consider value is a float.
  virtual void SetValueAsFloat(Int_t i, Int_t j, Float_t value) = 0;

  /** Set one value, for channel i, dimension j. Consider value is a float.
    Assume (i,j) are valid indices, i.e. do not check them.
    */
  virtual void SetValueAsFloatFast(Int_t i, Int_t j, Float_t value) = 0;

  /// Set one value, for channel i, dimension j. Consider value is an integer.
  virtual void SetValueAsInt(Int_t i, Int_t j, Int_t value) = 0;

  /// Same as above but w/o bound checkings.
  virtual void SetValueAsIntFast(Int_t i, Int_t j, Int_t value) = 0;

  /// The number of channels handled by this object.
  virtual Int_t Size() const = 0;

  /// Whether we can store double precision values   
  virtual Bool_t IsDoublePrecision() const { return kFALSE; }
  
  /** Retrieve the value for a given (channel,dim) as a double.
      Only ok if IsDoublePrecision() is kTRUE.
      (i,j) are checked to within boundaries
    */
  virtual Double_t ValueAsDouble(Int_t i, Int_t j=0) const;

  /** Retrieve the value for a given (channel,dim) as a double.
    Only ok if IsDoublePrecision() is kTRUE.
    Fast means there's no bound checking on (i,j)
    */
  virtual Double_t ValueAsDoubleFast(Int_t i, Int_t j=0) const;
  
  /** Retrieve the value for a given (channel,dim) as a float, with
      bound checking on (i,j).
    */
  virtual Float_t ValueAsFloat(Int_t i, Int_t j=0) const = 0;

  /// Same as above but without bound checking.
  virtual Float_t ValueAsFloatFast(Int_t i, Int_t j=0) const = 0;

  /** Retrieve the value for a given (channel,dim) as an integer.
      With bound checking.
    */
  virtual Int_t ValueAsInt(Int_t i, Int_t j=0) const = 0;

  /// Same as above but w/o bound checking.
  virtual Int_t ValueAsIntFast(Int_t i, Int_t j=0) const = 0;

  /// Unpack a value into a couple (a,b). Returns false if IsPacked()==kFALSE
  virtual Bool_t UnpackValue(Int_t /*value*/, Int_t& /*a*/, Int_t& /*b*/) const { return kFALSE; }
  
  /// Pack (a,b) as a single int. Returns false if IsPacked()==kFALSE
  virtual Bool_t PackValues(Int_t /*a*/, Int_t /*b*/, Int_t& /*packedValue*/) const { return kFALSE; }
  
  /// Return 1E38 as invalid float value
  static Float_t InvalidFloatValue() { return 1E38; }

protected:
    
  static UInt_t BuildUniqueID(Int_t id0, Int_t id1);
  static void DecodeUniqueID(UInt_t uniqueID, Int_t& id0, Int_t& id1);
  static Int_t ID0(UInt_t uniqueID);
  static Int_t ID1(UInt_t uniqueID);
  
  ClassDef(AliMUONVCalibParam,0) // Base class for a calibration data holder (usually for 64 channels)
};

#endif
