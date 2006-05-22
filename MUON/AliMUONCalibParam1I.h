/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUONCalibParam1I
/// \brief Implementation of AliMUONVCalibParam.
/// 
/// Handle the case of 1 integer parameter per channel.
///
/// \author Laurent Aphecetche

#ifndef ALIMUONCALIBPARAM1I_H
#define ALIMUONCALIBPARAM1I_H

#ifndef ALIMUONVCALIBPARAM_H
#  include "AliMUONVCalibParam.h"
#endif

class AliMUONCalibParam1I : public AliMUONVCalibParam
{
public:
  AliMUONCalibParam1I();
  AliMUONCalibParam1I(Int_t theSize, Int_t fillWithValue=0);
  AliMUONCalibParam1I(const AliMUONCalibParam1I& other);
  AliMUONCalibParam1I& operator=(const AliMUONCalibParam1I& other);
  virtual ~AliMUONCalibParam1I();

  /// The dimension of this object is 1, i.e. it's a vector.
  virtual Int_t Dimension() const { return 1; }
  
  virtual void Print(Option_t* opt="") const;
  
  /** Set the value located at i (j must be zero).
    * Note that the internal structure store ints, so the float is rounded
    * prior to be stored
    */
  virtual void SetValueAsFloat(Int_t i, Int_t j, Float_t value);
  
  /// Set the value located at i (j must be zero), as an int.
  virtual void SetValueAsInt(Int_t i, Int_t j, Int_t value);
  
  /// The number of values we hold.
  virtual Int_t Size() const { return fSize; }

  /// Retrieve the value located at i (j must remain zero), as a float.
  virtual Float_t ValueAsFloat(Int_t i, Int_t j=0) const;
  
  /// Retrieve the value located at i (j must remain zero).
  virtual Int_t ValueAsInt(Int_t i, Int_t j=0) const;

private:
  /// Copy *this to destination (used by copy ctor and operator=).
  void CopyTo(AliMUONCalibParam1I& destination) const;
  
private:
  Int_t fSize;    ///< The number of values we hold

  ///  The values array 
  Int_t* fValues; //[fSize] The values array
  
  ClassDef(AliMUONCalibParam1I,1) // Container for calibration parameters
};

#endif
