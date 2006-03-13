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
  virtual ~AliMUONCalibParam1I();

  virtual Int_t Dimension() const { return 1; }
  
  virtual void Print(Option_t* opt="") const;
  
  virtual void SetValueAsFloat(Int_t i, Int_t j, Float_t value);
  virtual void SetValueAsInt(Int_t i, Int_t j, Int_t value);
  
  virtual Int_t Size() const { return fSize; }

  virtual Float_t ValueAsFloat(Int_t i, Int_t j=0) const;
  virtual Int_t ValueAsInt(Int_t i, Int_t j=0) const;

 protected:
  AliMUONCalibParam1I(const AliMUONCalibParam1I& right);
  AliMUONCalibParam1I&  operator = (const AliMUONCalibParam1I& right);
 
 private:
  Int_t fSize;
  Int_t* fValues; //[fSize]
  
  ClassDef(AliMUONCalibParam1I,1) // Container for calibration parameters
};

#endif
