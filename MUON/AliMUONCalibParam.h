/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUONCalibParam
/// \brief Container for 2 floats (one mean and one sigma)
/// 
/// \author Laurent Aphecetche

#ifndef AliMUONCALIBPARAM_H
#define AliMUONCALIBPARAM_H

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONCalibParam : public TObject
{
public:
  AliMUONCalibParam(Float_t mean=0.0, Float_t sigma=0.0);
  virtual ~AliMUONCalibParam();
  
  Float_t Mean() const;
  
  virtual void Print(Option_t*) const;
  
  void Set(Float_t mean, Float_t sigma);
  
  Float_t Sigma() const;
  
private:
  Float_t fMean;
  Float_t fSigma;
  
  ClassDef(AliMUONCalibParam,1) // 
};

#endif
