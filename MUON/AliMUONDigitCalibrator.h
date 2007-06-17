/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONDigitCalibrator
/// \brief Class to calibrate the digits
/// 
//  Author Laurent Aphecetche

#ifndef ALIMUONDIGITCALIBRATOR_H
#define ALIMUONDIGITCALIBRATOR_H

#ifndef ROOT_TObject
#include "TObject.h"
#endif

class AliMUONCalibrationData;
class AliMUONLogger;
class AliMUONVStore;
class AliMUONVDigitStore;
class AliMUONVDigit;

class AliMUONDigitCalibrator : public TObject
{
public:
  AliMUONDigitCalibrator(const AliMUONCalibrationData& calib,
                         Bool_t createAndUseStatusMap=kTRUE);
  
  virtual ~AliMUONDigitCalibrator();
  
  virtual void Calibrate(AliMUONVDigitStore& digitStore);
    
private:    
    /// Not implemented
    AliMUONDigitCalibrator(const AliMUONDigitCalibrator& other);
    /// Not implemented
    AliMUONDigitCalibrator& operator=(const AliMUONDigitCalibrator& other);

    virtual void CalibrateDigit(AliMUONVDigit& digit);

private:
    const AliMUONCalibrationData& fCalibrationData; //!< Calibration data
    AliMUONVStore* fStatusMap; //!< Channel status map
    AliMUONLogger* fLogger; //!< to log repeated messages
    
  ClassDef(AliMUONDigitCalibrator,3) // Calibrate raw digit
};

#endif
