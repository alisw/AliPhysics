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
class AliMUONPadStatusMaker;
class AliMUONPadStatusMapMaker;

class AliMUONDigitCalibrator : public TObject
{
public:
  AliMUONDigitCalibrator(const AliMUONCalibrationData& calib);
  
  virtual ~AliMUONDigitCalibrator();
  
  virtual void Calibrate(AliMUONVDigitStore& digitStore);
    
private:    
    /// Not implemented
    AliMUONDigitCalibrator(const AliMUONDigitCalibrator& other);
    /// Not implemented
    AliMUONDigitCalibrator& operator=(const AliMUONDigitCalibrator& other);

    virtual void CalibrateDigit(AliMUONVDigit& digit);

private:
    AliMUONLogger* fLogger; //!< to log repeated messages
    AliMUONPadStatusMaker* fStatusMaker; //!< to build pad statuses
    AliMUONPadStatusMapMaker* fStatusMapMaker; //!< to build status map
    AliMUONVStore* fPedestals; //!< pedestal values
    AliMUONVStore* fGains; //!< gain values
    
  ClassDef(AliMUONDigitCalibrator,4) // Calibrate raw digit
};

#endif
