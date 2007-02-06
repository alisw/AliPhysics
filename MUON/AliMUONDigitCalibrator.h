/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONDigitCalibrator
/// \brief Class to calibrate the digits
/// 
/// \author Laurent Aphecetche

#ifndef ALIMUONDIGITCALIBRATOR_H
#define ALIMUONDIGITCALIBRATOR_H

#ifndef ROOT_TTask
#include "TTask.h"
#endif

class AliMUONCalibrationData;
class AliMUONData;
class AliMUONV2DStore;

class AliMUONDigitCalibrator : public TTask
{
public:
  AliMUONDigitCalibrator(AliMUONData* data, AliMUONCalibrationData* calib);
  virtual ~AliMUONDigitCalibrator();
  
  virtual void Exec(Option_t*);

private:    
    AliMUONDigitCalibrator(const AliMUONDigitCalibrator& other);
    AliMUONDigitCalibrator& operator=(const AliMUONDigitCalibrator& other);

    AliMUONData* fData;                       //!< MUON data 
    AliMUONCalibrationData* fCalibrationData; //!< Calibration data
    AliMUONV2DStore* fStatusMap; //!< Channel status map
    
  ClassDef(AliMUONDigitCalibrator,2) // Calibrate raw digit
};

#endif
