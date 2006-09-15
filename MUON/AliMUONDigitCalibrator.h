/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONDigitCalibrator
/// \brief Class to calibrate the digits
/// 
/// \author Laurent Aphecetche

#ifndef AliMUONDigitCalibrator_H
#define AliMUONDigitCalibrator_H

#ifndef ROOT_TTask
#include "TTask.h"
#endif

class AliMUONCalibrationData;
class AliMUONData;

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
  
  ClassDef(AliMUONDigitCalibrator,1) // Subtract pedestal from digit charge.
};

#endif
