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

#ifndef ROOT_TTask
#include "TTask.h"
#endif

class AliMUONCalibrationData;
class AliMUONData;
class AliMUONLogger;
class AliMUONV2DStore;

class AliMUONDigitCalibrator : public TTask
{
public:
  AliMUONDigitCalibrator(AliMUONData* data, 
                         AliMUONCalibrationData* calib,
                         Bool_t createAndUseStatusMap=kTRUE);
  virtual ~AliMUONDigitCalibrator();
  
  virtual void Exec(Option_t*);

private:    
    /// Not implemented
    AliMUONDigitCalibrator(const AliMUONDigitCalibrator& other);
    /// Not implemented
    AliMUONDigitCalibrator& operator=(const AliMUONDigitCalibrator& other);

    AliMUONData* fData;                       //!< MUON data 
    AliMUONCalibrationData* fCalibrationData; //!< Calibration data
    AliMUONV2DStore* fStatusMap; //!< Channel status map
    AliMUONLogger* fLogger; //!< to log repeated messages
    
  ClassDef(AliMUONDigitCalibrator,2) // Calibrate raw digit
};

#endif
