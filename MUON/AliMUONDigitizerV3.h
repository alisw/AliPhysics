/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup sim
/// \class AliMUONDigitizerV3
/// \brief Digitizer (from SDigit to Digit), performing digit de-calibration.
///
/// \author Laurent Aphecetche

#ifndef ALIMUONDIGITIZERV3_H
#define ALIMUONDIGITIZERV3_H

#ifndef ALIDIGITIZER_H
#include "AliDigitizer.h"
#endif

class AliMUONCalibrationData;
class AliMUONData;
class AliMUONDigit;
class TClonesArray;
class TString;

class AliMUONDigitizerV3 : public AliDigitizer
{
public:
  enum ETriggerCodeVersion
  {
    kTriggerDecision=-1,
    kTriggerElectronics
  };
  
  AliMUONDigitizerV3(AliRunDigitizer* manager=0, 
                     ETriggerCodeVersion=kTriggerDecision);
  virtual ~AliMUONDigitizerV3();

  virtual void Exec(Option_t* opt="");
  
  virtual Bool_t Init();

private:
    
  void AddOrUpdateDigit(TClonesArray& array, 
                        const AliMUONDigit& digit);
    
  void ApplyResponse();

  void ApplyResponseToDigit(AliMUONDigit& digit);

  Int_t FindDigitIndex(TClonesArray& array, const AliMUONDigit& digit);

  AliMUONData* GetDataAccess(const TString& folderName);

  Bool_t MergeDigits(const AliMUONDigit& src, AliMUONDigit& srcAndDest);

  void MergeWithSDigits(AliMUONData& outputData, const AliMUONData& inputData, 
                        Int_t mask);
  
private:
  Bool_t fIsInitialized; // are we initialized ?
  AliMUONData* fOutputData; //! pointer to access digits
  AliMUONCalibrationData* fCalibrationData; //! pointer to access calib parameters
  TTask* fTriggerProcessor; // pointer to the trigger part of the job
  ETriggerCodeVersion fTriggerCodeVersion; // which version of trigger job
  
  ClassDef(AliMUONDigitizerV3,2) // MUON Digitizer V3-2
};

#endif
