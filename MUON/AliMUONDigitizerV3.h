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

#ifndef ROOT_TStopwatch
#  include "TStopwatch.h"
#endif

class AliMUONCalibrationData;
class AliMUONData;
class AliMUONDigit;
class AliMUONTriggerEfficiencyCells;
class TClonesArray;
class TF1;
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
                     ETriggerCodeVersion=kTriggerDecision,
                     Bool_t generateNoisyDigits=kTRUE);
  virtual ~AliMUONDigitizerV3();

  virtual void Exec(Option_t* opt="");
  
  virtual Bool_t Init();

private:
    
  AliMUONDigitizerV3(const AliMUONDigitizerV3& other);
  AliMUONDigitizerV3& operator=(const AliMUONDigitizerV3& other);
  
  void AddOrUpdateDigit(TClonesArray& array, 
                        const AliMUONDigit& digit);
    
  void ApplyResponse();

  void ApplyResponseToTrackerDigit(AliMUONDigit& digit, Bool_t addNoise);
  void ApplyResponseToTriggerDigit(AliMUONDigit& digit, AliMUONData* data);

private:  
  AliMUONDigit* FindCorrespondingDigit(AliMUONDigit& digit,AliMUONData* data) const;
  
  Int_t FindDigitIndex(TClonesArray& array, const AliMUONDigit& digit) const;

  void GenerateNoisyDigits();
  void GenerateNoisyDigitsForOneCathode(Int_t detElemId, Int_t cathode);

  AliMUONData* GetDataAccess(const TString& folderName);

  Bool_t MergeDigits(const AliMUONDigit& src, AliMUONDigit& srcAndDest);

  void MergeWithSDigits(AliMUONData& outputData, const AliMUONData& inputData, 
                        Int_t mask);
  
private:
  Bool_t fIsInitialized; ///< are we initialized ?
  AliMUONData* fOutputData; //!< pointer to access digits
  AliMUONCalibrationData* fCalibrationData; //!< pointer to access calib parameters
  TTask* fTriggerProcessor; ///< pointer to the trigger part of the job
  ETriggerCodeVersion fTriggerCodeVersion; ///< which version of trigger job
  AliMUONTriggerEfficiencyCells* fTriggerEfficiency; ///< trigger efficiency map  
  mutable TStopwatch fFindDigitIndexTimer; //!< counting time spent in FindDigitIndex
  TStopwatch fGenerateNoisyDigitsTimer; //!< counting time spent in GenerateNoisyDigits()
  TStopwatch fExecTimer; //!< couting time spent in Exec()  
  TF1* fNoiseFunction; //!< function to randomly get signal above n*sigma_ped
  Bool_t fGenerateNoisyDigits; //!< whether or not we should generate noise-only digits for tracker
  static const Double_t fgkNSigmas; ///< \brief number of sigmas above ped to use 
  /// for noise-only digit generation and zero-suppression

  ClassDef(AliMUONDigitizerV3,3) // MUON Digitizer V3-3
};

#endif
