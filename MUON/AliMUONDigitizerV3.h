/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup sim
/// \class AliMUONDigitizerV3
/// \brief Digitizer (from SDigit to Digit), performing digit de-calibration.
///
//  Author Laurent Aphecetche

#ifndef ALIMUONDIGITIZERV3_H
#define ALIMUONDIGITIZERV3_H

#ifndef ALIDIGITIZER_H
#include "AliDigitizer.h"
#endif

#include "TArrayI.h"

class AliMUONCalibrationData;
class AliMUONVDigit;
class AliMUONLogger;
class TClonesArray;
class TF1;
class TString;
class AliMUONVDigitStore;
class AliLoader;
class AliMUONVTriggerStore;
class AliMUONTriggerElectronics;
class AliMUONVCalibParam;
class AliMUONRecoParam;
class AliMUONTriggerChamberEfficiency;
class AliMUONTriggerUtilities;

class AliMUONDigitizerV3 : public AliDigitizer
{
public:
  AliMUONDigitizerV3(AliRunDigitizer* manager=0, 
                     Int_t generateNoisyDigits=1);
  
  virtual ~AliMUONDigitizerV3();

  virtual void Exec(Option_t* opt="");
  
  virtual Bool_t Init();

  static Int_t DecalibrateTrackerDigit(const AliMUONVCalibParam& pedestals,
                                       const AliMUONVCalibParam& gains,
                                       Int_t channel,
                                       Float_t charge,
                                       Bool_t addNoise=kFALSE,
                                       Bool_t noiseOnly=kFALSE,
                                       const TString& calibrationMode="NOGAIN");
  
  /// Set calibration (and recoparam) data
  void SetCalibrationData(AliMUONCalibrationData* calibrationData, AliMUONRecoParam* recoParam);

  /// Set the number of sigmas for pedestal cut
  static void SetNSigmas(Double_t nsigmas=4.0) { fgNSigmas = nsigmas; }

private:
  /// Not implemented
  AliMUONDigitizerV3(const AliMUONDigitizerV3& other);
  /// Not implemented
  AliMUONDigitizerV3& operator=(const AliMUONDigitizerV3& other);
    
  void ApplyResponse(const AliMUONVDigitStore& store, AliMUONVDigitStore& filteredStore);

  void ApplyResponseToTrackerDigit(AliMUONVDigit& digit, Bool_t addNoise);
  void ApplyResponseToTriggerDigit(AliMUONVDigit& digit);

  AliLoader* GetLoader(const TString& foldername);
  
private:  

  void GenerateNoisyDigits(AliMUONVDigitStore& digitStore);
  void GenerateNoisyDigitsForOneCathode(AliMUONVDigitStore& digitStore, 
                                        Int_t detElemId, Int_t cathode);
  void GenerateNoisyDigitsForTrigger(AliMUONVDigitStore& digitStore);

  void MergeWithSDigits(AliMUONVDigitStore*& digitStore,
                        const AliMUONVDigitStore& input,
                        Int_t mask);
  
  static TF1* NoiseFunction();
  
  void CreateInputDigitStores();
  
  void BuildTriggerStatusMap();
  Int_t GetArrayIndex(Int_t cathode, Int_t trigCh, Int_t localCircuit);

private:
  Bool_t fIsInitialized; ///< are we initialized ?
  AliMUONCalibrationData* fCalibrationData; //!< pointer to access calib parameters
  AliMUONTriggerElectronics* fTriggerProcessor; ///< pointer to the trigger part of the job
  TF1* fNoiseFunctionTrig; //!< function to get noise disribution on trig. chambers
  Int_t fGenerateNoisyDigits; //!< whether or not we should generate noise-only digits for tracker (1) and trigger (2)
  static Double_t fgNSigmas; ///< \brief number of sigmas above ped to use 
  /// for noise-only digit generation and zero-suppression
  AliMUONLogger* fLogger; //!< to keep track of messages
  AliMUONVTriggerStore* fTriggerStore; //!< trigger objects
  AliMUONVDigitStore* fDigitStore; //!< temporary digits
  AliMUONVDigitStore* fOutputDigitStore; //!< digits we'll output to disk
  TObjArray* fInputDigitStores; //!< input digit stores (one per input file
  AliMUONRecoParam* fRecoParam; //!< reco params (to know how to decalibrate) (not owner)
  AliMUONTriggerChamberEfficiency* fTriggerEfficiency; //!< trigger efficiency map
  AliMUONTriggerUtilities* fTriggerUtilities; //!< Trigger utilities for masks
  TArrayI fEfficiencyResponse; //!< Local board efficiency response
  
  ClassDef(AliMUONDigitizerV3,11) // MUON Digitizer V3-11
};

#endif
