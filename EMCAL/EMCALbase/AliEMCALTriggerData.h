#ifndef ALIEMCALTRIGGERDATA_H
#define ALIEMCALTRIGGERDATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliEMCALTriggerTypes.h"

#include <TObject.h>
#include <TVector2.h>
#include <TClonesArray.h>

class AliEMCALTriggerData : public TObject 
{

/// \class AliEMCALTriggerTRUDCSConfig
/// \brief  Trigger data container
/// \ingroup EMCALbase
/// \author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
///
/// EMCal trigger data container: can be used independently 
/// of the data stream (simulation or raw data)
/// for transient storage of trigger data
public:
  /// \enum ETriggerDataMode_t
  /// \brief Type of the data
  enum ETriggerDataMode_t {
    kSimulation,            ///< Simulation
    kRawData                ///< Raw Data
  };

  /// \enum L0DecisionOrigin
  enum L0DecisionOrigin_t {
    kOriginRawData,         ///< Decision from raw stream
    kOriginRecalc           ///< Decision from recalculation based on the L0 amplitudes
  };

  /// \enum L0TriggerStatus
  /// \brief L0 trigger decision (fired or not fired)
  enum L0TriggerStatus_t {
    kNotFired,              ///< L0 trigger has been fired
    kFired                  ///< L0 trigger has not been fired
  };

  /// \brief Default constructor
  AliEMCALTriggerData();

  /// \brief Destructor	
  virtual      ~AliEMCALTriggerData();
  
  /// \brief Set the mode of the data handling
  /// \param mode Mode of the data handling (Simulation / Raw Data)
  virtual void  SetMode(ETriggerDataMode_t mode) { fMode = mode ; }

  /// \brief Get the mode of the data handling
  /// \return Mode of hte data handling (Simulation / Reconstruction)
  virtual ETriggerDataMode_t GetMode()    const  { return fMode ; }
  
  /// \brief Set the trigger status of a given TRU
  /// \param datamode Reconstruction mode providing the L0 decision
  /// \param truID TRU number for which to set the trigger status
  /// \param triggerstatus Trigger status (fired or not fired)
  virtual void  SetL0Trigger(L0DecisionOrigin_t datamode, Int_t truID,  L0TriggerStatus_t triggestatus)       { fL0Trigger[datamode][truID]    = triggestatus ; }

  /// \brief Get L0 status for a given TRU
  /// \param datamode Mode of the L0 statsus (raw stream or recalculated)
  /// \param truID Index of the TRU
  /// \param triggerstatus[out] L0 trigger status (fired or not fired)
  virtual void  GetL0Trigger( L0DecisionOrigin_t datamode, Int_t truID, L0TriggerStatus_t& triggestatus ) const {    triggestatus = fL0Trigger[datamode][truID] ; }

  /// \brief Get L0 status for a given TRU
  /// \param datamode Mode of the L0 statsus (raw stream or recalculated)
  /// \param truID Index of the TRU
  /// \return L0 trigger status (fired or not fired)
  virtual L0TriggerStatus_t GetL0Trigger( L0DecisionOrigin_t datamode, Int_t truID           ) const { return fL0Trigger[datamode][truID] ; }
  
  /// \brief Set L1 Gamma threshold
  /// \param thresholdID Number ID of the threshold (0 - EG1, 1 - EG2)
  /// \param thresholdvalue ADC value of the threshold
  virtual void  SetL1GammaThreshold(int thresholdID, int  thresholdvalue) { fL1GammaThreshold[thresholdID]    = thresholdvalue ; }

  /// \brief Set L1 Jet threshold
  /// \param thresholdID Number ID of the threshold (0 - EJ1, 1 - EJ2)
  /// \param thresholdvalue ADC value of the threshold
  virtual void  SetL1JetThreshold(  int thresholdID, int  thresholdvalue) {   fL1JetThreshold[thresholdID]    = thresholdvalue ; }

  virtual void  SetL1FrameMask(     Int_t  v)      {     fL1FrameMask        = v ; }            
  
  /// \brief Get the L1 Gamma ADC threshold for a given threshold ID
  /// \param thresholdID Index of the threshold (0 - EG1, 1 - EG2)
  /// \return Threshold value in ADC counts for given threshold ID
  virtual Int_t GetL1GammaThreshold(  int thresholdID) const { return fL1GammaThreshold[thresholdID] ; }

  /// \brief Get the L1 Jet ADC threshold for a given threshold ID
  /// \param thresholdID Index of the threshold (0 - EJ1, 1 - EJ2)
  /// \return Threshold value in ADC counts for given threshold ID
  virtual Int_t GetL1JetThreshold(    int thresholdID) const { return   fL1JetThreshold[thresholdID] ; }

  virtual Int_t GetL1FrameMask(            ) const { return      fL1FrameMask    ; }            
  
  /// \brief set the L1 V0 params for V0-multiplicity dependent threshold adjustment
  /// \param v0params Params for V0-dependent threshold calculation
  virtual void  SetL1V0(            Int_t* v0params)      { for (int i = 0; i <  2; i++) fL1V0[i]          = v0params[i] ; }

  virtual void  SetL1TriggerType(   Int_t* v)      { for (int i = 0; i < 15; i++) fL1TriggerType[i] = v[i] ; }
  
  /// \brief Get parameters for the V0-multiplicity dependent thresholds
  /// \param v0params[out] V0-dependent threholds settings
  virtual void  GetL1V0(         Int_t  v0params[]) const { for (int i = 0; i <  2; i++) v0params[i] = fL1V0[i]          ; }

  virtual void  GetL1TriggerType(Int_t  v[]) const { for (int i = 0; i < 15; i++) v[i] = fL1TriggerType[i] ; }
  
  /// \brief Mark L1 data as decoded / not decoded
  /// \param decoded If true the raw data has been decoded
  virtual void  SetL1DataDecoded(  Bool_t  decoded)       {    fL1DataDecoded = decoded ; } 

  /// \brief Mark L1 raw data as present / not present
  /// \param hasrawdata If true the the raw data is present
  virtual void  SetL1RawData(      Bool_t  hasrawdata)    {        fL1RawData = hasrawdata ; }

  /// \brief Set L1 median value
  /// \param median Median value
  virtual void  SetMedian(         Int_t  median)       {           fMedian = median ; }
  
  /// \brief Check if the L1 raw data was decoded
  /// \return Decoding status
  virtual Bool_t IsL1DataDecoded(          ) const { return fL1DataDecoded ; }

  /// \brief Check if the L1 raw data was present
  /// \return True if the L1 raw data was present
  virtual Int_t     HasL1RawData(          ) const { return     fL1RawData ; }

  /// \brief Get L1 median value
  /// \return L1 median value
  virtual Int_t        GetMedian(          ) const { return        fMedian ; }

  /// \brief Check whether the data handling is in simulation mode
  /// \return True if the data handling is in simulation mode, false otherwise 
  bool IsSimulationMode() const { return fMode == ETriggerDataMode_t::kSimulation; }

  /// \brief Check whether the data handling is in raw data mode
  /// \brief True if running in simulation mode, false otherwise
  bool IsRawDataMode() const { return fMode == ETriggerDataMode_t::kRawData; }

  /// \brief Set data handling to simulation mode
  void SetSimulationMode() { fMode = kSimulation; }

  /// \brief Set data handling to reconstruction mode
  void SetRawDataMode() { fMode = kRawData; }
  
  /// \brief Dump parameters value
  virtual void  Scan() const;

  /// \brief Reset
  virtual void  Reset();
  
private:

  enum {
    kMaxTRU = 52,
    kMaxL0Modes = 2,
    kNL1Thresholds = 2
  };
  
  AliEMCALTriggerData(           const AliEMCALTriggerData& rhs); // NOT implemented
  AliEMCALTriggerData& operator=(const AliEMCALTriggerData& rhs); // NOT implemented
  
  ETriggerDataMode_t  fMode;                                  ///< Simulation/Raw
  
  L0TriggerStatus_t   fL0Trigger[kMaxL0Modes][kMaxTRU];       ///< Triggering TRU
  
  Int_t        fL1GammaThreshold[kNL1Thresholds];             ///< L1-g threshold
  Int_t          fL1JetThreshold[kNL1Thresholds];             ///< L1-j threshold
  
  Int_t                    fL1V0[2];                          ///< V0 charges
  Int_t             fL1FrameMask;                             ///< Frame mask
  Int_t           fL1TriggerType[19];                         ///< Trigger type
  
  Bool_t          fL1DataDecoded;                             ///< Raw data decoded
  Bool_t              fL1RawData;                             ///< Raw data present
  Int_t                  fMedian;                             ///< Median
  
  /// \cond CLASSIMP
  ClassDef(AliEMCALTriggerData,3) ;
  /// \endcond
  
};

#endif // ALIEMCALTRIGGERDATA_H
