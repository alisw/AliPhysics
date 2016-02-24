#ifndef ALIMUONPADSTATUSMAKER_H
#define ALIMUONPADSTATUSMAKER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONPadStatusMaker
/// \brief Make a 2DStore of pad statuses, using different sources of information
/// 
//  Author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif
#ifndef ROOT_TVector2
#  include "TVector2.h"
#endif

class TExMap;
class AliMUONCalibrationData;
class AliMUONRecoParam;
class AliMUONVCalibParam;
class AliMUONVTrackerData;
class AliMUONVStore;

using std::ostream;

class AliMUONPadStatusMaker : public TObject
{
public:
  AliMUONPadStatusMaker(const AliMUONCalibrationData& calibData);
  virtual ~AliMUONPadStatusMaker();
  
  /// Get the reference to the calibrationdata object we use.
  const AliMUONCalibrationData& CalibrationData() const { return fkCalibrationData; }
  
  /** Get access to internal status store (for debug only, as it may not be complete,
    depending on whether you've already called PadStatus for all possible de,manu
    combinations or not...
    */
  AliMUONVStore* StatusStore() const { return fStatus; }

  AliMUONVCalibParam* PadStatus(Int_t detElemId, Int_t manuId) const;

  Int_t PadStatus(Int_t detElemId, Int_t manuId, Int_t manuChannel) const;
  
  AliMUONVStore* NeighboursStore() const;
  
  AliMUONVCalibParam* Neighbours(Int_t detElemId, Int_t manuId) const;
  
  static TString AsString(Int_t status);

  static TString AsCondition(Int_t status);
  
  /// Return HV threshold
  Double_t HVLimit(Int_t chamberId) const;
  
  /// Return Low and High threshold for pedestal mean
  TVector2 PedMeanLimits() const { return fPedMeanLimits; }
  /// Return Low and High threshold for pedestal sigma
  TVector2 PedSigmaLimits() const { return fPedSigmaLimits; }
  
  /// Set HV limit
  void SetHVLimit(Int_t chamberId, Double_t hv);

  /// Set Low and High threshold for pedestal mean
  void SetPedMeanLimits(float low, float high) { fPedMeanLimits.Set(low,high); }
  /// Set Low and High threshold for pedestal sigma 
  void SetPedSigmaLimits(float low, float high) { fPedSigmaLimits.Set(low,high); }
  	
  /// Set Low and High manu occupancy limits
  void SetManuOccupancyLimits(float low, float high) { fManuOccupancyLimits.Set(low,high); }
  /// Get manu occupancy limits
  TVector2 ManuOccupancyLimits() const { return fManuOccupancyLimits; }

  /// Set Low and High bus patch occupancy limits
  void SetBuspatchOccupancyLimits(float low, float high) { fBuspatchOccupancyLimits.Set(low,high); }
  /// Get bus patch occupancy limits
  TVector2 BuspatchOccupancyLimits() const { return fBuspatchOccupancyLimits; }

  /// Set Low and High DE occupancy limits
  void SetDEOccupancyLimits(float low, float high) { fDEOccupancyLimits.Set(low,high); }
  /// Get DE occupancy limits
  TVector2 DEOccupancyLimits() const { return fDEOccupancyLimits; }
  
  void SetLimits(const AliMUONRecoParam& recoParams);

  void Report(UInt_t mask);
  
  static Float_t SwitchValue(const TObjArray& dcsArray);
  
  Int_t LVStatus(Int_t detElemId) const;
  
  Int_t HVStatus(Int_t detElemId, Int_t manuId) const;
  
  Int_t OccupancyStatus(Int_t detElemId, Int_t manuId) const;

  static void DecodeStatus(Int_t status, Int_t& pedStatus, Int_t& hvStatus, 
                           Int_t&  lvStatus, Int_t& otherStatus);
  static Int_t BuildStatus(Int_t pedStatus, Int_t hvStatus, 
                           Int_t lvStatus, Int_t otherStatus);
private:
  /// Not implemented
  AliMUONPadStatusMaker(const AliMUONPadStatusMaker&);
  /// Not implemented
  AliMUONPadStatusMaker& operator=(const AliMUONPadStatusMaker&);

  
  AliMUONVCalibParam* ComputeStatus(Int_t detElemId, Int_t manuId) const;

  Bool_t HVSt12Status(Int_t detElemId, Int_t sector,
                      Bool_t& hvChannelTooLow,
                      Bool_t& hvChannelTooHigh,
                      Bool_t& hvChannelON) const;
  
  
  Bool_t HVSt345Status(Int_t detElemId, Int_t pcbIndex,
                       Bool_t& hvChannelTooLow,
                       Bool_t& hvChannelTooHigh,
                       Bool_t& hvChannelON,
                       Bool_t& hvSwitchON) const;

  void SetHVStatus(Int_t detElemId, Int_t index, Int_t status) const;

  Int_t CheckConfigConsistencyWithPedestalInformation(Int_t detElemId,Int_t manuId) const;

private:
  /// General status
  enum EGeneralStatus
  {
    kMissing = (1<<7)
  };

  /// Pedestal status
  enum EPedestalStatus
  {
    kPedOK = 0,
    kPedMeanZero = (1<<1),
    kPedMeanTooLow = (1<<2),
    kPedMeanTooHigh = (1<<3),
    kPedSigmaTooLow = (1<<4),
    kPedSigmaTooHigh = (1<<5),
    
    kPedMissing = kMissing // please always use last bit for meaning "missing"
  };
  
  /// LV status
  enum ELVStatus
  {
    kLVOK = 0,
    kLVTooLow = (1<<3),
    
    kLVMissing = kMissing // please always use last bit for meaning "missing"
  };
  
  /// HV Error
  enum EHVError 
  {
    kHVOK = 0,
    kHVError = (1<<0),
    kHVTooLow = (1<<1),
    kHVTooHigh = (1<<2), // no longer to be used
    kHVChannelOFF = (1<<3),
    kHVSwitchOFF = (1<<4),

    kHVMissing = kMissing // please always use last bit for meaning "missing"
  };
  
  /// Other
  enum EOccupancyStatus
  {
    kManuOccupancyTooLow = (1<<1),
    kManuOccupancyTooHigh = (1<<2),
    kBusPatchOccupancyTooLow = (1<<3),
    kBusPatchOccupancyTooHigh = (1<<4),
    kDEOccupancyTooLow = (1<<5),
    kDEOccupancyTooHigh = (1<<6),
    kBusPatchRemovedByPAR = (1<<7),
    
  };
  
  const AliMUONCalibrationData& fkCalibrationData; //!<! helper class to get data access (not owner)
  
  Double_t fHVLimit[10]; //!<! Low thresholds for HV

  TVector2 fPedMeanLimits; //!<! Low and High threshold for pedestal mean
  TVector2 fPedSigmaLimits; //!<! Low and High threshold for pedestal sigma
  
  TVector2 fManuOccupancyLimits; //!<! Low and High manu occupancy limits
  TVector2 fBuspatchOccupancyLimits; //!<! Low and High buspatch occupancy limits
  TVector2 fDEOccupancyLimits; //!<! Low and High DE occupancy limits
  
  AliMUONVStore* fStatus; //!<! statuses of the pads
  
  mutable TExMap* fHV; //!<! cache of hv statuses

  AliMUONVStore* fPedestals; //!<! pedestal values
  AliMUONVStore* fConfig; //!<! readout configuration
  
  AliMUONVTrackerData* fTrackerData; //!<! to get occupancies...
  
  ClassDef(AliMUONPadStatusMaker,0) // Creates pad statuses from ped,lv,hv
};

#endif
