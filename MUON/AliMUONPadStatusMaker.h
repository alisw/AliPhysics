#ifndef ALIMUONPADSTATUSMAKER_H
#define ALIMUONPADSTATUSMAKER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup reco
/// \class AliMUONPadStatusMaker
/// \brief Make a 2DStore of pad statuses, using different sources of information
/// 
/// \author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif
#ifndef ROOT_TVector2
#  include "TVector2.h"
#endif

class AliMUONCalibrationData;
class AliMUONV2DStore;
class TMap;

class AliMUONPadStatusMaker : public TObject
{
public:
  AliMUONPadStatusMaker(const AliMUONCalibrationData& calibData);
  virtual ~AliMUONPadStatusMaker();
      
  AliMUONV2DStore* MakeGainStatus(const AliMUONV2DStore& gainValues) const;

  AliMUONV2DStore* MakeHVStatus(const TMap& hvMap) const;

  AliMUONV2DStore* MakePedestalStatus(const AliMUONV2DStore& pedValues) const;

  AliMUONV2DStore* MakeStatus() const;

  TVector2 HVSt12Limits() const { return fHVSt12Limits; }
  TVector2 HVSt345Limits() const { return fHVSt345Limits; }
  
  TVector2 PedMeanLimits() const { return fPedMeanLimits; }
  TVector2 PedSigmaLimits() const { return fPedSigmaLimits; }
  
  void SetHVSt12Limits(float low, float high) { fHVSt12Limits.Set(low,high); }
  void SetHVSt345Limits(float low, float high) { fHVSt345Limits.Set(low,high); }

  void SetPedMeanLimits(float low, float high) { fPedMeanLimits.Set(low,high); }
  void SetPedSigmaLimits(float low, float high) { fPedSigmaLimits.Set(low,high); }
  
private:
  AliMUONPadStatusMaker(const AliMUONPadStatusMaker&);
  AliMUONPadStatusMaker& operator=(const AliMUONPadStatusMaker&);
  
private:
  
  AliMUONV2DStore* Combine(const AliMUONV2DStore& store1,
                           const AliMUONV2DStore& store2,
                           Int_t binShift) const;
  
  Bool_t GetSt12Status(const TMap& hvMap, Int_t detElemId, Int_t sector,
                         Bool_t& hvChannelTooLow,
                         Bool_t& hvChannelTooHigh,
                         Bool_t& hvChannelON) const;
  
  
  Bool_t GetSt345Status(const TMap& hvMap, Int_t detElemId, Int_t pcbIndex,
                        Bool_t& hvChannelTooLow,
                        Bool_t& hvChannelTooHigh,
                        Bool_t& hvChannelON,
                        Bool_t& hvSwitchON) const;
  
  void SetStatusSt12(AliMUONV2DStore& hvStatus, 
                     Int_t detElemId, Int_t sector, Int_t status) const;

  void SetStatusSt345(AliMUONV2DStore& hvStatus,
                      Int_t detElemId, Int_t pcbIndex, Int_t status) const;

  
private:
  
  enum EPedestalStatus
  {
    kPedOK = 0,
    kPedMeanZero = (1<<1),
    kPedMeanTooLow = (1<<2),
    kPedMeanTooHigh = (1<<3),
    kPedSigmaTooLow = (1<<4),
    kPedSigmaTooHigh = (1<<5)
  };
  
  enum EHVError 
  {
    kHVOK = 0,
    kHVError = (1<<0),
    kHVTooLow = (1<<1),
    kHVTooHigh = (1<<2),
    kHVChannelOFF = (1<<3),
    kHVSwitchOFF = (1<<4)
    
  };
  
  const AliMUONCalibrationData& fCalibrationData; //!< helper class to get data access (not owner)
  TVector2 fPedMeanLimits; //!< Low and High threshold for pedestal mean
  TVector2 fPedSigmaLimits; //!< Low and High threshold for pedestal sigma
  
  TVector2 fHVSt12Limits; //!< Low and High threshold for St12 HV
  TVector2 fHVSt345Limits; //!< Low and High threshold for St345 HV
  
  ClassDef(AliMUONPadStatusMaker,0) // Creates pad statuses from ped,gain,hv
};

#endif
