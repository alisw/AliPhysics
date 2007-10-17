/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// $Id$

//-----------------------------------------------------------------------------
/// \class AliMUONPadStatusMaker
///
/// Make a 2DStore of pad statuses, using different sources of information,
/// like pedestal values, gain values, and HV values.
///
/// \author Laurent Aphecetche
//-----------------------------------------------------------------------------

#include "AliMUONPadStatusMaker.h"

#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCodeTimer.h"
#include "AliDCSValue.h"
#include "AliLog.h"
#include "AliMUON2DMap.h"
#include "AliMUON2DStoreValidator.h"
#include "AliMUONCalibParamNI.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONStringIntMap.h"
#include "AliMUONVCalibParam.h"
#include "AliMpArea.h"
#include "AliMpArrayI.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"
#include "AliMpDetElement.h"
#include "AliMpExMap.h"
#include "AliMpHVNamer.h"
#include "AliMpIntPair.h"
#include "AliMpManuUID.h"
#include "AliMpMotifMap.h"
#include "AliMpMotifPosition.h"
#include "AliMpPCB.h"
#include "AliMpPad.h"
#include "AliMpSegmentation.h"
#include "AliMpSlat.h"
#include "AliMpSlatSegmentation.h"
#include "AliMpStationType.h"
#include "AliMpVPadIterator.h"
#include "AliMpVSegmentation.h"
#include <Riostream.h>
#include <TArrayI.h>
#include <TExMap.h>
#include <TMap.h>
#include <TString.h>

/// \cond CLASSIMP
ClassImp(AliMUONPadStatusMaker)
/// \endcond

//_____________________________________________________________________________
AliMUONPadStatusMaker::AliMUONPadStatusMaker(const AliMUONCalibrationData& calibData)
: fCalibrationData(calibData),
  fGainA0Limits(0,1E30),
  fGainA1Limits(-1E-30,1E30),
  fGainThresLimits(0,4095),
  fHVSt12Limits(0,5000),
  fHVSt345Limits(0,5000),
  fPedMeanLimits(0,4095),
  fPedSigmaLimits(0,4095),
  fStatus(new AliMUON2DMap(true)),
  fHV(new TExMap),
  fPedestals(calibData.Pedestals()),
  fGains(calibData.Gains())
{
   /// ctor
    AliInfo(Form("ped store %s gain store %s",
                 fPedestals->ClassName(),
                 fGains->ClassName()));
}

//_____________________________________________________________________________
AliMUONPadStatusMaker::~AliMUONPadStatusMaker()
{
  /// dtor.
  delete fStatus;
  delete fHV;
}

//_____________________________________________________________________________
TString
AliMUONPadStatusMaker::AsString(Int_t status)
{
  /// return a human readable version of the integer status
  Int_t pedStatus;
  Int_t gainStatus;
  Int_t hvStatus;
  
  DecodeStatus(status,pedStatus,hvStatus,gainStatus);
  
//  /// Gain status
//  enum EGainStatus
//  {
//    kGainOK = 0,
//    kGainA0TooLow = (1<<1),
//    kGainA0TooHigh = (1<<2),
//    kGainA1TooLow = (1<<3),
//    kGainA1TooHigh = (1<<4),
//    kGainThresTooLow = (1<<5),
//    kGainThresTooHigh = (1<<6),
//    
//    kGainMissing = kMissing // please always use last bit for meaning "missing"
//  };
//  
//  /// Pedestal status
//  enum EPedestalStatus
//  {
//    kPedOK = 0,
//    kPedMeanZero = (1<<1),
//    kPedMeanTooLow = (1<<2),
//    kPedMeanTooHigh = (1<<3),
//    kPedSigmaTooLow = (1<<4),
//    kPedSigmaTooHigh = (1<<5),
//    
//    kPedMissing = kMissing // please always use last bit for meaning "missing"
//  };
//  
  TString s("PED ");
  
  if ( pedStatus == 0 ) s+= " OK";
  if ( pedStatus & kPedMeanZero ) s += " Mean is Zero. ";
  if ( pedStatus & kPedMeanTooLow ) s += " Mean Too Low. ";
  if ( pedStatus & kPedMeanTooHigh ) s += " Mean Too High. ";
  if ( pedStatus & kPedSigmaTooLow ) s += " Sigma Too Low. ";
  if ( pedStatus & kPedSigmaTooHigh ) s += " Sigma Too High. ";
  if ( pedStatus & kPedMissing ) s += " is missing.";
  
//  /// HV Error
//  enum EHVError 
//  {
//    kHVOK = 0,
//    kHVError = (1<<0),
//    kHVTooLow = (1<<1),
//    kHVTooHigh = (1<<2),
//    kHVChannelOFF = (1<<3),
//    kHVSwitchOFF = (1<<4),
//    
//    kHVMissing = kMissing // please always use last bit for meaning "missing"
//  };

  return s;
}

//_____________________________________________________________________________
Int_t
AliMUONPadStatusMaker::BuildStatus(Int_t pedStatus, 
                                   Int_t hvStatus, 
                                   Int_t gainStatus)
{
  /// Build a complete status from specific parts (ped,hv,gain)
  
  return ( hvStatus & 0xFF ) | ( ( pedStatus & 0xFF ) << 8 ) | 
  ( ( gainStatus & 0xFF ) << 16 );
}

//_____________________________________________________________________________
void
AliMUONPadStatusMaker::DecodeStatus(Int_t status, 
                                    Int_t& pedStatus, 
                                    Int_t& hvStatus, 
                                    Int_t& gainStatus)
{
  /// Decode complete status into specific parts (ped,hv,gain)
  
  gainStatus = ( status & 0xFF0000 ) >> 16;
  pedStatus = ( status & 0xFF00 ) >> 8;
  hvStatus = (status & 0xFF);
}

//_____________________________________________________________________________
Bool_t 
AliMUONPadStatusMaker::HVSt12Status(Int_t detElemId, Int_t sector,
                                    Bool_t& hvChannelTooLow,
                                    Bool_t& hvChannelTooHigh,
                                    Bool_t& hvChannelON) const
{
  /// Get HV status for one HV sector of St12
  
  /// For a given PCB in a given DE, get the HV status (both the channel
  /// and the switch).
  /// Returns false if hv switch changed during the run.
  
  AliCodeTimerAuto("")
  
  Bool_t error = kFALSE;
  hvChannelTooLow = kFALSE;
  hvChannelTooHigh = kFALSE;
  hvChannelON = kTRUE;

  AliMpHVNamer hvNamer;
  
  TString hvChannel(hvNamer.DCSHVChannelName(detElemId,sector));
  
  TMap* hvMap = fCalibrationData.HV();
  TPair* hvPair = static_cast<TPair*>(hvMap->FindObject(hvChannel.Data()));
  if (!hvPair)
  {
    AliError(Form("Did not find expected alias (%s) for DE %d",
                  hvChannel.Data(),detElemId));  
    error = kTRUE;
  }
  else
  {
    TObjArray* values = static_cast<TObjArray*>(hvPair->Value());
    if (!values)
    {
      AliError(Form("Could not get values for alias %s",hvChannel.Data()));
      error = kTRUE;
    }
    else
    {
      // find out min and max value, and makes a cut
      Float_t hvMin(1E9);
      Float_t hvMax(0);
      TIter next(values);
      AliDCSValue* val;
      
      while ( ( val = static_cast<AliDCSValue*>(next()) ) )
      {
        Float_t hv = val->GetFloat();
        hvMin = TMath::Min(hv,hvMin);
        hvMax = TMath::Max(hv,hvMax);
      }
      
      float lowThreshold = fHVSt12Limits.X();
      float highThreshold = fHVSt12Limits.Y();
            
      if ( hvMin < lowThreshold ) hvChannelTooLow = kTRUE;
      if ( hvMax > highThreshold ) hvChannelTooHigh = kTRUE;
      if ( hvMin < 1 ) hvChannelON = kFALSE;
    }
  }
  
  return error;
}

//_____________________________________________________________________________
Bool_t 
AliMUONPadStatusMaker::HVSt345Status(Int_t detElemId, Int_t pcbIndex,
                                     Bool_t& hvChannelTooLow,
                                     Bool_t& hvChannelTooHigh,
                                     Bool_t& hvChannelON,
                                     Bool_t& hvSwitchON) const
{
  /// For a given PCB in a given DE, get the HV status (both the channel
  /// and the switch).
  /// Returns false if something goes wrong (in particular if 
  /// hv switch changed during the run).
  
  AliCodeTimerAuto("")
  
  Bool_t error = kFALSE;
  hvChannelTooLow = kFALSE;
  hvChannelTooHigh = kFALSE;
  hvSwitchON = kTRUE;
  hvChannelON = kTRUE;
  
  AliMpHVNamer hvNamer;
  
  TString hvChannel(hvNamer.DCSHVChannelName(detElemId));
  
  TMap* hvMap = fCalibrationData.HV();
  
  TPair* hvPair = static_cast<TPair*>(hvMap->FindObject(hvChannel.Data()));
  if (!hvPair)
  {
    AliError(Form("Did not find expected alias (%s) for DE %d",
                  hvChannel.Data(),detElemId));  
    error = kTRUE;
  }
  else
  {
    TObjArray* values = static_cast<TObjArray*>(hvPair->Value());
    if (!values)
    {
      AliError(Form("Could not get values for alias %s",hvChannel.Data()));
      error = kTRUE;
    }
    else
    {
      // find out min and max value, and makes a cut
      Float_t hvMin(1E9);
      Float_t hvMax(0);
      TIter next(values);
      AliDCSValue* val;
      
      while ( ( val = static_cast<AliDCSValue*>(next()) ) )
      {
        Float_t hv = val->GetFloat();
        hvMin = TMath::Min(hv,hvMin);
        hvMax = TMath::Max(hv,hvMax);
      }

      float lowThreshold = fHVSt345Limits.X();
      float highThreshold = fHVSt345Limits.Y();

      if ( hvMin < lowThreshold ) hvChannelTooLow = kTRUE;
      else if ( hvMax > highThreshold ) hvChannelTooHigh = kTRUE;
      if ( hvMin < 1 ) hvChannelON = kFALSE;
    }
  }
  
  TString hvSwitch(hvNamer.DCSHVSwitchName(detElemId,pcbIndex));
  TPair* switchPair = static_cast<TPair*>(hvMap->FindObject(hvSwitch.Data()));
  if (!switchPair)
  {
    AliError(Form("Did not find expected alias (%s) for DE %d PCB %d",
                  hvSwitch.Data(),detElemId,pcbIndex));
    error = kTRUE;
  }
  else
  {
    TObjArray* values = static_cast<TObjArray*>(switchPair->Value());
    if (!values)
    {    
      AliError(Form("Could not get values for alias %s",hvSwitch.Data()));
      error = kTRUE;
    }
    else
    {
      // we'll count the number of ON/OFF for this pad, to insure
      // consistency (i.e. if status changed during the run, we should
      // at least notify this fact ;-) and hope it's not the norm)
      Int_t nTrue(0);
      Int_t nFalse(0);
      TIter next(values);
      AliDCSValue* val;
      
      while ( ( val = static_cast<AliDCSValue*>(next()) ) )
      {
        if ( val->GetBool() )
        {
          ++nTrue;
        }
        else
        {
          ++nFalse;
        }
      }
      
      if ( (nTrue>0 && nFalse>0) )
      {
        AliWarning(Form("Status of HV Switch %s changed during this run nTrue=%d nFalse=%d! Will consider it OFF",
                        hvSwitch.Data(),nTrue,nFalse));
        error = kTRUE;
      }
      
      if ( nFalse ) hvSwitchON = kFALSE;
    }
  }
  return error;
}

//_____________________________________________________________________________
Int_t
AliMUONPadStatusMaker::HVStatus(Int_t detElemId, Int_t manuId) const
{
  /// Get HV status of one manu
  
  AliCodeTimerAuto("")
  
  if ( !fCalibrationData.HV() ) return kMissing;

  Long_t lint = fHV->GetValue(AliMpManuUID::BuildUniqueID(detElemId,manuId));
  
  if ( lint ) 
  {
    return (Int_t)(lint - 1);
  }

  Int_t status(0);
  
  AliMpHVNamer hvNamer;
  
  switch ( AliMpDEManager::GetStationType(detElemId) )
  {
    case AliMp::kStation1:
    case AliMp::kStation2:
    {
      int sector = hvNamer.ManuId2Sector(detElemId,manuId);
      if ( sector >= 0 ) 
      {
        Bool_t hvChannelTooLow, hvChannelTooHigh, hvChannelON;
        Bool_t error = HVSt12Status(detElemId,sector,
                                    hvChannelTooLow,
                                    hvChannelTooHigh,
                                    hvChannelON);
        if ( error ) status |= kHVError;
        if ( hvChannelTooLow ) status |= kHVTooLow;
        if ( hvChannelTooHigh ) status |= kHVTooHigh; 
        if ( !hvChannelON ) status |= kHVChannelOFF;
        // assign this status to all the other manus handled by the same HV channel
        SetHVStatus(detElemId,sector,status);
      }
    }
      break;
    case AliMp::kStation345:
    {
      int pcbIndex = hvNamer.ManuId2PCBIndex(detElemId,manuId);
      if ( pcbIndex >= 0 ) 
      {
        Bool_t hvChannelTooLow, hvChannelTooHigh, hvChannelON,hvSwitchON;
        Bool_t error = HVSt345Status(detElemId,pcbIndex,
                                     hvChannelTooLow,hvChannelTooHigh,
                                     hvChannelON,hvSwitchON);
        if ( error ) status |= kHVError;
        if ( hvChannelTooLow ) status |= kHVTooLow;
        if ( hvChannelTooHigh ) status |= kHVTooHigh; 
        if ( !hvSwitchON ) status |= kHVSwitchOFF; 
        if ( !hvChannelON) status |= kHVChannelOFF;
        // assign this status to all the other manus handled by the same HV channel
        SetHVStatus(detElemId,pcbIndex,status);
      }
    }
      break;
    default:
      break;
  }
  
  return status;
}

//_____________________________________________________________________________
AliMUONVCalibParam* 
AliMUONPadStatusMaker::Neighbours(Int_t detElemId, Int_t manuId) const
{
  /// Get the neighbours parameters for a given manu
  AliMUONVStore* neighbourStore = fCalibrationData.Neighbours();
  return static_cast<AliMUONVCalibParam*>(neighbourStore->FindObject(detElemId,manuId));
}

//_____________________________________________________________________________
AliMUONVStore* 
AliMUONPadStatusMaker::NeighboursStore() const
{
  /// Return the store containing all the neighbours
  return fCalibrationData.Neighbours();
}

//_____________________________________________________________________________
AliMUONVCalibParam*
AliMUONPadStatusMaker::ComputeStatus(Int_t detElemId, Int_t manuId) const
{
  /// Compute the status of a given manu, using all available information,
  /// i.e. pedestals, gains, and HV
  
//  AliCodeTimerAuto("")
  
//  AliCodeTimerStart("Param creation");
  AliMUONVCalibParam* param = new AliMUONCalibParamNI(1,AliMpConstants::ManuNofChannels(),detElemId,manuId,-1);
  fStatus->Add(param);
//  AliCodeTimerStop("Param creation");
  
//  AliCodeTimerStart("FindObject");
  AliMUONVCalibParam* pedestals = static_cast<AliMUONVCalibParam*>(fPedestals->FindObject(detElemId,manuId));

  AliMUONVCalibParam* gains = static_cast<AliMUONVCalibParam*>(fGains->FindObject(detElemId,manuId));
//  AliCodeTimerStop("FindObject");
  
  Int_t hvStatus = HVStatus(detElemId,manuId);

//  AliCodeTimerStart("Loop");
  
  for ( Int_t manuChannel = 0; manuChannel < param->Size(); ++manuChannel )
  {
    Int_t pedStatus(0);
    
    if (pedestals) 
    {
      Float_t pedMean = pedestals->ValueAsFloatFast(manuChannel,0);
      Float_t pedSigma = pedestals->ValueAsFloatFast(manuChannel,1);
      if ( pedMean < fPedMeanLimits.X() ) pedStatus |= kPedMeanTooLow;
      else if ( pedMean > fPedMeanLimits.Y() ) pedStatus |= kPedMeanTooHigh;
      if ( pedSigma < fPedSigmaLimits.X() ) pedStatus |= kPedSigmaTooLow;
      else if ( pedSigma > fPedSigmaLimits.Y() ) pedStatus |= kPedSigmaTooHigh;
      if ( pedMean == 0 ) pedStatus |= kPedMeanZero;
    }
    else
    {
      pedStatus = kPedMissing;
    }
    
    Int_t gainStatus(0);
  
    if ( gains ) 
    {
      Float_t a0 = gains->ValueAsFloatFast(manuChannel,0);
      Float_t a1 = gains->ValueAsFloatFast(manuChannel,1);
      Float_t thres = gains->ValueAsFloatFast(manuChannel,2);
  
      if ( a0 < fGainA0Limits.X() ) gainStatus |= kGainA0TooLow;
      else if ( a0 > fGainA0Limits.Y() ) gainStatus |= kGainA0TooHigh;
      if ( a1 < fGainA1Limits.X() ) gainStatus |= kGainA1TooLow;
      else if ( a1 > fGainA1Limits.Y() ) gainStatus |= kGainA1TooHigh;
      if ( thres < fGainThresLimits.X() ) gainStatus |= kGainThresTooLow;
      else if ( thres > fGainThresLimits.Y() ) gainStatus |= kGainThresTooHigh;
    }
    else
    {
      gainStatus = kGainMissing;
    }
        
    Int_t status = BuildStatus(pedStatus,hvStatus,gainStatus);
      
    param->SetValueAsIntFast(manuChannel,0,status);
  }
  
//  AliCodeTimerStop("Loop");
  
  return param;
}

//_____________________________________________________________________________
AliMUONVCalibParam* 
AliMUONPadStatusMaker::PadStatus(Int_t detElemId, Int_t manuId) const
{
  /// Get the status for a given channel
  
  AliMUONVCalibParam* param = static_cast<AliMUONVCalibParam*>(fStatus->FindObject(detElemId,manuId));
  if (!param)
  {
    // not already there, so compute it now
    AliCodeTimerAuto("ComputeStatus");
    param = ComputeStatus(detElemId,manuId);
  }
  return param;
}

//_____________________________________________________________________________
Int_t 
AliMUONPadStatusMaker::PadStatus(Int_t detElemId, Int_t manuId, Int_t manuChannel) const
{
  /// Get the status for a given channel
  
  AliMUONVCalibParam* param = static_cast<AliMUONVCalibParam*>(fStatus->FindObject(detElemId,manuId));
  if (!param)
  {
    // not already there, so compute it now
    param = ComputeStatus(detElemId,manuId);
  }
  return param->ValueAsInt(manuChannel,0);
}

//_____________________________________________________________________________
void
AliMUONPadStatusMaker::SetHVStatus(Int_t detElemId, Int_t index, Int_t status) const
{
  /// Assign status to all manus in a given HV "zone" (defined by index, meaning
  /// is different thing from St12 and St345)
  
  AliCodeTimerAuto("")
  
  AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
  
  const AliMpArrayI* manus = de->ManusForHV(index);
  
  for ( Int_t i = 0; i < manus->GetSize(); ++ i ) 
  {
    Int_t manuId = manus->GetValue(i);
    fHV->Add(AliMpManuUID::BuildUniqueID(detElemId,manuId),status + 1);
  }
}
