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

#include "AliQA.h"

#include "AliMUON2DMap.h"
#include "AliMUON2DStoreValidator.h"
#include "AliMUONCalibParamNI.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONStringIntMap.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVTrackerData.h"

#include "AliMpArea.h"
#include "AliMpArrayI.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
#include "AliMpDEManager.h"
#include "AliMpDetElement.h"
#include "AliMpDCSNamer.h"
#include "AliMpManuUID.h"

#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCodeTimer.h"
#include "AliDCSValue.h"
#include "AliLog.h"

#include <Riostream.h>
#include <TArrayI.h>
#include <TExMap.h>
#include <TFile.h>
#include <TKey.h>
#include <TMap.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>

/// \cond CLASSIMP
ClassImp(AliMUONPadStatusMaker)
/// \endcond

//_____________________________________________________________________________
AliMUONPadStatusMaker::AliMUONPadStatusMaker(const AliMUONCalibrationData& calibData)
: fkCalibrationData(calibData),
fGainA1Limits(0,1E30),
fGainA2Limits(-1E-30,1E30),
fGainThresLimits(0,4095),
fHVSt12Limits(0,5000),
fHVSt345Limits(0,5000),
fPedMeanLimits(0,4095),
fPedSigmaLimits(0,4095),
fManuOccupancyLimits(0,0.1),
fStatus(new AliMUON2DMap(true)),
fHV(new TExMap),
fPedestals(calibData.Pedestals()),
fGains(calibData.Gains()),
fTrackerData(0x0)
{
  /// ctor
  AliDebug(1,Form("ped store %s gain store %s",
                  fPedestals->ClassName(),
                  fGains->ClassName()));
  
  TString qaFileName(AliQA::GetQADataFileName("MUON",calibData.RunNumber()));
  
  // search the QA file in memory first.
  TFile* f = static_cast<TFile*>(gROOT->GetListOfFiles()->FindObject(qaFileName.Data()));

  if (!f)
  {
    // then tries to open it
    if ( gSystem->AccessPathName(qaFileName.Data()) == kFALSE ) 
    {
      f = TFile::Open(qaFileName.Data());
      if ( f )
      {
        AliDebug(1,Form("Got %s from disk",qaFileName.Data()));
      }
    }
  }
  else
  {
    AliDebug(1,Form("Got %s from memory",qaFileName.Data()));
  }
  
  if (f)
  {
    TDirectory* d = gDirectory;
    
    f->cd("MUON/Raws");
    
    TIter next(gDirectory->GetListOfKeys());
    TKey* key;
    
    while ( ( key = static_cast<TKey*>(next()) ) && !fTrackerData )
    {
      TString name(key->GetName());
      
      if ( name.Contains("CALZ") )
      {
        fTrackerData = dynamic_cast<AliMUONVTrackerData*>(key->ReadObj());
      }
      
    }
    
    gDirectory = d;
    
    if ( fTrackerData ) 
    {
      AliInfo(Form("Will make a cut on MANU occupancy from TrackerData=%s",fTrackerData->GetName()));
    }
    else
    {
      AliWarning(Form("Found a QA file = %s, but could not get the expected TrackerData in there... (probably not a serious problem though)",
                      f->GetName()));
    }
  }
  else
  {
    AliWarning("Did not find QA file, so will not use manu occupancy as a criteria");
  }
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
  Int_t otherStatus;
  
  DecodeStatus(status,pedStatus,hvStatus,gainStatus,otherStatus);
  
  TString s;
  
  if ( pedStatus & kPedMeanZero ) s += "& Ped Mean is Zero ";
  if ( pedStatus & kPedMeanTooLow ) s += "& Ped Mean Too Low ";
  if ( pedStatus & kPedMeanTooHigh ) s += "& Ped Mean Too High ";
  if ( pedStatus & kPedSigmaTooLow ) s += "& Ped Sigma Too Low ";
  if ( pedStatus & kPedSigmaTooHigh ) s += "& Ped Sigma Too High ";
  if ( pedStatus & kPedMissing ) s += "& Ped is missing ";
  
	if ( gainStatus & kGainA1TooLow ) s+="& Gain A1 is Too Low ";
	if ( gainStatus & kGainA1TooHigh ) s+="& Gain A1 is Too High ";
	if ( gainStatus & kGainA2TooLow ) s+="& Gain A2 is Too Low ";
	if ( gainStatus & kGainA2TooHigh ) s+="& Gain A2 is Too High ";
	if ( gainStatus & kGainThresTooLow ) s+="& Gain Thres is Too Low ";
	if ( gainStatus & kGainThresTooHigh ) s+="& Gain Thres is Too High ";
	if ( gainStatus & kGainMissing ) s+="& Gain is missing ";
	
	if ( hvStatus & kHVError ) s+="& HV is on error ";
	if ( hvStatus & kHVTooLow ) s+="& HV is Too Low ";
	if ( hvStatus & kHVTooHigh ) s+="& HV is Too High ";
	if ( hvStatus & kHVChannelOFF ) s+="& HV has channel OFF ";
	if ( hvStatus & kHVSwitchOFF ) s+="& HV has switch OFF ";
	if ( hvStatus & kHVMissing ) s+="& HV is missing ";

  if ( otherStatus & kManuOccupancyTooHigh ) s+="& manu occupancy too high ";
  if ( otherStatus & kManuOccupancyTooLow ) s+="& manu occupancy too low ";
  
  if ( s[0] == '&' ) s[0] = ' ';
  
  return s;
}

//_____________________________________________________________________________
TString
AliMUONPadStatusMaker::AsCondition(Int_t mask)
{
  /// return a human readable version of the mask's equivalent condition
  
  TString s(AsString(mask));
  
  s.ReplaceAll("&","|");
  
  return s;
}

//_____________________________________________________________________________
Int_t
AliMUONPadStatusMaker::BuildStatus(Int_t pedStatus, 
                                   Int_t hvStatus, 
                                   Int_t gainStatus,
                                   Int_t otherStatus)
{
  /// Build a complete status from specific parts (ped,hv,gain)
  
  return ( hvStatus & 0xFF ) | ( ( pedStatus & 0xFF ) << 8 ) | 
  ( ( gainStatus & 0xFF ) << 16 ) |
  ( ( otherStatus & 0xFF ) << 24 ) ;
}

//_____________________________________________________________________________
void
AliMUONPadStatusMaker::DecodeStatus(Int_t status, 
                                    Int_t& pedStatus, 
                                    Int_t& hvStatus, 
                                    Int_t& gainStatus,
                                    Int_t& otherStatus)
{
  /// Decode complete status into specific parts (ped,hv,gain)
  
  otherStatus = ( status & 0xFF000000 ) >> 24;
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

  AliMpDCSNamer hvNamer("TRACKER");
  
  TString hvChannel(hvNamer.DCSChannelName(detElemId,sector));
  
  TMap* hvMap = fkCalibrationData.HV();
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
  
  AliMpDCSNamer hvNamer("TRACKER");
  
  TString hvChannel(hvNamer.DCSChannelName(detElemId));
  
  TMap* hvMap = fkCalibrationData.HV();
  
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
  
  TString hvSwitch(hvNamer.DCSSwitchName(detElemId,pcbIndex));
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
  
  if ( !fkCalibrationData.HV() ) return kMissing;

  Long_t lint = fHV->GetValue(AliMpManuUID::BuildUniqueID(detElemId,manuId));
  
  if ( lint ) 
  {
    return (Int_t)(lint - 1);
  }

  Int_t status(0);
  
  AliMpDCSNamer hvNamer("TRACKER");
  
  switch ( AliMpDEManager::GetStationType(detElemId) )
  {
    case AliMp::kStation12:
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
  AliMUONVStore* neighbourStore = fkCalibrationData.Neighbours();
  return static_cast<AliMUONVCalibParam*>(neighbourStore->FindObject(detElemId,manuId));
}

//_____________________________________________________________________________
AliMUONVStore* 
AliMUONPadStatusMaker::NeighboursStore() const
{
  /// Return the store containing all the neighbours
  return fkCalibrationData.Neighbours();
}

//_____________________________________________________________________________
AliMUONVCalibParam*
AliMUONPadStatusMaker::ComputeStatus(Int_t detElemId, Int_t manuId) const
{
  /// Compute the status of a given manu, using all available information,
  /// i.e. pedestals, gains, and HV
  
  AliMUONVCalibParam* param = new AliMUONCalibParamNI(1,AliMpConstants::ManuNofChannels(),detElemId,manuId,-1);
  fStatus->Add(param);

  AliMUONVCalibParam* pedestals = static_cast<AliMUONVCalibParam*>(fPedestals->FindObject(detElemId,manuId));

  AliMUONVCalibParam* gains = static_cast<AliMUONVCalibParam*>(fGains->FindObject(detElemId,manuId));
  
  Int_t hvStatus = HVStatus(detElemId,manuId);

  Int_t otherStatus = OtherStatus(detElemId,manuId);
  
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
  
      if ( a0 < fGainA1Limits.X() ) gainStatus |= kGainA1TooLow;
      else if ( a0 > fGainA1Limits.Y() ) gainStatus |= kGainA1TooHigh;
      if ( a1 < fGainA2Limits.X() ) gainStatus |= kGainA2TooLow;
      else if ( a1 > fGainA2Limits.Y() ) gainStatus |= kGainA2TooHigh;
      if ( thres < fGainThresLimits.X() ) gainStatus |= kGainThresTooLow;
      else if ( thres > fGainThresLimits.Y() ) gainStatus |= kGainThresTooHigh;
    }
    else
    {
      gainStatus = kGainMissing;
    }
        
    Int_t status = BuildStatus(pedStatus,hvStatus,gainStatus,otherStatus);
      
    param->SetValueAsIntFast(manuChannel,0,status);
  }
  
  return param;
}

//_____________________________________________________________________________
Int_t 
AliMUONPadStatusMaker::OtherStatus(Int_t detElemId, Int_t manuId) const
{
  /// Get the "other" status for a given manu
  if ( fTrackerData ) 
  {
    Double_t occ = fTrackerData->Manu(detElemId,manuId,2);
    if ( occ < fManuOccupancyLimits.X() )
    {
      return kManuOccupancyTooLow;
    }
    if ( occ > fManuOccupancyLimits.Y() )
    {
      return kManuOccupancyTooHigh;
    }
  }
  return 0;
}

//_____________________________________________________________________________
AliMUONVCalibParam* 
AliMUONPadStatusMaker::PadStatus(Int_t detElemId, Int_t manuId) const
{
  /// Get the status container for a given manu
  
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
